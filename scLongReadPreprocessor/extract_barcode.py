#!/usr/bin/env python3
import os
import gzip
import argparse
import editdistance
import edlib
from tqdm import tqdm
from collections import Counter
from multiprocessing import Pool
import logging
import sys

# ------------------------
# BK-Tree for barcode correction
# ------------------------
class BKTree:
    def __init__(self, dist_fn):
        self.dist_fn = dist_fn
        self.tree = None

    def add(self, word):
        node = self.tree
        if node is None:
            self.tree = (word, {})
            return
        while True:
            root, children = node
            d = self.dist_fn(word, root)
            if d in children:
                node = children[d]
            else:
                children[d] = (word, {})
                break

    def search(self, word, max_dist):
        results = []
        node = self.tree
        if node is None:
            return results
        candidates = [self.tree]
        while candidates:
            root, children = candidates.pop()
            d = self.dist_fn(word, root)
            if d <= max_dist:
                results.append((root, d))
            candidates += [c for dist, c in children.items() if d - max_dist <= dist <= d + max_dist]
        return results

# ------------------------
# FASTQ Utilities
# ------------------------
def read_fastq_chunk(fq_handle, chunk_size=100000):
    chunk = []
    for _ in range(chunk_size):
        lines = [fq_handle.readline() for _ in range(4)]
        if not lines[0]:
            break
        # 检查是否完整4行
        if any(line == '' for line in lines):
            logging.warning("Incomplete FASTQ record detected, skipping.")
            break
        chunk.append(tuple(line.rstrip('\n') for line in lines))
    return chunk

# 写fastq，确保4行且格式正确
def write_fastq(path, records):
    try:
        with gzip.open(path, 'wt') as f:
            for r in records:
                if len(r) == 4:
                    f.write('\n'.join(r) + '\n')
                else:
                    logging.warning(f"Record skipped because it doesn't have 4 lines: {r}")
    except Exception as e:
        logging.error(f"Failed to write fastq file {path}: {e}")
        raise

# ------------------------
# polyA/polyT尾巴剪切功能
# ------------------------
def auto_trim_poly_tail(seq, qual, min_len=6):
    """
    自动检测并剪切末端polyT或polyA连续区，长度≥min_len
    """
    for base in ['T', 'A']:
        tail_len = 0
        for b in reversed(seq):
            if b == base:
                tail_len += 1
            else:
                break
        if tail_len >= min_len:
            return seq[:-tail_len], qual[:-tail_len]
    return seq, qual

# ------------------------
# 核心reads处理函数
# ------------------------
def process_reads(chunk, args, whitelist_set, bk_tree=None):
    filtered = []
    unmatched_tso = []
    invalid_bc = []
    bc_counter = Counter()

    total_len = args.bc_len + args.umi_len

    for header, seq, sep, qual in chunk:
        try:
            tso_result = edlib.align(args.tso, seq, mode="HW", task="locations", k=args.tso_max_dist)
            if tso_result['editDistance'] == -1 or not tso_result.get('locations'):
                unmatched_tso.append((header, seq, sep, qual))
                continue
            tso_start = tso_result['locations'][0][0]

            if tso_start < total_len:
                invalid_bc.append((header, seq, sep, qual))
                continue

            bc_umi_seq = seq[tso_start - total_len: tso_start]
            bc_seq = bc_umi_seq[:args.bc_len]
            umi_seq = bc_umi_seq[args.bc_len:]

            # Barcode纠错
            if args.use_bktree and bk_tree is not None:
                matches = bk_tree.search(bc_seq, args.bc_max_dist)
                if not matches:
                    corrected_bc = None
                else:
                    matches_sorted = sorted(matches, key=lambda x: x[1])
                    min_dist = matches_sorted[0][1]
                    first_match = [m for m in matches_sorted if m[1] == min_dist][0]
                    corrected_bc = first_match[0]
            else:
                matches = []
                bc_len = len(bc_seq)
                for bc in whitelist_set:
                    if abs(len(bc) - bc_len) <= args.bc_max_dist:
                        dist = editdistance.eval(bc_seq, bc)
                        if dist <= args.bc_max_dist:
                            matches.append((bc, dist))
                if not matches:
                    corrected_bc = None
                else:
                    matches_sorted = sorted(matches, key=lambda x: x[1])
                    min_dist = matches_sorted[0][1]
                    first_match = [m for m in matches_sorted if m[1] == min_dist][0]
                    corrected_bc = first_match[0]

            if corrected_bc is None:
                invalid_bc.append((header, seq, sep, qual))
                continue

            cDNA_seq = seq[tso_start + len(args.tso):]
            cDNA_qual = qual[tso_start + len(args.tso):]

            # 剪polyA/polyT尾巴，使用新参数
            cDNA_seq, cDNA_qual = auto_trim_poly_tail(cDNA_seq, cDNA_qual, min_len=args.poly_tail_min_len)

            if len(cDNA_seq) != len(cDNA_qual):
                logging.warning(f"Seq and qual length mismatch after trimming for read {header}")
                invalid_bc.append((header, seq, sep, qual))
                continue
            # 如果cDNA太短，也忽略
            if len(cDNA_seq) == 0:
                invalid_bc.append((header, seq, sep, qual))
                continue

            first_field = header.split()[0]  # 只保留第一个字段
            new_header = f"{first_field}\tCB:Z:{corrected_bc}\tUB:Z:{umi_seq}"
            filtered.append((new_header, cDNA_seq, sep, cDNA_qual))
            bc_counter[corrected_bc] += 1

        except Exception as e:
            logging.error(f"Error processing read {header}: {e}")
            continue

    return filtered, unmatched_tso, invalid_bc, bc_counter

# ------------------------
# 参数解析
# ------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Extract barcode + UMI with TSO correction and polyA/T trimming")
    parser.add_argument("--fastq", required=True, help="Input FASTQ file (gz)")
    parser.add_argument("--whitelist", required=True, help="Barcode whitelist")
    parser.add_argument("--out-prefix", default="output", help="Output file prefix")
    parser.add_argument("--out-dir", default=".", help="Output directory")
    parser.add_argument("--bc-len", type=int, required=True, help="Barcode length")
    parser.add_argument("--umi-len", type=int, required=True, help="UMI length")
    parser.add_argument("--tso", default="TTTCTTATATGGG", help="TSO sequence to match")
    parser.add_argument("--tso-max-dist", type=int, default=2, help="Max TSO edit distance")
    parser.add_argument("--bc-max-dist", type=int, default=1, help="Max barcode correction distance")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--chunk-size", type=int, default=100000, help="Reads per chunk")
    parser.add_argument("--use-bktree", action='store_true', help="Use BK-tree for barcode matching")
    parser.add_argument("--poly-tail-min-len", type=int, default=6,
                        help="Minimum length to consider for trimming polyA/T tail (default: 6)")  # 新增参数
    return parser.parse_args()

# ------------------------
# 全局变量初始化
# ------------------------
global_args = None
global_whitelist_set = None
global_bk_tree = None

def init_worker(args, whitelist_set, bk_tree):
    global global_args
    global global_whitelist_set
    global global_bk_tree
    global_args = args
    global_whitelist_set = whitelist_set
    global_bk_tree = bk_tree

def worker_process(chunk):
    return process_reads(chunk, global_args, global_whitelist_set, global_bk_tree)

# ------------------------
# 主程序入口
# ------------------------
def main():
    args = parse_args()

    try:
        os.makedirs(args.out_dir, exist_ok=True)
    except Exception as e:
        print(f"Error creating output directory {args.out_dir}: {e}", file=sys.stderr)
        sys.exit(1)

    log_path = os.path.join(args.out_dir, args.out_prefix + ".log")
    logging.basicConfig(filename=log_path, level=logging.INFO,
                        format='[%(asctime)s] %(levelname)s: %(message)s')

    try:
        with open(args.whitelist) as wl_file:
            whitelist = [line.strip() for line in wl_file if line.strip()]
        whitelist_set = set(whitelist)
    except Exception as e:
        logging.error(f"Failed to load whitelist file {args.whitelist}: {e}")
        print(f"Failed to load whitelist file {args.whitelist}", file=sys.stderr)
        sys.exit(1)

    logging.info(f"Loaded whitelist with {len(whitelist)} barcodes")

    if args.use_bktree:
        logging.info("Building BK-tree for whitelist...")
        bk_tree = BKTree(lambda x, y: editdistance.eval(x, y))
        for bc in whitelist:
            bk_tree.add(bc)
        logging.info("BK-tree build complete")
    else:
        bk_tree = None

    filtered_path = os.path.join(args.out_dir, args.out_prefix + ".filtered.fastq.gz")
    unmatched_path = os.path.join(args.out_dir, args.out_prefix + ".unmatched_tso.fastq.gz")
    invalid_path = os.path.join(args.out_dir, args.out_prefix + ".invalid_barcode.fastq.gz")
    barcodesummary_path = os.path.join(args.out_dir, args.out_prefix + ".barcode_summary.tsv")

    bc_counter = Counter()
    total_filtered = 0
    total_unmatched = 0
    total_invalid = 0

    try:
        fq_in = gzip.open(args.fastq, 'rt')
    except Exception as e:
        logging.error(f"Failed to open FASTQ file {args.fastq}: {e}")
        print(f"Failed to open FASTQ file {args.fastq}", file=sys.stderr)
        sys.exit(1)

    pool = Pool(processes=args.threads, initializer=init_worker, initargs=(args, whitelist_set, bk_tree))

    try:
        filtered_fq = gzip.open(filtered_path, 'wt')
        unmatched_fq = gzip.open(unmatched_path, 'wt')
        invalid_fq = gzip.open(invalid_path, 'wt')
    except Exception as e:
        logging.error(f"Failed to open output files: {e}")
        print(f"Failed to open output files", file=sys.stderr)
        sys.exit(1)

    total_reads = 0
    jobs = []

    try:
        while True:
            chunk = read_fastq_chunk(fq_in, args.chunk_size)
            if not chunk:
                break
            total_reads += len(chunk)
            jobs.append(pool.apply_async(worker_process, (chunk,)))
    except Exception as e:
        logging.error(f"Error reading input FASTQ: {e}")
        print("Error reading input FASTQ", file=sys.stderr)
        sys.exit(1)

    logging.info(f"Total reads to process (estimated): {total_reads}")
    try:
        with tqdm(total=total_reads, desc="Processing reads") as pbar:
            for job in jobs:
                filtered, unmatched, invalid, bc_c = job.get()
                total_filtered += len(filtered)
                total_unmatched += len(unmatched)
                total_invalid += len(invalid)
                bc_counter.update(bc_c)

                for r in filtered:
                    if len(r) == 4:
                        filtered_fq.write('\n'.join(r) + '\n')
                for r in unmatched:
                    if len(r) == 4:
                        unmatched_fq.write('\n'.join(r) + '\n')
                for r in invalid:
                    if len(r) == 4:
                        invalid_fq.write('\n'.join(r) + '\n')

                pbar.update(len(filtered)+len(unmatched)+len(invalid))
    except Exception as e:
        logging.error(f"Error during processing/writing output: {e}")
        print("Error during processing/writing output", file=sys.stderr)
        sys.exit(1)
    finally:
        fq_in.close()
        filtered_fq.close()
        unmatched_fq.close()
        invalid_fq.close()
        pool.close()
        pool.join()

    logging.info(f"Reads matched TSO and barcode: {total_filtered}")
    logging.info(f"Reads with unmatched TSO: {total_unmatched}")
    logging.info(f"Reads with invalid barcode or UMI length: {total_invalid}")
    logging.info(f"Sum reads in categories: {total_filtered + total_unmatched + total_invalid}")

    try:
        with open(barcodesummary_path, 'w') as f:
            f.write("Barcode\tCount\n")
            for bc, count in bc_counter.most_common():
                f.write(f"{bc}\t{count}\n")
    except Exception as e:
        logging.error(f"Failed to write barcode summary: {e}")
        print("Failed to write barcode summary", file=sys.stderr)

    logging.info("Finished successfully.")
    print(f"✅ All done! Outputs written with prefix '{args.out_prefix}' in folder '{args.out_dir}'")

if __name__ == "__main__":
    main()