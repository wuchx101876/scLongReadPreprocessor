#!/usr/bin/env python3
import pysam
from intervaltree import IntervalTree
from collections import defaultdict
import argparse

def load_gene_bed(gene_bed_file):
    """加载基因区间，建立染色体->IntervalTree映射"""
    gene_intervals = defaultdict(IntervalTree)
    with open(gene_bed_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split('\t')
            if len(cols) < 4:
                continue
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            gene_id = cols[3]
            gene_intervals[chrom].addi(start, end, gene_id)
    return gene_intervals

def get_read_genes(read, gene_intervals):
    """
    返回 read 比对区间(s)对应所有重叠gene的 set
    推荐用 get_blocks() 获取连续比对区间，避免简单的start-end不准确
    """
    if read.is_unmapped:
        return set()
    chrom = read.reference_name
    if chrom not in gene_intervals:
        return set()

    ivtree = gene_intervals[chrom]
    genes = set()
    try:
        blocks = read.get_blocks()  # list of (start,end)
    except:
        # 保险，返回空
        return genes

    for start, end in blocks:
        overlaps = ivtree.overlap(start, end)
        for iv in overlaps:
            genes.add(iv.data)
    return genes

def main(args):
    gene_intervals = load_gene_bed(args.gene_bed)
    bamfile = pysam.AlignmentFile(args.input_bam, "rb")

    # 三重字典：counts[cell][gene][UMI] = True， 用于去重
    counts = defaultdict(lambda: defaultdict(set))

    cnt_reads = 0
    cnt_annotated = 0

    for read in bamfile.fetch(until_eof=True):
        cnt_reads += 1
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        # 读出CB和UB标签
        try:
            cb = read.get_tag("CB")
            ub = read.get_tag("UB")
        except KeyError:
            # 没有条码或UMI不考虑
            continue

        genes = get_read_genes(read, gene_intervals)
        if not genes:
            continue

        # 这里假设read可能匹配多个基因，全部统统计入
        for gene in genes:
            counts[cb][gene].add(ub)
        cnt_annotated += 1

        if cnt_reads % 100000 == 0:
            print(f"处理reads: {cnt_reads}; 已注释基因reads: {cnt_annotated}")

    bamfile.close()

    # 输出count矩阵，格式: cell_gene_umi_counts.tsv
    # 行: cell barcode
    # 列: gene
    # 内容: UMI计数
    print("输出计数矩阵到文件:", args.output)
    # 收集所有gene和cell列表，用于构建矩阵
    all_cells = sorted(counts.keys())
    all_genes = sorted({gene for cell in counts for gene in counts[cell]})

    with open(args.output, "w") as out:
        # 写表头
        out.write("Cell\\Gene\t" + "\t".join(all_genes) + "\n")
        for cell in all_cells:
            out.write(cell)
            for gene in all_genes:
                umi_count = len(counts[cell][gene]) if gene in counts[cell] else 0
                out.write(f"\t{umi_count}")
            out.write("\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="三代单细胞测序BAM统计基因count矩阵")
    parser.add_argument("-i", "--input_bam", required=True, help="输入BAM文件，排序并索引")
    parser.add_argument("-g", "--gene_bed", required=True, help="基因注释BED文件，格式chr start end gene_id")
    parser.add_argument("-o", "--output", required=True, help="输出count矩阵TSV文件")
    args = parser.parse_args()

    main(args)