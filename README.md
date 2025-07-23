# scLongReadPreprocessor

**A Python tool for extracting and correcting barcodes and UMIs from third-generation single-cell sequencing data with TSO correction and polyA/T tail trimming.**

---

## Features

- Supports gzipped FASTQ input (3rd-gen scRNA-seq long reads)
- Extracts barcode and UMI sequences based on customizable barcode & UMI lengths
- Detects and corrects barcodes with fuzzy matching using BK-tree or edit distance
- Detects TSO (template switching oligo) sequence with user-defined max edit distance using `edlib`
- Automatically trims polyA or polyT tails from cDNA sequence ends
- Multi-threaded processing for fast performance on large datasets
- Outputs three FASTQ files for filtered, unmatched TSO, and invalid barcode reads
- Generates barcode count summary

---

## Installation

```bash
git clone https://github.com/yourusername/scLongReadPreprocessor.git
cd scLongReadPreprocessor
pip install -r requirements.txt  # install dependencies: edlib, editdistance, tqdm
```

## Usage

```
python extract_barcode.py \
  --fastq input_reads.fastq.gz \
  --whitelist barcodes.txt \
  --out-prefix sample1 \
  --out-dir ./results \
  --bc-len 16 \
  --umi-len 12 \
  --tso "TTTCTTATATGGG" \
  --tso-max-dist 2 \
  --bc-max-dist 1 \
  --threads 8 \
  --chunk-size 100000 \
  --use-bktree \
  --poly-tail-min-len 6
```

### Parameters

- `--fastq` : Input gzipped FASTQ file of long reads
- `--whitelist` : File containing barcode whitelist, one barcode per line
- `--bc-len` : Barcode length (int)
- `--umi-len` : UMI length (int)
- `--tso` : Template switching oligo (TSO) sequence (default: TTTCTTATATGGG)
- `--tso-max-dist` : Max edit distance allowed for TSO matching (default: 2)
- `--bc-max-dist` : Max edit distance allowed for barcode correction (default: 1)
- `--use-bktree` : Use BK-tree data structure for faster barcode fuzzy matching
- `--poly-tail-min-len` : Minimum length to trim polyA/T tails from sequence ends (default: 6)
- `--threads` : Number of CPU threads to run in parallel
- `--chunk-size` : Number of reads to process per chunk

## Output

- `{out-prefix}.filtered.fastq.gz` : Reads with matched TSO and corrected barcode+UMI, with modified header
- `{out-prefix}.unmatched_tso.fastq.gz` : Reads that failed TSO matching
- `{out-prefix}.invalid_barcode.fastq.gz` : Reads with invalid or unmatched barcodes/UMI
- `{out-prefix}.barcode_summary.tsv` : Tab-delimited summary of barcode counts in filtered reads

## Notes

- The script requires Python3 and following packages: `edlib`, `editdistance`, `tqdm`
- Designed for third-generation single-cell long read data (PacBio/ONT)
- The barcode whitelist file should include all expected barcode sequences (one per line)
- Adjust parameters such as TSO sequence and max distances to fit your experimental design



