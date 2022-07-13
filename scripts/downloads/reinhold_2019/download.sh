#!/bin/bash
# Remark: Make sure build_star_index.sh was run before starting this script
# bash scripts/downloads/reinhold_2019/build_star_index.sh

# bash scripts/downloads/reinhold_2019/download_meta.sh

# Download to downloads/pub/reinhold_2019/full as fastq.gz files, subfolder per SRR
xargs -L1 bash scripts/downloads/reinhold_2019/download_single.sh < \
scripts/downloads/reinhold_2019/files.txt

# Run scripts/downloads/reinhold_2019/run_trim_align_rna.sh on all SRR IDs
# Output files go into downloads/pub/reinhold_2019/HLA
mkdir downloads/pub/reinhold_2019/HLA

for SRR in $(ls downloads/pub/reinhold_2019/full); do
  scripts/downloads/reinhold_2019/run_trim_align_rna.sh "$SRR"
done