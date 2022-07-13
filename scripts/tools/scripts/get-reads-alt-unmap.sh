#!/bin/bash
# Adapted from version on xHLA Github
set -eu -o pipefail

[[ -z $@ ]] && echo \
"
Filter reads from a .bam file for HLA typing.\n
Usage: ./get-reads-alt-unmap.sh <input.bam> <output.bam>"

IN="$1"
OUT="$2"

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    NCORE=$(grep -c '^processor' /proc/cpuinfo)
    # POSIX compliant system memory query
    MEM=$(awk 'BEGIN{for (i=1; i<ARGC;i++)
       printf "%.0fGB\n", ARGV[i]}' $(echo "$(grep MemTotal /proc/meminfo | awk '{print $2}') * 0.95 / 1000 / 1000" | bc))
elif [[ "$OSTYPE" == darwin* ]]; then
    NCORE=$(sysctl -n hw.ncpu)
    MEM=8GB
fi

echo "NCORE=$NCORE MEM=$MEM"

BIN="`dirname \"$0\"`"

cleanup(){

  rm "${OUT}.1.fq" || true
  rm "${OUT}.2.fq" || true
  rm "${OUT}.full.bam" || true
  rm "${OUT}.full.bam.bai" || true

}

trap cleanup EXIT

# Compared to the original xHLA preprocessing script,
# we remove the code filtering by position of the reads.
# Other tools than xHLA might use regions outside chr6:29844528-33100696
# (as this range does not include all MHC-genes)
sambamba view \
  -f "bam" -h -p -l 0 -t $NCORE \
  "$IN" | sambamba sort -p -n -t $NCORE -o - /dev/stdin |
  samtools fastq -@ 32 /dev/stdin \
  -1 "${OUT}.1.fq" \
  -2 "${OUT}.2.fq" \
  -0 /dev/null -s /dev/null -n

# Note: reference genome was changed here to "GCA_000001405.15_GRCh38_no_alt_analysis_set"
# xHLA used a sliced version of a non-alt aware reference.
# As the origin of this file is unclear, we use the primary reference used by the GDC portal.
bwa mem -t $NCORE "/home/labgroups/ccgg/downloads/genomes/RefSeq/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" "${OUT}.1.fq" "${OUT}.2.fq" |
samtools view -b - |
sambamba sort -m $MEM -t $NCORE -o "${OUT}.full.bam" /dev/stdin

sambamba index "${OUT}.full.bam"

# Use the entire primary chromosome 6 to make sure nothing is missing
# (input files were sliced during download anyway)
samtools view -o "$OUT" -b "${OUT}.full.bam" chr6 '*'
sambamba index "$OUT"