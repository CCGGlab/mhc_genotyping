#!/bin/bash
if [ -z "$bam_input" ]; then
    echo "The bam_input must be defined" >&2
    exit 1    
fi

if [ -z "$sample_id" ]; then
    echo "The sample_id must be defined" >&2
    exit 1    
fi

if [ -z "$temp_folder" ]; then
    echo "A temp folder must be defined" >&2
    exit 1
fi

# FASTQ output files:
# Not all tools support compression.
# To enable compression: set extension to .fastq.gz
fastq1="${temp_folder}/${sample_id}_1.fastq"
fastq2="${temp_folder}/${sample_id}_2.fastq"

# Temporary files by samtools sort are too large to keep in memory.
# Store them on the disk instead.
samtools_tmp="$(mktemp -d -p temp/nosnapshot)"

# First step (sort by QNAME) is disk intensive.
# For two parallel BAM files we often have 100% disk I/O.
# As the disk usage is the bottleneck. Only use 4 cores (per file) for this step.
# Samtools FASTQ writes directly to tmpfs memory (unless available memory is limited).
# We can use more threads for this step.
# With compression: add -c 7 to samtools fastq
samtools sort -l 9 -T "$samtools_tmp" -@ 4 -n "${bam_input}" -o - |
samtools fastq -@ 32 /dev/stdin \
-1 "${fastq1}" \
-2 "${fastq2}" \
-0 /dev/null -s /dev/null -n

rm -r "$samtools_tmp"