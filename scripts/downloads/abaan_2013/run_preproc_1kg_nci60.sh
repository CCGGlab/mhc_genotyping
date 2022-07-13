#!/bin/bash
set -euxo pipefail

# First download FASTQ files to temp/nosnapshot/download_abaan
# using download.sh

align_fastq() {
    sample_id="$1"
    output="downloads/pub/abaan_2013/HLA/${sample_id}.bam"
    # Skip if the required output file has been created already.
    if [ -e "$output" ]; then
        return
    fi
    bash 'scripts/downloads/abaan_2013/align_fastq_dna_1kg.sh' \
        "temp/nosnapshot/download_abaan/${sample_id}" \
        "${sample_id}" \
        "fastq.gz" \
        "$output"
}
export -f align_fastq

echo 5 > /tmp/align_jobs
ls temp/nosnapshot/download_abaan | grep -o 'SRR[0-9]*' | sort -u |
parallel --memfree 50G \
-j /tmp/align_jobs \
--env 'align_fastq' \
--results 'temp/nosnapshot/logalign/{/.}_{= $_=yyyymmddhhmmss() =}' \
align_fastq '{}'
