#!/bin/bash
source scripts/tools/scripts/general_tool_preamble.sh

# TODO: When many tools need FASTQ files, this is redundant.
# Find a way to generate the FASTQ once and share these resources.
source scripts/tools/scripts/bam_to_fastq.sh

if [ ! -e "$fastq1" ] || [ ! -e "$fastq2" ]; then
    echo "FASTQ not created successfully. Exiting for the current sample"
    exit 1
fi