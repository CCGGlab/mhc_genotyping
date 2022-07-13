#!/bin/bash
data_type="RNA-Seq"
tool="seq2HLA"
source scripts/tools/scripts/bam2fastq_tool_preamble.sh

# Tool specific code
# ------------------
# Run tool
toolscript="${toolbase}/seq2HLA/seq2HLA.py"

containerize_cmd <<HERE
source "${conda_profile}"
conda activate mhc2_genotyping_seq2hla_3001

date '+%s' > "$resmon_folder/start.txt"
# -r: sample_id is the prefix of the filename. Seq2HLA will append underscore and a suffix.
time python2 "$toolscript" -1 "${fastq1}" -2 "${fastq2}" -r "${workdir}/${sample_id}"
date '+%s' > "$resmon_folder/end.txt"
HERE

# Move results from working directory to results folder
# Results: keep all files
mv "${workdir}"/* "${res_folder}/"