#!/bin/bash
tool="PHLAT"
source scripts/tools/scripts/bam2fastq_tool_preamble.sh

# Tool specific code
# ------------------
# Run tool
phlatdir="${toolbase}/PHLAT"
indexdir="$phlatdir/b2folder"
b2url="bowtie2"
toolscript="${phlatdir}/dist/PHLAT.py"

containerize_cmd <<HERE
# Initialize conda
source "${conda_profile}"
conda activate mhc2_genotyping_phlat

date '+%s' > "$resmon_folder/start.txt"

time python -O "$toolscript" -1 "$fastq1" -2 "$fastq2" -index $indexdir -b2url $b2url -orientation "--fr" -tag "$sample_id" -e $phlatdir -o "$workdir/" -tmp 0

date '+%s' > "$resmon_folder/end.txt"
HERE

# Move results from working directory to results folder
# Results: keep all files
mv "${workdir}"/*.sum "${res_folder}/"