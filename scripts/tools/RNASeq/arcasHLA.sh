#!/bin/bash
data_type="RNA-Seq"
tool="arcasHLA"
source scripts/tools/scripts/bam2fastq_tool_preamble.sh
if [ "$monitor_on" = true ]; then
    THREADS=1
else
    THREADS=96
fi

# Tool specific code
# ------------------
# Initialize conda
# Uses default Conda

toolscript="${toolbase}/arcasHLA/arcasHLA"
containerize_cmd <<HERE
source "${conda_profile}"
conda activate mhc2_genotyping_default

date '+%s' > "$resmon_folder/start.txt"
time "$toolscript" genotype "${fastq1}" "${fastq2}" -g A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1 -o "${workdir}/" -t $THREADS -v
date '+%s' > "$resmon_folder/end.txt"
HERE

# Move results from working directory to results folder
# Results: *.json (.genes.json) and (.genotype.json)
# We need to specify the output filename manually.
# ArcasHLA automatically truncates the filename which can lead to overwritten files.
mv "${workdir}"/*.genes.json "${res_folder}/${sample_id}.genes.json"
mv "${workdir}"/*.genotype.json "${res_folder}/${sample_id}.genotype.json"