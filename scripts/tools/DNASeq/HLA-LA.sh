#!/bin/bash
# Requirement: BAM file allowed to one of the reference genomes supported by HLA-LA
# (auto-detects used reference genome)
data_type="DNA-Seq"
tool="HLA-LA"
# HLA-LA works directly from the BAM file
source scripts/tools/scripts/bam_tool_preamble.sh

if [ "$monitor_on" = true ]; then
    THREADS=1
else
    THREADS=32
fi

# Tool specific code
# ------------------
toolfolder="${toolbase}/hla-la"
toolscript="${toolfolder}/src/HLA-LA.pl"

containerize_cmd <<HERE
# Initialize conda
source "${conda_profile}"
conda activate mhc2_genotyping_hlala

date '+%s' > "$resmon_folder/start.txt"

# Sample ID must consist only of alphanumeric characters
# -> Sample ID is set to "sample": below the files are renamed anyway
time "${toolscript}" --BAM "${bam_input}" --graph PRG_MHC_GRCh38_withIMGT --sampleID "sample" --maxThreads "$THREADS" --workingDir "${workdir}"

date '+%s' > "$resmon_folder/end.txt"
HERE

# Move results from working directory to results folder
# Results:
# For each sample with ID $mySampleID, the main output file is ../working/$mySampleID/hla/R1_bestguess_G.txt. Columns:
# (see documentation HLA-LA)
# Do not fail if one command fails
set +e
for i in "${workdir}"/*/hla/*.txt; do
    destfile="${sample_id}_$(basename "$i")"
    mv "$i" "${res_folder}/${destfile}"
done