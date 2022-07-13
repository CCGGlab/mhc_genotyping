#!/bin/bash
tool="HLA-HD"
source scripts/tools/scripts/bam2fastq_tool_preamble.sh

if [ "$monitor_on" = true ]; then
    THREADS=1
else
    THREADS=96
fi

# Tool specific code
# ------------------
# Initialize conda
# Nothing needs to be done. Script uses default mhc2_genotyping conda environment
toolfolder="${toolbase}/HLA-HD/hlahd.1.3.0"
toolbin="${toolfolder}/bin"

containerize_cmd <<HERE
source "${conda_profile}"
conda activate mhc2_genotyping_default

# Test if scripts in HLA-HD's bin folder are included in the PATH. If not, do so.
if [[ $PATH != *"${toolbin}"* ]]; then
    export PATH="${toolbin}":"$PATH"
fi

date '+%s' > "$resmon_folder/start.txt"
time hlahd.sh -t "$THREADS" -m 30 -c 0.95 -f "${toolfolder}/freq_data/" "${fastq1}" "${fastq2}" "${toolfolder}/HLA_gene.split.txt" "${toolfolder}/dictionary/" "${sample_id}" "${workdir}/"
date '+%s' > "$resmon_folder/end.txt"
HERE

# Move results from working directory to results folder
# Results: 
mv "${workdir}"/*/result/*.txt "${res_folder}/"