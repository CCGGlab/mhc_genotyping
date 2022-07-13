#!/bin/bash
tool="HLAminer"
source scripts/tools/scripts/bam2fastq_tool_preamble.sh

# Tool specific code
# ------------------
# Run tool
## We are obliged to change the working directory for HLAminer (to prevent clotting the current folder)
## -> Replace all defined paths by absolute paths
temp_folder="$(realpath "$temp_folder")"
workdir="$(realpath "$workdir")"
toolbase="$(realpath "$toolbase")"
fastq1="$(realpath "$fastq1")"
fastq2="$(realpath "$fastq2")"

current_toolbase="${toolbase}/HLAminer-1.4/HLAminer_v1.4"

# Tool subfolders
toolscript="${current_toolbase}/bin/HLAminer.pl"
database="${current_toolbase}/database"

# Database files
if (grep -iq 'dna' <<< "$data_type"); then
    reference_seq="${database}/HLA-I_II_GEN.fasta"
elif (grep -iq 'rna' <<< "$data_type"); then
    reference_seq="${database}/HLA-I_II_CDS.fasta"
fi

hla_nom_p="${database}/hla_nom_p.txt"

## Align FASTQ using BWA
sai1_output="${temp_folder}/aln_1.sai"
sai2_output="${temp_folder}/aln_2.sai"
sam_output="${temp_folder}/aln.sam"

containerize_cmd <<HERE
source "${conda_profile}"
conda activate mhc2_genotyping_hlaminer

# Alignment to HLA sequence database is considered part of the official pipeline
date '+%s' > "$resmon_folder/start.txt"

time (
## Create bwa index
if [ ! -f "${fastq1}.amb" ] || [ ! -f "${fastq2}.amb" ]; then
    echo "----- Running bwa index since index could not be found ----- " >&2
    bwa index "$fastq1"
    bwa index "$fastq2"
fi

# As suggested in bin/HPRAwgs_classI-II.sh and bin/HPRArnaseq_classI-II.sh
bwa aln -e 0 -o 0 "$reference_seq" "$fastq1" > "$sai1_output"
bwa aln -e 0 -o 0 "$reference_seq" "$fastq2" > "$sai2_output"
bwa sampe -o 1000 "$reference_seq" "$sai1_output" "$sai2_output" "$fastq1" "$fastq2" > "$sam_output"

## Predict HLA
# Set output folder of tool to workdir
cd "$workdir"

"$toolscript" -a "${sam_output}" -h "$reference_seq" -p "$hla_nom_p" -s 500 -l "${sample_id}"
)

# Go back to the previous folder
cd -

date '+%s' > "$resmon_folder/end.txt"
HERE

# Move results from working directory to results folder
# Results: keep all files
mv "${workdir}"/*.csv "${res_folder}/"
mv "${workdir}"/*.log "${res_folder}/"