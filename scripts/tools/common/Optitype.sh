#!/bin/bash
tool="Optitype"
source scripts/tools/scripts/bam2fastq_tool_preamble.sh

# Tool specific code
# ------------------
toolconda="${conda_base}/envs/mhc2_genotyping_optitype"

toolfolder="${toolconda}/bin"
toolscript="${toolfolder}/OptiTypePipeline.py"

tooldata="${toolfolder}/data"
reference="${tooldata}/hla_reference_rna.fasta"

# Run Optitype
containerize_cmd <<HERE
# Initialize conda
source "${conda_profile}"
conda activate mhc2_genotyping_optitype

date '+%s' > "$resmon_folder/start.txt"

time (
# Add "enumerate" parameter to make soft predictions (output a few suboptimal solutions as well)
if (grep -iq 'dna' <<< "$data_type"); then
    python "${toolscript}" \
    -i "${fastq1}" "${fastq2}" \
    --dna \
    -v \
    --outdir "${workdir}"
elif (grep -iq 'rna' <<< "$data_type"); then
    python "${toolscript}" \
    -i "${fastq1}" "${fastq2}" \
    --rna \
    -v \
    --outdir "${workdir}"
fi
)

date '+%s' > "$resmon_folder/end.txt"
HERE

# Move results (*.txt) from working directory to results folder
# The coverage plot (*.pdf) might be interesting to keep as well
while read fl
do
    name="$(basename "$fl")"
    cp "$fl" "${res_folder}/${sample_id}_${name}"
done < <(find "${workdir}" -name '*.pdf' -o -name '*.tsv')