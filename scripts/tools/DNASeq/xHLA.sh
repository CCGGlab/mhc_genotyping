#!/bin/bash
data_type="DNA-Seq"
tool="xHLA"
source scripts/tools/scripts/bam_tool_preamble.sh

# Tool specific code
# ------------------

# We need to change the working directory as the xHLA clutters the current folder with temorary files
# Ensure all paths are absolute paths
toolbase="$(realpath "${toolbase}")"
bam_input="$(realpath "${bam_input}")"
workdir="$(realpath "${workdir}")"

toolfolder="${toolbase}/xHLA/HLA"

# Change working directory
cd "${workdir}"
# The input BAM file for all tools is aligned according to 1000 genomes specifications
# -> contains ALT contigs
# xHLA does not support this. We apply the recommended preprocessing script.
# See: https://github.com/humanlongevity/HLA/wiki/BAMs-compatible-with-xHLA
"$curdir/scripts/tools/scripts/get-reads-alt-unmap.sh" "${bam_input}" 'output.alt-unmap.bam'

containerize_cmd <<HERE
date '+%s' > "$resmon_folder/start.txt"
time /opt/bin/run.py --sample_id "${sample_id}" --input_bam_path "${workdir}/output.alt-unmap.bam" --output_path "${workdir}"
date '+%s' > "$resmon_folder/end.txt"
HERE

cd -

# Move results from working directory to results folder
# Results: 
mv "${workdir}"/*.json "${res_folder}/"