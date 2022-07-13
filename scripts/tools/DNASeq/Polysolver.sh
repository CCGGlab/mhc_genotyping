#!/bin/bash
data_type="DNA-Seq"
tool="Polysolver"
source scripts/tools/scripts/bam_tool_preamble.sh

# Tool specific code
# ------------------

# We need to change the working directory as the tools below clutter the current folder with temorary files
# Ensure all paths are absolute paths
toolbase="$(realpath "${toolbase}")"
bam_input="$(realpath "${bam_input}")"
workdir="$(realpath "${workdir}")"

# Change working directory
cd "${workdir}"

# The input BAM file for all tools is aligned according to 1000 genomes specifications
# -> contains ALT contigs
# This region is not extracted by Polysolver and likely ignored.
# => Realign to non-ALT aware version
"$curdir/scripts/tools/scripts/get-reads-alt-unmap.sh" "${bam_input}" 'output.alt-unmap.bam'

# Polysolver expects contig names without "chr"
"$curdir/scripts/tools/scripts/reheader_chr.sh" 'output.alt-unmap.bam' 'output.alt-unmap-nochr.bam'
rm 'output.alt-unmap.bam'

# Run Polysolver
containerize_cmd <<HERE
source /home/polysolver/scripts/config.bash

date '+%s' > "$resmon_folder/start.txt"
time /bin/bash /home/polysolver/scripts/shell_call_hla_type 'output.alt-unmap-nochr.bam' Unknown 0 hg38 STDFQ 0 "${workdir}"
date '+%s' > "$resmon_folder/end.txt"
HERE
cd -

# Move results from working directory to results folder
for fl in "${workdir}"/*.txt
do
    mv "${fl}" "${res_folder}/${sample_id}_$(basename "$fl")"
done