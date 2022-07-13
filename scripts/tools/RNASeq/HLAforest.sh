#!/bin/bash
data_type="RNA-Seq"
tool="HLAforest"
source scripts/tools/scripts/bam2fastq_tool_preamble.sh

# Tool specific code
# ------------------
# Run tool
toolscript="${toolbase}/hlaforest/scripts/CallHaplotypesPE.sh"
# For whatever reason HLAForest refuses to write to our usual workdir
# Experimentally, this seems to be caused by the presence of a dot in the filename
# Strange this is not an issue for PHLAT that also uses bowtie...
# As a workaround we specify a symbolic link in /tmp
altworkdir="/tmp/hla_forest/$(date '+%N')"
workdir_abs="$(realpath "$workdir")"

# Important to not stop in case of a non-zero exit code.
# As soon as the calling for one gene (even one of the rare genes) fails, the exit code of HLAforest is non-zero.
# Make sure to perform the last step of the script: moving the result to the results folder.
set +e
containerize_cmd <<HERE
source "${conda_profile}"
conda activate mhc2_genotyping_HLAforest

set -x

# Workaround (see above)
mkdir -p /tmp/hla_forest
ln -s "$workdir_abs" "$altworkdir"

# Important to not stop in case of a non-zero exit code.
# As soon as the calling for one gene (even one of the rare genes) fails, the exit code of HLAforest is non-zero.
# Make sure to perform the last step of the script: moving the result to the results folder.
set +e

date '+%s' > "$resmon_folder/start.txt"

time "$toolscript" "$altworkdir" "$fastq1" "$fastq2"

date '+%s' > "$resmon_folder/end.txt"
HERE

# Move results from working directory to results folder
# Results: ending in haplotypes.txt
mv "${altworkdir}/haplotypes.txt" "${res_folder}/${sample_id}_haplotypes.txt"

# Explicitly remove the symbolic link
rm "${altworkdir}"