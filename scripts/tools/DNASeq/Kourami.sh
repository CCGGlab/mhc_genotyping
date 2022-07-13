#!/bin/bash
data_type="DNA-Seq"
tool="Kourami"
source scripts/tools/scripts/bam_tool_preamble.sh

# Tool specific code
# ------------------
# We need to change the working directory as the Kourami scripts write to the current folder
# Ensure all paths are absolute paths
toolbase="$(realpath "${toolbase}")"
bam_input="$(realpath "${bam_input}")"
workdir="$(realpath "${workdir}")"

toolfolder="${toolbase}/Kourami"
KOURAMI_JAR="${toolfolder}/target/Kourami.jar"
DATABASE="${toolfolder}/custom_db/3.43.0"

# Note: alignAndExtract_hs38DH.sh was used for 1000 genomes and NCI-60 datasets.
# TCGA does not have HLA contigs. Here the alignAndExtract_hs38 script was used.
if [ "$DATASET" = 'TCGA' ] ; then
    kourami_ref_script='alignAndExtract_hs38.sh'
else
    kourami_ref_script='alignAndExtract_hs38DH.sh'
fi

cd "$workdir"

containerize_cmd <<HERE
# Initialize conda
source "${conda_profile}"
conda activate mhc2_genotyping_kourami

date '+%s' > "$resmon_folder/start.txt"

time (
# Generate alignment to Kourami panel for each sample
"${toolfolder}/scripts/${kourami_ref_script}" -d $DATABASE "${sample_id}" "${bam_input}" 
# Use the new alignment to call Kourami
java -jar $KOURAMI_JAR -a -d $DATABASE -o "${workdir}/${sample_id}" "${sample_id}_on_KouramiPanel.bam"
)

date '+%s' > "$resmon_folder/end.txt"
HERE

# -o parameter: outputs different files in $workdir with sample_id as the base of the filename
# with a suffix determined by Kourami and as extensions: .result and .log

cd -

# Move results from working directory to results folder
# Results: 
mv "${workdir}"/*.log "${workdir}"/*.result "${res_folder}/"