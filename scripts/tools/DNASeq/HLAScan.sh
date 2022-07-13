#!/bin/bash
data_type="DNA-Seq"
tool="HLAScan"
# HLAScan can run in either FASTQ or BAM mode
# HLAScan in FASTQ mode appears to be extremely slow.
# -> BAM mode was used
# HLAScan supports at least chr6 and chr6 ALT
source scripts/tools/scripts/bam_tool_preamble.sh

# Tool specific code
# ------------------
toolfolder="${toolbase}/HLAScan"
toolscript="${toolfolder}/hla_scan_r_v2.1.4"

DATABASE="${toolfolder}/db/HLA-ALL.IMGT"
if [ "$monitor_on" = true ]; then
    THREADS=1
else
    THREADS=96
fi
REF_VERSION=38

# For some reason HLAScan always returns exit code 1...
# (Issue is known on GitHub: https://github.com/SyntekabioTools/HLAscan/issues/4)
# Disable "stop on error" in Bash
containerize_cmd <<HERE
# Initialize conda
source "${conda_profile}"
conda activate mhc2_genotyping_default

date '+%s' > "$resmon_folder/start.txt"

time (
set +e
for locus in A B C DPA1 DPB1 DQA1 DQB1 DRB1; do
    time "${toolscript}" -b "${bam_input}" -d $DATABASE -v "${REF_VERSION}" -t $THREADS -g "HLA-\${locus}" > "${workdir}/${sample_id}_\${locus}_result.txt"
    # fi
done
)

date '+%s' > "$resmon_folder/end.txt"
HERE

# Move result file
mv "${workdir}"/*_result.txt "${res_folder}/"