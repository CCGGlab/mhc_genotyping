#!/bin/bash
set -exo pipefail

toolpath="$(realpath "$1")"
outdir="$(realpath "$3")"
resmon_path="$(realpath "$4")"
logbase="$(realpath "$5")"
cpu_num=$(($6+10))
export cpu_num

# Useful for Kourami as another reference was used for TCGA
export DATASET='TCGA'

IFS=$'\t' read -r mapped unmapped sid <<< "$2"

tempfolder="$(mktemp -d -p temp/nosnapshot)"
cleanup() {
    if [ -n "$tempfolder" ]; then
        echo cleanup
        rm -r "$tempfolder"
    fi
}
trap cleanup EXIT

cd "$tempfolder"

# Merge mapped and unmapped
# merge the two files, store temporarily in TMP folder
samtools merge -@ 144 "${sid}.bam" "$mapped" "$unmapped"

# create index since most tools require it
samtools index -@ 144 "${sid}.bam"
cd -

toolname="$(basename "$toolpath")"
logdir="$logbase/$toolname"
mkdir -p "$logdir"
logfile="$logdir/${sid}_$(date +%s)"

bash "$toolpath" "${tempfolder}/${sid}.bam" "$outdir" "$resmon_path" > "$logfile" 2>&1
