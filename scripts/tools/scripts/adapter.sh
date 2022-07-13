#!/bin/bash
toolpath="$1"
file_path="$2"
outdir="$3"
resmon_path="$4"
logbase="$5"
cpu_num=$(($6+10))
export cpu_num

# Useful for Kourami as another reference was used for TCGA (no ALT contigs)
# For NCI60 the same reference as 1000 genomes was used
# as the samples were realigned to that reference directly after download
export DATASET='DEFAULT'

toolname="$(basename "$toolpath")"
file_basename="$(basename "$file_path")"

logdir="$logbase/$toolname"
mkdir -p "$logdir"
logfile="$logdir/${file_basename}_$(date +%s)"

bash "$toolpath" "$file_path" "$outdir" "$resmon_path" > "$logfile" 2>&1