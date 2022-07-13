#!/bin/bash
# Extract the entire chromosome 6 and unmapped reads
# No alternative contigs were present for Geuvadis data
region='chr6'
no_contig='*'

filename="$(basename "$1")"
output_file="$(sed 's/\.bam$//' <<< "$filename")_extracted.bam"

if [ -e "$output_file" ]; then
	echo File already exists: "$output_file" 2>&1
	exit 1
fi

# samtools creates a bunch of temporary files in the working directory
# -> move directly to output folder
cd "downloads/1kg/bam_geuvadis/"

wget "${1}.bai"
wget "$1"

samtools view -b -h "${filename}" "$region" "$no_contig" > "${output_file}"
sambamba index "${output_file}"

# Remove original files to save disk space
rm "${filename}"
rm "${filename}.bai"
