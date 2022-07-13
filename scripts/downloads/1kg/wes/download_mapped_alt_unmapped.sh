#!/bin/bash
# Extract the MHC region and HLA contigs for each bam file
region='chr6:28509970-33480727'
hla_contigs=$(<"scripts/downloads/1kg/wes/hla_contigs.txt")
no_contig='*'

filename="$(basename "$1")"
output_file="$(sed 's/\.cram$//' <<< "${filename}").bam"

# samtools creates a bunch of temporary files in the woring directory.
# Change working directory to prevent this.
cd "downloads/1kg/bam_hs38/merged_alt_mapped"

samtools view -b -h "$1" "$region" $hla_contigs "$no_contig" > "${output_file}"
sambamba index "${output_file}"
