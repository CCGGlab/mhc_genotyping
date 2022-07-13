#!/bin/bash
# Project accession number: PRJNA474665
# FASTQ file list obtained from
# https://www.ebi.ac.uk/ena/browser/view/PRJNA474665?show=reads
# -> file_list.txt

# Download files
mkdir temp/nosnapshot/download_abaan

xargs -L1 bash scripts/downloads/abaan_2013/download_single.sh < \
scripts/downloads/abaan_2013/file_list.txt

# Next step:
# scripts/downloads/abaan_2013/run_preproc_1kg_nci60.sh

# Cleanup temp/nosnapshot/download_abaan when successful