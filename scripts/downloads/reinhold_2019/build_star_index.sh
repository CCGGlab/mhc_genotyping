#!/bin/bash
# Download GRCh38.d1.vd1.fa from:
# https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files
# -> https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834

# Download gencode.v38.chr_patch_hapl_scaff.annotation.gtf
# -> https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz

# Use default parameters
STAR --runThreadN 64 \
--runMode genomeGenerate \
--genomeDir temp/STAR/index_rna \
--genomeFastaFiles downloads/genomes/GDC/GRCh38.d1.vd1.fa \
--sjdbGTFfile downloads/gencode/gencode.v38.chr_patch_hapl_scaff.annotation.gtf
