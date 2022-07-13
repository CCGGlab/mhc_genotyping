#!/bin/bash
token=$(</home/arne/.user/gdc.txt)
curl --header "X-Auth-Token: $token" 'https://api.gdc.cancer.gov/slicing/view/'"$1"'?region=unmapped' --output "/mnt/ccgg_data/tmp/arne/bam_slicing/TCGA_unmapped/${1}.bam"