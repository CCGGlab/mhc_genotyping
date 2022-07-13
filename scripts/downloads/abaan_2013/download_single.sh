#!/bin/bash

url="$1"
sample_id="$(basename "$url" | grep -oP 'SRR[0-9]{7}')"
echo "${sample_id}"

# Download file (wget)
# All FASTQs per sample in separate folder
mkdir "temp/nosnapshot/download_abaan/${sample_id}"
wget -P "temp/nosnapshot/download_abaan/${sample_id}" "${url}"