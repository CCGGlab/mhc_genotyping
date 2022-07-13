#!/bin/bash
# Remark:
# scripts/tools/scripts/create_tool_containers.sh must be run before this script
# HLA-LA, Optitype, Polysolver and xHLA will be installed
# when the Docker containers / Conda environments are created

project_folder="$PWD"
toolbase="${project_folder}/temp/toolTest"
cd "${toolbase}"

# arcasHLA
git clone -b 'v0.2.0' --single-branch 'https://github.com/RabadanLab/arcasHLA' arcasHLA

# HLAforest
git clone 'https://github.com/FNaveed786/hlaforest' hlaforest
cd hlaforest
git reset --hard '75edb464d27d265f45f380b401111e2f32c60d56'
cd -

# HLA*LA
# Tool itself is installed when its Conda environment is installed
# Extract required resources in the graphs folder
cd "${project_folder}/temp/conda/conda_base/envs/mhc2_genotyping_hlala/opt/hla-la/graphs"
wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz
tar -xvzf PRG_MHC_GRCh38_withIMGT.tar.gz
# Create link to Conda HLA-LA install directory to the toolbase folder
ln -s "${project_folder}/temp/conda/conda_base/envs/mhc2_genotyping_hlala/opt/hla-la" "${project_folder}/temp/toolbase/hla-la"

# HLAminer
cd "$toolbase"
wget 'https://github.com/bcgsc/HLAminer/raw/master/tarballs/HLAminer_1-4.tar.gz'
tar xvf 'HLAminer_1-4.tar.gz'
rm 'HLAminer_1-4.tar.gz'

# HLA-HD
mkdir 'HLA-HD'
# Version 1.3.0 was evaluated
# Should be requested at:
# https://www.genome.med.kyoto-u.ac.jp/HLA-HD/download-request/
# Extract to temp/toolbase/HLA-HD

# HLAScan
mkdir HLAScan
cd HLAScan
wget 'https://github.com/SyntekabioTools/HLAscan/releases/download/v2.1.4/hla_scan_r_v2.1.4'
wget 'https://github.com/SyntekabioTools/HLAscan/releases/download/v2.0.0/dataset.zip'
unzip dataset.zip
cd "$toolbase"

# HLA-VBSeq
mkdir HLA_VBSeq
cd HLA_VBSeq
wget 'http://nagasakilab.csml.org/hla/HLAVBSeq.jar'
wget 'http://nagasakilab.csml.org/hla/bamNameIndex.jar'
wget 'http://nagasakilab.csml.org/hla/parse_result.pl'
wget 'http://nagasakilab.csml.org/hla/call_hla_digits.py'
# Database
wget 'http://nagasakilab.csml.org/hla/hla_all_v2.fasta'
wget 'http://nagasakilab.csml.org/hla/Allelelist_v2.txt'
# Download SamToFastq.jar (from picard-tools-1.119.zip) to temp/toolbase/HLA_VBSeq
# (If link is broken: download from https://sourceforge.net/projects/picard/)
wget 'https://altushost-swe.dl.sourceforge.net/project/picard/picard-tools/1.119/picard-tools-1.119.zip'
unzip picard-tools-1.119.zip picard-tools-1.119/SamToFastq.jar
mv 'picard-tools-1.119/SamToFastq.jar' SamToFastq.jar
rmdir picard-tools-1.119
rm picard-tools-1.119.zip
cd "$toolbase"

git clone -b 'v0.9.6' --single-branch 'https://github.com/Kingsford-Group/kourami' Kourami
# Follow compilation instructions at https://github.com/Kingsford-Group/kourami/tree/v0.9.6:
cd Kourami
scripts/download_panel.sh
mvn install

# PHLAT
# PHLAT was obtained by contacting the authors (wen.fury@regeneron.com)
# Extract folder to temp/toolbase/PHLAT

# seq2HLA
git clone 'https://github.com/TRON-Bioinformatics/seq2HLA' seq2HLA
cd seq2HLA
git reset --hard '7e0e8f5a84ed18a7e9add0bb70c3a581ca3f2455'