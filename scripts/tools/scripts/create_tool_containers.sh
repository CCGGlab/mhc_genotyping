#!/bin/bash
# Please install Miniconda3 to temp/conda/conda_base before running this script
# Answer "no" to "Do you wish the installer to initialize Miniconda3"

# Initialize Conda environments
conda_base="$(realpath "temp/conda/conda_base")"
conda_profile="${conda_base}/etc/profile.d/conda.sh"
source "${conda_profile}"

yml_dir='scripts/tools/scripts/conda_environments'

for i in "$(find "$yml_dir" -type f)"; do
  conda env create -f "$i"
done

# Download and install Docker images
# Ubuntu 20.04 container is used for all tools that use Anaconda
docker pull ubuntu:20.04
docker pull humanlongevity/hla:0.0.0
docker pull sachet/polysolver:v4