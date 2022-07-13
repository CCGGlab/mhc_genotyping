#!/bin/bash
cleanup() { 
    echo "Cleanup at exit" 2>&1
    if [ -n "$temp_folder" ]; then
        # Note: we must use the absolute path here.
        # Necessary in case crash occurs in a script where we need to change the working dir.
        source '/home/arne/projects/mhc2_genotyping/scripts/tools/scripts/bam2fastq_tool_footer.sh'
    fi
}

trap cleanup EXIT
set -x

if [ $# -ne 3 ]; then
    echo Syntax: "$0" "bam_input" "res_folder" "resmon_base" >&2
    exit 1
fi

bam_input="$1"
res_base="$2"
resmon_base="$3"
sample_id="$(basename "${bam_input}" | sed 's/\.bam$//')"

if [ -z "$data_type" ]; then
    echo "data_type must be defined" >&2
    exit 1
fi

if [ -z "$tool" ]; then
    echo "tool must be defined" >&2
    exit 1
fi

# Create results folder
res_folder="${res_base}/${data_type}/${tool}"
mkdir -p "${res_folder}"

if ls "${res_folder}" | grep -q "${sample_id}"; then
    echo "Skipping ${sample_id} for ${tool} as the result file already exists." >&2
    exit 0
fi

# Load default mhc2_genotyping conda environment
conda_base="$(realpath "temp/conda/conda_base")"
# Used in containers: sets CONDARC and then executes the profile script
conda_profile="$(realpath "temp/conda/conda_profile.sh")"
source "$conda_profile"

conda_default_env="${conda_base}/envs/mhc2_genotyping_default"
conda activate mhc2_genotyping_default

# Debugging options
set -exo pipefail
# -e: exit immediately if any command fails
# (-u: reference to undeclared variables is an error)
# --> not chosen, gives problem with Conda activation
# -x: be verbose
# -o pipefail: one command in pipeline fails -> entire command fails

# Create temporary folder
temp_folder_base="$(mktemp -d -p temp/nosnapshot)"
temp_folder="${temp_folder_base}/${data_type}/${tool}/${sample_id}"
mkdir -p "${temp_folder}"

# Create working directory (output created by tool itself)
workdir="${temp_folder}/workdir"
mkdir -p "${workdir}"

# The scripts for the tools are currently stored in the downloads folder...
# Define this location as a variable so we can easily change this afterwards.
toolbase="temp/toolbase"

# Save current folder: tool specific scripts might change the working directory
curdir="$PWD"

# Resource monitoring
# ===================
# Disable resource monitoring when dummy folder is used
monitor_on=true
if (egrep -q '_dummy/?$' <<< "${resmon_base}"); then
    monitor_on=false
fi
export "$monitor_on"

# Create folder for resource monitoring
# Convert to absolute path (realpath) as we run this inside a Docker environment
resmon_folder="$(realpath -m "${resmon_base}/${data_type}/${tool}/${sample_id}")"
mkdir -p "${resmon_folder}"

# Folder for containerized scripts
if ! [ -d /tmp/benchmark_scripts ]; then
    mkdir /tmp/benchmark_scripts
fi
resmon_script_file="$(mktemp -p /tmp/benchmark_scripts/)"

# This approach with piping the script from a here document
# and reading with cat allows to execute an action for all scripts
# before or after the actual script
containerize_cmd() {
cat > "$resmon_script_file"
echo 'chown -R 1000:1002 "'"$workdir"'"' >> "$resmon_script_file"
"$curdir/scripts/resources/docker_monitor.sh" "$tool" "$resmon_script_file" "$resmon_folder" "$monitor_on"
}