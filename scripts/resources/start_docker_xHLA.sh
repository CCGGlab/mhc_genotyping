#!/bin/bash
set -x

if [ "$#" -ne 3 ]; then
    echo Syntax: "$0" 'docker_cid script_file monitor_on' >&2
    exit 1
fi

docker_cid="$1"
script_file="$2"
monitor_on="$3"

# Note:
# --cpuset-cpus: Limit the specific CPUs or cores a container can use
# - This option was only defined when resource consumption was monitored.
cpuset_param=""
if [ $monitor_on = true ]; then
    cpuset_param="--cpuset-cpus=$cpu_num"
fi

docker run \
$cpuset_param \
--cidfile "$docker_cid" \
--rm \
--entrypoint='' \
-e HOME="$HOME" \
-v /home:/home \
-v /tmp:/tmp \
-v /home/labgroups/ccgg:/home/labgroups/ccgg \
-w "$PWD" \
humanlongevity/hla:0.0.0 \
/bin/bash "$script_file"

# Create CID file in case of an error
# Allows the main monitor script to continue
echo CLOSED > "$docker_cid"
