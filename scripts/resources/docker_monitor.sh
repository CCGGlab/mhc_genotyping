#!/bin/bash
# Must be hard coded: we have no guarantee this script is ran from the project directory
projectdir='/home/arne/projects/mhc2_genotyping'

if [ "$#" -ne 4 ]; then
    echo Syntax: "$0" 'tool script_file resource_log_folder monitor_on' >&2
    exit 1
fi

tool="$1"
script_file="$2"
resource_log_folder="$3"
monitor_on="$4"

# Need to create the CID file in an unsafe manner...
# Docker does not allow this file to exist already
docker_cid="$(mktemp -u)"
docker_start_script="${projectdir}/scripts/resources/start_docker_${tool}.sh"
if ! [ -e "$docker_start_script" ]; then
    docker_start_script="${projectdir}/scripts/resources/start_docker_conda.sh"
fi

set -x
"$docker_start_script" \
"${docker_cid}" "$script_file" "$monitor_on" &
pid=$!

if [ "$monitor_on" = true ]; then

# Wait until CID file exists and is non empty before we can start monitoring
docker_id=''
while [ -z "$docker_id" ]; do
    if [ -e "$docker_cid" ]; then
        docker_id="$(<${docker_cid})"
    fi
    sleep 0.5
done

set +x
# This outputs the stats each 500 ms (see Docker source code)
# (without streaming this interval is not possible)
stdbuf -i0 -o0 -e0 docker stats \
--no-trunc \
--format '{{json .}}' "$docker_id" > "${resource_log_folder}/resources.txt" &
pidmon=$!
fi

# Wait until Docker exits, then stop monitoring
wait $pid
if [ "$monitor_on" = true ]; then
    kill $pidmon
fi

rm ${docker_cid}