#!/bin/bash

if [ -z "$temp_folder_base" ]; then
    echo "temp_folder_base was not defined" >&2
    exit 1    
fi

rm -r "$temp_folder_base"

# Remove tool specific script that is containerized
if [ -n "$resmon_script_file" ] && [ -e "$resmon_script_file" ]; then
    rm "$resmon_script_file"
fi