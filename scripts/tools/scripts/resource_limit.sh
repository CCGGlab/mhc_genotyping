#!/bin/bash
# Exit code:
# 0   Below limit. Start another job.
# 1   Over limit. Start no jobs.
# 2   Way over limit. Kill the youngest job.

# Based on functions defined in GNU parallel itself

# Adapted so io > 99 => kill job

cleanup() {
    if [ -n "$tempfile" ] && [ -e "$tempfile" ]; then
        rm "$tempfile"
    fi
}
trap cleanup EXIT

tempfile=$(mktemp)

io() {
    limit=$1
    io_file=$2
    # Do the measurement in the background
    (
        (
            LANG=C iostat -x 1 2 > $io_file;
        ) </dev/null >/dev/null
    );
    # if ($max >= 99) { exit(2); }
    perl -e '
    for(reverse <>) {
        /Device/ and last;
        /(\S+)$/ and $max = $max > $1 ? $max : $1;
    };
    exit ('$limit' < $max);
    ' $io_file;
};

mem() {
    # Original script overestimated the available memory
    # = total amount of memory in SwapCached, Cached, MemFree, Buffers
    # -> replace by MemAvailable
    limit=$1;
    # ^((Swap)?Cached|MemFree|Buffers):/{ sum += $2}

    # Largest part of ARC cached is available as well
    # Subtract 13 GB (minimum ARC size that cannot be reclaimed)
    arc_size=$(arc_summary -r | grep '^ *size ' | grep -o '[0-9]*$')
    if [ ${arc_size} -gt $((13*1024**3)) ]; then
        arc_size=$((arc_size - 13*1024**3))
    else
        arc_size=0
    fi
    
    awk '/^MemAvailable:/{sum = $2}
        END {
            available = sum * 1024 + '$arc_size';
            if (sum*1024 < '$limit'/2) { exit 2; }
            else { exit (sum*1024 < '$limit') }
        }' /proc/meminfo;
};

# Combine mem 20G and io 70 criteria
# Parameter 1 of "mem": size in bytes
if mem $((70 * 1024**3)); then
    io 70 $tempfile
    errcode=$?
    if [ $errcode -eq 1 ]; then
        date >> /tmp/ionostart
    fi
    if [ $errcode -eq 2 ]; then
        date >> /tmp/ioterm
    fi
    exit $errcode
else
    errcode=$?
    if [ $errcode -eq 2 ]; then
        date >> /tmp/memterm
    fi
    exit $errcode
fi
