#!/bin/bash
# Used by Polysolver pipeline: remove the "chr" from the contig names
samtools view -H "$1" | sed '/^@SQ/s/SN\:chr/SN\:/' | samtools reheader - "$1" > "$2"
samtools index "$2"