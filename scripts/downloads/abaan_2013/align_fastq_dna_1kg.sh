#!/bin/bash
set -x
cleanup() { 
    echo "Cleanup at exit" 2>&1
    if [ -e "$tempdir" ]; then
        rm -r "$tempdir"
    fi
}

trap cleanup EXIT
tempdir="$(mktemp -d -p temp/nosnapshot/)"
# tempdir="$(mktemp -d -p /mnt/zram0/fastq/)"

input_folder="$1"
sample_id="$2"
file_ext="$3"
OUT="$4"

F0="${input_folder}/${sample_id}.${file_ext}"
F1="${input_folder}/${sample_id}_1.${file_ext}"
F2="${input_folder}/${sample_id}_2.${file_ext}"

NCORE=$(grep -c '^processor' /proc/cpuinfo)
MEM=32GB
echo "NCORE=$NCORE MEM=$MEM"

# Reference used by 1000 genomes project
REFERENCE='downloads/genomes/hs38DH_1kg/GRCh38_full_analysis_set_plus_decoy_hla.fa'

# "As of now, ALT mapping is done in two separate steps: BWA-MEM mapping and postprocessing."
# prefix_hla_hit: prefix of output files containing sequences matching HLA genes [null]
# This creates separate FASTQ files e.g. $tempdir/${filename}_hla.HLA-A.fq
# Ignored as these files are not used in our pipeline.
alt_aligned_bam_file="${tempdir}/aligned.bam"

align_fastq () {
    if [ -n "$2" ]; then
        files=("$1" "$2")
    else
        files=("$1")
    fi
    output="$3"

    # The M parameter (markShorterSplits) is the only non-default parameter
    # → It was used in the 1000 genomes alignment pipeline
    # -1 parameter for samtools view (compression level 1) was not used.
    # We want to used default compression settings.
    bwa mem -t $NCORE -M -B 4 -O 6 -E 1 "$REFERENCE" "${files[@]}" > "${output}.aln.bam"
    echo "Alignment finished"

    cat "${output}.aln.bam" |
    k8 /home/labgroups/ccgg/tools/bwa/bwakit/bwa-postalt.js "${REFERENCE}.alt" |
    samtools view -bh -o "${output}.post.bam" -
    echo "Post-alignment finished"
    
    sambamba sort -m $MEM -t $NCORE -o "${output}" "${output}.post.bam"
    echo "Sort finished"
}

paired=''
single=''
if ! [ -e "$alt_aligned_bam_file" ]; then
    # Single end reads
    if [ -e "$F0" ]; then
        single="${tempdir}/single.bam"
        align_fastq "$F0" '' "$single"
    fi
    # Paired end reads
    if [ -e "$F1" ] && [ -e "$F2" ]; then
        paired="${tempdir}/paired.bam"
        align_fastq "$F1" "$F2" "$paired"
    fi
    # Combine single and paired
    if [ -n "$paired" ] && [ -n "$single" ]; then
        samtools merge "$alt_aligned_bam_file" "$single" "$paired"
    # Only files for paired-end reads exist [most common]
    elif [ -n "$paired" ]; then
        mv "$paired" "$alt_aligned_bam_file"
    # Only files for single-end reads exist
    elif [ -n "$single" ]; then
        mv "$single" "$alt_aligned_bam_file"
    else # None of these files exists. The input parameters are incorrect.
        echo "${input_folder}" "${sample_id}" "${file_ext}" 'is invalid' >&2
        exit 1
    fi
fi

# -> Step 2 and 3 from 1000 genomes alignment pipeline were skipped.
# The indel realignment tools are not longer a part of GATK v4.
# (excluded from recommended variant calling pipeline)
# Unclear whether this is still recommended for HLA genotyping.

# 4. Mark Duplicates using BioBamBam
# Go to tempdir to avoid bammarkduplicates creates temporary files in the project folder
alt_aligned_bam_file="$(realpath "$alt_aligned_bam_file")"
cd "${tempdir}"
bammarkduplicates I="$alt_aligned_bam_file" O="mdup.bam" index=1 rmdup=0
cd -

# Cleanup intermediate files
rm "${alt_aligned_bam_file}"

primary_contig='chr6:28509970-33480727'
alt_contigs=$(<"scripts/downloads/abaan_2013/alt_chr6_contigs.txt")
hla_contigs=$(<"scripts/downloads/abaan_2013/hla_contigs.txt")
no_contig='*'

samtools view -o "${OUT}" -b "${tempdir}/mdup.bam" $primary_contig $alt_contigs $hla_contigs "$no_contig"
sambamba index "${OUT}"
