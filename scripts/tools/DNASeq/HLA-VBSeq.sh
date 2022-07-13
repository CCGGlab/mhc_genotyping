#!/bin/bash
data_type="DNA-Seq"
tool="HLA-VBSeq"
source scripts/tools/scripts/bam_tool_preamble.sh

if [ "$monitor_on" = true ]; then
    THREADS=1
else
    THREADS=48
fi

# Tool specific code
# ------------------
set -x

# HLA-VBSeq pipeline needs to know the read length
# Derive this from the BAM file
set +e
echo "Calculate read length"
read_length=$(samtools view "${bam_input}" |
head -n 1000000 |
cut -f 10 |
perl -ne 'chomp;print length($_) . "\n"' |
sort |
uniq -c |
sed 's/^ *//' |
sort -k1n -t ' ' |
tail -n 1 |
cut -f2 -d ' ')

echo "Finished calc read length"
echo "Read length: $read_length"
set -e

# We need to change the working directory as the tools below clutter the current folder with temorary files
# Ensure all paths are absolute paths
toolbase="$(realpath "${toolbase}")"
bam_input="$(realpath "${bam_input}")"
workdir="$(realpath "${workdir}")"

# Change working directory
cd "${workdir}"

# HLA-VBSeq does not look for reads in ALT contigs and the specified coordinates
# in the documentation are for GRCh37:
# Realign to non-ALT aware version of GRCh37
"$curdir/scripts/tools/scripts/get-reads-alt-unmap-to-hg19.sh" "${bam_input}" 'output.alt-unmap.bam'

# Run HLA-VBSeq
toolfolder="${toolbase}/HLA_VBSeq"
DATABASE_DIR="${toolfolder}/hla_v2_db"
DATABASE_FILE="${DATABASE_DIR}/hla_all_v2.fasta"

JAR_DIR_BAM_NAME_INDEX="${toolfolder}/bamNameIndex"
JAR_DIR_HLA_VBSEQ="${toolfolder}/HLAVBSeq"

# Build class path
BAM_NAME_INDEX_CLASS_PATH_arr=(
    '.'
    'args4j-2.0.25.jar'
    'collections-0.1.4.jar'
    'commons-cli-1.2.jar'
    'commons-math3-3.0.jar'
    'hppc-0.5.2.jar'
    'log4j-1.2.17.jar'
    'picard-1.97.jar'
    'sam-1.97.jar'
)
BAM_NAME_INDEX_CLASS_PATH_arr2=()
for i in "${BAM_NAME_INDEX_CLASS_PATH_arr[@]}"; do
    path="$(realpath "${JAR_DIR_BAM_NAME_INDEX}/${i}")"
    BAM_NAME_INDEX_CLASS_PATH_arr2+=("$path")
done

BAM_NAME_INDEX_CLASS_PATH="$(printf '%s:' "${BAM_NAME_INDEX_CLASS_PATH_arr2[@]}" | sed 's/:$//')"

# Build HLAVBSeq Java class path
HLA_VBSEQ_CLASS_PATH_arr=(
    '.'
    'args4j-2.0.21.jar'
    'picard-1.74.jar'
    'commons-math3-3.0.jar'
    'sam-1.74.jar'
    'ssj.jar'
)

HLA_VBSEQ_CLASS_PATH_arr2=()
for i in "${HLA_VBSEQ_CLASS_PATH_arr[@]}"; do
    path="$(realpath "${JAR_DIR_HLA_VBSEQ}/${i}")"
    HLA_VBSEQ_CLASS_PATH_arr2+=("$path")
done

HLA_VBSEQ_CLASS_PATH="$(printf '%s:' "${HLA_VBSEQ_CLASS_PATH_arr2[@]}" | sed 's/:$//')"

containerize_cmd <<HERE
# Initialize conda
source "${conda_profile}"
conda activate mhc2_genotyping_vbseq_0122

# Stop on error
set -e
date '+%s' > "$resmon_folder/start.txt"

time (
# Extract a list of read names that were aligned to HLA loci
# Remark: $$1 in the instructions was replaced by $1
# The original code was suspicious and did not work!
samtools view output.alt-unmap.bam chr6:29907037-29915661 chr6:31319649-31326989 chr6:31234526-31241863 \
chr6:32914391-32922899 chr6:32900406-32910847 chr6:32969960-32979389 chr6:32778540-32786825 \
chr6:33030346-33050555 chr6:33041703-33059473 chr6:32603183-32613429 chr6:32707163-32716664 \
chr6:32625241-32636466 chr6:32721875-32733330 chr6:32405619-32414826 chr6:32544547-32559613 \
chr6:32518778-32554154 chr6:32483154-32559613 chr6:30455183-30463982 chr6:29689117-29699106 \
chr6:29792756-29800899 chr6:29793613-29978954 chr6:29855105-29979733 chr6:29892236-29899009 \
chr6:30225339-30236728 chr6:31369356-31385092 chr6:31460658-31480901 chr6:29766192-29772202 \
chr6:32810986-32823755 chr6:32779544-32808599 chr6:29756731-29767588 |
awk '{print \$1}' | sort | uniq > partial_reads.txt

# Build read name index and search read pairs and their sequences on HLA loci

java -Xmx32g -Xms32g -cp "$BAM_NAME_INDEX_CLASS_PATH" "org.eclipse.jdt.internal.jarinjarloader.JarRsrcLoader" index output.alt-unmap.bam --indexFile output.alt-unmap.bam.idx
java -cp "$BAM_NAME_INDEX_CLASS_PATH" "org.eclipse.jdt.internal.jarinjarloader.JarRsrcLoader" search output.alt-unmap.bam --name partial_reads.txt --output partial.sam
# Remark: VALIDATION_STRINGENCY=LENIENT
# Not provided in HLA-VBseq instructions, but the only way to get SamToFastq working without an error
# TODO: What causes SAMFormatException: SAM validation error: ERROR: Found 2 unpaired mates
java -jar "${toolfolder}/SamToFastq.jar" VALIDATION_STRINGENCY=LENIENT I=partial.sam F=partial_1.fastq F2=partial_2.fastq

# Extract unmapped reads
samtools view -bh -f 12 output.alt-unmap.bam > sorted_unmapped.bam
# Remark: VALIDATION_STRINGENCY=LENIENT
# Not provided in HLA-VBseq instructions, but the only way to get SamToFastq working without an error
java -jar "${toolfolder}/SamToFastq.jar" VALIDATION_STRINGENCY=LENIENT I=sorted_unmapped.bam F=unmapped_1.fastq F2=unmapped_2.fastq

# Combine reads in FASTQ format
cat partial_1.fastq unmapped_1.fastq > part_1.fastq
cat partial_2.fastq unmapped_2.fastq > part_2.fastq

# Alignment by BWA-MEM allowing multiple alignments for each read
# This is has to be done once: bwa index $DATABASE_DIR/$DATABASE_FILE
# -t    Number of threads
# -P    In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair.
# -L    Clipping penalty. When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied. Note that in
#       this case, the SAM AS tag reports the best SW score; clipping penalty is not deduced. If two numbers are provided, the first is for 5'-end clipping and second for 3'-end clipping. [5]
# -a    Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.
bwa mem -t $THREADS -P -L 10000 -a "$DATABASE_FILE" part_1.fastq part_2.fastq > part.sam

# Here, alpha_zero is a hyperparameter as described in the paper and we recommend to use 0.01.
java -cp "$HLA_VBSEQ_CLASS_PATH" org.eclipse.jdt.internal.jarinjarloader.JarRsrcLoader "$DATABASE_FILE" "part.sam" "result.txt" --alpha_zero 0.01 --is_paired

echo "--- Analyzing output ---"
python "${toolfolder}/call_hla_digits.py" -v result.txt -a "$DATABASE_DIR/Allelelist_v2.txt" -r $read_length -d 4 --ispaired > report.d4.txt
)

date '+%s' > "$resmon_folder/end.txt"
HERE

cd -

# Move results from working directory to results folder
for fl in "${workdir}"/*.txt
do
    mv "${fl}" "${res_folder}/${sample_id}_$(basename "$fl")"
done