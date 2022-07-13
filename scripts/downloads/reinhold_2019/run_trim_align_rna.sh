#!/bin/bash
set -exo pipefail
# Path to project Conda environment created using mhc2_genotyping.yml
condapath="$PWD/mhc2_genotyping/"

sample_id="$1"
output_folder="$PWD/downloads/pub/reinhold_2019/HLA"

ls downloads/pub/reinhold_2019/HLA | grep -q "${sample_id}" &&
{ echo File exists; exit 0; }

tempfolder="$(mktemp -d -p "$PWD/temp/nosnapshot/")"
cleanup() { 
    if [ -n "$tempfolder" ]; then
	echo cleanup
	rm -r "$tempfolder"
    fi
}
trap cleanup EXIT

# align to GRCh38.d1.vd1.fa using STAR
REFERENCE="$(realpath temp/STAR/index_rna)"

# Download samples
find "$PWD/downloads/pub/reinhold_2019/full/${sample_id}" -name '*.fastq.gz' -exec ln -s '{}' "$tempfolder/" ';'

cd "$tempfolder"

# Trim reads
trimmomatic \
PE \
-threads 8 \
"${sample_id}_1.fastq.gz" \
"${sample_id}_2.fastq.gz" \
"${sample_id}_1p.fastq.gz" \
"${sample_id}_1u.fastq.gz" \
"${sample_id}_2p.fastq.gz" \
"${sample_id}_2u.fastq.gz" \
ILLUMINACLIP:"$condapath/share/trimmomatic/adapters/NexteraPE-PE.fa":2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36

# Cleanup to save space
rm "${sample_id}_1.fastq.gz" "${sample_id}_2.fastq.gz"

F1="${sample_id}_1p.fastq.gz"
F2="${sample_id}_2p.fastq.gz"
FU1="${sample_id}_1u.fastq.gz"
FU2="${sample_id}_2u.fastq.gz"
#FU="${sample_id}_u.fastq.gz"

Fpaired=("$F1" "$F2")
Funpaired=("$FU1","$FU2")

for variable in paired unpaired; do
eval 'readfiles=("${F'$variable'[@]}")'
STAR \
--readFilesIn "${readfiles[@]}" \
--alignIntronMax 1000000 \
--alignIntronMin 20 \
--alignMatesGapMax 1000000 \
--alignSJDBoverhangMin 1 \
--alignSJoverhangMin 8 \
--alignSoftClipAtReferenceEnds Yes \
--chimJunctionOverhangMin 15 \
--chimMainSegmentMultNmax 1 \
--chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
--chimSegmentMin 15 \
--genomeDir "${REFERENCE}" \
--genomeLoad NoSharedMemory \
--limitSjdbInsertNsj 1200000 \
--outBAMcompression 10 \
--outFileNamePrefix "${sample_id}_${variable}" \
--outFilterIntronMotifs None \
--outFilterMatchNminOverLread 0.33 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--outFilterMultimapNmax 20 \
--outFilterScoreMinOverLread 0.33 \
--outFilterType BySJout \
--outSAMattributes NH HI AS nM NM ch \
--outSAMstrandField intronMotif \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--readFilesCommand zcat \
--runThreadN 8 \
--twopassMode Basic

if [ "$variable" = "paired" ]; then
	rm "${readfiles[@]}"
elif [ "$variable" = "unpaired" ]; then
        rm "$FU1" "$FU2" 
fi
done

samtools merge - "${sample_id}_pairedAligned.out.bam" "${sample_id}_unpairedAligned.out.bam" |
samtools sort -@ 8 -o "${sample_id}Aligned.out.bam" -
rm "${sample_id}_pairedAligned.out.bam" "${sample_id}_unpairedAligned.out.bam"
samtools index "${sample_id}Aligned.out.bam"

primary_contig='chr6:28509970-33480727'
no_contig='*'

samtools view -o "${sample_id}.bam" -b "${sample_id}Aligned.out.bam" $primary_contig "$no_contig"
samtools index "${sample_id}.bam"
rm "${sample_id}Aligned.out.bam"

cp "${sample_id}.bam" "${sample_id}.bam.bai" "$output_folder/"
