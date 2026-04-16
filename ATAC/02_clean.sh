#!/bin/bash
#SBATCH -J atac_clean
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=bigmem
#SBATCH --mem=200G
#SBATCH --mail-type=ALL

set -euo pipefail
source config/paths.sh

SAMPLE_ID="$1"
STAGE="$2"
BASE_DIR="${PROJECT_ROOT}/${STAGE}"

module load SAMtools
module load picard

BAM_DIR="${BASE_DIR}/Cleaned/${SAMPLE_ID}_out"
mkdir -p "$BAM_DIR"

samtools view -S -b "${BASE_DIR}/Aligned/${SAMPLE_ID}_alignment/${SAMPLE_ID}_bowtie2.sam" > "$BAM_DIR/${SAMPLE_ID}_bowtie2.bam"
samtools sort -o "$BAM_DIR/${SAMPLE_ID}_sorted.bam" "$BAM_DIR/${SAMPLE_ID}_bowtie2.bam"
samtools index "$BAM_DIR/${SAMPLE_ID}_sorted.bam"
samtools view -h "$BAM_DIR/${SAMPLE_ID}_sorted.bam" | grep -v chrM | samtools view -b - > "$BAM_DIR/${SAMPLE_ID}_rmChrM.bam"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    -I "$BAM_DIR/${SAMPLE_ID}_rmChrM.bam" \
    -O "$BAM_DIR/${SAMPLE_ID}_rmDup.bam" \
    -METRICS_FILE "$BAM_DIR/${SAMPLE_ID}_picard.dupMark.txt" \
    -REMOVE_DUPLICATES true
samtools view -b -q 10 "$BAM_DIR/${SAMPLE_ID}_rmDup.bam" > "$BAM_DIR/${SAMPLE_ID}_callpeak.bam"
