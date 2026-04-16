#!/bin/bash
#SBATCH -J atac_align
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=bigmem
#SBATCH --mem=200G
#SBATCH --mail-type=ALL

set -euo pipefail
source config/paths.sh

SAMPLE_ID="$1"
STAGE="$2"
REF_NAME="$3"

BASE_DIR="${PROJECT_ROOT}/${STAGE}"

if [[ "$REF_NAME" == "hmm10" ]]; then
    REF="$HUMAN_REF"
else
    REF="$CHIMP_REF"
fi

module load Bowtie2

mkdir -p "${BASE_DIR}/Aligned/${SAMPLE_ID}_alignment"

bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 10 \
    -x "$REF" \
    -1 "${BASE_DIR}/Trimmed/${SAMPLE_ID}_trimmed/${SAMPLE_ID}_R1.fastq" \
    -2 "${BASE_DIR}/Trimmed/${SAMPLE_ID}_trimmed/${SAMPLE_ID}_R2.fastq" \
    -S "${BASE_DIR}/Aligned/${SAMPLE_ID}_alignment/${SAMPLE_ID}_bowtie2.sam" \
    &> "${BASE_DIR}/Aligned/${SAMPLE_ID}_alignment/${SAMPLE_ID}_bowtie2.txt"
