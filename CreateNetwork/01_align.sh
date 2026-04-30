#!/bin/bash
#SBATCH -J align
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=bigmem
#SBATCH --mem=200G
#SBATCH --mail-type=ALL

set -euo pipefail

SAMPLE_ID="$1"
STAGE="$2"
BASE_DIR="/vast/palmer/scratch/noonan/yj345/${STAGE}"
MM10_REF="/gpfs/gibbs/project/noonan/yj345/mm10"

module load Bowtie2

mkdir -p "${BASE_DIR}/Aligned/${SAMPLE_ID}_alignment"

if [[ "$SAMPLE_ID" =~ ^(SRR4835044|SRR4835045|SRX3395791|SRR3950341|SRR3950342)$ ]]; then
    bowtie2 --local --very-sensitive --phred33 -p 10 \
        -x "$MM10_REF" \
        -U "${BASE_DIR}/RawFastq/${SAMPLE_ID}.fastq" \
        -S "${BASE_DIR}/Aligned/${SAMPLE_ID}_alignment/${SAMPLE_ID}_bowtie2.sam" \
        &> "${BASE_DIR}/Aligned/${SAMPLE_ID}_alignment/${SAMPLE_ID}_bowtie2.txt"

else
    bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 10 \
        -x "$MM10_REF" \
        -1 "${BASE_DIR}/Trimmed/${SAMPLE_ID}_trimmed/${SAMPLE_ID}_R1.fastq" \
        -2 "${BASE_DIR}/Trimmed/${SAMPLE_ID}_trimmed/${SAMPLE_ID}_R2.fastq" \
        -S "${BASE_DIR}/Aligned/${SAMPLE_ID}_alignment/${SAMPLE_ID}_bowtie2.sam" \
        &> "${BASE_DIR}/Aligned/${SAMPLE_ID}_alignment/${SAMPLE_ID}_bowtie2.txt"
fi
