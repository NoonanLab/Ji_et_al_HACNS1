#!/bin/bash
#SBATCH -J PIPVelocyto
#SBATCH -t 1-00:00:00
#SBATCH --nodes=10
#SBATCH --cpus-per-task=10
#SBATCH --mem=110G
#SBATCH --mail-type=ALL

set -euo pipefail

source config/velocyto_paths.sh

module load "$STAR_MODULE"

SAMPLE_SHEET="config/velocyto_samples.tsv"

tail -n +2 "$SAMPLE_SHEET" | while IFS=$'\t' read -r sample_id stage tissue root_path; do
    echo "${sample_id} STAR alignment begun:"
    date

    mkdir -p "${OUT_ROOT}/${tissue}/${stage}/${sample_id}_Velocyto"

    STAR --runThreadN "${THREADS}" \
        --genomeDir "${GENOME_DIR}/" \
        --readFilesPrefix "${root_path}/${sample_id}_out/barcoded_fastqs/" \
        --readFilesIn \
barcoded_1_R2.fastq.gz,barcoded_2_R2.fastq.gz,barcoded_3_R2.fastq.gz,barcoded_4_R2.fastq.gz,barcoded_5_R2.fastq.gz,barcoded_6_R2.fastq.gz \
barcoded_1_R1.fastq.gz,barcoded_2_R1.fastq.gz,barcoded_3_R1.fastq.gz,barcoded_4_R1.fastq.gz,barcoded_5_R1.fastq.gz,barcoded_6_R1.fastq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix "${OUT_ROOT}/${tissue}/${stage}/${sample_id}_Velocyto/" \
        --outSAMtype BAM Unsorted \
        --outSAMattributes NH HI AS nM jM jI CR CY UR UY GX GN \
        --outSAMprimaryFlag AllBestScore \
        --outSAMmultNmax 10 \
        --outBAMcompression 10 \
        --soloFeatures Gene GeneFull SJ Velocyto \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "${root_path}/${sample_id}_out/barcodes/barcode_whitelist.txt" \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 12 \
        --soloBarcodeReadLength 1 \
        --soloCBmatchWLtype 1MM_multi \
        --soloUMIdedup 1MM_All \
        --soloUMIfiltering - \
        --soloCellFilter CellRanger2.2

    echo "${sample_id} STAR alignment finished:"
    date
done
