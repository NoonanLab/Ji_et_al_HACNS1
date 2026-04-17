#!/bin/bash
#SBATCH -J PIPbarcode
#SBATCH -t 1-00:00:00
#SBATCH --nodes=10
#SBATCH --cpus-per-task=25
#SBATCH --mem=110G
#SBATCH --mail-type=ALL

set -euo pipefail

module load PIPseeker/2.1.4

SAMPLE_SHEET="config/samples.tsv"

tail -n +2 "$SAMPLE_SHEET" | while IFS=$'\t' read -r sample_id stage tissue root_path fastq_dir; do
    echo "Initiating PIPseeker barcode for ${sample_id}:"
    date

    mkdir -p "${root_path}/${sample_id}_out"

    pipseeker barcode \
        --threads 25 \
        --fastq "${fastq_dir}/." \
        --output-path "${root_path}/${sample_id}_out" \
        --chemistry v4 \
        --verbosity 2

    echo "Finished PIPseeker barcode for ${sample_id}"
    date
done
