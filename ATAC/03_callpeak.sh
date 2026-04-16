#!/bin/bash
#SBATCH -J atac_peak
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=bigmem
#SBATCH --mem=200G
#SBATCH --mail-type=ALL

set -euo pipefail
source config/paths.sh

GROUP="$1"
STAGE="$2"
BASE_DIR="${PROJECT_ROOT}/${STAGE}"

module load MACS2

mkdir -p "${BASE_DIR}/Peaks/${GROUP}_peak"

for rep in R1 R2 merged; do
    mkdir -p "${BASE_DIR}/Peaks/${GROUP}_peak/${rep}"
    macs2 callpeak \
        -t "${BASE_DIR}/Cleaned/${GROUP}${rep}_out/${GROUP}${rep}_callpeak.bam" \
        -f BAMPE -n "${GROUP}${rep}_macs2_peak" -g mm -q 0.1 \
        --outdir "${BASE_DIR}/Peaks/${GROUP}_peak/${rep}" --keep-dup all
done
