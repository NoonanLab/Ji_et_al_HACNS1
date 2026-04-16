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
module load SAMtools
module load Kent_tools
module load BEDTools

echo "$GROUP Peak Calling begun:" && date
mkdir -p "${BASE_DIR}/Peaks/${GROUP}_peak"
mkdir -p "${BASE_DIR}/Final"

# Peak calling for R1, R2, and merged
for rep in R1 R2 merged; do
    mkdir -p "${BASE_DIR}/Peaks/${GROUP}_peak/${rep}"
    macs2 callpeak \
        -t "${BASE_DIR}/Cleaned/${GROUP}${rep}_out/${GROUP}${rep}_callpeak.bam" \
        -f BAMPE \
        -n "${GROUP}${rep}_macs2_peak" \
        -g mm \
        -q 0.1 \
        --outdir "${BASE_DIR}/Peaks/${GROUP}_peak/${rep}" \
        --keep-dup all
done
echo "$GROUP Peak Calling finished:" && date

# Conversion to BigWig
for rep in R1 R2 merged; do
    SAMPLE="${GROUP}${rep}"
    echo "$SAMPLE Format Conversion begun:" && date

    mkdir -p "${BASE_DIR}/Format/${SAMPLE}_format_bw"

    samtools sort \
        -o "${BASE_DIR}/Format/${SAMPLE}_format_bw/${SAMPLE}_sorted.bam" \
        "${BASE_DIR}/Cleaned/${SAMPLE}_out/${SAMPLE}_callpeak.bam"

    samtools index "${BASE_DIR}/Format/${SAMPLE}_format_bw/${SAMPLE}_sorted.bam"

    bedtools genomecov \
        -ibam "${BASE_DIR}/Format/${SAMPLE}_format_bw/${SAMPLE}_sorted.bam" \
        -bg -scale 1 \
        > "${BASE_DIR}/Format/${SAMPLE}_format_bw/${SAMPLE}.bedgraph"

    grep -f "$CONTIGS_MOUSE" \
        "${BASE_DIR}/Format/${SAMPLE}_format_bw/${SAMPLE}.bedgraph" \
        > "${BASE_DIR}/Format/${SAMPLE}_format_bw/${SAMPLE}_cleaned.bedgraph"

    bedGraphToBigWig \
        "${BASE_DIR}/Format/${SAMPLE}_format_bw/${SAMPLE}_cleaned.bedgraph" \
        "$GENOME_SIZE" \
        "${BASE_DIR}/Final/${SAMPLE}.bw"

    echo "$SAMPLE Format Conversion finished:" && date
done

# Convert Peaks to BigBed
echo "$GROUP Peak Format Conversion begun:" && date
for rep in R1 R2 merged; do
    grep -f "$CONTIGS_MOUSE" \
        "${BASE_DIR}/Peaks/${GROUP}_peak/${rep}/${GROUP}${rep}_macs2_peak_peaks.narrowPeak" \
        > "${BASE_DIR}/Peaks/${GROUP}_peak/${rep}/${GROUP}${rep}.narrowPeak"

    bedToBigBed \
        -type=bed4+6 \
        "${BASE_DIR}/Peaks/${GROUP}_peak/${rep}/${GROUP}${rep}.narrowPeak" \
        "$GENOME_SIZE" \
        "${BASE_DIR}/Final/${GROUP}${rep}_peaks.bb"
done
echo "$GROUP Peak Format Conversion finished:" && date
