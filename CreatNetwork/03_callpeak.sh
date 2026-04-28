#!/bin/bash
#SBATCH -J peak
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=bigmem
#SBATCH --mem=200G
#SBATCH --mail-type=ALL

set -euo pipefail

GROUP="$1"
STAGE="$2"
BASE_DIR="/vast/palmer/scratch/noonan/yj345/${STAGE}"
GENOME_SIZE="/gpfs/gibbs/pi/noonan/yj345/genome/mm10.chrom.sizes"
CONTIGS_MOUSE="/vast/palmer/scratch/noonan/yj345/contigs_mouse.txt"

module load MACS2
module load SAMtools
module load Kent_tools
module load BEDTools

mkdir -p "${BASE_DIR}/Peaks/${GROUP}_peak"
mkdir -p "${BASE_DIR}/Final"

# Decide whether this group is single-end ChIP-seq or paired-end ATAC-style data
if [[ "$GROUP" =~ K27ac$ ]]; then
    MACS_FORMAT="BAM"
else
    MACS_FORMAT="BAMPE"
fi

# Build list of available members
AVAILABLE_REPS=()

for rep in R1 R2; do
    if [[ -f "${BASE_DIR}/Cleaned/${GROUP}${rep}_out/${GROUP}${rep}_callpeak.bam" ]]; then
        AVAILABLE_REPS+=("${rep}")
    fi
done

# Merge only if both R1 and R2 exist
if [[ -f "${BASE_DIR}/Cleaned/${GROUP}R1_out/${GROUP}R1_callpeak.bam" && -f "${BASE_DIR}/Cleaned/${GROUP}R2_out/${GROUP}R2_callpeak.bam" ]]; then
    mkdir -p "${BASE_DIR}/Cleaned/${GROUP}merged_out"
    samtools merge "${BASE_DIR}/Cleaned/${GROUP}merged_out/${GROUP}merged_callpeak.bam" \
        "${BASE_DIR}/Cleaned/${GROUP}R1_out/${GROUP}R1_callpeak.bam" \
        "${BASE_DIR}/Cleaned/${GROUP}R2_out/${GROUP}R2_callpeak.bam"
    samtools index "${BASE_DIR}/Cleaned/${GROUP}merged_out/${GROUP}merged_callpeak.bam"
    AVAILABLE_REPS+=("merged")
fi

# Peak calling
for rep in "${AVAILABLE_REPS[@]}"; do
    mkdir -p "${BASE_DIR}/Peaks/${GROUP}_peak/${rep}"

    macs2 callpeak \
        -t "${BASE_DIR}/Cleaned/${GROUP}${rep}_out/${GROUP}${rep}_callpeak.bam" \
        -f "${MACS_FORMAT}" \
        -n "${GROUP}${rep}_macs2_peak" \
        -g mm \
        -q 0.1 \
        --outdir "${BASE_DIR}/Peaks/${GROUP}_peak/${rep}" \
        --keep-dup all
done

# BigWig conversion
for rep in "${AVAILABLE_REPS[@]}"; do
    SAMPLE="${GROUP}${rep}"
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
done

# BigBed conversion for MACS2 peaks
for rep in "${AVAILABLE_REPS[@]}"; do
    grep -f "$CONTIGS_MOUSE" \
        "${BASE_DIR}/Peaks/${GROUP}_peak/${rep}/${GROUP}${rep}_macs2_peak_peaks.narrowPeak" \
        > "${BASE_DIR}/Peaks/${GROUP}_peak/${rep}/${GROUP}${rep}.narrowPeak"

    bedToBigBed \
        -type=bed4+6 \
        "${BASE_DIR}/Peaks/${GROUP}_peak/${rep}/${GROUP}${rep}.narrowPeak" \
        "$GENOME_SIZE" \
        "${BASE_DIR}/Final/${GROUP}${rep}_peaks.bb"
done
