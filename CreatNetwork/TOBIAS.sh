#!/bin/bash
#SBATCH -J Footprint
#SBATCH -t 1-00:00:00
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=500G
#SBATCH --mail-type=ALL

module load miniconda
conda activate TOBIAS

GENOME="/gpfs/gibbs/project/noonan/yj345/makeHumanizedFasta/mm10.fa"
MOTIFS="/gpfs/gibbs/project/noonan/yj345/mm10_fasta_gtf/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt"
ORIGIN="/vast/palmer/scratch/noonan/yj345/ATAC/motif2gene_mouse.tsv"
BASE="/vast/palmer/scratch/noonan/yj345/ATAC"

for tissue in FL HL P1 P2; do

    echo "Running TOBIAS for ${tissue}"
  
    cd ${BASE}/${tissue}

    PEAKS="${tissue}_E10_Merged_peaks_annotated_withEnh.bed"

    HUMAN="H${tissue}E10merged"
    CHIMP="C${tissue}E10merged"

    TOBIAS ATACorrect \
        --bam ${HUMAN}/${HUMAN}_callpeak.bam \
        --genome ${GENOME} \
        --peaks ${PEAKS} \
        --outdir ${HUMAN}/ATACorrect \
        --cores 10

    TOBIAS FootprintScores \
        --signal ${HUMAN}/ATACorrect/*_corrected.bw \
        --regions ${PEAKS} \
        --output ${HUMAN}/${HUMAN}_footprints.bw \
        --cores 10

    TOBIAS ATACorrect \
        --bam ${CHIMP}/${CHIMP}_callpeak.bam \
        --genome ${GENOME} \
        --peaks ${PEAKS} \
        --outdir ${CHIMP}/ATACorrect \
        --cores 10

    TOBIAS FootprintScores \
        --signal ${CHIMP}/ATACorrect/*_corrected.bw \
        --regions ${PEAKS} \
        --output ${CHIMP}/${CHIMP}_footprints.bw \
        --cores 10

    TOBIAS BINDetect \
        --motifs ${MOTIFS} \
        --signals ${HUMAN}/${HUMAN}_footprints.bw ${CHIMP}/${CHIMP}_footprints.bw \
        --genome ${GENOME} \
        --peaks ${PEAKS} \
        --outdir BINDetect_out \
        --cond_names Human Chimp \
        --cores 10

    TOBIAS CreateNetwork \
        --TFBS BINDetect_out/*/beds/*_bound.bed \
        --origin ${ORIGIN} \
        --outdir CreateNetwork_out

done
