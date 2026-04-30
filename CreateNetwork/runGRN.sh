#!/bin/bash
#SBATCH -J runGRN
#SBATCH -t 1-00:00:00
#SBATCH --partition=day
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --array=1-160

module load pySCENIC/0.12.1-20240311-foss-2022b

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

BASE_DIR="/vast/palmer/scratch/noonan/yj345/snRNA/grn"
TISSUES=(FL HL PA1 PA2)
N_SUB=40

TASK_ID=${SLURM_ARRAY_TASK_ID}

TISSUE_INDEX=$(( (TASK_ID - 1) / N_SUB ))
SUB_ID=$(( (TASK_ID - 1) % N_SUB + 1 ))

TISSUE=${TISSUES[$TISSUE_INDEX]}
ID_PAD=$(printf "%02d" "${SUB_ID}")

SUB_DIR="${BASE_DIR}/${TISSUE}/sub"
OUT_DIR="${BASE_DIR}/${TISSUE}/adj_sub"
TF_LIST="${BASE_DIR}/${TISSUE}/tf_${TISSUE}.txt"

mkdir -p "${OUT_DIR}"

IN_LOOM="${SUB_DIR}/${TISSUE}_E10_magic_sub_${ID_PAD}_3k.loom"
OUT_TSV="${OUT_DIR}/adj_sub_${SUB_ID}.tsv"

echo "Running pySCENIC GRN for tissue: ${TISSUE}"

pyscenic grn \
    --num_workers 10 \
    --output "${OUT_TSV}" \
    --method grnboost2 \
    "${IN_LOOM}" \
    "${TF_LIST}"
