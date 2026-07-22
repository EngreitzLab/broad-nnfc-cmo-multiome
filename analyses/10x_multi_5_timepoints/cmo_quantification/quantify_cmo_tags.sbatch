#!/bin/bash
#SBATCH -J cmo-quantification
#SBATCH -p engreitz
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH -o /oak/stanford/groups/engreitz/Users/emattei/git/broad-nnfc-cmo-multiome/log/cmo_quant_%j.log
#SBATCH -e /oak/stanford/groups/engreitz/Users/emattei/git/broad-nnfc-cmo-multiome/log/cmo_quant_%j.err

set -euo pipefail

source /home/users/emattei/miniforge3/etc/profile.d/conda.sh
conda activate /oak/stanford/groups/engreitz/Users/emattei/git/broad-nnfc-cmo-multiome/env/nnfc-cmo-multiome

SCRIPT="/oak/stanford/groups/engreitz/Users/emattei/git/broad-nnfc-cmo-multiome/analyses/10x_multi_5_timepoints/cmo_quantification/quantify_cmo_tags.py"
CHANNEL="${1:-}"  # optional: channel1 or channel2; omit to run both
WORK_DIR="${L_SCRATCH}/cmo_quant${CHANNEL:+_${CHANNEL}}"

CHANNEL_ARG=""
if [ -n "${CHANNEL}" ]; then
    CHANNEL_ARG="--channels ${CHANNEL}"
fi

python "${SCRIPT}" \
    --work-dir "${WORK_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 8G \
    ${CHANNEL_ARG}
