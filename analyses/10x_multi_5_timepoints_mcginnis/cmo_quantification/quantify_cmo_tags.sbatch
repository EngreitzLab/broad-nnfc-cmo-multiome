#!/bin/bash
#SBATCH -J cmo-quantification-mcginnis
#SBATCH -p normal
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH -o /oak/stanford/groups/engreitz/Users/emattei/git/broad-nnfc-cmo-multiome/log/cmo_quant_mcginnis_%j.log
#SBATCH -e /oak/stanford/groups/engreitz/Users/emattei/git/broad-nnfc-cmo-multiome/log/cmo_quant_mcginnis_%j.err

set -euo pipefail

source /home/users/emattei/miniforge3/etc/profile.d/conda.sh
conda activate /oak/stanford/groups/engreitz/Users/emattei/git/broad-nnfc-cmo-multiome/env/nnfc-cmo-multiome

SCRIPT="/oak/stanford/groups/engreitz/Users/emattei/git/broad-nnfc-cmo-multiome/analyses/10x_multi_5_timepoints_mcginnis/cmo_quantification/quantify_cmo_tags.py"
WORK_DIR="${L_SCRATCH}/cmo_quant_mcginnis"

python "${SCRIPT}" \
    --work-dir "${WORK_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 8G
