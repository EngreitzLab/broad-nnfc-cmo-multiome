#!/bin/bash
#SBATCH -J conda-env-build
#SBATCH -p normal
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH -o /oak/stanford/groups/engreitz/Users/emattei/git/broad-nnfc-cmo-multiome/log/build_%j.log
#SBATCH -e /oak/stanford/groups/engreitz/Users/emattei/git/broad-nnfc-cmo-multiome/log/build_%j.err

set -euo pipefail

source /home/users/emattei/miniforge3/etc/profile.d/conda.sh

PROJECT_DIR="/oak/stanford/groups/engreitz/Users/emattei/git/broad-nnfc-cmo-multiome"
ENV_NAME="nnfc-cmo-multiome"
ENV_PREFIX="$PROJECT_DIR/env/$ENV_NAME"
TARGET_DIR="$PROJECT_DIR/config/conda"

# Use L_SCRATCH for fast package cache only
export CONDA_PKGS_DIRS="$L_SCRATCH/conda_pkgs"
mkdir -p "$ENV_PREFIX" "$CONDA_PKGS_DIRS"

echo "[$(date)] Building conda env at $ENV_PREFIX"

conda create \
    --prefix "$ENV_PREFIX" \
    -c bioconda -c conda-forge \
    --yes \
    python=3.12 \
    splitcode \
    "kb-python=0.29.1"

conda activate "$ENV_PREFIX"

echo "[$(date)] Installing pip packages"
pip install "seqspec==0.3.0" igvf-utils marimo

echo "[$(date)] Exporting frozen environment spec"
conda env export --prefix "$ENV_PREFIX" | gzip > "$TARGET_DIR/environment.yaml.gz"

echo "[$(date)] Done. Environment: $ENV_PREFIX"
echo "[$(date)] Frozen spec: $TARGET_DIR/environment.yaml.gz"
