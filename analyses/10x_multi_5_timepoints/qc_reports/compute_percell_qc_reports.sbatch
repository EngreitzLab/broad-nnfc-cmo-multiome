#!/bin/bash
#SBATCH --job-name=percell_qc_reports
#SBATCH --array=0-4
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%A_%a.log
# Per-cell QC report for each pseudobulk (one array task per timepoint): ATAC metrics with
# compute_atac_percell_qc.R, then merged with the RNA metrics by build_percell_qc_report.py.
#
# Set before submitting (sbatch --export=ALL ...):
#   PSEUDOBULK_DIR   directory with RNA/rna_<day>.h5ad and ATAC/<day>.fragments.tsv.gz
#   GENE_BOUNDS_BED  scE2G resources/genome_annotations/CollapsedGeneBounds.hg38.intGENCODEv43.bed
#   OUT_DIR          output directory
#   PYTHON           python with h5py and pandas (default: python)
#   WORK_DIR         ArchR working directory, e.g. node-local scratch (default: OUT_DIR)
#   QC_SCRIPTS_DIR   this folder (default: <submit dir>/analyses/10x_multi_5_timepoints/qc_reports)
# R with ArchR, BSgenome.Hsapiens.UCSC.hg38 and data.table must be available (e.g. module load +
# R_LIBS_USER).
set -euo pipefail
: "${PSEUDOBULK_DIR:?}" "${GENE_BOUNDS_BED:?}" "${OUT_DIR:?}"
PYTHON=${PYTHON:-python}
# Slurm runs a spooled copy of this script, so locate the scripts from the submit directory
# (submit from the repository root) unless QC_SCRIPTS_DIR is set.
HERE=${QC_SCRIPTS_DIR:-${SLURM_SUBMIT_DIR}/analyses/10x_multi_5_timepoints/qc_reports}
DAYS=(d0 d1 d2 d3 d4)
DAY=${DAYS[$SLURM_ARRAY_TASK_ID]}
mkdir -p "$OUT_DIR"

Rscript "$HERE/compute_atac_percell_qc.R" "$DAY" "$PSEUDOBULK_DIR/ATAC/${DAY}.fragments.tsv.gz" \
    "$GENE_BOUNDS_BED" "$OUT_DIR/atac_percell_qc_${DAY}.tsv.gz" "${WORK_DIR:-$OUT_DIR}/archr_work_${DAY}"
"$PYTHON" "$HERE/build_percell_qc_report.py" "$DAY" "$PSEUDOBULK_DIR/RNA/rna_${DAY}.h5ad" \
    "$OUT_DIR/atac_percell_qc_${DAY}.tsv.gz" "$OUT_DIR/percell_qc_report_${DAY}.tsv.gz"
