# Per-cell QC reports

One per-cell quality report per pseudobulk (timepoint d0-d4), describing exactly the cells in that
pseudobulk. It runs after `merge/merge_channels.py`, on its outputs: the per-timepoint RNA h5ad and
the per-timepoint ATAC fragments (sorted, bgzip-compressed and tabix-indexed).

| Script | Description |
|---|---|
| `compute_atac_percell_qc.R` | ATAC metrics per barcode from the pseudobulk fragments: fragment count and TSS enrichment (ArchR, with `minTSS=0`, `minFrags=0` so every barcode is kept) and nucleosome signal (ratio of 147-294 bp to <147 bp fragments) |
| `build_percell_qc_report.py` | Adds the RNA metrics stored in the pseudobulk h5ad `.obs` (computed by `qc/channel{1,2}_qc.py`) and writes the report; fails unless the RNA and ATAC tables hold exactly the same barcodes |
| `compute_percell_qc_reports.sbatch` | Runs both for the 5 timepoints (Slurm array 0-4) |
| `FILE_SPEC_PERCELL_QC_REPORT.txt` | Column definitions of the report |

## Inputs (per timepoint)

- `rna_<day>.h5ad`: pseudobulk RNA matrix from `merge/merge_channels.py`
- `<day>.fragments.tsv.gz`: pseudobulk ATAC fragments from `merge/merge_channels.py`
- Gene bounds BED: `resources/genome_annotations/CollapsedGeneBounds.hg38.intGENCODEv43.bed` from
  [scE2G](https://github.com/EngreitzLab/scE2G) (used for TSS positions)

## Outputs (per timepoint)

- `atac_percell_qc_<day>.tsv.gz`: barcode, atac_fragment_count, atac_tss_enrichment, atac_nucleosome_signal
- `percell_qc_report_<day>.tsv.gz`: barcode, the 3 ATAC columns, rna_umi_count, rna_gene_count,
  rna_pct_mito, rna_scrublet_doublet_score, pass_strict_qc (see `FILE_SPEC_PERCELL_QC_REPORT.txt`)

## Usage

```bash
export PSEUDOBULK_DIR=/path/to/Pseudobulks     # contains RNA/rna_<day>.h5ad and ATAC/<day>.fragments.tsv.gz
export GENE_BOUNDS_BED=/path/to/CollapsedGeneBounds.hg38.intGENCODEv43.bed
export OUT_DIR=/path/to/qc_reports
export PYTHON=/path/to/python                  # with h5py and pandas
# R with ArchR, BSgenome.Hsapiens.UCSC.hg38 and data.table on the library path
sbatch --export=ALL analyses/10x_multi_5_timepoints/qc_reports/compute_percell_qc_reports.sbatch
```

Or for one timepoint:

```bash
Rscript compute_atac_percell_qc.R d0 $PSEUDOBULK_DIR/ATAC/d0.fragments.tsv.gz $GENE_BOUNDS_BED $OUT_DIR/atac_percell_qc_d0.tsv.gz
python build_percell_qc_report.py d0 $PSEUDOBULK_DIR/RNA/rna_d0.h5ad $OUT_DIR/atac_percell_qc_d0.tsv.gz $OUT_DIR/percell_qc_report_d0.tsv.gz
```

## Dependencies

Versions used: R 4.4.2 with ArchR 1.0.3, BSgenome.Hsapiens.UCSC.hg38 1.4.5 and data.table 1.18.2.1;
Python 3.10 with pandas 2.3.3 and h5py 3.15.1.
