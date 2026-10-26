"""Build the per-cell quality report of one pseudobulk (timepoint): RNA per-cell QC metrics from
the pseudobulk RNA h5ad (.obs, computed by qc/channel{1,2}_qc.py) merged with the ATAC per-cell QC
metrics from compute_atac_percell_qc.R, one row per cell. Fails if the two inputs do not hold
exactly the same barcodes.

Usage:
    python build_percell_qc_report.py <day> <rna_<day>.h5ad> <atac_percell_qc_<day>.tsv.gz> <out.tsv.gz>

Output columns: barcode, atac_fragment_count, atac_tss_enrichment, atac_nucleosome_signal,
rna_umi_count, rna_gene_count, rna_pct_mito, rna_scrublet_doublet_score, pass_strict_qc
(described in FILE_SPEC_PERCELL_QC_REPORT.txt).

Requires h5py and pandas.
"""

import sys

import h5py
import pandas as pd

if len(sys.argv) != 5:
    sys.exit(__doc__)
day, h5ad_path, atac_qc_path, out_path = sys.argv[1:5]

f = h5py.File(h5ad_path, "r")
obs = f["obs"]
rna = pd.DataFrame(
    {
        "barcode": [x.decode() for x in obs["barcode"][:]],
        "rna_umi_count": obs["total_counts"][:],
        "rna_gene_count": obs["n_genes_by_counts"][:],
        "rna_pct_mito": obs["pct_counts_mt"][:],
        "rna_scrublet_doublet_score": obs["scrublet_doublet_score"][:],
        "pass_strict_qc": obs["pass_strict_qc"][:],
    }
)

atac = pd.read_csv(atac_qc_path, sep="\t")

merged = rna.merge(atac, on="barcode", how="inner")
assert len(merged) == len(rna) == len(atac), (
    f"[{day}] expected exact barcode match, got rna={len(rna)} atac={len(atac)} merged={len(merged)}"
)

cols = [
    "barcode",
    "atac_fragment_count",
    "atac_tss_enrichment",
    "atac_nucleosome_signal",
    "rna_umi_count",
    "rna_gene_count",
    "rna_pct_mito",
    "rna_scrublet_doublet_score",
    "pass_strict_qc",
]
merged = merged[cols]
merged.to_csv(out_path, sep="\t", index=False, compression="gzip")
print(f"[{day}] wrote {len(merged)} cells -> {out_path}")
