# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "altair>=6.1.0",
#     "igraph>=1.0.0",
#     "ipython>=9.13.0",
#     "marimo>=0.23.3",
#     "scanpy[scrublet]>=1.12.1",
#     "snapatac2>=2.9.0",
# ]
# ///

import marimo

__generated_with = "0.23.13"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Introduction

    This notebook processes the 10x multiome data for the endothelial differentiation time course: 5 timepoints (d0-d4), 5 biological replicates each. Samples were multiplexed using the MULTI-seq technique with CMO (Cell Multiplexing Oligo) barcodes. Data was processed with the IGVF pipeline using kallisto and the GENCODE v43 annotation, and CMO quantification was performed with the `kite` workflow from the `kallisto-bustools` suite.

    The goals of this notebook are to:

    1. Perform quality control filtering on the RNA data
    2. Run CMO hash classification (mimicking `Seurat::HTODemux`) on QC-passing cells
    3. Assign each cell barcode to a CMO (and therefore to a timepoint/replicate)
    4. Process the corresponding ATAC data with `snapatac2`
    """)
    return


@app.cell(hide_code=True)
def load_gene_metadata_intro(mo):
    mo.md(r"""
    ## Loading gene metadata

    The transcription start sites (TSSs) for the protein-coding genes annotated in GENCODE v43 were curated using the MANE annotations. The source file stores a 500bp window centered on each TSS (-250/+249); we collapse that down to the single-bp TSS position itself.

    Columns loaded:

    - `chr` - chromosome
    - `start` / `end` - single-bp TSS position (BED 0-based, half-open)
    - `gene_symbol` - gene symbol
    - `score` - unused placeholder column from the source BED file
    - `strand` - `+` or `-`
    - `gene_id` - Ensembl gene ID
    - `gene_type` - GENCODE gene biotype (all `protein_coding` here)
    """)
    return


@app.cell(hide_code=True)
def load_gene_metadata(pd, project_root):
    # The bed file stores a 500bp window around each gene's TSS (-250/+249,
    # BED 0-based half-open). The TSS itself, not "the first base of the
    # window", the base the window is centered on, is the window's midpoint.
    _gene_metadata_fnp = "annotations/gencode.v43.protein_coding.TSS500bp.bed"
    gene_metadata_df = pd.read_csv(
        project_root / _gene_metadata_fnp,
        sep="\t",
        header=0,
        names=["chr", "start", "end", "gene_symbol", "score", "strand", "gene_id", "gene_type"]
    )
    _tss_pos = (gene_metadata_df["start"] + gene_metadata_df["end"]) // 2
    gene_metadata_df["start"] = _tss_pos
    gene_metadata_df["end"] = _tss_pos + 1
    gene_metadata_df
    return (gene_metadata_df,)


@app.cell(hide_code=True)
def get_igvf_gencode_intro(mo):
    mo.md(r"""
    ## Downloading the GENCODE GTF (IGVF-hosted)
    The full GENCODE v43 GTF (used for the ATAC TSSE metric below) is downloaded from IGVF's reference-file API if not already present locally, so the notebook doesn't depend on a symlinked annotation file outside this repo.
    """)
    return


@app.cell
def get_igvf_gencode(project_root):
    import urllib.request

    igvf_gencode_gtf_path = project_root / "annotations/IGVFFI9573KOZR.gtf.gz"
    _igvf_gencode_url = "https://api.data.igvf.org/reference-files/IGVFFI9573KOZR/@@download/IGVFFI9573KOZR.gtf.gz"

    if not igvf_gencode_gtf_path.exists():
        igvf_gencode_gtf_path.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(_igvf_gencode_url, igvf_gencode_gtf_path)

    # Show only the path relative to the repo, not the full local filesystem path
    igvf_gencode_gtf_path.relative_to(project_root)
    return (igvf_gencode_gtf_path,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Loading the counts in scanpy

    We are going to load the total counts as produced by `kallisto-bustool` and annotate the genes with the metadata loaded the cell above. The default matrix stored in `.X` is the total counts: the sum of `ambigous`, `mature`, and `nascent` matrices.
    """)
    return


@app.cell
def load_h5_counts(ad, gene_metadata_df, project_root):
    _h5ad_fnp = "data/h5ad/10x_5timepoints_channel1.h5ad"
    adata = ad.read_h5ad(project_root / _h5ad_fnp)

    # -------------------------
    # Add gene metadata to var
    # -------------------------
    _original_var_names = adata.var_names.copy()

    adata.var["gene_id_base"] = _original_var_names.str.replace(r"\.\d+$", "", regex=True)

    adata.var = adata.var.merge(
        gene_metadata_df[
            ["gene_id", "gene_symbol", "gene_type", "chr", "start", "end", "strand"]
        ],
        left_on="gene_id_base",
        right_on="gene_id",
        how="left",
        sort=False
    )

    adata.var.index = _original_var_names
    adata.var.index.name = None

    # -------------------------
    # QC gene flags
    # -------------------------
    _gene_symbol = adata.var["gene_symbol"].fillna("")

    adata.var["mt"] = _gene_symbol.str.startswith("MT-")
    adata.var["ribo"] = _gene_symbol.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = _gene_symbol.str.contains(r"^HB(?!P)", regex=True)


    adata
    return (adata,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # RNA: Quality control
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Compute QC metrics
    """)
    return


@app.cell
def _(adata, sc):
    # calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Knee plot
    """)
    return


@app.cell
def plot_knee_plot(adata, alt, mo, np, pd):
    # Total UMI counts per barcode (from QC metrics, computed above). Drop
    # zero-count barcodes since they're undefined on a log-scale y-axis.
    _counts_sorted = np.sort(adata.obs["total_counts"].to_numpy())[::-1]
    _counts_sorted = _counts_sorted[_counts_sorted > 0]
    knee_full_counts = _counts_sorted
    _ranks = np.arange(1, len(_counts_sorted) + 1)

    # Log-spaced subsample for a light, smooth curve (full barcode set is >600k points)
    _n_points = 2000
    _log_idx = np.unique(np.logspace(0, np.log10(len(_ranks) - 1), _n_points).astype(int))
    _knee_df = pd.DataFrame({"rank": _ranks[_log_idx], "n_umis": _counts_sorted[_log_idx]})

    knee_umi_selection = alt.selection_interval(encodings=["y"], value={"y": [500, 10_000]})

    _knee_chart = alt.Chart(_knee_df).mark_circle(size=15, color="black").encode(
        x=alt.X("rank:Q", scale=alt.Scale(type="log"), title="Cell rank"),
        y=alt.Y("n_umis:Q", scale=alt.Scale(type="log"), title="Total UMIs"),
        tooltip=[alt.Tooltip("rank:Q", title="Rank"), alt.Tooltip("n_umis:Q", title="# UMIs", format=",")],
    ).add_params(knee_umi_selection).properties(width=550, height=400, title="Knee plot")

    knee_chart_ui = mo.ui.altair_chart(_knee_chart, chart_selection=False, legend_selection=False)
    knee_chart_ui
    return knee_chart_ui, knee_full_counts


@app.cell(hide_code=True)
def knee_plot_brush_count(knee_chart_ui, knee_full_counts, mo):
    # Read the brush's y-range (falls back to the default 500-10,000 before any
    # interaction, since the interval's initial "value" isn't reported back until
    # the user actually drags it) and count how many barcodes it covers, using the
    # full (non-downsampled) counts array for an accurate number.
    _selections = knee_chart_ui.selections
    if _selections:
        _lo, _hi = next(iter(_selections.values()))["n_umis"]
    else:
        _lo, _hi = 500, 10_000

    _n_in_range = int(((knee_full_counts >= _lo) & (knee_full_counts <= _hi)).sum())

    mo.md(f"**{_n_in_range:,} barcodes** fall between **{_lo:,.0f}** and **{_hi:,.0f}** total UMIs (drag the brush above to adjust).")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Violin plots for different metrics
    """)
    return


@app.cell
def basic_qc_violin_plots(adata, mo, sc):
    _metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"]
    _labels = ["Number of genes", "Total UMI counts", "% mitochondrial", "% ribosomal", "% hemoglobin"]

    _panels = []
    for _metric, _label in zip(_metrics, _labels):
        _ax = sc.pl.violin(adata, _metric, jitter=0.4, ylabel=_label, show=False)
        _ax.set_title(_label)
        _panels.append(_ax.figure)

    mo.carousel(_panels)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### UMIs from protein-coding genes
    """)
    return


@app.cell(hide_code=True)
def pct_protein_coding_dist(adata, np, sc):
    _pc_mask = (adata.var["gene_type"] == "protein_coding").to_numpy()
    _pc_counts = np.asarray(adata.X[:, _pc_mask].sum(axis=1)).ravel()
    _total_counts = adata.obs["total_counts"].to_numpy()

    adata.obs["pct_counts_pc"] = np.divide(
        _pc_counts, _total_counts, out=np.zeros_like(_pc_counts), where=_total_counts > 0
    ) * 100

    sc.pl.violin(adata, "pct_counts_pc", jitter=0.4, ylabel="% UMIs from protein-coding genes")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Flag low quality barcodes
    """)
    return


@app.cell(hide_code=True)
def mask_filter_cells(adata, mo, pd):
    # QC pass mask mirroring the cell-level filters applied below
    _max_counts_cutoff = 10_000

    adata.obs["pass_min_umi_filter"] = adata.obs["total_counts"] >= 500
    adata.obs["pass_max_umi_filter"] = adata.obs["total_counts"] <= _max_counts_cutoff
    adata.obs["pass_min_gene_filter"] = adata.obs["n_genes_by_counts"] >= 200

    _steps = [
        ("min_counts >= 500", adata.obs["pass_min_umi_filter"]),
        (f"max_counts <= {_max_counts_cutoff}", adata.obs["pass_max_umi_filter"]),
        ("min_genes >= 200", adata.obs["pass_min_gene_filter"]),
    ]

    _remaining_mask = pd.Series(True, index=adata.obs.index)
    _rows = []
    for _label, _step_mask in _steps:
        _before_all = int(_remaining_mask.sum())
        _remaining_mask &= _step_mask
        _after_all = int(_remaining_mask.sum())
        _lost_all = _before_all - _after_all
        _rows.append(
            f"| `{_label}` | {_before_all:,} → {_after_all:,} "
            f"(lost {_lost_all:,}, {_lost_all / _before_all:.1%}) |"
        )

    adata.obs["pass_lvl_1_qc"] = _remaining_mask

    mo.md(f"""
    **QC filter breakdown** (min_counts=500, max_counts={_max_counts_cutoff:,}, min_genes=200)

    | Filter step | All barcodes |
    |---|---|
    {chr(10).join(_rows)}

    **Net:** {int(adata.obs["pass_lvl_1_qc"].sum()):,} of {adata.n_obs:,} barcodes survive this QC filter.
    """)
    return


@app.cell(hide_code=True)
def checking_mito_content(adata, plt, sns):
    _df = adata.obs.loc[adata.obs["pass_lvl_1_qc"], ["total_counts", "pct_counts_mt"]].copy()

    _g = sns.JointGrid(data=_df, x="total_counts", y="pct_counts_mt", height=6)
    _scatter = _g.ax_joint.scatter(
        _df["total_counts"], _df["pct_counts_mt"],
        c=_df["pct_counts_mt"], cmap="viridis", s=5, alpha=0.5,
    )

    _g.ax_marg_x.hist(_df["total_counts"], bins=100, color="gray", alpha=0.7)
    _g.ax_marg_y.hist(_df["pct_counts_mt"], bins=100, orientation="horizontal", color="gray", alpha=0.7)

    for _y in (5, 15):
        _g.ax_joint.axhline(y=_y, color="red", linestyle="--")
        _g.ax_joint.annotate(
            f"{_y}%",
            xy=(_df["total_counts"].max(), _y),
            xytext=(20, 5),
            textcoords="offset points",
            va="bottom",
            ha="right",
            color="red",
            fontweight="bold",
        )

    _g.ax_joint.axvline(x=500, color="blue", linestyle="--")
    _g.ax_joint.annotate(
        "500 UMI",
        xy=(500, _df["pct_counts_mt"].max()),
        xytext=(5, -5),
        textcoords="offset points",
        va="top",
        ha="left",
        color="blue",
        fontweight="bold",
    )

    _g.ax_joint.set_xlabel("Total UMI counts")
    _g.ax_joint.set_ylabel("% mitochondrial counts")
    _g.figure.suptitle("QC-passing barcodes: total counts vs. % mitochondrial", y=1.02)
    _g.figure.text(0.5, 0.96, f"n = {len(_df):,} barcodes", ha="center", fontsize=9, color="dimgray")

    # Dedicated axes for a horizontal colorbar, placed below the joint plot
    # without stealing space from it (keeps marginal alignment intact).
    _joint_pos = _g.ax_joint.get_position()
    _cax = _g.figure.add_axes([_joint_pos.x0, _joint_pos.y0 - 0.1, _joint_pos.width, 0.03])
    _cbar = _g.figure.colorbar(_scatter, cax=_cax, orientation="horizontal")
    _cbar.set_label("% mitochondrial counts")

    plt.show()
    return


@app.cell(hide_code=True)
def qc_round1_summary(adata, mo):
    mo.md(f"""
    ### QC round 1 summary

    Starting from {adata.n_obs:,} raw barcodes, the level-1 filters (`min_counts >= 500`, `max_counts <= 10,000`, `min_genes >= 200`) leave **{int(adata.obs["pass_lvl_1_qc"].sum()):,} barcodes** (`pass_lvl_1_qc`):

    - `pass_min_umi_filter`: {int(adata.obs["pass_min_umi_filter"].sum()):,} pass (most of the loss here is empty droplets/background)
    - `pass_max_umi_filter`: {int(adata.obs["pass_max_umi_filter"].sum()):,} pass
    - `pass_min_gene_filter`: {int(adata.obs["pass_min_gene_filter"].sum()):,} pass

    Among barcodes passing these filters, %mitochondrial content is still elevated (median {adata.obs.loc[adata.obs["pass_lvl_1_qc"], "pct_counts_mt"].median():.1f}%, {adata.obs.loc[adata.obs["pass_lvl_1_qc"], "pct_counts_mt"].gt(15).mean():.1%} above 15%) despite this being nuclei input, which should show close to 0% mito. No mito-based filtering is applied yet at this stage; see the note below.
    """)
    return


@app.cell(hide_code=True)
def notes_on_mito_content(mo):
    mo.md(r"""
    The 10X multi-ome protocol requires nuclei in input, thus we would expect mitochondrial percentage to be close to 0. We see a lot of barcodes with high mitochondrial content. It could be due to incomplete nuclei isolation or ambient/cytoplasmic contamination. For now I am not filtering and check if other metric downstream will clarify if we want to filter them or not.
    """)
    return


@app.cell
def filter_cells(adata):
    adata_flt = adata[adata.obs["pass_lvl_1_qc"]].copy()
    adata_flt.layers["counts"] = adata_flt.X.copy()
    adata_flt
    return (adata_flt,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # RNA downstream analyses
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Normalization

    We use `scclr`'s PFlog (shifted centered log-ratio) normalization instead of the conventional `log1p(CP10K)` approach. PFlog jointly stabilizes technical variance, normalizes for sequencing depth, and preserves within-cell gene ranking by calibrating the log-transform pseudocount from the data's own overdispersion and CLR-centering each cell, rather than using a fixed round-number pseudocount and scale factor. See Booeshaghi et al., ["Normalization for sampled count data"](https://www.biorxiv.org/content/10.1101/2022.05.06.490859) for details.
    """)
    return


@app.cell
def _(adata_flt, scclr):
    scclr.pp.pflog(adata_flt, target="auto")
    adata_flt.uns["pflog"]
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Find highly variable genes
    """)
    return


@app.cell
def _(adata_flt, sc):
    sc.pp.highly_variable_genes(adata_flt, layer="pflog", n_top_genes=2000)
    return


@app.cell
def _(adata_flt, sc):
    sc.pl.highly_variable_genes(adata_flt)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## PCA
    """)
    return


@app.cell
def _(adata_flt, mo, scclr):
    _ncomps = 50
    _ncv = 2*_ncomps + 1
    scclr.tl.pca(adata_flt, n_comps=_ncomps, ncv=_ncv)
    mo.md(f"Computed PCA using {_ncomps} components and {_ncv} Lanczos vectors")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Plot variance ratio
    """)
    return


@app.cell
def _(adata_flt, sc):
    sc.pl.pca_variance_ratio(adata_flt, n_pcs=15, log=False)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Plot PCA
    """)
    return


@app.cell
def pca_axis_setup(adata_flt, pca_axis_title):
    # ------------ Change the PC numbers here--------------------#
    _pc_number_x_axis = 0
    _pc_number_y_axis = 1
    # -----------------------------------------------------------#

    pca_x_title = pca_axis_title(adata_flt, _pc_number_x_axis)
    pca_y_title = pca_axis_title(adata_flt, _pc_number_y_axis)
    return pca_x_title, pca_y_title


@app.cell(hide_code=True)
def pca_pc1_pc2_colored(adata_flt, sc):
    sc.pl.pca(
        adata_flt,
        color=["pct_counts_mt", "pct_counts_pc", "total_counts"],
        dimensions=[(0, 1)],
        ncols=1,
        size=5,
    )
    return


@app.cell(hide_code=True)
def pca_pc3_pc4_colored(adata_flt, sc):
    sc.pl.pca(
        adata_flt,
        color=["pct_counts_mt", "pct_counts_pc", "total_counts"],
        dimensions=[(2, 3)],
        ncols=1,
        size=5,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Clustering
    """)
    return


@app.cell
def compute_knn_neighbors(adata_flt, sc):
    sc.pp.neighbors(adata_flt)
    return


@app.cell
def leiden_clustering(adata_flt, sc):
    # Using the igraph implementation and a fixed number of iterations can be significantly faster,
    # especially for larger datasets
    sc.tl.leiden(adata_flt, flavor="igraph", resolution=0.5, n_iterations=2)
    return


@app.cell
def _(adata_flt, sc):
    sc.pl.pca(
        adata_flt,
        color=["leiden"],
        dimensions=[(0, 1)],
        ncols=1,
        size=5,
    )
    return


@app.cell
def _(adata_flt, sc):
    sc.pl.pca(
        adata_flt,
        color=["leiden"],
        dimensions=[(2, 3)],
        ncols=1,
        size=5,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## UMAP
    """)
    return


@app.cell
def compute_umap(adata_flt, sc):
    sc.tl.umap(adata_flt)
    return


@app.cell
def plot_umap(adata_flt, sc):
    sc.pl.umap(adata_flt, color=["leiden"])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Doublet detection

    We are running doublet detection using `scrublet`. The `expected doublet rate` has been set using the observations made [here](https://www.biostars.org/p/9467937/).
    ```py
    _expected_doublet_rate = adata_flt.n_obs / 1000 * 0.008
    ```
    """)
    return


@app.cell(hide_code=True)
def scrublet_doublet_detection(ad, adata_flt, sc):
    # Run on raw counts (stashed in .layers["counts"] before normalization),
    # regardless of where this cell sits relative to the normalization step.

    # This formula comes from the observations made in 
    # https://www.biostars.org/p/9467937/
    # The specific 10X technical manuals are loaded there
    _expected_doublet_rate = adata_flt.n_obs / 1000 * 0.008

    _raw = ad.AnnData(X=adata_flt.layers["counts"], obs=adata_flt.obs[[]])
    sc.pp.scrublet(_raw, expected_doublet_rate=_expected_doublet_rate)

    adata_flt.obs["doublet_score"] = _raw.obs["doublet_score"].to_numpy()
    adata_flt.obs["predicted_doublet"] = _raw.obs["predicted_doublet"].to_numpy()

    sc.pl.scrublet_score_distribution(_raw)
    return


@app.cell
def doublets_stats(adata_flt):
    adata_flt.obs["predicted_doublet"].value_counts()
    return


@app.cell(hide_code=True)
def doublets_per_leiden_pct(adata_flt, okabe_ito_palette, pd, plt):
    _singlet = (~adata_flt.obs["predicted_doublet"]).rename("singlet")
    _doublet_crosstab = pd.crosstab(adata_flt.obs["leiden"], _singlet)
    _cluster_order = _doublet_crosstab[True].sort_values(ascending=False).index

    _doublet_crosstab_pct = _doublet_crosstab.div(_doublet_crosstab.sum(axis=1), axis=0) * 100
    _doublet_crosstab_pct = _doublet_crosstab_pct.loc[_cluster_order]
    _doublet_crosstab_pct.columns = _doublet_crosstab_pct.columns.astype(str)

    _doublet_crosstab_pct.plot(
        kind="barh", stacked=False,
        color={"True": okabe_ito_palette[0], "False": okabe_ito_palette[6]},
        figsize=(8, 5),
    )
    plt.ylabel("Leiden cluster")
    plt.xlabel("% of barcodes")
    plt.title("Singlets vs. doublets per Leiden cluster")
    plt.legend(title="singlet", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.gcf()
    return


@app.cell
def _(adata_flt, okabe_ito_palette, sc):
    sc.pl.pca(
        adata_flt,
        color=["predicted_doublet"],
        dimensions=[(0, 1)],
        palette={"False": okabe_ito_palette[0], "True": okabe_ito_palette[6]},
        ncols=1,
        size=5,
    )
    return


@app.cell
def qc_doublets_umap(adata_flt, sc):
    sc.pl.umap(
        adata_flt,
        color=["doublet_score"],
        wspace=0.5,
        ncols=2,
    )
    return


@app.cell(hide_code=True)
def cmo_introduction(mo):
    mo.md(r"""
    # CMO hash classification

    Classify each barcode's CMO (Cell Multiplexing Oligo) identity using an approach that mirrors `Seurat::HTODemux`. Briefly, CLR normalize the CMO counts per cell, then call a CMO "positive" if its normalized signal exceeds a per-CMO threshold (95th percentile). Barcodes positive for more than one CMO are called doublets; barcodes positive for none are negatives.

    Classification is restricted to barcodes that already pass the RNA-based QC filtering (`adata_flt`), so the CMO calls reflect real cells rather than empty droplets or debris.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Load CMO counts and perform assignment
    """)
    return


@app.cell
def initialize_cmo(ad, adata_flt, mmread, pd, project_root):
    # path to the cmo counts
    _cmo_counts_path = project_root / "data/cmo_counts/channel1/counts_unfiltered"

    # the counts are in matrix market format
    _cmo_mat = mmread(_cmo_counts_path / "cells_x_features.mtx").tocsr()
    _cmo_barcodes = [
        bc + "_10x_5timepoints_channel1"
        for bc in (_cmo_counts_path / "cells_x_features.barcodes.txt").read_text().splitlines()
    ]
    _cmo_genes = (_cmo_counts_path / "cells_x_features.genes.txt").read_text().splitlines()
    _cmo_gene_names = (_cmo_counts_path / "cells_x_features.genes.names.txt").read_text().splitlines()

    # creating the adata for CMO counts
    adata_cmo = ad.AnnData(
        X=_cmo_mat,
        obs=pd.DataFrame(index=_cmo_barcodes),
        var=pd.DataFrame(
            {"gene_name": _cmo_gene_names},
            index=_cmo_genes,
        ),
    )

    # Filter adata_cmo data to cells in adata_flt
    _cells_in_use = adata_flt.obs_names.intersection(adata_cmo.obs_names.to_list())
    adata_cmo = adata_cmo[_cells_in_use, :]

    # Map each CMO to its timepoint: CMO01-05 = d0, CMO06-10 = d1, ..., CMO21-25 = d4
    cmo_to_timepoint = {f"CMO{i:02d}": f"d{(i - 1) // 5}" for i in range(1, 26)}

    adata_cmo
    return adata_cmo, cmo_to_timepoint


@app.cell
def cmo_clr_normalization(adata_cmo, np):
    _cmo = adata_cmo.X
    _cmo = _cmo.toarray() if hasattr(_cmo, "toarray") else np.asarray(_cmo)

    # 1. Compute Seurat-style geometric mean per cell across the 25 features
    # Seurat sums log1p of non-zero values, divides by total features (25), then takes exp
    _log_counts = np.log1p(np.where(_cmo > 0, _cmo, 0))
    _gm = np.exp(np.sum(_log_counts, axis=1) / 25)

    # 2. Divide raw counts by the geometric mean, then apply log1p
    clr = np.log1p(_cmo / _gm[:, None])
    return (clr,)


@app.cell
def find_cmo_thresholds(clr, np):
    # which quantile to use as cutoff
    positive_quantile_threshold = 0.95
    # get the per tag cutoff
    thresholds = np.quantile(clr, positive_quantile_threshold, axis=0)
    # apply the filter
    positive = clr > thresholds
    n_pos = positive.sum(axis=1)
    return n_pos, positive, thresholds


@app.cell
def assign_cmo_tag(
    adata_cmo,
    adata_flt,
    clr,
    cmo_to_timepoint,
    n_pos,
    np,
    positive,
):
    # Assign each barcode to its strongest-signal CMO (used for singlets and doublets alike),
    # and classify Negative/Singlet/Doublet by how many CMOs cleared their threshold.
    _top_cmo_idx = np.argmax(clr, axis=1)
    _top_cmo_tag = adata_cmo.var["gene_name"].values[_top_cmo_idx]

    adata_flt.obs["cmo_tag_scanpy"] = np.where(n_pos == 0, "Negative", _top_cmo_tag)
    adata_flt.obs["cmo_status_scanpy"] = np.where(
        n_pos == 0, "Negative", np.where(n_pos == 1, "Singlet", "Doublet")
    )
    adata_flt.obs["cmo_positive_tags_scanpy"] = [
        ",".join(adata_cmo.var["gene_name"].values[row]) if row.any() else "Negative"
        for row in positive
    ]

    adata_flt.obs["timepoint_scanpy"] = (
        adata_flt.obs["cmo_positive_tags_scanpy"]
        .map(cmo_to_timepoint)
        .fillna(adata_flt.obs["cmo_status_scanpy"])
    )
    return


@app.cell
def _(adata_flt):
    adata_flt.obs["cmo_status_scanpy"].value_counts()
    return


@app.cell(hide_code=True)
def cmo_threshold_vs_assignment(
    adata_cmo,
    adata_flt,
    cmo_to_timepoint,
    ec_diff_palette,
    pd,
    plt,
    thresholds,
):
    # Per-CMO detection threshold vs. how often that CMO actually wins the argmax
    # assignment. A CMO with a low threshold AND a low assigned count suggests weak/
    # inefficient staining (background and positive signal are both compressed), rather
    # than a contamination or misclassification issue.
    _cmo_summary = pd.DataFrame({
        "threshold": thresholds,
        "n_assigned": adata_flt.obs["cmo_tag_scanpy"].value_counts().reindex(adata_cmo.var["gene_name"]).to_numpy(),
    }, index=adata_cmo.var["gene_name"])
    _cmo_summary["timepoint"] = _cmo_summary.index.map(cmo_to_timepoint)
    _colors = _cmo_summary["timepoint"].map(ec_diff_palette)

    plt.figure(figsize=(7, 6))
    plt.scatter(_cmo_summary["threshold"], _cmo_summary["n_assigned"], c=_colors, s=80, edgecolor="black")
    for _cmo, _row in _cmo_summary.iterrows():
        plt.annotate(_cmo, (_row["threshold"], _row["n_assigned"]), fontsize=7, xytext=(3, 3), textcoords="offset points")
    plt.xlabel("CLR detection threshold (95th percentile)")
    plt.ylabel("Number of barcodes assigned to this CMO")
    plt.title("Per-CMO threshold vs. assignment count")
    plt.tight_layout()
    plt.gcf()
    return


@app.cell(hide_code=True)
def plot_cmo_assignment_counts(
    adata_cmo,
    adata_flt,
    cmo_to_timepoint,
    ec_diff_palette,
    plt,
):
    # Number of barcodes assigned (hash_ID) to each CMO
    _cmo_order = list(adata_cmo.var["gene_name"]) + ["Negative"]
    _counts = adata_flt.obs["cmo_tag_scanpy"].value_counts().reindex(_cmo_order)
    _colors = [ec_diff_palette[cmo_to_timepoint.get(c, "Unassigned")] for c in _counts.index]

    plt.figure(figsize=(6, 8))
    plt.barh(_counts.index, _counts.values, color=_colors)
    plt.ylabel("CMO")
    plt.xlabel("Number of barcodes")
    plt.title("Barcodes assigned per CMO")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.gcf()
    return


@app.cell(hide_code=True)
def plot_timepoint_assignment_counts(
    adata_flt,
    cmo_to_timepoint,
    ec_diff_palette,
    plt,
):
    # Number of barcodes assigned per timepoint (CMOs summed within each timepoint block)
    _timepoint_order = ["d0", "d1", "d2", "d3", "d4", "Negative"]
    _hash_to_timepoint = adata_flt.obs["cmo_tag_scanpy"].map(cmo_to_timepoint).fillna("Negative")
    _timepoint_counts = _hash_to_timepoint.value_counts().reindex(_timepoint_order)
    _colors = [ec_diff_palette.get(t, ec_diff_palette["Unassigned"]) for t in _timepoint_counts.index]

    plt.figure(figsize=(6, 4))
    plt.barh(_timepoint_counts.index, _timepoint_counts.values, color=_colors)
    plt.ylabel("Timepoint")
    plt.xlabel("Number of barcodes")
    plt.title("Barcodes assigned per timepoint")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.gcf()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Comparison with `scrublet` results
    """)
    return


@app.cell(hide_code=True)
def cmo_doublets_per_leiden_pct(adata_flt, okabe_ito_palette, pd, plt):
    # Same cluster ordering (by Scrublet singlet count) applied to both panels for direct comparison
    _singlet_scrublet = (~adata_flt.obs["predicted_doublet"]).rename("singlet")
    _scrublet_crosstab = pd.crosstab(adata_flt.obs["leiden"], _singlet_scrublet)
    _cluster_order = _scrublet_crosstab[True].sort_values(ascending=False).index
    _scrublet_pct = _scrublet_crosstab.div(_scrublet_crosstab.sum(axis=1), axis=0) * 100
    _scrublet_pct = _scrublet_pct.loc[_cluster_order]
    _scrublet_pct.columns = _scrublet_pct.columns.astype(str)

    _singlet_cmo = (adata_flt.obs["cmo_status_scanpy"] == "Singlet").rename("singlet")
    _cmo_crosstab = pd.crosstab(adata_flt.obs["leiden"], _singlet_cmo)
    _cmo_pct = _cmo_crosstab.div(_cmo_crosstab.sum(axis=1), axis=0) * 100
    _cmo_pct = _cmo_pct.loc[_cluster_order]
    _cmo_pct.columns = _cmo_pct.columns.astype(str)

    _fig, (_ax1, _ax2) = plt.subplots(1, 2, figsize=(14, 5), sharey=True, sharex=True)

    _scrublet_pct.plot(
        kind="barh", stacked=False,
        color={"True": okabe_ito_palette[0], "False": okabe_ito_palette[6]},
        ax=_ax1, legend=False,
    )
    _ax1.set_title("Scrublet: singlet vs. doublet")
    _ax1.set_xlabel("% of barcodes")
    _ax1.set_ylabel("Leiden cluster")
    _ax1.set_xlim(0, 100)
    _ax1.invert_yaxis()

    _cmo_pct.plot(
        kind="barh", stacked=False,
        color={"True": okabe_ito_palette[0], "False": okabe_ito_palette[6]},
        ax=_ax2, legend=False,
    )
    _ax2.set_title("CMO hashing: singlet vs. non-singlet")
    _ax2.set_xlabel("% of barcodes")
    _ax2.set_xlim(0, 100)

    _handles = [
        plt.Rectangle((0, 0), 1, 1, color=okabe_ito_palette[0]),
        plt.Rectangle((0, 0), 1, 1, color=okabe_ito_palette[6]),
    ]
    _fig.legend(_handles, ["Singlet", "Not singlet"], loc="upper center", bbox_to_anchor=(0.5, 1.05), ncol=2)
    _fig.suptitle("Clusters ordered by Scrublet singlet count (shared across both panels)", y=1.0, fontsize=9)
    plt.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def scrublet_score_per_cmo(adata_flt, plt, sc):
    # Approximate threshold as the midpoint between the highest non-doublet score
    # and the lowest predicted-doublet score (Scrublet doesn't expose it directly here).
    _scrublet_threshold = (
        adata_flt.obs.loc[adata_flt.obs["predicted_doublet"], "doublet_score"].min()
        + adata_flt.obs.loc[~adata_flt.obs["predicted_doublet"], "doublet_score"].max()
    ) / 2

    _ax = sc.pl.violin(
        adata_flt,
        keys="doublet_score",
        groupby="cmo_status_scanpy",
        order=["Singlet", "Doublet", "Negative"],
        ylabel="Scrublet doublet score",
        palette=["lightgray", "lightgray", "lightgray"],
        show=False,
    )
    _ax.axhline(_scrublet_threshold, color="red", linestyle="--", label=f"Scrublet threshold ({_scrublet_threshold:.2f})")
    _ax.legend()
    plt.gcf()
    return


@app.cell(hide_code=True)
def umap_by_cmo_hashing_altair(adata_flt, alt, mo, okabe_ito_palette, pd):
    _umap_df = pd.DataFrame(
        adata_flt.obsm["X_umap"],
        columns=["UMAP1", "UMAP2"],
        index=adata_flt.obs_names,
    )
    _umap_df["cmo_status_scanpy"] = adata_flt.obs["cmo_status_scanpy"].to_numpy()

    _umap_selection = alt.selection_point(fields=["cmo_status_scanpy"], bind="legend")

    _umap_chart = alt.Chart(_umap_df).mark_circle(size=10).encode(
        x="UMAP1:Q",
        y="UMAP2:Q",
        color=alt.Color(
            "cmo_status_scanpy:N",
            scale=alt.Scale(
                domain=["Singlet", "Doublet", "Negative"],
                range=[okabe_ito_palette[0], okabe_ito_palette[6], okabe_ito_palette[8]],
            ),
        ),
        opacity=alt.condition(_umap_selection, alt.value(0.6), alt.value(0)),
        tooltip=["cmo_status_scanpy"],
    ).properties(title="UMAP colored by CMO hashing classification", width=500, height=500).add_params(_umap_selection)

    mo.ui.altair_chart(_umap_chart, chart_selection=False)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Final round of quality control
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Quality control summary statistics per cluster
    """)
    return


@app.cell(hide_code=True)
def leiden_doublet_overview(adata_flt):
    # Per-cluster overview: does the CMO-hashing doublet rate track the Scrublet doublet rate?
    leiden_doublet_summary = adata_flt.obs.groupby("leiden", observed=True).agg(
        n_cells=("leiden", "size"),
        pct_cmo_singlet=("cmo_status_scanpy", lambda s: (s == "Singlet").mean() * 100),
        pct_cmo_doublet=("cmo_status_scanpy", lambda s: (s == "Doublet").mean() * 100),
        pct_cmo_negative=("cmo_status_scanpy", lambda s: (s == "Negative").mean() * 100),
        pct_scrublet_doublet=("predicted_doublet", lambda s: s.mean() * 100),
        median_mito=("pct_counts_mt", "median"),
        median_counts=("total_counts", "median"),
    ).round(1)
    leiden_doublet_summary.sort_values("pct_cmo_doublet", ascending=False)
    return (leiden_doublet_summary,)


@app.cell(hide_code=True)
def doublet_cluster_investigation_intro(mo):
    mo.md(r"""
    ## Doublet-dominated cluster investigation

    Both methods agree that clusters 2, 7, 8, and 11 are heavily doublet-dominated (CMO hashing 61-77%, Scrublet 76-99%). Cells below dig into cluster 8 (the highest CMO-doublet-rate cluster) as a representative case: (1) whether its CMO "Singlet" calls are threshold-sensitive borderline cases, and (2) whether "Negative" barcodes in the singlet-dominant clusters look transcriptionally like real cells worth rescuing.
    """)
    return


@app.cell(hide_code=True)
def cluster8_threshold_sensitivity(adata_flt, clr, np, pd):
    # Are cluster 8's "Singlet" calls robust, or do they flip to Doublet under a
    # slightly stricter (higher-quantile) threshold? If the CMO thresholds are the
    # problem, singlet counts here should collapse quickly as the quantile rises.
    _cluster8_mask = (adata_flt.obs["leiden"] == "8").to_numpy()
    _clr_8 = clr[_cluster8_mask]

    _rows = []
    for _q in [0.90, 0.95, 0.975, 0.99]:
        _t = np.quantile(clr, _q, axis=0)
        _pos = _clr_8 > _t
        _n_pos = _pos.sum(axis=1)
        _rows.append({
            "quantile": _q,
            "n_singlet": int((_n_pos == 1).sum()),
            "n_doublet": int((_n_pos > 1).sum()),
            "n_negative": int((_n_pos == 0).sum()),
        })

    pd.DataFrame(_rows).set_index("quantile")
    return


@app.cell(hide_code=True)
def cluster8_singlet_margin(adata_flt, clr, np, plt, thresholds):
    # For cluster 8's current "Singlet" calls, how close is the 2nd-strongest CMO to
    # also crossing its own threshold? A small (near-zero) gap means the call is
    # borderline; a large gap means the singlet call is robust regardless of quantile.
    _cluster8_mask = (adata_flt.obs["leiden"] == "8").to_numpy()
    _singlet_8_mask = _cluster8_mask & (adata_flt.obs["cmo_status_scanpy"] == "Singlet").to_numpy()
    _clr_singlet_8 = clr[_singlet_8_mask]
    _margin_to_threshold = -np.sort(-(_clr_singlet_8 - thresholds), axis=1)[:, 1]  # 2nd-highest (value - threshold)

    plt.figure(figsize=(6, 4))
    plt.hist(_margin_to_threshold, bins=40)
    plt.axvline(0, color="red", linestyle="--", label="threshold (0 = right at cutoff)")
    plt.xlabel("2nd-strongest CMO: CLR value minus its threshold")
    plt.ylabel("Number of cluster-8 singlets")
    plt.title("How close is the 2nd CMO to also being called positive?")
    plt.legend()
    plt.tight_layout()
    plt.gcf()
    return


@app.cell(hide_code=True)
def doublet_cluster_conclusion(mo):
    mo.md(r"""
    ### Doublet-dominated clusters: conclusions

    1. **Clusters 2, 7, 8, and 11 are pure doublet clusters, trust Scrublet over CMO here.**
       Both methods agree these four clusters are almost entirely doublets: CMO hashing calls 61-77% doublet, Scrublet independently calls 76-99% doublet.

    2. **The residual CMO "singlets" in these clusters are not robust.**
       Using cluster 8 as the representative case: sweeping the CMO detection quantile from 0.90 to 0.99 makes its singlet count swing wildly (8 to 705 to 1,613 to 830), with no stable quantile. Of the current singlet calls, 25.4% have their 2nd-strongest CMO within 0.1 CLR units of also clearing threshold. As before, this looks like CMO hashing's structural blindness to same-CMO (same-sample) doublets rather than a tunable threshold problem, so we drop these clusters wholesale rather than keep their CMO-labeled singlets.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Cluster 9
    """)
    return


@app.cell(hide_code=True)
def cluster9_mito_filter_impact(adata_flt, mo):
    # How big is cluster 9, and what happens if we retroactively apply a 15% mito cutoff?
    _mask9 = (adata_flt.obs["leiden"] == "9")
    _mito_pass = adata_flt.obs["pct_counts_mt"] <= 15

    _n_total = adata_flt.n_obs
    _n_cluster9 = int(_mask9.sum())
    _n_cluster9_pass_mito = int((_mask9 & _mito_pass).sum())

    _n_other = _n_total - _n_cluster9
    _n_other_pass_mito = int((~_mask9 & _mito_pass).sum())

    mo.md(f"""
    **Cluster 9 size:** {_n_cluster9:,} of {_n_total:,} cells ({_n_cluster9 / _n_total:.1%} of the whole dataset)

    **Effect of a retroactive `pct_counts_mt <= 15` filter:**

    | Group | Before | After | Lost | Lost % |
    |---|---|---|---|---|
    | Cluster 9 | {_n_cluster9:,} | {_n_cluster9_pass_mito:,} | {_n_cluster9 - _n_cluster9_pass_mito:,} | {(_n_cluster9 - _n_cluster9_pass_mito) / _n_cluster9:.1%} |
    | Rest of dataset | {_n_other:,} | {_n_other_pass_mito:,} | {_n_other - _n_other_pass_mito:,} | {(_n_other - _n_other_pass_mito) / _n_other:.1%} |
    | **Total** | {_n_total:,} | {_n_cluster9_pass_mito + _n_other_pass_mito:,} | {_n_total - (_n_cluster9_pass_mito + _n_other_pass_mito):,} | {(_n_total - (_n_cluster9_pass_mito + _n_other_pass_mito)) / _n_total:.1%} |

    Cluster 9 accounts for {(_n_cluster9 - _n_cluster9_pass_mito) / (_n_total - (_n_cluster9_pass_mito + _n_other_pass_mito)):.1%} of all cells that a 15% mito cutoff would remove dataset-wide, despite being only {_n_cluster9 / _n_total:.1%} of the dataset, i.e. a mito filter would disproportionately clear out cluster 9.
    """)
    return


@app.cell(hide_code=True)
def negative_tag_rescue(mo):
    mo.md(r"""
    ## Rescue Negative tags

    **"Negative" barcodes in the singlet-dominant clusters look worth rescuing, and the rescued tags are independently supported by the sub-threshold CMO signal.**

    Restricted to clusters that actually survive `pass_strict_qc` (i.e. excluding both the doublet-dominated clusters and the debris cluster 9), Negatives vs. Singlets: %mito 4.8% vs. 6.2% (Negatives actually lower), Scrublet doublet_score 0.192 vs. 0.195 (essentially identical). Comparable-to-better QC profile supports treating them as real cells with failed CMO staining rather than debris.

       Going further (`rescue_tag_cmo_signal_check`): for each rescued Negative, we checked whether its single closest-to-threshold CMO (the one nearest to, but not over, its cutoff) belongs to the same timepoint as the tag we assigned via cluster consensus. Match rate is 48.9% overall, vs. a 15.9% null baseline (1,000-permutation test) if the rescue tag were random -- roughly 3.1x enrichment (permutation p < 0.001, binomial p ~ 2e-305). This varies by timepoint (d3 92.5%, d0 52.4%, d1/d2/d4 38-46%), but every timepoint sits well above the random baseline. Matched barcodes also have a smaller median gap-to-threshold than mismatched ones, meaning genuine near-miss CMO signal really does track the assigned rescue tag rather than being noise.
    """)
    return


@app.cell(hide_code=True)
def negative_vs_singlet_qc(
    adata_flt,
    doublet_dominated_clusters,
    high_mito_cluster,
    leiden_doublet_summary,
    sc,
):
    # Are "Negative" barcodes that land in largely-Singlet clusters real cells that
    # simply failed CMO hash detection (similar QC profile to their cluster's Singlets),
    # or lower-quality debris (worse QC profile)? Restrict to clusters that actually
    # survive pass_strict_qc (not just "CMO doublet rate < 50%", which would
    # incorrectly include the debris cluster excluded via the separate mito rule).
    _singlet_dominant_clusters = [c for c in leiden_doublet_summary.index if c not in doublet_dominated_clusters and c != high_mito_cluster]
    _compare_mask = adata_flt.obs["leiden"].isin(_singlet_dominant_clusters) & adata_flt.obs["cmo_status_scanpy"].isin(["Singlet", "Negative"])

    sc.pl.violin(
        adata_flt[_compare_mask.to_numpy()],
        keys=["pct_counts_mt", "doublet_score"],
        groupby="cmo_status_scanpy",
        rotation=0,
        multi_panel=True,
    )
    return


@app.cell(hide_code=True)
def rescue_tag_cmo_signal_check(
    adata_cmo,
    adata_flt,
    clr,
    cmo_to_timepoint,
    doublet_dominated_clusters,
    high_mito_cluster,
    leiden_doublet_summary,
    mo,
    np,
    pd,
    thresholds,
):
    # Deeper validation of the rescue heuristic: for each rescued "Negative" barcode
    # (restricted to clusters that actually survive pass_strict_qc -- NOT just
    # "pct_cmo_doublet < 50", which incorrectly includes the debris cluster excluded
    # via the separate mito rule), find its single closest-to-threshold CMO (the
    # smallest gap between CLR value and that CMO's own threshold, even though none
    # cleared it) and check whether that CMO's timepoint matches the tag we assigned
    # via cluster consensus. A high match rate means the sub-threshold CMO signal
    # independently supports the rescue, not just cluster popularity.
    from scipy import stats

    _good_clusters = [c for c in leiden_doublet_summary.index if c not in doublet_dominated_clusters and c != high_mito_cluster]
    _neg_mask = (adata_flt.obs["cmo_status_scanpy"] == "Negative") & adata_flt.obs["leiden"].isin(_good_clusters)

    _clr_neg = clr[_neg_mask.to_numpy()]
    _gap = _clr_neg - thresholds
    _best_cmo_idx = np.argmax(_gap, axis=1)
    _best_cmo = adata_cmo.var["gene_name"].to_numpy()[_best_cmo_idx]
    _best_cmo_timepoint = pd.Series(_best_cmo).map(cmo_to_timepoint).to_numpy()
    _rescued_tag = adata_flt.obs.loc[_neg_mask, "rescued_cmo_tag"].to_numpy()
    _match = _best_cmo_timepoint == _rescued_tag
    _best_gap = _gap[np.arange(len(_gap)), _best_cmo_idx]
    _n = len(_match)
    _k = int(_match.sum())

    # Proper null: 1000-permutation empirical distribution, plus a binomial test
    # against that null rate (rather than eyeballing a handful of shuffles).
    _rng = np.random.default_rng(0)
    _null_rates = np.array([(_best_cmo_timepoint == _rng.permutation(_rescued_tag)).mean() for _ in range(1000)])
    _null_rate = _null_rates.mean()
    _p_permutation = (_null_rates >= _match.mean()).mean()
    _p_binomial = stats.binomtest(_k, _n, _null_rate, alternative="greater").pvalue

    _by_tag = pd.DataFrame({"rescued_tag": _rescued_tag, "match": _match}).groupby("rescued_tag")["match"].agg(["mean", "count"])
    _by_tag = _by_tag.rename(columns={"mean": "pct_match", "count": "n_barcodes"})
    _by_tag["pct_match"] = (_by_tag["pct_match"] * 100).round(1)

    mo.md(f"""
    **Rescue-tag validation:** for {_n:,} rescued Negative barcodes (restricted to
    clusters that actually survive `pass_strict_qc`), the CMO closest to (but not over) its
    threshold matches the assigned rescue tag **{_match.mean():.1%}** of the time
    ({_k:,} of {_n:,}), vs. a **{_null_rate:.1%}** null (1,000-permutation mean, SD
    {_null_rates.std():.1%}) if the rescue tag were random
    (~{_match.mean() / _null_rate:.1f}x enrichment over chance).

    Statistical test: 0 of 1,000 permutations reached the observed rate (permutation
    p < 0.001); a binomial test against the null rate gives p = {_p_binomial:.2e}.
    With n in the thousands this significance is expected -- the 3x-ish enrichment
    is the number that matters for practical significance, not the p-value itself.

    As a complementary check (not relying on the shuffle null): matched barcodes have a
    smaller median gap-to-threshold ({np.median(_best_gap[_match]):.2f}) than mismatched ones
    ({np.median(_best_gap[~_match]):.2f}) -- when there's a genuine near-miss CMO signal, it
    tends to agree with the rescue tag; when the signal is weak/noisy, agreement is closer to chance.

    {_by_tag.to_markdown()}
    """)
    return


@app.cell(hide_code=True)
def rescue_tag_gap_distributions(
    adata_cmo,
    adata_flt,
    clr,
    cmo_to_timepoint,
    doublet_dominated_clusters,
    high_mito_cluster,
    leiden_doublet_summary,
    np,
    okabe_ito_palette,
    pd,
    plt,
    thresholds,
):
    # Distribution of gap-to-threshold for the rescue-tag validation's best CMO,
    # split by whether it matched the assigned rescue tag or not.
    _good_clusters = [c for c in leiden_doublet_summary.index if c not in doublet_dominated_clusters and c != high_mito_cluster]
    _neg_mask = (adata_flt.obs["cmo_status_scanpy"] == "Negative") & adata_flt.obs["leiden"].isin(_good_clusters)

    _clr_neg = clr[_neg_mask.to_numpy()]
    _gap = _clr_neg - thresholds
    _best_cmo_idx = np.argmax(_gap, axis=1)
    _best_cmo = adata_cmo.var["gene_name"].to_numpy()[_best_cmo_idx]
    _best_cmo_timepoint = pd.Series(_best_cmo).map(cmo_to_timepoint).to_numpy()
    _rescued_tag = adata_flt.obs.loc[_neg_mask, "rescued_cmo_tag"].to_numpy()
    _match = _best_cmo_timepoint == _rescued_tag
    _best_gap = _gap[np.arange(len(_gap)), _best_cmo_idx]

    plt.figure(figsize=(7, 5))
    plt.hist(_best_gap[_match], bins=40, alpha=0.6, density=True, color=okabe_ito_palette[5], label=f"Matched (n={_match.sum():,})")
    plt.hist(_best_gap[~_match], bins=40, alpha=0.6, density=True, color=okabe_ito_palette[6], label=f"Mismatched (n={(~_match).sum():,})")
    plt.axvline(0, color="black", linestyle="--", label="threshold (0 = right at cutoff)")
    plt.xlabel("Best CMO: CLR value minus its threshold")
    plt.ylabel("Density")
    plt.title("Gap-to-threshold: matched vs. mismatched rescue tags")
    plt.legend()
    plt.tight_layout()
    plt.gcf()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Final number of barcodes per QC flag
    """)
    return


@app.cell(hide_code=True)
def apply_strict_qc_and_rescue(adata_flt, leiden_doublet_summary, np):
    # --- Consolidated QC decision ---------------------------------------------
    # 1. Doublet-dominated clusters: both CMO hashing AND Scrublet independently
    #    call >50% doublet -> drop the whole cluster.
    # 2. Cluster 9: distinct high-mitochondrial / damaged-nuclei cluster -> drop
    #    the whole cluster regardless of its own per-cell mito value.
    # 3. General mito filter: pct_counts_mt > 15% -> drop, dataset-wide.
    # 4. Scrublet-predicted doublets -> drop.
    # 5. CMO-hashing "Doublet" calls -> drop, even where Scrublet disagrees.
    _doublet_cluster_mask = leiden_doublet_summary["pct_cmo_doublet"].gt(50) & leiden_doublet_summary["pct_scrublet_doublet"].gt(50)
    doublet_dominated_clusters = leiden_doublet_summary.index[_doublet_cluster_mask].tolist()
    high_mito_cluster = "9"

    _is_doublet_cluster = adata_flt.obs["leiden"].isin(doublet_dominated_clusters)
    _is_high_mito_cluster = adata_flt.obs["leiden"] == high_mito_cluster
    _is_high_mito = adata_flt.obs["pct_counts_mt"] > 15
    _is_scrublet_doublet = adata_flt.obs["predicted_doublet"]
    _is_cmo_doublet = adata_flt.obs["cmo_status_scanpy"] == "Doublet"

    adata_flt.obs["qc_exclude_reason"] = np.select(
        [_is_doublet_cluster, _is_high_mito_cluster, _is_scrublet_doublet, _is_cmo_doublet, _is_high_mito],
        ["doublet_cluster", "high_mito_cluster", "scrublet_doublet", "cmo_doublet", "high_mito"],
        default="",
    )
    adata_flt.obs["pass_strict_qc"] = adata_flt.obs["qc_exclude_reason"] == ""

    # --- Rescue barcodes with no CMO tag ("Negative") using their cluster's
    # consensus timepoint (the majority vote among that cluster's own CMO singlets) ---
    cluster_consensus_tag = (
        adata_flt.obs.loc[adata_flt.obs["cmo_status_scanpy"] == "Singlet"]
        .groupby("leiden", observed=True)["timepoint_scanpy"]
        .agg(lambda s: s.mode().iat[0])
    )

    adata_flt.obs["rescued_cmo_tag"] = np.select(
        [adata_flt.obs["cmo_status_scanpy"] == "Singlet", adata_flt.obs["cmo_status_scanpy"] == "Negative"],
        [adata_flt.obs["timepoint_scanpy"], adata_flt.obs["leiden"].map(cluster_consensus_tag)],
        default=None,
    )

    adata_flt.obs["qc_exclude_reason"].value_counts()
    return cluster_consensus_tag, doublet_dominated_clusters, high_mito_cluster


@app.cell
def _(adata_flt):
    adata_flt
    return


@app.cell(hide_code=True)
def pca_by_rescued_timepoint(adata_flt, ec_diff_palette, sc):
    _adata_qc_view = adata_flt[adata_flt.obs["pass_strict_qc"]]
    _palette = {**ec_diff_palette, "None": "#4D4D4D"}
    sc.pl.pca(
        _adata_qc_view,
        color="rescued_cmo_tag",
        palette=_palette,
        na_color="#4D4D4D",
        dimensions=[(0, 1)],
        size=5,
        title="PCA colored by timepoint (pass_strict_qc only)",
    )
    return


@app.cell(hide_code=True)
def umap_by_rescued_timepoint(adata_flt, ec_diff_palette, sc):
    _adata_qc_view = adata_flt[adata_flt.obs["pass_strict_qc"]]
    sc.pl.umap(
        _adata_qc_view,
        color="rescued_cmo_tag",
        palette=ec_diff_palette,
        na_color="#4D4D4D",
        size=5,
        title="UMAP colored by timepoint (pass_strict_qc only)",
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Saving it to file
    """)
    return


@app.cell(hide_code=True)
def write_qc_annotations_tsv(adata_flt, project_root):
    _outfile = project_root / "results/channel1/adata_flt_qc_annotations.tsv"
    _outfile.parent.mkdir(parents=True, exist_ok=True)
    adata_flt.obs.to_csv(_outfile, sep="\t")

    # Show only the path relative to the repo, not the full local filesystem path
    _outfile.relative_to(project_root)
    return


@app.cell(hide_code=True)
def cluster9_conclusion(mo):
    mo.md(r"""
    ## Cluster 9: confirmed debris / damaged nuclei, not a real cell type

    Cluster 9 (1,807 cells, 5.2% of the dataset) is excluded from `pass_strict_qc` based on its mitochondrial fraction: median %mito 23.6% vs. 6.6% for the rest of the dataset, with 88.3% of the cluster exceeding the 15% mito cutoff regardless of CMO status (Doublet 77.5%, Singlet 83.3%, Negative 93.9% individually). See `cluster9_mito_filter_impact` for the dataset-wide impact of a retroactive mito filter.

    **Marker genes** (`cluster9_marker_genes`, Wilcoxon vs. rest, on the `pflog` layer) confirm this directly: the top of the list is dominated by mitochondrially-encoded genes (`MT-CO1/2/3`, `MT-ND1-5`, `MT-ND4L`, `MT-CYB`, `MT-ATP6`) plus a heat-shock/ER-stress signature (`HSP90B1`, `HSP90AB1`, `HSPA5`, `CANX`) and generic ribosomal genes. No cell-type-specific marker or lineage transcription factor appears near the top.

    **Interpretation:** this is the same damaged/dying-nuclei signature we found before the normalization change (previously cluster 7), just under a different arbitrary Leiden cluster number. Excluding cluster 9 wholesale remains well justified.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Cluster 9 marker gene analysis
    """)
    return


@app.cell(hide_code=True)
def cluster9_marker_genes(adata_flt, sc):
    # Cluster 9 marker genes -- if this is debris/damaged nuclei rather than a real
    # cell type, we'd expect either no coherent marker signature, or a signature
    # dominated by stress/mito-adjacent genes rather than specific biology. Uses
    # the full adata_flt since cluster 9 is entirely excluded from pass_strict_qc.
    # Uses the pflog layer (PFlog-normalized) rather than .X, which is raw counts.
    _c9_view = adata_flt.copy()
    sc.tl.rank_genes_groups(_c9_view, groupby="leiden", groups=["9"], reference="rest", method="wilcoxon", layer="pflog", use_raw=False)

    cluster9_markers = sc.get.rank_genes_groups_df(_c9_view, group="9")
    cluster9_markers["gene_symbol"] = cluster9_markers["names"].map(adata_flt.var["gene_symbol"])
    cluster9_markers.head(25)
    return


@app.cell(hide_code=True)
def cluster4_conclusion(mo):
    mo.md(r"""
    ## Cluster 4: a lagging, tip-cell-like transitional population (reproduces under PFlog)

    Cluster 4 (522 cells) is CMO-hashing-labeled mostly d3/d4 (297 d3, 94 d4 among its 401 singlets). Its own-tag-vs-cluster-consensus mismatch rate is 25.9%, in line with the 24.5% "non-majority" (non-d3) fraction of its singlets, and its QC metrics remain unremarkable (total_counts 1,125 vs. 1,450 rest-of-dataset; %mito 8.3% vs. 6.9%; doublet_score 0.28 vs. 0.25, both well under the levels seen in the doublet-dominated clusters).

    **Marker genes** (`cluster4_marker_genes`, Wilcoxon vs. rest, on the `pflog` layer) reproduce the signatures: axon-guidance/cell-motility genes (`KIF26B`, `UNC5C`, `ROBO2`, `SLIT3`, `SEMA6D`, `NRP2`, the last a classic endothelial tip-cell marker), plus adhesion/ECM genes (`PCDH7`, `FN1`, `ALCAM`, `HMCN1`) and progenitor-associated genes (`MLLT3`, `CDK6`).

    **Interpretation:** this cluster reproducing almost identically (composition, QC profile, and marker genes) under a completely different normalization scheme (PFlog vs. the original log1p) is strong evidence this is real biology rather than a normalization artifact: a migratory, tip-cell-like subpopulation that transcriptionally lags the bulk d3/d4 differentiation trajectory. Not excluded; kept as-is.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Tag mismatch analysis
    """)
    return


@app.cell(hide_code=True)
def cluster_tag_mismatch_check(
    adata_flt,
    cluster_consensus_tag,
    doublet_dominated_clusters,
    high_mito_cluster,
):
    # Sanity check on the rescue heuristic: among barcodes that DO have a direct CMO
    # singlet call, how often does that call disagree with what we'd have assigned
    # them via their cluster's consensus (i.e. the same rule used to rescue Negatives)?
    _singlets = adata_flt.obs.loc[adata_flt.obs["cmo_status_scanpy"] == "Singlet"].copy()
    _singlets["cluster_consensus_tag"] = _singlets["leiden"].map(cluster_consensus_tag)
    _mismatch = _singlets["timepoint_scanpy"] != _singlets["cluster_consensus_tag"]

    _n_mismatch = int(_mismatch.sum())
    _n_singlets = len(_singlets)

    _by_cluster = (
        _singlets.assign(mismatch=_mismatch)
        .groupby("leiden", observed=True)["mismatch"]
        .agg(["sum", "count"])
        .rename(columns={"sum": "n_mismatch", "count": "n_singlets"})
    )
    _by_cluster["pct_mismatch"] = (_by_cluster["n_mismatch"] / _by_cluster["n_singlets"] * 100).round(1)
    _by_cluster["marked_for_removal"] = [
        (cl in doublet_dominated_clusters) or (cl == high_mito_cluster) for cl in _by_cluster.index
    ]

    print(f"{_n_mismatch} of {_n_singlets} singlets ({_n_mismatch / _n_singlets:.1%}) have an own-tag that disagrees with their cluster's consensus tag.")
    _by_cluster.sort_values("pct_mismatch", ascending=False)
    return


@app.cell(hide_code=True)
def cluster4_investigation(adata_flt, pd):
    # Cluster 4 flagged earlier for a 24.4% singlet/cluster-consensus mismatch rate,
    # despite sitting below the 50% doublet cutoff for both methods (CMO 24.1%,
    # Scrublet 48.5%). Check whether that mismatch reflects genuine mixed-timepoint
    # biology (a d3/d4 transition state) rather than a technical/doublet problem.
    _mask4 = adata_flt.obs["leiden"] == "4"
    _sub4 = adata_flt.obs.loc[_mask4]

    print(_sub4["cmo_status_scanpy"].value_counts())
    print()
    print("Singlet timepoint composition:")
    print(_sub4.loc[_sub4["cmo_status_scanpy"] == "Singlet", "timepoint_scanpy"].value_counts())
    print()
    print("Cluster 4 median QC vs. rest of dataset:")
    print(pd.DataFrame({
        "cluster_4": _sub4[["total_counts", "n_genes_by_counts", "pct_counts_mt", "doublet_score"]].median(),
        "rest_of_dataset": adata_flt.obs.loc[~_mask4, ["total_counts", "n_genes_by_counts", "pct_counts_mt", "doublet_score"]].median(),
    }))
    return


@app.cell(hide_code=True)
def cluster4_pca_highlight(
    adata_flt,
    alt,
    mo,
    np,
    okabe_ito_palette,
    pca_dataframe,
    pca_scatter,
    pca_x_title,
    pca_y_title,
):
    _qc_view = adata_flt[adata_flt.obs["pass_strict_qc"]]
    _pca_cluster_df = pca_dataframe(_qc_view, 0, 1, ["leiden"])
    _pca_cluster_df["is_cluster4"] = np.where(_pca_cluster_df["leiden"] == "4", "Cluster 4", "Other")

    _cluster4_chart = pca_scatter(
        _pca_cluster_df, pca_x_title, pca_y_title, "is_cluster4", "N",
        "PCA (pass_strict_qc only): cluster 4 highlighted",
        color_scale=alt.Scale(domain=["Cluster 4", "Other"], range=[okabe_ito_palette[6], okabe_ito_palette[8]]),
    )

    mo.ui.altair_chart(_cluster4_chart, legend_selection=["is_cluster4"])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Marker genes analysis
    """)
    return


@app.cell(hide_code=True)
def cluster4_marker_genes(adata_flt, sc):
    # Marker genes for cluster 4 vs. the rest, restricted to QC-passing cells.
    # Uses the pflog layer (PFlog-normalized), since .X is still raw counts --
    # scclr stores the normalized matrix in .layers rather than overwriting .X.
    _qc_view = adata_flt[adata_flt.obs["pass_strict_qc"]].copy()
    sc.tl.rank_genes_groups(_qc_view, groupby="leiden", groups=["4"], reference="rest", method="wilcoxon", layer="pflog", use_raw=False)

    cluster4_markers = sc.get.rank_genes_groups_df(_qc_view, group="4")
    cluster4_markers["gene_symbol"] = cluster4_markers["names"].map(adata_flt.var["gene_symbol"])
    cluster4_markers.head(25)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # ATAC data processing with snapatac2
    """)
    return


@app.cell
def _():
    import snapatac2 as snap

    return (snap,)


@app.cell
def _(pd, project_root):
    chrom_dict = pd.read_csv(project_root / "annotations/GRCh38_EBV.chrom.sizes.no.alt.tsv", sep="\t", header=None, names=["chr", "size"])
    chrom_dict = chrom_dict.set_index("chr")["size"].to_dict()
    list(chrom_dict.items())[1:5]
    return (chrom_dict,)


@app.cell
def _(pd, project_root):
    # Read the list of barcodes passing strict RNA-based QC and use it as the
    # snapatac2 import whitelist -- we trust the RNA-side filtering (doublet
    # clusters, high mito, Scrublet, CMO doublets) and don't re-filter on ATAC
    # metrics; we still compute and plot them for visibility.
    _qc_annotations = pd.read_csv(
        project_root / "results/channel1/adata_flt_qc_annotations.tsv", sep="\t", index_col=0
    )
    assigned_cells = _qc_annotations.index[_qc_annotations["pass_strict_qc"]].to_list()
    len(assigned_cells)
    return (assigned_cells,)


@app.cell
def _(assigned_cells, chrom_dict, igvf_gencode_gtf_path, project_root, snap):
    atac_adata = snap.pp.import_fragments(
        project_root / "data/fragments/10x_5timepoints_channel1.fragments.tsv.gz",
        sorted_by_barcode=False,
        chrom_sizes=chrom_dict,
        whitelist=assigned_cells,
    )
    snap.metrics.tsse(atac_adata, igvf_gencode_gtf_path)
    atac_adata
    return (atac_adata,)


@app.cell
def _(atac_adata):
    atac_adata.obs["tsse"]
    return


@app.cell
def _(atac_adata, mo):
    df = atac_adata.obs[["n_fragment", "tsse"]].copy()

    min_frag = mo.ui.slider(
        start=0,
        stop=int(df["n_fragment"].max()),
        value=1000,
        step=500,
        label="Min fragments"
    )

    min_tsse = mo.ui.slider(
        start=0,
        stop=float(df["tsse"].max()),
        value=5.0,
        step=0.5,
        label="Min TSSE"
    )

    mo.vstack([min_frag, min_tsse])
    return df, min_frag, min_tsse


@app.cell
def _(alt, df, min_frag, min_tsse, mo, okabe_ito_palette):
    # plot unique fragments vs tsse
    alt.renderers.enable("html")

    # QC flag
    plot_df = df.copy()
    plot_df["pass_qc"] = (
        (plot_df["n_fragment"] >= min_frag.value) &
        (plot_df["tsse"] >= min_tsse.value) &
        (plot_df["n_fragment"] <= 15_000)
    ).astype(str)

    # counts
    print(plot_df["pass_qc"].value_counts())

    _pass_qc_scale = alt.Scale(domain=["True", "False"], range=[okabe_ito_palette[1], okabe_ito_palette[5]])

    # base scatter
    base = alt.Chart(plot_df).encode(
        x=alt.X("n_fragment:Q", scale=alt.Scale(type="log"), title="Unique fragments"),
        y=alt.Y("tsse:Q", title="TSSE"),
        color=alt.Color("pass_qc:N", title="Pass QC", scale=_pass_qc_scale)
    )

    scatter = base.mark_circle(size=10, opacity=0.3)

    # marginal distributions
    top_hist = alt.Chart(plot_df).mark_bar(opacity=0.5).encode(
        x=alt.X("n_fragment:Q", scale=alt.Scale(type="log")),
        y=alt.Y("count()"),
        color=alt.Color("pass_qc:N", scale=_pass_qc_scale, legend=None)
    ).properties(height=100)

    right_hist = alt.Chart(plot_df).mark_bar(opacity=0.5).encode(
        y=alt.Y("tsse:Q"),
        x=alt.X("count()"),
        color=alt.Color("pass_qc:N", scale=_pass_qc_scale, legend=None)
    ).properties(width=100)

    # layout
    mo.ui.altair_chart((top_hist & (scatter | right_hist)).resolve_legend(color="shared"))
    return


@app.cell
def _(adata_flt, atac_adata):
    # Use the fuller rescued_cmo_tag (CMO singlet call + cluster-rescued negatives)
    # rather than the raw timepoint_scanpy, which would just read "Negative" for
    # rescued cells even though we're keeping them.
    atac_adata.obs["timepoint_scanpy"] = adata_flt.obs.loc[atac_adata.obs_names, "rescued_cmo_tag"].astype("category")
    return


@app.cell
def _(atac_adata):
    atac_adata.obs["timepoint_scanpy"].value_counts()
    return


@app.cell
def _(mo):
    # Intentionally NOT calling snap.pp.filter_cells here: the barcode set was
    # already restricted to pass_strict_qc cells at import time (see the whitelist
    # above), so we trust the RNA-based QC rather than layering an additional
    # ATAC-metric-based filter on top. The QC plots below are for visibility only.
    mo.md("""
    Skipping ATAC-based cell filtering -- relying on the RNA `pass_strict_qc` whitelist used at import.
    """)
    return


@app.cell
def _(atac_adata, snap):
    snap.pp.add_tile_matrix(atac_adata)
    return


@app.cell
def _(atac_adata, snap):
    snap.pp.select_features(atac_adata, n_features=250000)
    return


@app.cell
def _(atac_adata, snap):
    snap.tl.spectral(atac_adata)
    snap.tl.umap(atac_adata)
    snap.pp.knn(atac_adata)
    snap.tl.leiden(atac_adata)
    return


@app.cell(hide_code=True)
def atac_umap_helpers(alt, pd):
    def atac_umap_dataframe(adata, color_cols):
        _coords = adata.obsm["X_umap"]
        return pd.DataFrame(
            {"UMAP1": _coords[:, 0], "UMAP2": _coords[:, 1], **{c: adata.obs[c].to_numpy() for c in color_cols}},
            index=adata.obs_names,
        )


    def atac_umap_scatter(df, color, color_type, title, color_scale=None):
        _color = alt.Color(f"{color}:{color_type}", title=color)
        if color_scale is not None:
            _color = _color.scale(color_scale)

        return alt.Chart(df).mark_circle(size=10, opacity=0.6).encode(
            x=alt.X("UMAP1:Q"),
            y=alt.Y("UMAP2:Q"),
            color=_color,
            tooltip=[color],
        ).properties(title=title, width=500, height=500)

    return atac_umap_dataframe, atac_umap_scatter


@app.cell
def _(atac_adata, atac_umap_dataframe, atac_umap_scatter, mo):
    _leiden_df = atac_umap_dataframe(atac_adata, ["leiden"])
    _leiden_chart = atac_umap_scatter(_leiden_df, "leiden", "N", "ATAC UMAP colored by Leiden cluster")

    mo.ui.altair_chart(_leiden_chart, legend_selection=["leiden"])
    return


@app.cell
def _(atac_adata):
    print(atac_adata.obs["timepoint_scanpy"].cat.categories.tolist())
    return


@app.cell
def _(atac_adata, pd):
    pd.crosstab(atac_adata.obs["leiden"], atac_adata.obs["timepoint_scanpy"])
    return


@app.cell
def _(
    alt,
    atac_adata,
    atac_umap_dataframe,
    atac_umap_scatter,
    ec_diff_palette,
    mo,
):
    _timepoint_df = atac_umap_dataframe(atac_adata, ["timepoint_scanpy"])
    _present_timepoints = [t for t in ec_diff_palette if t in _timepoint_df["timepoint_scanpy"].unique()]
    _timepoint_scale = alt.Scale(domain=_present_timepoints, range=[ec_diff_palette[t] for t in _present_timepoints])

    _timepoint_chart = atac_umap_scatter(
        _timepoint_df, "timepoint_scanpy", "N", "ATAC UMAP colored by timepoint", color_scale=_timepoint_scale
    )

    mo.ui.altair_chart(_timepoint_chart, legend_selection=["timepoint_scanpy"])
    return


@app.cell(hide_code=True)
def rna_cluster4_on_atac_umap(
    adata_flt,
    alt,
    atac_adata,
    atac_umap_dataframe,
    atac_umap_scatter,
    mo,
    np,
    okabe_ito_palette,
):
    _atac_cluster4_df = atac_umap_dataframe(atac_adata, [])
    _atac_cluster4_df["is_rna_cluster4"] = np.where(
        adata_flt.obs.loc[atac_adata.obs_names, "leiden"].to_numpy() == "4", "RNA cluster 4", "Other"
    )

    _atac_cluster4_chart = atac_umap_scatter(
        _atac_cluster4_df, "is_rna_cluster4", "N", "ATAC UMAP: RNA cluster 4 barcodes highlighted",
        color_scale=alt.Scale(domain=["RNA cluster 4", "Other"], range=[okabe_ito_palette[6], okabe_ito_palette[8]]),
    )

    mo.ui.altair_chart(_atac_cluster4_chart, legend_selection=["is_rna_cluster4"])
    return


@app.cell
def _(atac_adata):
    atac_adata.obs["timepoint_scanpy"].value_counts()
    return


@app.cell
def imports():
    import altair as alt
    import anndata as ad
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    import scanpy as sc
    import seaborn as sns
    from dotenv import find_dotenv, load_dotenv
    from pathlib import Path
    from scipy.io import mmread

    # get the project root path using the '.env' file.
    _env_path = find_dotenv(usecwd=True)
    project_root = Path(_env_path).parent


    return ad, alt, mmread, np, pd, plt, project_root, sc, sns


@app.cell
def _():
    import scclr

    # PFlog (shifted centered log-ratio) normalization instead of log1p(CP10K):
    # https://www.biorxiv.org/content/10.1101/2022.05.06.490859
    # jointly stabilizes technical variance, normalizes for sequencing depth, and
    # preserves within-cell gene ranking (monotonicity) via a data-calibrated
    # pseudocount and CLR centering, rather than a fixed round-number pseudocount.
    return (scclr,)


@app.cell
def _():
    # color paletter for the endothelial differentiation
    ec_diff_palette = {
      "d0": "#C6C7C7",
      "d1": "#A8B1D6",
      "d2": "#EBBC9E",
      "d3": "#FBC1C3",
      "d4": "#F7999C",
      "Unassigned": "#4D4D4D",
    }

    # Okabe-Ito colorblind-safe 9-color palette, used for all non-timepoint categorical plots
    okabe_ito_palette = [
      "#000000",  # black
      "#E69F00",  # orange
      "#56B4E9",  # sky blue
      "#009E73",  # bluish green
      "#F0E442",  # yellow
      "#0072B2",  # blue
      "#D55E00",  # vermillion
      "#CC79A7",  # reddish purple
      "#999999",  # gray
    ]
    return ec_diff_palette, okabe_ito_palette


@app.cell(hide_code=True)
def pca_helpers(alt, pd):
    def pca_axis_title(adata, pc_idx):
        _pct = adata.uns["pca"]["variance_ratio"][pc_idx] * 100
        return f"PC{pc_idx + 1} ({_pct:.1f}% var.)"


    def pca_dataframe(adata, x_pc, y_pc, color_cols):
        _coords = adata.obsm["X_pca"][:, [x_pc, y_pc]]
        return pd.DataFrame(
            {"x": _coords[:, 0], "y": _coords[:, 1], **{c: adata.obs[c].to_numpy() for c in color_cols}}
        )


    def pca_scatter(df, x_title, y_title, color, color_type, title, color_scale=None, hide_on_deselect=False, x_domain=None, y_domain=None):
        _color = alt.Color(
            f"{color}:{color_type}",
            title=color,
            legend=alt.Legend(labelFontSize=13, titleFontSize=14, symbolSize=100),
        )
        if color_scale is not None:
            _color = _color.scale(color_scale)

        _x_scale = alt.Scale(domain=x_domain) if x_domain is not None else alt.Undefined
        _y_scale = alt.Scale(domain=y_domain) if y_domain is not None else alt.Undefined

        _chart = alt.Chart(df).mark_circle(size=20, opacity=1).encode(
            x=alt.X("x:Q", title=x_title, scale=_x_scale),
            y=alt.Y("y:Q", title=y_title, scale=_y_scale),
            color=_color,
            tooltip=[color],
        ).properties(title=title, width=350, height=350)

        if hide_on_deselect:
            _selection = alt.selection_point(fields=[color], bind="legend")
            _chart = _chart.encode(
                opacity=alt.condition(_selection, alt.value(1), alt.value(0))
            ).add_params(_selection)

        return _chart

    return pca_axis_title, pca_dataframe, pca_scatter


if __name__ == "__main__":
    app.run()
