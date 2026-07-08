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


    _env_path = find_dotenv(usecwd=True)
    project_root = Path(_env_path).parent

    # color paletter for the endothelial differentiation
    ec_diff_palette = {
      "d0": "#C6C7C7",
      "d1": "#A8B1D6",
      "d2": "#EBBC9E",
      "d3": "#FBC1C3",
      "d4": "#F7999C",
      "Unassigned": "#4D4D4D",
    }
    return ad, alt, ec_diff_palette, mmread, np, pd, plt, project_root, sc, sns


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Introduction

    This notebook processes the 10x multiome data for the endothelial differentiation time course: 5 timepoints (d0-d4), 5 biological replicates each. Samples were multiplexed using the MULTI-seq technique with CMO (Cell Multiplexing Oligo) barcodes. Data was processed with the IGVF pipeline using kallisto and the GENCODE v43 annotation, and CMO quantification was performed with the `kite` workflow from the `kallisto-bustools` suite.

    The goals of this notebook are to:

    1. Perform quality control filtering on the RNA data
    2. Run CMO hash classification (mimicking `Seurat::HTODemux`) on QC-passing cells
    3. Assign each cell barcode to a CMO (and therefore to a timepoint/replicate)
    4. Process the corresponding ATAC data with `snapatac2`
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Loading old CMO assignments for QC purposes
    """)
    return


@app.cell
def load_old_cmo(pd, project_root):
    # relative path to the previously annotated CMO assignments
    _old_cmo_assignments_fnp = "data/10x_multi_5_timepoints_original_barcode_cmo_assignments_from_dulguun.dedup.tsv"
    # load the assignments
    _old_cmo_assignments = pd.read_csv(
        project_root / _old_cmo_assignments_fnp,
        sep="\t",
        header=None,
        names=["barcode", "annotation"]
    )

    # keep only channel 1 explicitly
    _old_cmo_assignments_ch1 = _old_cmo_assignments[
        _old_cmo_assignments["barcode"].str.endswith("-1_1")
    ].copy()

    # build cleaned barcode
    _old_cmo_assignments_ch1["cleaned_barcode"] = (
        _old_cmo_assignments_ch1["barcode"]
        .str.replace(r"-1_1$", "", regex=True)
        + "_10x_5timepoints_channel1"
    )

    _old_cmo_assignments_ch1.reset_index(drop=True, inplace=True)
    # create final dataframe for downstream analysis
    old_assignments_df = (
        _old_cmo_assignments_ch1
        .loc[:, ["cleaned_barcode", "annotation"]]
        .rename(columns={
            "cleaned_barcode": "cell_barcode",
            "annotation": "cmo_old_annotation"
        })
        .set_index("cell_barcode")
    )
    old_assignments_df["cmo_old_annotation"].value_counts()
    return (old_assignments_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Loading the TSS annotations
    The transcription start sites(TSSs) for the protein-coding genes annotated in GENCODE v43 were curate using the MANE annotations. The annotation file includes the -250/+249 region around the TSS for a total of 500bp.
    """)
    return


@app.cell
def load_gene_metadata(pd, project_root):
    _gene_metadata_fnp = "annotations/gencode.v43.protein_coding.TSS500bp.bed"
    gene_metadata_df = pd.read_csv(
        project_root / _gene_metadata_fnp,
        sep="\t",
        header=0,
        names=["chr", "start", "end", "gene_symbol", "score", "strand", "gene_id", "gene_type"]
    )
    gene_metadata_df
    return (gene_metadata_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Loading the counts in scanpy

    We are going to load the total counts as produced by `kallisto-bustool` and annotate the genes with the metadata loaded the cell above. We are going to add to the barcodes the old CMO assignments for QC purposes.
    """)
    return


@app.cell
def load_h5_counts(
    ad,
    ec_diff_palette,
    gene_metadata_df,
    old_assignments_df,
    project_root,
    sc,
):
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

    adata.obs = adata.obs.merge(
        old_assignments_df,
        left_index=True,
        right_index=True,
        how="left"
    )
    adata.obs["cmo_old_annotation"] = adata.obs["cmo_old_annotation"].fillna("Unassigned")

    # Transform to categories
    adata.obs["cmo_old_annotation"] = adata.obs["cmo_old_annotation"].astype("category")

    _categories = adata.obs["cmo_old_annotation"].cat.categories
    adata.uns["cmo_old_annotation_colors"] = [
        ec_diff_palette[c] for c in _categories
    ]


    # -------------------------
    # QC gene flags
    # -------------------------
    _gene_symbol = adata.var["gene_symbol"].fillna("")

    adata.var["mt"] = _gene_symbol.str.startswith("MT-")
    adata.var["ribo"] = _gene_symbol.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = _gene_symbol.str.contains(r"^HB(?!P)", regex=True)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

    adata
    return (adata,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### QC cDNA library
    """)
    return


@app.cell
def basic_qc_violin_plots(adata, sc):
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
        jitter=0.4,
        multi_panel=False,
    )
    return


@app.cell
def plot_knee_plot(adata, np, plt):
    # Plot knee plot for n_genes_by_counts
    # total counts per cell
    counts = adata.obs.loc[:, "total_counts"].values

    # sort descending
    counts_sorted = np.sort(counts)[::-1]

    # rank
    ranks = np.arange(1, len(counts_sorted) + 1)

    # plot
    plt.figure(figsize=(5,5))
    plt.loglog(ranks, counts_sorted)
    plt.xlabel("Cell rank")
    plt.ylabel("Total counts")
    plt.title("Knee plot")
    plt.show()
    return


@app.cell
def old_cmo_qc_setup(adata, np):
    adata.obs["old_cmo_group"] = np.where(
        adata.obs["cmo_old_annotation"] != "Unassigned",
        "Has old CMO assignment",
        "Unassigned (rest)",
    )
    return


@app.cell(hide_code=True)
def old_cmo_qc_total_counts(adata, sc):
    sc.pl.violin(
        adata,
        "total_counts",
        groupby="old_cmo_group",
        jitter=0.4,
        rotation=45,
        ylabel="Total UMI counts",
    )
    return


@app.cell(hide_code=True)
def old_cmo_qc_n_genes_by_counts(adata, sc):
    sc.pl.violin(
        adata,
        "n_genes_by_counts",
        groupby="old_cmo_group",
        jitter=0.4,
        rotation=45,
        ylabel="Number of genes detected",
    )
    return


@app.cell(hide_code=True)
def old_cmo_qc_pct_counts_mt(adata, sc):
    sc.pl.violin(
        adata,
        "pct_counts_mt",
        groupby="old_cmo_group",
        jitter=0.4,
        rotation=45,
        ylabel="% mitochondrial counts",
    )
    return


@app.cell(hide_code=True)
def old_cmo_qc_pct_counts_ribo(adata, sc):
    sc.pl.violin(
        adata,
        "pct_counts_ribo",
        groupby="old_cmo_group",
        jitter=0.4,
        rotation=45,
        ylabel="% ribosomal counts",
    )
    return


@app.cell(hide_code=True)
def old_cmo_qc_pct_counts_hb(adata, sc):
    sc.pl.violin(
        adata,
        "pct_counts_hb",
        groupby="old_cmo_group",
        jitter=0.4,
        rotation=45,
        ylabel="% hemoglobin counts",
    )
    return


@app.cell
def extra_qc_scatter(adata, sc):
    sc.pl.scatter(adata, "total_counts", "pct_counts_mt", color="old_cmo_group")
    sc.pl.scatter(adata, "total_counts", "pct_counts_ribo", color="old_cmo_group")
    sc.pl.scatter(adata, "total_counts", "pct_counts_hb", color="old_cmo_group")
    return


@app.cell(hide_code=True)
def old_cmo_total_counts_dist(adata, plt):
    adata.obs.loc[adata.obs["cmo_old_annotation"] != "Unassigned", "total_counts"].plot.hist(
        bins=100,
        title="Total UMI counts for barcodes with an old CMO assignment",
    )
    plt.xlabel("Total UMI counts")
    plt.gcf()
    return


@app.cell(hide_code=True)
def old_cmo_n_genes_dist(adata, plt):
    adata.obs.loc[adata.obs["cmo_old_annotation"] != "Unassigned", "n_genes_by_counts"].plot.hist(
        bins=100,
        title="Number of genes detected for barcodes with an old CMO assignment",
    )
    plt.xlabel("Number of genes detected")
    plt.gcf()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Filter cells
    """)
    return


@app.cell
def mask_filter_cells(adata, mo, pd):
    # QC pass mask mirroring the cell-level filters applied below
    #_max_counts_cutoff = adata.obs["total_counts"].quantile(0.999)
    _max_counts_cutoff = 10_000

    _has_old_assignment = adata.obs["cmo_old_annotation"] != "Unassigned"

    _steps = [
        ("min_counts >= 500", adata.obs["total_counts"] >= 500),
        (f"max_counts <= {_max_counts_cutoff}", adata.obs["total_counts"] <= _max_counts_cutoff),
        ("min_genes >= 200", adata.obs["n_genes_by_counts"] >= 200),
    ]

    _remaining_mask = pd.Series(True, index=adata.obs.index)
    _rows = []
    for _label, _step_mask in _steps:
        _before_all = int(_remaining_mask.sum())
        _before_old = int((_remaining_mask & _has_old_assignment).sum())
        _remaining_mask &= _step_mask
        _after_all = int(_remaining_mask.sum())
        _after_old = int((_remaining_mask & _has_old_assignment).sum())
        _lost_all = _before_all - _after_all
        _lost_old = _before_old - _after_old
        _rows.append(
            f"| `{_label}` | {_before_all:,} → {_after_all:,} "
            f"(lost {_lost_all:,}, {_lost_all / _before_all:.1%}) | "
            f"{_before_old:,} → {_after_old:,} "
            f"(lost {_lost_old:,}, {_lost_old / _before_old:.1%}) |"
        )

    adata.obs["pass_qc"] = _remaining_mask

    _n_old_assignment_surviving = int((_has_old_assignment & adata.obs["pass_qc"]).sum())
    _n_old_assignment_total = int(_has_old_assignment.sum())

    mo.md(f"""
    **QC filter breakdown** (min_counts=500, max_counts={_max_counts_cutoff:,}, min_genes=200)

    | Filter step | All barcodes | Old-assignment barcodes |
    |---|---|---|
    {chr(10).join(_rows)}

    **Net:** {_n_old_assignment_surviving:,} of {_n_old_assignment_total:,} barcodes with an old CMO assignment
    ({_n_old_assignment_surviving / _n_old_assignment_total:.1%}) would survive this QC filter.
    """)
    return


@app.cell
def checking_mito_content(adata, np, plt, sns):
    _df = adata.obs.loc[adata.obs["pass_qc"], ["total_counts", "pct_counts_mt", "cmo_old_annotation"]].copy()
    _df["has_assignment"] = np.where(_df["cmo_old_annotation"] != "Unassigned", "Has old CMO assignment", "Unassigned (rest)")

    _g = sns.JointGrid(data=_df, x="total_counts", y="pct_counts_mt", height=6)
    _scatter = _g.ax_joint.scatter(
        _df["total_counts"], _df["pct_counts_mt"],
        c=_df["pct_counts_mt"], cmap="viridis", s=5, alpha=0.5,
    )

    _groups = {
        "Unassigned (rest)": "gray",
        "Has old CMO assignment": "orange",
    }
    for _label, _color in _groups.items():
        _sub = _df.loc[_df["has_assignment"] == _label]
        _g.ax_marg_x.hist(_sub["total_counts"], bins=100, color=_color, alpha=0.5, density=True, label=_label)
        _g.ax_marg_y.hist(_sub["pct_counts_mt"], bins=100, orientation="horizontal", color=_color, alpha=0.5, density=True, label=_label)

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

    _handles, _labels = _g.ax_marg_x.get_legend_handles_labels()
    _g.figure.legend(_handles, _labels, loc="lower center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize=8)

    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The 10X multi-ome protocol requires nuclei in input, thus we would expect mitochondrial percentage to be close to 0. We see a lot of barcodes with high mitochondrial content. It could be due to incomplete nuclei isolation or ambient/cytoplasmic contamination. For now I am not filtering and check if other metric downstream will clarify if we want to filter them or not.
    """)
    return


@app.cell
def filter_cells(adata):
    adata_flt = adata[adata.obs["pass_qc"]].copy()
    adata_flt.layers["counts"] = adata_flt.X.copy()
    adata_flt
    return (adata_flt,)


@app.cell
def _(adata_flt, sc):
    # Normalizing to median total counts
    sc.pp.normalize_total(adata_flt)
    # Logarithmize the data
    sc.pp.log1p(adata_flt)
    return


@app.cell
def _(adata_flt, sc):
    sc.pp.highly_variable_genes(adata_flt, n_top_genes=2000)
    return


@app.cell
def _(adata_flt, sc):
    sc.pl.highly_variable_genes(adata_flt)
    return


@app.cell
def _(adata_flt, sc):
    sc.tl.pca(adata_flt)
    return


@app.cell
def _(adata_flt, sc):
    sc.pl.pca_variance_ratio(adata_flt, n_pcs=15, log=True)
    return


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


@app.cell
def qc_pca_before_new_assignment(
    adata_flt,
    alt,
    ec_diff_palette,
    mo,
    pca_axis_title,
    pca_dataframe,
    pca_scatter,
):
    # ------------ Change the PC numbers here--------------------#
    _pc_number_x_axis = 0
    _pc_number_y_axis = 1
    # -----------------------------------------------------------#

    pca_x_title = pca_axis_title(adata_flt, _pc_number_x_axis)
    pca_y_title = pca_axis_title(adata_flt, _pc_number_y_axis)
    pca_df = pca_dataframe(adata_flt, _pc_number_x_axis, _pc_number_y_axis, ["cmo_old_annotation", "pct_counts_mt"])

    pca_x_domain = [pca_df["x"].min(), pca_df["x"].max()]
    pca_y_domain = [pca_df["y"].min(), pca_df["y"].max()]

    _counts = pca_df["cmo_old_annotation"].value_counts()
    _label_map = {k: f"{k} ({_counts[k]:,})" for k in _counts.index}
    pca_df["cmo_old_annotation_label"] = pca_df["cmo_old_annotation"].map(_label_map)

    _cmo_color_scale = alt.Scale(
        domain=[_label_map[k] for k in ec_diff_palette if k in _label_map],
        range=[v for k, v in ec_diff_palette.items() if k in _label_map],
    )

    _annotation_chart = pca_scatter(
        pca_df, pca_x_title, pca_y_title, "cmo_old_annotation_label", "N",
        "PCA colored by old CMO assignment (click legend to filter)",
        color_scale=_cmo_color_scale,
        hide_on_deselect=True,
        x_domain=pca_x_domain,
        y_domain=pca_y_domain,
    )

    pca_old_label_selection = mo.ui.altair_chart(_annotation_chart, chart_selection=False)
    pca_old_label_selection
    return (
        pca_df,
        pca_old_label_selection,
        pca_x_domain,
        pca_x_title,
        pca_y_domain,
        pca_y_title,
    )


@app.cell(hide_code=True)
def qc_pca_mito_before_new_assignment(
    mo,
    pca_df,
    pca_old_label_selection,
    pca_scatter,
    pca_x_domain,
    pca_x_title,
    pca_y_domain,
    pca_y_title,
):
    _selected_df = pca_old_label_selection.value if len(pca_old_label_selection.value) > 0 else pca_df

    _mito_chart = pca_scatter(
        _selected_df, pca_x_title, pca_y_title, "pct_counts_mt", "Q",
        "PCA colored by % mitochondrial counts",
        x_domain=pca_x_domain,
        y_domain=pca_y_domain,
    )

    mo.ui.altair_chart(_mito_chart, chart_selection=False, legend_selection=False)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Leiden clustering
    """)
    return


@app.cell
def compute_knn_neighbors(adata_flt, sc):
    sc.pp.neighbors(adata_flt)
    return


@app.cell
def compute_umap(adata_flt, sc):
    sc.tl.umap(adata_flt)
    return


@app.cell
def leiden_clustering(adata_flt, sc):
    # Using the igraph implementation and a fixed number of iterations can be significantly faster,
    # especially for larger datasets
    sc.tl.leiden(adata_flt, flavor="igraph", n_iterations=2)
    return


@app.cell
def plot_umap(adata_flt, sc):
    sc.pl.umap(adata_flt, color=["leiden"])
    return


@app.cell(hide_code=True)
def leiden_cmo_barplot(adata_flt, ec_diff_palette, pd, plt):
    _crosstab = pd.crosstab(adata_flt.obs["leiden"], adata_flt.obs["cmo_old_annotation"])
    _crosstab = _crosstab[[c for c in ec_diff_palette if c in _crosstab.columns and c != "Unassigned"]]

    # Order clusters: first the cluster with most d0, then (of those left) most d1, etc.
    _remaining_clusters = list(_crosstab.index)
    _cluster_order = []
    for _tp in _crosstab.columns:
        _best_cluster = _crosstab.loc[_remaining_clusters, _tp].idxmax()
        _cluster_order.append(_best_cluster)
        _remaining_clusters.remove(_best_cluster)
    _cluster_order.extend(_remaining_clusters)

    _crosstab = _crosstab.loc[_cluster_order]

    _crosstab.plot(kind="bar", stacked=False, color=ec_diff_palette, figsize=(8, 5))
    plt.xlabel("Leiden cluster")
    plt.ylabel("Number of barcodes")
    plt.title("CMO assignment composition per Leiden cluster")
    plt.legend(title="cmo_old_annotation", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.gcf()
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


@app.cell
def qc_doublets_umap(adata_flt, sc):
    sc.pl.umap(
        adata_flt[~adata_flt.obs["predicted_doublet"]],
        color=["doublet_score"],
        wspace=0.5,
        ncols=2,
    )
    return


@app.cell(hide_code=True)
def doublets_per_leiden_pct(adata_flt, pd, plt):
    _singlet = (~adata_flt.obs["predicted_doublet"]).rename("singlet")
    _doublet_crosstab = pd.crosstab(adata_flt.obs["leiden"], _singlet)
    _cluster_order = _doublet_crosstab[True].sort_values(ascending=False).index

    _doublet_crosstab_pct = _doublet_crosstab.div(_doublet_crosstab.sum(axis=1), axis=0) * 100
    _doublet_crosstab_pct = _doublet_crosstab_pct.loc[_cluster_order]
    _doublet_crosstab_pct.columns = _doublet_crosstab_pct.columns.astype(str)

    _doublet_crosstab_pct.plot(kind="barh", stacked=False, color={"True": "black", "False": "red"}, figsize=(8, 5))
    plt.ylabel("Leiden cluster")
    plt.xlabel("% of barcodes")
    plt.title("Singlets vs. doublets per Leiden cluster")
    plt.legend(title="singlet", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.gcf()
    return


@app.cell(hide_code=True)
def cmo_introduction(mo):
    mo.md(r"""
    ## CMO hash classification

    Classify each barcode's CMO (Cell Multiplexing Oligo) identity using an approach that mirrors `Seurat::HTODemux`: CLR: normalize the CMO counts per cell, then call a CMO "positive" if its normalized signal exceeds a per-CMO threshold (95th percentile). Barcodes positive for more than one CMO are called doublets; barcodes positive for none are negatives.

    Classification is restricted to barcodes that already pass the RNA-based QC filtering (`adata_flt`), so the CMO calls reflect real cells rather than empty droplets or debris.
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
    return n_pos, positive, positive_quantile_threshold, thresholds


@app.cell
def plot_cmo_thresholds(
    adata_cmo,
    plt,
    positive_quantile_threshold,
    thresholds,
):
    # Bar plot of the per-CMO detection threshold (x threshold percentile of CLR signal)
    plt.figure(figsize=(6, 8))
    plt.barh(adata_cmo.var["gene_name"], thresholds)
    plt.ylabel("CMO")
    plt.xlabel("CLR threshold")
    plt.title(f"CLR thresholds at {positive_quantile_threshold} quantile")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.gcf()
    return


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
def cmo_doublets_per_leiden_pct(adata_flt, pd, plt):
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

    _scrublet_pct.plot(kind="barh", stacked=False, color={"True": "black", "False": "red"}, ax=_ax1, legend=False)
    _ax1.set_title("Scrublet: singlet vs. doublet")
    _ax1.set_xlabel("% of barcodes")
    _ax1.set_ylabel("Leiden cluster")
    _ax1.set_xlim(0, 100)
    _ax1.invert_yaxis()

    _cmo_pct.plot(kind="barh", stacked=False, color={"True": "black", "False": "red"}, ax=_ax2, legend=False)
    _ax2.set_title("CMO hashing: singlet vs. non-singlet")
    _ax2.set_xlabel("% of barcodes")
    _ax2.set_xlim(0, 100)

    _handles = [plt.Rectangle((0, 0), 1, 1, color="black"), plt.Rectangle((0, 0), 1, 1, color="red")]
    _fig.legend(_handles, ["Singlet", "Not singlet"], loc="upper center", bbox_to_anchor=(0.5, 1.05), ncol=2)
    _fig.suptitle("Clusters ordered by Scrublet singlet count (shared across both panels)", y=1.0, fontsize=9)
    plt.tight_layout()
    _fig
    return


@app.cell
def _(adata_flt):
    adata_flt.obs["cmo_status_scanpy"].value_counts()
    return


@app.cell(hide_code=True)
def umap_by_cmo_hashing_altair(adata_flt, alt, mo, pd):
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
            scale=alt.Scale(domain=["Singlet", "Doublet", "Negative"], range=["black", "red", "lightgray"]),
        ),
        opacity=alt.condition(_umap_selection, alt.value(0.6), alt.value(0)),
        tooltip=["cmo_status_scanpy"],
    ).properties(title="UMAP colored by CMO hashing classification", width=500, height=500).add_params(_umap_selection)

    mo.ui.altair_chart(_umap_chart, chart_selection=False)
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
def _(mo):
    mo.md(r"""
    ### Cluster 11 conclusions
    1. Cluster 11 is essentially a pure doublet cluster — trust Scrublet over CMO here.
    The per-cluster overview shows five clusters (1, 3, 10, 11, 12) where CMO hashing calls 67–78% doublet and Scrublet independently calls 97.5–99.7% doublet — two completely different methods (oligo counting vs. transcriptome complexity) converging on "almost everything here is a doublet." That agreement is strong evidence this isn't a fluke of one method.

    2. The residual 20% "singlets" in cluster 11 are not robust — this does look like a threshold problem, but a structural one, not a "pick a better number" one.
    - Sweeping the quantile from 0.90 → 0.99 makes cluster 11's singlet count swing wildly: 10 → 411 → 902 → 443. There is no stable quantile where the singlet fraction settles down — it's whipsawing, which means whatever number you pick is largely arbitrary for this cluster.
    - Of the current 411 "singlet" calls, 27% have their 2nd-strongest CMO within 0.1 CLR units of also clearing threshold — i.e., a quarter of them are one small nudge from being called doublets too.
    - My read: this isn't really a "wrong threshold" problem — it's that CMO hashing is structurally blind to doublets formed from two cells carrying the same CMO tag (same-sample doublets). Those look identical to a true singlet to the hashing assay no matter what threshold you pick, but Scrublet (which looks at transcriptome complexity, not oligo identity) still catches them. Given Scrublet's near-100% call rate here, I'd treat the whole cluster as doublet-contaminated and drop it — including the CMO "singlets" — rather than cherry-picking survivors. I'd apply the same logic to clusters 1, 3, 10, and 12.

    3. The Negatives look worth rescuing, at least in singlet-dominant clusters.
    Restricted to clusters where CMO doublet-calling is <50% (i.e., Negatives sitting among otherwise-confident cells), Negatives vs. Singlets: total_counts 880 vs. 1,176, genes 629 vs. 836, %mito 8.4% vs. 6.8%, Scrublet doublet_score 0.146 vs. 0.192 (lower, not higher). So Negatives are modestly lower-complexity but comfortably above your QC cutoffs, and — importantly — Scrublet doesn't think they're doublets either. That combination (real-cell-range QC + low doublet score + clustering with confident cells) points to "failed CMO staining" rather than "junk," so I'd keep them rather than discard.

    One outlier worth a separate look: cluster 7 has an unusually high 36.8% Negative rate (vs. single digits/low teens elsewhere) — that's worth checking on its own; it might be a cell type that stains poorly for CMOs rather than a QC problem.

    Bottom line recommendation: drop clusters 1/3/10/11/12 wholesale (trust Scrublet's near-unanimous call over CMO hashing's structurally-blind-spot singlets), keep Negatives elsewhere, and take a closer look at cluster 7 specifically before deciding on its negatives.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Cluster 7 conclusions
    - Cluster 7's overall median %mito is 19.5%, vs. 6.5% for the rest of the dataset — nearly 3x higher. This isn't just the Negatives; the whole cluster runs hot on mito.
    - Breaking it down by CMO status: 86.2% of cluster 7's Negatives, but also 51.6% of its Singlets and 51.0% of its Doublets, exceed the 15% mito cutoff we discussed earlier as a debris/damage signal.
    - Cluster 7's Negatives specifically: median %mito 23.7%, total_counts 655, n_genes 316 — all worse than cluster 7's own Singlets (15.4% mito, 820 counts, 536 genes), and their CMO signal isn't borderline either (median gap to threshold is larger than for negatives elsewhere in the dataset, i.e. these are confidently, not marginally, negative).

    Conclusion for cluster 7: this looks like a genuine high-mitochondrial / damaged-nuclei cluster rather than a "real cell type that just failed CMO staining." Even its Singlet-labeled cells are majority high-mito. That's the opposite situation from the other singlet-dominant clusters, where Negatives looked comparably healthy to their Singlets. So I wouldn't rescue cluster 7's Negatives the way I'd suggest for the rest of the dataset — instead, I'd treat cluster 7 as a candidate for exclusion (or at least apply the 15% mito filter retroactively), since roughly half the cluster fails that bar regardless of CMO status.
    """)
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
    ok so overall we said, we remove the clusters with more than 50% doublets, both scrublet and cmo agrees
      about those(1/3/10/11/12). Cluster 7 is removed because of high mito and this calls for a general filter
      at 15% mito. It seems like we can trust scrublet calls so we will remove all the predicted doublet. We
      will remove all the negatives as per cmo tag and finally what's left we believe it to be real and we are
      going to rescue the cells without a tag by assigning the one of the cluster.


    12.8% of singlets (2,697 of 21,117) have an own-tag that disagrees with their cluster's consensus tag, and it also tells us the rescue heuristic for Negatives will be wrong roughly 1 in 8 times. Let's look at which clusters drive that mismatch rate.

    That mismatch is almost entirely concentrated in the clusters you're already excluding — 46–58% mismatch in clusters 1, 3, 7, 10, 11, 12 (which makes sense: doublet-dominated/damaged clusters have mixed timepoint composition, so a majority vote there is close to a coin flip). Within the clusters you're actually keeping, the rescue tag is far more trustworthy: 1.3–7.5% mismatch for clusters 0, 2, 5, 6, 8, 9.

    One thing worth flagging: cluster 4 sits at 24.4% mismatch despite not being on your exclusion list (its doublet rates were 24.1% CMO / 48.5% Scrublet — both under the 50% cutoff, so it wasn't flagged). That's noticeably worse than the other kept clusters and might deserve a second look, though I'll leave it in for now since you didn't ask me to change the exclusion criteria.
    """)
    return


app._unparsable_cell(
    r"""
    _outfile = _project_root / "results/channel1/filtered_barcodes.txt")
    _outfile.write_text("\n".join(adata.obs_names[adata.obs["doublet_score"] < 0.15].to_list()))
    """,
    column=None, disabled=True, hide_code=False, name="_"
)


@app.cell(hide_code=True)
def cluster11_investigation_intro(mo):
    mo.md(r"""
    ## Cluster 11 doublet investigation

    Both methods flag cluster 11 as heavily doublet-dominated: Scrublet calls ~100%, CMO hashing calls ~75%. Cells below dig into (1) whether the CMO "Singlet" calls in that cluster are threshold-sensitive borderline cases, and (2) whether "Negative" barcodes that cluster confidently with assigned cells look transcriptionally like real cells worth rescuing.
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
    ).round(1)
    leiden_doublet_summary.sort_values("pct_cmo_doublet", ascending=False)
    return (leiden_doublet_summary,)


@app.cell(hide_code=True)
def cluster11_threshold_sensitivity(adata_flt, clr, np, pd):
    # Are cluster 11's "Singlet" calls robust, or do they flip to Doublet under a
    # slightly stricter (higher-quantile) threshold? If the CMO thresholds are the
    # problem, singlet counts here should collapse quickly as the quantile rises.
    _cluster11_mask = (adata_flt.obs["leiden"] == "11").to_numpy()
    _clr_11 = clr[_cluster11_mask]

    _rows = []
    for _q in [0.90, 0.95, 0.975, 0.99]:
        _t = np.quantile(clr, _q, axis=0)
        _pos = _clr_11 > _t
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
def cluster11_singlet_margin(adata_flt, clr, np, plt, thresholds):
    # For cluster 11's current "Singlet" calls, how close is the 2nd-strongest CMO to
    # also crossing its own threshold? A small (near-zero) gap means the call is
    # borderline; a large gap means the singlet call is robust regardless of quantile.
    _cluster11_mask = (adata_flt.obs["leiden"] == "11").to_numpy()
    _singlet_11_mask = _cluster11_mask & (adata_flt.obs["cmo_status_scanpy"] == "Singlet").to_numpy()
    _clr_singlet_11 = clr[_singlet_11_mask]
    _margin_to_threshold = -np.sort(-(_clr_singlet_11 - thresholds), axis=1)[:, 1]  # 2nd-highest (value - threshold)

    plt.figure(figsize=(6, 4))
    plt.hist(_margin_to_threshold, bins=40)
    plt.axvline(0, color="red", linestyle="--", label="threshold (0 = right at cutoff)")
    plt.xlabel("2nd-strongest CMO: CLR value minus its threshold")
    plt.ylabel("Number of cluster-11 singlets")
    plt.title("How close is the 2nd CMO to also being called positive?")
    plt.legend()
    plt.tight_layout()
    plt.gcf()
    return


@app.cell(hide_code=True)
def negative_vs_singlet_qc(adata_flt, leiden_doublet_summary, sc):
    # Are "Negative" barcodes that land in largely-Singlet clusters real cells that
    # simply failed CMO hash detection (similar QC profile to their cluster's Singlets),
    # or lower-quality debris (worse QC profile)? Restrict to clusters where CMO hashing
    # is NOT doublet-dominated, i.e. where Negatives are clustering with confidently
    # assigned cells rather than with other ambiguous/doublet barcodes.
    _singlet_dominant_clusters = leiden_doublet_summary.loc[leiden_doublet_summary["pct_cmo_doublet"] < 50].index
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
def cluster7_mito_filter_impact(adata_flt, mo):
    # How big is cluster 7, and what happens if we retroactively apply a 15% mito cutoff?
    _mask7 = (adata_flt.obs["leiden"] == "7")
    _mito_pass = adata_flt.obs["pct_counts_mt"] <= 15

    _n_total = adata_flt.n_obs
    _n_cluster7 = int(_mask7.sum())
    _n_cluster7_pass_mito = int((_mask7 & _mito_pass).sum())

    _n_other = _n_total - _n_cluster7
    _n_other_pass_mito = int((~_mask7 & _mito_pass).sum())

    mo.md(f"""
    **Cluster 7 size:** {_n_cluster7:,} of {_n_total:,} cells ({_n_cluster7 / _n_total:.1%} of the whole dataset)

    **Effect of a retroactive `pct_counts_mt <= 15` filter:**

    | Group | Before | After | Lost | Lost % |
    |---|---|---|---|---|
    | Cluster 7 | {_n_cluster7:,} | {_n_cluster7_pass_mito:,} | {_n_cluster7 - _n_cluster7_pass_mito:,} | {(_n_cluster7 - _n_cluster7_pass_mito) / _n_cluster7:.1%} |
    | Rest of dataset | {_n_other:,} | {_n_other_pass_mito:,} | {_n_other - _n_other_pass_mito:,} | {(_n_other - _n_other_pass_mito) / _n_other:.1%} |
    | **Total** | {_n_total:,} | {_n_cluster7_pass_mito + _n_other_pass_mito:,} | {_n_total - (_n_cluster7_pass_mito + _n_other_pass_mito):,} | {(_n_total - (_n_cluster7_pass_mito + _n_other_pass_mito)) / _n_total:.1%} |

    Cluster 7 accounts for {(_n_cluster7 - _n_cluster7_pass_mito) / (_n_total - (_n_cluster7_pass_mito + _n_other_pass_mito)):.1%} of all cells that a 15% mito cutoff would remove dataset-wide, despite being only {_n_cluster7 / _n_total:.1%} of the dataset — i.e. a mito filter would disproportionately clear out cluster 7.
    """)
    return


@app.cell(hide_code=True)
def apply_strict_qc_and_rescue(adata_flt, leiden_doublet_summary, np):
    # --- Consolidated QC decision ---------------------------------------------
    # 1. Doublet-dominated clusters: both CMO hashing AND Scrublet independently
    #    call >50% doublet -> drop the whole cluster.
    # 2. Cluster 7: distinct high-mitochondrial / damaged-nuclei cluster -> drop
    #    the whole cluster regardless of its own per-cell mito value.
    # 3. General mito filter: pct_counts_mt > 15% -> drop, dataset-wide.
    # 4. Scrublet-predicted doublets -> drop.
    # 5. CMO-hashing "Doublet" calls -> drop, even where Scrublet disagrees.
    _doublet_cluster_mask = leiden_doublet_summary["pct_cmo_doublet"].gt(50) & leiden_doublet_summary["pct_scrublet_doublet"].gt(50)
    doublet_dominated_clusters = leiden_doublet_summary.index[_doublet_cluster_mask].tolist()
    high_mito_cluster = "7"

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
    return (cluster_consensus_tag,)


@app.cell(hide_code=True)
def cluster_tag_mismatch_check(adata_flt, cluster_consensus_tag):
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

    print(f"{_n_mismatch} of {_n_singlets} singlets ({_n_mismatch / _n_singlets:.1%}) have an own-tag that disagrees with their cluster's consensus tag.")
    _by_cluster.sort_values("pct_mismatch", ascending=False)
    return


@app.cell(hide_code=True)
def write_qc_annotations_tsv(adata_flt, project_root):
    _outfile = project_root / "results/channel1/adata_flt_qc_annotations.tsv"
    _outfile.parent.mkdir(parents=True, exist_ok=True)
    adata_flt.obs.to_csv(_outfile, sep="\t")
    _outfile
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
def _(mo):
    mo.md(r"""
    ### Cluster 4 conclusions

    - Timepoint composition of cluster 4's singlets: 75.5% d3, 18.3% d4, and a small tail of d2/d1/d0 (6%). That 24.5% "non-majority" fraction lines up almost exactly with the 24.4% mismatch rate flagged earlier — so the mismatch isn't noise, it's explained by cluster 4 being a genuine mixed d3/d4 transition population, not a technical artifact.
    -
    - QC metrics are unremarkable: total_counts and n_genes are essentially identical to the rest of the dataset; %mito is a bit elevated (8.2% vs 6.9%) but nowhere near cluster 7's territory (19.5%); doublet_score is somewhat higher (0.31 vs 0.25, right around the Scrublet threshold) but CMO's doublet rate here (24.1%) is well under the 50% cutoff used for the other excluded clusters.

    Conclusion: cluster 4 doesn't look like a contamination/doublet problem — it looks like a real transitional cell state spanning the d3→d4 boundary. I wouldn't exclude it; the elevated mismatch rate is just the expected consequence of a majority-vote heuristic applied to a genuinely mixed population. No further action needed here unless you want to split cluster 4's d3/d4 cells apart some other way (e.g. re-clustering at higher resolution) rather than treating it as one group.
    """)
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
def _(assigned_cells, chrom_dict, project_root, snap):
    atac_adata = snap.pp.import_fragments(
        project_root / "data/fragments/10x_5timepoints_channel1.fragments.tsv.gz",
        sorted_by_barcode=False,
        chrom_sizes=chrom_dict,
        whitelist=assigned_cells,
    )
    snap.metrics.tsse(atac_adata, project_root / "annotations/gencode.v43.primary_assembly.annotation.gtf.gz")
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
def _(alt, df, min_frag, min_tsse, mo):
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

    # base scatter
    base = alt.Chart(plot_df).encode(
        x=alt.X("n_fragment:Q", scale=alt.Scale(type="log"), title="Unique fragments"),
        y=alt.Y("tsse:Q", title="TSSE"),
        color=alt.Color("pass_qc:N", title="Pass QC")
    )

    scatter = base.mark_circle(size=10, opacity=0.3)

    # marginal distributions
    top_hist = alt.Chart(plot_df).mark_bar(opacity=0.5).encode(
        x=alt.X("n_fragment:Q", scale=alt.Scale(type="log")),
        y=alt.Y("count()"),
        color=alt.Color("pass_qc:N", legend=None)
    ).properties(height=100)

    right_hist = alt.Chart(plot_df).mark_bar(opacity=0.5).encode(
        y=alt.Y("tsse:Q"),
        x=alt.X("count()"),
        color=alt.Color("pass_qc:N", legend=None)
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


@app.cell
def _(atac_adata):
    atac_adata.obs["timepoint_scanpy"].value_counts()
    return


if __name__ == "__main__":
    app.run()
