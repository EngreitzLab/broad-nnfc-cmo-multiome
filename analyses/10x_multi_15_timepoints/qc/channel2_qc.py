import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Introduction

    This notebook processes the 10x multiome data (**channel2**) for the endothelial differentiation time course: 15 timepoints (d0-d5, with arterial and venous branches from d2.5 onward), 2 to 3 biological replicates per timepoint. Samples were multiplexed using the MULTI-seq technique with 23 CMO (Cell Multiplexing Oligo) barcodes; unlike the 5-timepoints dataset, CMO barcodes are reused across channels for different timepoint/replicate combinations, so a CMO id alone does not identify a timepoint (see `metadata/10X_multi_15_timepoints_CMO_design.tsv` for the per-channel mapping). GEX and CMO reads were already demultiplexed into separate FASTQs by the sequencer, so `splitcode` was not needed upstream. Data was processed with kallisto/bustools (GENCODE v43 annotation) for RNA, and CMO quantification was performed with the `kite` workflow -- see `analyses/10x_multi_15_timepoints/README.md` for the full processing pipeline.

    The goals of this notebook are to:

    1. Perform quality control filtering on the RNA data
    2. Run CMO hash classification (mimicking `Seurat::HTODemux`) on QC-passing cells
    3. Assign each cell barcode to a CMO (and therefore to a timepoint/replicate)
    4. Process the corresponding ATAC data with `snapatac2`

    This dataset is not yet on the IGVF portal; all data is read locally from `data/10X_multiome_15_timepoints/channel2/`.
    """)
    return


@app.cell(hide_code=True)
def imports():
    import altair as alt
    import anndata as ad
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    import scanpy as sc
    import seaborn as sns
    import sys
    import gzip
    import re

    from dotenv import find_dotenv
    from pathlib import Path

    # get the project root path using the '.env' file.
    _env_path = find_dotenv(usecwd=True)
    project_root = Path(_env_path).parent

    sys.path.insert(0, str(project_root / "scripts" / "python"))
    from utils import gtf_to_gene_metadata

    # PFlog (shifted centered log-ratio) normalization instead of log1p(CP10K):
    # https://www.biorxiv.org/content/10.1101/2022.05.06.490859
    # jointly stabilizes technical variance, normalizes for sequencing depth, and
    # preserves within-cell gene ranking (monotonicity) via a data-calibrated
    # pseudocount and CLR centering, rather than a fixed round-number pseudocount.
    import scclr

    # VegaFusion pre-aggregates chart data in Python before sending it to the
    # browser, raising Altair's default 5,000-row embed limit. This is enabled
    # notebook-wide (rather than toggled per-cell) since alt.data_transformers is
    # global module state, not a per-chart setting.
    alt.data_transformers.enable("vegafusion")

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
    def pca_axis_title(adata, pc_idx):
        _pct = adata.uns["pca"]["variance_ratio"][pc_idx] * 100
        return f"PC{pc_idx + 1} ({_pct:.1f}% var.)"

    return (
        Path,
        ad,
        alt,
        gtf_to_gene_metadata,
        gzip,
        np,
        okabe_ito_palette,
        pca_axis_title,
        pd,
        plt,
        project_root,
        re,
        sc,
        scclr,
        sns,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Imports and palettes
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Local data paths

    This dataset is not yet on the IGVF portal, so there is no accession-based download step here. The RNA counts (kallisto/bustools h5ad) and ATAC fragments were produced locally by the pipeline in `analyses/10x_multi_15_timepoints/README.md` and already live under `data/10X_multiome_15_timepoints/channel2/`.
    """)
    return


@app.cell
def _(Path, project_root):
    ch2_data_root_path = Path(project_root / "data/10X_multiome_15_timepoints/channel2/")

    ch2_outdir_root = Path(project_root / "results/10X_multiome_15_timepoints/channel2")
    return ch2_data_root_path, ch2_outdir_root


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Loading gene metadata

    Per-gene metadata parsed directly from the GENCODE v43 GTF via `gtf_to_gene_metadata` (`scripts/python/utils.py`), covering all annotated genes.

    Columns loaded:

    - `chr` - chromosome
    - `start` / `end` - gene body coordinates (BED 0-based, half-open)
    - `gene_symbol` - gene symbol
    - `strand` - `+` or `-`
    - `gene_id` - Ensembl gene ID (version suffix stripped)
    - `gene_type` - GENCODE gene biotype
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Gene metadata (GENCODE GTF)

    Reuses the GENCODE v43 GTF already downloaded from IGVF for the 5-timepoints analysis (`annotations/IGVFFI9573KOZR.gtf.gz`), shared across analyses in this repo rather than re-downloaded here.
    """)
    return


@app.cell
def get_igvf_gencode(gtf_to_gene_metadata, project_root):
    igvf_gencode_gtf_path = project_root / "annotations/IGVFFI9573KOZR.gtf.gz"

    # Show only the path relative to the repo, not the full local filesystem path
    igvf_gencode_gtf_path.relative_to(project_root)

    gene_metadata_df = gtf_to_gene_metadata(igvf_gencode_gtf_path)
    gene_metadata_df
    return gene_metadata_df, igvf_gencode_gtf_path


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Loading the counts in scanpy

    We are going to load the total counts as produced by `kallisto-bustool` and annotate the genes with the metadata loaded the cell above. The default matrix stored in `.X` is the total counts: the sum of `ambigous`, `mature`, and `nascent` matrices.
    """)
    return


@app.cell
def load_h5_counts(
    ad,
    ch2_data_root_path,
    gene_metadata_df,
    gzip,
    igvf_gencode_gtf_path,
    re,
):
    try:
        _h5ad_path = ch2_data_root_path / "h5ad/10x_15timepoints_channel2.h5ad"
        adata = ad.read_h5ad(_h5ad_path)

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
        _gene_symbol = adata.var["gene_symbol"].fillna("").astype(str)
        _pc_mask = (adata.var["gene_type"] == "protein_coding").to_numpy()

        # The protein-coding-only TSS metadata merge above leaves mito rRNAs/tRNAs
        # without a gene_symbol, so a plain "MT-" symbol prefix match misses them.
        # Use the GENCODE v43 GTF (loaded above) as ground truth instead: every
        # gene actually annotated on chrM, not just the ones with a MANE symbol.
        _chrm_gene_ids = set()
        with gzip.open(igvf_gencode_gtf_path, "rt") as _gtf:
            for _line in _gtf:
                if _line.startswith("#"):
                    continue
                _fields = _line.split("\t")
                if _fields[0] == "chrM" and _fields[2] == "gene":
                    _gene_id_match = re.search(r'gene_id "([^"]+)"', _fields[8])
                    _chrm_gene_ids.add(_gene_id_match.group(1).split(".")[0])

        adata.var["mt"] = _gene_symbol.str.startswith("MT-") | adata.var["gene_id_base"].isin(_chrm_gene_ids)
        adata.var["ribo"] = _gene_symbol.str.startswith(("RPS", "RPL"))
        adata.var["hb"] = _gene_symbol.str.contains(r"^HB(?!P)", regex=True)
        # filtering to protein-coding genes not overlapping the above categories
        adata.var["pc_flt"] = _pc_mask & ~adata.var["mt"] & ~adata.var["ribo"] & ~adata.var["hb"]

        # variable to maintain marimo reactivity
        adata_ready = True
    except:
        adata_ready = False

    adata
    return adata, adata_ready


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
def _(adata, adata_ready, mo, sc):
    mo.stop(not adata_ready)

    # calculate QC metrics
    def _run_calculate_qc_metrics(adata):
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb", "pc_flt"], inplace=True, log1p=True)
        return True

    # Captured so downstream cells can reference this name and get a real,
    # trackable marimo dependency edge instead of just sharing a ref to "adata"
    # -- see the marimo-pair race-condition note.
    qc_metrics_computed = _run_calculate_qc_metrics(adata)

    mo.md(f"**DONE**: Quality control metrics calculation finished")
    return (qc_metrics_computed,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Knee plot
    """)
    return


@app.cell
def plot_knee_plot(adata, alt, mo, np, okabe_ito_palette, pd):
    # Total UMI counts per barcode (from QC metrics, computed above). Drop
    # zero-count barcodes since they are undefined on a log-scale y-axis.
    _counts_sorted = np.sort(adata.obs["total_counts"].to_numpy())[::-1]
    _counts_sorted = _counts_sorted[_counts_sorted > 0]
    knee_full_counts = _counts_sorted
    _ranks = np.arange(1, len(_counts_sorted) + 1)

    # Log-spaced subsample for a light, smooth curve (full barcode set is >600k points)
    _n_points = 2000
    _log_idx = np.unique(np.logspace(0, np.log10(len(_ranks) - 1), _n_points).astype(int))
    _knee_df = pd.DataFrame({"rank": _ranks[_log_idx], "n_umis": _counts_sorted[_log_idx]})

    knee_umi_selection = alt.selection_interval(encodings=["y"], value={"y": [1000, 10_000]})

    _knee_chart = alt.Chart(_knee_df).mark_circle(size=15, color=okabe_ito_palette[0]).encode(
        x=alt.X("rank:Q", scale=alt.Scale(type="log"), title="Cell rank"),
        y=alt.Y("n_umis:Q", scale=alt.Scale(type="log"), title="Total UMIs"),
        tooltip=[alt.Tooltip("rank:Q", title="Rank"), alt.Tooltip("n_umis:Q", title="# UMIs", format=",")],
    ).add_params(knee_umi_selection).properties(width=550, height=400, title="Knee plot")

    # mo.ui.altair_chart calls chart.to_json() without format="vega", which
    # conflicts with the vegafusion transformer enabled globally in `imports`
    # (it requires that format explicitly). This chart's data is already a
    # ~2000-point downsample, so vegafusion is not needed here anyway --
    # temporarily drop back to the default transformer just for this widget.
    with alt.data_transformers.enable("default"):
        knee_chart_ui = mo.ui.altair_chart(_knee_chart, chart_selection=False, legend_selection=False)
    knee_chart_ui
    return knee_chart_ui, knee_full_counts


@app.cell(hide_code=True)
def knee_plot_brush_count(knee_chart_ui, knee_full_counts, mo):
    # Read the brush's y-range (falls back to the default 1000-10,000 before any
    # interaction, since the interval's initial "value" is not reported back until
    # the user actually drags it) and count how many barcodes it covers, using the
    # full (non-downsampled) counts array for an accurate number.
    _selections = knee_chart_ui.selections
    if _selections:
        _lo, _hi = next(iter(_selections.values()))["n_umis"]
    else:
        _lo, _hi = 1000, 10_000

    _n_in_range = int(((knee_full_counts >= _lo) & (knee_full_counts <= _hi)).sum())

    mo.md(f"**{_n_in_range:,} barcodes** fall between **{_lo:,.0f}** and **{_hi:,.0f}** total UMIs (drag the brush above to adjust).")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Flag low quality barcodes
    """)
    return


@app.cell(hide_code=True)
def mask_filter_cells(adata, mo, np, pd):
    # QC pass mask mirroring the cell-level filters applied below
    _max_counts_cutoff = 15_000
    _min_counts_cutoff = 1000

    adata.obs["pass_min_umi_filter"] = adata.obs["total_counts"] >= _min_counts_cutoff
    adata.obs["pass_max_umi_filter"] = adata.obs["total_counts"] <= _max_counts_cutoff
    adata.obs["pass_min_gene_filter"] = adata.obs["n_genes_by_counts"] >= 200

    pre_mito_qc_mask = (
        adata.obs["pass_min_umi_filter"] & adata.obs["pass_max_umi_filter"] & adata.obs["pass_min_gene_filter"]
    )
    _mito_vals = adata.obs.loc[pre_mito_qc_mask, "pct_counts_mt"].to_numpy()
    _mito_median = np.median(_mito_vals)
    _mito_mad = np.median(np.abs(_mito_vals - _mito_median))
    mito_mad_cutoff = _mito_median + 3 * _mito_mad  # kept for reference/comparison only
    mito_mad_cutoff_2_5 = _mito_median + 2.5 * _mito_mad
    adata.obs["pass_mito_mad_filter"] = adata.obs["pct_counts_mt"] <= mito_mad_cutoff
    adata.obs["pass_mito_mad_filter_2_5"] = adata.obs["pct_counts_mt"] <= mito_mad_cutoff_2_5

    _steps = [
        (f"min_counts >= {_min_counts_cutoff}", adata.obs["pass_min_umi_filter"]),
        (f"max_counts <= {_max_counts_cutoff}", adata.obs["pass_max_umi_filter"]),
        ("min_genes >= 200", adata.obs["pass_min_gene_filter"]),
        (f"pct_counts_mt <= {mito_mad_cutoff_2_5:.1f} (median + 2.5 MADs, computed on the population passing the other lenient filters)", adata.obs["pass_mito_mad_filter_2_5"]),
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

    # Exposed as its own top-level name (not just an adata.obs mutation) so
    # downstream cells get a real, trackable marimo dependency edge instead of
    # just sharing a ref to "adata" -- see the marimo-pair race-condition note.
    pass_lenient_qc_mask = _remaining_mask
    adata.obs["pass_lenient_qc"] = pass_lenient_qc_mask

    mo.md(f"""
    **QC filter breakdown** (min_counts={_min_counts_cutoff:,}, max_counts={_max_counts_cutoff:,}, min_genes=200, mito_mad_cutoff_2_5={mito_mad_cutoff_2_5:.1f}% [3-MAD reference cutoff: {mito_mad_cutoff:.1f}%])

    | Filter step | All barcodes |
    |---|---|
    {chr(10).join(_rows)}

    **Net:** {int(pass_lenient_qc_mask.sum()):,} of {adata.n_obs:,} barcodes survive this QC filter.
    """)
    return mito_mad_cutoff_2_5, pass_lenient_qc_mask, pre_mito_qc_mask


@app.cell(hide_code=True)
def checking_mito_content(
    adata,
    mito_mad_cutoff_2_5,
    okabe_ito_palette,
    plt,
    pre_mito_qc_mask,
    sns,
):
    pre_mito_qc_mask  # ran after the pre-mito lenient QC mask (min_counts/max_counts/min_genes) was computed
    mito_mad_cutoff_2_5  # ran after the 2.5-MAD mito cutoff was computed

    # Plotted on the pre-mito population (counts + genes filters only), not
    # pass_lenient_qc, so this stays a complete picture of what the mito filter
    # below actually removes, rather than only showing survivors.
    _df = adata.obs.loc[pre_mito_qc_mask, ["total_counts", "pct_counts_mt"]].copy()

    _g = sns.JointGrid(data=_df, x="total_counts", y="pct_counts_mt", height=6)
    _scatter = _g.ax_joint.scatter(
        _df["total_counts"], _df["pct_counts_mt"],
        c=_df["pct_counts_mt"], cmap="cividis", s=5, alpha=0.5,
    )

    _g.ax_marg_x.hist(_df["total_counts"], bins=100, color="gray", alpha=0.7)
    _g.ax_marg_y.hist(_df["pct_counts_mt"], bins=100, orientation="horizontal", color="gray", alpha=0.7)

    _g.ax_joint.axvline(x=1000, color="blue", linestyle="--")
    _g.ax_joint.annotate(
        "1000 UMI",
        xy=(1000, _df["pct_counts_mt"].max()),
        xytext=(5, -5),
        textcoords="offset points",
        va="top",
        ha="left",
        color="blue",
        fontweight="bold",
    )

    # 2.5-MAD mito cutoff actually applied in mask_filter_cells, drawn here for
    # reference against the full (pre-mito-filter) population.
    _g.ax_joint.axhline(y=mito_mad_cutoff_2_5, color=okabe_ito_palette[0], linestyle="--")
    _g.ax_joint.annotate(
        f"{mito_mad_cutoff_2_5:.1f}% (2.5 MAD, applied cutoff)",
        xy=(_df["total_counts"].max(), mito_mad_cutoff_2_5),
        xytext=(20, 5),
        textcoords="offset points",
        va="bottom",
        ha="right",
        color=okabe_ito_palette[0],
        fontweight="bold",
    )

    _g.ax_joint.set_xlabel("Total UMI counts")
    _g.ax_joint.set_ylabel("% mitochondrial counts")
    _g.figure.suptitle("Counts/genes-filtered barcodes: total counts vs. % mitochondrial", y=1.02)
    _g.figure.text(0.5, 0.96, f"n = {len(_df):,} barcodes", ha="center", fontsize=9, color="dimgray")

    # Dedicated axes for a horizontal colorbar, placed below the joint plot
    # without stealing space from it (keeps marginal alignment intact).
    _joint_pos = _g.ax_joint.get_position()
    _cax = _g.figure.add_axes([_joint_pos.x0, _joint_pos.y0 - 0.1, _joint_pos.width, 0.03])
    _cbar = _g.figure.colorbar(_scatter, cax=_cax, orientation="horizontal")
    _cbar.set_label("% mitochondrial counts")
    # Colorbar swatch should read at full opacity even though the scatter points
    # are drawn with alpha=0.5 -- otherwise the legend inherits that
    # transparency and looks washed out relative to the actual color scale.
    _cbar.solids.set_alpha(1)

    plt.show()
    return


@app.cell(hide_code=True)
def qc_round1_summary(adata, mito_mad_cutoff_2_5, mo, pre_mito_qc_mask):
    mo.md(f"""
    ### Lenient QC summary

    Starting from {adata.n_obs:,} raw barcodes, the level-1 filters (`min_counts >= 1000`, `max_counts <= 15,000`, `min_genes >= 200`) leave **{int(pre_mito_qc_mask.sum()):,} barcodes**:

    - `pass_min_umi_filter`: {int(adata.obs["pass_min_umi_filter"].sum()):,} pass (most of the loss here is empty droplets/background)
    - `pass_max_umi_filter`: {int(adata.obs["pass_max_umi_filter"].sum()):,} pass
    - `pass_min_gene_filter`: {int(adata.obs["pass_min_gene_filter"].sum()):,} pass

    Among barcodes passing these filters, %mitochondrial content is still elevated (median {adata.obs.loc[pre_mito_qc_mask, "pct_counts_mt"].median():.1f}%, {adata.obs.loc[pre_mito_qc_mask, "pct_counts_mt"].gt(15).mean():.1%} above 15%) despite this being nuclei input, which should show close to 0% mito. See the note below.

    A mito filter (`pct_counts_mt <= {mito_mad_cutoff_2_5:.1f}%`, median + 2.5 MADs of the level-1-passing population) is applied on top of the level-1 filters, leaving **{int(adata.obs["pass_lenient_qc"].sum()):,} barcodes** (`pass_lenient_qc`).
    """)
    return


@app.cell(hide_code=True)
def notes_on_mito_content(mo):
    mo.md(r"""
    The 10X multi-ome protocol requires nuclei in input, thus we would expect mitochondrial percentage to be close to 0. We see a lot of barcodes with high mitochondrial content. It could be due to incomplete nuclei isolation or ambient/cytoplasmic contamination. We filter on a median + 2.5 MADs mito cutoff (computed on the level-1-passing population) rather than a fixed threshold, to adapt to this dataset's own distribution.
    """)
    return


@app.cell
def filter_cells(adata, pass_lenient_qc_mask):
    pass_lenient_qc_mask  # ran after the lenient QC mask (incl. mito filter) was computed

    adata_flt = adata[adata.obs["pass_lenient_qc"]].copy()
    adata_flt.layers["counts"] = adata_flt.X.copy()
    adata_flt
    return (adata_flt,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Violin plots for different metrics
    """)
    return


@app.cell
def qc_violin_metric_dropdown_ui(mo, qc_metrics_computed):
    mo.stop(not qc_metrics_computed)

    # Defined here, displayed together with the plot below (basic_qc_violin_plots),
    # since its .value cannot be read in the same cell that creates it.
    _qc_metric_options = {
        "Number of genes": "n_genes_by_counts",
        "Total UMI counts": "total_counts",
        "% mitochondrial": "pct_counts_mt",
        "% ribosomal": "pct_counts_ribo",
        "% hemoglobin": "pct_counts_hb",
        "% protein coding": "pct_counts_pc_flt",
    }
    qc_violin_metric_dropdown = mo.ui.dropdown(
        options=_qc_metric_options,
        value="Total UMI counts",
        label="QC metric",
    )
    return (qc_violin_metric_dropdown,)


@app.cell
def basic_qc_violin_plots(
    adata_flt,
    alt,
    mo,
    np,
    okabe_ito_palette,
    pd,
    qc_metrics_computed,
    qc_violin_metric_dropdown,
):
    mo.stop(not qc_metrics_computed)

    from scipy.stats import gaussian_kde as _gaussian_kde_qc_violin

    # KDE on a fixed grid rather than Altair's transform_density: VegaFusion
    # cannot pre-aggregate transform_density, so it would embed the full raw
    # dataframe and trip the output-too-large limit.
    _metric = qc_violin_metric_dropdown.value
    _label = qc_violin_metric_dropdown.selected_key

    _values = adata_flt.obs[_metric].to_numpy()
    _grid = np.linspace(_values.min(), _values.max(), 200)
    _density = _gaussian_kde_qc_violin(_values)(_grid)
    _density_df = pd.DataFrame({_metric: _grid, "density": _density})

    # Mirror the density manually (x0/x2) instead of Altair's stack="center":
    # with a per-row groupby of size 1, Vega-Lite's stack transform did not
    # produce a symmetric shape here, and the combined x scale (quantitative on
    # both channels) also flipped the tick mark to a vertical orientation stuck
    # at the domain edge instead of a centered horizontal dash.
    _density_df["_x0"] = -_density_df["density"] / 2
    _density_df["_x1"] = _density_df["density"] / 2

    _violin = alt.Chart(_density_df).mark_area(orient="horizontal", color=okabe_ito_palette[8], opacity=0.7).encode(
        y=alt.Y(f"{_metric}:Q", title=_label),
        x=alt.X("_x0:Q", title=None, axis=None),
        x2="_x1:Q",
    )

    _median = np.median(_values)
    _median_tick = alt.Chart(pd.DataFrame({"y": [_median], "_x0": [0.0]})).mark_tick(
        color=okabe_ito_palette[0], thickness=2, size=80, orient="horizontal",
    ).encode(y=alt.Y("y:Q"), x=alt.X("_x0:Q"))

    _chart = (_violin + _median_tick).properties(
        title=f"{_label} (n={len(_values):,} barcodes, median={_median:,.1f})",
        width=250, height=400,
    ).configure_view(strokeWidth=0)

    mo.vstack([qc_violin_metric_dropdown, _chart])
    return


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
    try:
        scclr.pp.pflog(adata_flt, target="auto")
        is_data_normalized = True
    except:
        is_data_normalized = False

    adata_flt.uns["pflog"]
    return (is_data_normalized,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Find highly variable genes
    """)
    return


@app.cell
def _(adata_flt, is_data_normalized, mo, sc):
    mo.stop(not is_data_normalized)  # ran after PFlog normalization

    def _run_hvg(adata):
        sc.pp.highly_variable_genes(adata, layer="pflog", n_top_genes=2000)
        return True

    hvg_computed = _run_hvg(adata_flt)
    return (hvg_computed,)


@app.cell
def _(adata_flt, hvg_computed, sc):
    hvg_computed  # ran after highly variable genes were computed
    sc.pl.highly_variable_genes(adata_flt)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## PCA
    """)
    return


@app.cell
def _(adata_flt, hvg_computed, mo, scclr):
    mo.stop(not hvg_computed)

    def _run_pca(adata, n_comps, ncv):
        scclr.tl.pca(adata, n_comps=n_comps, ncv=ncv)
        return True

    _ncomps = 50
    _ncv = 2 * _ncomps + 1
    pca_computed = _run_pca(adata_flt, _ncomps, _ncv)
    mo.md(f"Computed PCA using {_ncomps} components and {_ncv} Lanczos vectors")
    return (pca_computed,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Plot variance ratio
    """)
    return


@app.cell
def _(adata_flt, mo, pca_computed, sc):
    mo.stop(not pca_computed) # ran after PCA
    sc.pl.pca_variance_ratio(adata_flt, n_pcs=15, log=False)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Plot PCA
    """)
    return


@app.cell
def pca_axis_setup(adata_flt, mo, pca_axis_title, pca_computed):
    mo.stop(not pca_computed)  # ran after PCA

    # ------------ Change the PC numbers here--------------------#
    _pc_number_x_axis = 0
    _pc_number_y_axis = 1
    # -----------------------------------------------------------#

    pca_x_title = pca_axis_title(adata_flt, _pc_number_x_axis)
    pca_y_title = pca_axis_title(adata_flt, _pc_number_y_axis)
    return


@app.cell(hide_code=True)
def pca_pc1_pc2_colored(adata_flt, pca_computed, sc):
    pca_computed  # ran after PCA
    sc.pl.pca(
        adata_flt,
        color=["pct_counts_mt", "pct_counts_pc_flt", "total_counts"],
        dimensions=[(0, 1)],
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
def compute_knn_neighbors(adata_flt, mo, pca_computed, sc):
    mo.stop(not pca_computed) # ran after PCA

    def _run_neighbors(adata):
        sc.pp.neighbors(adata, random_state=0)
        return True

    neighbors_computed = _run_neighbors(adata_flt)
    return (neighbors_computed,)


@app.cell
def leiden_clustering(adata_flt, mo, neighbors_computed, sc):
    mo.stop(not neighbors_computed)  # ran after the neighbor graph was computed

    # Using the igraph implementation and a fixed number of iterations can be significantly faster,
    # especially for larger datasets. random_state is pinned explicitly (rather than
    # relying on scanpy's default) so cluster identities are reproducible run-to-run;
    # note this only guards against algorithm-internal nondeterminism, not against
    # genuine changes to the input data (QC, CMO assignment, etc.), which will still
    # legitimately shift which cells land in which cluster.
    def _run_leiden(adata):
        sc.tl.leiden(adata, flavor="igraph", resolution=0.5, n_iterations=2, random_state=0)
        return True

    leiden_computed = _run_leiden(adata_flt)
    return (leiden_computed,)


@app.cell
def mito_pct_per_leiden_cluster(
    adata_flt,
    alt,
    leiden_computed,
    okabe_ito_palette,
):
    leiden_computed  # ran after leiden clustering

    # % mito per leiden cluster, to spot a mito-driven cluster before it gets a
    # dedicated deep-dive (see the cluster investigation sections below).
    _cluster_order = sorted(adata_flt.obs["leiden"].cat.categories, key=int)

    _box = alt.Chart(adata_flt.obs[["leiden", "pct_counts_mt"]]).mark_boxplot(
        color=okabe_ito_palette[5], size=25,
    ).encode(
        x=alt.X("leiden:N", title="Leiden cluster", sort=_cluster_order),
        y=alt.Y("pct_counts_mt:Q", title="% mitochondrial counts"),
    )

    _box.properties(
        title="% mitochondrial counts per leiden cluster",
        width=650, height=350,
    ).configure_view(strokeWidth=0)
    return


@app.cell
def _(adata_flt, leiden_computed, pca_computed, sc):
    pca_computed  # ran after PCA
    leiden_computed  # ran after leiden clustering
    sc.pl.pca(
        adata_flt,
        color=["leiden"],
        dimensions=[(0, 1)],
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
def compute_umap(adata_flt, mo, neighbors_computed, sc):
    mo.stop(not neighbors_computed)   # ran after the neighbor graph was computed

    def _run_umap(adata):
        sc.tl.umap(adata)
        return True

    umap_computed = _run_umap(adata_flt)
    return (umap_computed,)


@app.cell
def plot_umap(adata_flt, leiden_computed, sc, umap_computed):
    umap_computed  # ran after UMAP was computed
    leiden_computed  # ran after leiden clustering
    sc.pl.umap(adata_flt, color=["leiden", "pct_counts_mt"], cmap="cividis")
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

    adata_flt.obs["scrublet_doublet_score"] = _raw.obs["doublet_score"].to_numpy()
    adata_flt.obs["scrublet_predicted_doublet"] = _raw.obs["predicted_doublet"].to_numpy()

    # Exposed for the interactive threshold plot below
    scrublet_scores_sim = _raw.uns["scrublet"]["doublet_scores_sim"]
    scrublet_auto_threshold = _raw.uns["scrublet"].get("threshold")
    return scrublet_auto_threshold, scrublet_scores_sim


@app.cell(hide_code=True)
def scrublet_threshold_number_ui(mo, scrublet_auto_threshold):
    # Numeric threshold input, initialized to Scrublet's own automatic threshold.
    # Displayed together with the plots below (scrublet_threshold_interactive_plot),
    # not here -- this cell only defines it, since its .value cannot be read in the
    # same cell that creates it.
    scrublet_threshold_number = mo.ui.number(
        start=0.0, stop=1.0, step=0.01,
        value=float(scrublet_auto_threshold) if scrublet_auto_threshold is not None else 0.5,
        label="Doublet score threshold",
    )
    return (scrublet_threshold_number,)


@app.cell(hide_code=True)
def scrublet_threshold_interactive_plot(
    adata_flt,
    alt,
    leiden_computed,
    mo,
    np,
    okabe_ito_palette,
    pd,
    scrublet_auto_threshold,
    scrublet_scores_sim,
    scrublet_threshold_number,
):
    leiden_computed  # ran after leiden clustering

    # Raw per-cell scores, no manual numpy pre-binning -- VegaFusion pre-aggregates
    # in Python, so Altair's own bin transform can run directly on the full
    # dataset without hitting the row-embed limit.
    _obs_df = pd.DataFrame({"doublet_score": adata_flt.obs["scrublet_doublet_score"].to_numpy()})
    _sim_df = pd.DataFrame({"doublet_score": np.asarray(scrublet_scores_sim)})

    _threshold_df = pd.DataFrame({"x": [scrublet_threshold_number.value]})
    _threshold_rule = alt.Chart(_threshold_df).mark_rule(color="red", strokeDash=[4, 4]).encode(x="x:Q")

    _obs_chart = alt.layer(
        alt.Chart(_obs_df).mark_bar(color=okabe_ito_palette[0], opacity=0.8).encode(
            x=alt.X("doublet_score:Q", bin=alt.Bin(maxbins=80), title="Doublet score"),
            y=alt.Y("count():Q", scale=alt.Scale(type="log"), title="Observed cells (log)", axis=alt.Axis(grid=False)),
        ),
        _threshold_rule,
    ).properties(title="Observed transcriptomes", width=340, height=280).configure_view(strokeWidth=0)

    _sim_chart = alt.layer(
        alt.Chart(_sim_df).mark_bar(color=okabe_ito_palette[0], opacity=0.8).encode(
            x=alt.X("doublet_score:Q", bin=alt.Bin(maxbins=80), title="Doublet score"),
            y=alt.Y("count():Q", title="Simulated doublets"),
        ),
        _threshold_rule,
    ).properties(title="Simulated doublets", width=340, height=280).configure_view(strokeWidth=0)

    # --- Per-leiden-cluster singlet/doublet breakdown, as an Altair chart so it
    # can sit alongside the two histograms above.
    _is_doublet_at_threshold = pd.Series(
        adata_flt.obs["scrublet_doublet_score"].to_numpy() >= scrublet_threshold_number.value,
        index=adata_flt.obs_names,
    )
    _singlet = (~_is_doublet_at_threshold).rename("singlet")
    _doublet_crosstab = pd.crosstab(adata_flt.obs["leiden"].astype(str), _singlet)

    _doublet_crosstab_pct = _doublet_crosstab.div(_doublet_crosstab.sum(axis=1), axis=0) * 100
    _cluster_order = _doublet_crosstab_pct[True].sort_values(ascending=False).index.tolist()
    _singlet_counts_by_cluster = _doublet_crosstab[True].to_dict()

    _leiden_df = _doublet_crosstab_pct.reset_index().melt(
        id_vars="leiden", value_vars=[True, False], var_name="singlet", value_name="pct"
    )
    _leiden_df["singlet"] = _leiden_df["singlet"].astype(str)
    _leiden_df["n_singlets"] = _leiden_df["leiden"].map(_singlet_counts_by_cluster).astype(int)
    _leiden_df["label"] = np.where(
        _leiden_df["singlet"] == "True", "n=" + _leiden_df["n_singlets"].map("{:,}".format), ""
    )

    _total_singlets = int(_doublet_crosstab[True].sum())

    _leiden_bars = alt.Chart(_leiden_df).mark_bar(opacity=0.9).encode(
        y=alt.Y("leiden:N", sort=_cluster_order, title="Leiden cluster (sorted by % of singlets)"),
        x=alt.X("pct:Q", title="% of barcodes"),
        yOffset=alt.YOffset("singlet:N", sort=["True", "False"]),
        color=alt.Color(
            "singlet:N",
            sort=["True", "False"],
            scale=alt.Scale(domain=["True", "False"], range=[okabe_ito_palette[3], okabe_ito_palette[8]]),
            legend=alt.Legend(title="singlet"),
        ),
    )

    _leiden_labels = alt.Chart(_leiden_df).mark_text(align="left", dx=3, fontSize=8).encode(
        y=alt.Y("leiden:N", sort=_cluster_order),
        x=alt.X("pct:Q"),
        yOffset=alt.YOffset("singlet:N", sort=["True", "False"]),
        text="label:N",
    )

    _leiden_chart = (_leiden_bars + _leiden_labels).properties(
        title=["Scrublet singlets vs. doublets per Leiden cluster", f"n = {_total_singlets:,} singlets total"],
        width=340, height=280,
    ).configure_view(strokeWidth=0)

    _n_singlet = int((adata_flt.obs["scrublet_doublet_score"].to_numpy() < scrublet_threshold_number.value).sum())
    _pct_singlet = _n_singlet / adata_flt.n_obs

    mo.vstack([
        scrublet_threshold_number,
        mo.md(
            f"**{_n_singlet:,} of {adata_flt.n_obs:,} cells ({_pct_singlet:.1%})** "
            f"are singlets at threshold **{scrublet_threshold_number.value:.2f}** "
            f"(automatic threshold was {scrublet_auto_threshold:.2f})."
        ),
        mo.hstack([_obs_chart, _sim_chart, _leiden_chart], justify="start"),
    ])
    return


@app.cell
def qc_doublets_umap(adata_flt, okabe_ito_palette, sc):
    # Left: singlet vs. doublet call (Scrublet's automatic threshold). Right:
    # the underlying continuous doublet score, for the same UMAP layout.
    adata_flt.obs["scrublet_singlet"] = ~adata_flt.obs["scrublet_predicted_doublet"]

    sc.pl.umap(
        adata_flt,
        color=["scrublet_singlet", "scrublet_doublet_score"],
        palette={"True": okabe_ito_palette[3], "False": okabe_ito_palette[8]},
        cmap="cividis",
        ncols=2,
        size=5,
    )
    return


@app.cell(hide_code=True)
def cmo_introduction(mo):
    mo.md(r"""
    # CMO hash classification

    Classify each barcode's CMO (Cell Multiplexing Oligo) identity using an approach that mirrors `Seurat::HTODemux`. Briefly, CLR normalize the CMO counts per cell, then call a CMO "positive" if its normalized signal exceeds a per-CMO threshold (95th percentile). Barcodes positive for more than one CMO are called doublets; barcodes positive for none are negatives.

    Classification is restricted to barcodes that already pass the RNA-based QC filtering (`adata_flt`), so the CMO calls reflect real cells rather than empty droplets or debris.

    Unlike the 5-timepoints dataset, CMO ids are reused across channels for different timepoint/replicate combinations, so the CMO -> timepoint mapping is read from `metadata/10X_multi_15_timepoints_CMO_design.tsv` (filtered to this channel) rather than computed from a formula.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Load CMO counts and perform assignment
    """)
    return


@app.cell
def initialize_cmo(ad, adata_flt, ch2_data_root_path, pd, project_root):
    # path to the cmo counts
    _cmo_counts_path = ch2_data_root_path / "cmo_counts/adata.h5ad"

    adata_cmo = ad.read_h5ad(_cmo_counts_path)

    # CMO barcodes need the same per-sample suffix as adata_flt's RNA barcodes to match up
    _barcode_suffix = "_" + adata_flt.obs_names[0].split("_", 1)[1]
    adata_cmo.obs_names = adata_cmo.obs_names + _barcode_suffix

    adata_cmo.var["gene_name"] = adata_cmo.var_names

    # Filter adata_cmo data to cells in adata_flt
    _cells_in_use = adata_flt.obs_names.intersection(adata_cmo.obs_names.to_list())
    adata_cmo = adata_cmo[_cells_in_use, :]

    # Map each CMO to its timepoint for this channel. Unlike the 5-timepoints
    # dataset, CMO ids are reused across channels for different timepoints, so
    # this mapping is channel-specific and comes from the design table rather
    # than an arithmetic rule.
    _cmo_design = pd.read_csv(project_root / "metadata/10X_multi_15_timepoints_CMO_design.tsv", sep="	")
    _cmo_design_ch2 = _cmo_design[_cmo_design["Channel"] == "channel2"]
    cmo_to_timepoint = dict(zip(_cmo_design_ch2["CMO"], _cmo_design_ch2["Timepoint"]))

    # The KITE quantification panel covers all 23 CMOs regardless of channel, but
    # not every CMO is actually used in every channel (e.g. CMO23 is not part of
    # channel2's design -- see the README). Kept in the panel deliberately (not
    # dropped): a cell whose strongest signal is a CMO not in this channel's
    # design is itself a useful QC signal (contamination / index hopping), and is
    # flagged via cmo_tag_out_of_panel in assign_cmo_tag rather than hidden here.

    adata_cmo
    return adata_cmo, cmo_to_timepoint


@app.cell
def cmo_clr_normalization(adata_cmo, np):
    _cmo = adata_cmo.X
    _cmo = _cmo.toarray() if hasattr(_cmo, "toarray") else np.asarray(_cmo)

    # 1. Compute Seurat-style geometric mean per cell across all CMO features
    # Seurat sums log1p of non-zero values, divides by the total number of
    # features, then takes exp. Generalized to adata_cmo.n_vars (this
    # channel's actual CMO count -- channel2 uses only 22 of the 23 CMOs,
    # CMO23 was not used here) instead of a hardcoded feature count.
    _log_counts = np.log1p(np.where(_cmo > 0, _cmo, 0))
    _gm = np.exp(np.sum(_log_counts, axis=1) / adata_cmo.n_vars)

    # 2. Divide raw counts by the geometric mean, then apply log1p
    clr = np.log1p(_cmo / _gm[:, None])
    return (clr,)


@app.cell
def find_cmo_thresholds(clr, np):
    # For channel1, all 23 CMOs showed clean, similarly-shaped bimodal CLR
    # distributions under a per-tag GMM diagnostic (this dataset is much less
    # noisy than the mcginnis dataset that motivated the per-tag GMM approach
    # there), so a simple global CLR quantile threshold was used instead.
    # Re-check with the same per-tag diagnostic for channel2's own CMOs before
    # trusting this threshold here.
    positive_quantile_threshold = 0.95
    thresholds = np.quantile(clr, positive_quantile_threshold, axis=0)
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
    thresholds,
):
    def _run_assign_cmo_tag(adata, clr, thresholds, positive, n_pos, adata_cmo, cmo_to_timepoint):
        # Assign each barcode to its strongest-signal CMO (used for singlets and doublets alike),
        # and classify Negative/Singlet/Doublet by how many CMOs cleared their threshold.
        _top_cmo_idx = np.argmax(clr, axis=1)
        _top_cmo_tag = adata_cmo.var["gene_name"].values[_top_cmo_idx]

        adata.obs["cmo_tag_scanpy"] = np.where(n_pos == 0, "Negative", _top_cmo_tag)
        adata.obs["cmo_status_scanpy"] = np.where(
            n_pos == 0, "Negative", np.where(n_pos == 1, "Singlet", "Doublet")
        )
        adata.obs["cmo_positive_tags_scanpy"] = [
            ",".join(adata_cmo.var["gene_name"].values[row]) if row.any() else "Negative"
            for row in positive
        ]
        adata.obs["timepoint_scanpy"] = (
            adata.obs["cmo_positive_tags_scanpy"]
            .map(cmo_to_timepoint)
            .fillna(adata.obs["cmo_status_scanpy"])
        )

        # The KITE panel quantifies all 23 CMOs regardless of channel, but not
        # every CMO is used in every channel's design (e.g. CMO23 is not part of
        # channel2 -- see the README). A cell whose strongest signal is a CMO not
        # in this channel's design is itself a useful QC signal (contamination /
        # index hopping), so it is kept visible here rather than dropped from the
        # panel -- but it must not be silently treated as a normal Singlet/Doublet
        # downstream. Flagged separately; see cmo_out_of_panel_report below.
        adata.obs["cmo_tag_out_of_panel"] = (n_pos > 0) & ~adata.obs["cmo_tag_scanpy"].isin(cmo_to_timepoint)
        return True

    # Captured so downstream cells can reference this name and get a real, trackable
    # marimo dependency edge instead of just sharing a ref to "adata_flt" -- see the
    # marimo-pair race-condition note.
    cmo_assignment_computed = _run_assign_cmo_tag(adata_flt, clr, thresholds, positive, n_pos, adata_cmo, cmo_to_timepoint)
    return (cmo_assignment_computed,)


@app.cell(hide_code=True)
def cmo_out_of_panel_report(adata_flt, cmo_assignment_computed, mo):
    cmo_assignment_computed  # ran after CMO tags were assigned

    _n_out_of_panel = int(adata_flt.obs["cmo_tag_out_of_panel"].sum())
    _out_of_panel_breakdown = (
        adata_flt.obs.loc[adata_flt.obs["cmo_tag_out_of_panel"], "cmo_tag_scanpy"]
        .value_counts()
    )

    mo.md(f"""
    **{_n_out_of_panel:,} barcodes** have a strongest-signal CMO that is not part of this channel's design (not in `cmo_to_timepoint`):

    {_out_of_panel_breakdown.to_frame("n_barcodes").to_markdown() if _n_out_of_panel else "*(none)*"}

    These are flagged via `cmo_tag_out_of_panel` and excluded from `pass_strict_qc` -- not used for downstream analyses -- rather than dropped from the CMO panel itself.
    """)
    return


@app.cell
def timepoint_palette(cmo_to_timepoint, np, okabe_ito_palette, re):
    def _timepoint_sort_key(label):
        if label == "Negative":
            return (float("inf"), 0)
        _match = re.match(r"d([\d.]+)(?:_(arterial|venous))?$", label)
        _day = float(_match.group(1))
        _branch = _match.group(2)
        _branch_rank = {None: 0, "arterial": 1, "venous": 2}[_branch]
        return (_day, _branch_rank)

    timepoint_order = sorted(set(cmo_to_timepoint.values()), key=_timepoint_sort_key)

    # Pre-branch timepoints (no arterial/venous split yet) get a neutral grey ramp;
    # arterial and venous branches (from d2.5 onward) get distinct warm/cool hues
    # so the two lineages read apart at a glance -- arterial=red (oxygenated),
    # venous=blue (deoxygenated), matching the physiological convention. The grey
    # ramp is kept away from black (unlike the 0-1 full colormap range) so it
    # never gets confused with the dedicated black used for Negative below, but
    # spread wide within that range (0.1-0.75) since grey has no hue to lean on --
    # a narrow lightness range made adjacent pre-branch stages hard to tell apart.
    _pre_branch = [t for t in timepoint_order if "_" not in t]
    _arterial = [t for t in timepoint_order if t.endswith("_arterial")]
    _venous = [t for t in timepoint_order if t.endswith("_venous")]

    from matplotlib.colors import to_hex
    from matplotlib import colormaps as _colormaps

    def _ramp(cmap_name, n, lo=0.35, hi=0.85):
        _cmap = _colormaps[cmap_name]
        return [to_hex(_cmap(_t)) for _t in np.linspace(lo, hi, n)]

    ec_diff_palette = {
        **dict(zip(_pre_branch, _ramp("Greys", len(_pre_branch), lo=0.1, hi=0.75))),
        **dict(zip(_arterial, _ramp("Reds", len(_arterial)))),
        **dict(zip(_venous, _ramp("Blues", len(_venous)))),
        "Negative": okabe_ito_palette[0],
    }
    timepoint_order = timepoint_order + ["Negative"]

    ec_diff_palette
    return ec_diff_palette, timepoint_order


@app.cell(hide_code=True)
def cmo_scatter_multiselect_ui(mo, timepoint_order):
    # Filters which timepoints appear in the threshold-vs-assignment scatter below.
    # Displayed together with the plots (in cmo_threshold_vs_assignment), not here
    # -- this cell only defines it, since its .value cannot be read in the same
    # cell that creates it.
    _timepoints_only = [t for t in timepoint_order if t != "Negative"]
    timepoint_multiselect = mo.ui.multiselect(
        options=_timepoints_only,
        value=_timepoints_only,
        label="Timepoints to show in scatter",
    )
    return (timepoint_multiselect,)


@app.cell(hide_code=True)
def cmo_threshold_vs_assignment(
    adata_cmo,
    adata_flt,
    alt,
    cmo_assignment_computed,
    cmo_to_timepoint,
    ec_diff_palette,
    mo,
    okabe_ito_palette,
    pd,
    thresholds,
    timepoint_multiselect,
    timepoint_order,
):
    cmo_assignment_computed  # ran after CMO tags were assigned

    # --- Per-CMO threshold vs. assignment count -------------------------------
    # A CMO with a low threshold AND a low assigned count suggests weak/
    # inefficient staining (background and positive signal are both compressed),
    # rather than a contamination or misclassification issue.
    _cmo_summary = pd.DataFrame({
        "threshold": thresholds,
        "n_assigned": adata_flt.obs["cmo_tag_scanpy"].value_counts().reindex(adata_cmo.var["gene_name"]).to_numpy(),
    }, index=adata_cmo.var["gene_name"])
    _cmo_summary["timepoint"] = _cmo_summary.index.map(cmo_to_timepoint)
    _cmo_summary = _cmo_summary.reset_index().rename(columns={"gene_name": "cmo"})
    _cmo_summary_filtered = _cmo_summary[_cmo_summary["timepoint"].isin(timepoint_multiselect.value)]

    _present_timepoints = [t for t in ec_diff_palette if t in _cmo_summary["timepoint"].unique()]
    _timepoint_scale = alt.Scale(domain=_present_timepoints, range=[ec_diff_palette[t] for t in _present_timepoints])

    # Fixed axis domains from the FULL (unfiltered) data, with a little padding, so
    # the scatter's scale does not rescale as the timepoint selection changes.
    _x_pad = (_cmo_summary["threshold"].max() - _cmo_summary["threshold"].min()) * 0.1
    _y_pad = (_cmo_summary["n_assigned"].max() - _cmo_summary["n_assigned"].min()) * 0.1
    _x_domain = [_cmo_summary["threshold"].min() - _x_pad, _cmo_summary["threshold"].max() + _x_pad]
    _y_domain = [_cmo_summary["n_assigned"].min() - _y_pad, _cmo_summary["n_assigned"].max() + _y_pad]

    _scatter_points = alt.Chart(_cmo_summary_filtered).mark_circle(size=100, opacity=0.9, stroke=okabe_ito_palette[0], strokeWidth=0.5).encode(
        x=alt.X("threshold:Q", title="CLR detection threshold (95th percentile)", scale=alt.Scale(domain=_x_domain)),
        y=alt.Y("n_assigned:Q", title="Number of barcodes assigned to this CMO", scale=alt.Scale(domain=_y_domain)),
        color=alt.Color("timepoint:N", scale=_timepoint_scale),
        tooltip=["cmo", "threshold", "n_assigned", "timepoint"],
    )
    _scatter_labels = alt.Chart(_cmo_summary_filtered).mark_text(align="left", dx=5, dy=-5, fontSize=8).encode(
        x=alt.X("threshold:Q", scale=alt.Scale(domain=_x_domain)),
        y=alt.Y("n_assigned:Q", scale=alt.Scale(domain=_y_domain)),
        text="cmo:N",
    )
    _scatter_chart = (_scatter_points + _scatter_labels).properties(
        title="Per-CMO threshold vs. assignment count", width=380, height=340,
    ).configure_view(strokeWidth=0)

    # --- Number of barcodes assigned per timepoint (CMOs summed within each timepoint block) ---
    _hash_to_timepoint = adata_flt.obs["cmo_tag_scanpy"].map(cmo_to_timepoint).fillna("Negative")
    _timepoint_counts = _hash_to_timepoint.value_counts().reindex(timepoint_order).reset_index()
    _timepoint_counts.columns = ["timepoint", "n_barcodes"]
    _total_cells = int(_timepoint_counts["n_barcodes"].sum())

    _bar_colors = {t: ec_diff_palette[t] for t in timepoint_order}
    _bar_chart = alt.Chart(_timepoint_counts).mark_bar().encode(
        y=alt.Y("timepoint:N", sort=timepoint_order, title="Timepoint"),
        x=alt.X("n_barcodes:Q", title="Number of barcodes"),
        color=alt.Color(
            "timepoint:N",
            scale=alt.Scale(domain=list(_bar_colors.keys()), range=list(_bar_colors.values())),
            legend=None,
        ),
    )
    _bar_labels = alt.Chart(_timepoint_counts).mark_text(align="left", dx=3, fontSize=9).encode(
        y=alt.Y("timepoint:N", sort=timepoint_order),
        x=alt.X("n_barcodes:Q"),
        text=alt.Text("n_barcodes:Q", format=","),
    )
    _bar_final = (_bar_chart + _bar_labels).properties(
        title=["Barcodes assigned per timepoint", f"n = {_total_cells:,} cells total"],
        width=380, height=340,
    ).configure_view(strokeWidth=0)

    # --- Number of barcodes per CMO-hashing status (Singlet/Doublet/Negative) ---
    _status_order = ["Singlet", "Doublet", "Negative"]
    _status_colors = {"Singlet": okabe_ito_palette[3], "Doublet": okabe_ito_palette[8], "Negative": okabe_ito_palette[0]}
    _status_counts = adata_flt.obs["cmo_status_scanpy"].value_counts().reindex(_status_order).reset_index()
    _status_counts.columns = ["cmo_status", "n_barcodes"]
    _total_status = int(_status_counts["n_barcodes"].sum())

    _status_chart = alt.Chart(_status_counts).mark_bar().encode(
        y=alt.Y("cmo_status:N", sort=_status_order, title="CMO hashing status"),
        x=alt.X("n_barcodes:Q", title="Number of barcodes"),
        color=alt.Color(
            "cmo_status:N",
            scale=alt.Scale(domain=_status_order, range=[_status_colors[s] for s in _status_order]),
            legend=None,
        ),
    )
    _status_labels = alt.Chart(_status_counts).mark_text(align="left", dx=3, fontSize=9).encode(
        y=alt.Y("cmo_status:N", sort=_status_order),
        x=alt.X("n_barcodes:Q"),
        text=alt.Text("n_barcodes:Q", format=","),
    )
    _status_final = (_status_chart + _status_labels).properties(
        title=["Barcodes per CMO hashing status", f"n = {_total_status:,} cells total"],
        width=380, height=340,
    ).configure_view(strokeWidth=0)

    mo.vstack([
        timepoint_multiselect,
        mo.hstack([_scatter_chart, _bar_final, _status_final], justify="start"),
    ])
    return


@app.cell(hide_code=True)
def cmo_doublet_pair_confounding(
    adata_cmo,
    alt,
    cmo_to_timepoint,
    mo,
    n_pos,
    np,
    okabe_ito_palette,
    pd,
    positive,
    timepoint_order,
):
    # Are CMO doublets consistently the same two CMOs getting confounded together
    # (e.g. adjacent wells cross-contaminating), or is which two CMOs co-occur
    # essentially random given how often each CMO individually gets falsely
    # triggered? Restricted to doublets with exactly two positive CMOs (n_pos == 2),
    # since a clean pairwise co-occurrence test needs exactly one pair per barcode --
    # barcodes with 3+ positive CMOs do not have a single well-defined pair.
    from statsmodels.stats.multitest import multipletests

    _cmo_names = adata_cmo.var["gene_name"].to_numpy()
    _doublet_mask = (n_pos == 2)
    _positive_doublets = positive[_doublet_mask]
    _pair_idx_all = np.array([np.flatnonzero(row) for row in _positive_doublets])

    # Exclude pairs involving a CMO not in this channel's design (e.g. CMO23 in
    # channel2) -- these have no defined timepoint to compare against, and are
    # already flagged/excluded from pass_strict_qc via cmo_tag_out_of_panel (see
    # cmo_out_of_panel_report above), so they are dropped from this pairwise
    # analysis too rather than crashing on the timepoint lookup below.
    _in_panel_pair_mask = np.array([
        _cmo_names[_p[0]] in cmo_to_timepoint and _cmo_names[_p[1]] in cmo_to_timepoint
        for _p in _pair_idx_all
    ]) if len(_pair_idx_all) else np.array([], dtype=bool)
    _n_out_of_panel_pairs = int((~_in_panel_pair_mask).sum())
    _pair_idx = _pair_idx_all[_in_panel_pair_mask]
    _n_doublets = len(_pair_idx)

    _n_all_doublets = int((n_pos >= 2).sum())
    _n_multi_excluded = _n_all_doublets - len(_pair_idx_all)
    _pct_multi_excluded = _n_multi_excluded / _n_all_doublets

    _observed_pairs = pd.Series(
        [tuple(sorted(_cmo_names[p])) for p in _pair_idx]
    ).value_counts()

    # Null model: pool every CMO "hit" event across all doublets (2 per doublet),
    # then randomly re-pair them. This preserves each CMO's overall rate of being
    # falsely triggered, but destroys any preferential pairing between specific CMOs.
    _hit_pool = _cmo_names[_pair_idx.ravel()]
    _rng = np.random.default_rng(0)
    _n_perm = 1000

    _perm_counts = {}
    _top5_share_null = np.empty(_n_perm)
    for _i in range(_n_perm):
        _shuffled = _rng.permutation(_hit_pool)
        _perm_pairs = pd.Series(
            [tuple(sorted(pair)) for pair in _shuffled.reshape(-1, 2)]
        ).value_counts()
        _top5_share_null[_i] = _perm_pairs.head(5).sum() / _n_doublets
        for _pair, _count in _perm_pairs.items():
            _perm_counts.setdefault(_pair, []).append(_count)

    _rows = []
    for _pair, _obs_count in _observed_pairs.items():
        _null = np.array(_perm_counts.get(_pair, []))
        _null = np.pad(_null, (0, _n_perm - len(_null)))
        _p = (_null >= _obs_count).mean()
        _rows.append({
            "cmo_pair": " + ".join(_pair),
            "observed": int(_obs_count),
            "null_mean": round(_null.mean(), 2),
            "null_95th_pct": np.quantile(_null, 0.95),
            "p_value": _p,
        })

    cmo_pair_confounding_table = pd.DataFrame(_rows).sort_values("observed", ascending=False).reset_index(drop=True)

    # Multiple-testing correction across all distinct pairs tested (Benjamini-Hochberg
    # FDR, less conservative than Bonferroni and standard for many simultaneous tests).
    # Permutation p-values are floored at 1/n_perm so BH does not choke on exact zeros.
    _p_floored = cmo_pair_confounding_table["p_value"].clip(lower=1 / _n_perm)
    _reject_bh, _p_adj_bh, _, _ = multipletests(_p_floored, alpha=0.20, method="fdr_bh")
    cmo_pair_confounding_table["p_adj_bh"] = _p_adj_bh
    cmo_pair_confounding_table["significant_bh"] = _reject_bh

    # Global summary: do the top pairs account for more of the doublets than chance
    # would predict, given the same per-CMO false-trigger rates?
    _top5_share_obs = _observed_pairs.head(5).sum() / _n_doublets
    _p_top5 = (_top5_share_null >= _top5_share_obs).mean()

    _n_significant_raw = int((cmo_pair_confounding_table["p_value"] < 0.05).sum())
    _n_significant_bh = int(cmo_pair_confounding_table["significant_bh"].sum())

    _summary_md = mo.md(f"""
    **CMO doublet pairing: {_n_doublets:,} doublets with exactly two positive CMOs,
    {_observed_pairs.shape[0]} distinct CMO pairs observed.**

    {_n_multi_excluded:,} of {_n_all_doublets:,} doublets ({_pct_multi_excluded:.1%}) have
    3 or more positive CMOs and are excluded from this pairwise analysis, since they
    do not reduce to a single CMO pair. Of the remaining exactly-two-CMO doublets,
    {_n_out_of_panel_pairs:,} involve a CMO not in this channel's design (see
    `cmo_out_of_panel_report` above) and are also excluded here.

    The top 5 CMO pairs account for **{_top5_share_obs:.1%}** of all exactly-two-CMO doublets, vs. a
    **{_top5_share_null.mean():.1%}** null mean (1,000-permutation shuffle of which CMO-hit
    events get paired together, preserving each CMO's own false-trigger rate) -- permutation
    p = {_p_top5:.3f}.

    Across the {cmo_pair_confounding_table.shape[0]} distinct pairs tested, {_n_significant_raw}
    have an uncorrected permutation p < 0.05; after Benjamini-Hochberg FDR correction
    (20% FDR), **{_n_significant_bh}** remain significant.
    """)

    _table_md = mo.md(cmo_pair_confounding_table.head(10).to_markdown(index=False))

    # --- Same pairs, grouped by timepoint: does confounding happen mostly within
    # the same timepoint (adjacent CMOs in the same block) or across timepoints? ---
    _pair_timepoints = np.array([
        [cmo_to_timepoint[_cmo_names[p[0]]], cmo_to_timepoint[_cmo_names[p[1]]]]
        for p in _pair_idx
    ])

    # This channel's own timepoint labels, not the 5-timepoints dataset's fixed
    # d0-d4 -- see the timepoint_palette cell.
    _timepoint_order = [t for t in timepoint_order if t != "Negative"]
    _confusion = pd.DataFrame(0, index=_timepoint_order, columns=_timepoint_order)
    for _ta, _tb in _pair_timepoints:
        _confusion.loc[_ta, _tb] += 1
        if _ta != _tb:
            _confusion.loc[_tb, _ta] += 1

    _confusion_long = _confusion.reset_index().rename(columns={"index": "timepoint_a"}).melt(
        id_vars="timepoint_a", var_name="timepoint_b", value_name="count"
    )

    _confusion_heatmap = alt.Chart(_confusion_long).mark_rect().encode(
        x=alt.X("timepoint_b:N", sort=_timepoint_order, title="Timepoint"),
        y=alt.Y("timepoint_a:N", sort=_timepoint_order, title="Timepoint"),
        color=alt.Color("count:Q", scale=alt.Scale(scheme="cividis"), title="Doublets"),
        tooltip=["timepoint_a", "timepoint_b", "count"],
    )
    _confusion_labels = alt.Chart(_confusion_long).mark_text(fontSize=11).encode(
        x=alt.X("timepoint_b:N", sort=_timepoint_order),
        y=alt.Y("timepoint_a:N", sort=_timepoint_order),
        text=alt.Text("count:Q", format=","),
        color=alt.condition(alt.datum.count > _confusion_long["count"].max() / 2, alt.value(okabe_ito_palette[0]), alt.value("white")),
    )
    _confusion_chart = (_confusion_heatmap + _confusion_labels).properties(
        title="Doublet pairs by timepoint (raw count)", width=280, height=280,
    ).configure_view(strokeWidth=0)

    # Same matrix, normalized by how many possible CMO pairs exist in each cell --
    # timepoints with only 1 CMO in this channel have zero within-timepoint pairs
    # possible, and same-timepoint vs. cross-timepoint cells can have very
    # different numbers of possible CMO pairs, so raw counts alone make
    # higher-CMO-count timepoints look inflated just from having more combinations.
    # Restricted to in-panel CMOs (cmo_to_timepoint keys), since a CMO not used in
    # this channel's design has no timepoint to count towards.
    _cmos_per_timepoint = {
        t: [c for c in _cmo_names if c in cmo_to_timepoint and cmo_to_timepoint[c] == t] for t in _timepoint_order
    }
    _n_possible = pd.DataFrame(0.0, index=_timepoint_order, columns=_timepoint_order)
    for _ta in _timepoint_order:
        for _tb in _timepoint_order:
            _na, _nb = len(_cmos_per_timepoint[_ta]), len(_cmos_per_timepoint[_tb])
            _n_possible.loc[_ta, _tb] = (_na * (_na - 1) / 2) if _ta == _tb else (_na * _nb)

    _rate = _confusion / _n_possible.replace(0, np.nan)
    _rate_long = _rate.reset_index().rename(columns={"index": "timepoint_a"}).melt(
        id_vars="timepoint_a", var_name="timepoint_b", value_name="rate"
    )

    _rate_heatmap = alt.Chart(_rate_long).mark_rect().encode(
        x=alt.X("timepoint_b:N", sort=_timepoint_order, title="Timepoint"),
        y=alt.Y("timepoint_a:N", sort=_timepoint_order, title="Timepoint"),
        color=alt.Color("rate:Q", scale=alt.Scale(scheme="cividis"), title="Doublets / possible pair"),
        tooltip=["timepoint_a", "timepoint_b", alt.Tooltip("rate:Q", format=".1f")],
    )
    _rate_labels = alt.Chart(_rate_long).mark_text(fontSize=11).encode(
        x=alt.X("timepoint_b:N", sort=_timepoint_order),
        y=alt.Y("timepoint_a:N", sort=_timepoint_order),
        text=alt.Text("rate:Q", format=".1f"),
        color=alt.condition(alt.datum.rate > _rate_long["rate"].max() / 2, alt.value(okabe_ito_palette[0]), alt.value("white")),
    )
    _rate_chart = (_rate_heatmap + _rate_labels).properties(
        title="Doublet pairs by timepoint (per possible CMO pair)", width=280, height=280,
    ).configure_view(strokeWidth=0)

    def _tight_row(*items):
        # mo.hstack with widths=None adds no wrapper/flex styling around children,
        # so a block-level markdown table (which wants to be as wide as its
        # container) and the chart just fill the row between them -- no slack left
        # for justify-content to redistribute. Build the flex row by hand instead.
        _items_html = "".join(
            f'<div style="flex: 0 0 auto;">{mo.as_html(it).text}</div>' for it in items
        )
        return mo.Html(f'<div style="display:flex; justify-content:flex-start; gap:1rem;">{_items_html}</div>')

    mo.vstack([
        _summary_md,
        _tight_row(_table_md, _confusion_chart, _rate_chart),
    ])
    return


@app.cell(hide_code=True)
def plot_cmo_assignment_counts(
    adata_cmo,
    adata_flt,
    cmo_assignment_computed,
    cmo_to_timepoint,
    ec_diff_palette,
    plt,
):
    cmo_assignment_computed  # ran after CMO tags were assigned

    # Number of barcodes assigned (hash_ID) to each CMO
    _cmo_order = list(adata_cmo.var["gene_name"]) + ["Negative"]
    _counts = adata_flt.obs["cmo_tag_scanpy"].value_counts().reindex(_cmo_order)
    _colors = [ec_diff_palette[cmo_to_timepoint.get(c, "Negative")] for c in _counts.index]

    plt.figure(figsize=(10, 5))
    plt.bar(_counts.index, _counts.values, color=_colors)
    plt.xlabel("CMO")
    plt.ylabel("Number of barcodes")
    plt.title("Barcodes assigned per CMO")
    plt.xticks(rotation=90)
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
def cmo_doublets_per_leiden_pct(
    adata_flt,
    alt,
    cmo_assignment_computed,
    leiden_computed,
    mo,
    np,
    okabe_ito_palette,
    pd,
):
    cmo_assignment_computed  # ran after CMO tags were assigned
    leiden_computed  # ran after leiden clustering

    # Same cluster ordering (by Scrublet singlet %) applied to both panels for direct comparison
    _singlet_scrublet = (~adata_flt.obs["scrublet_predicted_doublet"]).rename("singlet")
    _scrublet_crosstab = pd.crosstab(adata_flt.obs["leiden"].astype(str), _singlet_scrublet)
    _scrublet_pct = _scrublet_crosstab.div(_scrublet_crosstab.sum(axis=1), axis=0) * 100
    _cluster_order = _scrublet_pct[True].sort_values(ascending=False).index.tolist()

    # --- Scrublet panel: binary singlet vs. doublet ---
    _long_scrublet = _scrublet_pct.reset_index().melt(id_vars="leiden", value_vars=[True, False], var_name="singlet", value_name="pct")
    _long_scrublet["singlet"] = _long_scrublet["singlet"].astype(str)
    _long_scrublet["opacity"] = 1.0
    _scrublet_counts_by_cluster = _scrublet_crosstab[True].to_dict()
    _long_scrublet["n"] = _long_scrublet["leiden"].map(_scrublet_counts_by_cluster).astype(int)
    _long_scrublet["label"] = np.where(_long_scrublet["singlet"] == "True", "n=" + _long_scrublet["n"].map("{:,}".format), "")

    # Fake, invisible third slot so bar thickness/spacing lines up row-for-row with
    # the 3-category CMO panel (which has one more group per cluster).
    _pad_rows = pd.DataFrame({
        "leiden": _cluster_order, "singlet": "pad", "pct": 0.0, "n": 0, "label": "", "opacity": 0.0,
    })
    _long_scrublet = pd.concat([_long_scrublet, _pad_rows], ignore_index=True)

    _scrublet_bars = alt.Chart(_long_scrublet).mark_bar().encode(
        y=alt.Y("leiden:N", sort=_cluster_order, title="Leiden cluster (sorted by Scrublet singlet %)"),
        x=alt.X("pct:Q", title="% of barcodes", scale=alt.Scale(domain=[0, 100])),
        yOffset=alt.YOffset("singlet:N", sort=["True", "False", "pad"]),
        color=alt.Color(
            "singlet:N",
            sort=["True", "False"],
            scale=alt.Scale(domain=["True", "False"], range=[okabe_ito_palette[3], okabe_ito_palette[8]]),
            legend=alt.Legend(title="singlet"),
        ),
        opacity=alt.Opacity("opacity:Q", legend=None),
    )
    _scrublet_labels = alt.Chart(_long_scrublet).mark_text(align="left", dx=3, fontSize=8).encode(
        y=alt.Y("leiden:N", sort=_cluster_order),
        x=alt.X("pct:Q"),
        yOffset=alt.YOffset("singlet:N", sort=["True", "False", "pad"]),
        text="label:N",
    )
    _scrublet_chart = (_scrublet_bars + _scrublet_labels).properties(
        title="Scrublet: singlet vs. doublet", width=380, height=340,
    ).configure_view(strokeWidth=0)

    # --- CMO hashing panel: full 3-way status (Singlet / Doublet / Negative) --
    # "not singlet" is not one thing for CMO hashing -- it is either a Doublet call
    # or a Negative (no CMO cleared threshold), which are different failure modes.
    _status_order = ["Singlet", "Doublet", "Negative"]
    _status_colors = {"Singlet": okabe_ito_palette[3], "Doublet": okabe_ito_palette[8], "Negative": okabe_ito_palette[0]}
    _cmo_crosstab = pd.crosstab(adata_flt.obs["leiden"].astype(str), adata_flt.obs["cmo_status_scanpy"])
    _cmo_crosstab = _cmo_crosstab.reindex(columns=_status_order, fill_value=0)
    _cmo_pct = _cmo_crosstab.div(_cmo_crosstab.sum(axis=1), axis=0) * 100

    _long_cmo = _cmo_pct.reset_index().melt(id_vars="leiden", value_vars=_status_order, var_name="cmo_status", value_name="pct")
    _cmo_counts_by_cluster = _cmo_crosstab["Singlet"].to_dict()
    _long_cmo["n"] = _long_cmo["leiden"].map(_cmo_counts_by_cluster).astype(int)
    _long_cmo["label"] = np.where(_long_cmo["cmo_status"] == "Singlet", "n=" + _long_cmo["n"].map("{:,}".format), "")

    _cmo_bars = alt.Chart(_long_cmo).mark_bar(opacity=0.9).encode(
        y=alt.Y("leiden:N", sort=_cluster_order, title=None),
        x=alt.X("pct:Q", title="% of barcodes", scale=alt.Scale(domain=[0, 100])),
        yOffset=alt.YOffset("cmo_status:N", sort=_status_order),
        color=alt.Color(
            "cmo_status:N",
            sort=_status_order,
            scale=alt.Scale(domain=_status_order, range=[_status_colors[s] for s in _status_order]),
            legend=alt.Legend(title="CMO status"),
        ),
    )
    _cmo_labels = alt.Chart(_long_cmo).mark_text(align="left", dx=3, fontSize=8).encode(
        y=alt.Y("leiden:N", sort=_cluster_order),
        x=alt.X("pct:Q"),
        yOffset=alt.YOffset("cmo_status:N", sort=_status_order),
        text="label:N",
    )
    _cmo_chart = (_cmo_bars + _cmo_labels).properties(
        title="CMO hashing: singlet / doublet / negative", width=380, height=340,
    ).configure_view(strokeWidth=0)

    mo.vstack([
        mo.md("*Clusters ordered by Scrublet singlet % (shared across both panels).*"),
        mo.hstack([_scrublet_chart, _cmo_chart], justify="start"),
    ])
    return


@app.cell(hide_code=True)
def scrublet_score_per_cmo(
    adata_flt,
    cmo_assignment_computed,
    okabe_ito_palette,
    plt,
    sc,
    sns,
):
    cmo_assignment_computed  # ran after CMO tags were assigned

    # Approximate threshold as the midpoint between the highest non-doublet score
    # and the lowest predicted-doublet score (Scrublet does not expose it directly here).
    _scrublet_threshold = (
        adata_flt.obs.loc[adata_flt.obs["scrublet_predicted_doublet"], "scrublet_doublet_score"].min()
        + adata_flt.obs.loc[~adata_flt.obs["scrublet_predicted_doublet"], "scrublet_doublet_score"].max()
    ) / 2

    _status_order = ["Singlet", "Doublet", "Negative"]
    _status_palette = [okabe_ito_palette[3], okabe_ito_palette[8], okabe_ito_palette[0]]

    _fig, (_ax1, _ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    # Left: violin (distribution shape). Points are off here -- scanpy's stripplot
    # dots are always black, invisible against the black Negative violin -- so the
    # actual per-cell points are shown separately on the right instead.
    sc.pl.violin(
        adata_flt,
        keys="scrublet_doublet_score",
        groupby="cmo_status_scanpy",
        order=_status_order,
        ylabel="Scrublet doublet score",
        palette=_status_palette,
        stripplot=False,
        show=False,
        ax=_ax1,
    )
    _ax1.axhline(_scrublet_threshold, color="red", linestyle="--", label=f"Scrublet threshold ({_scrublet_threshold:.2f})")
    _ax1.legend()
    _ax1.set_title("Distribution")

    # Right: just the points, jittered, colored per status (not the scanpy default
    # black) so the Negative group's points are actually visible.
    sns.stripplot(
        data=adata_flt.obs,
        x="cmo_status_scanpy",
        y="scrublet_doublet_score",
        order=_status_order,
        hue="cmo_status_scanpy",
        hue_order=_status_order,
        palette=_status_palette,
        size=2,
        alpha=0.5,
        jitter=0.3,
        legend=False,
        ax=_ax2,
    )
    _ax2.axhline(_scrublet_threshold, color="red", linestyle="--")
    _ax2.set_ylabel("")
    _ax2.set_xlabel("")
    _ax2.set_title("Individual cells")

    plt.tight_layout()
    _fig
    return


@app.cell
def cmo_status_multiselect_ui(mo):
    # Filters which CMO-hashing status categories appear in the two UMAPs below.
    # Displayed together with the plots (in umap_by_cmo_hashing_altair), not here --
    # this cell only defines it, since its .value cannot be read in the same cell
    # that creates it.
    cmo_status_multiselect = mo.ui.multiselect(
        options=["Singlet", "Doublet", "Negative"],
        value=["Singlet", "Doublet", "Negative"],
        label="CMO hashing status to show",
    )
    return (cmo_status_multiselect,)


@app.cell(hide_code=True)
def umap_by_cmo_hashing_altair(
    adata_flt,
    alt,
    cmo_assignment_computed,
    cmo_status_multiselect,
    mo,
    np,
    okabe_ito_palette,
    pd,
    umap_computed,
):
    cmo_assignment_computed  # ran after CMO tags were assigned
    umap_computed  # ran after UMAP was computed

    _status_order = ["Singlet", "Doublet", "Negative"]
    _status_colors = {"Singlet": okabe_ito_palette[3], "Doublet": okabe_ito_palette[8], "Negative": okabe_ito_palette[0]}

    _mask = adata_flt.obs["cmo_status_scanpy"].isin(cmo_status_multiselect.value).to_numpy()

    _umap_df = pd.DataFrame(adata_flt.obsm["X_umap"][_mask], columns=["UMAP1", "UMAP2"])
    _umap_df["cmo_status"] = adata_flt.obs["cmo_status_scanpy"].to_numpy()[_mask]
    _umap_df["scrublet_status"] = np.where(
        adata_flt.obs["scrublet_predicted_doublet"].to_numpy()[_mask], "Doublet", "Singlet"
    )

    # Subsample for the plot itself -- two full-size raw scatter specs in one cell
    # output is too large for marimo to render (no aggregating transform here for
    # VegaFusion to pre-reduce, unlike the binned histograms). Both panels share
    # the same subsampled points so they stay directly comparable.
    _max_points = 8000
    _n_before_subsample = len(_umap_df)
    if _n_before_subsample > _max_points:
        _umap_df = _umap_df.sample(n=_max_points, random_state=0)

    # Fixed axis domains from the FULL (unfiltered) UMAP, with a little padding, so
    # the scale does not rescale as the status selection changes.
    _full_umap = adata_flt.obsm["X_umap"]
    _x_pad = (_full_umap[:, 0].max() - _full_umap[:, 0].min()) * 0.05
    _y_pad = (_full_umap[:, 1].max() - _full_umap[:, 1].min()) * 0.05
    _x_domain = [_full_umap[:, 0].min() - _x_pad, _full_umap[:, 0].max() + _x_pad]
    _y_domain = [_full_umap[:, 1].min() - _y_pad, _full_umap[:, 1].max() + _y_pad]

    _scrublet_umap = alt.Chart(_umap_df).mark_circle(size=10, opacity=0.6).encode(
        x=alt.X("UMAP1:Q", scale=alt.Scale(domain=_x_domain)),
        y=alt.Y("UMAP2:Q", scale=alt.Scale(domain=_y_domain)),
        color=alt.Color(
            "scrublet_status:N",
            scale=alt.Scale(domain=["Singlet", "Doublet"], range=[_status_colors["Singlet"], _status_colors["Doublet"]]),
            legend=alt.Legend(title="Scrublet"),
        ),
        tooltip=["scrublet_status"],
    ).properties(title="UMAP colored by Scrublet", width=380, height=380).configure_view(strokeWidth=0)

    _cmo_umap = alt.Chart(_umap_df).mark_circle(size=10, opacity=0.6).encode(
        x=alt.X("UMAP1:Q", scale=alt.Scale(domain=_x_domain)),
        y=alt.Y("UMAP2:Q", scale=alt.Scale(domain=_y_domain)),
        color=alt.Color(
            "cmo_status:N",
            scale=alt.Scale(domain=_status_order, range=[_status_colors[s] for s in _status_order]),
            legend=alt.Legend(title="CMO hashing"),
        ),
        tooltip=["cmo_status"],
    ).properties(title="UMAP colored by CMO hashing", width=380, height=380).configure_view(strokeWidth=0)

    def _tight_row(*items):
        # mo.hstack with widths=None adds no wrapper/flex styling around children,
        # so block-level chart divs just fill the row (no slack left for
        # justify-content to redistribute). Build the flex row by hand instead.
        _items_html = "".join(
            f'<div style="flex: 0 0 auto;">{mo.as_html(it).text}</div>' for it in items
        )
        return mo.Html(f'<div style="display:flex; justify-content:flex-start; gap:1rem;">{_items_html}</div>')

    _caption = mo.md(
        f"*Showing {len(_umap_df):,} of {_n_before_subsample:,} matching cells "
        f"(subsampled to keep the plot size manageable).*"
        if _n_before_subsample > _max_points else ""
    )

    mo.vstack([
        cmo_status_multiselect,
        _caption,
        _tight_row(_scrublet_umap, _cmo_umap),
    ])
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
    ## Summary of quality control metrics per cluster
    """)
    return


@app.cell(hide_code=True)
def leiden_doublet_overview(
    adata_flt,
    cmo_assignment_computed,
    leiden_computed,
):
    cmo_assignment_computed  # ran after CMO tags were assigned
    leiden_computed  # ran after leiden clustering

    # Per-cluster overview: does the CMO-hashing doublet rate track the Scrublet doublet rate?
    leiden_doublet_summary = adata_flt.obs.groupby("leiden", observed=True).agg(
        n_cells=("leiden", "size"),
        pct_cmo_singlet=("cmo_status_scanpy", lambda s: (s == "Singlet").mean() * 100),
        pct_cmo_doublet=("cmo_status_scanpy", lambda s: (s == "Doublet").mean() * 100),
        pct_cmo_negative=("cmo_status_scanpy", lambda s: (s == "Negative").mean() * 100),
        pct_scrublet_doublet=("scrublet_predicted_doublet", lambda s: s.mean() * 100),
        median_mito=("pct_counts_mt", "median"),
        median_counts=("total_counts", "median"),
    ).round(1)

    leiden_doublet_summary[
        ["n_cells", "pct_scrublet_doublet", "pct_cmo_doublet", "median_mito"]
    ].sort_values("pct_cmo_doublet", ascending=False)
    return (leiden_doublet_summary,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Filtering strategy: conclusions
    """)
    return


@app.cell
def strict_qc_cascade_intro(mo):
    mo.md(r"""
    ## Consolidated strict QC decision

    Three independent exclusion criteria combine into `pass_strict_qc`: doublet-dominated clusters, Scrublet-predicted doublets, and CMO hashing status != Singlet (excludes both CMO Doublet and CMO Negative calls -- there is no dedicated negative-tag rescue analysis for this dataset, so Negatives are excluded outright rather than kept unfiltered). The mito filtering itself was already applied upstream, pre-clustering, in `mask_filter_cells` (median + 2.5 MADs); no additional cluster- or cell-level mito exclusion is applied here. The breakdown below applies the three criteria as a cascade, ordered from least to most restrictive (by how many cells each filter alone would keep), the same style used for the lenient QC filter earlier in `mask_filter_cells`.
    """)
    return


@app.cell(hide_code=True)
def strict_qc_cascade(
    adata_flt,
    cmo_assignment_computed,
    leiden_doublet_summary,
    mo,
    np,
    pd,
):
    cmo_assignment_computed  # ran after CMO tags were assigned

    def _run_strict_qc(adata, doublet_dominated_clusters):
        _is_doublet_cluster = adata.obs["leiden"].isin(doublet_dominated_clusters)
        _is_scrublet_doublet = adata.obs["scrublet_predicted_doublet"]
        _is_cmo_doublet = adata.obs["cmo_status_scanpy"] == "Doublet"
        _is_cmo_negative = adata.obs["cmo_status_scanpy"] == "Negative"

        adata.obs["qc_exclude_reason"] = np.select(
            [_is_doublet_cluster, _is_scrublet_doublet, _is_cmo_doublet, _is_cmo_negative],
            ["doublet_cluster", "scrublet_doublet", "cmo_doublet", "cmo_negative"],
            default="",
        )
        adata.obs["pass_not_doublet_cluster"] = ~_is_doublet_cluster
        adata.obs["pass_not_scrublet_doublet"] = ~_is_scrublet_doublet
        # Not rescuing Negatives (no dedicated rescue analysis for this dataset), so
        # the CMO step must require Singlet outright -- "!= Doublet" alone would
        # silently let Negatives through. Also excludes cmo_tag_out_of_panel cells
        # (top signal is a CMO not in this channel's design, e.g. CMO23 in
        # channel2) -- see cmo_out_of_panel_report above.
        adata.obs["pass_cmo_singlet"] = (adata.obs["cmo_status_scanpy"] == "Singlet") & ~adata.obs["cmo_tag_out_of_panel"]
        return True

    # I am being more strict here
    _doublet_cluster_mask = leiden_doublet_summary["pct_cmo_doublet"].gt(40) | leiden_doublet_summary["pct_scrublet_doublet"].gt(40)
    doublet_dominated_clusters = leiden_doublet_summary.index[_doublet_cluster_mask].tolist()

    # Picked programmatically (not hardcoded) so the write-up and charts referencing
    # "the" doublet cluster / "the" well-behaved cluster stay correct even if Leiden
    # renumbers clusters on a future re-run.
    representative_doublet_cluster = leiden_doublet_summary.loc[doublet_dominated_clusters, "pct_cmo_doublet"].idxmax()
    representative_well_behaved_cluster = leiden_doublet_summary["pct_cmo_doublet"].idxmin()

    _ok = _run_strict_qc(adata_flt, doublet_dominated_clusters)

    # Cascade breakdown, like mask_filter_cells, steps ordered from least to most
    # restrictive (by how many cells each filter alone would keep) so the table
    # reads as a natural narrowing funnel.
    _n_total = adata_flt.n_obs
    _steps = [
        (f"not in doublet-dominated cluster ({', '.join(doublet_dominated_clusters)})", adata_flt.obs["pass_not_doublet_cluster"]),
        ("CMO hashing status == Singlet (and top CMO in panel)", adata_flt.obs["pass_cmo_singlet"]),
        ("Scrublet not predicted_doublet", adata_flt.obs["pass_not_scrublet_doublet"]),
    ]
    _steps = sorted(_steps, key=lambda s: -int(s[1].sum()))

    _remaining_mask = pd.Series(True, index=adata_flt.obs.index)
    _rows = []
    for _label, _step_mask in _steps:
        _before = int(_remaining_mask.sum())
        _remaining_mask &= _step_mask
        _after = int(_remaining_mask.sum())
        _lost = _before - _after
        _rows.append(
            f"| `{_label}` | {_before:,} | {_after:,} | {_lost:,} | {_lost / _before:.1%} |"
        )

    # Exposed as its own top-level name, and written back to adata_flt.obs, so
    # downstream cells get a real, trackable marimo dependency edge instead of
    # just sharing a ref to "adata_flt". See the marimo-pair race-condition note.
    pass_strict_qc_mask = _remaining_mask
    adata_flt.obs["pass_strict_qc"] = pass_strict_qc_mask
    strict_qc_computed = _ok

    mo.md(f"""
    **Strict QC filter breakdown** (steps ordered from least to most restrictive):

    | Filter step | Before | After | Lost | Lost % |
    |---|---|---|---|---|
    {chr(10).join(_rows)}

    **Net:** {int(pass_strict_qc_mask.sum()):,} of {_n_total:,} barcodes ({pass_strict_qc_mask.sum() / _n_total:.1%}) survive strict QC (`pass_strict_qc`).
    """)
    return (strict_qc_computed,)


@app.cell
def filter_strict_qc(adata_flt, strict_qc_computed):
    strict_qc_computed  # ran after strict QC flags were computed

    adata_qc = adata_flt[adata_flt.obs["pass_strict_qc"]].copy()
    adata_qc
    return (adata_qc,)


@app.cell(hide_code=True)
def pca_by_timepoint(
    adata_qc,
    ec_diff_palette,
    pca_axis_title,
    pca_computed,
    plt,
    strict_qc_computed,
    timepoint_order,
):
    strict_qc_computed  # ran after strict QC flags were computed
    pca_computed  # ran after PCA

    # Plotted directly with matplotlib (not sc.pl.pca) so the 15-category legend
    # gets full control over layout, drawn outside the axes -- scanpy's default
    # right-margin legend was getting clipped with this many categories.
    _pca_coords = adata_qc.obsm["X_pca"][:, [0, 1]]
    _timepoints = adata_qc.obs["timepoint_scanpy"].to_numpy()

    _fig, _ax = plt.subplots(figsize=(7, 6))
    for _tp in timepoint_order:
        if _tp == "Negative":
            continue
        _mask = _timepoints == _tp
        if not _mask.any():
            continue
        _ax.scatter(_pca_coords[_mask, 0], _pca_coords[_mask, 1], s=5, color=ec_diff_palette[_tp], label=_tp)

    _ax.set_xlabel(pca_axis_title(adata_qc, 0))
    _ax.set_ylabel(pca_axis_title(adata_qc, 1))
    _ax.set_title("PCA colored by timepoint (pass_strict_qc only)")
    _ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, markerscale=2, frameon=False)

    plt.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def umap_by_timepoint(
    adata_qc,
    ec_diff_palette,
    plt,
    sc,
    strict_qc_computed,
    timepoint_order,
    umap_computed,
):
    strict_qc_computed  # ran after strict QC flags were computed
    umap_computed  # ran after UMAP was computed

    _umap_coords = adata_qc.obsm["X_umap"]
    _timepoints = adata_qc.obs["timepoint_scanpy"].to_numpy()

    _fig, _axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: timepoint, drawn directly with matplotlib (not sc.pl.umap) so the
    # 15-category legend gets full control over layout, drawn outside the axes --
    # scanpy's default right-margin legend was getting clipped with this many
    # categories.
    for _tp in timepoint_order:
        if _tp == "Negative":
            continue
        _mask = _timepoints == _tp
        if not _mask.any():
            continue
        _axes[0].scatter(_umap_coords[_mask, 0], _umap_coords[_mask, 1], s=5, color=ec_diff_palette[_tp], label=_tp)
    _axes[0].set_xlabel("UMAP1")
    _axes[0].set_ylabel("UMAP2")
    _axes[0].set_title("UMAP colored by timepoint")
    _axes[0].legend(loc="center left", bbox_to_anchor=(1.15, 0.5), fontsize=8, markerscale=2, frameon=False)

    # Right: %mito is continuous, so scanpy's usual colorbar (not a category
    # legend) is not affected by the same clipping issue -- plotted the same way
    # for a consistent look.
    _mito_scatter = _axes[1].scatter(
        _umap_coords[:, 0], _umap_coords[:, 1], s=5, c=adata_qc.obs["pct_counts_mt"], cmap="cividis",
    )
    _axes[1].set_xlabel("UMAP1")
    _axes[1].set_ylabel("UMAP2")
    _axes[1].set_title("UMAP colored by %mito")
    _fig.colorbar(_mito_scatter, ax=_axes[1], label="% mitochondrial counts")

    plt.tight_layout()
    _fig

    # Leiden clusters shown separately -- ec_diff_palette only covers timepoint
    # labels, not the Leiden cluster IDs, so this uses scanpy's own default
    # categorical palette instead.
    sc.pl.umap(
        adata_qc,
        color="leiden",
        size=5,
        title="UMAP colored by Leiden cluster",
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Saving it to file
    """)
    return


@app.cell(hide_code=True)
def write_qc_annotations_tsv(
    adata_flt,
    ch2_outdir_root,
    project_root,
    strict_qc_computed,
):
    strict_qc_computed  # ran after strict QC flags were computed

    _outfile = ch2_outdir_root / "adata_flt_qc_annotations.tsv"
    _outfile.parent.mkdir(parents=True, exist_ok=True)
    adata_flt.obs.to_csv(_outfile, sep="	")

    # Show only the path relative to the repo, not the full local filesystem path
    _outfile.relative_to(project_root)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # ATAC data processing with snapatac2
    """)
    return


@app.cell
def import_snapatac2():
    import snapatac2 as snap

    return (snap,)


@app.cell
def load_chrom_sizes(pd, project_root):
    chrom_dict = pd.read_csv(project_root / "annotations/GRCh38_EBV.chrom.sizes.no.alt.tsv", sep="	", header=None, names=["chr", "size"])
    chrom_dict = chrom_dict.set_index("chr")["size"].to_dict()
    return (chrom_dict,)


@app.cell(hide_code=True)
def atac_import_intro(mo):
    mo.md(r"""
    ### Import ATAC fragments

    We import fragments restricted to the RNA-side `pass_strict_qc` whitelist rather than re-filtering on ATAC metrics separately. We trust the RNA-based QC (doublet clusters, Scrublet, CMO doublets/negatives) as the primary cell call, and only compute ATAC QC metrics for visibility.
    """)
    return


@app.cell
def load_strict_qc_whitelist(ch2_outdir_root, pd):
    # See atac_import_intro above.
    _qc_annotations = pd.read_csv(
        ch2_outdir_root / "adata_flt_qc_annotations.tsv", sep="	", index_col=0
    )
    assigned_cells = _qc_annotations.index[_qc_annotations["pass_strict_qc"]].to_list()
    len(assigned_cells)
    return (assigned_cells,)


@app.cell
def import_atac_fragments(
    assigned_cells,
    ch2_data_root_path,
    chrom_dict,
    igvf_gencode_gtf_path,
    snap,
):
    atac_adata = snap.pp.import_fragments(
        ch2_data_root_path / "fragments/10x_15timepoints_channel2.fragments.tsv.gz",
        sorted_by_barcode=False,
        chrom_sizes=chrom_dict,
        whitelist=assigned_cells,
    )
    snap.metrics.tsse(atac_adata, igvf_gencode_gtf_path)
    atac_adata
    return (atac_adata,)


@app.cell
def atac_tsse_preview(atac_adata):
    atac_adata.obs["tsse"]
    return


@app.cell(hide_code=True)
def atac_qc_plot_intro(mo):
    mo.md(r"""
    ### ATAC QC, for visibility only

    Unique fragments vs. TSSE, with fixed thresholds (`atac_qc_thresholds`) shown for reference. We are not applying this as an additional filter, see the note below, since the barcode set already inherits the RNA `pass_strict_qc` whitelist at import.
    """)
    return


@app.cell
def atac_qc_thresholds(atac_adata):
    df = atac_adata.obs[["n_fragment", "tsse"]].copy()

    # Fixed thresholds for visibility only, see the note below: we do not actually
    # filter on these, the ATAC barcode set already inherits the RNA pass_strict_qc
    # whitelist at import time.
    min_frag = 500
    min_tsse = 5.0
    return df, min_frag, min_tsse


@app.cell
def atac_fragment_tsse_plot(alt, df, min_frag, min_tsse, okabe_ito_palette):
    # See atac_qc_plot_intro above.
    alt.renderers.enable("html")

    plot_df = df.copy()
    plot_df["pass_qc"] = (
        (plot_df["n_fragment"] >= min_frag) &
        (plot_df["tsse"] >= min_tsse) &
        (plot_df["n_fragment"] <= 15_000)
    ).astype(str)

    print(plot_df["pass_qc"].value_counts())

    _pass_qc_scale = alt.Scale(domain=["True", "False"], range=[okabe_ito_palette[1], okabe_ito_palette[5]])

    base = alt.Chart(plot_df).encode(
        x=alt.X("n_fragment:Q", scale=alt.Scale(type="log"), title="Unique fragments"),
        y=alt.Y("tsse:Q", title="TSSE"),
        color=alt.Color("pass_qc:N", title="Pass QC", scale=_pass_qc_scale)
    )

    scatter = base.mark_circle(size=10, opacity=0.3)

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

    (top_hist & (scatter | right_hist)).resolve_legend(color="shared")
    return (plot_df,)


@app.cell
def map_timepoint_to_atac(
    adata_flt,
    atac_adata,
    atac_spectral_umap_leiden_computed,
):
    # Forces a deterministic write order into atac_adata.obs: both this cell
    # and atac_spectral_umap_leiden mutate atac_adata.obs in place (adding
    # "timepoint_scanpy" and "leiden" respectively) with no data dependency
    # between them, so without this, which column lands first -- and thus
    # atac_adata.obs's final column order -- depends on arbitrary execution
    # order rather than being reproducible.
    atac_spectral_umap_leiden_computed  # ran after spectral embedding, UMAP, KNN, and Leiden clustering

    def _map_timepoint_to_atac(atac_adata, adata_flt):
        # Every whitelisted barcode is a CMO Singlet (pass_strict_qc requires it,
        # with no negative-tag rescue for this dataset), so timepoint_scanpy is
        # already a real timepoint value for all of them -- no rescue needed here.
        atac_adata.obs["timepoint_scanpy"] = adata_flt.obs.loc[atac_adata.obs_names, "timepoint_scanpy"].astype("category")
        return True

    timepoint_mapped_to_atac = _map_timepoint_to_atac(atac_adata, adata_flt)
    return (timepoint_mapped_to_atac,)


@app.cell
def atac_timepoint_counts(atac_adata, timepoint_mapped_to_atac):
    timepoint_mapped_to_atac  # ran after timepoint_scanpy was mapped onto atac_adata

    atac_adata.obs["timepoint_scanpy"].value_counts()
    return


@app.cell
def atac_no_filtering_note(mo):
    mo.md("""
    Intentionally not calling `snap.pp.filter_cells` here, the barcode set already reflects the RNA `pass_strict_qc` whitelist used at import (see `atac_import_intro` above), so we trust the RNA-based QC rather than layering an additional ATAC-metric-based filter on top. The QC plot above is for visibility only.
    """)
    return


@app.cell(hide_code=True)
def atac_processing_intro(mo):
    mo.md(r"""
    ### ATAC dimensionality reduction and clustering

    Standard snapATAC2 processing on the tile matrix: feature selection, spectral embedding, UMAP, KNN graph, and Leiden clustering.
    """)
    return


@app.cell
def atac_add_tile_matrix(atac_adata, snap):
    def _run_add_tile_matrix(adata):
        snap.pp.add_tile_matrix(adata)
        return True

    tile_matrix_added = _run_add_tile_matrix(atac_adata)
    return (tile_matrix_added,)


@app.cell
def atac_select_features(atac_adata, snap, tile_matrix_added):
    tile_matrix_added  # ran after the tile matrix was added

    def _run_select_features(adata, n_features):
        snap.pp.select_features(adata, n_features=n_features)
        return True

    features_selected = _run_select_features(atac_adata, 250_000)
    return (features_selected,)


@app.cell
def atac_spectral_umap_leiden(atac_adata, features_selected, snap):
    features_selected  # ran after feature selection

    def _run_atac_spectral_umap_leiden(adata):
        snap.tl.spectral(adata)
        snap.tl.umap(adata)
        snap.pp.knn(adata)
        snap.tl.leiden(adata)
        return True

    atac_spectral_umap_leiden_computed = _run_atac_spectral_umap_leiden(atac_adata)
    return (atac_spectral_umap_leiden_computed,)


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
def atac_umap_intro(mo):
    mo.md(r"""
    ### ATAC UMAP

    Colored by timepoint, Leiden cluster, and RNA %mito, using the helper functions above.
    """)
    return


@app.cell
def atac_umap_by_timepoint_and_leiden(
    atac_adata,
    atac_spectral_umap_leiden_computed,
    ec_diff_palette,
    plt,
    sc,
    timepoint_mapped_to_atac,
    timepoint_order,
):
    timepoint_mapped_to_atac  # ran after timepoint_scanpy was mapped onto atac_adata
    atac_spectral_umap_leiden_computed  # ran after spectral embedding, UMAP, KNN, and Leiden clustering

    # Timepoint plotted directly with matplotlib (not sc.pl.umap) so the
    # 15-category legend gets full control over layout, drawn outside the axes --
    # scanpy's default right-margin legend gets clipped with this many categories
    # (see the RNA umap_by_timepoint cell above).
    _umap_coords = atac_adata.obsm["X_umap"]
    _timepoints = atac_adata.obs["timepoint_scanpy"].to_numpy()

    _fig, _ax = plt.subplots(figsize=(8, 6))
    for _tp in timepoint_order:
        if _tp == "Negative":
            continue
        _mask = _timepoints == _tp
        if not _mask.any():
            continue
        _ax.scatter(_umap_coords[_mask, 0], _umap_coords[_mask, 1], s=5, color=ec_diff_palette[_tp], label=_tp)
    _ax.set_xlabel("UMAP1")
    _ax.set_ylabel("UMAP2")
    _ax.set_title("ATAC UMAP colored by timepoint")
    _ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, markerscale=2, frameon=False)
    plt.tight_layout()
    _fig

    # Leiden cluster count is small enough that scanpy's own legend handles it fine.
    sc.pl.umap(
        atac_adata,
        color="leiden",
        size=5,
        title="ATAC UMAP colored by Leiden cluster",
    )
    return


@app.cell(hide_code=True)
def atac_umap_and_confusion_vs_rna_cluster(
    adata_qc,
    alt,
    atac_adata,
    atac_spectral_umap_leiden_computed,
    atac_umap_dataframe,
    atac_umap_scatter,
    mo,
    okabe_ito_palette,
    pd,
    strict_qc_computed,
):
    atac_spectral_umap_leiden_computed  # ran after ATAC spectral embedding, UMAP, and Leiden clustering
    strict_qc_computed  # ran after strict QC flags were computed (adata_qc)

    # atac_adata's barcodes are a subset of adata_qc's (every ATAC barcode has an
    # RNA counterpart here), so no explicit intersection is needed, but the count
    # is still worth stating explicitly rather than assuming it.
    _shared_barcodes = atac_adata.obs_names.intersection(adata_qc.obs_names)
    _n_shared = len(_shared_barcodes)
    _rna_cluster_on_atac = adata_qc.obs.loc[_shared_barcodes, "leiden"]

    _rna_categories = sorted(_rna_cluster_on_atac.astype(str).unique(), key=int)
    _rna_palette = {_cl: okabe_ito_palette[_i % len(okabe_ito_palette)] for _i, _cl in enumerate(_rna_categories)}

    _umap_df = atac_umap_dataframe(atac_adata[_shared_barcodes], [])
    _umap_df["rna_leiden"] = _rna_cluster_on_atac.astype(str).to_numpy()

    _umap_chart = atac_umap_scatter(
        _umap_df, "rna_leiden", "N", "ATAC UMAP colored by RNA cluster",
        color_scale=alt.Scale(domain=list(_rna_palette), range=list(_rna_palette.values())),
    )

    # Confusion matrix: how ATAC leiden clusters (computed on the spectral
    # embedding) map onto RNA leiden clusters (computed independently on gene
    # expression). Cluster IDs are arbitrary and not comparable by number across
    # modalities, so this crosstab is the only way to track correspondence.
    _atac_cluster_on_shared = atac_adata.obs.loc[_shared_barcodes, "leiden"].astype(str)
    _confusion = pd.crosstab(_atac_cluster_on_shared.rename("atac_cluster"), _rna_cluster_on_atac.astype(str).rename("rna_cluster"))

    _atac_order = sorted(_confusion.index, key=int)
    _rna_order = _confusion.idxmax(axis=0).sort_values(key=lambda s: s.astype(int)).index.tolist()
    _atac_rank = {v: i for i, v in enumerate(_atac_order)}
    _rna_rank = {v: i for i, v in enumerate(_rna_order)}

    _confusion_long = _confusion.reset_index().melt(id_vars="atac_cluster", var_name="rna_cluster", value_name="n_cells")
    _confusion_long = _confusion_long[_confusion_long["n_cells"] > 0].copy()
    # Altair's sort=<list> shorthand can silently fail under VegaFusion; precomputing
    # the rank as an explicit data column and sorting via EncodingSortField
    # sidesteps it.
    _confusion_long["atac_rank"] = _confusion_long["atac_cluster"].map(_atac_rank)
    _confusion_long["rna_rank"] = _confusion_long["rna_cluster"].map(_rna_rank)

    _heatmap = alt.Chart(_confusion_long).mark_rect().encode(
        x=alt.X("rna_cluster:N", title="RNA cluster", sort=alt.EncodingSortField(field="rna_rank", op="min")),
        y=alt.Y("atac_cluster:N", title="ATAC cluster", sort=alt.EncodingSortField(field="atac_rank", op="min")),
        color=alt.Color("n_cells:Q", title="Cells", scale=alt.Scale(scheme="cividis")),
        tooltip=["atac_cluster", "rna_cluster", "n_cells"],
    )
    _labels = alt.Chart(_confusion_long).mark_text(fontSize=9).encode(
        x=alt.X("rna_cluster:N", sort=alt.EncodingSortField(field="rna_rank", op="min")),
        y=alt.Y("atac_cluster:N", sort=alt.EncodingSortField(field="atac_rank", op="min")),
        text="n_cells:Q",
        color=alt.condition(alt.datum.n_cells > _confusion_long["n_cells"].max() / 2, alt.value(okabe_ito_palette[0]), alt.value("white")),
    )

    _confusion_chart = (_heatmap + _labels).properties(
        title="ATAC vs. RNA leiden cluster assignment (cell counts, columns ordered to show correspondence)",
        width=450, height=450,
    ).configure_view(strokeWidth=0)

    def _tight_row(*items):
        # mo.hstack with widths=None adds no wrapper/flex styling around children,
        # so block-level chart divs just fill the row (no slack left for
        # justify-content to redistribute). Build the flex row by hand instead.
        _items_html = "".join(
            f'<div style="flex: 0 0 auto;">{mo.as_html(it).text}</div>' for it in items
        )
        return mo.Html(f'<div style="display:flex; justify-content:flex-start; gap:1rem;">{_items_html}</div>')

    mo.vstack([
        mo.md(f"**{_n_shared:,} barcodes shared between `atac_adata` and `adata_qc`** (of {atac_adata.n_obs:,} ATAC and {adata_qc.n_obs:,} RNA barcodes); both panels below are restricted to this shared set."),
        _tight_row(_umap_chart, _confusion_chart),
    ])
    return


@app.cell
def atac_umap_by_rna_mito(
    adata_flt,
    alt,
    atac_adata,
    atac_spectral_umap_leiden_computed,
    atac_umap_dataframe,
    atac_umap_scatter,
):
    atac_spectral_umap_leiden_computed  # ran after spectral embedding, UMAP, KNN, and Leiden clustering

    _mito_df = atac_umap_dataframe(atac_adata, [])
    _mito_df["pct_counts_mt"] = adata_flt.obs.loc[atac_adata.obs_names, "pct_counts_mt"].to_numpy()

    _mito_chart = atac_umap_scatter(
        _mito_df, "pct_counts_mt", "Q", "ATAC UMAP colored by RNA %mito", color_scale=alt.Scale(scheme="cividis"),
    )

    _mito_chart
    return


@app.cell(hide_code=True)
def joint_analyses_header(mo):
    mo.md(r"""
    ## Joint analyses (RNA + ATAC)
    """)
    return


@app.cell
def ncount_atac_vs_rna_counts(
    adata_flt,
    alt,
    atac_adata,
    okabe_ito_palette,
    pd,
    plot_df,
):
    # Do the barcodes with low RNA quality (low total_counts) also show low ATAC
    # quality (few unique fragments)? Uses the same pass_qc flag from the plot
    # above (ATAC fragment/TSSE thresholds) to color points.
    _scatter_df = pd.DataFrame({
        "n_fragment": atac_adata.obs["n_fragment"].to_numpy(),
        "total_counts_rna": adata_flt.obs.loc[atac_adata.obs_names, "total_counts"].to_numpy(),
        "pass_qc": plot_df["pass_qc"].to_numpy(),
    })

    _pass_qc_scale = alt.Scale(domain=["True", "False"], range=[okabe_ito_palette[1], okabe_ito_palette[5]])

    _ncount_chart = alt.Chart(_scatter_df).mark_circle(size=10, opacity=0.3).encode(
        x=alt.X(
            "n_fragment:Q", title="Unique ATAC fragments (nCount_ATAC)",
            scale=alt.Scale(type="log", domain=[10, _scatter_df["n_fragment"].max()]),
        ),
        y=alt.Y(
            "total_counts_rna:Q", title="Total RNA counts",
            scale=alt.Scale(type="log", domain=[10, _scatter_df["total_counts_rna"].max()]),
        ),
        color=alt.Color("pass_qc:N", title="Pass ATAC QC", scale=_pass_qc_scale),
    ).properties(title="nCount_ATAC vs. total RNA counts", width=500, height=500)

    _ncount_chart
    return


@app.cell
def saving_h5ad(
    adata_qc,
    atac_adata,
    atac_spectral_umap_leiden_computed,
    ch2_outdir_root,
    project_root,
    strict_qc_computed,
):
    strict_qc_computed  # ran after strict QC flags were computed (RNA)
    atac_spectral_umap_leiden_computed  # ran after spectral embedding, UMAP, KNN, and Leiden clustering (ATAC)

    # Saving the RNA and ATAC anndata objects.
    _rna_outfile = ch2_outdir_root / "rna_adata.h5ad"
    _atac_outfile = ch2_outdir_root / "atac_adata.h5ad"
    _rna_outfile.parent.mkdir(parents=True, exist_ok=True)

    adata_qc.write_h5ad(_rna_outfile)
    atac_adata.write_h5ad(_atac_outfile)

    # Show only the paths relative to the repo, not the full local filesystem path
    [_rna_outfile.relative_to(project_root), _atac_outfile.relative_to(project_root)]
    return


if __name__ == "__main__":
    app.run()
