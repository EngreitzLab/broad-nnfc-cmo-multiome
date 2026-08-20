# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "altair==6.2.2",
#     "anndata==0.13.2",
#     "igraph>=1.0.0",
#     "igvf-utils==3.1.1",
#     "ipython>=9.13.0",
#     "marimo>=0.23.3",
#     "matplotlib==3.11.1",
#     "numpy==2.4.6",
#     "pandas==2.3.3",
#     "python-dotenv==1.2.2",
#     "scanpy[scrublet]==1.12.2",
#     "scclr==0.1.0",
#     "scikit-learn==1.9.0",
#     "scipy==1.18.0",
#     "seaborn==0.13.2",
#     "snapatac2==2.9.0",
#     "statsmodels==0.14.6",
#     "tabulate==0.10.0",
#     "vegafusion==2.0.3",
#     "vl-convert-python==1.9.0.post1",
# ]
#
# [tool.uv.sources]
# scclr = { git = "https://github.com/cleartools/scclr.git" }
# ///

import marimo

__generated_with = "0.23.14"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Introduction

    This notebook processes the 10x multiome data (**produced by Chris McGinnis**) for the endothelial differentiation time course: 5 timepoints (d0-d4), 5 biological replicates each. Samples were multiplexed using the MULTI-seq technique with CMO (Cell Multiplexing Oligo) barcodes. Data was processed with the IGVF pipeline using kallisto and the GENCODE v43 annotation, and CMO quantification was performed with the `kite` workflow from the `kallisto-bustools` suite. Single channel (no channel1/channel2 split).

    The goals of this notebook are to:

    1. Perform quality control filtering on the RNA data
    2. Run CMO hash classification (mimicking `Seurat::HTODemux`) on QC-passing cells
    3. Assign each cell barcode to a CMO (and therefore to a timepoint/replicate)
    4. Process the corresponding ATAC data with `snapatac2`

    IGVF portal accessions

    - Analysis set : [IGVFDS1612ZNCA](https://data.igvf.org/analysis-sets/IGVFDS1612ZNCA/)
    - Gene count matrix: [IGVFFI5144GDQB](https://data.igvf.org/matrix-files/IGVFFI5144GDQB/)
    - Fragment file: [IGVFFI6722CCZW](https://data.igvf.org/tabular-files/IGVFFI6722CCZW/)
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

    from dotenv import find_dotenv, load_dotenv
    from igvf_utils.connection import Connection
    from pathlib import Path
    from scipy.io import mmread

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

    return (
        Connection,
        Path,
        ad,
        alt,
        ec_diff_palette,
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
    # Imports and palettes
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Download data from the portal

    To be able to download the files you need to create a `.env` file in the repository root and add the lines
    ```
    IGVF_API_KEY=<your-api-key>
    IGVF_SECRET_KEY=<your-api-secret>
    ```
    """)
    return


@app.cell
def _(Connection, Path, project_root):
    data_root_path = Path(project_root / "data/10x_multi_5_timepoints_mcginnis/")

    outdir_root = Path(project_root /"results/10x_multi_5_timepoints_mcginnis")


    # 1. Initialize the connection targeting the production portal
    # It automatically picks up IGVF_API_KEY and IGVF_SECRET_KEY from the environment
    conn = Connection(igvf_mode="prod")

    # 2. Files to download: (accession, destination subfolder)
    _files_to_download = [
        ("IGVFFI5144GDQB", data_root_path / "rna/h5ad"),
        ("IGVFFI6722CCZW", data_root_path / "atac/fragments"),
    ]

    for _file_accession, _dest_dir in _files_to_download:
        # Make sure the destination folder exists before writing into it
        _dest_dir.mkdir(parents=True, exist_ok=True)

        # Skip re-downloading if a file for this accession is already there (name
        # matched via the portal's own href, since the served filename isn't known
        # until conn.download() reads it off the response's Content-Disposition header)
        _expected_name = conn.get(rec_ids=_file_accession)["href"].rsplit("/", 1)[-1]
        _expected_path = _dest_dir / _expected_name
        if _expected_path.exists():
            print(f"{_expected_path} already exists, skipping download.")
            continue

        print(f"Downloading {_file_accession} to {_dest_dir}...")
        _downloaded_path = conn.download(rec_id=_file_accession, directory=str(_dest_dir))
        print(f"Successfully downloaded: {_downloaded_path}")
    return data_root_path, outdir_root


@app.cell(hide_code=True)
def load_gene_metadata_intro(mo):
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
def get_igvf_gencode_intro(mo):
    mo.md(r"""
    ## Downloading the GENCODE GTF (IGVF-hosted)
    The full GENCODE v43 GTF (used for the ATAC TSSE metric below) is downloaded from IGVF's reference-file API if not already present locally, so the notebook doesn't depend on a symlinked annotation file outside this repo.
    """)
    return


@app.cell
def get_igvf_gencode(Connection, gtf_to_gene_metadata, project_root):
    igvf_gencode_gtf_path = project_root / "annotations/IGVFFI9573KOZR.gtf.gz"
    igvf_gencode_gtf_path.parent.mkdir(parents=True, exist_ok=True)

    if not igvf_gencode_gtf_path.exists():
        _conn = Connection(igvf_mode="prod")
        _conn.download(rec_id="IGVFFI9573KOZR", directory=str(igvf_gencode_gtf_path.parent))

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
    data_root_path,
    gene_metadata_df,
    gzip,
    igvf_gencode_gtf_path,
    re,
):
    try:
        _h5ad_path = data_root_path / "rna/h5ad/IGVFFI5144GDQB.h5ad"
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
        # Use the GENCODE v43 GTF (downloaded above) as ground truth instead: every
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
def plot_knee_plot(adata, alt, mo, np, pd, qc_metrics_computed):
    mo.stop( not qc_metrics_computed)

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

    knee_umi_selection = alt.selection_interval(encodings=["y"], value={"y": [1500, 30_000]})

    _knee_chart = alt.Chart(_knee_df).mark_circle(size=15, color="black").encode(
        x=alt.X("rank:Q", scale=alt.Scale(type="log"), title="Cell rank"),
        y=alt.Y("n_umis:Q", scale=alt.Scale(type="log"), title="Total UMIs"),
        tooltip=[alt.Tooltip("rank:Q", title="Rank"), alt.Tooltip("n_umis:Q", title="# UMIs", format=",")],
    ).add_params(knee_umi_selection).properties(width=550, height=400, title="Knee plot")

    # mo.ui.altair_chart calls chart.to_json() without format="vega", which
    # conflicts with the vegafusion transformer enabled globally in `imports`
    # (it requires that format explicitly). This chart's data is already a
    # ~2000-point downsample, so vegafusion isn't needed here anyway --
    # temporarily drop back to the default transformer just for this widget.
    with alt.data_transformers.enable("default"):
        knee_chart_ui = mo.ui.altair_chart(_knee_chart, chart_selection=False, legend_selection=False)
    knee_chart_ui
    return knee_chart_ui, knee_full_counts


@app.cell(hide_code=True)
def knee_plot_brush_count(knee_chart_ui, knee_full_counts, mo):
    # Read the brush's y-range (falls back to the default 1000-10,000 before any
    # interaction, since the interval's initial "value" isn't reported back until
    # the user actually drags it) and count how many barcodes it covers, using the
    # full (non-downsampled) counts array for an accurate number.
    _selections = knee_chart_ui.selections
    if _selections:
        _lo, _hi = next(iter(_selections.values()))["n_umis"]
    else:
        _lo, _hi = 1500, 30_000

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
def mask_filter_cells(adata, mo, np, pd, qc_metrics_computed):
    mo.stop(not qc_metrics_computed)

    # QC pass mask mirroring the cell-level filters applied below.
    # min_counts_cutoff set to 1500 based on the external final_anno_igvf0.tsv
    # annotation: its own curated real cells extend down to 578 UMIs (5th
    # percentile 2,445), so the previous 5000 floor was excluding a real,
    # externally-validated chunk of cells (~4,864 barcodes). 1500 sits below the
    # 10th percentile of that population (2,980) while still well above the
    # empty-droplet noise floor.
    _max_counts_cutoff = 30_000
    _min_counts_cutoff = 1_500

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
    _g.ax_joint.axhline(y=mito_mad_cutoff_2_5, color="black", linestyle="--")
    _g.ax_joint.annotate(
        f"{mito_mad_cutoff_2_5:.1f}% (2.5 MAD, applied cutoff)",
        xy=(_df["total_counts"].max(), mito_mad_cutoff_2_5),
        xytext=(20, 5),
        textcoords="offset points",
        va="bottom",
        ha="right",
        color="black",
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

    Starting from {adata.n_obs:,} raw barcodes, the level-1 filters (`min_counts >= 1000`, `max_counts <= 10,000`, `min_genes >= 200`) leave **{int(pre_mito_qc_mask.sum()):,} barcodes**:

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
    # since its .value can'\''t be read in the same cell that creates it.
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

    # KDE on a fixed grid rather than Altair'\''s transform_density: VegaFusion
    # cannot pre-aggregate transform_density, so it would embed the full raw
    # dataframe and trip the output-too-large limit.
    _metric = qc_violin_metric_dropdown.value
    _label = qc_violin_metric_dropdown.selected_key

    _values = adata_flt.obs[_metric].to_numpy()
    _grid = np.linspace(_values.min(), _values.max(), 200)
    _density = _gaussian_kde_qc_violin(_values)(_grid)
    _density_df = pd.DataFrame({_metric: _grid, "density": _density})

    # Mirror the density manually (x0/x2) instead of Altair'\''s stack="center":
    # with a per-row groupby of size 1, Vega-Lite'\''s stack transform did not
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


@app.cell
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
        sc.tl.leiden(adata, flavor="igraph", resolution=0.4, n_iterations=2, random_state=0)
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
    # dedicated deep-dive (see the cluster 2 sections below).
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
    # not here -- this cell only defines it, since its .value can't be read in the
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
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Load CMO counts and perform assignment
    """)
    return


@app.cell
def initialize_cmo(ad, adata_flt, data_root_path):
    # path to the cmo counts
    _cmo_counts_path = data_root_path / "cmo_counts/adata.h5ad"

    adata_cmo = ad.read_h5ad(_cmo_counts_path)

    # CMO barcodes need the same per-sample suffix as adata_flt's RNA barcodes to match up
    _barcode_suffix = "_" + adata_flt.obs_names[0].split("_", 1)[1]
    adata_cmo.obs_names = adata_cmo.obs_names + _barcode_suffix

    adata_cmo.var["gene_name"] = adata_cmo.var_names

    # Filter adata_cmo data to cells in adata_flt
    _cells_in_use = adata_flt.obs_names.intersection(adata_cmo.obs_names.to_list())
    adata_cmo = adata_cmo[_cells_in_use, :]

    # wtc11_endo_d1_5, wtc11_endo_d2_3, and wtc11_endo_d4_2 have no real signal in
    # this dataset -- confirmed directly against the raw R3 FASTQ (zero
    # occurrences in a 5M-read sample) and matching final_anno_igvf0.tsv (these
    # 3 tags are also entirely absent from its own sample_multiseq column). Those
    # wells simply have no data in this experiment. Drop them here rather than
    # special-casing them in every downstream cell (geometric mean, thresholds,
    # assignment) -- their near-zero counts are just background noise, and
    # percentile-based thresholding would otherwise always manufacture some
    # spurious "positive" calls for them regardless. See CMO_ISSUE.txt.
    _dead_cmo_tags = ["wtc11_endo_d1_5", "wtc11_endo_d2_3", "wtc11_endo_d4_2"]
    adata_cmo = adata_cmo[:, ~adata_cmo.var["gene_name"].isin(_dead_cmo_tags)].copy()

    # CMO tag names here aren't "CMO1".."CMO25" like channel1/channel2; they're
    # descriptive sample_multiseq names (e.g. "wtc11_endo_d3_2" = timepoint d3,
    # replicate 2), and the timepoint is always the third underscore-separated
    # field. Confirmed against data/10x_multi_5_timepoints_mcginnis/cmo_counts/
    # final_anno_igvf0.tsv's own sample_multiseq column, which has the same
    # names prefixed with "ansuman-satpathy:". Note the replicate count isn't
    # uniform across timepoints here (3 for d0, 5 for d1-d4, 20 CMOs remaining
    # after dropping the 3 dead tags above), unlike channel1/channel2's uniform
    # 5-per-timepoint/25-CMO panel.
    cmo_to_timepoint = {name: name.split("_")[2] for name in adata_cmo.var_names}

    adata_cmo
    return adata_cmo, cmo_to_timepoint


@app.cell
def _(adata_cmo, np, okabe_ito_palette, pd, plt):
    # 1. Calculate total counts per CMO
    col_sums = np.ravel(adata_cmo.X.sum(axis=0))

    # 2. Pair with CMO names and sort values
    cmo_counts = pd.Series(col_sums, index=adata_cmo.var_names).sort_values(ascending=True)

    # 3. Plot horizontal bar chart
    plt.figure(figsize=(10, max(4, len(cmo_counts) * 0.35)))  # Scales height based on number of CMOs
    plt.barh(cmo_counts.index, cmo_counts.values, color=okabe_ito_palette[0], edgecolor='black')

    plt.title('Total Counts per CMO', fontsize=14)
    plt.xlabel('Total Counts', fontsize=12)
    plt.ylabel('CMO Name', fontsize=12)
    plt.grid(axis='x', linestyle='--', alpha=0.7)

    plt.tight_layout()
    plt.show()
    return


@app.cell
def cmo_clr_normalization(adata_cmo, np):
    _cmo = adata_cmo.X
    _cmo = _cmo.toarray() if hasattr(_cmo, "toarray") else np.asarray(_cmo)

    # 1. Compute Seurat-style geometric mean per cell across all CMO features
    # Seurat sums log1p of non-zero values, divides by total features, then takes exp.
    # Divide by the actual panel size (_cmo.shape[1]) -- 20 CMOs after
    # initialize_cmo drops the 3 dead tags (wtc11_endo_d1_5, d2_3, d4_2; see
    # CMO_ISSUE.txt), not channel1/channel2's uniform 25.
    _log_counts = np.log1p(np.where(_cmo > 0, _cmo, 0))
    _gm = np.exp(np.sum(_log_counts, axis=1) / _cmo.shape[1])

    # 2. Divide raw counts by the geometric mean, then apply log1p
    clr = np.log1p(_cmo / _gm[:, None])
    return (clr,)


@app.cell
def find_cmo_thresholds(clr, np):
    # Per-tag Gaussian Mixture Model on each tag's own CLR distribution, using a
    # shared posterior-probability cutoff instead of (a) a single global CLR
    # quantile (previously 0.95) shared across all 20 CMOs, or (b) each tag's own
    # individually-optimal split (posterior > 0.5). (b) validates every tag
    # one-vs-rest in isolation and ignores that a cell only counts as a true
    # Singlet if exactly one of 20 simultaneous per-tag tests fires: at
    # posterior > 0.5 that inflated the doublet rate so much (n_pos>=2 for
    # roughly half of adata_flt) that overlap with final_anno_igvf0.tsv barely
    # improved over the old approach (70.0%). Sweeping one SHARED posterior
    # cutoff across all tags (still tag-adaptive -- each tag keeps its own GMM
    # shape) found a real joint optimum around 0.9999: 78.2% IGVF0 overlap, 81.0%
    # concordance -- beating both the old global-quantile peak (75.1% overlap /
    # 81.4% concordance at quantile 0.94) and the naive per-tag GMM default
    # (70.0% / 79.8%).
    positive_posterior_cutoff = 0.9999

    from sklearn.mixture import GaussianMixture

    def _gmm_fit_column(values):
        _gmm = GaussianMixture(n_components=2, random_state=0, n_init=3)
        _gmm.fit(values.reshape(-1, 1))
        _high = int(np.argmax(_gmm.means_.flatten()))
        _probs_high = _gmm.predict_proba(values.reshape(-1, 1))[:, _high]
        _is_positive = _probs_high > positive_posterior_cutoff

        # Grid-crossing approximation of the decision boundary at the same
        # cutoff, for display only (cmo_threshold_vs_assignment) -- classification
        # itself uses _is_positive directly, not this reconstructed value.
        _grid = np.linspace(values.min(), values.max(), 2000)
        _grid_probs_high = _gmm.predict_proba(_grid.reshape(-1, 1))[:, _high]
        _crossed = np.flatnonzero(_grid_probs_high >= positive_posterior_cutoff)
        _boundary = _grid[_crossed[0]] if len(_crossed) else values.max()
        return _is_positive, _boundary

    _gmm_results = [_gmm_fit_column(clr[:, _col]) for _col in range(clr.shape[1])]
    positive = np.column_stack([_r[0] for _r in _gmm_results])
    thresholds = np.array([_r[1] for _r in _gmm_results])
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
        return True

    # Captured so downstream cells can reference this name and get a real, trackable
    # marimo dependency edge instead of just sharing a ref to "adata_flt" -- see the
    # marimo-pair race-condition note.
    cmo_assignment_computed = _run_assign_cmo_tag(adata_flt, clr, thresholds, positive, n_pos, adata_cmo, cmo_to_timepoint)
    return (cmo_assignment_computed,)


@app.cell
def cmo_scatter_multiselect_ui(cmo_to_timepoint, mo):
    # Filters which timepoints appear in the threshold-vs-assignment scatter below.
    # Displayed together with the plots (in cmo_threshold_vs_assignment), not here
    # -- this cell only defines it, since its .value can't be read in the same
    # cell that creates it.
    timepoint_multiselect = mo.ui.multiselect(
        options=list(cmo_to_timepoint.values()),
        value=sorted(set(cmo_to_timepoint.values())),
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
    # the scatter's scale doesn't rescale as the timepoint selection changes.
    _x_pad = (_cmo_summary["threshold"].max() - _cmo_summary["threshold"].min()) * 0.1
    _y_pad = (_cmo_summary["n_assigned"].max() - _cmo_summary["n_assigned"].min()) * 0.1
    _x_domain = [_cmo_summary["threshold"].min() - _x_pad, _cmo_summary["threshold"].max() + _x_pad]
    _y_domain = [_cmo_summary["n_assigned"].min() - _y_pad, _cmo_summary["n_assigned"].max() + _y_pad]

    _scatter_points = alt.Chart(_cmo_summary_filtered).mark_circle(size=100, opacity=0.9, stroke="black", strokeWidth=0.5).encode(
        x=alt.X("threshold:Q", title="CLR detection threshold (per-tag GMM boundary)", scale=alt.Scale(domain=_x_domain)),
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
    _timepoint_order = ["d0", "d1", "d2", "d3", "d4", "Negative"]
    _hash_to_timepoint = adata_flt.obs["cmo_tag_scanpy"].map(cmo_to_timepoint).fillna("Negative")
    _timepoint_counts = _hash_to_timepoint.value_counts().reindex(_timepoint_order).reset_index()
    _timepoint_counts.columns = ["timepoint", "n_barcodes"]
    _total_cells = int(_timepoint_counts["n_barcodes"].sum())

    _bar_colors = {t: ec_diff_palette.get(t, ec_diff_palette["Unassigned"]) for t in _timepoint_order}
    _bar_chart = alt.Chart(_timepoint_counts).mark_bar().encode(
        y=alt.Y("timepoint:N", sort=_timepoint_order, title="Timepoint"),
        x=alt.X("n_barcodes:Q", title="Number of barcodes"),
        color=alt.Color(
            "timepoint:N",
            scale=alt.Scale(domain=list(_bar_colors.keys()), range=list(_bar_colors.values())),
            legend=None,
        ),
    )
    _bar_labels = alt.Chart(_timepoint_counts).mark_text(align="left", dx=3, fontSize=9).encode(
        y=alt.Y("timepoint:N", sort=_timepoint_order),
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
    pd,
    positive,
):
    # Are CMO doublets consistently the same two CMOs getting confounded together
    # (e.g. adjacent wells cross-contaminating), or is which two CMOs co-occur
    # essentially random given how often each CMO individually gets falsely
    # triggered? Restricted to doublets with exactly two positive CMOs (n_pos == 2),
    # since a clean pairwise co-occurrence test needs exactly one pair per barcode --
    # barcodes with 3+ positive CMOs don't have a single well-defined pair.
    from statsmodels.stats.multitest import multipletests

    _cmo_names = adata_cmo.var["gene_name"].to_numpy()
    _doublet_mask = (n_pos == 2)
    _positive_doublets = positive[_doublet_mask]
    _pair_idx = np.array([np.flatnonzero(row) for row in _positive_doublets])
    _n_doublets = len(_pair_idx)

    _n_all_doublets = int((n_pos >= 2).sum())
    _n_multi_excluded = _n_all_doublets - _n_doublets
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
    # Permutation p-values are floored at 1/n_perm so BH doesn't choke on exact zeros.
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
    don't reduce to a single CMO pair.

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

    _timepoint_order = ["d0", "d1", "d2", "d3", "d4"]
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
        color=alt.condition(alt.datum.count > _confusion_long["count"].max() / 2, alt.value("black"), alt.value("white")),
    )
    _confusion_chart = (_confusion_heatmap + _confusion_labels).properties(
        title="Doublet pairs by timepoint (raw count)", width=280, height=280,
    ).configure_view(strokeWidth=0)

    # Same matrix, normalized by how many possible CMO pairs exist in each cell --
    # same-timepoint cells have C(n,2) possible pairs, where n is that
    # timepoint's own CMO count (not uniform across timepoints here: 3 for d0,
    # 5 for d1-d4), while cross-timepoint cells have n_a * n_b, so raw counts
    # alone make cross-timepoint pairing look inflated just from having more
    # combinations, and from d0's smaller panel giving it fewer possible pairs.
    _cmos_per_timepoint = {
        t: [c for c in adata_cmo.var["gene_name"] if cmo_to_timepoint[c] == t] for t in _timepoint_order
    }
    _n_possible = pd.DataFrame(0.0, index=_timepoint_order, columns=_timepoint_order)
    for _ta in _timepoint_order:
        for _tb in _timepoint_order:
            _na, _nb = len(_cmos_per_timepoint[_ta]), len(_cmos_per_timepoint[_tb])
            _n_possible.loc[_ta, _tb] = (_na * (_na - 1) / 2) if _ta == _tb else (_na * _nb)

    _rate = _confusion / _n_possible
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
        color=alt.condition(alt.datum.rate > _rate_long["rate"].max() / 2, alt.value("black"), alt.value("white")),
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
    np,
    pd,
    plt,
):
    cmo_assignment_computed  # ran after CMO tags were assigned

    # Number of barcodes assigned (hash_ID) to each CMO. Uses the threshold-based
    # call (cmo_positive_tags_scanpy, restricted to Singlets) rather than the raw
    # argmax tag (cmo_tag_scanpy) -- wtc11_endo_d1_4 has an abnormally elevated
    # background across nearly every cell in this dataset, so it wins the argmax
    # comparison for huge numbers of barcodes that don\'t actually clear its
    # threshold, making cmo_tag_scanpy misleadingly show ~0 barcodes for CMOs
    # like wtc11_endo_d0_1 that are actually correctly Singlet-called via
    # thresholding. Restricted to real CMO tags only -- Negative and Doublet
    # aren\'t attributable to one CMO, so they\'re dropped from this per-CMO view
    # rather than shown as their own bars.
    _cmo_order = list(adata_cmo.var["gene_name"])
    _display_tag = np.where(
        adata_flt.obs["cmo_status_scanpy"] == "Singlet",
        adata_flt.obs["cmo_positive_tags_scanpy"],
        adata_flt.obs["cmo_status_scanpy"],
    )
    _counts = pd.Series(_display_tag).value_counts().reindex(_cmo_order)
    _colors = [ec_diff_palette[cmo_to_timepoint.get(c, "Unassigned")] for c in _counts.index]

    plt.figure(figsize=(10, 5))
    plt.bar(_counts.index, _counts.values, color=_colors)
    plt.xlabel("CMO")
    plt.ylabel("Number of barcodes")
    plt.title("Barcodes assigned per CMO (Singlets only)")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.gcf()
    return


@app.cell(hide_code=True)
def igvf_dataset0_intro(mo):
    mo.md(r"""
    ### External reference: IGVF Dataset0 annotation

    Before running this notebook's own CMO hash classification, it's worth checking an independent reference: `data/10x_multi_5_timepoints_mcginnis/cmo_counts/final_anno_igvf0.tsv`, Chris McGinnis's own pre-computed MULTI-seq sample assignment and cell-type annotation for this dataset ("IGVF Dataset0"). This is useful context given the CMO hash classification below currently produces a heavily skewed result (see CMO_ISSUE.txt) -- this file shows what a healthy, roughly balanced classification looks like for comparison.
    """)
    return


@app.cell(hide_code=True)
def igvf_dataset0_load(data_root_path, pd):
    igvf_dataset0_anno_path = data_root_path / "cmo_counts/final_anno_igvf0.tsv"
    igvf_dataset0_anno = pd.read_csv(igvf_dataset0_anno_path, sep="\t")

    # Match against adata_flt's obs_names, which carry a "_<sample accession>"
    # suffix instead of cellBC's "_CharacterizationMcGinnis_Dataset0" suffix.
    igvf_dataset0_anno["barcode_raw"] = igvf_dataset0_anno["cellBC"].str.split("_", n=1).str[0]
    igvf_dataset0_anno["sample_multiseq_clean"] = igvf_dataset0_anno["sample_multiseq"].str.replace("ansuman-satpathy:", "", regex=False)
    igvf_dataset0_anno["timepoint"] = igvf_dataset0_anno["sample_multiseq_clean"].str.split("_").str[2]
    return (igvf_dataset0_anno,)


@app.cell(hide_code=True)
def igvf_dataset0_cmo_assignment_counts(
    adata_cmo,
    cmo_to_timepoint,
    ec_diff_palette,
    igvf_dataset0_anno,
    plt,
):
    # Same style as plot_cmo_assignment_counts below, but using the external
    # igvf_dataset0_anno (sample_multiseq) assignment instead of this notebook's
    # own CMO hash classification -- a reference point for what a healthy,
    # roughly even per-CMO distribution looks like.
    _cmo_order = list(adata_cmo.var["gene_name"])
    _counts = igvf_dataset0_anno["sample_multiseq_clean"].value_counts().reindex(_cmo_order)
    _colors = [ec_diff_palette[cmo_to_timepoint.get(c, "Unassigned")] for c in _counts.index]

    plt.figure(figsize=(10, 5))
    plt.bar(_counts.index, _counts.values, color=_colors)
    plt.xlabel("CMO (sample_multiseq)")
    plt.ylabel("Number of barcodes")
    plt.title("IGVF Dataset0 annotation: barcodes assigned per CMO")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.gcf()
    return


@app.cell(hide_code=True)
def igvf_dataset0_confusion_matrix(
    adata_flt,
    alt,
    cmo_assignment_computed,
    igvf_dataset0_anno,
    mo,
    np,
    pd,
):
    cmo_assignment_computed  # ran after CMO tags were assigned

    # Confusion matrix: external ground truth (igvf_dataset0's own
    # sample_multiseq timepoint) vs. this notebook's own RAW CMO assignment
    # (timepoint_scanpy). Uses timepoint_scanpy, not rescued_cmo_tag --
    # rescued_cmo_tag folds in RNA-cluster-consensus info for rescued Negatives,
    # which isn\'t a clean test of CMO-quantification agreement on its own.
    # timepoint_scanpy already has explicit "Negative"/"Doublet" string values
    # (not NaN) for non-Singlet calls, so nothing gets silently dropped by
    # pd.crosstab.
    _adata_flt_raw_to_full = pd.Series(adata_flt.obs_names, index=adata_flt.obs_names.str.split("_", n=1).str[0])
    _matched = igvf_dataset0_anno[igvf_dataset0_anno["barcode_raw"].isin(_adata_flt_raw_to_full.index)].copy()
    _matched["obs_name"] = _adata_flt_raw_to_full.loc[_matched["barcode_raw"]].to_numpy()

    _timepoint_order = ["d0", "d1", "d2", "d3", "d4"]
    _sub_obs = adata_flt.obs.loc[_matched["obs_name"]]

    _cm = pd.crosstab(_matched["timepoint"].to_numpy(), _sub_obs["timepoint_scanpy"].astype(str).to_numpy())
    _col_order = _timepoint_order + ["Negative", "Doublet"]
    _cm = _cm.reindex(index=_timepoint_order, columns=_col_order, fill_value=0)

    # Match rate over the FULL matched population: Negative/Doublet-bucketed
    # barcodes count as non-matches, so this is an honest denominator.
    _match_rate = np.diag(_cm[_timepoint_order].to_numpy()).sum() / _cm.to_numpy().sum()

    _cm_long = _cm.reset_index().rename(columns={"row_0": "igvf_dataset0_timepoint"}).melt(
        id_vars="igvf_dataset0_timepoint", var_name="cmo_assignment", value_name="n_barcodes"
    )

    _heatmap = alt.Chart(_cm_long).mark_rect().encode(
        x=alt.X("cmo_assignment:N", sort=_col_order, title="this notebook\'s own CMO assignment (timepoint_scanpy)"),
        y=alt.Y("igvf_dataset0_timepoint:N", sort=_timepoint_order, title="IGVF Dataset0 timepoint"),
        color=alt.Color("n_barcodes:Q", scale=alt.Scale(scheme="cividis"), title="Barcodes"),
        tooltip=["igvf_dataset0_timepoint", "cmo_assignment", "n_barcodes"],
    )
    _labels = alt.Chart(_cm_long).mark_text(fontSize=11).encode(
        x=alt.X("cmo_assignment:N", sort=_col_order),
        y=alt.Y("igvf_dataset0_timepoint:N", sort=_timepoint_order),
        text=alt.Text("n_barcodes:Q", format=","),
        color=alt.condition(alt.datum.n_barcodes > _cm_long["n_barcodes"].max() / 2, alt.value("black"), alt.value("white")),
    )
    _chart = (_heatmap + _labels).properties(
        title=f"IGVF Dataset0 timepoint vs. this notebook\'s raw CMO assignment (match rate: {_match_rate:.1%}, n={_cm.to_numpy().sum():,})",
        width=460, height=340,
    ).configure_view(strokeWidth=0)

    mo.vstack([_chart, mo.md(_cm.to_markdown())])
    return


@app.cell(hide_code=True)
def igvf_dataset0_cmo_signal_check(
    adata_cmo,
    adata_flt,
    alt,
    clr,
    igvf_dataset0_anno,
    mo,
    np,
    pd,
    thresholds,
):
    # Does this notebook's own CLR-normalized CMO signal for wtc11_endo_d0_1
    # actually separate the barcodes the external IGVF Dataset0 annotation
    # assigns to that exact sample from everyone else?
    from scipy.stats import gaussian_kde as _gaussian_kde_igvf0

    _target_cmo = "wtc11_endo_d0_1"
    _col_idx = list(adata_cmo.var["gene_name"]).index(_target_cmo)

    _adata_flt_raw_to_full = pd.Series(adata_flt.obs_names, index=adata_flt.obs_names.str.split("_", n=1).str[0])
    _target_barcodes_raw = igvf_dataset0_anno.loc[
        igvf_dataset0_anno["sample_multiseq_clean"] == _target_cmo, "barcode_raw"
    ]
    _matched_full = _adata_flt_raw_to_full.reindex(_target_barcodes_raw).dropna()
    _in_group_mask = adata_cmo.obs_names.isin(_matched_full.to_numpy())

    _vals_in = clr[_in_group_mask, _col_idx]
    _vals_out = clr[~_in_group_mask, _col_idx]

    _grid = np.linspace(0, max(_vals_in.max(), _vals_out.max()), 300)
    _density_df = pd.DataFrame({
        "clr_value": np.tile(_grid, 2),
        "density": np.concatenate([
            _gaussian_kde_igvf0(_vals_in)(_grid),
            _gaussian_kde_igvf0(_vals_out)(_grid),
        ]),
        "group": [f"IGVF0 {_target_cmo} (n={int(_in_group_mask.sum()):,})"] * len(_grid)
            + [f"Rest (n={int((~_in_group_mask).sum()):,})"] * len(_grid),
    })

    _density_chart = alt.Chart(_density_df).mark_line(interpolate="monotone").encode(
        x=alt.X("clr_value:Q", title=f"CLR-normalized {_target_cmo} signal"),
        y=alt.Y("density:Q", title="Density"),
        color=alt.Color("group:N", title=None),
    )
    _threshold_rule = alt.Chart(pd.DataFrame({"x": [thresholds[_col_idx]]})).mark_rule(
        strokeDash=[4, 4], color="black"
    ).encode(x="x:Q")

    _chart = (_density_chart + _threshold_rule).properties(
        title=f"CLR-normalized {_target_cmo} signal: IGVF0-assigned vs. rest (dashed = this notebook's own detection threshold)",
        width=550, height=350,
    ).configure_view(strokeWidth=0)

    mo.vstack([
        _chart,
        mo.md(
            f"**{_target_cmo}** -- threshold = {thresholds[_col_idx]:.3f} | "
            f"IGVF0-assigned: median = {np.median(_vals_in):.3f}, mean = {_vals_in.mean():.3f} | "
            f"rest: median = {np.median(_vals_out):.3f}, mean = {_vals_out.mean():.3f}"
        ),
    ])
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
    # "not singlet" isn't one thing for CMO hashing -- it's either a Doublet call
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
    # and the lowest predicted-doublet score (Scrublet doesn't expose it directly here).
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
    # this cell only defines it, since its .value can't be read in the same cell
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

    # Subsample for the plot itself -- two full-size (~32k-row) raw scatter specs
    # in one cell output is too large for marimo to render (no aggregating
    # transform here for VegaFusion to pre-reduce, unlike the binned histograms).
    # Both panels share the same subsampled points so they stay directly comparable.
    _max_points = 8000
    _n_before_subsample = len(_umap_df)
    if _n_before_subsample > _max_points:
        _umap_df = _umap_df.sample(n=_max_points, random_state=0)

    # Fixed axis domains from the FULL (unfiltered) UMAP, with a little padding, so
    # the scale doesn't rescale as the status selection changes.
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
def final_round_intro(mo):
    mo.md(r"""
    # Final round of quality control
    """)
    return


@app.cell(hide_code=True)
def create_qc_passing_adata(
    adata_flt,
    alt,
    cmo_assignment_computed,
    ec_diff_palette,
    igvf_dataset0_anno,
    pd,
):
    cmo_assignment_computed  # ran after CMO tags were assigned

    # Direct per-cell QC filter only, no cluster-based exclusion here (clustering
    # happens after this filter, not before it): keep barcodes that are both a
    # Scrublet singlet and a CMO singlet.
    _is_scrublet_singlet = ~adata_flt.obs["scrublet_predicted_doublet"]
    _is_cmo_singlet = adata_flt.obs["cmo_status_scanpy"] == "Singlet"
    _qc_pass_mask = _is_scrublet_singlet & _is_cmo_singlet

    adata_clean = adata_flt[_qc_pass_mask].copy()

    _n_total = adata_flt.n_obs
    _n_pass = adata_clean.n_obs
    print(f"Scrublet singlet: {int(_is_scrublet_singlet.sum()):,} / {_n_total:,}")
    print(f"CMO singlet: {int(_is_cmo_singlet.sum()):,} / {_n_total:,}")
    print(f"Both (kept in adata_clean): {_n_pass:,} / {_n_total:,} ({_n_pass / _n_total:.1%})")

    _timepoint_order = ["d0", "d1", "d2", "d3", "d4"]
    _tag_counts = adata_clean.obs["timepoint_scanpy"].value_counts().reindex(_timepoint_order).reset_index()
    _tag_counts.columns = ["timepoint", "n_barcodes"]

    _bar_colors = {t: ec_diff_palette.get(t, ec_diff_palette["Unassigned"]) for t in _timepoint_order}
    _bar_chart = alt.Chart(_tag_counts).mark_bar().encode(
        y=alt.Y("timepoint:N", sort=_timepoint_order, title="Timepoint"),
        x=alt.X("n_barcodes:Q", title="Number of barcodes"),
        color=alt.Color(
            "timepoint:N",
            scale=alt.Scale(domain=list(_bar_colors), range=list(_bar_colors.values())),
            legend=None,
        ),
    )
    _bar_labels = alt.Chart(_tag_counts).mark_text(align="left", dx=3, fontSize=9).encode(
        y=alt.Y("timepoint:N", sort=_timepoint_order),
        x=alt.X("n_barcodes:Q"),
        text=alt.Text("n_barcodes:Q", format=","),
    )
    _timepoint_chart = (_bar_chart + _bar_labels).properties(
        title=["Barcodes per timepoint (adata_clean: Scrublet singlet + CMO singlet)", f"n = {_n_pass:,} barcodes total"],
        width=480, height=300,
    ).configure_view(strokeWidth=0)

    # Concordance with IGVF0's independent annotation, restricted to adata_clean.
    _raw_to_full = pd.Series(adata_clean.obs_names, index=adata_clean.obs_names.str.split("_", n=1).str[0])
    _matched = igvf_dataset0_anno[igvf_dataset0_anno["barcode_raw"].isin(_raw_to_full.index)].copy()
    _matched["obs_name"] = _raw_to_full.loc[_matched["barcode_raw"]].to_numpy()
    _matched["our_timepoint"] = adata_clean.obs.loc[_matched["obs_name"], "timepoint_scanpy"].to_numpy()
    _agree = _matched["timepoint"] == _matched["our_timepoint"]

    _n_igvf0_total = len(igvf_dataset0_anno)
    _n_igvf0_matched = len(_matched)
    print(f"IGVF0 barcodes in common with adata_clean: {_n_igvf0_matched:,} / {_n_igvf0_total:,} ({_n_igvf0_matched / _n_igvf0_total:.1%})")

    # Why the rest were lost: not present in adata_flt at all (failed the lenient
    # QC / counts cutoff, or no barcode overlap), vs. present but excluded by the
    # Scrublet-singlet / CMO-singlet filter.
    _raw_to_flt = pd.Series(adata_flt.obs_names, index=adata_flt.obs_names.str.split("_", n=1).str[0])
    _in_adata_flt = igvf_dataset0_anno["barcode_raw"].isin(_raw_to_flt.index)
    _n_not_in_flt = int((~_in_adata_flt).sum())

    _present = igvf_dataset0_anno[_in_adata_flt].copy()
    _present["obs_name"] = _raw_to_flt.loc[_present["barcode_raw"]].to_numpy()
    _present_obs = adata_flt.obs.loc[_present["obs_name"]]
    _present["scrublet_doublet"] = _present_obs["scrublet_predicted_doublet"].to_numpy()
    _present["cmo_singlet"] = (_present_obs["cmo_status_scanpy"] == "Singlet").to_numpy()
    _lost_present = _present[~(~_present["scrublet_doublet"] & _present["cmo_singlet"])]

    print(f"\nOf the {_n_igvf0_total - _n_igvf0_matched:,} lost barcodes:")
    print(f"  Not in adata_flt at all (failed lenient QC / counts cutoff, or no overlap): {_n_not_in_flt:,}")
    print(f"  In adata_flt but Scrublet-predicted doublet (CMO singlet otherwise): {int((_lost_present['scrublet_doublet'] & _lost_present['cmo_singlet']).sum()):,}")
    print(f"  In adata_flt but CMO not singlet (Scrublet singlet otherwise): {int((~_lost_present['scrublet_doublet'] & ~_lost_present['cmo_singlet']).sum()):,}")
    print(f"  In adata_flt but failed both filters: {int((_lost_present['scrublet_doublet'] & ~_lost_present['cmo_singlet']).sum()):,}")

    _igvf0_concordance = pd.crosstab(
        _matched["timepoint"], _matched["our_timepoint"]
    ).reindex(index=_timepoint_order, columns=_timepoint_order, fill_value=0)

    print(f"IGVF0 concordance: {int(_agree.sum()):,} / {len(_matched):,} matched barcodes agree on timepoint ({_agree.mean():.1%}).")
    print(_igvf0_concordance.to_string())

    _timepoint_chart

    return (adata_clean,)


@app.cell
def cluster_adata_clean(adata_clean, sc, scclr):
    def _run_clustering_pipeline(adata):
        scclr.pp.pflog(adata, target="auto")
        sc.pp.highly_variable_genes(adata, layer="pflog", n_top_genes=2000)
        _ncomps = 50
        _ncv = 2 * _ncomps + 1
        scclr.tl.pca(adata, n_comps=_ncomps, ncv=_ncv)
        sc.pp.neighbors(adata, random_state=0)
        sc.tl.leiden(adata, flavor="igraph", resolution=0.25, n_iterations=2, random_state=0)
        sc.tl.umap(adata)
        return True

    adata_clean_clustering_computed = _run_clustering_pipeline(adata_clean)

    return (adata_clean_clustering_computed,)


@app.cell
def umap_adata_clean(
    adata_clean,
    adata_clean_clustering_computed,
    ec_diff_palette,
    sc,
):
    adata_clean_clustering_computed  # ran after normalization, HVG, PCA, neighbors, leiden, and UMAP

    sc.pl.umap(
        adata_clean,
        color=["leiden", "pct_counts_mt"],
        cmap="cividis",
        ncols=2,
        size=5,
        title=["UMAP colored by leiden cluster", "UMAP colored by %mito"],
    )

    # Separate call: timepoint_scanpy needs the timepoint-specific ec_diff_palette,
    # not scanpy's default categorical palette. adata_clean is already restricted
    # to CMO singlets, so timepoint_scanpy is a clean d0-d4 label here (no
    # Negative/Doublet values to worry about).
    sc.pl.umap(
        adata_clean,
        color="timepoint_scanpy",
        palette=ec_diff_palette,
        size=5,
        title="UMAP colored by timepoint",
    )

    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Summary of quality control metrics per cluster
    """)
    return


@app.cell
def cluster_summary_adata_clean(
    adata_clean,
    adata_clean_clustering_computed,
    cmo_assignment_computed,
):
    cmo_assignment_computed  # ran after CMO tags were assigned
    adata_clean_clustering_computed  # ran after leiden clustering on adata_clean

    # Per-cluster overview: composition of the new clusters computed directly on
    # the already-filtered adata_clean (Scrublet-singlet + CMO-singlet
    # survivors). No more doublet-rate columns here since adata_clean has none by
    # construction -- this is now a plain per-cluster QC + timepoint-composition
    # summary, used to spot confounders in the clustering itself.
    cluster_summary_adata_clean = adata_clean.obs.groupby("leiden", observed=True).agg(
        n_cells=("leiden", "size"),
        median_mito=("pct_counts_mt", "median"),
        median_counts=("total_counts", "median"),
    ).round(1)

    _timepoint_order = ["d0", "d1", "d2", "d3", "d4"]
    _timepoint_pct = (
        adata_clean.obs.groupby("leiden", observed=True)["timepoint_scanpy"]
        .value_counts(normalize=True)
        .unstack(fill_value=0)
        .reindex(columns=_timepoint_order, fill_value=0) * 100
    ).round(1)
    _timepoint_pct.columns = [f"pct_{_tp}" for _tp in _timepoint_pct.columns]

    cluster_summary_adata_clean = cluster_summary_adata_clean.join(_timepoint_pct)

    cluster_summary_adata_clean[
        ["n_cells", "median_counts", "median_mito", *_timepoint_pct.columns]
    ].sort_values("n_cells", ascending=False)

    return


@app.cell(hide_code=True)
def mixed_timepoint_clusters_investigation(
    adata_clean,
    adata_clean_clustering_computed,
    mo,
    pd,
    sc,
):
    adata_clean_clustering_computed  # ran after leiden clustering on adata_clean

    # Clusters 3 and 6 are the least timepoint-pure clusters in
    # cluster_summary_adata_clean (3: 16.6/9.9/18.7/42.9/12.0% across d0-d4; 6:
    # 6.3/54.4/27.2/11.3/0.8%, also the lowest median_counts of any cluster).
    # Check whether either reflects a proliferation/cell-cycle-driven clustering
    # artifact (cells clustering by cycle state rather than lineage/timepoint)
    # instead of a distinct identity, the same check run on the previous
    # clustering's mixed clusters.
    _view = adata_clean.copy()
    sc.tl.rank_genes_groups(_view, groupby="leiden", groups=["3", "6"], reference="rest", method="wilcoxon", layer="pflog", use_raw=False)

    _top_n = 15
    _marker_dfs = {}
    for _grp in ["3", "6"]:
        _df = sc.get.rank_genes_groups_df(_view, group=_grp).reset_index(drop=True)
        _df["gene_symbol"] = _df["names"].map(adata_clean.var["gene_symbol"])
        _df["rank"] = _df.index + 1
        _marker_dfs[_grp] = _df

    mixed_timepoint_top_markers = pd.DataFrame({_grp: _df.head(_top_n)["gene_symbol"].tolist() for _grp, _df in _marker_dfs.items()})
    mixed_timepoint_top_markers.index = range(1, _top_n + 1)
    mixed_timepoint_top_markers.index.name = "rank"

    # Canonical S-phase and G2/M proliferation markers.
    _cell_cycle_genes = [
        "MKI67", "TOP2A", "CENPF", "MCM2", "MCM4", "MCM6", "CLSPN", "DTL",
        "PCNA", "CCNB1", "CCNB2", "CDK1", "AURKA", "AURKB", "BUB1", "TYMS",
    ]
    _n_genes = adata_clean.n_vars

    _rows = []
    for _gene in _cell_cycle_genes:
        _row = {"gene": _gene}
        for _grp, _df in _marker_dfs.items():
            _match = _df.loc[_df["gene_symbol"] == _gene]
            if len(_match):
                _row[f"cluster_{_grp}_rank"] = f"{int(_match['rank'].iloc[0]):,} / {_n_genes:,}"
                _row[f"cluster_{_grp}_logfc"] = round(float(_match["logfoldchanges"].iloc[0]), 2)
            else:
                _row[f"cluster_{_grp}_rank"] = "not found"
                _row[f"cluster_{_grp}_logfc"] = None
        _rows.append(_row)

    mixed_timepoint_cell_cycle_ranks = pd.DataFrame(_rows).set_index("gene")

    mo.vstack([
        mo.md("**Top 15 markers, clusters 3 and 6 vs. rest (Wilcoxon, `pflog` layer):**"),
        mixed_timepoint_top_markers,
        mo.md(
            "**Rank and logFC of canonical proliferation/cell-cycle genes in each cluster\'s full ranking "
            f"(out of {_n_genes:,} genes; low rank + positive logFC would support a cell-cycle-driven cluster):**"
        ),
        mixed_timepoint_cell_cycle_ranks,
    ])

    return


@app.cell
def same_timepoint_mito_split_investigation(
    adata_clean,
    adata_clean_clustering_computed,
    mo,
    pd,
    sc,
):
    adata_clean_clustering_computed  # ran after leiden clustering on adata_clean

    # Within the same timepoint, cluster_summary_adata_clean shows a striking
    # mito split: d2 has cluster 1 (0.9% median mito) vs cluster 2 (14.3%), and
    # d4 has cluster 7 (1.4%) vs cluster 4 (13.9%). Compare each low-mito vs
    # high-mito pair directly (not vs "rest") to see whether this is genuine
    # biological heterogeneity within the timepoint or a residual technical
    # (stress/dying-nuclei) confound.
    _pairs = {"d2": ("1", "2"), "d4": ("7", "4")}  # timepoint: (low-mito cluster, high-mito cluster)

    _marker_dfs = {}
    for _tp, (_lo, _hi) in _pairs.items():
        _view = adata_clean[adata_clean.obs["leiden"].isin([_lo, _hi])].copy()
        sc.tl.rank_genes_groups(_view, groupby="leiden", groups=[_hi], reference=_lo, method="wilcoxon", layer="pflog", use_raw=False)
        _df = sc.get.rank_genes_groups_df(_view, group=_hi).reset_index(drop=True)
        _df["gene_symbol"] = _df["names"].map(adata_clean.var["gene_symbol"])
        _df["rank"] = _df.index + 1
        _marker_dfs[_tp] = _df

    _top_n = 15
    mito_split_top_markers = pd.DataFrame({
        f"{_tp} (cl.{_hi} vs cl.{_lo})": _marker_dfs[_tp].head(_top_n)["gene_symbol"].tolist()
        for _tp, (_lo, _hi) in _pairs.items()
    })
    mito_split_top_markers.index = range(1, _top_n + 1)
    mito_split_top_markers.index.name = "rank"

    # Mito genes dominate the raw top-15 almost tautologically -- these clusters
    # were picked BY %mito, so of course mito-encoded genes separate them. The
    # real question is what distinguishes them beyond that, so also show the top
    # markers with mito-encoded and ribosomal genes excluded.
    mito_split_top_markers_non_mito = pd.DataFrame({
        f"{_tp} (cl.{_hi} vs cl.{_lo})": (
            _marker_dfs[_tp][~_marker_dfs[_tp]["gene_symbol"].str.startswith(("MT-", "RPL", "RPS"), na=False)]
            .head(_top_n)["gene_symbol"].tolist()
        )
        for _tp, (_lo, _hi) in _pairs.items()
    })
    mito_split_top_markers_non_mito.index = range(1, _top_n + 1)
    mito_split_top_markers_non_mito.index.name = "rank"

    # Mito-encoded and heat-shock/ER-stress chaperone genes, the same panel used
    # for the pre-filter clustering's earlier high-mito investigation.
    _stress_genes = [
        "MT-ND5", "MT-CO3", "MT-ND2", "MT-ND1", "MT-ND4", "MT-CO2", "MT-ND3",
        "MT-ATP6", "MT-CO1", "MT-CYB", "MT-ND4L", "HSP90AB1", "HSP90B1", "HSPA5",
        "CANX", "RPLP1",
    ]
    _n_genes = adata_clean.n_vars

    _rows = []
    for _gene in _stress_genes:
        _row = {"gene": _gene}
        for _tp, _df in _marker_dfs.items():
            _match = _df.loc[_df["gene_symbol"] == _gene]
            if len(_match):
                _row[f"{_tp}_rank"] = f"{int(_match['rank'].iloc[0]):,} / {_n_genes:,}"
                _row[f"{_tp}_logfc"] = round(float(_match["logfoldchanges"].iloc[0]), 2)
            else:
                _row[f"{_tp}_rank"] = "not found"
                _row[f"{_tp}_logfc"] = None
        _rows.append(_row)

    mito_split_stress_gene_ranks = pd.DataFrame(_rows).set_index("gene")

    mo.vstack([
        mo.md("**Top 15 markers, high-mito cluster vs. low-mito cluster within the same timepoint (Wilcoxon, `pflog` layer):**"),
        mito_split_top_markers,
        mo.md(
            "**Same comparison, mito-encoded and ribosomal genes excluded** (those dominate the raw "
            "list almost tautologically, since these clusters were picked by %mito in the first place -- "
            "this shows what else distinguishes them):"
        ),
        mito_split_top_markers_non_mito,
        mo.md(
            "**Rank and logFC of mito-encoded and heat-shock/ER-stress chaperone genes "
            f"(out of {_n_genes:,} genes; low rank + positive logFC in the high-mito cluster would support a stress/dying-nuclei explanation):**"
        ),
        mito_split_stress_gene_ranks,
    ])

    return


@app.cell
def cluster_replicate_composition_check(
    adata_clean,
    adata_clean_clustering_computed,
    pd,
):
    adata_clean_clustering_computed  # ran after leiden clustering on adata_clean

    # If the mito split within a timepoint were purely intrinsic biology, both
    # clusters should draw roughly evenly from all of that timepoint's
    # biological replicate wells. A strong skew toward specific wells instead
    # points at a well-specific technical effect (e.g. differential ambient
    # contamination from uneven lysis/cell health at harvest), since replicates
    # are supposed to be interchangeable copies of the same condition.
    from scipy.stats import chi2_contingency

    _pairs = {"d2": ("1", "2", ["wtc11_endo_d2_1", "wtc11_endo_d2_2", "wtc11_endo_d2_4", "wtc11_endo_d2_5"]),
              "d4": ("7", "4", ["wtc11_endo_d4_1", "wtc11_endo_d4_3", "wtc11_endo_d4_4", "wtc11_endo_d4_5"])}

    _results = {}
    for _tp, (_lo, _hi, _tags) in _pairs.items():
        _sub = adata_clean.obs[adata_clean.obs["leiden"].isin([_lo, _hi])]
        _sub = _sub[_sub["cmo_positive_tags_scanpy"].isin(_tags)]
        _ct = pd.crosstab(_sub["leiden"], _sub["cmo_positive_tags_scanpy"])
        _chi2, _p, _dof, _ = chi2_contingency(_ct)
        _pct = (_ct.div(_ct.sum(axis=1), axis=0) * 100).round(1)
        _results[_tp] = {"counts": _ct, "pct": _pct, "chi2": _chi2, "p": _p}
        print(f"{_tp} (cluster {_lo} vs {_hi}): chi2={_chi2:.2f}, p={_p:.2e}")
        print(_pct)
        print()

    replicate_composition_by_cluster = _results

    return


@app.cell(hide_code=True)
def cluster_confounder_conclusion(mo):
    mo.md(r"""
    ### Cluster confounders in the new clustering: cell-cycle artifact, real transitional state, and real (not stressed) endothelial biology

    Re-clustering `adata_clean` after the GMM-calibrated CMO threshold change produced different clusters, so the earlier "clusters 2 and 3" investigation no longer targeted the right ones. Two independent patterns stood out in `cluster_summary_adata_clean`: mixed-timepoint clusters (3, 6), and a striking within-timepoint mito split (clusters 1 vs 2 in d2, 7 vs 4 in d4).

    **Mixed-timepoint clusters (`mixed_timepoint_clusters_investigation`): same pattern as before, just renumbered.** Cluster 6 is a cell-cycle/S-phase artifact -- `CENPF` (rank 1 of 62,757), `MKI67` (rank 2), `DTL` (rank 5), `CLSPN` (rank 30), `MCM4` (rank 64) all strongly elevated, while G2/M mitotic genes (`CCNB1`, `CCNB2`, `CDK1`, `AURKA`, `AURKB`, `BUB1`) are depleted. Cluster 3 shows no cell-cycle signature at all and is instead defined by axon-guidance/neuronal-migration genes (`UNC5C`, `SLIT3`, `KIF26B`, `RALYL`, `PCDH7`, `ROBO2`) -- the same tip-cell-like transitional identity found previously.

    **Same-timepoint mito split (`same_timepoint_mito_split_investigation`): initially looked like a stress artifact, but checking the actual identity markers overturns that.** The raw top markers are dominated by mito-encoded genes in both pairs (d2: cluster 2 vs 1; d4: cluster 4 vs 7) -- expected and not meaningful on its own, since these clusters were picked BY %mito in the first place. Heat-shock/ER-stress chaperone genes (`HSP90AB1`, `HSP90B1`, `HSPA5`, `CANX`, `RPLP1`) are also elevated in both pairs, which looked like a stress signature at first.

    But excluding mito-encoded and ribosomal genes reveals real, coherent cell-identity biology underneath. **Cluster 4 (d4)** is marked by unambiguous endothelial genes: `KDR` (VEGFR2), `CD93`, `PLVAP`, `ESM1`, `PLXND1`, `GJA1` -- and `SERPINH1` (HSP47, a collagen-specific chaperone) alongside real ECM genes explains the "ER stress" signal as the protein-folding load of active collagen/ECM secretion, not death. This is the same conclusion the pre-filter clustering already reached for its own moderately-elevated-mito clusters 6 and 7 (real metabolic/secretory endothelial biology, not debris), which were correctly kept. **Cluster 2 (d2)** shows a related pattern: `ROBO2`, `UNC5B`, `UNC5C`, `PCDH10` -- the same Slit-Robo/protocadherin axon-guidance family already established as real tip-cell-like transitional biology for cluster 3.

    **Conclusion:** neither cluster 2 nor cluster 4 should be excluded. Both show real, coherent, differentiation-relevant biology once the tautological mito-gene signal is set aside, matching this dataset's own precedent from the pre-filter investigation. Practically, excluding them would also have been a large and disproportionate cut: cluster 4 alone is 1,565 of roughly 2,144 d4 cells (about 73% of the timepoint), and cluster 2 is 2,872 of `adata_clean`'s 12,418 cells (23.1%) -- losing that much of specific timepoints on the strength of an incomplete marker check (mito/stress genes only) would have been a mistake.
    """)

    return


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

    `pass_strict_qc` marks the same barcodes kept in `adata_clean`: a Scrublet singlet AND a CMO-hashing Singlet (see `create_qc_passing_adata`). Unlike channel1/channel2, no cluster-based exclusion criterion is applied here -- clustering happens after this filter rather than before it (see `final_round_intro`), so there is no doublet-dominated-cluster concept to apply at this stage.
    """)

    return


@app.cell(hide_code=True)
def strict_qc_cascade(
    adata_clean,
    adata_clean_clustering_computed,
    adata_flt,
    mo,
):
    adata_clean_clustering_computed  # ran after adata_clean was built and clustered

    # pass_strict_qc marks the same barcodes kept in adata_clean, so downstream
    # consumers that need this as a column on adata_flt.obs (write_qc_annotations_tsv,
    # the ATAC whitelist) stay in sync with adata_clean by construction rather than
    # recomputing the filter a second time.
    adata_flt.obs["pass_strict_qc"] = adata_flt.obs_names.isin(adata_clean.obs_names)
    strict_qc_computed = True

    _n_total = adata_flt.n_obs
    _n_pass = int(adata_flt.obs["pass_strict_qc"].sum())

    mo.md(f"""
    **Net:** {_n_pass:,} of {_n_total:,} lenient-QC barcodes ({_n_pass / _n_total:.1%}) survive strict QC (`pass_strict_qc`), matching `adata_clean`.
    """)

    return (strict_qc_computed,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Saving it to file
    """)
    return


@app.cell(hide_code=True)
def write_qc_annotations_tsv(
    adata_flt,
    outdir_root,
    project_root,
    strict_qc_computed,
):
    strict_qc_computed  # ran after pass_strict_qc was set

    _outfile = outdir_root / "adata_flt_qc_annotations.tsv"
    _outfile.parent.mkdir(parents=True, exist_ok=True)
    adata_flt.obs.to_csv(_outfile, sep="\t")

    # Show only the path relative to the repo, not the full local filesystem path
    _outfile.relative_to(project_root)

    return


@app.cell(hide_code=True)
def filtering_strategy_conclusion(
    adata_clean,
    adata_clean_clustering_computed,
    adata_flt,
    igvf_dataset0_anno,
    mo,
    pd,
    strict_qc_computed,
):
    strict_qc_computed  # ran after pass_strict_qc was set
    adata_clean_clustering_computed  # ran after leiden clustering on adata_clean

    _n_lenient = adata_flt.n_obs
    _n_pass = adata_clean.n_obs
    _pct_pass = _n_pass / _n_lenient
    _n_scrublet_doublet = int(adata_flt.obs["scrublet_predicted_doublet"].sum())
    _n_cmo_not_singlet = int((adata_flt.obs["cmo_status_scanpy"] != "Singlet").sum())

    _raw_to_clean = pd.Series(adata_clean.obs_names, index=adata_clean.obs_names.str.split("_", n=1).str[0])
    _igvf0_in_clean = igvf_dataset0_anno["barcode_raw"].isin(_raw_to_clean.index)
    _pct_igvf0_overlap = _igvf0_in_clean.mean()

    _matched = igvf_dataset0_anno[_igvf0_in_clean].copy()
    _matched["obs_name"] = _raw_to_clean.loc[_matched["barcode_raw"]].to_numpy()
    _matched["our_timepoint"] = adata_clean.obs.loc[_matched["obs_name"], "timepoint_scanpy"].to_numpy()
    _pct_igvf0_agree = (_matched["timepoint"] == _matched["our_timepoint"]).mean()

    _n_clusters = adata_clean.obs["leiden"].nunique()

    mo.md(f"""
    ### Filtering strategy: conclusions

    Starting from {_n_lenient:,} barcodes surviving the lenient QC filter (median + 2.5 MADs mito cutoff, applied pre-clustering in `mask_filter_cells`), strict QC keeps only barcodes that are both a Scrublet singlet and a CMO-hashing Singlet (`create_qc_passing_adata`), leaving **{_n_pass:,} of {_n_lenient:,} barcodes ({_pct_pass:.1%})** as `adata_clean` / `pass_strict_qc`. This drops {_n_scrublet_doublet:,} Scrublet-predicted doublets and {_n_cmo_not_singlet:,} barcodes that are CMO-hashing Negative or Doublet (with overlap between the two). Unlike this dataset's own earlier approach, no cluster-based exclusion criterion is applied at this stage: clustering happens after this filter rather than before it, so there is no doublet-dominated-cluster concept left to apply, by construction.

    This direct per-cell filter already agrees well with the external `final_anno_igvf0.tsv` annotation (`create_qc_passing_adata`): {_pct_igvf0_overlap:.1%} of IGVF0's barcodes are present in `adata_clean`, and of those, {_pct_igvf0_agree:.1%} agree on timepoint. Since every barcode in `adata_clean` is a direct CMO Singlet by construction, `timepoint_scanpy` is already a real, direct d0-d4 label for every survivor -- there is no cluster-consensus rescue step for this dataset, unlike channel1/channel2 (that strategy was tried and retired; see the rescue-strategy discussion above).

    Clustering `adata_clean` from scratch (`cluster_adata_clean`) gives {_n_clusters} clusters, five of which are each over 90% a single timepoint (`doublet_overview_adata_clean`). The remaining two, clusters 2 and 3, have mixed timepoint composition but for different reasons (`clusters_2_3_conclusion`): cluster 2 is a cell-cycle/S-phase clustering artifact, not a distinct identity, while cluster 3 is real, timepoint-spanning transitional biology. Neither is excluded from `adata_clean` at this stage.
    """)

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
    chrom_dict = pd.read_csv(project_root / "annotations/GRCh38_EBV.chrom.sizes.no.alt.tsv", sep="\t", header=None, names=["chr", "size"])
    chrom_dict = chrom_dict.set_index("chr")["size"].to_dict()
    return (chrom_dict,)


@app.cell
def atac_import_intro(mo):
    mo.md(r"""
    ### Import ATAC fragments

    We import fragments restricted to the RNA-side `pass_strict_qc` whitelist rather than re-filtering on ATAC metrics separately. We trust the RNA-based QC (doublet clusters, high mito, Scrublet, CMO doublets) as the primary cell call, and only compute ATAC QC metrics for visibility.
    """)
    return


@app.cell
def load_strict_qc_whitelist(outdir_root, pd):
    # See atac_import_intro above.
    _qc_annotations = pd.read_csv(
        outdir_root / "adata_flt_qc_annotations.tsv", sep="\t", index_col=0
    )
    assigned_cells = _qc_annotations.index[_qc_annotations["pass_strict_qc"]].to_list()
    len(assigned_cells)
    return (assigned_cells,)


@app.cell
def import_atac_fragments(
    assigned_cells,
    chrom_dict,
    data_root_path,
    igvf_gencode_gtf_path,
    snap,
):
    atac_adata = snap.pp.import_fragments(
        data_root_path / "atac/fragments/IGVFFI6722CCZW.bed.gz",
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


@app.cell
def atac_qc_plot_intro(mo):
    mo.md(r"""
    ### ATAC QC, for visibility only

    Unique fragments vs. TSSE, with fixed thresholds (`atac_qc_thresholds`) shown for reference. We are not applying this as an additional filter, see the note below, since the barcode set already inherits the RNA `pass_strict_qc` whitelist at import.
    """)
    return


@app.cell
def atac_qc_thresholds(atac_adata):
    df = atac_adata.obs[["n_fragment", "tsse"]].copy()

    # Fixed thresholds for visibility only, see the note below: we don't actually
    # filter on these, the ATAC barcode set already inherits the RNA pass_strict_qc
    # whitelist at import time.
    min_frag = 1000
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
        (plot_df["n_fragment"] <= 70_000)
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
def atac_no_filtering_note(mo):
    mo.md("""
    Intentionally not calling `snap.pp.filter_cells` here, the barcode set already reflects the RNA `pass_strict_qc` whitelist used at import (see `atac_import_intro` above), so we trust the RNA-based QC rather than layering an additional ATAC-metric-based filter on top. The QC plot above is for visibility only.
    """)
    return


@app.cell
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


@app.cell
def map_rescued_cmo_tag_to_atac(
    adata_flt,
    atac_adata,
    atac_spectral_umap_leiden_computed,
    cmo_assignment_computed,
):
    cmo_assignment_computed  # ran after CMO tags were assigned, for timepoint_scanpy below
    # Forces a deterministic write order into atac_adata.obs: both this cell
    # and atac_spectral_umap_leiden mutate atac_adata.obs in place (adding
    # "rescued_cmo_tag" and "leiden" respectively) with no data dependency
    # between them, so without this, which column lands first -- and thus
    # atac_adata.obs's final column order -- depends on arbitrary execution
    # order rather than being reproducible.
    atac_spectral_umap_leiden_computed  # ran after spectral embedding, UMAP, KNN, and Leiden clustering

    def _map_rescued_cmo_tag_to_atac(atac_adata, adata_flt):
        # No cluster-consensus rescue for this dataset (rescue was retired, see
        # rescue_negative_tags_intro history) -- every barcode here already
        # passed pass_strict_qc, i.e. is a direct CMO Singlet, so timepoint_scanpy
        # is already a real, direct d0-d4 label for all of them.
        # remove_unused_categories() matters here: timepoint_scanpy's dtype
        # carries "Negative"/"Doublet" categories from the full adata_flt, even
        # though zero ATAC barcodes have those values -- sc.pl.umap builds its
        # palette lookup from ALL declared categories, not just observed ones,
        # and KeyErrors on any category missing from ec_diff_palette otherwise.
        atac_adata.obs["rescued_cmo_tag"] = (
            adata_flt.obs.loc[atac_adata.obs_names, "timepoint_scanpy"]
            .astype("category")
            .cat.remove_unused_categories()
        )
        return True

    rescued_cmo_tag_mapped = _map_rescued_cmo_tag_to_atac(atac_adata, adata_flt)

    return (rescued_cmo_tag_mapped,)


@app.cell
def atac_timepoint_counts(atac_adata, rescued_cmo_tag_mapped):
    rescued_cmo_tag_mapped  # ran after rescued_cmo_tag was mapped onto atac_adata

    atac_adata.obs["rescued_cmo_tag"].value_counts()
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
    okabe_ito_palette,
    rescued_cmo_tag_mapped,
    sc,
):
    rescued_cmo_tag_mapped  # ran after rescued_cmo_tag was mapped onto atac_adata
    atac_spectral_umap_leiden_computed  # ran after spectral embedding, UMAP, KNN, and Leiden clustering

    # sc.pl.umap applies one shared palette dict across all color= panels, so build
    # a single dict covering both namespaces: ec_diff_palette for timepoint, and
    # Okabe-Ito (cycled) for leiden cluster IDs. Use .unique() rather than
    # .cat.categories, since snap.tl.leiden does not always store "leiden" as a
    # categorical dtype (depends on the number of clusters found).
    _leiden_categories = sorted(atac_adata.obs["leiden"].astype(str).unique())
    _leiden_palette = {_cl: okabe_ito_palette[_i % len(okabe_ito_palette)] for _i, _cl in enumerate(_leiden_categories)}
    _atac_umap_palette = {**ec_diff_palette, **_leiden_palette}

    sc.pl.umap(
        atac_adata,
        color=["rescued_cmo_tag", "leiden"],
        palette=_atac_umap_palette,
        ncols=2,
        size=5,
        title=["ATAC UMAP colored by timepoint", "ATAC UMAP colored by Leiden cluster"],
    )
    return


@app.cell(hide_code=True)
def atac_umap_and_confusion_vs_rna_cluster(
    adata_clean,
    adata_clean_clustering_computed,
    alt,
    atac_adata,
    atac_spectral_umap_leiden_computed,
    atac_umap_dataframe,
    atac_umap_scatter,
    mo,
    okabe_ito_palette,
    pd,
):
    atac_spectral_umap_leiden_computed  # ran after ATAC spectral embedding, UMAP, and Leiden clustering
    adata_clean_clustering_computed  # ran after leiden clustering on adata_clean

    # atac_adata's barcodes are a subset of adata_clean's (every ATAC barcode has an
    # RNA counterpart here), so no explicit intersection is needed, but the count
    # is still worth stating explicitly rather than assuming it.
    _shared_barcodes = atac_adata.obs_names.intersection(adata_clean.obs_names)
    _n_shared = len(_shared_barcodes)
    _rna_cluster_on_atac = adata_clean.obs.loc[_shared_barcodes, "leiden"]

    _rna_categories = sorted(_rna_cluster_on_atac.astype(str).unique(), key=int)
    _rna_palette = {_cl: okabe_ito_palette[_i % len(okabe_ito_palette)] for _i, _cl in enumerate(_rna_categories)}

    _umap_df = atac_umap_dataframe(atac_adata[_shared_barcodes], [])
    _umap_df["rna_leiden"] = _rna_cluster_on_atac.astype(str).to_numpy()

    _umap_chart = atac_umap_scatter(
        _umap_df, "rna_leiden", "N", "ATAC UMAP colored by RNA cluster",
        color_scale=alt.Scale(domain=list(_rna_palette), range=list(_rna_palette.values())),
    )

    # Confusion matrix: how ATAC leiden clusters (computed on the spectral
    # embedding) map onto RNA leiden clusters (computed independently on
    # post-filter gene expression). Cluster IDs are arbitrary and not comparable
    # by number across modalities, so this crosstab is the only way to track
    # correspondence.
    _atac_cluster_on_shared = atac_adata.obs.loc[_shared_barcodes, "leiden"].astype(str)
    _confusion = pd.crosstab(_atac_cluster_on_shared.rename("atac_cluster"), _rna_cluster_on_atac.astype(str).rename("rna_cluster"))

    _atac_order = sorted(_confusion.index, key=int)
    _rna_order = _confusion.idxmax(axis=0).sort_values(key=lambda s: s.astype(int)).index.tolist()
    _atac_rank = {v: i for i, v in enumerate(_atac_order)}
    _rna_rank = {v: i for i, v in enumerate(_rna_order)}

    _confusion_long = _confusion.reset_index().melt(id_vars="atac_cluster", var_name="rna_cluster", value_name="n_cells")
    _confusion_long = _confusion_long[_confusion_long["n_cells"] > 0].copy()
    # Altair's sort=<list> shorthand can silently fail under VegaFusion; precomputing
    # the rank as an explicit data column and sorting via EncodingSortField sidesteps it.
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
        color=alt.condition(alt.datum.n_cells > _confusion_long["n_cells"].max() / 2, alt.value("black"), alt.value("white")),
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
        mo.md(f"**{_n_shared:,} barcodes shared between `atac_adata` and `adata_clean`** (of {atac_adata.n_obs:,} ATAC and {adata_clean.n_obs:,} RNA barcodes); both panels below are restricted to this shared set."),
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


@app.cell
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
    adata_clean,
    adata_clean_clustering_computed,
    atac_adata,
    atac_spectral_umap_leiden_computed,
    outdir_root,
    project_root,
):
    adata_clean_clustering_computed  # ran after normalization, HVG, PCA, neighbors, leiden, and UMAP (RNA)
    atac_spectral_umap_leiden_computed  # ran after spectral embedding, UMAP, KNN, and Leiden clustering (ATAC)

    # Saving the RNA and ATAC anndata objects.
    _rna_outfile = outdir_root / "rna_adata.h5ad"
    _atac_outfile = outdir_root / "atac_adata.h5ad"
    _rna_outfile.parent.mkdir(parents=True, exist_ok=True)

    adata_clean.write_h5ad(_rna_outfile)
    atac_adata.write_h5ad(_atac_outfile)

    # Show only the paths relative to the repo, not the full local filesystem path
    [_rna_outfile.relative_to(project_root), _atac_outfile.relative_to(project_root)]
    return


if __name__ == "__main__":
    app.run()
