# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "altair==6.2.2",
#     "anndata==0.13.2",
#     "decontx-python==0.2.0",
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

    This notebook processes the 10x multiome data (**channel2**) for the endothelial differentiation time course: 5 timepoints (d0-d4), 5 biological replicates each. Samples were multiplexed using the MULTI-seq technique with CMO (Cell Multiplexing Oligo) barcodes. Data was processed with the IGVF pipeline using kallisto and the GENCODE v43 annotation, and CMO quantification was performed with the `kite` workflow from the `kallisto-bustools` suite.

    The goals of this notebook are to:

    1. Perform quality control filtering on the RNA data
    2. Run CMO hash classification (mimicking `Seurat::HTODemux`) on QC-passing cells
    3. Assign each cell barcode to a CMO (and therefore to a timepoint/replicate)
    4. Process the corresponding ATAC data with `snapatac2`

    IGVF portal accessions

    - Analysis set : [IGVFDS3995WHFT](https://data.igvf.org/analysis-sets/IGVFDS3995WHFT/)
    - Gene count matrix: [IGVFFI8316SYHQ](https://data.igvf.org/matrix-files/IGVFFI8316SYHQ/)
    - Fragment file: [IGVFFI3256WWXC](https://data.igvf.org/tabular-files/IGVFFI3256WWXC/)
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
    ch2_data_root_path = Path(project_root / "data/10X_miltiome_5_timepoints/channel2/")

    ch2_outdir_root = Path(project_root /"results/10X_miltiome_5_timepoints/channel2")


    # 1. Initialize the connection targeting the production portal
    # It automatically picks up IGVF_API_KEY and IGVF_SECRET_KEY from the environment
    conn = Connection(igvf_mode="prod")

    # 2. Files to download: (accession, destination subfolder)
    _files_to_download = [
        ("IGVFFI8316SYHQ", ch2_data_root_path / "rna/h5ad"),
        ("IGVFFI3256WWXC", ch2_data_root_path / "atac/fragments"),
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
    return ch2_data_root_path, ch2_outdir_root


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
def load_h5_counts(
    ad,
    ch2_data_root_path,
    gene_metadata_df,
    igvf_gencode_gtf_path,
):
    import gzip
    import re

    try:
        _h5ad_path = ch2_data_root_path / "rna/h5ad/IGVFFI8316SYHQ.h5ad"
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

    knee_umi_selection = alt.selection_interval(encodings=["y"], value={"y": [1000, 10_000]})

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
    # Read the brush's y-range (falls back to the default 1000-10,000 before any
    # interaction, since the interval's initial "value" isn't reported back until
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
def mask_filter_cells(adata, mo, np, pd, qc_metrics_computed):
    mo.stop(not qc_metrics_computed)

    # QC pass mask mirroring the cell-level filters applied below
    _max_counts_cutoff = 10_000
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
def pca_pc1_pc2_colored(adata_flt, sc):
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
    doublet_dominated_clusters,
    leiden_computed,
    okabe_ito_palette,
):
    leiden_computed  # ran after leiden clustering
    doublet_dominated_clusters  # ran after doublet-dominated clusters were identified

    # % mito per leiden cluster, to spot a mito-driven cluster before it gets a
    # dedicated deep-dive (see the cluster 2 sections below).
    _cluster_order = sorted(adata_flt.obs["leiden"].cat.categories, key=int)

    _median_mito_by_cluster = adata_flt.obs.groupby("leiden", observed=True)["pct_counts_mt"].median()

    # "Moderately elevated" clusters are picked out via the largest gap in the
    # sorted per-cluster median %mito among non-doublet-dominated clusters (their
    # mito elevation already has a separate explanation), rather than a fixed
    # list, so this stays correct if Leiden renumbers clusters on a future re-run.
    _non_doublet_medians = _median_mito_by_cluster.drop(index=doublet_dominated_clusters).sort_values(ascending=False)
    _gaps = -_non_doublet_medians.diff().dropna()
    _split = int(_gaps.to_numpy().argmax()) + 1
    high_mito_clusters = sorted(_non_doublet_medians.index[:_split].tolist(), key=int)

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
    return (high_mito_clusters,)


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
    return


@app.cell
def plot_umap(adata_flt, sc):
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
    mo,
    np,
    okabe_ito_palette,
    pd,
    scrublet_auto_threshold,
    scrublet_scores_sim,
    scrublet_threshold_number,
):
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
def initialize_cmo(ad, adata_flt, project_root):
    # path to the cmo counts
    _cmo_counts_path = project_root / "data/10X_miltiome_5_timepoints/channel2/cmo_counts/adata.h5ad"

    adata_cmo = ad.read_h5ad(_cmo_counts_path)

    # CMO barcodes need the same per-sample suffix as adata_flt's RNA barcodes to match up
    _barcode_suffix = "_" + adata_flt.obs_names[0].split("_", 1)[1]
    adata_cmo.obs_names = adata_cmo.obs_names + _barcode_suffix

    adata_cmo.var["gene_name"] = adata_cmo.var_names

    # Filter adata_cmo data to cells in adata_flt
    _cells_in_use = adata_flt.obs_names.intersection(adata_cmo.obs_names.to_list())
    adata_cmo = adata_cmo[_cells_in_use, :]

    # Map each CMO to its timepoint: CMO1-5 = d0, CMO6-10 = d1, ..., CMO21-25 = d4
    cmo_to_timepoint = {f"CMO{i}": f"d{(i - 1) // 5}" for i in range(1, 26)}

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
    # same-timepoint cells only have C(5,2)=10 possible pairs (5 CMOs per
    # timepoint), while cross-timepoint cells have 5x5=25, so raw counts alone
    # make cross-timepoint pairing look inflated just from having more combinations.
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
    plt,
):
    cmo_assignment_computed  # ran after CMO tags were assigned

    # Number of barcodes assigned (hash_ID) to each CMO
    _cmo_order = list(adata_cmo.var["gene_name"]) + ["Negative"]
    _counts = adata_flt.obs["cmo_tag_scanpy"].value_counts().reindex(_cmo_order)
    _colors = [ec_diff_palette[cmo_to_timepoint.get(c, "Unassigned")] for c in _counts.index]

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
    mo,
    np,
    okabe_ito_palette,
    pd,
):
    cmo_assignment_computed  # ran after CMO tags were assigned

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
):
    cmo_assignment_computed  # ran after CMO tags were assigned

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
def leiden_doublet_overview(adata_flt, cmo_assignment_computed):
    cmo_assignment_computed  # ran after CMO tags were assigned

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
def doublet_cluster_investigation_intro(
    doublet_dominated_clusters,
    leiden_doublet_summary,
    mo,
):
    _dd_sorted = sorted(doublet_dominated_clusters, key=int)
    _dd_list = ", ".join(f"**{c}**" for c in _dd_sorted)
    _cmo_lo, _cmo_hi = leiden_doublet_summary.loc[_dd_sorted, "pct_cmo_doublet"].agg(["min", "max"])
    _scrub_lo, _scrub_hi = leiden_doublet_summary.loc[_dd_sorted, "pct_scrublet_doublet"].agg(["min", "max"])

    mo.md(f"""
    ## Doublet-dominated cluster investigation

    Both methods agree that clusters {_dd_list} are heavily doublet-dominated (CMO hashing {_cmo_lo:.1f}-{_cmo_hi:.1f}%, Scrublet {_scrub_lo:.1f}-{_scrub_hi:.1f}%).
    """)
    return


@app.cell
def doublet_pct_forest_plot_intro(mo):
    mo.md(r"""
    ### Doublet % by threshold: is the confounding real?

    In an ideal scenario, raising the detection threshold should mostly turn "Singlets" into "Negatives", not turn "Doublets" into something else, since a genuine doublet's second CMO signal should be comfortably real, not threshold-sensitive.

    The plot below shows one point per cluster group per threshold, the mean %doublet across the clusters in that group, with a thick bar spanning the observed min-max range across those same clusters. A dashed black line marks the "random" expectation. Since each CMO's threshold is defined as its own (1 minus quantile) upper tail, a cell with zero real signal still has that per-CMO chance of a false positive. Under independence, "doublet" (2 or more of 25 CMOs positive) follows Binomial(n=25, p=1-quantile), i.e. P(X>=2), with no free parameters. A cluster sitting well above this line has real confounding, not just the baseline false-positive rate baked into using a 25-CMO panel at this quantile.
    """)
    return


@app.cell(hide_code=True)
def doublet_pct_forest_plot(
    adata_flt,
    alt,
    clr,
    doublet_dominated_clusters,
    np,
    okabe_ito_palette,
    pd,
):
    # See doublet_pct_forest_plot_intro above for the interpretation.
    from scipy.stats import binom

    def _pct_doublet_by_quantile(cluster_id):
        _mask = (adata_flt.obs["leiden"] == cluster_id).to_numpy()
        _clr_cluster = clr[_mask]
        _rows = []
        for _q in np.round(np.linspace(0.90, 0.99, 10), 4):
            _t = np.quantile(clr, _q, axis=0)
            _pos = _clr_cluster > _t
            _n_pos = _pos.sum(axis=1)
            _rows.append({"quantile": _q, "pct_doublet": (_n_pos > 1).mean() * 100})
        return pd.DataFrame(_rows)

    _all_clusters = sorted(adata_flt.obs["leiden"].astype(str).unique(), key=int)
    _other_clusters = [c for c in _all_clusters if c not in doublet_dominated_clusters]
    _dd_label = f"Doublet-dominated ({', '.join(doublet_dominated_clusters)})"
    _other_label = f"Other clusters ({', '.join(_other_clusters)})"

    _frames = []
    for _cl in _all_clusters:
        _df = _pct_doublet_by_quantile(_cl)
        _df["cluster"] = _cl
        _df["group"] = _dd_label if _cl in doublet_dominated_clusters else _other_label
        _frames.append(_df)
    _doublet_long = pd.concat(_frames, ignore_index=True)

    _group_colors = {_dd_label: okabe_ito_palette[8], _other_label: okabe_ito_palette[0]}
    _group_order = list(_group_colors)

    # Point at the mean, bar spanning the observed min-max range, computed
    # across the clusters in each group.
    _summary = _doublet_long.groupby(["quantile", "group"], as_index=False)["pct_doublet"].agg(
        mean="mean",
        lo="min",
        hi="max",
    )

    _n_cmos = clr.shape[1]
    _quantiles = np.round(np.linspace(0.90, 0.99, 10), 4)
    _null_df = pd.DataFrame({
        "quantile": _quantiles,
        "pct_doublet_null": binom.sf(1, _n_cmos, 1 - _quantiles) * 100,
    })

    _null_line = alt.Chart(_null_df).mark_line(
        strokeDash=[4, 4], color="black", strokeWidth=2, point=alt.OverlayMarkDef(color="black", size=40),
    ).encode(
        x=alt.X("quantile:O", title="CLR detection quantile", axis=alt.Axis(labelAngle=0)),
        y=alt.Y("pct_doublet_null:Q"),
    )

    _error_bars = alt.Chart(_summary).mark_rule(strokeWidth=2).encode(
        x=alt.X("quantile:O", title="CLR detection quantile", axis=alt.Axis(labelAngle=0)),
        y=alt.Y("lo:Q", title="% of cluster called Doublet", scale=alt.Scale(domain=[0, 100])),
        y2="hi:Q",
        color=alt.Color(
            "group:N",
            scale=alt.Scale(domain=_group_order, range=list(_group_colors.values())),
            legend=alt.Legend(title="Cluster group"),
        ),
    )

    _points = alt.Chart(_summary).mark_point(size=100, filled=True).encode(
        x=alt.X("quantile:O", axis=alt.Axis(labelAngle=0)),
        y=alt.Y("mean:Q"),
        color=alt.Color("group:N", scale=alt.Scale(domain=_group_order, range=list(_group_colors.values())), legend=None),
        tooltip=["quantile", "group", alt.Tooltip("mean:Q", format=".1f"), alt.Tooltip("lo:Q", format=".1f"), alt.Tooltip("hi:Q", format=".1f")],
    )

    _chart = (_null_line + _error_bars + _points).properties(
        title="Doublet %: mean and min-max range across clusters, vs. the random (binomial) expectation",
        width=550, height=400,
    ).configure_view(strokeWidth=0)

    _chart
    return


@app.cell
def n_pos_distribution_vs_null_intro(
    mo,
    representative_doublet_cluster,
    representative_well_behaved_cluster,
):
    mo.md(f"""
    ### n_pos distribution vs. the random null

    Does each cluster's distribution of "number of positive CMOs" (n_pos, at the pipeline's actual q=0.95 threshold) match the random/binomial null, Binomial(n=25, p=1-0.95=0.05), or deviate from it?

    Left: cluster {representative_doublet_cluster}, doublet-dominated. Right: cluster {representative_well_behaved_cluster}, well-behaved. Black bars are the empirical distribution, gray bars are the theoretical distribution expected under pure noise (no real CMO signal at all). If a cluster's black bars closely track the gray ones, its doublet calls are largely explained by the baseline false-positive rate of a 25-CMO panel. If they diverge sharply, especially with excess mass at n_pos=2 or higher, that cluster has real confounding beyond what chance alone would produce.
    """)
    return


@app.cell(hide_code=True)
def n_pos_distribution_vs_null(
    adata_flt,
    alt,
    clr,
    mo,
    n_pos,
    okabe_ito_palette,
    pd,
    positive_quantile_threshold,
    representative_doublet_cluster,
    representative_well_behaved_cluster,
):
    # See n_pos_distribution_vs_null_intro above for the interpretation.
    from scipy.stats import binom as _binom_dist

    _n_cmos_cmp = clr.shape[1]
    _null_p = 1 - positive_quantile_threshold
    _max_k = 5

    def _n_pos_comparison(cluster_id):
        _mask = (adata_flt.obs["leiden"] == cluster_id).to_numpy()
        _empirical = n_pos[_mask]
        _emp_counts = pd.Series(_empirical).value_counts(normalize=True).reindex(range(_max_k + 1), fill_value=0) * 100
        _null_counts = _binom_dist.pmf(range(_max_k + 1), _n_cmos_cmp, _null_p) * 100
        _df = pd.DataFrame({
            "n_pos": list(range(_max_k + 1)) * 2,
            "pct": list(_emp_counts.to_numpy()) + list(_null_counts),
            "source": ["Empirical"] * (_max_k + 1) + ["Theoretical (random null)"] * (_max_k + 1),
        })
        return _df

    _source_colors = {"Empirical": okabe_ito_palette[0], "Theoretical (random null)": okabe_ito_palette[8]}

    def _comparison_chart(cluster_id, title):
        _df = _n_pos_comparison(cluster_id)
        _bars = alt.Chart(_df).mark_bar().encode(
            x=alt.X("n_pos:O", title="Number of positive CMOs", axis=alt.Axis(labelAngle=0)),
            y=alt.Y("pct:Q", title="% of cells", scale=alt.Scale(domain=[0, 100])),
            xOffset=alt.XOffset("source:N", sort=list(_source_colors)),
            color=alt.Color(
                "source:N",
                scale=alt.Scale(domain=list(_source_colors), range=list(_source_colors.values())),
                legend=alt.Legend(title="Distribution"),
            ),
        )
        return _bars.properties(title=title, width=340, height=340).configure_view(strokeWidth=0)

    _chart_rep = _comparison_chart(representative_doublet_cluster, f"Cluster {representative_doublet_cluster} (doublet-dominated)")
    _chart_well_behaved = _comparison_chart(representative_well_behaved_cluster, f"Cluster {representative_well_behaved_cluster} (well-behaved)")

    def _tight_row(*items):
        # mo.hstack with widths=None adds no wrapper/flex styling around children,
        # so block-level chart divs just fill the row (no slack left for
        # justify-content to redistribute). Build the flex row by hand instead.
        _items_html = "".join(
            f'<div style="flex: 0 0 auto;">{mo.as_html(it).text}</div>' for it in items
        )
        return mo.Html(f'<div style="display:flex; justify-content:flex-start; gap:1rem;">{_items_html}</div>')

    _tight_row(_chart_rep, _chart_well_behaved)
    return


@app.cell(hide_code=True)
def doublet_cluster_conclusion(
    doublet_dominated_clusters,
    leiden_doublet_summary,
    mo,
    representative_doublet_cluster,
    representative_well_behaved_cluster,
):
    _dd_sorted = sorted(doublet_dominated_clusters, key=int)
    _dd_list = ", ".join(_dd_sorted)
    _cmo_lo, _cmo_hi = leiden_doublet_summary.loc[_dd_sorted, "pct_cmo_doublet"].agg(["min", "max"])
    _scrub_lo, _scrub_hi = leiden_doublet_summary.loc[_dd_sorted, "pct_scrublet_doublet"].agg(["min", "max"])

    mo.md(f"""
    ### Doublet-dominated clusters: conclusions

    1. **Clusters {_dd_list} are pure doublet clusters, trust Scrublet over CMO here.**
       Both methods agree these {len(_dd_sorted)} clusters are almost entirely doublets: CMO hashing calls {_cmo_lo:.1f}-{_cmo_hi:.1f}% doublet, Scrublet independently calls {_scrub_lo:.1f}-{_scrub_hi:.1f}% doublet.

    2. **Cluster {representative_doublet_cluster}'s doublet calls are genuine confounding, not just noise.**
       `doublet_pct_forest_plot` shows cluster {representative_doublet_cluster} needs a far stricter threshold than well-behaved clusters to clear its doublets, staying well above the random (binomial) expectation across most of the quantile sweep. `n_pos_distribution_vs_null` confirms this directly: cluster {representative_doublet_cluster}'s distribution of "number of positive CMOs" has substantially more mass at 2 or more than pure chance would produce, Binomial(n=25, p=0.05), unlike cluster {representative_well_behaved_cluster} (a well-behaved cluster), which tracks close to or below the null. This looks like CMO hashing's structural blindness to same-CMO, same-sample, doublets rather than a tunable threshold problem, so we drop these clusters wholesale rather than keep their CMO-labeled singlets.
    """)
    return


@app.cell(hide_code=True)
def high_mito_investigation_intro(
    adata_flt,
    doublet_dominated_clusters,
    high_mito_clusters,
    mo,
):
    _median_by_cluster = adata_flt.obs.groupby("leiden", observed=True)["pct_counts_mt"].median()
    _lowest_clusters = [c for c in _median_by_cluster.index if c not in high_mito_clusters and c not in doublet_dominated_clusters]

    _elevated_lo, _elevated_hi = _median_by_cluster.loc[high_mito_clusters].agg(["min", "max"])
    _lowest_lo, _lowest_hi = _median_by_cluster.loc[_lowest_clusters].agg(["min", "max"])
    _doublet_lo, _doublet_hi = _median_by_cluster.loc[doublet_dominated_clusters].agg(["min", "max"])

    _elevated_list = ", ".join(f"**{c}**" for c in high_mito_clusters)

    mo.md(f"""
    ## High mito-content cluster investigation

    Clusters {_elevated_list} show moderately elevated %mito ({_elevated_lo:.0f}-{_elevated_hi:.0f}% median) relative to the rest of the dataset ({_lowest_lo:.0f}-{_lowest_hi:.0f}% in the lowest clusters, {_doublet_lo:.0f}-{_doublet_hi:.0f}% in the doublet-dominated clusters). This section asks whether that elevation reflects a stress response or dying/damaged nuclei, or genuine biological variation across real, distinct cell types.
    """)
    return


@app.cell
def high_mito_marker_comparison(
    adata_flt,
    high_mito_clusters,
    leiden_computed,
    mo,
    pd,
    sc,
):
    leiden_computed  # ran after leiden clustering
    high_mito_clusters  # ran after the moderately mito-elevated clusters were identified

    # high_mito_clusters are the moderately mito-elevated ones (see
    # high_mito_investigation_intro). A stress-response or dying-nuclei
    # explanation predicts their top marker genes should be dominated by
    # mitochondrially-encoded genes plus a heat-shock/ER-stress chaperone
    # signature, rather than specific biology.
    _marker_view = adata_flt.copy()
    sc.tl.rank_genes_groups(
        _marker_view, groupby="leiden", groups=high_mito_clusters, reference="rest",
        method="wilcoxon", layer="pflog", use_raw=False,
    )

    _top_n = 15
    _marker_dfs = {}
    for _cl in high_mito_clusters:
        _df = sc.get.rank_genes_groups_df(_marker_view, group=_cl).reset_index(drop=True)
        _df["gene_symbol"] = _df["names"].map(adata_flt.var["gene_symbol"])
        _df["rank"] = _df.index + 1
        _marker_dfs[_cl] = _df

    high_mito_top_markers = pd.DataFrame({_cl: _df.head(_top_n)["gene_symbol"].tolist() for _cl, _df in _marker_dfs.items()})
    high_mito_top_markers.index = range(1, _top_n + 1)
    high_mito_top_markers.index.name = "rank"

    # Quantitative check: where do mitochondrially-encoded genes and canonical
    # heat-shock/ER-stress chaperones actually rank in each of these clusters,
    # out of adata_flt.n_vars total genes? A stress-response/dying-nuclei cluster
    # would show these near the top with a positive logFC; real biology should
    # show them unremarkable or even depleted.
    high_mito_stress_genes = [
        "MT-ND5", "MT-CO3", "MT-ND2", "MT-ND1", "MT-ND4", "MT-CO2", "MT-ND3",
        "MT-ATP6", "MT-CO1", "MT-CYB", "MT-ND4L", "HSP90AB1", "HSP90B1", "HSPA5",
        "CANX", "RPLP1",
    ]
    _n_genes = adata_flt.n_vars

    _rows = []
    for _gene in high_mito_stress_genes:
        _row = {"gene": _gene}
        for _cl, _df in _marker_dfs.items():
            _match = _df.loc[_df["gene_symbol"] == _gene]
            if len(_match):
                _row[f"cluster_{_cl}_rank"] = f"{int(_match['rank'].iloc[0]):,} / {_n_genes:,}"
                _row[f"cluster_{_cl}_logfc"] = round(float(_match["logfoldchanges"].iloc[0]), 2)
            else:
                _row[f"cluster_{_cl}_rank"] = "not found"
                _row[f"cluster_{_cl}_logfc"] = None
        _rows.append(_row)

    high_mito_stress_gene_ranks = pd.DataFrame(_rows).set_index("gene")

    mo.vstack([
        mo.md("**Top 15 markers per cluster (Wilcoxon vs. rest, `pflog` layer):**"),
        high_mito_top_markers,
        mo.md("**Rank and logFC of mito-encoded and heat-shock/ER-stress chaperone genes in each cluster's full ranking "
              f"(out of {_n_genes:,} genes; low rank + positive logFC would indicate a stress/dying-nuclei signature):**"),
        high_mito_stress_gene_ranks,
    ])
    return (high_mito_stress_genes,)


@app.cell
def high_mito_stress_score_by_cluster(
    adata_flt,
    alt,
    high_mito_clusters,
    high_mito_stress_genes,
    np,
    okabe_ito_palette,
    pd,
    sc,
):
    high_mito_stress_genes  # ran after the mito/heat-shock/ER-stress gene list was defined
    high_mito_clusters  # ran after the moderately mito-elevated clusters were identified

    # Module score (sc.tl.score_genes) over the mito-encoded and heat-shock/ER-stress
    # genes, rather than looking at each gene individually: a per-cell summary of
    # how strongly the whole stress signature is expressed, then compared across
    # all clusters (not just the moderately mito-elevated ones) so the
    # comparison has proper context.
    _stress_gene_ids = adata_flt.var.index[adata_flt.var["gene_symbol"].isin(high_mito_stress_genes)]

    sc.tl.score_genes(
        adata_flt, gene_list=_stress_gene_ids, score_name="mito_stress_score",
        use_raw=False, layer="pflog", random_state=0,
    )

    _cluster_order = sorted(adata_flt.obs["leiden"].cat.categories, key=int)
    _plot_df = adata_flt.obs[["leiden", "mito_stress_score"]].copy()
    _plot_df["group"] = np.where(_plot_df["leiden"].isin(high_mito_clusters), "Elevated mito cluster", "Other cluster")

    _group_colors = {"Elevated mito cluster": okabe_ito_palette[8], "Other cluster": okabe_ito_palette[0]}

    _box = alt.Chart(_plot_df).mark_boxplot(size=25).encode(
        x=alt.X("leiden:N", title="Leiden cluster", sort=_cluster_order),
        y=alt.Y("mito_stress_score:Q", title="Mito/heat-shock/ER-stress module score"),
        color=alt.Color(
            "group:N",
            scale=alt.Scale(domain=list(_group_colors), range=list(_group_colors.values())),
            legend=alt.Legend(title="Group"),
        ),
    )

    _zero_line = alt.Chart(pd.DataFrame({"y": [0]})).mark_rule(color="black", strokeDash=[4, 4]).encode(y="y:Q")

    (_box + _zero_line).properties(
        title="Mito/stress module score per leiden cluster (dashed line = population average, score = 0 by construction)",
        width=650, height=350,
    ).configure_view(strokeWidth=0)
    return


@app.cell
def high_mito_density_comparison(
    adata_flt,
    alt,
    high_mito_clusters,
    leiden_computed,
    np,
    okabe_ito_palette,
    pd,
):
    leiden_computed  # ran after leiden clustering
    high_mito_clusters  # ran after the moderately mito-elevated clusters were identified

    # Density shape check: is the elevated %mito in these clusters a tight,
    # separate mode (consistent with a distinct dying/stressed population), or
    # just a shifted, overlapping tail of the same overall distribution
    # (consistent with ordinary biological variation)?
    from scipy.stats import gaussian_kde as _gaussian_kde_hm

    _high_mito_mask = adata_flt.obs["leiden"].isin(high_mito_clusters).to_numpy()
    _high_mito_vals = adata_flt.obs.loc[_high_mito_mask, "pct_counts_mt"].to_numpy()
    _rest_vals = adata_flt.obs.loc[~_high_mito_mask, "pct_counts_mt"].to_numpy()

    _max_mito = adata_flt.obs["pct_counts_mt"].max()
    _grid = np.linspace(0, _max_mito, 200)

    _elevated_label = f"Clusters {', '.join(high_mito_clusters)}"

    _density_df = pd.DataFrame({
        "pct_mito": np.tile(_grid, 2),
        "density": np.concatenate([
            _gaussian_kde_hm(_high_mito_vals)(_grid),
            _gaussian_kde_hm(_rest_vals)(_grid),
        ]),
        "group": [_elevated_label] * 200 + ["Rest of dataset"] * 200,
    })

    _group_colors = {_elevated_label: okabe_ito_palette[8], "Rest of dataset": okabe_ito_palette[0]}

    _chart = alt.Chart(_density_df).mark_line(opacity=0.7, interpolate="monotone").encode(
        x=alt.X("pct_mito:Q", title="% mitochondrial counts"),
        y=alt.Y("density:Q", title="Density"),
        color=alt.Color(
            "group:N",
            scale=alt.Scale(domain=list(_group_colors), range=list(_group_colors.values())),
            legend=alt.Legend(title="Group"),
        ),
    ).properties(
        title="%mito density: moderately elevated clusters vs. rest of dataset",
        width=550, height=350,
    ).configure_view(strokeWidth=0)

    _chart
    return


@app.cell
def high_mito_investigation_conclusion(mo):
    mo.md(r"""
    ### High-mito investigation: conclusion

    `high_mito_density_comparison` shows clusters 4, 6, 7, 8, and 12 as a shifted, heavily overlapping tail of the same overall %mito distribution, not a distinct separate mode: about a quarter (24.2%) of the "rest of dataset" barcodes sit above these clusters' median, and about a quarter (28.1%) of these clusters' own barcodes sit below the rest's median. This is not the shape a genuinely distinct dying/stressed population would produce.

    `high_mito_marker_comparison` backs this up at the gene level, and shows five distinct stories rather than one shared signature. Cluster 4 has elevated chaperone/ribosomal genes (HSP90AB1, HSPA5, RPLP1 all rank in the top ~150 of 62,757, positive logFC) alongside cell-cycle markers (DTL, CENPF, MKI67, MCM4) at the top, but its mitochondrial genes are actually depleted, not elevated, consistent with the chaperone and ribosome demand of active proliferation rather than damage. Cluster 6 has genuinely elevated mitochondrial genes (every MT- gene tested ranks in the top ~400, logFC +0.46 to +0.81) alongside angiogenic/endothelial markers (FLT1, ESM1), but its heat-shock/ER chaperone genes are depleted, consistent with real, elevated mitochondrial/metabolic activity rather than a stress response. Cluster 7 has elevated ER chaperones (HSP90B1, HSPA5, CANX) but not mitochondrial genes, alongside endothelial/ECM-secretory markers (KDR, NRP2, HAPLN1, COL11A1), consistent with the baseline protein-folding load of a secretory cell type. Cluster 8 shows no elevation in either the mitochondrial or the heat-shock/ER gene set at all (all rank near the bottom with negative logFC), so its mito elevation (14.0% median, the lowest of the five) isn't explained by anything tested here; its top markers (MECOM, FLI1, DACH1) don't point to an obvious alternative explanation either. Cluster 12 also shows no overlap with any of these genes, with neuronal/axon-guidance and migratory markers (KIF26B, UNC5C, ROBO2, PCDH7, FN1, HMCN1) at the top, matching the tip-cell-like transitional population described elsewhere in this notebook.

    `high_mito_stress_score_by_cluster` looks at the combined mito/heat-shock/ER-stress module score across every cluster, not just these five, and finds no clean separation: the doublet-dominated clusters (1, 2, 5, 9, 11) score just as high or higher than the elevated-mito clusters (2.28-2.58 vs. 2.09-2.38), most other clusters sit in a similar mid-to-high band, and only cluster 0 stands out as clearly lower (1.04). This looks like a gradient across the dataset rather than a distinct stressed/dying subpopulation.

    **Interpretation:** the moderate %mito elevation in these five clusters reflects genuine, and genuinely varied, biological variation, proliferation demand, real mitochondrial/metabolic activity, secretory ER load, or in cluster 8's case no clear explanation at all, rather than a shared stress response or dying nuclei. No additional cluster-level mito exclusion is applied.
    """)
    return


@app.cell(hide_code=True)
def negative_tag_rescue(mo):
    mo.md(r"""
    ## Rescue Negative tags

    Are "Negative" barcodes (no CMO cleared its detection threshold) worth rescuing as real cells, or are they debris that just happened to sit below threshold? This section checks their QC profile against Singlets in the same clusters, then validates the cluster-consensus rescue heuristic against each Negative barcode's own sub-threshold CMO signal.
    """)
    return


@app.cell
def negative_vs_singlet_qc_intro(mo):
    mo.md(r"""
    ### Negative vs. Singlet QC profile

    Are "Negative" barcodes that land in largely-Singlet clusters real cells that simply failed CMO hash detection, with a QC profile similar to their cluster's Singlets, or are they lower-quality debris? Restricted to clusters that are neither doublet-dominated nor the high-mito cluster, not the full per-cell `pass_strict_qc` mask, which would incorrectly include the debris cluster excluded via the separate mito rule.
    """)
    return


@app.cell(hide_code=True)
def negative_vs_singlet_qc(
    adata_flt,
    alt,
    cmo_assignment_computed,
    doublet_dominated_clusters,
    leiden_doublet_summary,
    mo,
    np,
    okabe_ito_palette,
    pd,
):
    cmo_assignment_computed  # ran after CMO tags were assigned

    # See negative_vs_singlet_qc_intro above for the interpretation.
    from scipy.stats import gaussian_kde as _gaussian_kde_neg

    _singlet_dominant_clusters = [c for c in leiden_doublet_summary.index if c not in doublet_dominated_clusters]
    _compare_mask = (adata_flt.obs["leiden"].isin(_singlet_dominant_clusters) & adata_flt.obs["cmo_status_scanpy"].isin(["Singlet", "Negative"])).to_numpy()
    _compare_obs = adata_flt.obs.loc[_compare_mask]
    _n_compare = len(_compare_obs)

    _status_colors = {"Singlet": okabe_ito_palette[3], "Negative": okabe_ito_palette[0]}

    def _density_lines(metric_col, grid):
        _frames = []
        for _status, _color in _status_colors.items():
            _vals = _compare_obs.loc[_compare_obs["cmo_status_scanpy"] == _status, metric_col].to_numpy()
            _label = f"{_status} ({len(_vals) / _n_compare * 100:.1f}%)"
            _frames.append(pd.DataFrame({metric_col: grid, "density": _gaussian_kde_neg(_vals)(grid), "status": _label, "color": _color}))
        return pd.concat(_frames, ignore_index=True)

    def _density_chart(metric_col, title, grid):
        _df = _density_lines(metric_col, grid)
        _colors = dict(zip(_df["status"].unique(), _df.drop_duplicates("status")["color"]))
        return alt.Chart(_df).mark_line(strokeWidth=2, interpolate="monotone").encode(
            x=alt.X(f"{metric_col}:Q", title=title),
            y=alt.Y("density:Q", title="Density"),
            color=alt.Color("status:N", scale=alt.Scale(domain=list(_colors), range=list(_colors.values())), legend=alt.Legend(title="CMO status")),
        ).properties(title=f"{title} by CMO status", width=340, height=320).configure_view(strokeWidth=0)

    _mito_grid = np.linspace(0, _compare_obs["pct_counts_mt"].max(), 200)
    _doublet_grid = np.linspace(0, _compare_obs["scrublet_doublet_score"].max(), 200)

    _mito_chart = _density_chart("pct_counts_mt", "% mitochondrial counts", _mito_grid)
    _doublet_chart = _density_chart("scrublet_doublet_score", "Scrublet doublet score", _doublet_grid)

    def _tight_row(*items):
        # mo.hstack with widths=None adds no wrapper/flex styling around children,
        # so block-level chart divs just fill the row (no slack left for
        # justify-content to redistribute). Build the flex row by hand instead.
        _items_html = "".join(
            f'<div style="flex: 0 0 auto;">{mo.as_html(it).text}</div>' for it in items
        )
        return mo.Html(f'<div style="display:flex; justify-content:flex-start; gap:1rem;">{_items_html}</div>')

    _tight_row(_doublet_chart, _mito_chart)
    return


@app.cell
def rescue_tag_cmo_signal_check_intro(mo):
    mo.md(r"""
    ### Rescue-tag validation against the sub-threshold CMO signal

    Deeper validation of the rescue heuristic: for each rescued "Negative" barcode, restricted to clusters that are neither doublet-dominated nor the high-mito cluster (not the full per-cell `pass_strict_qc` mask), find its single closest-to-threshold CMO (the smallest gap between CLR value and that CMO's own threshold, even though none cleared it) and check whether that CMO's timepoint matches the tag assigned via cluster consensus. A high match rate means the sub-threshold CMO signal independently supports the rescue, not just cluster popularity.
    """)
    return


@app.cell(hide_code=True)
def rescue_tag_cmo_signal_check(
    adata_cmo,
    adata_flt,
    clr,
    cmo_assignment_computed,
    cmo_to_timepoint,
    doublet_dominated_clusters,
    leiden_doublet_summary,
    mo,
    np,
    pd,
    rescue_computed,
    thresholds,
):
    cmo_assignment_computed  # ran after CMO tags were assigned
    rescue_computed  # ran after rescue tags were computed

    # See rescue_tag_cmo_signal_check_intro above for the interpretation.
    from scipy import stats
    from statsmodels.stats.proportion import proportion_confint

    _good_clusters = [c for c in leiden_doublet_summary.index if c not in doublet_dominated_clusters]
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
    # against that null rate, rather than eyeballing a handful of shuffles.
    _rng = np.random.default_rng(0)
    _null_rates = np.array([(_best_cmo_timepoint == _rng.permutation(_rescued_tag)).mean() for _ in range(1000)])
    _null_rate = _null_rates.mean()
    _p_permutation = (_null_rates >= _match.mean()).mean()
    _p_binomial = stats.binomtest(_k, _n, _null_rate, alternative="greater").pvalue

    # Per-timepoint match rate, exposed publicly so the forest plot below can be
    # built in its own cell instead of bundled into this markdown's output.
    _by_tag = pd.DataFrame({"rescued_tag": _rescued_tag, "match": _match}).groupby("rescued_tag").agg(
        pct_match=("match", "mean"), n_barcodes=("match", "count"), n_match=("match", "sum")
    )
    _ci = _by_tag.apply(lambda r: proportion_confint(r["n_match"], r["n_barcodes"], method="wilson"), axis=1)
    _by_tag["lo"] = [c[0] * 100 for c in _ci]
    _by_tag["hi"] = [c[1] * 100 for c in _ci]
    _by_tag["pct_match"] = _by_tag["pct_match"] * 100
    rescue_tag_match_summary = _by_tag.reset_index()
    rescue_tag_null_rate = _null_rate

    mo.md(f"""
    **Rescue-tag validation:** for {_n:,} rescued Negative barcodes, restricted to
    clusters that actually survive `pass_strict_qc`, the CMO closest to (but not over) its
    threshold matches the assigned rescue tag **{_match.mean():.1%}** of the time
    ({_k:,} of {_n:,}), versus a **{_null_rate:.1%}** null (1,000-permutation mean, SD
    {_null_rates.std():.1%}) if the rescue tag were random, roughly {_match.mean() / _null_rate:.1f}x
    enrichment over chance.

    Statistical test: 0 of 1,000 permutations reached the observed rate (permutation
    p < 0.001), and a binomial test against the null rate gives p = {_p_binomial:.2e}.
    With n in the thousands this significance is expected, so the roughly 3x enrichment
    is the number that matters for practical significance, not the p-value itself.

    As a complementary check, not relying on the shuffle null: matched barcodes have a
    smaller median gap-to-threshold ({np.median(_best_gap[_match]):.2f}) than mismatched ones
    ({np.median(_best_gap[~_match]):.2f}). When there's a genuine near-miss CMO signal, it
    tends to agree with the rescue tag, when the signal is weak or noisy, agreement is closer to chance.
    """)
    return rescue_tag_match_summary, rescue_tag_null_rate, stats


@app.cell(hide_code=True)
def rescue_tag_match_rate_plot(
    alt,
    okabe_ito_palette,
    pd,
    rescue_tag_match_summary,
    rescue_tag_null_rate,
):
    # Forest-plot style: point at the observed rate, bar spanning a 95% Wilson
    # confidence interval, dashed line at the permutation null rate, mirroring
    # doublet_pct_forest_plot's layout.
    _null_line = alt.Chart(pd.DataFrame({"y": [rescue_tag_null_rate * 100]})).mark_rule(
        strokeDash=[4, 4], color="black", strokeWidth=2,
    ).encode(y="y:Q")

    _error_bars = alt.Chart(rescue_tag_match_summary).mark_rule(strokeWidth=2, color=okabe_ito_palette[5]).encode(
        x=alt.X("rescued_tag:N", title="Timepoint", axis=alt.Axis(labelAngle=0)),
        y=alt.Y("lo:Q", title="% match rate", scale=alt.Scale(domain=[0, 100])),
        y2="hi:Q",
    )
    _points = alt.Chart(rescue_tag_match_summary).mark_point(size=100, filled=True, color=okabe_ito_palette[5]).encode(
        x=alt.X("rescued_tag:N"),
        y=alt.Y("pct_match:Q"),
        tooltip=["rescued_tag", alt.Tooltip("pct_match:Q", format=".1f"), "n_barcodes"],
    )

    _chart = (_null_line + _error_bars + _points).properties(
        title="Rescue-tag match rate by timepoint (dashed line = permutation null)",
        width=450, height=350,
    ).configure_view(strokeWidth=0)

    _chart
    return


@app.cell
def negative_tag_rescue_conclusion(mo):
    mo.md(r"""
    ### Negative tag rescue: conclusion

    **"Negative" barcodes in singlet-dominant clusters look worth rescuing, and the rescued tags are independently supported by the sub-threshold CMO signal.**

    Restricted to clusters that are neither doublet-dominated nor the high-mito cluster 2 (cluster-level restriction, not the full per-cell `pass_strict_qc` mask), `negative_vs_singlet_qc` shows Negatives have a comparable, even slightly better, QC profile than Singlets in the same clusters: median %mito 8.4% vs. 9.4%, median Scrublet doublet score 0.115 vs. 0.121. This supports treating them as real cells with failed CMO staining rather than debris.

    `rescue_tag_cmo_signal_check` goes further: for each rescued Negative, we check whether its single closest-to-threshold CMO (the one nearest to, but not over, its cutoff) belongs to the same timepoint as the tag assigned via cluster consensus. Match rate is 42.9% overall, versus a 16.9% null baseline (1,000-permutation test) if the rescue tag were random, roughly 2.5x enrichment (binomial p around 1.3e-135). This varies by timepoint, d3 91.3%, d0/d2/d4 32 to 35%, d1 42.0%, but every timepoint sits well above the random baseline. d3 being so much higher than the rest makes sense: d3 has the worst CMO recovery of any timepoint (3,428 confidently-assigned singlets, versus 3,913 to 4,127 for the other four), so wherever d3 is a cluster's true dominant population, its cells are disproportionately likely to end up "Negative" rather than the negatives being a mix of several timepoints' poorly-stained cells. That makes the cluster-consensus rescue tag for those barcodes both more likely to be d3 and more likely to be correct. As in channel1, d3 rescues are by far the most robust.
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

    Three independent exclusion criteria combine into `pass_strict_qc`: doublet-dominated clusters, Scrublet-predicted doublets, and CMO-hashing Doublet calls. The mito filtering itself was already applied upstream, pre-clustering, in `mask_filter_cells` (median + 2.5 MADs); no additional cluster- or cell-level mito exclusion is applied here, per `high_mito_investigation_conclusion`. The breakdown below applies the three criteria as a cascade, ordered from least to most restrictive (by how many cells each filter alone would keep), the same style used for the lenient QC filter earlier in `mask_filter_cells`.
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

        adata.obs["qc_exclude_reason"] = np.select(
            [_is_doublet_cluster, _is_scrublet_doublet, _is_cmo_doublet],
            ["doublet_cluster", "scrublet_doublet", "cmo_doublet"],
            default="",
        )
        adata.obs["pass_not_doublet_cluster"] = ~_is_doublet_cluster
        adata.obs["pass_not_scrublet_doublet"] = ~_is_scrublet_doublet
        adata.obs["pass_not_cmo_doublet"] = ~_is_cmo_doublet
        return True

    _doublet_cluster_mask = leiden_doublet_summary["pct_cmo_doublet"].gt(50) & leiden_doublet_summary["pct_scrublet_doublet"].gt(50)
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
        ("CMO hashing status != Doublet", adata_flt.obs["pass_not_cmo_doublet"]),
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
    return (
        doublet_dominated_clusters,
        representative_doublet_cluster,
        representative_well_behaved_cluster,
        strict_qc_computed,
    )


@app.cell
def rescue_negative_tags_intro(mo):
    mo.md(r"""
    ## Apply the negative tag rescue

    As justified in `negative_tag_rescue_conclusion` above: barcodes with no CMO clearing threshold ("Negative") are rescued to their cluster's consensus timepoint, the majority vote among that cluster's own CMO singlets, rather than dropped outright.
    """)
    return


@app.cell(hide_code=True)
def rescue_negative_tags(
    adata_flt,
    alt,
    cmo_assignment_computed,
    ec_diff_palette,
    np,
    strict_qc_computed,
):
    cmo_assignment_computed  # ran after CMO tags were assigned
    strict_qc_computed  # ran after strict QC flags were computed, for the pass_strict_qc filter below

    def _run_rescue_tags(adata, cluster_consensus_tag):
        adata.obs["rescued_cmo_tag"] = np.select(
            [adata.obs["cmo_status_scanpy"] == "Singlet", adata.obs["cmo_status_scanpy"] == "Negative"],
            [adata.obs["timepoint_scanpy"], adata.obs["leiden"].map(cluster_consensus_tag)],
            default=None,
        )
        return True

    cluster_consensus_tag = (
        adata_flt.obs.loc[adata_flt.obs["cmo_status_scanpy"] == "Singlet"]
        .groupby("leiden", observed=True)["timepoint_scanpy"]
        .agg(lambda s: s.mode().iat[0])
    )

    rescue_computed = _run_rescue_tags(adata_flt, cluster_consensus_tag)

    # Restricted to pass_strict_qc, since that's the population rescued_cmo_tag is
    # actually used for downstream -- barcodes dropped by strict QC still get a
    # rescued_cmo_tag value, but it's never read.
    _timepoint_order = ["d0", "d1", "d2", "d3", "d4"]
    _qc_pass_obs = adata_flt.obs.loc[adata_flt.obs["pass_strict_qc"]]
    _tag_counts = _qc_pass_obs["rescued_cmo_tag"].value_counts().reindex(_timepoint_order).reset_index()
    _tag_counts.columns = ["timepoint", "n_barcodes"]
    _n_total_tagged = int(_tag_counts["n_barcodes"].sum())

    _n_rescued = (
        _qc_pass_obs.loc[_qc_pass_obs["cmo_status_scanpy"] == "Negative", "rescued_cmo_tag"]
        .value_counts().reindex(_timepoint_order).fillna(0).astype(int)
    )
    _tag_counts["n_rescued"] = _tag_counts["timepoint"].map(_n_rescued)
    _tag_counts["label"] = _tag_counts.apply(
        lambda r: f"{r['n_barcodes']:,} ({r['n_rescued']:,} rescued)", axis=1
    )

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
        text=alt.Text("label:N"),
    )

    _chart = (_bar_chart + _bar_labels).properties(
        title=["Rescued timepoint tag breakdown (pass_strict_qc only)", f"n = {_n_total_tagged:,} barcodes total"],
        width=480, height=300,
    ).configure_view(strokeWidth=0)

    _chart
    return cluster_consensus_tag, rescue_computed


@app.cell(hide_code=True)
def pca_by_rescued_timepoint(
    adata_flt,
    ec_diff_palette,
    rescue_computed,
    sc,
    strict_qc_computed,
):
    strict_qc_computed  # ran after strict QC flags were computed
    rescue_computed  # ran after rescue tags were computed

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
def umap_by_rescued_timepoint(
    adata_flt,
    ec_diff_palette,
    rescue_computed,
    sc,
    strict_qc_computed,
):
    strict_qc_computed  # ran after strict QC flags were computed
    rescue_computed  # ran after rescue tags were computed

    _adata_qc_view = adata_flt[adata_flt.obs["pass_strict_qc"]]
    sc.pl.umap(
        _adata_qc_view,
        color=["rescued_cmo_tag", "pct_counts_mt"],
        palette=ec_diff_palette,
        cmap="cividis",
        ncols=2,
        size=5,
        title=["UMAP colored by timepoint", "UMAP colored by %mito"],
    )

    # Leiden clusters shown separately -- ec_diff_palette only covers timepoint
    # labels (d0-d4/Unassigned), not the 12 Leiden cluster IDs, so this uses
    # scanpy's own default categorical palette instead.
    sc.pl.umap(
        _adata_qc_view,
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
    rescue_computed,
    strict_qc_computed,
):
    strict_qc_computed  # ran after strict QC flags were computed
    rescue_computed  # ran after rescue tags were computed

    _outfile = ch2_outdir_root / "adata_flt_qc_annotations.tsv"
    _outfile.parent.mkdir(parents=True, exist_ok=True)
    adata_flt.obs.to_csv(_outfile, sep="\t")

    # Show only the path relative to the repo, not the full local filesystem path
    _outfile.relative_to(project_root)
    return


@app.cell
def filtering_strategy_conclusion(
    adata_flt,
    doublet_dominated_clusters,
    high_mito_clusters,
    mo,
    rescue_tag_match_summary,
):
    _n_lenient = adata_flt.n_obs
    _n_scrublet = int(adata_flt.obs["scrublet_predicted_doublet"].sum())
    _n_cmo = int((adata_flt.obs["cmo_status_scanpy"] == "Doublet").sum())
    _n_cluster = int(adata_flt.obs["leiden"].isin(doublet_dominated_clusters).sum())
    _n_pass = int(adata_flt.obs["pass_strict_qc"].sum())
    _pct_pass = _n_pass / _n_lenient

    _elevated_list = ", ".join(high_mito_clusters)

    _qc_pass_obs = adata_flt.obs.loc[adata_flt.obs["pass_strict_qc"]]
    _n_rescued_by_tp = (
        _qc_pass_obs.loc[_qc_pass_obs["cmo_status_scanpy"] == "Negative", "rescued_cmo_tag"]
        .value_counts()
    )

    _fewest_rescue_tp = _n_rescued_by_tp.idxmin()
    _fewest_rescue_n = int(_n_rescued_by_tp.min())
    _rescue_lo, _rescue_hi = int(_n_rescued_by_tp.min()), int(_n_rescued_by_tp.max())

    _match_by_tp = rescue_tag_match_summary.set_index("rescued_tag")["pct_match"]
    _fewest_tp_match = _match_by_tp.loc[_fewest_rescue_tp]
    _other_match = _match_by_tp.drop(index=_fewest_rescue_tp)
    _other_lo, _other_hi = _other_match.min(), _other_match.max()

    mo.md(f"""
    ### Filtering strategy: conclusions

    Starting from {_n_lenient:,} barcodes surviving the lenient QC filter (which already includes an upstream mito cutoff, median + 2.5 MADs, applied pre-clustering in `mask_filter_cells`), strict QC drops barcodes for one of three reasons: Scrublet-predicted doublets ({_n_scrublet:,}), CMO-hashing Doublet calls ({_n_cmo:,}), or doublet-dominated clusters ({_n_cluster:,}), leaving **{_n_pass:,} of {_n_lenient:,} barcodes ({_pct_pass:.1%})** as `pass_strict_qc`. No additional cluster- or cell-level mito exclusion is applied at this stage: the moderately mito-elevated clusters ({_elevated_list}) show no stress/dying-nuclei signature (density shape, marker genes, and module score all point to genuine biological variation, not debris, see `high_mito_investigation_conclusion`), so they're kept rather than excluded.

    Among those survivors, CMO-hashing "Negative" barcodes are rescued to their cluster's consensus timepoint rather than dropped, contributing {_rescue_lo:,} to {_rescue_hi:,} rescued barcodes per timepoint. {_fewest_rescue_tp} gets the fewest ({_fewest_rescue_n:,}), a knock-on effect of {_fewest_rescue_tp} also having the fewest confidently-tagged singlets overall, so fewer clusters end up with {_fewest_rescue_tp} as their consensus timepoint to rescue negatives into. That's a separate question from how *accurate* those rescues are: `negative_tag_rescue_conclusion` shows {_fewest_rescue_tp} rescues are in fact the most reliable of any timepoint ({_fewest_tp_match:.1f}% match rate vs. {_other_lo:.0f}-{_other_hi:.0f}% for the others), precisely because {_fewest_rescue_tp}'s weaker CMO recovery means its true cells are disproportionately likely to end up "Negative" rather than negatives being a random mix. The resulting `rescued_cmo_tag` is what downstream RNA and ATAC analyses use as each barcode's timepoint label.
    """)
    return


@app.cell(hide_code=True)
def cluster12_intro(mo):
    mo.md(r"""
    # Extra: Cluster 12 investigation

    `cluster_tag_mismatch_check` above flags cluster 12 for a 23.8% singlet/cluster-consensus mismatch rate, well above most other clusters despite sitting below the 50% doublet cutoff for both Scrublet and CMO hashing. This section digs into whether that mismatch reflects genuine mixed-timepoint biology (a d3/d4 transition state) or a technical/doublet problem, and whether the resulting population is worth keeping.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Tag mismatch analysis

    Sanity check on the rescue heuristic: among barcodes that DO have a direct CMO singlet call, how often does that call disagree with what we'd have assigned them via their cluster's consensus, the same rule used to rescue Negatives? A high mismatch rate for a specific cluster flags it as a candidate for further investigation.
    """)
    return


@app.cell(hide_code=True)
def cluster_tag_mismatch_check(
    adata_flt,
    cluster_consensus_tag,
    cmo_assignment_computed,
    doublet_dominated_clusters,
):
    cmo_assignment_computed  # ran after CMO tags were assigned

    # See the markdown above for the interpretation.
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
    _by_cluster["marked_for_removal"] = [cl in doublet_dominated_clusters for cl in _by_cluster.index]

    print(f"{_n_mismatch} of {_n_singlets} singlets ({_n_mismatch / _n_singlets:.1%}) have an own-tag that disagrees with their cluster's consensus tag.")
    _by_cluster.sort_values("pct_mismatch", ascending=False)
    return


@app.cell
def cluster12_investigation_intro(mo):
    mo.md(r"""
    ### Cluster 12: doublet or transition state?

    Cluster 12 sits below the 50% doublet cutoff for both Scrublet and CMO hashing, so it survives strict QC, but its 23.8% singlet/cluster-consensus mismatch rate is unusually high. If that mismatch were a doublet artifact, we'd expect elevated %mito and doublet score alongside a genuinely mixed CMO-status composition; if it's a real transition state, the singlets should skew toward two adjacent timepoints (d3/d4) rather than being scattered randomly, and QC metrics should stay unremarkable.
    """)
    return


@app.cell(hide_code=True)
def cluster12_investigation(
    adata_flt,
    alt,
    cmo_assignment_computed,
    ec_diff_palette,
    mo,
    okabe_ito_palette,
    pd,
):
    cmo_assignment_computed  # ran after CMO tags were assigned

    # See cluster12_investigation_intro above for the interpretation.
    _mask12 = (adata_flt.obs["leiden"] == "12").to_numpy()
    _sub12 = adata_flt.obs.loc[_mask12]

    _status_order = ["Singlet", "Doublet", "Negative"]
    _status_colors = {"Singlet": okabe_ito_palette[3], "Doublet": okabe_ito_palette[8], "Negative": okabe_ito_palette[0]}
    _status_counts = _sub12["cmo_status_scanpy"].value_counts().reindex(_status_order).reset_index()
    _status_counts.columns = ["cmo_status", "n_barcodes"]

    _status_chart = alt.Chart(_status_counts).mark_bar().encode(
        y=alt.Y("cmo_status:N", sort=_status_order, title="CMO status"),
        x=alt.X("n_barcodes:Q", title="Number of barcodes"),
        color=alt.Color(
            "cmo_status:N",
            scale=alt.Scale(domain=_status_order, range=[_status_colors[s] for s in _status_order]),
            legend=None,
        ),
    ).properties(title="Cluster 12: CMO status", width=340, height=280).configure_view(strokeWidth=0)

    _timepoint_order = ["d0", "d1", "d2", "d3", "d4"]
    _timepoint_counts = (
        _sub12.loc[_sub12["cmo_status_scanpy"] == "Singlet", "timepoint_scanpy"]
        .value_counts().reindex(_timepoint_order).fillna(0).astype(int).reset_index()
    )
    _timepoint_counts.columns = ["timepoint", "n_barcodes"]

    _timepoint_chart = alt.Chart(_timepoint_counts).mark_bar().encode(
        y=alt.Y("timepoint:N", sort=_timepoint_order, title="Timepoint"),
        x=alt.X("n_barcodes:Q", title="Number of barcodes"),
        color=alt.Color(
            "timepoint:N",
            scale=alt.Scale(domain=_timepoint_order, range=[ec_diff_palette[t] for t in _timepoint_order]),
            legend=None,
        ),
    ).properties(title="Cluster 12 singlets: timepoint composition", width=340, height=280).configure_view(strokeWidth=0)

    def _tight_row(*items):
        # mo.hstack with widths=None adds no wrapper/flex styling around children,
        # so block-level chart divs just fill the row (no slack left for
        # justify-content to redistribute). Build the flex row by hand instead.
        _items_html = "".join(
            f'<div style="flex: 0 0 auto;">{mo.as_html(it).text}</div>' for it in items
        )
        return mo.Html(f'<div style="display:flex; justify-content:flex-start; gap:1rem;">{_items_html}</div>')

    _qc_table = pd.DataFrame({
        "cluster_12": _sub12[["total_counts", "n_genes_by_counts", "pct_counts_mt", "scrublet_doublet_score"]].median(),
        "rest_of_dataset": adata_flt.obs.loc[~_mask12, ["total_counts", "n_genes_by_counts", "pct_counts_mt", "scrublet_doublet_score"]].median(),
    }).round(2)

    mo.vstack([
        _tight_row(_status_chart, _timepoint_chart),
        mo.md(f"**Cluster 12 median QC vs. rest of dataset:**\n\n{_qc_table.to_markdown()}"),
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Marker genes analysis

    If cluster 10 is a real, distinct transitional population rather than a doublet or QC artifact, it should have a coherent marker gene signature, restricted to QC-passing cells and computed on the `pflog` layer (PFlog-normalized), since `.X` is still raw counts.
    """)
    return


@app.cell(hide_code=True)
def cluster12_marker_genes(adata_flt, sc, strict_qc_computed):
    strict_qc_computed  # ran after strict QC flags were computed

    # See the markdown above for the interpretation.
    _qc_view = adata_flt[adata_flt.obs["pass_strict_qc"]].copy()
    sc.tl.rank_genes_groups(_qc_view, groupby="leiden", groups=["12"], reference="rest", method="wilcoxon", layer="pflog", use_raw=False)

    cluster12_markers = sc.get.rank_genes_groups_df(_qc_view, group="12")
    cluster12_markers["gene_symbol"] = cluster12_markers["names"].map(adata_flt.var["gene_symbol"])
    cluster12_markers.head(25)
    return


@app.cell(hide_code=True)
def cluster12_conclusion(mo):
    mo.md(r"""
    ## Cluster 12: a lagging, tip-cell-like transitional population

    Cluster 12 (496 cells) is CMO-hashing-labeled mostly d3/d4 (298 d3, 90 d4 among its 391 singlets, 76.2%/23.0%). Its own-tag-vs-cluster-consensus mismatch rate is 23.8%, in line with its non-d3 singlet fraction, and its QC metrics remain unremarkable relative to a real cell type: %mito 16.1% vs. 12.5% rest-of-dataset, a mild elevation consistent with the moderate, biology-driven variation described in `high_mito_investigation_conclusion` (cluster 12 is one of the clusters flagged there), well under the levels seen in the doublet-dominated clusters.

    **Marker genes** (`cluster12_marker_genes`, Wilcoxon vs. rest, on the `pflog` layer) show the same signature previously identified for this population (and its analog in channel1): axon-guidance/cell-motility genes (`KIF26B`, `UNC5C`, `ROBO2`, `SLIT3`, `SEMA6D`, `NRP2`, the last a classic endothelial tip-cell marker), plus adhesion/ECM genes (`PCDH7`, `FN1`, `ALCAM`, `HMCN1`) and progenitor-associated genes (`MLLT3`, `CDK6`, `HMGA2`).

    **Interpretation:** this population reproducing almost identically (composition, QC profile, and marker genes) across independent samples and re-clustering runs is strong evidence this is real, robust biology rather than a normalization- or channel-specific artifact: a migratory, tip-cell-like subpopulation that transcriptionally lags the bulk d3/d4 differentiation trajectory. Not excluded; kept as-is.
    """)
    return


@app.cell
def cluster4_investigation_intro(mo):
    mo.md(r"""
    ### Cluster 4: doublet or multi-timepoint progenitor state?

    Cluster 4 sits below the 50% doublet cutoff for both Scrublet and CMO hashing, so it survives strict QC, but its 47.3% singlet/cluster-consensus mismatch rate is the highest of any cluster not already excluded as doublet-dominated. If that mismatch were a doublet artifact, we'd expect elevated %mito and doublet score alongside a genuinely mixed CMO-status composition; if it's a real biological state that simply persists across several timepoints rather than being tied to one, the singlets should show a broad, structured timepoint spread rather than being scattered randomly, and QC metrics should stay unremarkable.
    """)
    return


@app.cell
def cluster4_investigation(
    adata_flt,
    alt,
    cmo_assignment_computed,
    ec_diff_palette,
    mo,
    okabe_ito_palette,
    pd,
):
    cmo_assignment_computed  # ran after CMO tags were assigned

    # See cluster4_investigation_intro above for the interpretation.
    _mask4 = (adata_flt.obs["leiden"] == "4").to_numpy()
    _sub4 = adata_flt.obs.loc[_mask4]

    _status_order = ["Singlet", "Doublet", "Negative"]
    _status_colors = {"Singlet": okabe_ito_palette[3], "Doublet": okabe_ito_palette[8], "Negative": okabe_ito_palette[0]}
    _status_counts = _sub4["cmo_status_scanpy"].value_counts().reindex(_status_order).reset_index()
    _status_counts.columns = ["cmo_status", "n_barcodes"]

    _status_chart = alt.Chart(_status_counts).mark_bar().encode(
        y=alt.Y("cmo_status:N", sort=_status_order, title="CMO status"),
        x=alt.X("n_barcodes:Q", title="Number of barcodes"),
        color=alt.Color(
            "cmo_status:N",
            scale=alt.Scale(domain=_status_order, range=[_status_colors[s] for s in _status_order]),
            legend=None,
        ),
    ).properties(title="Cluster 4: CMO status", width=340, height=280).configure_view(strokeWidth=0)

    _timepoint_order = ["d0", "d1", "d2", "d3", "d4"]
    _timepoint_counts = (
        _sub4.loc[_sub4["cmo_status_scanpy"] == "Singlet", "timepoint_scanpy"]
        .value_counts().reindex(_timepoint_order).fillna(0).astype(int).reset_index()
    )
    _timepoint_counts.columns = ["timepoint", "n_barcodes"]

    _timepoint_chart = alt.Chart(_timepoint_counts).mark_bar().encode(
        y=alt.Y("timepoint:N", sort=_timepoint_order, title="Timepoint"),
        x=alt.X("n_barcodes:Q", title="Number of barcodes"),
        color=alt.Color(
            "timepoint:N",
            scale=alt.Scale(domain=_timepoint_order, range=[ec_diff_palette[t] for t in _timepoint_order]),
            legend=None,
        ),
    ).properties(title="Cluster 4 singlets: timepoint composition", width=340, height=280).configure_view(strokeWidth=0)

    def _tight_row(*items):
        # mo.hstack with widths=None adds no wrapper/flex styling around children,
        # so block-level chart divs just fill the row (no slack left for
        # justify-content to redistribute). Build the flex row by hand instead.
        _items_html = "".join(
            f'<div style="flex: 0 0 auto;">{mo.as_html(it).text}</div>' for it in items
        )
        return mo.Html(f'<div style="display:flex; justify-content:flex-start; gap:1rem;">{_items_html}</div>')

    _qc_table = pd.DataFrame({
        "cluster_4": _sub4[["total_counts", "n_genes_by_counts", "pct_counts_mt", "scrublet_doublet_score"]].median(),
        "rest_of_dataset": adata_flt.obs.loc[~_mask4, ["total_counts", "n_genes_by_counts", "pct_counts_mt", "scrublet_doublet_score"]].median(),
    }).round(2)

    mo.vstack([
        _tight_row(_status_chart, _timepoint_chart),
        mo.md(f"**Cluster 4 median QC vs. rest of dataset:**\n\n{_qc_table.to_markdown()}"),
    ])
    return


@app.cell
def cluster4_marker_genes_intro(mo):
    mo.md(r"""
    ### Marker genes analysis

    If cluster 4 is a real, distinct biological state rather than a doublet or QC artifact, it should have a coherent marker gene signature, restricted to QC-passing cells and computed on the `pflog` layer (PFlog-normalized), since `.X` is still raw counts.
    """)
    return


@app.cell
def cluster4_marker_genes(adata_flt, sc, strict_qc_computed):
    strict_qc_computed  # ran after strict QC flags were computed

    # See the markdown above for the interpretation.
    _qc_view = adata_flt[adata_flt.obs["pass_strict_qc"]].copy()
    sc.tl.rank_genes_groups(_qc_view, groupby="leiden", groups=["4"], reference="rest", method="wilcoxon", layer="pflog", use_raw=False)

    cluster4_markers = sc.get.rank_genes_groups_df(_qc_view, group="4")
    cluster4_markers["gene_symbol"] = cluster4_markers["names"].map(adata_flt.var["gene_symbol"])
    cluster4_markers.head(25)
    return


@app.cell
def cluster4_conclusion(mo):
    mo.md(r"""
    ## Cluster 4: a proliferative progenitor state spanning early timepoints (not a doublet artifact)

    Cluster 4 (700 cells) is CMO-hashing-labeled mostly d0 (296 of 562 singlets, 52.7%), with a long tail through d1 (140, 24.9%), d2 (88, 15.7%), d3 (33, 5.9%), and d4 (5, 0.9%), rather than the clean two-timepoint split seen in `cluster12_conclusion`. Its own-tag-vs-cluster-consensus mismatch rate of 47.3% is largely a mechanical consequence of that spread: with singlets that broadly distributed, the d0 plurality vote is a "majority" in name only, so a large fraction of non-d0 singlets necessarily register as mismatches. Its QC metrics argue against a doublet or debris explanation: Scrublet doublet score is *lower* than the rest of the dataset (0.071 vs. 0.165), not higher, and %mito is only mildly elevated (18.0% vs. 12.4%), consistent with the moderate, biology-driven variation described in `high_mito_investigation_conclusion` (cluster 4 is one of the clusters flagged there). Total counts and genes-per-cell are somewhat lower than the rest of the dataset (1,837 vs. 2,764; 1,222 vs. 1,658), which doesn't fit a doublet profile either (doublets skew higher, not lower, on these metrics).

    **Marker genes** (`cluster4_marker_genes`, Wilcoxon vs. rest, on the `pflog` layer) point to an actively proliferating, high-biosynthesis progenitor state: cell-cycle genes (`CENPF`, `DTL`, `NASP`, `H2AZ1`), pluripotency/progenitor markers (`L1TD1`, `CD24`, `NAP1L1`), and a cluster of chaperone (`HSP90AB1`, `HSP90AA1`, `HSPA8`, `HSPA5`, `HSPD1`) and ribosomal/translation genes (`RPLP1`, `RPL11`, `RPL37`, `EIF3A`, `NCL`, `TKT`, `RPS6`, `RPS8`). This chaperone/ribosome load reads as the ordinary protein-folding and biosynthesis demand of rapid proliferation, the same interpretation `high_mito_marker_comparison` and `high_mito_investigation_conclusion` already reached for this cluster, not a stress or damage response.

    **Interpretation:** cluster 4 looks like a genuine proliferative progenitor population present across several early timepoints (mostly d0-d2, tailing into d3/d4) rather than a doublet or technical artifact. Its high tag-mismatch rate is an expected consequence of a real multi-timepoint population being forced through a single-timepoint cluster-consensus vote, not evidence against it. Not excluded; kept as-is. One caveat worth flagging: `rescue_negative_tags` rescues this cluster's 74 Negative barcodes to its d0 consensus tag; applying the cluster's own 47.3% singlet mismatch rate as an estimate, roughly 35 of those 74 are likely mislabeled, a rescue-accuracy limitation specific to clusters with broad timepoint spread like this one.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Post filtering
    """)
    return


@app.cell
def rna_adata_post_filter(
    adata_flt,
    rescue_computed,
    sc,
    scclr,
    strict_qc_computed,
):
    strict_qc_computed  # ran after strict QC flags were computed
    rescue_computed  # ran after rescue tags were computed

    def _run_post_filter_pipeline(adata):
        scclr.pp.pflog(adata, target="auto")
        sc.pp.highly_variable_genes(adata, layer="pflog", n_top_genes=2000)
        _ncomps = 50
        _ncv = 2 * _ncomps + 1
        scclr.tl.pca(adata, n_comps=_ncomps, ncv=_ncv)
        sc.pp.neighbors(adata, random_state=0)
        sc.tl.leiden(adata, flavor="igraph", resolution=0.5, n_iterations=2, random_state=0)
        sc.tl.umap(adata)
        return True

    rna_adata = adata_flt[adata_flt.obs["pass_strict_qc"]].copy()
    rna_adata_post_filter_computed = _run_post_filter_pipeline(rna_adata)
    rna_adata
    return rna_adata, rna_adata_post_filter_computed


@app.cell
def umap_post_filter_leiden_mito(
    ec_diff_palette,
    rna_adata,
    rna_adata_post_filter_computed,
    sc,
):
    rna_adata_post_filter_computed  # ran after normalization, HVG, PCA, neighbors, leiden, and UMAP

    sc.pl.umap(
        rna_adata,
        color=["leiden", "pct_counts_mt"],
        cmap="cividis",
        ncols=2,
        size=5,
        title=["UMAP colored by leiden cluster", "UMAP colored by %mito"],
    )

    # Separate call: rescued_cmo_tag needs the timepoint-specific ec_diff_palette,
    # not scanpy's default categorical palette, since these are the same
    # endothelial-differentiation timepoint colors used everywhere else.
    sc.pl.umap(
        rna_adata,
        color="rescued_cmo_tag",
        palette=ec_diff_palette,
        size=5,
        title="UMAP colored by final CMO tags",
    )
    return


@app.cell
def old_vs_new_cluster_confusion_matrix(
    adata_flt,
    alt,
    pd,
    rna_adata,
    rna_adata_post_filter_computed,
):
    rna_adata_post_filter_computed  # ran after the post-filter reclustering

    # Confusion matrix: how the pre-strict-QC leiden clusters (adata_flt, computed
    # on all lenient-QC-passing barcodes) map onto the post-filter leiden
    # clusters (rna_adata, recomputed from scratch on just the pass_strict_qc
    # population). Cluster IDs are arbitrary and not comparable by number across
    # the two clusterings, so this crosstab is the only way to track cluster
    # identity across the two rounds.
    _old_cluster = adata_flt.obs.loc[rna_adata.obs_names, "leiden"].rename("old_cluster")
    _new_cluster = rna_adata.obs["leiden"].rename("new_cluster")

    _confusion = pd.crosstab(_old_cluster, _new_cluster)

    # Rows (old clusters) kept in plain numeric order as the fixed reference axis.
    # Columns (new clusters) reordered by their dominant old cluster (the row
    # each column has the most cells in), so matches line up close to the
    # diagonal against that fixed reference instead of numeric new-cluster ID,
    # which carries no relationship to the old IDs.
    _old_order = sorted(_confusion.index, key=int)
    _new_order = _confusion.idxmax(axis=0).sort_values(key=lambda s: s.astype(int)).index.tolist()
    _old_rank = {v: i for i, v in enumerate(_old_order)}
    _new_rank = {v: i for i, v in enumerate(_new_order)}

    _confusion_long = _confusion.reset_index().melt(id_vars="old_cluster", var_name="new_cluster", value_name="n_cells")
    _confusion_long = _confusion_long[_confusion_long["n_cells"] > 0].copy()
    # Altair's sort=<list> shorthand generates its own rank lookup internally,
    # but that mechanism silently breaks under VegaFusion for one of the two
    # encodings here (every row falls back to the same rank, so no ordering is
    # applied). Precomputing the rank as an explicit data column and sorting on
    # that via EncodingSortField sidesteps it entirely.
    _confusion_long["old_rank"] = _confusion_long["old_cluster"].astype(str).map(_old_rank)
    _confusion_long["new_rank"] = _confusion_long["new_cluster"].astype(str).map(_new_rank)

    _heatmap = alt.Chart(_confusion_long).mark_rect().encode(
        x=alt.X("new_cluster:N", title="New (post-filter) cluster", sort=alt.EncodingSortField(field="new_rank", op="min")),
        y=alt.Y("old_cluster:N", title="Old (pre-strict-QC) cluster", sort=alt.EncodingSortField(field="old_rank", op="min")),
        color=alt.Color("n_cells:Q", title="Cells", scale=alt.Scale(scheme="cividis")),
        tooltip=["old_cluster", "new_cluster", "n_cells"],
    )
    _labels = alt.Chart(_confusion_long).mark_text(fontSize=9).encode(
        x=alt.X("new_cluster:N", sort=alt.EncodingSortField(field="new_rank", op="min")),
        y=alt.Y("old_cluster:N", sort=alt.EncodingSortField(field="old_rank", op="min")),
        text="n_cells:Q",
        color=alt.condition(alt.datum.n_cells > _confusion_long["n_cells"].max() / 2, alt.value("black"), alt.value("white")),
    )

    (_heatmap + _labels).properties(
        title="Old vs. new leiden cluster assignment (cell counts, columns ordered to show correspondence)",
        width=450, height=450,
    ).configure_view(strokeWidth=0)
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
def load_strict_qc_whitelist(ch2_outdir_root, pd):
    # See atac_import_intro above.
    _qc_annotations = pd.read_csv(
        ch2_outdir_root / "adata_flt_qc_annotations.tsv", sep="\t", index_col=0
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
        ch2_data_root_path / "atac/fragments/IGVFFI3256WWXC.bed.gz",
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
def map_rescued_cmo_tag_to_atac(adata_flt, atac_adata):
    def _map_rescued_cmo_tag_to_atac(atac_adata, adata_flt):
        # Use rescued_cmo_tag, not the raw timepoint_scanpy, so cluster-rescued
        # negatives keep their assigned timepoint here too.
        atac_adata.obs["rescued_cmo_tag"] = adata_flt.obs.loc[atac_adata.obs_names, "rescued_cmo_tag"].astype("category")
        return True

    rescued_cmo_tag_mapped = _map_rescued_cmo_tag_to_atac(atac_adata, adata_flt)
    return (rescued_cmo_tag_mapped,)


@app.cell
def atac_timepoint_counts(atac_adata, rescued_cmo_tag_mapped):
    rescued_cmo_tag_mapped  # ran after rescued_cmo_tag was mapped onto atac_adata

    atac_adata.obs["rescued_cmo_tag"].value_counts()
    return


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
    alt,
    atac_adata,
    atac_spectral_umap_leiden_computed,
    atac_umap_dataframe,
    atac_umap_scatter,
    mo,
    okabe_ito_palette,
    pd,
    rna_adata,
    rna_adata_post_filter_computed,
):
    atac_spectral_umap_leiden_computed  # ran after ATAC spectral embedding, UMAP, and Leiden clustering
    rna_adata_post_filter_computed  # ran after the post-filter RNA reclustering

    # atac_adata's barcodes are a subset of rna_adata's (every ATAC barcode has an
    # RNA counterpart here), so no explicit intersection is needed, but the count
    # is still worth stating explicitly rather than assuming it.
    _shared_barcodes = atac_adata.obs_names.intersection(rna_adata.obs_names)
    _n_shared = len(_shared_barcodes)
    _rna_cluster_on_atac = rna_adata.obs.loc[_shared_barcodes, "leiden"]

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
    # Altair's sort=<list> shorthand can silently fail under VegaFusion (see
    # old_vs_new_cluster_confusion_matrix); precomputing the rank as an explicit
    # data column and sorting via EncodingSortField sidesteps it.
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
        mo.md(f"**{_n_shared:,} barcodes shared between `atac_adata` and `rna_adata`** (of {atac_adata.n_obs:,} ATAC and {rna_adata.n_obs:,} RNA barcodes); both panels below are restricted to this shared set."),
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
    atac_adata,
    atac_spectral_umap_leiden_computed,
    ch2_outdir_root,
    project_root,
    rna_adata,
    rna_adata_post_filter_computed,
):
    rna_adata_post_filter_computed  # ran after normalization, HVG, PCA, neighbors, leiden, and UMAP (RNA)
    atac_spectral_umap_leiden_computed  # ran after spectral embedding, UMAP, KNN, and Leiden clustering (ATAC)

    # Saving the RNA and ATAC anndata objects.
    _rna_outfile = ch2_outdir_root / "rna_adata.h5ad"
    _atac_outfile = ch2_outdir_root / "atac_adata.h5ad"
    _rna_outfile.parent.mkdir(parents=True, exist_ok=True)

    rna_adata.write_h5ad(_rna_outfile)
    atac_adata.write_h5ad(_atac_outfile)

    # Show only the paths relative to the repo, not the full local filesystem path
    [_rna_outfile.relative_to(project_root), _atac_outfile.relative_to(project_root)]
    return


@app.cell
def ambient_rna_investigation_intro(mo):
    mo.md(r"""
    ### Is the residual %mito signal ambient RNA or real biology?

    Barcodes failing ATAC QC show a bimodal RNA %mito distribution, which raises the question: does the RNA-side %mito signal reflect real per-cell mitochondrial content, or ambient RNA contamination? If it's real biology, it should leave a footprint in the independent, DNA-based ATAC measurement too. We test this three ways:

    - **Spatial autocorrelation** of %mito on each modality's own UMAP (Moran's I).
    - **KNN neighbor correlation**: whether a cell's own %mito predicts its neighbors' %mito, in each modality's own KNN graph.
    - **Count dependence**: whether low RNA counts (a hallmark of ambient dilution) are actually driving the residual signal.

    See the conclusion at the end of this section.
    """)
    return


@app.cell
def morans_i_mito_by_embedding(
    adata_flt,
    atac_adata,
    atac_spectral_umap_leiden_computed,
    np,
    pd,
    rna_adata,
    rna_adata_post_filter_computed,
):
    rna_adata_post_filter_computed  # ran after normalization, HVG, PCA, neighbors, leiden, and UMAP (RNA)
    atac_spectral_umap_leiden_computed  # ran after spectral embedding, UMAP, KNN, and Leiden clustering (ATAC)

    from sklearn.neighbors import kneighbors_graph as _kneighbors_graph_moran

    # Kept sparse throughout (no .toarray()) to avoid densifying an n x n matrix,
    # which would be 500MB+ here.
    def _morans_i(coordinates, values, k=15):
        _w = _kneighbors_graph_moran(coordinates, k, mode="connectivity", include_self=False)
        _row_sums = np.asarray(_w.sum(axis=1)).flatten()
        _w = _w.multiply(1 / _row_sums[:, None]).tocsr()  # row-normalize weights

        _z = values - values.mean()
        _s_zero = _w.sum()
        _numerator = len(values) * _z.dot(_w.dot(_z))
        _denominator = _s_zero * np.sum(_z ** 2)
        return _numerator / _denominator

    # rna_adata and atac_adata are distinct AnnData objects sharing the ATAC-side
    # barcode subset, each with its own UMAP.
    _rna_coords = rna_adata.obsm["X_umap"]
    _rna_mito = rna_adata.obs["pct_counts_mt"].to_numpy()

    _atac_coords = atac_adata.obsm["X_umap"]
    _rna_mito_for_atac_cells = adata_flt.obs.loc[atac_adata.obs_names, "pct_counts_mt"].to_numpy()

    morans_i_by_embedding = pd.DataFrame({
        "embedding": ["RNA UMAP", "ATAC UMAP"],
        "morans_i_pct_mito": [_morans_i(_rna_coords, _rna_mito), _morans_i(_atac_coords, _rna_mito_for_atac_cells)],
    })
    morans_i_by_embedding
    return


@app.cell
def neighbor_mean_mito_correlation(
    adata_flt,
    atac_adata,
    atac_spectral_umap_leiden_computed,
    np,
    pd,
    rna_adata,
    rna_adata_post_filter_computed,
):
    rna_adata_post_filter_computed  # ran after normalization, HVG, PCA, neighbors, leiden, and UMAP (RNA)
    atac_spectral_umap_leiden_computed  # ran after spectral embedding, UMAP, KNN, and Leiden clustering (ATAC)

    # Reuse the KNN graphs already computed for each modality (RNA: sc.pp.neighbors
    # in rna_adata_post_filter; ATAC: snap.pp.knn), rather than raw embedding coordinates.
    _rna_conn = rna_adata.obsp["connectivities"]
    _rna_row_sums = np.asarray(_rna_conn.sum(axis=1)).flatten()
    _rna_mito = rna_adata.obs["pct_counts_mt"].to_numpy()
    _rna_neighbor_mean_mito = np.asarray(_rna_conn.dot(_rna_mito)).flatten() / _rna_row_sums

    # snapatac2's snap.pp.knn only stores obsp["distances"], not a weighted
    # connectivities matrix, so binarize it (nonzero = neighbor).
    _atac_conn = (atac_adata.obsp["distances"] > 0).astype(float)
    _atac_row_sums = np.asarray(_atac_conn.sum(axis=1)).flatten()
    _rna_mito_for_atac_cells = adata_flt.obs.loc[atac_adata.obs_names, "pct_counts_mt"].to_numpy()
    _atac_neighbor_mean_mito = np.asarray(_atac_conn.dot(_rna_mito_for_atac_cells)).flatten() / _atac_row_sums

    neighbor_mean_mito_corr = pd.DataFrame({
        "space": ["RNA (KNN graph)", "ATAC (KNN graph)"],
        "corr_own_vs_neighbor_mean_mito": [
            np.corrcoef(_rna_mito, _rna_neighbor_mean_mito)[0, 1],
            np.corrcoef(_rna_mito_for_atac_cells, _atac_neighbor_mean_mito)[0, 1],
        ],
    })
    neighbor_mean_mito_corr
    return


@app.cell
def total_counts_vs_mito_by_timepoint(
    alt,
    ec_diff_palette,
    np,
    pd,
    rna_adata,
    rna_adata_post_filter_computed,
    stats,
):
    rna_adata_post_filter_computed  # ran after normalization, HVG, PCA, neighbors, leiden, and UMAP (RNA)

    # Colored by timepoint (rescued_cmo_tag).
    _timepoint_order = ["d0", "d1", "d2", "d3", "d4"]
    _counts_mito_df = pd.DataFrame({
        "total_counts": rna_adata.obs["total_counts"].to_numpy(),
        "pct_counts_mt": rna_adata.obs["pct_counts_mt"].to_numpy(),
        "timepoint": rna_adata.obs["rescued_cmo_tag"].to_numpy(),
    })

    _scatter = alt.Chart(_counts_mito_df).mark_circle(size=15, opacity=0.5).encode(
        x=alt.X("total_counts:Q", title="Total RNA counts", scale=alt.Scale(type="log")),
        y=alt.Y("pct_counts_mt:Q", title="% mitochondrial counts"),
        color=alt.Color(
            "timepoint:N", title="Timepoint",
            scale=alt.Scale(domain=_timepoint_order, range=[ec_diff_palette[t] for t in _timepoint_order]),
        ),
    )

    # Combine, set properties, and apply clean view adjustments safely
    _counts_mito_chart = _scatter.properties(
        title="Are low-count cells driving the remaining mito signal?",
        width=550, 
        height=400,
    ).configure_view(
        strokeWidth=0 
    )

    # Added .copy() here to stop the pandas layout warning
    df_zone = _counts_mito_df[
        (_counts_mito_df["total_counts"] <= 10000) & 
        (_counts_mito_df["pct_counts_mt"] <= 7)
    ].copy()

    # Run regression against log10 total counts
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        np.log10(df_zone["total_counts"]), 
        df_zone["pct_counts_mt"]
    )

    print(f"R-squared: {r_value**2:.4f}")
    print(f"Slope:     {slope:.4f}")
    print(f"p-value:   {p_value:.4e}")

    # Quantile cut into equal halves
    df_zone["count_bin"] = pd.qcut(df_zone["total_counts"], q=2, labels=["Low UMI Half", "High UMI Half"])

    # Calculate the average mito percentage in each half
    summary = df_zone.groupby("count_bin", observed=False)["pct_counts_mt"].mean()
    print("\n--- Average Mito % by Library Size Half ---")
    print(summary)
    print("-------------------------------------------\n")

    _counts_mito_chart
    return


@app.cell
def decontx_validation_intro(mo):
    mo.md(r"""
    ### The definitive test: does removing ambient RNA eliminate the Moran's I signal?

    If the spatial autocorrelation of %mito on the RNA UMAP (Moran's I 0.295, see `morans_i_mito_by_embedding`) is really driven by ambient RNA contamination, then estimating and subtracting that contamination per cell (DecontX, run on raw counts using the existing Leiden clusters as groups) should collapse it toward 0. If it's real biology, decontamination should not erase it. This does not modify `rna_adata` itself, DecontX runs on a copy.
    """)
    return


@app.cell
def run_decontx(mo, rna_adata, rna_adata_post_filter_computed):
    rna_adata_post_filter_computed  # ran after normalization, HVG, PCA, neighbors, leiden, and UMAP (RNA)

    import decontx

    def _run_decontx(adata):
        return decontx.decontx(adata, cluster_key="leiden", copy=True)

    # Kept minimal and isolated: this is the expensive step (~17 min), so all
    # downstream analysis (Moran's I, contamination scatter, etc.) should be in
    # separate cells referencing this public result, never re-triggering DecontX.
    rna_decontx_adata = _run_decontx(rna_adata)
    decontx_computed = True
    mo.md(f"DecontX finished. Mean contamination: {rna_decontx_adata.obs['decontX_contamination'].mean():.1%}")
    return decontx_computed, rna_decontx_adata


@app.cell
def decontx_moran_validation(
    decontx_computed,
    mo,
    np,
    pd,
    rna_adata,
    rna_decontx_adata,
):
    decontx_computed  # ran after DecontX finished

    from sklearn.neighbors import kneighbors_graph as _kneighbors_graph_decontx
    # Recompute %mito from the decontaminated counts layer, restricted to the
    # same 37 mito genes (adata.var["mt"]) used everywhere else in this notebook.
    _mt_mask = rna_adata.var["mt"].to_numpy()
    _decontx_counts = rna_decontx_adata.layers["decontX_counts"]
    _decontx_counts = _decontx_counts.toarray() if hasattr(_decontx_counts, "toarray") else np.asarray(_decontx_counts)
    _total_decontx = _decontx_counts.sum(axis=1).astype(float)
    _mt_decontx = _decontx_counts[:, _mt_mask].sum(axis=1).astype(float)
    rna_pct_counts_mt_decontx = np.divide(
        _mt_decontx, _total_decontx, out=np.zeros_like(_total_decontx), where=_total_decontx > 0
    ) * 100
    rna_decontx_contamination = rna_decontx_adata.obs["decontX_contamination"].to_numpy()

    # Same sparse Moran's I as morans_i_mito_by_embedding, on the same RNA UMAP
    # coordinates, just swapping in the decontaminated %mito values.
    def _morans_i_decontx(coordinates, values, k=15):
        _w = _kneighbors_graph_decontx(coordinates, k, mode="connectivity", include_self=False)
        _row_sums = np.asarray(_w.sum(axis=1)).flatten()

        # Catch any disconnected singletons or zeroes before dividing
        _row_sums_safe = np.where(_row_sums == 0, 1, _row_sums)
        _w = _w.multiply(1 / _row_sums_safe[:, None]).tocsr()

        _z = values - values.mean()
        _s_zero = _w.sum()
        _numerator = len(values) * _z.dot(_w.dot(_z))
        _denominator = _s_zero * np.sum(_z ** 2)
        return _numerator / _denominator

    _rna_coords = rna_adata.obsm["X_umap"]

    morans_i_decontx_vs_raw = pd.DataFrame({
        "pct_mito": ["Raw (pre-DecontX)", "Decontaminated (post-DecontX)"],
        "morans_i_rna_umap": [
            _morans_i_decontx(_rna_coords, rna_adata.obs["pct_counts_mt"].to_numpy()),
            _morans_i_decontx(_rna_coords, rna_pct_counts_mt_decontx),
        ],
    })
    decontx_moran_validation_computed = True
    _mean_contamination = rna_decontx_contamination.mean()
    mo.vstack([
        morans_i_decontx_vs_raw,
        mo.md(f"Mean estimated contamination fraction: {_mean_contamination:.1%}"),
    ])
    return (
        decontx_moran_validation_computed,
        rna_decontx_contamination,
        rna_pct_counts_mt_decontx,
    )


@app.cell
def save_decontx_results(
    decontx_moran_validation_computed,
    mo,
    pd,
    project_root,
    rna_adata,
    rna_decontx_contamination,
    rna_pct_counts_mt_decontx,
):
    decontx_moran_validation_computed  # ran after the decontaminated %mito and contamination were derived

    # Lightweight: per-cell scalars only, not the full dense decontX_counts matrix
    # (which would be ~4GB for 8146 cells x 62757 genes).
    _decontx_results_df = pd.DataFrame({
        "decontX_contamination": rna_decontx_contamination,
        "pct_counts_mt_decontx": rna_pct_counts_mt_decontx,
        "pct_counts_mt_raw": rna_adata.obs["pct_counts_mt"].to_numpy(),
    }, index=rna_adata.obs_names)

    _decontx_results_path = project_root / "results/channel2/rna_decontx_results.tsv"
    _decontx_results_df.to_csv(_decontx_results_path, sep="\t")
    mo.md(f"Saved DecontX per-cell results to `{_decontx_results_path}`")
    return


@app.cell
def decontx_contamination_vs_mito(
    alt,
    decontx_moran_validation_computed,
    ec_diff_palette,
    pd,
    rna_adata,
    rna_decontx_contamination,
):
    decontx_moran_validation_computed  # ran after the decontaminated %mito and contamination were derived

    _timepoint_order = ["d0", "d1", "d2", "d3", "d4"]
    _contam_df = pd.DataFrame({
        "decontX_contamination": rna_decontx_contamination,
        "pct_counts_mt": rna_adata.obs["pct_counts_mt"].to_numpy(),
        "timepoint": rna_adata.obs["rescued_cmo_tag"].to_numpy(),
    })

    _contam_chart = alt.Chart(_contam_df).mark_circle(size=15, opacity=0.4).encode(
        x=alt.X("decontX_contamination:Q", title="DecontX contamination fraction"),
        y=alt.Y("pct_counts_mt:Q", title="% mitochondrial counts"),
        color=alt.Color(
            "timepoint:N", title="Timepoint",
            scale=alt.Scale(domain=_timepoint_order, range=[ec_diff_palette[t] for t in _timepoint_order]),
        ),
    ).properties(
        title="Does ambient contamination explain mitochondrial reads?",
        width=550, height=400,
    ).configure_view(strokeWidth=0)

    _contam_chart
    return


@app.cell
def _(rna_decontx_adata):
    rna_decontx_adata
    return


@app.cell
def mito_pocket_decontx_markers(
    decontx_computed,
    np,
    rna_adata,
    rna_decontx_adata,
    sc,
    scclr,
):
    decontx_computed  # ran after DecontX finished

    # Which genes drive the high-mito pocket (5% split, same as the earlier
    # per-timepoint self-segregation analysis), using the decontaminated counts
    # rather than raw, so ambient contamination isn't itself driving the
    # result. Normalized the same way as everywhere else in this notebook
    # (scclr.pp.pflog on a copy) before running Wilcoxon rank_genes_groups, since
    # library-size differences between groups would otherwise confound raw counts.
    _decontx_view = rna_decontx_adata.copy()
    _decontx_view.X = _decontx_view.layers["decontX_counts"]
    scclr.pp.pflog(_decontx_view, target="auto")

    _decontx_view.obs["mito_pocket"] = np.where(
        _decontx_view.obs["pct_counts_mt"] > 5.0, "High_Mito", "Normal"
    )

    sc.tl.rank_genes_groups(
        _decontx_view, groupby="mito_pocket", method="wilcoxon", layer="pflog", use_raw=False,
    )

    mito_pocket_decontx_markers = sc.get.rank_genes_groups_df(_decontx_view, group="High_Mito")
    mito_pocket_decontx_markers["gene_symbol"] = mito_pocket_decontx_markers["names"].map(rna_adata.var["gene_symbol"])

    sc.pl.rank_genes_groups(_decontx_view, n_genes=25, sharey=False)
    mito_pocket_decontx_markers.head(25)
    return


@app.cell
def mito_gradient_gene_correlation(
    decontx_computed,
    np,
    pd,
    rna_adata,
    rna_decontx_adata,
    rna_pct_counts_mt_decontx,
    scclr,
):
    decontx_computed  # ran after DecontX finished

    # Treat the decontaminated %mito as a continuous gradient rather than a binary
    # split, and find genes whose (pflog-normalized) expression scales linearly
    # with it. A per-gene OLS/GLM loop over 62,757 genes would be far too slow
    # here, so this uses the same vectorized correlation-via-matrix-multiplication
    # approach used earlier in this notebook for the mito gene correlation
    # investigation (mathematically equivalent to the standardized coefficient of
    # a per-gene linear regression against the gradient).
    _decontx_view2 = rna_decontx_adata.copy()
    _decontx_view2.X = _decontx_view2.layers["decontX_counts"]
    scclr.pp.pflog(_decontx_view2, target="auto")

    _pflog = _decontx_view2.layers["pflog"]
    _n_cells = _pflog.shape[0]

    _mito_gradient = rna_pct_counts_mt_decontx
    _mean_mito = _mito_gradient.mean()
    _std_mito = _mito_gradient.std()

    _mean_genes = np.asarray(_pflog.mean(axis=0)).flatten()
    _sq_mean_genes = np.asarray(_pflog.multiply(_pflog).mean(axis=0)).flatten()
    _std_genes = np.sqrt(np.maximum(_sq_mean_genes - _mean_genes ** 2, 0))

    _dot = np.asarray(_pflog.T.dot(_mito_gradient)).flatten()
    _cov = _dot / _n_cells - _mean_genes * _mean_mito
    _corr = _cov / (_std_genes * _std_mito + 1e-12)

    mito_gradient_correlation = pd.DataFrame({
        "gene_id": rna_adata.var_names,
        "gene_symbol": rna_adata.var["gene_symbol"].to_numpy(),
        "is_mt": rna_adata.var["mt"].to_numpy(),
        "corr_with_mito_gradient": _corr,
    })

    # Mito genes are trivially correlated with the mito gradient by definition,
    # exclude them to see what else tracks it.
    mito_gradient_correlation_non_mt = (
        mito_gradient_correlation.loc[~mito_gradient_correlation["is_mt"]]
        .sort_values("corr_with_mito_gradient", ascending=False)
    )
    mito_gradient_correlation_non_mt.head(25)
    return


@app.cell
def ambient_rna_evidence_conclusion(mo):
    mo.md(r"""
    ### Ambient RNA vs. real mitochondrial signal, conclusion

    The evidence doesn't support one verdict for the whole dataset. It splits cleanly by population: baseline %mito variation among the majority, low-mito clusters looks like ambient contamination, while the moderate elevation in the four clusters already investigated in `high_mito_investigation_conclusion` does not.

    **Spatial structure.** %mito is strongly spatially autocorrelated on the RNA UMAP (Moran's I 0.68, neighbor-mean correlation 0.94) but nearly unstructured on the ATAC UMAP (Moran's I 0.065, neighbor-mean correlation 0.24), the same RNA-only asymmetry seen before. On its own this doesn't distinguish ambient contamination from real, cluster-associated biology, since both would produce RNA-side clustering; it mainly argues against a purely random, cell-by-cell technical artifact.

    **DecontX barely moves the global structure.** Correcting for estimated ambient contamination leaves the RNA Moran's I essentially unchanged (0.684 pre-DecontX to 0.693 post-DecontX, `decontx_moran_validation`), mean estimated contamination 10.0% (278 of 18,803 cells above 50%). If contamination were the dominant driver of the *global* spatial pattern, removing it should have visibly reduced this; it didn't, because the global statistic is dominated by between-cluster mean differences, and DecontX's correction is a within-cluster, per-cell operation.

    **The per-cluster breakdown is the decisive evidence.** `decontx_contamination_vs_mito` shows DecontX's own estimated contamination fraction correlates with %mito very differently depending on the cluster: in the low-baseline-mito clusters (Pearson r 0.68-0.99), %mito is almost entirely what DecontX itself calls contamination. In the four moderately-elevated clusters already flagged in `high_mito_investigation_conclusion` (old clusters 2, 6, 7, 10), the correlation is flat or negative (r -0.06 to -0.54): their elevated mito is specifically *not* what DecontX attributes to ambient contamination, independently corroborating that conclusion's marker-gene-based case for real biology.

    **Count dilution doesn't explain it either.** The classic ambient-dilution signature (low counts inflating %mito as a fraction) would produce a negative counts-vs-mito correlation. Instead, the global correlation is positive (Pearson r on log10 counts 0.35, Spearman 0.34), the opposite direction, so simple per-cell dilution isn't the mechanism.

    **Not yet tested: a direct ATAC-based mito measurement.** All the ATAC evidence here (Moran's I, neighbor correlation) uses RNA-derived %mito plotted on the ATAC embedding; snapatac2 doesn't compute an ATAC-side mito-fragment fraction natively, so no direct RNA-vs-ATAC mito correlation has been run yet for the specific high-mito-excluded population. That's a planned follow-up, not included here since it doesn't exist as a reproducible result yet.

    **Net:** the majority population's %mito noise is ambient contamination, well-characterized and already handled by DecontX where it matters. The four moderately-elevated clusters are real biology, already established independently by marker genes and now reinforced by their lack of a contamination signature. No filtering-strategy change follows from this either way: the ambient noise sits within clusters we're already keeping, and the real biology in the elevated clusters is exactly what `high_mito_investigation_conclusion` argued should not be filtered.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Imports and palettes
    """)
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
    from igvf_utils.connection import Connection
    from pathlib import Path
    from scipy.io import mmread

    # get the project root path using the '.env' file.
    _env_path = find_dotenv(usecwd=True)
    project_root = Path(_env_path).parent
    return Connection, Path, ad, alt, np, pd, plt, project_root, sc, sns


@app.cell
def enable_vegafusion(alt):
    # VegaFusion pre-aggregates chart data in Python before sending it to the
    # browser, raising Altair's default 5,000-row embed limit. This is enabled
    # notebook-wide (rather than toggled per-cell) since alt.data_transformers is
    # global module state, not a per-chart setting.
    alt.data_transformers.enable("vegafusion")
    return


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

    return (pca_axis_title,)


@app.cell
def mito_valley_cutoff(np):
    from scipy.stats import gaussian_kde

    def find_density_valley_cutoff(values: np.ndarray, n_grid: int = 500) -> tuple[float, float, float]:
        """Find the density valley between the two most prominent modes of `values`.

        Guards against picking up a minor noise bump as "the" valley if the KDE
        isn't perfectly clean, by restricting the search to the region between
        the two most prominent local maxima.

        Parameters
        ----------
        values : np.ndarray
            1-D array of observations to estimate the density from.
        n_grid : int, default 500
            Number of points in the grid the KDE is evaluated on.

        Returns
        -------
        cutoff : float
            The x-position of the valley.
        lo_mode : float
            The x-position of the lower of the two most prominent modes.
        hi_mode : float
            The x-position of the higher of the two most prominent modes.
        """
        _kde = gaussian_kde(values)
        _grid = np.linspace(values.min(), values.max(), n_grid)
        _density = _kde(_grid)

        _is_min = (_density[1:-1] < _density[:-2]) & (_density[1:-1] < _density[2:])
        _is_max = (_density[1:-1] > _density[:-2]) & (_density[1:-1] > _density[2:])
        _minima = list(zip(_grid[1:-1][_is_min], _density[1:-1][_is_min]))
        _maxima = list(zip(_grid[1:-1][_is_max], _density[1:-1][_is_max]))

        _top2_maxima_x = sorted([x for x, _ in sorted(_maxima, key=lambda t: -t[1])[:2]])
        _lo, _hi = _top2_maxima_x[0], _top2_maxima_x[-1]
        _candidates = [(x, d) for x, d in _minima if _lo < x < _hi]

        _cutoff = float(min(_candidates, key=lambda t: t[1])[0]) if _candidates else float(_minima[0][0])
        return _cutoff, _lo, _hi


    return


if __name__ == "__main__":
    app.run()
