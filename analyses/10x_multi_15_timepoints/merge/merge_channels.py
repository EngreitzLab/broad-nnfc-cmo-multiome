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
    # Merging channel1 and channel2

    This notebook loads the strict-QC-filtered RNA and ATAC AnnData objects produced by `channel1_qc.py` and `channel2_qc.py`, concatenates them into a single merged RNA object and a single merged ATAC object, and plots summary statistics to sanity-check the two channels are comparable (including checking for a batch effect between them). It is the starting point for trajectory analysis (RNA and ATAC), since this is a differentiation time course rather than a set of discrete populations.
    """)
    return


@app.cell(hide_code=True)
def imports():
    import altair as alt
    import anndata as ad
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import re
    import scanpy as sc
    import scclr
    import seaborn as sns
    import snapatac2 as snap
    from dotenv import find_dotenv
    from pathlib import Path
    alt.data_transformers.enable("vegafusion")

    # get the project root path using the '.env' file.
    _env_path = find_dotenv(usecwd=True)
    project_root = Path(_env_path).parent
    return Path, ad, alt, np, pd, plt, project_root, re, sc, scclr, snap, sns


@app.cell(hide_code=True)
def _(np, pd, project_root, re):
    # Okabe-Ito colorblind-safe palette (channel1 / channel2 borrow two of its
    # colors), matching the convention used in channel1_qc.py / channel2_qc.py.
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
    channel_palette = {"channel1": okabe_ito_palette[5], "channel2": okabe_ito_palette[6]}

    # Timepoint palette, same construction as timepoint_palette in channel1_qc.py /
    # channel2_qc.py: pre-branch stages get a neutral grey ramp, arterial/venous
    # branches (from d2.5 onward) get distinct warm/cool hues. Derived from the
    # full design table (all channels) rather than a single channel's
    # cmo_to_timepoint, since both channels share the same 15 timepoint labels.
    _cmo_design = pd.read_csv(project_root / "metadata/10X_multi_15_timepoints_CMO_design.tsv", sep="	")

    def _timepoint_sort_key(label):
        _match = re.match(r"d([\d.]+)(?:_(arterial|venous))?$", label)
        _day = float(_match.group(1))
        _branch = _match.group(2)
        _branch_rank = {None: 0, "arterial": 1, "venous": 2}[_branch]
        return (_day, _branch_rank)

    timepoint_order = sorted(_cmo_design["Timepoint"].unique(), key=_timepoint_sort_key)

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
    }
    return channel_palette, ec_diff_palette, okabe_ito_palette, timepoint_order


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Loading the per-channel RNA and ATAC results
    """)
    return


@app.cell
def load_rna_channels(ad, project_root):
    rna_channel1 = ad.read_h5ad(project_root / "results/10X_multiome_15_timepoints/channel1/rna_adata.h5ad")
    rna_channel2 = ad.read_h5ad(project_root / "results/10X_multiome_15_timepoints/channel2/rna_adata.h5ad")

    [rna_channel1.n_obs, rna_channel2.n_obs]
    return rna_channel1, rna_channel2


@app.cell
def load_atac_channels(project_root, snap):
    # Read with snap.read (not a generic h5ad reader) so the snapatac2-specific
    # internals (obsm["fragment_paired"], uns["reference_sequences"], etc.) come
    # back as the native snapatac2 AnnData type, not a plain anndata.AnnData.
    atac_channel1 = snap.read(project_root / "results/10X_multiome_15_timepoints/channel1/atac_adata.h5ad")
    atac_channel2 = snap.read(project_root / "results/10X_multiome_15_timepoints/channel2/atac_adata.h5ad")

    # snapatac2's Rust-side vstack (used by both snap.concat and AnnDataSet)
    # requires identical obs column *order* across components, not just identical
    # column names -- harmonize channel2 to channel1's order defensively, even
    # though channel2_qc.py mirrors channel1_qc.py's cells verbatim so the order
    # should already match.
    _col_order = atac_channel1.obs[:].columns
    atac_channel2.obs = atac_channel2.obs[:].select(_col_order)

    assert atac_channel1.obs[:].columns == atac_channel2.obs[:].columns

    [atac_channel1.n_obs, atac_channel2.n_obs]
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Concatenating the two channels
    """)
    return


@app.cell
def concat_rna(ad, rna_channel1, rna_channel2):
    # Barcodes already carry a unique per-channel suffix (e.g.
    # _10x_15timepoints_channel1 vs. _10x_15timepoints_channel2), and both
    # channels share an identical gene set (same GENCODE v43 annotation), so a
    # plain inner-join concat is safe here -- no obs_names_make_unique needed.
    rna_merged = ad.concat(
        {"channel1": rna_channel1, "channel2": rna_channel2},
        label="channel",
        join="inner",
        # var columns (gene_symbol, gene_type, gene_id, ...) are identical across
        # channels (same GENCODE v43 annotation) but ad.concat drops all
        # non-var_names var columns by default -- merge="same" keeps any column
        # that matches exactly across the concatenated objects.
        merge="same",
    )
    print(rna_merged.obs["timepoint_scanpy"].value_counts())
    rna_merged
    return (rna_merged,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Total cells per timepoint (both channels combined)
    """)
    return


@app.cell(hide_code=True)
def total_cells_per_timepoint_plot(
    alt,
    ec_diff_palette,
    rna_merged,
    timepoint_order,
):
    _counts = (
        rna_merged.obs["timepoint_scanpy"].value_counts()
        .reindex(timepoint_order).fillna(0).astype(int)
        .rename_axis("timepoint").reset_index(name="n_cells")
    )

    _bars = alt.Chart(_counts).mark_bar().encode(
        y=alt.Y("timepoint:N", sort=timepoint_order, title="Timepoint"),
        x=alt.X("n_cells:Q", title="Number of cells"),
        color=alt.Color(
            "timepoint:N",
            scale=alt.Scale(domain=list(ec_diff_palette), range=list(ec_diff_palette.values())),
            legend=None,
        ),
        tooltip=["timepoint", "n_cells"],
    )
    _labels = alt.Chart(_counts).mark_text(align="left", dx=3, fontSize=9).encode(
        y=alt.Y("timepoint:N", sort=timepoint_order),
        x=alt.X("n_cells:Q"),
        text=alt.Text("n_cells:Q", format=","),
    )

    (_bars + _labels).properties(
        title=f"Total cells assigned per timepoint (n={int(_counts['n_cells'].sum()):,}, channel1 + channel2)",
        width=550, height=400,
    ).configure_view(strokeWidth=0)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Joint normalization and clustering

    Same PFlog normalization / HVG / PCA / Leiden / UMAP pipeline as channel1_qc.py and channel2_qc.py, run once on the merged object so timepoint and channel effects can be compared directly in a shared embedding.
    """)
    return


@app.cell
def normalize_rna_merged(rna_merged, scclr):
    scclr.pp.pflog(rna_merged, target="auto")
    rna_merged.uns["pflog"]
    return


@app.cell
def hvg_rna_merged(rna_merged, sc):
    sc.pp.highly_variable_genes(rna_merged, layer="pflog", n_top_genes=2000)
    sc.pl.highly_variable_genes(rna_merged)
    return


@app.cell
def pca_rna_merged(mo, rna_merged, scclr):
    _ncomps = 50
    _ncv = 2 * _ncomps + 1
    scclr.tl.pca(rna_merged, n_comps=_ncomps, ncv=_ncv)
    mo.md(f"Computed PCA using {_ncomps} components and {_ncv} Lanczos vectors")
    return


@app.cell
def cluster_umap_rna_merged(
    ec_diff_palette,
    plt,
    rna_merged,
    sc,
    timepoint_order,
):
    sc.pp.neighbors(rna_merged, random_state=0)
    # Kept at the same default resolution used in channel1_qc.py / channel2_qc.py
    # (0.5) -- this is a differentiation time course (down to 0.5-day timepoints),
    # not discrete populations, so Leiden clusters here are only a coarse QC-level
    # view (batch-effect check, doublet-dominated clusters), not the primary tool
    # for describing the trajectory. Chasing a cluster count to match the 15
    # timepoints was tried and does not make sense for a continuous process --
    # see the trajectory analysis instead.
    sc.tl.leiden(rna_merged, flavor="igraph", resolution=0.5, n_iterations=2, random_state=0)
    sc.tl.umap(rna_merged)
    sc.pl.umap(rna_merged, color=["leiden", "channel"], ncols=2, size=5)

    # Timepoint plotted directly with matplotlib (not sc.pl.umap) so the
    # 15-category legend gets full control over layout, drawn outside the axes --
    # scanpy's default right-margin legend gets clipped with this many categories
    # (see channel1_qc.py / channel2_qc.py's umap_by_timepoint cells).
    _umap_coords = rna_merged.obsm["X_umap"]
    _timepoints = rna_merged.obs["timepoint_scanpy"].to_numpy()

    _fig, _ax = plt.subplots(figsize=(8, 6))
    for _tp in timepoint_order:
        _mask = _timepoints == _tp
        if not _mask.any():
            continue
        _ax.scatter(_umap_coords[_mask, 0], _umap_coords[_mask, 1], s=5, color=ec_diff_palette[_tp], label=_tp)
    _ax.set_xlabel("UMAP1")
    _ax.set_ylabel("UMAP2")
    _ax.set_title("UMAP colored by timepoint")
    _ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, markerscale=2, frameon=False)
    plt.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Cluster vs. timepoint composition

    Checking how well the (now higher-resolution) Leiden clusters line up with the 15 real timepoint labels, before re-checking channel representation per cluster.
    """)
    return


@app.cell(hide_code=True)
def cluster_timepoint_composition(
    alt,
    okabe_ito_palette,
    pd,
    rna_merged,
    timepoint_order,
):
    _tp_by_cluster = pd.crosstab(rna_merged.obs["leiden"].astype(str), rna_merged.obs["timepoint_scanpy"])
    _cluster_order = sorted(_tp_by_cluster.index, key=int)
    _dominant_tp = _tp_by_cluster.idxmax(axis=1)
    _dominant_pct = (_tp_by_cluster.max(axis=1) / _tp_by_cluster.sum(axis=1) * 100)

    _long = _tp_by_cluster.reset_index().melt(id_vars="leiden", var_name="timepoint", value_name="n_cells")
    _long = _long[_long["n_cells"] > 0]

    _heatmap = alt.Chart(_long).mark_rect().encode(
        x=alt.X("timepoint:N", title="Timepoint", sort=timepoint_order),
        y=alt.Y("leiden:N", title="Leiden cluster", sort=_cluster_order),
        color=alt.Color("n_cells:Q", title="Cells", scale=alt.Scale(scheme="cividis")),
        tooltip=["leiden", "timepoint", "n_cells"],
    )
    _labels_layer = alt.Chart(_long).mark_text(fontSize=9).encode(
        x=alt.X("timepoint:N", sort=timepoint_order),
        y=alt.Y("leiden:N", sort=_cluster_order),
        text=alt.Text("n_cells:Q", format=","),
        color=alt.condition(alt.datum.n_cells > _long["n_cells"].max() / 2, alt.value("white"), alt.value(okabe_ito_palette[0])),
    )

    _chart = (_heatmap + _labels_layer).properties(
        title="Leiden cluster vs. timepoint (cell counts)", width=650, height=400,
    ).configure_view(strokeWidth=0)

    print(pd.DataFrame({"dominant_timepoint": _dominant_tp, "pct_of_cluster": _dominant_pct.round(1)}).loc[_cluster_order])
    _chart
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Channel representation per cluster

    If channel1 and channel2 are well integrated (no strong batch effect), each Leiden cluster should draw from both channels in roughly the proportions they contribute overall, rather than clusters being channel-exclusive.
    """)
    return


@app.cell(hide_code=True)
def channel_representation_per_cluster(alt, channel_palette, pd, rna_merged):
    _crosstab = pd.crosstab(rna_merged.obs["leiden"].astype(str), rna_merged.obs["channel"])
    _pct = _crosstab.div(_crosstab.sum(axis=1), axis=0) * 100
    _cluster_order = sorted(_crosstab.index, key=int)

    _overall_pct = rna_merged.obs["channel"].value_counts(normalize=True) * 100

    _long = _pct.reset_index().melt(id_vars="leiden", var_name="channel", value_name="pct")
    _long["n_cells"] = _long.apply(lambda r: _crosstab.loc[r["leiden"], r["channel"]], axis=1)

    _bars = alt.Chart(_long).mark_bar().encode(
        y=alt.Y("leiden:N", sort=_cluster_order, title="Leiden cluster"),
        x=alt.X("pct:Q", title="% of cluster", scale=alt.Scale(domain=[0, 100])),
        xOffset=alt.XOffset("channel:N"),
        color=alt.Color(
            "channel:N",
            scale=alt.Scale(domain=list(channel_palette), range=list(channel_palette.values())),
        ),
        tooltip=["leiden", "channel", "n_cells", alt.Tooltip("pct:Q", format=".1f")],
    )

    # Reference lines at each channel's overall share, so a cluster's bars can be
    # compared against the expected split rather than eyeballed against 50/50.
    _ref_df = _overall_pct.rename_axis("channel").reset_index(name="pct")
    _ref_lines = alt.Chart(_ref_df).mark_rule(strokeDash=[4, 4]).encode(
        x="pct:Q",
        color=alt.Color("channel:N", scale=alt.Scale(domain=list(channel_palette), range=list(channel_palette.values())), legend=None),
    )

    (_bars + _ref_lines).properties(
        title="Channel representation per Leiden cluster (dashed lines: overall channel share)",
        width=550, height=400,
    ).configure_view(strokeWidth=0)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Trajectory analysis (RNA)

    This is a continuous differentiation time course (down to 0.5-day resolution), with one known branch point around d2.5 into arterial and venous fates -- Leiden clusters above are only a coarse QC view, not the primary description of this biology.
    """)
    return


@app.cell
def assign_branch_label(np, okabe_ito_palette, pd, rna_merged):
    # Pre-branch / Arterial / Venous, derived from the timepoint_scanpy suffix
    # (same convention as the _arterial / _venous split used to build ec_diff_palette).
    _tp = rna_merged.obs["timepoint_scanpy"].to_numpy()
    _branch = np.select(
        [np.char.endswith(_tp.astype(str), "_arterial"), np.char.endswith(_tp.astype(str), "_venous")],
        ["Arterial", "Venous"],
        default="Pre-branch",
    )

    rna_merged.obs["branch"] = pd.Categorical(_branch)
    branch_palette = {"Pre-branch": okabe_ito_palette[8], "Arterial": okabe_ito_palette[6], "Venous": okabe_ito_palette[5]}

    rna_merged.obs["branch"].value_counts()
    return (branch_palette,)


@app.cell
def compute_trajectory_neighbors(rna_merged, sc):
    sc.pp.neighbors(rna_merged, random_state=0)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Trajectory analysis via optimal transport (moscot)

    Both DPT and Palantir, run on the same PCA/diffusion structure, gave essentially the same result (~0.25 day-correlation, terminal states not clearly maximal) -- the bottleneck was the unsupervised embedding, not the pseudotime algorithm. This dataset has real, known timepoints, so instead of inferring time from expression similarity alone, `moscot`'s `TemporalProblem` computes optimal-transport couplings directly *between* consecutive real timepoints (Waddington-OT style). The arterial/venous branch is not given explicitly here; it should emerge from the transport cost (expression similarity) naturally splitting the two lineages apart once they diverge.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Does sequencing depth differ systematically across timepoints?

    Before dropping PC1 (the depth axis) from moscot's cost embedding, checking whether `total_counts` actually shifts meaningfully from one real timepoint to the next -- if it is just noise scattered independently of timepoint, it would average out in the OT coupling and not bias it much in practice.
    """)
    return


@app.cell(hide_code=True)
def check_depth_by_timepoint(
    ec_diff_palette,
    plt,
    rna_merged,
    sns,
    timepoint_order,
):
    _summary = rna_merged.obs.groupby("timepoint_scanpy", observed=True)["total_counts"].agg(["mean", "median", "std", "count"]).reindex(timepoint_order)
    print(_summary.round(1))
    print()

    from scipy.stats import f_oneway
    _groups = [g["total_counts"].to_numpy() for _, g in rna_merged.obs.groupby("timepoint_scanpy", observed=True)]
    _f_stat, _p_val = f_oneway(*_groups)
    print(f"One-way ANOVA across timepoints: F={_f_stat:.1f}, p={_p_val:.2e}")
    print(f"Between-timepoint mean range: {_summary['mean'].min():.0f} - {_summary['mean'].max():.0f} (overall std ~{rna_merged.obs['total_counts'].std():.0f})")

    _fig, _ax = plt.subplots(figsize=(10, 5))
    sns.boxplot(
        data=rna_merged.obs, x="timepoint_scanpy", y="total_counts", order=timepoint_order,
        hue="timepoint_scanpy", palette=ec_diff_palette, legend=False, ax=_ax,
    )
    _ax.set_yscale("log")
    _ax.set_xlabel(None)
    _ax.set_xticklabels(_ax.get_xticklabels(), rotation=90)
    plt.tight_layout()
    _fig
    return


@app.cell
def prepare_moscot_inputs(rna_merged):
    # PC1 correlates 0.93/0.97 with total_counts/n_genes_by_counts (technical
    # depth, not biology). PC3 correlates 0.70 with pct_counts_mt and PC4
    # correlates 0.37 with pct_counts_mt / 0.24 with channel -- channel2 runs
    # systematically hotter on mito than channel1 across nearly every timepoint
    # (a batch effect, not a per-sample defect -- see the d3.5_arterial
    # investigation above). All three would bias the optimal-transport cost the
    # same way they biased the diffusion map; dropped from the representation
    # used as moscot's joint_attr.
    _keep_pcs = [i for i in range(rna_merged.obsm["X_pca"].shape[1]) if i not in (0, 2, 3)]
    rna_merged.obsm["X_pca_cleaned"] = rna_merged.obsm["X_pca"][:, _keep_pcs]

    # moscot's sequential policy needs a numeric, orderable time key. Arterial vs.
    # venous is deliberately not encoded here -- the transport cost (expression
    # similarity in X_pca_cleaned) should separate the two lineages on its own
    # once they diverge, without being told about the branch upfront.
    rna_merged.obs["day"] = (
        rna_merged.obs["timepoint_scanpy"].astype(str).str.extract(r"d([\d.]+)")[0].astype(float)
    )
    sorted(rna_merged.obs["day"].unique())
    return


@app.cell
def prepare_temporal_problem(rna_merged):
    from moscot.problems.time import TemporalProblem

    tp = TemporalProblem(rna_merged)
    tp = tp.prepare(time_key="day", joint_attr="X_pca_cleaned")
    tp
    return (tp,)


@app.cell
def solve_temporal_problem(tp):
    # epsilon=1e-3 (moscot's raw default) left several of the 10 pairwise
    # problems not converged -- too small an entropic-regularization strength
    # for this cost scale, making the Sinkhorn iterations too peaked/unstable.
    # Bumped to 1e-2, a more typical starting point.
    tp_solved = tp.solve(epsilon=1e-2, tau_a=1.0, tau_b=1.0)
    tp_solved
    return (tp_solved,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Branch (arterial/venous) transition matrices and trajectory map

    `branch` (Pre-branch / Arterial / Venous, from `assign_branch_label`) was never given to moscot as a lineage constraint -- only real timepoints and expression similarity were used to solve the transport plan. If the couplings make biological sense, Pre-branch cells should stay Pre-branch through ~d2, then split into Arterial and Venous with each branch staying internally consistent (Arterial should not transition back to Venous, or vice versa) from d2.5 onward.
    """)
    return


@app.cell
def compute_branch_transitions(rna_merged, tp_solved):
    _days = sorted(rna_merged.obs["day"].unique())
    _day_pairs = list(zip(_days, _days[1:]))

    branch_transitions = {}
    for _src, _tgt in _day_pairs:
        branch_transitions[(_src, _tgt)] = tp_solved.cell_transition(
            _src, _tgt, source_groups="branch", target_groups="branch", forward=True, normalize=True,
        )

    for (_src, _tgt), _mat in branch_transitions.items():
        print(f"--- day {_src} -> {_tgt} ---")
        print(_mat.round(3))
        print()
    return


@app.cell
def plot_branch_sankey(rna_merged, tp_solved):
    import moscot.plotting as mpl

    tp_solved.sankey(
        source=0.0, target=5.0, source_groups="branch", target_groups="branch",
        normalize=True, key_added="sankey",
    )
    mpl.sankey(rna_merged, key="sankey", title="Branch composition across real timepoints")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Does arterial mass that "leaks" into the venous branch come back?

    Pushing `d3_arterial` mass forward to `d3.5`, splitting the pushed distribution by which branch label the landing d3.5 cells carry, then pushing just the venous-landing portion forward again to `d4` -- to tell apart transient wobble (drifts back to Arterial) from a real, sustained fate switch (stays Venous).
    """)
    return


@app.cell
def trace_arterial_leak_mass(np, tp_solved):
    # Step 1: push d3 Arterial-labeled mass forward to d3.5. push()'s output is
    # aligned to that sub-problem's own target subset (problems[...].adata_tgt),
    # not the full rna_merged -- use that subset's own obs, not global masks.
    _mass_at_35 = np.asarray(tp_solved.push(source=3.0, target=3.5, data="branch", subset="Arterial", key_added=None))
    _tgt_35 = tp_solved.problems[(3.0, 3.5)].adata_tgt
    print("push output shape:", _mass_at_35.shape, "target subset size:", _tgt_35.n_obs)

    _d35_arterial_mask = (_tgt_35.obs["branch"] == "Arterial").to_numpy()
    _d35_venous_mask = (_tgt_35.obs["branch"] == "Venous").to_numpy()

    _frac_to_arterial = _mass_at_35[_d35_arterial_mask].sum() / _mass_at_35.sum()
    _frac_to_venous = _mass_at_35[_d35_venous_mask].sum() / _mass_at_35.sum()
    print(f"d3_arterial mass landing at d3.5: {_frac_to_arterial:.1%} on Arterial-labeled cells, {_frac_to_venous:.1%} on Venous-labeled cells")

    # Step 2: take only the portion that landed on Venous-labeled d3.5 cells, push forward to d4
    # (moscot's push() returns a JAX array, which is immutable -- np.asarray above
    # already converted _mass_at_35, so a plain masked assignment works here.)
    _leaked_mass_at_35 = _mass_at_35.copy()
    _leaked_mass_at_35[~_d35_venous_mask] = 0.0

    _mass_at_4_from_leaked = np.asarray(tp_solved.push(source=3.5, target=4.0, data=_leaked_mass_at_35, key_added=None))
    _tgt_4 = tp_solved.problems[(3.5, 4.0)].adata_tgt

    _d4_arterial_mask = (_tgt_4.obs["branch"] == "Arterial").to_numpy()
    _d4_venous_mask = (_tgt_4.obs["branch"] == "Venous").to_numpy()

    _frac_back_to_arterial = _mass_at_4_from_leaked[_d4_arterial_mask].sum() / _mass_at_4_from_leaked.sum()
    _frac_stays_venous = _mass_at_4_from_leaked[_d4_venous_mask].sum() / _mass_at_4_from_leaked.sum()
    print(f"Of the d3_arterial mass that leaked onto Venous-labeled d3.5 cells, by d4: {_frac_back_to_arterial:.1%} back on Arterial-labeled cells, {_frac_stays_venous:.1%} stays on Venous-labeled cells")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Is the arterial-to-venous leak concentrated in one CMO/replicate?

    Breaking down the d3_arterial source population by its own `cmo_tag_scanpy` (its real biological replicate), then pushing each replicate's mass forward separately to check whether one specific CMO/tube accounts for most of the leak, vs. it being spread evenly across replicates (which would point to biology/embedding rather than a single mislabeled sample).
    """)
    return


@app.cell
def check_arterial_leak_by_cmo(np, pd, tp_solved):
    _src = tp_solved.problems[(3.0, 3.5)].adata_src
    _arterial_src_mask = (_src.obs["branch"] == "Arterial").to_numpy()
    _cmo_counts = _src.obs.loc[_arterial_src_mask, "cmo_tag_scanpy"].value_counts()
    print("d3_arterial source CMO composition:")
    print(_cmo_counts)
    print()
    print("channel breakdown:")
    print(pd.crosstab(_src.obs.loc[_arterial_src_mask, "cmo_tag_scanpy"], _src.obs.loc[_arterial_src_mask, "channel"]))
    print()

    _tgt_35 = tp_solved.problems[(3.0, 3.5)].adata_tgt
    _d35_venous_mask = (_tgt_35.obs["branch"] == "Venous").to_numpy()

    print("Leak fraction (pushed mass landing on Venous-labeled d3.5 cells) by source CMO:")
    for _cmo in _cmo_counts.index:
        if _cmo_counts[_cmo] < 20:
            continue  # skip tiny/noise CMOs, not real replicates
        _mass = np.asarray(tp_solved.push(source=3.0, target=3.5, data="cmo_tag_scanpy", subset=_cmo, key_added=None))
        _leak_frac = _mass[_d35_venous_mask].sum() / _mass.sum()
        print(f"  {_cmo} (n={_cmo_counts[_cmo]}): {_leak_frac:.1%}")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Are high-leak arterial cells transcriptionally more like real venous cells?

    Per-source-cell leak fraction (fraction of each d3_arterial cell's own transport mass landing on Venous-labeled d3.5 cells), then arterial/venous marker scores (derived from d4, a well-separated timepoint past the branch) compared across low-leak arterial, high-leak arterial, and real d3_venous cells -- to see whether the high-leak cells already look transcriptionally venous, or just less confidently arterial.
    """)
    return


@app.cell
def compute_leak_fraction_per_cell(np, pd, tp_solved):
    _tm = np.asarray(tp_solved.problems[(3.0, 3.5)].solution.transport_matrix)
    _src = tp_solved.problems[(3.0, 3.5)].adata_src
    _tgt_35 = tp_solved.problems[(3.0, 3.5)].adata_tgt

    _arterial_src_mask = (_src.obs["branch"] == "Arterial").to_numpy()
    _venous_tgt_mask = (_tgt_35.obs["branch"] == "Venous").to_numpy()

    _row_sums = _tm.sum(axis=1)
    _leak_fraction_all = _tm[:, _venous_tgt_mask].sum(axis=1) / _row_sums

    d3_arterial_leak_fraction = pd.Series(_leak_fraction_all[_arterial_src_mask], index=_src.obs_names[_arterial_src_mask])
    d3_arterial_leak_fraction.describe()
    return (d3_arterial_leak_fraction,)


@app.cell
def derive_branch_markers(mo, rna_merged, sc):
    _d4_mask = (rna_merged.obs["day"] == 4.0) & (rna_merged.obs["branch"].isin(["Arterial", "Venous"]))
    _d4_sub = rna_merged[_d4_mask].copy()
    sc.tl.rank_genes_groups(_d4_sub, groupby="branch", groups=["Arterial", "Venous"], layer="pflog", method="wilcoxon")

    arterial_markers = sc.get.rank_genes_groups_df(_d4_sub, group="Arterial").sort_values("scores", ascending=False)["names"].head(30).tolist()
    venous_markers = sc.get.rank_genes_groups_df(_d4_sub, group="Venous").sort_values("scores", ascending=False)["names"].head(30).tolist()

    mo.md(f"""
    Top arterial markers (d4, by Wilcoxon score): {", ".join(arterial_markers[:10])}

    Top venous markers (d4, by Wilcoxon score): {", ".join(venous_markers[:10])}
    """)
    return arterial_markers, venous_markers


@app.cell(hide_code=True)
def compare_leak_groups_to_real_venous(
    arterial_markers,
    d3_arterial_leak_fraction,
    okabe_ito_palette,
    pd,
    plt,
    rna_merged,
    sc,
    sns,
    venous_markers,
):
    sc.tl.score_genes(rna_merged, arterial_markers, score_name="arterial_marker_score", layer="pflog")
    sc.tl.score_genes(rna_merged, venous_markers, score_name="venous_marker_score", layer="pflog")

    _median_leak = d3_arterial_leak_fraction.median()
    _low_leak_ids = d3_arterial_leak_fraction[d3_arterial_leak_fraction <= _median_leak].index
    _high_leak_ids = d3_arterial_leak_fraction[d3_arterial_leak_fraction > _median_leak].index
    _real_venous_ids = rna_merged.obs_names[rna_merged.obs["timepoint_scanpy"] == "d3_venous"]

    _groups = {
        "Arterial (low leak)": _low_leak_ids,
        "Arterial (high leak)": _high_leak_ids,
        "Venous (real, d3_venous)": _real_venous_ids,
    }

    _rows = []
    _plot_rows = []
    for _label, _ids in _groups.items():
        _sub = rna_merged.obs.loc[_ids]
        _rows.append({
            "group": _label,
            "n": len(_ids),
            "median_arterial_score": round(float(_sub["arterial_marker_score"].median()), 3),
            "median_venous_score": round(float(_sub["venous_marker_score"].median()), 3),
        })
        for _metric, _col in [("arterial_marker_score", "arterial_marker_score"), ("venous_marker_score", "venous_marker_score")]:
            _plot_rows.append(pd.DataFrame({"group": _label, "metric": _metric, "value": _sub[_col].to_numpy()}))

    branch_score_summary = pd.DataFrame(_rows)
    print(branch_score_summary)

    _plot_df = pd.concat(_plot_rows, ignore_index=True)
    _fig, _axes = plt.subplots(1, 2, figsize=(12, 5))
    for _ax, _metric in zip(_axes, ["arterial_marker_score", "venous_marker_score"]):
        sns.boxplot(
            data=_plot_df[_plot_df["metric"] == _metric], x="group", y="value",
            hue="group", palette=[okabe_ito_palette[6], okabe_ito_palette[3], okabe_ito_palette[5]], legend=False, ax=_ax,
        )
        _ax.set_title(_metric)
        _ax.set_xlabel(None)
        _ax.set_xticklabels(_ax.get_xticklabels(), rotation=20, ha="right")
    plt.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Marker-based validation of the full trajectory

    Independent check using canonical, literature-known endothelial differentiation genes (not derived from this dataset): pan-endothelial/progenitor (`ETV2`, `PECAM1`, `CDH5`), arterial (`GJA5`, `EFNB2`, `DLL4`, `HEY1`, `HEY2`, `NOTCH4`), venous (`NR2F2`, `NRP2`, `EPHB4`). Expected pattern: `ETV2` (an early, transient progenitor TF) should decrease over time; `PECAM1`/`CDH5` (mature pan-EC identity) should increase and stabilize; arterial markers should stay low in Pre-branch/Venous and rise specifically in the Arterial branch after ~d2.5; venous markers should stay low in Pre-branch/Arterial and rise specifically in the Venous branch.
    """)
    return


@app.cell
def compute_literature_marker_expression(np, pd, rna_merged):
    pan_ec_markers = ["ETV2", "PECAM1", "CDH5"]
    arterial_lit_markers = ["GJA5", "EFNB2", "DLL4", "HEY1", "HEY2", "NOTCH4"]
    venous_lit_markers = ["NR2F2", "NRP2", "EPHB4"]
    all_lit_markers = pan_ec_markers + arterial_lit_markers + venous_lit_markers

    _gene_to_var = rna_merged.var.index[rna_merged.var["gene_symbol"].isin(all_lit_markers)]
    _symbol_by_var = rna_merged.var.loc[_gene_to_var, "gene_symbol"]

    _expr = rna_merged.layers["pflog"][:, rna_merged.var_names.get_indexer(_gene_to_var)]
    _expr = _expr.toarray() if hasattr(_expr, "toarray") else np.asarray(_expr)
    _expr_df = pd.DataFrame(_expr, columns=_symbol_by_var.to_numpy(), index=rna_merged.obs_names)
    _expr_df["day"] = rna_merged.obs["day"].to_numpy()
    _expr_df["branch"] = rna_merged.obs["branch"].to_numpy()

    marker_expr_by_day_branch = (
        _expr_df.melt(id_vars=["day", "branch"], var_name="gene", value_name="pflog_expr")
        .groupby(["gene", "day", "branch"], observed=True)["pflog_expr"]
        .mean()
        .reset_index()
    )
    marker_expr_by_day_branch.head()
    return all_lit_markers, marker_expr_by_day_branch


@app.cell(hide_code=True)
def plot_literature_marker_trends(
    all_lit_markers,
    branch_palette,
    marker_expr_by_day_branch,
    okabe_ito_palette,
    plt,
):
    _fig, _axes = plt.subplots(3, 4, figsize=(18, 11), sharex=True)
    _axes = _axes.flatten()

    for _ax, _gene in zip(_axes, all_lit_markers):
        _sub = marker_expr_by_day_branch[marker_expr_by_day_branch["gene"] == _gene]
        for _br, _color in branch_palette.items():
            _br_sub = _sub[_sub["branch"] == _br].sort_values("day")
            if _br_sub.empty:
                continue
            _ax.plot(_br_sub["day"], _br_sub["pflog_expr"], marker="o", markersize=3, color=_color, label=_br)
        _ax.set_title(_gene, fontsize=10)
        _ax.axvline(2.5, color=okabe_ito_palette[0], linestyle="--", linewidth=0.8)

    for _ax in _axes[len(all_lit_markers):]:
        _ax.axis("off")

    _axes[0].legend(fontsize=8)
    _fig.supxlabel("Real day")
    _fig.supylabel("Mean PFlog expression")
    _fig.suptitle("Canonical marker expression by day and branch (dashed line: branch point ~d2.5)", y=1.0)
    plt.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### CDH5 / PECAM1 are dropout-dominated, not smoothly declining

    The "mean expression declines at the terminals" pattern seen in the literature-marker validation above is actually a change in *detection rate*: both genes have a median (and even 10th percentile) of exactly 0 at every single real timepoint. Correlation with `pct_counts_mt` (0.11) and `total_counts` (-0.04) is weak, and the detection-rate pattern does not track sequencing depth in the direction a pure technical dropout would predict (depth is actually *higher* at d0 and d4.5-d5, exactly where CDH5 detection is worst) -- so this is not simply "fewer reads, more missed transcripts."
    """)
    return


@app.cell(hide_code=True)
def check_cdh5_pecam1_dropout(np, pd, rna_merged):
    _cdh5_var_id = rna_merged.var.index[rna_merged.var["gene_symbol"] == "CDH5"][0]
    _pecam1_var_id = rna_merged.var.index[rna_merged.var["gene_symbol"] == "PECAM1"][0]

    def _get_expr(var_id):
        _e = rna_merged.layers["pflog"][:, rna_merged.var_names.get_loc(var_id)]
        return np.asarray(_e.todense()).ravel() if hasattr(_e, "todense") else np.asarray(_e).ravel()

    _dropout_df = pd.DataFrame({
        "day": rna_merged.obs["day"].to_numpy(),
        "cdh5": _get_expr(_cdh5_var_id),
        "pecam1": _get_expr(_pecam1_var_id),
    })

    cdh5_pecam1_dropout_by_day = _dropout_df.groupby("day").agg(
        cdh5_frac_zero=("cdh5", lambda s: (s == 0).mean()),
        cdh5_p90=("cdh5", lambda s: s.quantile(0.90)),
        pecam1_frac_zero=("pecam1", lambda s: (s == 0).mean()),
        pecam1_p90=("pecam1", lambda s: s.quantile(0.90)),
    ).round(3)
    cdh5_pecam1_dropout_by_day
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Non-endothelial contamination check

    CDH5/PECAM1 detection rate rises then falls across the trajectory (peaking ~d2.5-3.5, dropping to ~10-12% of cells by d4.5-d5). Checking whether CDH5-negative cells are just noisy/low-quality EC cells, or a genuinely distinct, non-endothelial population.
    """)
    return


@app.cell
def derive_ec_vs_non_ec_markers(mo, rna_merged, sc):
    rna_merged.obs["cdh5_positive"] = _get_expr(_cdh5_var_id) > 0

    # QC comparison: CDH5-negative cells are not lower-depth/higher-mito (which
    # would suggest plain dropout noise) -- if anything, slightly the opposite.
    print("QC metrics by CDH5 detection status:")
    print(rna_merged.obs.groupby("cdh5_positive")[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].median().round(1))
    print()

    # Restricted to late timepoints (day >= 4), where detection is worst, to
    # characterize what actually distinguishes CDH5-negative cells there.
    _late_mask = (rna_merged.obs["day"].astype(float) >= 4.0).to_numpy()
    _late_sub = rna_merged[_late_mask].copy()
    _late_sub.obs["cdh5_positive"] = _late_sub.obs["cdh5_positive"].astype(str).astype("category")

    sc.tl.rank_genes_groups(_late_sub, groupby="cdh5_positive", groups=["True"], layer="pflog", method="wilcoxon")
    _de_res = sc.get.rank_genes_groups_df(_late_sub, group="True").sort_values("scores", ascending=False)

    _var_symbols = rna_merged.var["gene_symbol"]
    ec_identity_markers = _de_res.head(20)["names"].map(_var_symbols).dropna().tolist()
    non_ec_markers = _de_res.tail(20)["names"].map(_var_symbols).dropna().tolist()

    mo.md(f"""
    **Top genes up in CDH5-positive cells** (day >= 4, Wilcoxon): {", ".join(ec_identity_markers[:12])} -- canonical pan-endothelial identity genes (`KDR`/`FLI1`/`FLT1`/`ERG` are among the most established EC master regulators in the literature).

    **Top genes up in CDH5-negative cells**: {", ".join(non_ec_markers[:14])} -- a neuronal/mesenchymal-like signature (`ROBO2`/`UNC5C`/`RIMS2` axon guidance/synaptic genes, `VCAN` mesenchymal ECM, `LEF1` non-EC progenitor/mesenchymal), not "the same EC cells with lower CDH5."
    """)
    return ec_identity_markers, non_ec_markers


@app.cell
def score_ec_vs_non_ec_identity(
    ec_identity_markers,
    non_ec_markers,
    rna_merged,
    sc,
):
    # score_genes needs var_names (Ensembl IDs here), not gene symbols directly.
    # ec_identity_markers / non_ec_markers come from derive_ec_vs_non_ec_markers above.
    _ec_var_ids = rna_merged.var.index[rna_merged.var["gene_symbol"].isin(ec_identity_markers)]
    _non_ec_var_ids = rna_merged.var.index[rna_merged.var["gene_symbol"].isin(non_ec_markers)]

    sc.tl.score_genes(rna_merged, _ec_var_ids, score_name="ec_identity_score", layer="pflog")
    sc.tl.score_genes(rna_merged, _non_ec_var_ids, score_name="non_ec_score", layer="pflog")

    rna_merged.obs["is_non_ec_like"] = rna_merged.obs["non_ec_score"] > rna_merged.obs["ec_identity_score"]
    rna_merged.obs["is_non_ec_like"].value_counts()
    return


@app.cell(hide_code=True)
def plot_non_ec_fraction(branch_palette, okabe_ito_palette, plt, rna_merged):

    _contam_by_day_branch = (
        rna_merged.obs.groupby(["day", "branch"], observed=True)["is_non_ec_like"]
        .agg(["mean", "size"])
        .rename(columns={"mean": "frac_non_ec", "size": "n_cells"})
        .reset_index()
    )
    print(_contam_by_day_branch.to_string(index=False))

    _fig, _ax = plt.subplots(figsize=(9, 5))
    for _br, _color in branch_palette.items():
        _sub = _contam_by_day_branch[_contam_by_day_branch["branch"] == _br].sort_values("day")
        if _sub.empty:
            continue
        _ax.plot(_sub["day"], _sub["frac_non_ec"] * 100, marker="o", color=_color, label=_br)
    _ax.axvline(2.5, color=okabe_ito_palette[0], linestyle="--", linewidth=0.8)
    _ax.set_xlabel("Real day")
    _ax.set_ylabel("% non-EC-like cells")
    _ax.set_title("Non-endothelial-like cell fraction across the trajectory")
    _ax.legend()
    plt.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Is the d4+ non-EC spike a smooth decline, or a specific replicate failure?

    Breaking the non-EC fraction down by real biological replicate (`cmo_tag_scanpy` x `channel`) for every timepoint from d4 onward.
    """)
    return


@app.cell(hide_code=True)
def check_late_replicate_non_ec_breakdown(pd, rna_merged):
    _late_timepoints = ["d4_arterial", "d4_venous", "d4.5_venous", "d5_venous"]
    _rows = []
    for _tp in _late_timepoints:
        _sub = rna_merged.obs[rna_merged.obs["timepoint_scanpy"] == _tp]
        _ct = pd.crosstab([_sub["cmo_tag_scanpy"], _sub["channel"]], _sub["is_non_ec_like"])
        _ct = _ct[_ct.sum(axis=1) >= 20]  # drop tiny/noise CMO assignments, not real replicates
        for (_cmo, _channel), _row in _ct.iterrows():
            _n = int(_row.sum())
            _rows.append({
                "timepoint": _tp, "cmo": _cmo, "channel": _channel, "n_cells": _n,
                "pct_non_ec": round(float(_row.get(True, 0) / _n * 100), 1),
            })

    late_replicate_non_ec_breakdown = pd.DataFrame(_rows)
    late_replicate_non_ec_breakdown
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Summary: non-endothelial contamination

    - **d0-d1 (75-88% non-EC) is expected, not contamination**: these are pre-EC-commitment progenitor cells, which should not yet express CDH5/KDR/FLI1/ERG.
    - **d2.5-d3.5 Arterial is the cleanest window in the whole trajectory** (3-7% non-EC), matching the strong, unambiguous arterial marker validation found earlier.
    - **Venous is consistently worse than Arterial even in that same clean window** (23-28% vs 3-7%) -- both branches are meant to be fully differentiated ECs (arterial vs. venous is a subtype distinction layered on a shared EC identity, not EC-vs-non-EC), so this gap is itself a real signal, not expected biology.
    - **The d4+ spike (57-79% non-EC) is not a smooth biological decline -- it is a specific, near-total replicate failure concentrated in channel2**: every single channel2 replicate from d4 onward is 96-100% non-EC, with no exceptions, while channel1 replicates over the same window show a milder (30-60%) decline. The same CMO tag (CMO18) is 99.1% clean in channel1 and 96.4% failed in channel2 at d4_arterial -- ruling out a CMO/labeling explanation in favor of a channel2-specific technical failure (bad nuclei prep, sample degradation, etc.) affecting its d4+ samples specifically.
    - channel1's own milder late decline (30-60%) is a separate, smaller effect still worth attention, but is not on the same scale as channel2's near-total failure.

    **Practical implication:** channel2 data from d4 onward should likely be excluded from the trajectory analysis rather than fixed via further transport-plan tuning; channel1's later timepoints warrant a closer look but are not obviously unusable.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Independent check: projecting 5-timepoints d4 cells onto this trajectory

    Two possible explanations for the d4+ non-endothelial spike: (1) a channel2-specific technical failure, or (2) an intrinsic feature of this differentiation protocol (real EC identity loss/instability by this stage, independent of which experiment or channel). To tell them apart: the completely independent 5-timepoints dataset's own day-4 cells (a different set of channels, a different sequencing run, expected to be genuine ECs by that experiment's own design) are projected onto this 15-timepoints trajectory's PCA space and assigned a predicted day via kNN. If they land in a tight window, the failure here is more likely technical (channel2-specific); if they spread broadly across the whole trajectory, day-4-ish instability may be a real feature of the protocol itself, independent of channel.
    """)
    return


@app.cell
def project_5tp_d4_onto_trajectory(ad, np, rna_merged):
    d4_5tp = ad.read_h5ad("/Users/emattei/GitHub/broad-nnfc-cmo-multiome/results/10X_multiome_5_timepoints/merged/rna/rna_d4.h5ad")

    # Match genes by symbol: our var_names are Ensembl IDs with a gene_symbol
    # column; the 5tp file's var_names are already gene symbols. Restrict to
    # symbols that map 1:1 in our data (a small number of Ensembl IDs share a
    # legacy symbol -- see channel1_qc.py's save_rna_per_timepoint note -- those
    # are dropped rather than arbitrarily picking one).
    _symbol_counts = rna_merged.var["gene_symbol"].value_counts()
    _unique_symbols = set(_symbol_counts[_symbol_counts == 1].index)
    _shared_symbols = [g for g in d4_5tp.var_names if g in _unique_symbols]
    print(f"{len(_shared_symbols)} of {d4_5tp.n_vars} 5tp genes matched 1:1 to this dataset's genes")

    _our_var_ids = rna_merged.var.index[rna_merged.var["gene_symbol"].isin(_shared_symbols)]
    _our_gene_order = rna_merged.var.loc[_our_var_ids, "gene_symbol"]
    # Reorder our var ids to match _shared_symbols exactly
    _our_var_ids_ordered = _our_gene_order.reset_index().set_index("gene_symbol").loc[_shared_symbols, "index"]

    _pcs_shared = rna_merged.varm["PCs"][rna_merged.var_names.get_indexer(_our_var_ids_ordered), :]

    _new_pflog = d4_5tp[:, _shared_symbols].layers["pflog"]
    _new_pflog = np.asarray(_new_pflog.todense()) if hasattr(_new_pflog, "todense") else np.asarray(_new_pflog)
    _new_center = d4_5tp.obs["pflog_center"].to_numpy()

    # Same shifted-CLR projection scclr.tl.pca uses (subtract each cell's own
    # pflog_center, then project through the reference PCA loadings), restricted
    # to the shared gene subset.
    _5tp_d4_pca = (_new_pflog - _new_center[:, None]) @ _pcs_shared
    projected_5tp_d4_pca_cleaned = _5tp_d4_pca[:, [i for i in range(50) if i not in (0, 2, 3)]]
    projected_5tp_d4_pca_cleaned.shape
    return (projected_5tp_d4_pca_cleaned,)


@app.cell(hide_code=True)
def predict_day_for_projected_5tp_d4(
    np,
    okabe_ito_palette,
    pd,
    plt,
    projected_5tp_d4_pca_cleaned,
    rna_merged,
):
    from sklearn.neighbors import NearestNeighbors

    _nn = NearestNeighbors(n_neighbors=15).fit(rna_merged.obsm["X_pca_cleaned"])
    _dists, _idx = _nn.kneighbors(projected_5tp_d4_pca_cleaned)

    _ref_day = rna_merged.obs["day"].to_numpy()
    _ref_branch = rna_merged.obs["branch"].to_numpy()

    # Distance-weighted kNN regression for day, majority vote for branch.
    _weights = 1.0 / (_dists + 1e-6)
    _predicted_day = (_ref_day[_idx] * _weights).sum(axis=1) / _weights.sum(axis=1)

    def _majority_branch(row_idx):
        _vals, _counts = np.unique(_ref_branch[row_idx], return_counts=True)
        return _vals[np.argmax(_counts)]

    _predicted_branch = np.array([_majority_branch(row) for row in _idx])

    projected_5tp_d4 = pd.DataFrame({"predicted_day": _predicted_day, "predicted_branch": _predicted_branch})

    _day_range = projected_5tp_d4["predicted_day"].max() - projected_5tp_d4["predicted_day"].min()
    _day_iqr = projected_5tp_d4["predicted_day"].quantile(0.75) - projected_5tp_d4["predicted_day"].quantile(0.25)

    print(projected_5tp_d4["predicted_day"].describe())
    print()
    print(f"Range: {_day_range:.2f} days")
    print(f"IQR: {_day_iqr:.2f} days")
    print()
    print("predicted branch breakdown:")
    print(projected_5tp_d4["predicted_branch"].value_counts())

    _fig, _ax = plt.subplots(figsize=(9, 4))
    _ax.hist(projected_5tp_d4["predicted_day"], bins=np.arange(-0.25, 5.5, 0.5), color=okabe_ito_palette[5], edgecolor="white")
    _ax.set_xlabel("Predicted day (kNN onto this trajectory)")
    _ax.set_ylabel("Number of 5tp d4 cells")
    _ax.set_title("Where do independent, known-EC 5tp d4 cells land on this trajectory?")
    plt.tight_layout()
    _fig
    return (NearestNeighbors,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Conclusion: technical, not protocol-intrinsic

    The independent 5-timepoints dataset's own day-4 cells (a separate experiment, separate channels, expected to be genuine ECs by that experiment's own design) were projected onto this 15-timepoints trajectory via a shared-gene PCA projection (matching genes by symbol, using this dataset's own PCA loadings) followed by distance-weighted kNN regression against `day` and `branch`.

    **Result:** they land in an extremely tight cluster, not a spread -- range 1.40 days (2.31-3.71), IQR 0.23 days, with 85% of cells (4,289/5,038) falling in a single 0.5-day bin, and predicted branch 99.9% Arterial (5,031/5,038).

    If day-4 EC identity were intrinsically unstable as a feature of this differentiation protocol, an independent d4 population from a different experiment should have spread broadly across the trajectory when projected here. Instead it lands in a small, coherent, single-branch region -- consistent with a stable, well-defined transcriptional state. (Landing nearer "day 3" than "day 4," and reading as Arterial rather than split, likely reflects a timing/kinetics offset between the two separate protocols and Arterial being the cleaner, better-defined branch throughout this analysis -- neither undermines the main conclusion.)

    **This favors the technical explanation over the protocol-intrinsic one.** The d4+ non-endothelial spike found in this dataset -- channel2's near-total failure (96-100% non-EC in every replicate) and channel1's milder decline (30-60%) -- is more likely a processing/technical issue specific to this experiment's later timepoints than an inherent limitation of the differentiation protocol itself.

    **Next step (not yet applied):** exclude channel2 from day >= 4 before re-running the trajectory; channel1's milder late decline may warrant separate scrutiny before deciding whether it also needs exclusion.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Cleaning: removing technical artifacts before re-deriving pseudotime

    Excluding (1) all channel2 cells at day >= 4 (the near-total, wholesale technical failure found above -- even its EC-passing cells are suspect given how compromised that batch looks overall), and (2) any cell at day >= 3.5 flagged `is_non_ec_like` regardless of channel (catches channel1's milder late decline too). The non-EC filter starts at 3.5, not right at the branch point (2.5), since the non-EC fraction is still low and not clearly abnormal at d2.5-3.0 (3-28%) and only becomes a clear signal from d3.5 onward, sharply so by d4 (57-80%). Earlier/pre-branch cells are not filtered by EC score at all -- scoring non-EC there is expected biology (pre-differentiation progenitors), not contamination.
    """)
    return


@app.cell
def filter_technical_artifacts(rna_merged):
    _is_bad_channel2_late = (rna_merged.obs["channel"] == "channel2") & (rna_merged.obs["day"].astype(float) >= 4.0)
    _is_non_ec_from_3_5_onward = (rna_merged.obs["day"].astype(float) >= 3.5) & rna_merged.obs["is_non_ec_like"]

    pass_ec_qc_mask = ~(_is_bad_channel2_late | _is_non_ec_from_3_5_onward)
    rna_clean = rna_merged[pass_ec_qc_mask].copy()

    print(f"{rna_clean.n_obs:,} of {rna_merged.n_obs:,} cells kept ({rna_clean.n_obs / rna_merged.n_obs:.1%})")
    print(f"  excluded as bad channel2-late: {int(_is_bad_channel2_late.sum()):,}")
    print(f"  excluded as non-EC from day 3.5 onward: {int(_is_non_ec_from_3_5_onward.sum()):,}")
    print(f"  (overlap between the two: {int((_is_bad_channel2_late & _is_non_ec_from_3_5_onward).sum()):,})")
    return (rna_clean,)


@app.cell(hide_code=True)
def cleaning_report_by_day_branch_channel(pd, rna_clean, rna_merged):
    _before = rna_merged.obs.groupby(["day", "branch", "channel"], observed=True).size().rename("n_before")
    _after = rna_clean.obs.groupby(["day", "branch", "channel"], observed=True).size().rename("n_after")

    cleaning_report = pd.concat([_before, _after], axis=1).fillna(0).astype(int)
    cleaning_report["pct_kept"] = (cleaning_report["n_after"] / cleaning_report["n_before"] * 100).round(1)
    cleaning_report
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Re-embedding the cleaned population

    Reusing the existing `pflog` layer/centering (a per-cell transform, unaffected by which cells are kept), but recomputing HVGs and PCA fresh on `rna_clean` -- removing ~half the population (mostly the non-EC contamination) will genuinely shift which genes are highly variable and what the leading components represent.
    """)
    return


@app.cell
def reembed_clean_pca(mo, rna_clean, sc, scclr):
    sc.pp.highly_variable_genes(rna_clean, layer="pflog", n_top_genes=2000)

    _ncomps = 50
    _ncv = 2 * _ncomps + 1
    scclr.tl.pca(rna_clean, n_comps=_ncomps, ncv=_ncv)
    mo.md(f"Computed PCA on rna_clean using {_ncomps} components and {_ncv} Lanczos vectors")
    return


@app.cell(hide_code=True)
def check_clean_pca_correlations(np, pd, rna_clean):
    _pca = rna_clean.obsm["X_pca"][:, :10]
    _tc = rna_clean.obs["total_counts"].to_numpy()
    _mt = rna_clean.obs["pct_counts_mt"].to_numpy()
    _ch = (rna_clean.obs["channel"] == "channel2").to_numpy().astype(float)
    _day = rna_clean.obs["day"].to_numpy()

    _rows = []
    for _i in range(10):
        _pc = _pca[:, _i]
        _rows.append({
            "PC": _i + 1,
            "variance_pct": round(float(rna_clean.uns["pca"]["variance_ratio"][_i] * 100), 2),
            "vs_total_counts": round(float(np.corrcoef(_pc, _tc)[0, 1]), 3),
            "vs_pct_counts_mt": round(float(np.corrcoef(_pc, _mt)[0, 1]), 3),
            "vs_channel": round(float(np.corrcoef(_pc, _ch)[0, 1]), 3),
            "vs_day": round(float(np.corrcoef(_pc, _day)[0, 1]), 3),
        })
    pd.DataFrame(_rows).set_index("PC")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Re-deriving pseudotime with Palantir on the cleaned data

    Using all 50 PCs (no components dropped this time -- depth is entangled with real day in the cleaned population, not purely technical, see the check above). Same palantir workflow as before: diffusion maps, root selection (most extreme d0 cell along the DC that best separates d0 from the rest), terminal states (d4_arterial and d5_venous, same logic), then run_palantir.
    """)
    return


@app.cell
def compute_clean_neighbors(rna_clean, sc):
    sc.pp.neighbors(rna_clean, random_state=0)
    return


@app.cell
def _():
    import palantir

    return (palantir,)


@app.cell
def run_clean_diffusion_maps(mo, palantir, rna_clean):
    _dm_res = palantir.utils.run_diffusion_maps(rna_clean, n_components=10, knn=30)
    mo.md(f"Diffusion maps computed: {rna_clean.obsm['DM_EigenVectors'].shape[1]} components.")
    return


@app.cell
def select_clean_palantir_root(mo, np, rna_clean):
    # Find the diffusion component that best separates d0 (the known root
    # population) from the rest, then pick the most extreme d0 cell along that
    # component -- same logic as the earlier (pre-cleaning) root selection.
    _dcs = rna_clean.obsm["DM_EigenVectors"]
    _d0_mask = (rna_clean.obs["timepoint_scanpy"] == "d0").to_numpy()

    _separation = np.abs(_dcs[_d0_mask, 1:].mean(axis=0) - _dcs[~_d0_mask, 1:].mean(axis=0))
    _root_dc = int(np.argmax(_separation)) + 1
    _direction = np.sign(_dcs[_d0_mask, _root_dc].mean() - _dcs[~_d0_mask, _root_dc].mean())

    _d0_indices = np.flatnonzero(_d0_mask)
    _root_local_idx = np.argmax(_direction * _dcs[_d0_mask, _root_dc])
    _root_idx = int(_d0_indices[_root_local_idx])
    palantir_root_cell_clean = rna_clean.obs_names[_root_idx]

    mo.md(
        f"Root cell: `{palantir_root_cell_clean}` (timepoint `{rna_clean.obs['timepoint_scanpy'].iloc[_root_idx]}`), "
        f"the most extreme d0 cell along DC{_root_dc}."
    )
    return (palantir_root_cell_clean,)


@app.cell
def select_clean_palantir_terminal_states(mo, np, pd, rna_clean):
    def _pick_extreme_cell(group_mask, dcs):
        _separation = np.abs(dcs[group_mask, 1:].mean(axis=0) - dcs[~group_mask, 1:].mean(axis=0))
        _dc = int(np.argmax(_separation)) + 1
        _direction = np.sign(dcs[group_mask, _dc].mean() - dcs[~group_mask, _dc].mean())
        _indices = np.flatnonzero(group_mask)
        _local_idx = np.argmax(_direction * dcs[group_mask, _dc])
        return int(_indices[_local_idx]), _dc

    _dcs = rna_clean.obsm["DM_EigenVectors"]
    _arterial_mask = (rna_clean.obs["timepoint_scanpy"] == "d4_arterial").to_numpy()
    _venous_mask = (rna_clean.obs["timepoint_scanpy"] == "d5_venous").to_numpy()

    _arterial_idx, _arterial_dc = _pick_extreme_cell(_arterial_mask, _dcs)
    _venous_idx, _venous_dc = _pick_extreme_cell(_venous_mask, _dcs)

    palantir_terminal_states_clean = pd.Series({
        rna_clean.obs_names[_arterial_idx]: "Arterial",
        rna_clean.obs_names[_venous_idx]: "Venous",
    })

    mo.md(f"""
    Terminal states:

    - Arterial: `{rna_clean.obs_names[_arterial_idx]}` (most extreme d4_arterial cell along DC{_arterial_dc})
    - Venous: `{rna_clean.obs_names[_venous_idx]}` (most extreme d5_venous cell along DC{_venous_dc})
    """)
    return (palantir_terminal_states_clean,)


@app.cell
def clean_palantir_multiscale_space(mo, palantir, rna_clean):
    palantir.utils.determine_multiscale_space(rna_clean)
    mo.md(f"Multiscale space computed: {rna_clean.obsm['DM_EigenVectors_multiscaled'].shape[1]} components.")
    return


@app.cell
def run_clean_palantir(
    mo,
    palantir,
    palantir_root_cell_clean,
    palantir_terminal_states_clean,
    rna_clean,
):
    palantir.core.run_palantir(
        rna_clean,
        early_cell=palantir_root_cell_clean,
        terminal_states=palantir_terminal_states_clean,
        knn=30,
    )
    mo.md("Palantir pseudotime, entropy, and fate probabilities computed on the cleaned population.")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## QC and pseudotime plots (cleaned data)
    """)
    return


@app.cell
def compute_clean_umap(mo, rna_clean, sc):
    # Overwrites the stale X_umap inherited from rna_merged (via the .copy()
    # slice) with a fresh embedding computed on rna_clean's own neighbor graph.
    sc.tl.umap(rna_clean)
    mo.md("UMAP computed on the cleaned population's own neighbor graph.")
    return


@app.cell(hide_code=True)
def plot_qc_by_day_clean(plt, rna_clean, sns):
    _qc_metrics = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    _fig, _axes = plt.subplots(1, len(_qc_metrics), figsize=(5 * len(_qc_metrics), 4.5))
    for _ax, _metric in zip(_axes, _qc_metrics):
        sns.boxplot(
            data=rna_clean.obs, x="day", y=_metric,
            hue="day", palette="cividis", legend=False, ax=_ax,
        )
        _ax.set_title(_metric)
        _ax.set_xlabel("Real day")
    plt.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def plot_clean_umap_panels(
    branch_palette,
    ec_diff_palette,
    plt,
    rna_clean,
    timepoint_order,
):
    _umap_coords = rna_clean.obsm["X_umap"]

    _fig, _axes = plt.subplots(1, 2, figsize=(14, 6))

    _sc_pt = _axes[0].scatter(_umap_coords[:, 0], _umap_coords[:, 1], s=5, c=rna_clean.obs["palantir_pseudotime"], cmap="cividis")
    _axes[0].set_xlabel("UMAP1")
    _axes[0].set_ylabel("UMAP2")
    _axes[0].set_title("Colored by palantir_pseudotime")
    _fig.colorbar(_sc_pt, ax=_axes[0], label="pseudotime")

    for _br, _color in branch_palette.items():
        _mask = (rna_clean.obs["branch"] == _br).to_numpy()
        _axes[1].scatter(_umap_coords[_mask, 0], _umap_coords[_mask, 1], s=5, color=_color, label=_br)
    _axes[1].set_xlabel("UMAP1")
    _axes[1].set_ylabel("UMAP2")
    _axes[1].set_title("Colored by branch")
    _axes[1].legend()

    plt.tight_layout()
    _fig

    # Timepoint panel separately, with its own external legend (15 categories --
    # see the earlier notebooks' umap_by_timepoint cells for why this needs a
    # manual legend instead of scanpy's default).
    _fig2, _ax2 = plt.subplots(figsize=(8, 6))
    for _tp in timepoint_order:
        if _tp == "Negative":
            continue
        _mask = (rna_clean.obs["timepoint_scanpy"] == _tp).to_numpy()
        if not _mask.any():
            continue
        _ax2.scatter(_umap_coords[_mask, 0], _umap_coords[_mask, 1], s=5, color=ec_diff_palette[_tp], label=_tp)
    _ax2.set_xlabel("UMAP1")
    _ax2.set_ylabel("UMAP2")
    _ax2.set_title("Colored by timepoint")
    _ax2.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, markerscale=2, frameon=False)
    plt.tight_layout()
    _fig2
    return


@app.cell(hide_code=True)
def plot_clean_pseudotime_by_timepoint(
    ec_diff_palette,
    plt,
    rna_clean,
    sns,
    timepoint_order,
):
    _fig, _ax = plt.subplots(figsize=(10, 5))
    sns.boxplot(
        data=rna_clean.obs, x="timepoint_scanpy", y="palantir_pseudotime",
        order=[t for t in timepoint_order if t != "Negative"],
        hue="timepoint_scanpy", palette=ec_diff_palette, legend=False, ax=_ax,
    )
    _ax.set_xlabel(None)
    _ax.set_xticklabels(_ax.get_xticklabels(), rotation=90)
    _ax.set_title("palantir_pseudotime by real timepoint (cleaned population)")
    plt.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Validating the new pseudotime against the full 5-timepoints experiment (d0-d4)

    The 5-timepoints dataset spans the same early window (d0 through d4 EC) as an independent experiment. Loading all five of its per-timepoint files, projecting them into `rna_clean`'s PCA space (all 50 PCs, matching what the pseudotime was actually computed on -- no components dropped this time), then transferring `palantir_pseudotime` via kNN and checking it against the 5tp's own known real day labels.
    """)
    return


@app.cell
def load_5tp_all_timepoints(Path, ad):
    _5tp_dir = Path("/Users/emattei/GitHub/broad-nnfc-cmo-multiome/results/10X_multiome_5_timepoints/merged/rna")
    _5tp_files = {0.0: "rna_d0.h5ad", 1.0: "rna_d1.h5ad", 2.0: "rna_d2.h5ad", 3.0: "rna_d3.h5ad", 4.0: "rna_d4.h5ad"}

    _5tp_adatas = []
    for _day_val, _fname in _5tp_files.items():
        _a = ad.read_h5ad(_5tp_dir / _fname)
        _a.obs["day_5tp"] = _day_val
        _5tp_adatas.append(_a)

    five_tp_all = ad.concat(_5tp_adatas, join="inner")
    five_tp_all.obs["day_5tp"].value_counts()
    return (five_tp_all,)


@app.cell
def project_5tp_all_onto_clean_pca(five_tp_all, np, rna_clean):
    # Same gene-matching + shifted-CLR projection approach as the earlier d4-only
    # check, but onto rna_clean's PCA (all 50 PCs -- no components dropped this
    # time, since depth turned out to be entangled with real biology here).
    _symbol_counts = rna_clean.var["gene_symbol"].value_counts()
    _unique_symbols = set(_symbol_counts[_symbol_counts == 1].index)
    _shared_symbols = [g for g in five_tp_all.var_names if g in _unique_symbols]
    print(f"{len(_shared_symbols)} of {five_tp_all.n_vars} 5tp genes matched 1:1 to this dataset's genes")

    _our_var_ids = rna_clean.var.index[rna_clean.var["gene_symbol"].isin(_shared_symbols)]
    _our_gene_order = rna_clean.var.loc[_our_var_ids, "gene_symbol"]
    _our_var_ids_ordered = _our_gene_order.reset_index().set_index("gene_symbol").loc[_shared_symbols, "index"]

    _pcs_shared = rna_clean.varm["PCs"][rna_clean.var_names.get_indexer(_our_var_ids_ordered), :]

    _new_pflog = five_tp_all[:, _shared_symbols].layers["pflog"]
    _new_pflog = np.asarray(_new_pflog.todense()) if hasattr(_new_pflog, "todense") else np.asarray(_new_pflog)
    _new_center = five_tp_all.obs["pflog_center"].to_numpy()

    projected_5tp_all_pca = (_new_pflog - _new_center[:, None]) @ _pcs_shared
    projected_5tp_all_pca.shape
    return (projected_5tp_all_pca,)


@app.cell
def predict_pseudotime_for_5tp_all(
    NearestNeighbors,
    five_tp_all,
    np,
    projected_5tp_all_pca,
    rna_clean,
):
    _nn = NearestNeighbors(n_neighbors=15).fit(rna_clean.obsm["X_pca"])
    _dists, _idx = _nn.kneighbors(projected_5tp_all_pca)

    _ref_pt = rna_clean.obs["palantir_pseudotime"].to_numpy()
    _weights = 1.0 / (_dists + 1e-6)
    _predicted_pt = (_ref_pt[_idx] * _weights).sum(axis=1) / _weights.sum(axis=1)

    five_tp_all.obs["predicted_pseudotime"] = _predicted_pt

    _corr = np.corrcoef(five_tp_all.obs["predicted_pseudotime"], five_tp_all.obs["day_5tp"])[0, 1]
    print(f"predicted_pseudotime vs real 5tp day correlation: {_corr:.3f}")
    print()
    print(five_tp_all.obs.groupby("day_5tp")["predicted_pseudotime"].describe().round(3))
    return


@app.cell(hide_code=True)
def plot_5tp_all_pseudotime_validation(five_tp_all, plt, sns):
    _fig, _ax = plt.subplots(figsize=(8, 5))
    sns.boxplot(
        data=five_tp_all.obs, x="day_5tp", y="predicted_pseudotime",
        hue="day_5tp", palette="cividis", legend=False, ax=_ax,
    )
    _ax.set_xlabel("Real day (5-timepoints dataset)")
    _ax.set_ylabel("Predicted pseudotime (kNN transfer)")
    _ax.set_title("Independent validation: 5tp d0-d4 cells vs. this pseudotime")
    plt.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Conclusion: this pseudotime holds up under independent validation

    Full 5-timepoints d0-d4 cells, projected into `rna_clean`'s PCA space (all 50 PCs) and given a `palantir_pseudotime` value via distance-weighted kNN transfer from their nearest neighbors in this dataset, correlate **0.783** with their own known, ground-truth real day labels -- with medians increasing monotonically at every step:

    | 5tp real day | median predicted pseudotime |
    |---|---|
    | 0 | 0.053 |
    | 1 | 0.679 |
    | 2 | 0.792 |
    | 3 | 0.842 |
    | 4 | 0.846 |

    (The large d0->d1 jump followed by progressively smaller increments toward d4 is a typical, expected feature of diffusion-based pseudotime -- more resolution early, saturating later -- not a flaw.)

    Combined with the internal validation on this dataset itself (day-correlation 0.588, both branch terminals landing at their branch's pseudotime maximum, after removing the channel2 d4+ technical failure and the non-EC contamination from day >= 3.5), this pseudotime (`rna_clean.obs["palantir_pseudotime"]`) is solid enough to build on -- a substantial improvement over every earlier attempt this session (DPT: ~0.25 day-correlation; Palantir pre-cleaning: 0.249).

    **Next step (not yet started):** bin cells into ~5,000-cell pseudotime windows (per branch) for pseudobulk analysis, so groups reflect real differentiation progress rather than diluting signal by grouping strictly on nominal real timepoint.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Checkpoint: saving progress for tomorrow

    Saving `rna_clean` (with `palantir_pseudotime`, `branch`, `day`, and the EC-identity scores all attached) so the next session can reload directly instead of re-running the full pipeline (moscot solve, diffusion maps, Palantir) from a cold kernel restart.
    """)
    return


@app.cell
def save_checkpoint(project_root, rna_clean):
    _checkpoint_dir = project_root / "results/10X_multiome_15_timepoints/merged"
    _checkpoint_dir.mkdir(parents=True, exist_ok=True)
    _checkpoint_path = _checkpoint_dir / "rna_clean_with_pseudotime.h5ad"

    # moscot_results (sankey/cell_transition plotting data, a nested dict of
    # DataFrames) was inherited via the .copy() slice from rna_merged and is not
    # relevant to rna_clean -- anndata's h5ad writer chokes on its nested,
    # ragged structure ("inhomogeneous shape"). Drop it before saving.
    if "moscot_results" in rna_clean.uns:
        del rna_clean.uns["moscot_results"]

    rna_clean.write_h5ad(_checkpoint_path)

    # Show only the path relative to the repo, not the full local filesystem path
    _checkpoint_path.relative_to(project_root)
    return


if __name__ == "__main__":
    app.run()
