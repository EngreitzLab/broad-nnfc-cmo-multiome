# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "altair==6.2.2",
#     "anndata==0.13.2",
#     "marimo>=0.23.3",
#     "matplotlib==3.11.1",
#     "numpy==2.5.1",
#     "pandas==2.3.3",
#     "python-dotenv==1.2.2",
#     "seaborn==0.13.2",
#     "snapatac2==2.9.0",
#     "vegafusion==2.0.3",
#     "vl-convert-python==1.9.0.post1",
# ]
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
    # Merging channel1 and channel2

    This notebook loads the strict-QC-filtered RNA and ATAC AnnData objects produced by `channel1_qc.py` and `channel2_qc.py`, concatenates them into a single merged RNA object and a single merged ATAC object, plots summary statistics to sanity-check the two channels are comparable, and then splits the merged data by timepoint (`rescued_cmo_tag`), writing one RNA h5ad and one ATAC fragments file per timepoint (d0-d4) to `results/10X_multiome_5_timepoints/merged/`.
    """)
    return


@app.cell
def imports():
    import altair as alt
    import anndata as ad
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import snapatac2 as snap
    from dotenv import find_dotenv, load_dotenv
    from pathlib import Path
    alt.data_transformers.enable("vegafusion")

    # get the project root path using the '.env' file.
    _env_path = find_dotenv(usecwd=True)
    project_root = Path(_env_path).parent
    return ad, alt, pd, plt, project_root, snap, sns


@app.cell
def _():
    # Endothelial-differentiation timepoint palette, and Okabe-Ito colorblind-safe
    # palette (channel1 / channel2 borrow two of its colors), matching the
    # conventions used in channel1_qc.py / channel2_qc.py.
    ec_diff_palette = {
        "d0": "#C6C7C7",
        "d1": "#A8B1D6",
        "d2": "#EBBC9E",
        "d3": "#FBC1C3",
        "d4": "#F7999C",
        "Unassigned": "#4D4D4D",
    }

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
    return (channel_palette,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Loading the per-channel RNA and ATAC results
    """)
    return


@app.cell
def load_rna_channels(ad, project_root):
    rna_channel1 = ad.read_h5ad(project_root / "results/10X_multiome_5_timepoints/channel1/rna_adata.h5ad")
    rna_channel2 = ad.read_h5ad(project_root / "results/10X_multiome_5_timepoints/channel2/rna_adata.h5ad")

    [rna_channel1.n_obs, rna_channel2.n_obs]
    return rna_channel1, rna_channel2


@app.cell
def load_atac_channels(project_root, snap):
    # Read with snap.read (not a generic h5ad reader) so the snapatac2-specific
    # internals (obsm["fragment_paired"], uns["reference_sequences"], etc.) come
    # back as the native snapatac2 AnnData type, not a plain anndata.AnnData.
    atac_channel1 = snap.read(project_root / "results/10X_multiome_5_timepoints/channel1/atac_adata.h5ad")
    atac_channel2 = snap.read(project_root / "results/10X_multiome_5_timepoints/channel2/atac_adata.h5ad")

    # channel1 and channel2 have "leiden" and "rescued_cmo_tag" in swapped obs
    # column order. snapatac2's Rust-side vstack (used by both snap.concat and
    # AnnDataSet) requires identical column *order* across components, not just
    # identical column names, so mismatched order otherwise fails with
    # "unable to vstack, column names don't match". Harmonize channel2 to
    # channel1's column order. Per-column deletion (del atac_channel2.obs[col])
    # isn't implemented in this snapatac2 version, but reassigning the whole obs
    # to a reordered DataFrame at once works.
    _col_order = atac_channel1.obs[:].columns
    atac_channel2.obs = atac_channel2.obs[:].select(_col_order)

    assert atac_channel1.obs[:].columns == atac_channel2.obs[:].columns

    [atac_channel1.n_obs, atac_channel2.n_obs]
    return atac_channel1, atac_channel2


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Concatenating the two channels
    """)
    return


@app.cell
def concat_rna(ad, rna_channel1, rna_channel2):
    # Barcodes already carry a unique per-sample suffix (e.g. _IGVFSM6578YOFR
    # vs. _IGVFSM7208ZGNY), and both channels share an identical gene set (same
    # GENCODE v43 annotation), so a plain inner-join concat is safe here.
    rna_merged = ad.concat(
        {"channel1": rna_channel1, "channel2": rna_channel2},
        label="channel",
        join="inner",
        # var columns (gene_symbol, gene_type, gene_id, ...) are identical
        # across channels (same GENCODE v43 annotation) but ad.concat drops
        # all non-var_names var columns by default -- merge="same" keeps any
        # column that matches exactly across the concatenated objects, which
        # save_rna_per_timepoint below relies on (gene_type, gene_symbol).
        merge="same",
    )
    rna_merged
    return (rna_merged,)


@app.cell
def concat_atac(atac_channel1, atac_channel2, project_root, snap):
    # AnnDataSet (not snap.concat / generic anndata.concat): a lazy reference to
    # the two component files rather than a physical vstack, so it only requires
    # matching var_names (confirmed identical between channels) and sidesteps the
    # obs-column-*order* mismatch that broke snap.concat here (channel1 has
    # rescued_cmo_tag before leiden, channel2 has them swapped). It also natively
    # supports snap.ex.export_fragments and split_obs_by for the per-timepoint
    # outputs below.
    _dataset_path = project_root / "results/10X_multiome_5_timepoints/merged/atac_dataset.h5ads"
    _dataset_path.parent.mkdir(parents=True, exist_ok=True)
    # H5Fcreate can't truncate a path that's still open from a prior run of this
    # cell (stale fd in this same kernel process); unlinking first is safe on
    # POSIX since it only drops the directory entry, not any handle already open
    # on the old inode.
    _dataset_path.unlink(missing_ok=True)

    atac_merged = snap.AnnDataSet(
        adatas=[("channel1", atac_channel1), ("channel2", atac_channel2)],
        filename=_dataset_path,
        add_key="channel",
    )
    atac_merged
    return (atac_merged,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Summary statistics: are the two channels comparable?
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Cell counts per channel x timepoint
    """)
    return


@app.cell
def composition_plot(alt, channel_palette, rna_merged):
    _timepoint_order = ["d0", "d1", "d2", "d3", "d4"]
    _comp_df = (
        rna_merged.obs.groupby(["channel", "rescued_cmo_tag"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )

    alt.Chart(_comp_df).mark_bar().encode(
        x=alt.X("rescued_cmo_tag:N", title="Timepoint", sort=_timepoint_order),
        y=alt.Y("n_cells:Q", title="Number of cells"),
        xOffset=alt.XOffset("channel:N"),
        color=alt.Color(
            "channel:N",
            scale=alt.Scale(domain=list(channel_palette), range=list(channel_palette.values())),
        ),
        tooltip=["channel", "rescued_cmo_tag", "n_cells"],
    ).properties(
        title="Cell counts per channel x timepoint (RNA, pass_strict_qc)",
        width=500, height=350,
    ).configure_view(strokeWidth=0)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### RNA QC metrics by channel
    """)
    return


@app.cell
def rna_qc_by_channel(channel_palette, plt, rna_merged, sns):
    _rna_qc_metrics = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    _rna_qc_long = rna_merged.obs[["channel"] + _rna_qc_metrics].melt(
        id_vars="channel", var_name="metric", value_name="value"
    )

    _fig, _axes = plt.subplots(1, len(_rna_qc_metrics), figsize=(4 * len(_rna_qc_metrics), 4))
    for _ax, _metric in zip(_axes, _rna_qc_metrics):
        sns.boxplot(
            data=_rna_qc_long[_rna_qc_long["metric"] == _metric],
            x="channel",
            y="value",
            hue="channel",
            palette=channel_palette,
            legend=False,
            ax=_ax,
        )
        _ax.set_title(_metric)
        _ax.set_xlabel(None)
        _ax.set_ylabel(None)
    plt.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### ATAC QC metrics by channel
    """)
    return


@app.cell
def atac_qc_by_channel(atac_merged, channel_palette, plt, sns):
    _atac_qc_metrics = ["n_fragment", "tsse", "frac_mito"]
    _atac_obs = atac_merged.adatas.obs[:].select(_atac_qc_metrics).to_pandas()
    _atac_obs["channel"] = atac_merged.obs[:]["channel"].to_pandas()
    _atac_qc_long = _atac_obs[["channel"] + _atac_qc_metrics].melt(
        id_vars="channel", var_name="metric", value_name="value"
    )

    _fig, _axes = plt.subplots(1, len(_atac_qc_metrics), figsize=(4 * len(_atac_qc_metrics), 4))
    for _ax, _metric in zip(_axes, _atac_qc_metrics):
        sns.boxplot(
            data=_atac_qc_long[_atac_qc_long["metric"] == _metric],
            x="channel",
            y="value",
            hue="channel",
            palette=channel_palette,
            legend=False,
            ax=_ax,
        )
        _ax.set_title(_metric)
        _ax.set_xlabel(None)
        _ax.set_ylabel(None)
    plt.tight_layout()
    _fig
    return


@app.cell
def compute_shared_barcodes(atac_merged, pd, rna_merged):
    # Barcodes passing strict QC in both modalities. ATAC QC (fragment count,
    # TSSE) and RNA QC (mito content, doublet calls) are independent filters, so
    # a strict-QC RNA barcode isn't guaranteed to also be a strict-QC ATAC
    # barcode and vice versa. Both per-timepoint exports below (RNA h5ad and
    # ATAC fragments) are restricted to this intersection.
    shared_barcodes = pd.Index(rna_merged.obs_names).intersection(atac_merged.obs_names)
    len(shared_barcodes)
    return (shared_barcodes,)


@app.cell(hide_code=True)
def cell_counts_intro(mo):
    mo.md(r"""
    ### Cell counts: RNA, ATAC, and overlap per timepoint
    """)
    return


@app.cell
def cell_counts_per_timepoint(atac_merged, pd, rna_merged, shared_barcodes):
    _timepoint_order = ["d0", "d1", "d2", "d3", "d4"]

    _rna_tp = rna_merged.obs["rescued_cmo_tag"]
    _rna_tp.index = rna_merged.obs_names

    _atac_tp = pd.Series(
        atac_merged.adatas.obs[:]["rescued_cmo_tag"].to_numpy(),
        index=atac_merged.obs_names,
    )

    _rna_counts = _rna_tp.value_counts().reindex(_timepoint_order).fillna(0).astype(int)
    _atac_counts = _atac_tp.value_counts().reindex(_timepoint_order).fillna(0).astype(int)

    _both_tp = _rna_tp.loc[shared_barcodes]
    _both_counts = _both_tp.value_counts().reindex(_timepoint_order).fillna(0).astype(int)

    cell_counts_per_timepoint = pd.DataFrame({
        "RNA": _rna_counts,
        "ATAC": _atac_counts,
        "Both (RNA and ATAC)": _both_counts,
    })
    cell_counts_per_timepoint
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Splitting the merged data by timepoint
    """)
    return


@app.cell
def save_rna_per_timepoint(project_root, rna_merged, shared_barcodes):
    _timepoint_order = ["d0", "d1", "d2", "d3", "d4"]
    _rna_outdir = project_root / "results/10X_multiome_5_timepoints/merged/rna"
    _rna_outdir.mkdir(parents=True, exist_ok=True)

    _rna_shared = rna_merged[shared_barcodes].copy()

    # scE2G requires gene symbols as var_names, not Ensembl gene IDs. Restrict
    # to protein_coding genes first: gene_symbol collisions are pervasive
    # across the full annotation (~1,500, mostly rRNA/pseudogene loci sharing
    # a name) but drop to just 30 once limited to protein-coding genes. Most of
    # those remaining 30 are PAR (pseudoautosomal region) genes, annotated
    # twice under an X and a "_PAR_Y" gene_id for the same symbol (e.g. CD99,
    # ASMT) -- drop the _PAR_Y copy rather than suffix it, since it's a
    # duplicate annotation of the same gene, not a distinct locus.
    # var_names_make_unique handles the small remainder that isn't PAR-related
    # (e.g. HERC3, MATR3: two distinct Ensembl IDs sharing a legacy symbol)
    # via -1/-2 suffixes. The original Ensembl ID stays available in
    # var["gene_id"].
    _rna_shared = _rna_shared[:, _rna_shared.var["gene_type"] == "protein_coding"].copy()
    _rna_shared = _rna_shared[:, ~_rna_shared.var["gene_id"].str.endswith("_PAR_Y")].copy()
    # .to_numpy() strips the pandas Series name ("gene_symbol") so the new
    # index doesn't inherit it -- an index named "gene_symbol" alongside the
    # still-present "gene_symbol" column (whose values differ from the index
    # after make_unique's -1/-2 suffixes) fails to write to h5ad.
    _rna_shared.var_names = _rna_shared.var["gene_symbol"].astype(str).to_numpy()
    _rna_shared.var_names_make_unique()

    rna_timepoint_paths = {}
    for _tp in _timepoint_order:
        _tp_mask = (_rna_shared.obs["rescued_cmo_tag"] == _tp).to_numpy()
        _tp_path = _rna_outdir / f"rna_{_tp}.h5ad"
        _rna_shared[_tp_mask].copy().write_h5ad(_tp_path)
        rna_timepoint_paths[_tp] = _tp_path

    rna_timepoint_paths
    return


@app.cell
def save_atac_fragments_per_timepoint(
    atac_merged,
    pd,
    project_root,
    shared_barcodes,
    snap,
):
    _atac_fragments_outdir = project_root / "results/10X_multiome_5_timepoints/merged/atac_fragments"
    _atac_fragments_outdir.mkdir(parents=True, exist_ok=True)

    _timepoint_order = ["d0", "d1", "d2", "d3", "d4"]

    # groupby as a string resolves via atac_merged.obs[groupby], but atac_merged's
    # own obs table only carries the AnnDataSet add_key ("channel"); the per-cell
    # rescued_cmo_tag column lives in the lazy stacked view atac_merged.adatas.obs
    # instead. Cells outside shared_barcodes are relabeled to an excluded group
    # and dropped via `selections`, since export_fragments has no direct
    # per-cell subsetting.
    _is_shared = pd.Index(atac_merged.obs_names).isin(shared_barcodes)
    _group_labels = atac_merged.adatas.obs[:]["rescued_cmo_tag"].to_list()
    _group_labels = [g if keep else "excluded_not_in_rna" for g, keep in zip(_group_labels, _is_shared)]

    atac_fragment_paths = snap.ex.export_fragments(
        atac_merged,
        groupby=_group_labels,
        selections=_timepoint_order,
        out_dir=str(_atac_fragments_outdir),
        suffix=".fragments.tsv.gz",
        compression="gzip",
    )
    atac_fragment_paths
    return (atac_fragment_paths,)


@app.cell
def close_atac_handles(
    atac_channel1,
    atac_channel2,
    atac_fragment_paths,
    atac_merged,
    rna_channel1,
    rna_channel2,
):
    # Close the backed snapatac2 file handles once all downstream work is done,
    # so the underlying h5ad/h5ads files aren't left locked (snap.read defaults
    # to backed="r+", and the AnnDataSet is itself backed by its own file).
    # rna_channel1/rna_channel2 are loaded fully in-memory (no backed=), so
    # there's no file handle to close, but they're redundant once rna_merged
    # exists -- del them too to free the memory.
    atac_fragment_paths  # ran after fragments were exported

    atac_channel1.close()
    atac_channel2.close()
    atac_merged.close()
    del rna_channel1, rna_channel2
    return


if __name__ == "__main__":
    app.run()
