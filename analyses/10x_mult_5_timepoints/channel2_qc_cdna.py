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

__generated_with = "0.23.5"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import anndata as ad
    from pathlib import Path
    from scipy.io import mmread
    import matplotlib.pyplot as plt


    ec_diff_palette = {
      "d0": "#C6C7C7",
      "d1": "#A8B1D6",
      "d2": "#EBBC9E",
      "d3": "#FBC1C3",
      "d4": "#F7999C",
      "Unassigned": "#BEBEBE",
    }
    return Path, ad, mmread, np, pd, plt, sc


@app.cell
def _(pd):
    old_adata_cmo_assignments = pd.read_csv(
        "/Users/emattei/GitHub/broad-nnfc-cmo-multiome/data/barcode_mapping_multiome_dulguun.dedup.tsv",
        sep="\t",
        header=None,
        names=["barcode", "annotation"]
    )

    # keep only channel 2 explicitly
    old_adata_cmo_assignments_ch2 = old_adata_cmo_assignments[
        old_adata_cmo_assignments["barcode"].str.endswith("-1_2")
    ].copy()

    # build cleaned barcode
    old_adata_cmo_assignments_ch2["cleaned_barcode"] = (
        old_adata_cmo_assignments_ch2["barcode"]
        .str.replace(r"-1_2$", "", regex=True)
        + "_10x_5timepoints_channel2"
    )

    old_adata_cmo_assignments_ch2.reset_index(drop=True, inplace=True)
    old_assignments_df = (
        old_adata_cmo_assignments_ch2
        .loc[:, ["cleaned_barcode", "annotation"]]
        .rename(columns={
            "cleaned_barcode": "cell_barcode",
            "annotation": "old_annotation"
        })
        .set_index("cell_barcode")
    )
    old_assignments_df["old_annotation"].value_counts()
    return (old_assignments_df,)


@app.cell
def _(pd):
    gene_metadata_path = "../annotations/gencode.v43.protein_coding.TSS500bp.bed"
    gene_metadata_df = pd.read_csv(
        gene_metadata_path,
        sep="\t",
        header=0,
        names=["chr", "start", "end", "gene_symbol", "score", "strand", "gene_id", "gene_type"]
    )
    gene_metadata_df
    return (gene_metadata_df,)


@app.cell
def _(ad, gene_metadata_df, old_assignments_df, sc):
    h5ad_path = "../data/h5ad/10x_5timepoints_channel2.h5ad"
    adata = ad.read_h5ad(h5ad_path)

    # -------------------------
    # Add gene metadata to var
    # -------------------------
    original_var_names = adata.var_names.copy()

    adata.var["gene_id_base"] = original_var_names.str.replace(r"\.\d+$", "", regex=True)

    adata.var = adata.var.merge(
        gene_metadata_df[
            ["gene_id", "gene_symbol", "gene_type", "chr", "start", "end", "strand"]
        ],
        left_on="gene_id_base",
        right_on="gene_id",
        how="left",
        sort=False
    )

    adata.var.index = original_var_names
    adata.var.index.name = None

    adata.obs = adata.obs.merge(
        old_assignments_df,
        left_index=True,
        right_index=True,
        how="left"
    )
    adata.obs["old_annotation"] = adata.obs["old_annotation"].fillna("Unassigned")

    # -------------------------
    # QC gene flags
    # -------------------------
    gene_symbol = adata.var["gene_symbol"].fillna("")

    adata.var["mt"] = gene_symbol.str.startswith("MT-")
    adata.var["ribo"] = gene_symbol.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = gene_symbol.str.contains(r"^HB(?!P)", regex=True)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

    adata
    return (adata,)


@app.cell
def _(adata, sc):
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
        jitter=0.4,
        multi_panel=True,
    )
    return


@app.cell
def _(adata, np, plt):

    # Plot knee plot for n_genes_by_counts
    # total counts per cell
    counts = adata.obs.loc[adata.obs["old_annotation"] != "Unassigned", "total_counts"].values

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
def _(adata, sc):
    sc.pl.scatter(adata, "total_counts", "pct_counts_mt", color="pct_counts_mt")
    sc.pl.scatter(adata, "total_counts", "pct_counts_ribo", color="pct_counts_ribo")
    sc.pl.scatter(adata, "total_counts", "pct_counts_hb", color="pct_counts_hb")
    return


@app.cell
def _(adata, sc):
    # Compute quantile total counts
    #total_counts_95th_percentile = adata.obs["total_counts"].quantile(0.95)
    sc.pp.filter_cells(adata, min_counts=1000)
    sc.pp.filter_cells(adata, max_counts=5000)
    sc.pp.filter_cells(adata, min_genes=500)
    sc.pp.filter_genes(adata, min_cells=3)
    return


@app.cell
def _(adata, plt, sc):
    ax = sc.pl.scatter(
        adata,
        x="total_counts",
        y="pct_counts_mt",
        color="pct_counts_mt",
        show=False  # important
    )

    # horizontal line (e.g. mitochondrial cutoff)
    ax.axhline(y=8, color="red", linestyle="--")

    # vertical line (e.g. UMI cutoff)
    ax.axvline(x=1000, color="blue", linestyle="--")

    plt.show()
    return


@app.cell
def _(adata):
    adata
    return


@app.cell
def _(adata, sc):
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
        jitter=0.4,
        multi_panel=True,
    )
    return


@app.cell
def _(adata):
    raw = adata.copy()
    return (raw,)


@app.cell
def _(adata, sc):
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)
    return


@app.cell
def _(adata, sc):
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    return


@app.cell
def _(adata, sc):
    sc.pl.highly_variable_genes(adata)
    return


@app.cell
def _(adata, sc):
    sc.tl.pca(adata)
    return


@app.cell
def _(adata, sc):
    sc.pl.pca_variance_ratio(adata, n_pcs=15, log=True)
    return


@app.cell
def _(adata):
    adata.obs["old_annotation"] = adata.obs["old_annotation"].astype("category")
    return


@app.cell
def _(adata):
    cats = adata.obs["old_annotation"].cat.categories

    dpalette = {
      "d0": "#C6C7C7",
      "d1": "#A8B1D6",
      "d2": "#EBBC9E",
      "d3": "#FBC1C3",
      "d4": "#F7999C",
      "Unassigned": "#000000",
    }

    adata.uns["old_annotation_colors"] = [
        dpalette[c] for c in cats
    ]
    return


@app.cell
def _(adata, sc):
    sc.pl.pca(
        adata,
        color=["old_annotation", "pct_counts_mt", "pct_counts_mt"],
        dimensions=[(0, 1), (0, 1), (0, 1)],
        ncols=2,
        size=2,
    )
    return


@app.cell
def _(adata, sc):
    sc.pl.pca(
        adata,
        color=["old_annotation", "old_annotation", "old_annotation", "old_annotation"],
        dimensions=[(0, 1), (2, 3), (4, 5), (5,6)],
        ncols=2,
        size=2,
    )
    return


@app.cell
def _(adata, sc):
    sc.pp.neighbors(adata)
    return


@app.cell
def _(adata, sc):
    sc.tl.umap(adata)
    return


@app.cell
def _(adata):
    sum(adata.obs["old_annotation"] != "Unassigned")
    return


@app.cell
def _(adata, sc):
    sc.pl.umap(
        adata[adata.obs["old_annotation"] != "Unassigned"],
        color="old_annotation",
        # Setting a smaller point size to get prevent overlap
        #size=10,
    )
    return


@app.cell
def _(adata, sc):
    sc.pl.umap(
        adata,
        color="old_annotation",
        # Setting a smaller point size to get prevent overlap
        size=10,
    )
    return


@app.cell
def _(adata, sc):
    # Using the igraph implementation and a fixed number of iterations can be significantly faster,
    # especially for larger datasets
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
    return


@app.cell
def _(adata, sc):
    sc.pl.umap(adata, color=["leiden"])
    return


@app.cell
def _(adata):
    adata
    return


@app.cell
def _(adata, sc):
    sc.pl.umap(
        adata,
        color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
        wspace=0.5,
        ncols=2,
    )
    return


@app.cell
def _(adata, pd):
    pd.crosstab(adata.obs["leiden"], adata.obs["old_annotation"])
    return


@app.cell
def _(raw, sc):
    sc.pp.scrublet(
        raw,
        expected_doublet_rate=0.4,
        threshold=0.6
    )
    return


@app.cell
def _(raw, sc):
    sc.pl.scrublet_score_distribution(
        raw)
    return


@app.cell
def _(raw):
    sum(raw.obs["doublet_score"] < 0.35)
    return


@app.cell
def _(adata):
    adata
    return


@app.cell
def _(adata, raw, sc):
    adata.obs["doublet_score"] = raw.obs["doublet_score"]
    sc.pl.umap(
        adata[adata.obs["doublet_score"] < 0.35],
        color=["doublet_score"],
        wspace=0.5,
        ncols=2,
    )
    return


@app.cell
def _(adata, pd):
    pd.crosstab(adata.obs["leiden"], adata.obs["doublet_score"] < 0.35)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## CMO Hash classification

    Run the adata_cmo classification mimicking Seurat::HTOdemux.
    Only the cells passing scanpy filtering are included in the classification. This is because the adata_cmo classification is used to filter out low-quality cells, and including low-quality cells in the classification may lead to inaccurate results. By filtering out low-quality cells before running the adata_cmo classification, we can ensure that the classification is performed on high-quality cells, which can lead to more accurate results and improved computational efficiency.
    """)
    return


@app.cell
def _(Path, mmread):
    cmo_counts_path = Path("../data/cmo_counts/channel2/counts_unfiltered")

    mat = mmread(cmo_counts_path / "cells_x_features.mtx").tocsr()

    barcodes = [
        bc + "_10x_5timepoints_channel2"
        for bc in (cmo_counts_path / "cells_x_features.barcodes.txt").read_text().splitlines()
    ]
    genes = (cmo_counts_path / "cells_x_features.genes.txt").read_text().splitlines()
    gene_names = (cmo_counts_path / "cells_x_features.genes.names.txt").read_text().splitlines()
    return barcodes, gene_names, genes, mat


@app.cell
def _(Path, adata):
    Path("../results/channel2/filtered_barcodes.txt").write_text("\n".join(adata.obs_names[adata.obs["doublet_score"] < 0.35].to_list()))
    return


@app.cell
def _(ad, adata, barcodes, gene_names, genes, mat, pd):
    adata_cmo = ad.AnnData(
        X=mat,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(
            {"gene_name": gene_names},
            index=genes,
        ),
    )
    # Filter adata_cmo data to cells in adata
    cells_in_use = adata.obs_names.intersection(adata_cmo.obs_names.to_list())[adata.obs["doublet_score"] < 0.35]
    adata_cmo = adata_cmo[cells_in_use, :]
    adata_flt = adata[cells_in_use, :]
    adata_cmo
    return adata_cmo, adata_flt, cells_in_use


@app.cell
def _(adata_cmo, np):
    cmo = adata_cmo.X
    cmo = cmo.toarray() if hasattr(cmo, "toarray") else np.asarray(cmo)

    # CLR normalization per cell, like Seurat margin=2
    gm = np.exp(np.mean(np.log1p(cmo), axis=1))
    clr = np.log1p(cmo / gm[:, None])
    return (clr,)


@app.cell
def _(clr, np):
    # threshold per adata_cmo
    positive_quantile = 0.97
    thresholds = np.quantile(clr, positive_quantile, axis=0)
    return positive_quantile, thresholds


@app.cell
def _(clr, thresholds):
    #positive = clr > thresholds
    positive = clr > thresholds
    n_pos = positive.sum(axis=1)
    return n_pos, positive


@app.cell
def _(plt, positive_quantile, thresholds):
    # Plot histogram of thresholds
    plt.figure(figsize=(5,5))
    plt.hist(thresholds, bins=50)
    plt.xlabel("CLR threshold")
    plt.ylabel("Number of features")
    plt.title(f"CLR thresholds at {positive_quantile} quantile")
    plt.show()
    return


@app.cell
def _(adata_cmo, adata_flt, clr, n_pos, np, pd, positive):
    # assigned adata_cmo = strongest adata_cmo among positives, or max adata_cmo if doublet
    top_idx = np.argmax(clr, axis=1)
    top_adata_cmo = adata_cmo.var["gene_name"].values[top_idx]

    hash_ID = np.where(n_pos == 0, "Negative", top_adata_cmo)

    global_class = np.where(
        n_pos == 0,
        "Negative",
        np.where(n_pos == 1, "Singlet", "Doublet")
    )

    adata_flt.obs["hash_ID_scanpy"] = hash_ID
    adata_flt.obs["cmo_classification_global_scanpy"] = global_class

    # optional: exact positive adata_cmo list per cell
    adata_flt.obs["cmo_positive_tags_scanpy"] = [
        ",".join(adata_cmo.var["gene_name"].values[row]) if row.any() else "Negative"
        for row in positive
    ]

    cmo_to_timepoint = {
        f"CMO{i:02d}": f"d{(i - 1) // 5}"
        for i in range(1, 26)
    }

    adata_flt.obs["timepoint_scanpy"] = adata_flt.obs["cmo_positive_tags_scanpy"].map(cmo_to_timepoint)
    adata_flt.obs["timepoint_scanpy"] = adata_flt.obs["timepoint_scanpy"].fillna(
        adata_flt.obs["cmo_classification_global_scanpy"]
    )

    pd.crosstab(
        adata_flt.obs["cmo_classification_global_scanpy"],
        adata_flt.obs["cmo_positive_tags_scanpy"]
    )
    return (cmo_to_timepoint,)


@app.cell
def _(adata_flt, cmo_to_timepoint, pd):
    adata_flt.obs["timepoint_scanpy"] = adata_flt.obs["cmo_positive_tags_scanpy"].map(cmo_to_timepoint)
    adata_flt.obs["timepoint_scanpy"] = adata_flt.obs["timepoint_scanpy"].fillna(
        adata_flt.obs["cmo_classification_global_scanpy"]
    )
    adata_flt.obs["timepoint_scanpy"] = pd.Categorical(
        adata_flt.obs["timepoint_scanpy"],
        categories=["d0", "d1", "d2", "d3", "d4", "Doublet", "Negative"],
        ordered=True
    )
    return


@app.cell
def _(adata_flt):
    sum(adata_flt.obs["timepoint_scanpy"]=="Negative")
    return


@app.cell
def _(adata_flt):
    sum((adata_flt.obs["timepoint_scanpy"] != "Doublet") & (adata_flt.obs["timepoint_scanpy"] != "Negative"))
    return


@app.cell
def _(adata_flt, sc):
    sc.pl.umap(
        adata_flt[(adata_flt.obs["timepoint_scanpy"] != "Doublet") & (adata_flt.obs["timepoint_scanpy"] != "Negative")],
        color=["timepoint_scanpy"],
        wspace=0.5,
        ncols=2,
    )
    return


@app.cell
def _(adata_flt, pd):
    pd.crosstab(adata_flt.obs["timepoint_scanpy"], adata_flt.obs["old_annotation"])
    return


@app.cell
def _(adata_flt, pd):
    pd.crosstab(adata_flt.obs["leiden"], adata_flt.obs["timepoint_scanpy"])
    return


@app.cell
def _(adata_flt):
    sum((adata_flt.obs["timepoint_scanpy"] != "Doublet") & (adata_flt.obs["timepoint_scanpy"] != "Negative"))
    return


@app.cell
def _(adata_flt):
    adata_flt.obs["old_annotation"].value_counts()
    return


@app.cell
def _(adata_flt, sc):
    # adata_flt plot boxplot of doublet_scores by cmo_classification_global_scanpy
    sc.pl.violin(adata_flt[adata_flt.obs["cmo_classification_global_scanpy"] == "Singlet"], keys="doublet_score", groupby="timepoint_scanpy")
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
    import altair as alt

    return alt, snap


@app.cell
def _(pd):
    chrom_dict = pd.read_csv("../annotations/GRCh38_EBV.chrom.sizes.no.alt.tsv", sep="\t", header=None, names=["chr", "size"])
    chrom_dict = chrom_dict.set_index("chr")["size"].to_dict()
    list(chrom_dict.items())[1:5]
    return (chrom_dict,)


@app.cell
def _(cells_in_use, chrom_dict, snap):
    atac_adata = snap.pp.import_fragments(
        "../data/fragments/10x_5timepoints_channel2.fragments.tsv.gz",
        sorted_by_barcode=False,
        chrom_sizes=chrom_dict,
        whitelist=cells_in_use,
        min_num_fragments=1000,
    )
    #mask = atac_adata.obs_names.isin(cells_in_use)
    #atac_adata = atac_adata[mask, :].copy()
    snap.metrics.tsse(atac_adata, "../annotations/gencode.v43.primary_assembly.annotation.gtf.gz")
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
        (plot_df["n_fragment"] <= 30_000)
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
def _(atac_adata, snap):
    snap.pp.filter_cells(atac_adata, min_counts=1000, min_tsse=7.5, max_counts=20_000)
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


@app.cell
def _(atac_adata, snap):
    snap.pl.umap(atac_adata, color='leiden', interactive=False, height=500)
    return


@app.cell
def _(adata_flt, atac_adata, cells_in_use):
    common_cells = atac_adata.obs_names.intersection(cells_in_use)

    atac_adata.obs["timepoint_scanpy"] = adata_flt.obs.loc[common_cells, "timepoint_scanpy"].astype("category")
    return


@app.cell
def _(atac_adata):
    timepoint_palette = {
        "d0": "#C6C7C7",
        "d1": "#A8B1D6",
        "d2": "#EBBC9E",
        "d3": "#FBC1C3",
        "d4": "#F7999C",
        "Doublet": "#7F7F7F",
        "Negative": "#BEBEBE",
        "Unassigned": "#000000",
    }
    cats_atac = atac_adata.obs["timepoint_scanpy"].cat.categories
    atac_adata.obs["timepoint_color"] = (
        atac_adata.obs["timepoint_scanpy"]
        .map(timepoint_palette)
    )
    return (timepoint_palette,)


@app.cell
def _(atac_adata):
    print(atac_adata.obs["timepoint_scanpy"].cat.categories.tolist())
    return


@app.cell
def _(atac_adata, snap):
    snap.pl.umap(
        atac_adata,
        color="timepoint_scanpy",
        interactive=False,
        height=500
    )
    return


@app.cell
def _(atac_adata, snap, timepoint_palette):

    def snap_umap_with_palette(
        adata,
        color,
        palette,
        interactive=False,
        height=500,
        width=None,
        **kwargs,
    ):
        fig = snap.pl.umap(
            adata,
            color=color,
            interactive=interactive,
            height=height,
            width=width,
            show=False,
            **kwargs,
        )

        for trace in fig.data:
            label = str(trace.name)
            if label in palette:
                trace.marker.color = palette[label]

        fig.show()
        return fig

    fig = snap_umap_with_palette(
        atac_adata,
        color="timepoint_scanpy",
        palette=timepoint_palette,
        interactive=False,
        height=500,
    )
    for trace in fig.data:
        label = str(trace.name)
        trace.marker.color = timepoint_palette.get(label, trace.marker.color)
    return


@app.cell
def _(atac_adata):
    atac_adata.obs["timepoint_scanpy"].value_counts()
    return


if __name__ == "__main__":
    app.run()
