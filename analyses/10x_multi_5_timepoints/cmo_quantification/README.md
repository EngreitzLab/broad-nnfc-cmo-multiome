# CMO Quantification

Quantifies Cell Multiplexing Oligo (CMO) tags per cell barcode for both channels of the 10x multiome 5 timepoints experiment using the [kb-python KITE workflow](https://www.biorxiv.org/content/10.1101/2021.03.11.435036v1), orchestrated via the [`igvf-kite-cmo`](https://github.com/IGVF/atomic-workflows/tree/main/modules/igvf-kite-cmo) atomic workflow.

## Overview

The pipeline runs entirely from IGVF accessions — no local data files are required. It performs the following steps in order:

1. **CMO barcodes** — downloads `IGVFFI9308MNHO` (CMO barcode CSV) and converts it to the two-column TSV format expected by KITE (`run_kite prepare-barcodes`)
2. **KITE index** — builds the kallisto index from the CMO barcodes (`run_kite index`)
3. **Barcode onlist** — downloads `IGVFFI8751YQRY` (10x Genomics multiome cell barcode whitelist)
4. **Read format** — derived per channel from the seqspec YAML in `/oak/stanford/groups/engreitz/Projects/EC_Screen/Data/10x_5_timepoints/seqspecs/` via `seqspec index -m tag -s file -t kb`
5. **Quantification** — downloads channel FASTQs and runs `kb count` for channel 1 and channel 2 independently (`run_kite quantify`)

Downloads and intermediate files (index, FASTQs) are staged to `--work-dir` (use `$L_SCRATCH` on Slurm — wiped after the job). Results are written to `/oak/stanford/groups/engreitz/Projects/EC_Screen/Data/10x_5_timepoints/CMO_counts/`.

## IGVF accessions

| Resource | Accession | Description |
|---|---|---|
| Channel 1 FASTQs | `IGVFDS2186TMQE` | Auxiliary set — CMO R1 + R2 for channel 1 |
| Channel 2 FASTQs | `IGVFDS0370XXMJ` | Auxiliary set — CMO R1 + R2 for channel 2 |
| CMO barcodes | `IGVFFI9308MNHO` | MULTI-seq barcode CSV (`multiseq_bc`, `CMO ID` columns) |
| Cell barcode onlist | `IGVFFI8751YQRY` | 10x Genomics multiome barcode whitelist |

## Usage

### Slurm (recommended)

```bash
sbatch analyses/10x_multi_5_timepoints/cmo_quantification/quantify_cmo_tags.sbatch
```

Resources: 8 CPUs, 32 GB RAM, 8 h walltime on `engreitz` partition.

### Interactive (testing only)

```bash
source /home/users/emattei/miniforge3/etc/profile.d/conda.sh
conda activate /oak/stanford/groups/engreitz/Users/emattei/git/broad-nnfc-cmo-multiome/env/nnfc-cmo-multiome

python analyses/10x_multi_5_timepoints/cmo_quantification/quantify_cmo_tags.py \
    --work-dir /tmp/cmo_quant_test \
    --channels channel1
```

### Options

```
--ch1-accession     IGVF auxiliary set for channel 1  (default: IGVFDS2186TMQE)
--ch2-accession     IGVF auxiliary set for channel 2  (default: IGVFDS0370XXMJ)
--cmo-accession     IGVF tabular file for CMO barcodes (default: IGVFFI9308MNHO)
--onlist-accession  IGVF tabular file for barcode whitelist (default: IGVFFI8751YQRY)
--seqspec-dir       Directory with per-channel seqspec YAML.gz files (default: /oak/.../seqspecs/)
--channels          Subset of channels to run: channel1, channel2 (default: both)
--work-dir          Staging directory for downloads and index (required)
--output-dir        Results directory (default: results/10x_multi_5_timepoints/cmo_quantification/)
--threads           Threads for kb count (default: 8)
--memory            Memory for bustools sort (default: 4G)
```

## Outputs

```
/oak/stanford/groups/engreitz/Projects/EC_Screen/Data/10x_5_timepoints/CMO_counts/
├── channel1/
│   ├── counts_unfiltered/
│   │   ├── adata.h5ad          <- unfiltered CMO counts (AnnData)
│   │   ├── cells_x_genes.barcodes.txt
│   │   └── cells_x_genes.genes.txt
│   └── ...                     <- full kb count output
└── channel2/
    └── ...
```

## Authentication

IGVF downloads require credentials in a `.env` file at the project root:

```
IGVF_API_KEY=<your_key>
IGVF_SECRET_KEY=<your_secret>
```

## Dependencies

- [`igvf-kite-cmo`](https://github.com/IGVF/atomic-workflows/tree/main/modules/igvf-kite-cmo) atomic workflow (sibling repo at `../atomic-workflows/`)
- `kb-python==0.29.1`, `seqspec==0.3.0`, `igvf-utils` — available in the project conda env (`env/nnfc-cmo-multiome/`)
- Helper functions in `scripts/python/utils.py`
