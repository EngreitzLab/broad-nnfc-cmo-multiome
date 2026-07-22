# CMO Quantification — McGinnis 5 Timepoints

Quantifies Cell Multiplexing Oligo (CMO) tags per cell barcode for the McGinnis 10x multiome 5 timepoints dataset using the [kb-python KITE workflow](https://www.biorxiv.org/content/10.1101/2021.03.11.435036v1), orchestrated via the [`igvf-kite-cmo`](https://github.com/IGVF/atomic-workflows/tree/main/modules/igvf-kite-cmo) atomic workflow.

IGVF analysis set: [IGVFDS1612ZNCA](https://data.igvf.org/analysis-sets/IGVFDS1612ZNCA/)

## Overview

Single channel, three-read CMO library (R1 = cell barcode + UMI, R2 = CMO tag, R3 = cDNA). All inputs are downloaded from IGVF — no local data files required.

1. **CMO barcodes** — downloads `IGVFFI5955PKRW` and converts to KITE TSV (`barcode` → sequence, `sample description` → name)
2. **KITE index** — builds the kallisto index from the CMO barcodes
3. **Barcode onlist** — downloads `IGVFFI8751YQRY` (10x multiome cell barcode whitelist) and decompresses
4. **Seqspec** — downloads `IGVFFI6462KGWB` from IGVF and runs `seqspec index -m tag -s file -t kb` to derive the kb read format
5. **FASTQs** — downloads R1, R2, R3 from auxiliary set `IGVFDS4882ESVG`
6. **Quantification** — runs `kb count` with R1, R2, R3 in order

## IGVF accessions

| Resource | Accession | Description |
|---|---|---|
| CMO auxiliary set | `IGVFDS4882ESVG` | R1 (barcode+UMI), R2 (CMO tag), R3 (cDNA), I1 |
| CMO barcodes | `IGVFFI5955PKRW` | Barcode-to-sample map (`barcode`, `sample description` columns) |
| Cell barcode onlist | `IGVFFI8751YQRY` | 10x multiome barcode whitelist |
| CMO seqspec | `IGVFFI6462KGWB` | seqspec YAML for the CMO library |

## Usage

### Slurm (recommended)

```bash
sbatch analyses/10x_multi_5_timepoints_mcginnis/cmo_quantification/quantify_cmo_tags.sbatch
```

Resources: 8 CPUs, 32 GB RAM, 8 h walltime on `normal` partition.

### Interactive (testing only)

```bash
source /home/users/emattei/miniforge3/etc/profile.d/conda.sh
conda activate /oak/stanford/groups/engreitz/Users/emattei/git/broad-nnfc-cmo-multiome/env/nnfc-cmo-multiome

python analyses/10x_multi_5_timepoints_mcginnis/cmo_quantification/quantify_cmo_tags.py \
    --work-dir /tmp/cmo_quant_mcginnis_test
```

## Outputs

```
/oak/stanford/groups/engreitz/Projects/EC_Screen/Data/10x_5_timepoints_mcginnis/CMO_counts/
└── cmo/
    ├── counts_unfiltered/
    │   ├── adata.h5ad
    │   ├── cells_x_features.barcodes.txt
    │   └── cells_x_features.genes.txt
    └── ...
```

## Authentication

IGVF downloads require credentials in a `.env` file at the project root:

```
IGVF_API_KEY=<your_key>
IGVF_SECRET_KEY=<your_secret>
```

## Dependencies

- [`igvf-kite-cmo`](https://github.com/IGVF/atomic-workflows/tree/main/modules/igvf-kite-cmo) atomic workflow
- `kb-python`, `seqspec`, `igvf-utils` — available in the project conda env
- Helper functions in `scripts/python/utils.py`
