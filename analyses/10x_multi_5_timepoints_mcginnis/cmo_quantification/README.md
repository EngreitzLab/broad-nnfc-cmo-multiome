# CMO Quantification — McGinnis 5 Timepoints

Quantifies Cell Multiplexing Oligo (CMO) tags per cell barcode for the McGinnis 10x multiome 5 timepoints dataset using the [kb-python KITE workflow](https://www.biorxiv.org/content/10.1101/2021.03.11.435036v1), orchestrated via the [`igvf-kite-cmo`](https://github.com/IGVF/atomic-workflows/tree/main/modules/igvf-kite-cmo) atomic workflow.

IGVF analysis set: [IGVFDS1612ZNCA](https://data.igvf.org/analysis-sets/IGVFDS1612ZNCA/)

## Overview

Single channel. The CMO library has three reads (R1 = cell barcode + UMI, R2 = MULTI-seq CMO tag, R3 = cDNA) but only R1 and R2 are needed for CMO quantification. All inputs are downloaded from IGVF — no local data files required.

1. **CMO barcodes** — downloads `IGVFFI5955PKRW` and converts to KITE TSV (`barcode` → sequence, `sample description` → name)
2. **KITE index** — builds the kallisto index from the CMO barcodes
3. **Barcode onlist** — downloads `IGVFFI8751YQRY` (10x multiome cell barcode whitelist) and decompresses
4. **FASTQs** — downloads R1 and R2 from auxiliary set `IGVFDS4882ESVG`
5. **Quantification** — runs `kb count` with R1 and R2 using hardcoded read format `0,0,16:0,16,28:1,0,8`

## Read format

The kb read format is hardcoded as `0,0,16:0,16,28:1,0,8` (file0=R1: cell barcode bp 0–16, file0=R1: UMI bp 16–28, file1=R2: CMO tag bp 0–8). This was derived manually by inspecting the seqspec YAML (`IGVFFI6462KGWB`).

The `seqspec index -m tag -s file -t kb` tool cannot be used for this library because it produces a syntactically invalid format string — it emits six comma-separated values in the first segment instead of the required three, resulting in a string like `0,0,16,2,0,8:0,16,28:` that `kb count` rejects. The root cause is that the seqspec YAML for this library describes the CMO tag region in a way the tool does not handle correctly. The hardcoded format was validated against the seqspec source and confirmed correct.

## IGVF accessions

| Resource | Accession | Description |
|---|---|---|
| CMO auxiliary set | `IGVFDS4882ESVG` | R1 (barcode+UMI), R2 (CMO tag), R3 (cDNA, unused), I1 (unused) |
| CMO barcodes | `IGVFFI5955PKRW` | Barcode-to-sample map (`barcode`, `sample description` columns) |
| Cell barcode onlist | `IGVFFI8751YQRY` | 10x multiome barcode whitelist |
| CMO seqspec | `IGVFFI6462KGWB` | seqspec YAML (reference only — read format is hardcoded, see above) |

## Usage

### Slurm (recommended)

```bash
sbatch analyses/10x_multi_5_timepoints_mcginnis/cmo_quantification/quantify_cmo_tags.sbatch
```

Resources: 8 CPUs, 32 GB RAM, 8 h walltime on `engreitz` partition.

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
