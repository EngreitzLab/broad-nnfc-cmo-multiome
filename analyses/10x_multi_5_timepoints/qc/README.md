# QC

Quality control and demultiplexing of the 10x multiome 5 timepoints data (RNA + ATAC).

## Notebooks

| Notebook | Description |
|---|---|
| `channel1_qc.py` | QC filtering, CMO hash demultiplexing, and ATAC processing for channel 1 |
| `channel2_qc.py` | QC filtering, CMO hash demultiplexing, and ATAC processing for channel 2 |

Both are [marimo](https://marimo.io) notebooks with inline dependencies (PEP 723).

## Usage

```bash
marimo edit analyses/10x_multi_5_timepoints/qc/channel1_qc.py
marimo edit analyses/10x_multi_5_timepoints/qc/channel2_qc.py
```

## Inputs

- GEX counts: TBD
- CMO counts: `/oak/stanford/groups/engreitz/Projects/EC_Screen/Data/10x_5_timepoints/CMO_counts/channel{1,2}/counts_unfiltered/adata.h5ad`
- ATAC fragments: TBD

## Outputs

TBD
