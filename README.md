# Endothelial cell differentiation — CMO multiome processing

Processing pipelines for the 10x multiome EC differentiation datasets, covering FASTQ reconstruction, CMO/GEX demultiplexing, and CMO tag quantification.

## Datasets

### [10x multiome — 5 timepoints](analyses/10x_multi_5_timepoints/README.md)

GEX and CMO reads arrive mixed in the same FASTQs and must be separated with `splitcode` before quantification. CMO quantification is automated end-to-end via IGVF accessions.

| Channel | IGVF analysis set |
|---|---|
| Channel 1 | [IGVFDS5477BPOI](https://data.igvf.org/analysis-sets/IGVFDS5477BPOI) |
| Channel 2 | [IGVFDS3995WHFT](https://data.igvf.org/analysis-sets/IGVFDS3995WHFT) |

### [10x multiome — 15 timepoints](analyses/10x_multi_15_timepoints/README.md)

GEX, CMO, and ATAC reads arrive pre-demultiplexed (separate sample indices). No `splitcode` step required.

## Repository layout

```
analyses/
  10x_multi_5_timepoints/      # 5 tp pipeline (manual steps 1–14 + automated CMO quant)
    cmo_quantification/        # end-to-end CMO quantification from IGVF accessions
  10x_multi_15_timepoints/     # 15 tp pipeline scripts
  create_seqspecs/             # marimo notebook for generating seqspec YAMLs
config/
  conda/                       # conda environment definition and build sbatch
  splitcode/                   # splitcode tag configs for 5 tp and 15 tp datasets
metadata/                      # CMO design sheets (TSV)
scripts/
  bash/                        # exploratory / benchmarking scripts
  python/                      # shared Python utilities (IGVF auth, download, KITE helpers)
templates/                     # Jinja2 templates for seqspec YAML generation
```

## Environment

All tools (`kb-python`, `splitcode`, `seqspec`, `igvf-utils`) are in the project conda env. Build once:

```bash
sbatch config/conda/build_env.sbatch
```

Activate for interactive use:

```bash
source /home/users/emattei/miniforge3/etc/profile.d/conda.sh
conda activate /oak/stanford/groups/engreitz/Users/emattei/git/broad-nnfc-cmo-multiome/env/nnfc-cmo-multiome
```

## Authentication

IGVF downloads require a `.env` file at the repo root:

```
IGVF_API_KEY=<your_key>
IGVF_SECRET_KEY=<your_secret>
```
