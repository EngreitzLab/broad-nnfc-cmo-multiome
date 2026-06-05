# Endothelial cell differentiation CMO library processing

This repository contains the commands used to reconstruct per-channel GEX and CMO FASTQs from the mixed 10x multiome / CMO 5 timepoints sequencing data.

Summary of the steps:
1. Concatenate the FASTQs
2. Separate CMO from GEX with `splitcode`
3. Trim CMO R1 file
4. Quantify CMO tags per barcode using the `kite` workflow from `kb`



10X Genomics multi-ome 5 timepoints:

[README.md](analyses/10x_multi_5_timepoints/README.md)

10X Genomics multi-ome 15 timepoints:

[README.md](analyses/10x_multi_15_timepoints/README.md)
