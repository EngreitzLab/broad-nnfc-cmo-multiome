#!/bin/bash
#SBATCH --job-name=kb_multiome_15tp
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --partition=normal,engreitz
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.log

# conda create --prefix /home/groups/engreitz/Users/opushkar/.conda/envs/kb-splitcode -c conda-forge -c bioconda splitcode kb-python
source "/home/groups/engreitz/Software/anaconda3/etc/profile.d/conda.sh"
conda activate /home/groups/engreitz/Users/opushkar/.conda/envs/kb-splitcode

GEX_CMO_RAW="/oak/stanford/groups/engreitz/Users/munger/EC_Profiling/221212_Venous_GEX_CMO/sequence_processing_gex/outs/fastq_path/H7WWTBGXN"
ATAC_RAW="/oak/stanford/groups/engreitz/Users/munger/EC_Profiling/221214_Venous_ATAC/221214_Venous_ATAC_atac/outs/fastq_path/H7Y2CBGXN"

output_path="${OAK}/Projects/EC_Screen/Data/10x_15_timepoints"
raw_fastq_path="${output_path}/Raw_FASTQs"

kb ref \
    --workflow kite \
    -i ${output_path}/Kite/cmo.idx \
    -g ${output_path}/Kite/cmo_t2g.txt \
    -f1 ${output_path}/Kite/cmo_cdna.fa \
    ${output_path}/Metadata/cmo_barcodes_for_kite.tsv

kb count \
    --workflow kite \
    -i ${output_path}/Kite/cmo.idx \
    -g ${output_path}/Kite/cmo_t2g.txt \
    -x 0,0,16:0,16,28:1,0,15 \
    -w ${output_path}/Metadata/IGVFFI8751YQRY.tsv \
    --h5ad \
    -o ${output_path}/CMO_counts/channel1 \
    ${output_path}/Processed_FASTQs/channel1_CMO_R1.fastq.gz \
    ${output_path}/Processed_FASTQs/channel1_CMO_R2.fastq.gz


kb count \
    --workflow kite \
    -i ${output_path}/Kite/cmo.idx \
    -g ${output_path}/Kite/cmo_t2g.txt \
    -x 0,0,16:0,16,28:1,0,15 \
    -w ${output_path}/Metadata/IGVFFI8751YQRY.tsv \
    --h5ad \
    -o ${output_path}/CMO_counts/channel2 \
    ${output_path}/Processed_FASTQs/channel2_CMO_R1.fastq.gz \
    ${output_path}/Processed_FASTQs/channel2_CMO_R2.fastq.gz