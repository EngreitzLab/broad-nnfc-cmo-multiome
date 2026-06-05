#!/bin/bash
#SBATCH --job-name=multiome_15tp
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
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

mkdir -p ${raw_fastq_path}/GEX1
mkdir -p ${raw_fastq_path}/GEX2
mkdir -p ${raw_fastq_path}/CMO1
mkdir -p ${raw_fastq_path}/CMO2
mkdir -p ${raw_fastq_path}/ATAC1
mkdir -p ${raw_fastq_path}/ATAC2

mkdir -p ${output_path}/Metadata
mkdir -p ${output_path}/Kite
mkdir -p ${output_path}/CMO_counts/channel1
mkdir -p ${output_path}/CMO_counts/channel2
mkdir -p ${output_path}/Processed_FASTQs

# Concatenate GEX:
# Channel 1: GEX_1 fastqs
cat ${GEX_CMO_RAW}/GEX_1*_R1_*.fastq.gz \
  > ${raw_fastq_path}/GEX1/GEX1_concatenated_R1.fastq.gz

cat ${GEX_CMO_RAW}/GEX_1*_R2_*.fastq.gz \
  > ${raw_fastq_path}/GEX1/GEX1_concatenated_R2.fastq.gz

# Channel 2: GEX_2 fastqs
cat ${GEX_CMO_RAW}/GEX_2*_R1_*.fastq.gz \
  > ${raw_fastq_path}/GEX2/GEX2_concatenated_R1.fastq.gz
cat ${GEX_CMO_RAW}/GEX_2*_R2_*.fastq.gz \
  > ${raw_fastq_path}/GEX2/GEX2_concatenated_R2.fastq.gz


# Concatenate CMO:
# Channel 1: CMO_1 fastqs
cat ${GEX_CMO_RAW}/CMO_1*_R1_*.fastq.gz \
  > ${raw_fastq_path}/CMO1/CMO1_concatenated_R1.fastq.gz

cat ${GEX_CMO_RAW}/CMO_1*_R2_*.fastq.gz \
  > ${raw_fastq_path}/CMO1/CMO1_concatenated_R2.fastq.gz

# Channel 2: CMO_2 fastqs
cat ${GEX_CMO_RAW}/CMO_2*_R1_*.fastq.gz \
  > ${raw_fastq_path}/CMO2/CMO2_concatenated_R1.fastq.gz
cat ${GEX_CMO_RAW}/CMO_2*_R2_*.fastq.gz \
  > ${raw_fastq_path}/CMO2/CMO2_concatenated_R2.fastq.gz

# Concatenate ATAC:
# Channel 1: ATAC_1 fastqs
cat ${ATAC_RAW}/ATAC_1/ATAC_1*_R1_*.fastq.gz \
  > ${raw_fastq_path}/ATAC1/ATAC1_concatenated_R1.fastq.gz

cat ${ATAC_RAW}/ATAC_1/ATAC_1*_R2_*.fastq.gz \
  > ${raw_fastq_path}/ATAC1/ATAC1_concatenated_R2.fastq.gz

cat ${ATAC_RAW}/ATAC_1/ATAC_1*_R3_*.fastq.gz \
  > ${raw_fastq_path}/ATAC1/ATAC1_concatenated_R3.fastq.gz

# Channel 2: ATAC_2 fastqs
cat ${ATAC_RAW}/ATAC_2/ATAC_2*_R1_*.fastq.gz \
  > ${raw_fastq_path}/ATAC2/ATAC2_concatenated_R1.fastq.gz

cat ${ATAC_RAW}/ATAC_2/ATAC_2*_R2_*.fastq.gz \
  > ${raw_fastq_path}/ATAC2/ATAC2_concatenated_R2.fastq.gz

cat ${ATAC_RAW}/ATAC_2/ATAC_2*_R3_*.fastq.gz \
  > ${raw_fastq_path}/ATAC2/ATAC2_concatenated_R3.fastq.gz

# Calculate n reads in concatenated outputs
for f in \
  ${raw_fastq_path}/GEX1/GEX1_concatenated_R1.fastq.gz \
  ${raw_fastq_path}/GEX2/GEX2_concatenated_R1.fastq.gz \
  ${raw_fastq_path}/CMO1/CMO1_concatenated_R1.fastq.gz \
  ${raw_fastq_path}/CMO2/CMO2_concatenated_R1.fastq.gz \
  ${raw_fastq_path}/ATAC1/ATAC1_concatenated_R1.fastq.gz \
  ${raw_fastq_path}/ATAC2/ATAC2_concatenated_R1.fastq.gz;
do
  echo "$f"
  gzip -dc "$f" | awk 'END {print NR/4}'
done

# Explore CMO barcode distribution
gzip -dc ${raw_fastq_path}/CMO1/CMO1_concatenated_R2.fastq.gz \
  | awk 'NR % 4 == 2 {print substr($0,1,15)}' \
  | sort | uniq -c | sort -nr | head -30

gzip -dc ${raw_fastq_path}/CMO2/CMO2_concatenated_R2.fastq.gz \
  | awk 'NR % 4 == 2 {print substr($0,1,15)}' \
  | sort | uniq -c | sort -nr | head -30

# Check CMO read structure
paste \
  <(gzip -dc ${raw_fastq_path}/CMO1/CMO1_concatenated_R1.fastq.gz | awk 'NR % 4 == 1 || NR % 4 == 2 {print; if (++c==20) exit}') \
  <(gzip -dc ${raw_fastq_path}/CMO1/CMO1_concatenated_R2.fastq.gz | awk 'NR % 4 == 1 || NR % 4 == 2 {print; if (++c==20) exit}')

paste \
  <(gzip -dc ${raw_fastq_path}/CMO2/CMO2_concatenated_R1.fastq.gz | awk 'NR % 4 == 1 || NR % 4 == 2 {print; if (++c==20) exit}') \
  <(gzip -dc ${raw_fastq_path}/CMO2/CMO2_concatenated_R2.fastq.gz | awk 'NR % 4 == 1 || NR % 4 == 2 {print; if (++c==20) exit}')

# Trim CMO R1 to last 28 bp
# Channel 1:
gzip -dc ${raw_fastq_path}/CMO1/CMO1_concatenated_R1.fastq.gz \
  | awk 'NR % 4 == 1 {h=$0}
         NR % 4 == 2 {s=substr($0, length($0)-27, 28)}
         NR % 4 == 3 {p=$0}
         NR % 4 == 0 {
         q=substr($0, length($0)-27, 28)
         print h "\n" s "\n" p "\n" q
        }' \
  | gzip > ${raw_fastq_path}/CMO1/CMO1_concatenated_R1.trimmed.fastq.gz

# Channel 2:
gzip -dc ${raw_fastq_path}/CMO2/CMO2_concatenated_R1.fastq.gz \
  | awk 'NR % 4 == 1 {h=$0}
         NR % 4 == 2 {s=substr($0, length($0)-27, 28)}
         NR % 4 == 3 {p=$0}
         NR % 4 == 0 {
           q=substr($0, length($0)-27, 28)
           print h "\n" s "\n" p "\n" q
         }' \
  | gzip > ${raw_fastq_path}/CMO2/CMO2_concatenated_R1.trimmed.fastq.gz

# Validate length
# Channel 1:
gzip -dc ${raw_fastq_path}/CMO1/CMO1_concatenated_R1.trimmed.fastq.gz \
  | awk 'NR % 4 == 2 {print length($0), $0; if (++n==20) exit}'

# Channel 2:
gzip -dc ${raw_fastq_path}/CMO2/CMO2_concatenated_R1.trimmed.fastq.gz \
  | awk 'NR % 4 == 2 {print length($0), $0; if (++n==20) exit}'


# Symlink
# GEX channel 1
ln -s ${raw_fastq_path}/GEX1/GEX1_concatenated_R1.fastq.gz \
  ${output_path}/Processed_FASTQs/channel1_GEX_R1.fastq.gz

ln -s ${raw_fastq_path}/GEX1/GEX1_concatenated_R2.fastq.gz \
  ${output_path}/Processed_FASTQs/channel1_GEX_R2.fastq.gz

# GEX channel 2
ln -s ${raw_fastq_path}/GEX2/GEX2_concatenated_R1.fastq.gz \
  ${output_path}/Processed_FASTQs/channel2_GEX_R1.fastq.gz

ln -s ${raw_fastq_path}/GEX2/GEX2_concatenated_R2.fastq.gz \
  ${output_path}/Processed_FASTQs/channel2_GEX_R2.fastq.gz

# CMO (trimmed R1) channel 1
ln -s ${raw_fastq_path}/CMO1/CMO1_concatenated_R1.trimmed.fastq.gz \
  ${output_path}/Processed_FASTQs/channel1_CMO_R1.fastq.gz

ln -s ${raw_fastq_path}/CMO1/CMO1_concatenated_R2.fastq.gz \
  ${output_path}/Processed_FASTQs/channel1_CMO_R2.fastq.gz

# CMO (trimmed R1) channel 2
ln -s ${raw_fastq_path}/CMO2/CMO2_concatenated_R1.trimmed.fastq.gz \
  ${output_path}/Processed_FASTQs/channel2_CMO_R1.fastq.gz
  
ln -s ${raw_fastq_path}/CMO2/CMO2_concatenated_R2.fastq.gz \
  ${output_path}/Processed_FASTQs/channel2_CMO_R2.fastq.gz

# ATAC channel 1
ln -s ${raw_fastq_path}/ATAC1/ATAC1_concatenated_R1.fastq.gz \
  ${output_path}/Processed_FASTQs/channel1_ATAC_R1.fastq.gz

ln -s ${raw_fastq_path}/ATAC1/ATAC1_concatenated_R2.fastq.gz \
  ${output_path}/Processed_FASTQs/channel1_ATAC_R2.fastq.gz

ln -s ${raw_fastq_path}/ATAC1/ATAC1_concatenated_R3.fastq.gz \
  ${output_path}/Processed_FASTQs/channel1_ATAC_R3.fastq.gz

# ATAC channel 2
ln -s ${raw_fastq_path}/ATAC2/ATAC2_concatenated_R1.fastq.gz \
  ${output_path}/Processed_FASTQs/channel2_ATAC_R1.fastq.gz

ln -s ${raw_fastq_path}/ATAC2/ATAC2_concatenated_R2.fastq.gz \
  ${output_path}/Processed_FASTQs/channel2_ATAC_R2.fastq.gz

ln -s ${raw_fastq_path}/ATAC2/ATAC2_concatenated_R3.fastq.gz \
  ${output_path}/Processed_FASTQs/channel2_ATAC_R3.fastq.gz