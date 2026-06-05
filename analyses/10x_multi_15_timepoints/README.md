# CMO Multiome Processing Pipeline — 15 Timepoints

This README records the commands used to process the 10x multiome / CMO 15 timepoints sequencing data.

Unlike the 5 timepoints dataset, GEX and CMO reads are already demultiplexed into separate FASTQ files
by the sequencer (separate sample indices on the same flowcell). Therefore, **splitcode is not required**
to separate CMO from GEX reads. The pipeline proceeds directly from concatenation to CMO trimming and
quantification.

Additionally, this dataset includes ATAC FASTQs from a separate sequencing run.

## Library structure

```text
10X Lane 1 (channel 1):
  GEX: GEX_1_* files   -> channel1_GEX
  CMO: CMO_1_* files   -> channel1_CMO
  ATAC: ATAC_1_* files -> channel1_ATAC

10X Lane 2 (channel 2):
  GEX: GEX_2_* files   -> channel2_GEX
  CMO: CMO_2_* files   -> channel2_CMO
  ATAC: ATAC_2_* files -> channel2_ATAC
```

Do not merge channel1 and channel2 cell barcodes before downstream processing.

## Raw FASTQ locations

GEX and CMO FASTQs (flowcell H7WWTBGXN):
```
/oak/stanford/groups/engreitz/Users/munger/EC_Profiling/221212_Venous_GEX_CMO/sequence_processing_gex/outs/fastq_path/H7WWTBGXN/
```

ATAC FASTQs (flowcell H7Y2CBGXN):
```
/oak/stanford/groups/engreitz/Users/munger/EC_Profiling/221214_Venous_ATAC/221214_Venous_ATAC_atac/outs/fastq_path/H7Y2CBGXN/ATAC_1/
/oak/stanford/groups/engreitz/Users/munger/EC_Profiling/221214_Venous_ATAC/221214_Venous_ATAC_atac/outs/fastq_path/H7Y2CBGXN/ATAC_2/
```

A full list of all FASTQ paths is available at:
```
/oak/stanford/groups/engreitz/Users/munger/EC_Profiling/221212_Venous_GEX_CMO/sequence_processing_gex/outs/fastq_path/fastq_file_paths.csv
```

---

## Environment setup

```bash
conda create --prefix /home/groups/engreitz/Users/opushkar/.conda/envs/kb-splitcode -c conda-forge -c bioconda splitcode kb-python
source "/home/groups/engreitz/Software/anaconda3/etc/profile.d/conda.sh"
conda activate /home/groups/engreitz/Users/opushkar/.conda/envs/kb-splitcode
```

---

## 1. Create directories

```bash
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
```

---

## 2. Concatenate GEX_1 FASTQs (lanes L001–L004)

Channel 1 GEX reads are in files with prefix `GEX_1`.

```bash
cat ${GEX_CMO_RAW}/GEX_1*_R1_*.fastq.gz \
  > ${raw_fastq_path}/GEX1/GEX1_concatenated_R1.fastq.gz

cat ${GEX_CMO_RAW}/GEX_1*_R2_*.fastq.gz \
  > ${raw_fastq_path}/GEX1/GEX1_concatenated_R2.fastq.gz

```

---

## 3. Concatenate GEX_2 FASTQs (lanes L001–L004)

Channel 2 GEX reads are in files with prefix `GEX_2`.

```bash
cat ${GEX_CMO_RAW}/GEX_2*_R1_*.fastq.gz \
  > ${raw_fastq_path}/GEX2/GEX2_concatenated_R1.fastq.gz
cat ${GEX_CMO_RAW}/GEX_2*_R2_*.fastq.gz \
  > ${raw_fastq_path}/GEX2/GEX2_concatenated_R2.fastq.gz
```

---

## 4. Concatenate CMO_1 FASTQs (lanes L001–L004)

Channel 1 CMO reads are in files with prefix `CMO_1`.

```bash
cat ${GEX_CMO_RAW}/CMO_1*_R1_*.fastq.gz \
  > ${raw_fastq_path}/CMO1/CMO1_concatenated_R1.fastq.gz

cat ${GEX_CMO_RAW}/CMO_1*_R2_*.fastq.gz \
  > ${raw_fastq_path}/CMO1/CMO1_concatenated_R2.fastq.gz
```

---

## 5. Concatenate CMO_2 FASTQs (lanes L001–L004)

Channel 2 CMO reads are in files with prefix `CMO_2`.

```bash
cat ${GEX_CMO_RAW}/CMO_2*_R1_*.fastq.gz \
  > ${raw_fastq_path}/CMO2/CMO2_concatenated_R1.fastq.gz
cat ${GEX_CMO_RAW}/CMO_2*_R2_*.fastq.gz \
  > ${raw_fastq_path}/CMO2/CMO2_concatenated_R2.fastq.gz
```

---

## 6. Concatenate ATAC_1 FASTQs (lanes L001–L004)

10x multiome ATAC files use R1 (genomic read 1), R2 (barcode index), R3 (genomic read 2) naming.

```bash
cat ${ATAC_RAW}/ATAC_1/ATAC_1*_R1_*.fastq.gz \
  > ${raw_fastq_path}/ATAC1/ATAC1_concatenated_R1.fastq.gz

cat ${ATAC_RAW}/ATAC_1/ATAC_1*_R2_*.fastq.gz \
  > ${raw_fastq_path}/ATAC1/ATAC1_concatenated_R2.fastq.gz

cat ${ATAC_RAW}/ATAC_1/ATAC_1*_R3_*.fastq.gz \
  > ${raw_fastq_path}/ATAC1/ATAC1_concatenated_R3.fastq.gz
```

---

## 7. Concatenate ATAC_2 FASTQs (lanes L001–L004)

```bash
cat ${ATAC_RAW}/ATAC_2/ATAC_2*_R1_*.fastq.gz \
  > ${raw_fastq_path}/ATAC2/ATAC2_concatenated_R1.fastq.gz

cat ${ATAC_RAW}/ATAC_2/ATAC_2*_R2_*.fastq.gz \
  > ${raw_fastq_path}/ATAC2/ATAC2_concatenated_R2.fastq.gz

cat ${ATAC_RAW}/ATAC_2/ATAC_2*_R3_*.fastq.gz \
  > ${raw_fastq_path}/ATAC2/ATAC2_concatenated_R3.fastq.gz
```

---

## 8. Count reads in concatenated outputs

```bash
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
```

---

## 9. Inspect CMO barcode distributions

Question: Are the top sequences in the CMO FASTQs the expected CMO barcode sequences?
Expectation: Top sequences should match the known 15 timepoints CMO barcode sequences.

Channel 1 CMO:
All 23 CMOs are present in top 30
```bash
gzip -dc ${raw_fastq_path}/CMO1/CMO1_concatenated_R2.fastq.gz \
  | awk 'NR % 4 == 2 {print substr($0,1,15)}' \
  | sort | uniq -c | sort -nr | head -30
```

Channel 2 CMO:
22 exact matches in top 30; CMO23 is not used in channel 2

```bash
gzip -dc ${raw_fastq_path}/CMO2/CMO2_concatenated_R2.fastq.gz \
  | awk 'NR % 4 == 2 {print substr($0,1,15)}' \
  | sort | uniq -c | sort -nr | head -30
```

---

## 10. Inspect CMO R1/R2 read structure

Channel 1:
```bash
paste \
  <(gzip -dc ${raw_fastq_path}/CMO1/CMO1_concatenated_R1.fastq.gz | awk 'NR % 4 == 1 || NR % 4 == 2 {print; if (++c==20) exit}') \
  <(gzip -dc ${raw_fastq_path}/CMO1/CMO1_concatenated_R2.fastq.gz | awk 'NR % 4 == 1 || NR % 4 == 2 {print; if (++c==20) exit}')
```

Channel 2:
```bash
paste \
  <(gzip -dc ${raw_fastq_path}/CMO2/CMO2_concatenated_R1.fastq.gz | awk 'NR % 4 == 1 || NR % 4 == 2 {print; if (++c==20) exit}') \
  <(gzip -dc ${raw_fastq_path}/CMO2/CMO2_concatenated_R2.fastq.gz | awk 'NR % 4 == 1 || NR % 4 == 2 {print; if (++c==20) exit}')
```

Expected structure:

```text
CMO R1: variable leading A-run + 10x cell barcode + UMI + ...
CMO R2: CMO barcode + polyA/readthrough

Looks good
```

---

## 11. Trim CMO R1 to last 28 bp

The `kb kite` workflow expects R1 reads with cell barcode (16 bp) + UMI (12 bp) starting at position 0.
Trimming to the last 28 bp removes leading A-runs and other technical sequences.

Channel 1:

```bash
gzip -dc ${raw_fastq_path}/CMO1/CMO1_concatenated_R1.fastq.gz \
  | awk 'NR % 4 == 1 {h=$0}
         NR % 4 == 2 {s=substr($0, length($0)-27, 28)}
         NR % 4 == 3 {p=$0}
         NR % 4 == 0 {
         q=substr($0, length($0)-27, 28)
         print h "\n" s "\n" p "\n" q
        }' \
  | gzip > ${raw_fastq_path}/CMO1/CMO1_concatenated_R1.trimmed.fastq.gz
```

Channel 2:

```bash
gzip -dc ${raw_fastq_path}/CMO2/CMO2_concatenated_R1.fastq.gz \
  | awk 'NR % 4 == 1 {h=$0}
         NR % 4 == 2 {s=substr($0, length($0)-27, 28)}
         NR % 4 == 3 {p=$0}
         NR % 4 == 0 {
           q=substr($0, length($0)-27, 28)
           print h "\n" s "\n" p "\n" q
         }' \
  | gzip > ${raw_fastq_path}/CMO2/CMO2_concatenated_R1.trimmed.fastq.gz
```

---

## 12. Validate trimmed CMO R1 length

Channel 1:

```bash
gzip -dc ${OAK}/Projects/EC_Screen/Data/10x_15_timepoints/Raw_FASTQs/CMO1/CMO1_concatenated_R1.trimmed.fastq.gz \
  | awk 'NR % 4 == 2 {print length($0), $0; if (++n==20) exit}'
```

Channel 2:

```bash
gzip -dc ${OAK}/Projects/EC_Screen/Data/10x_15_timepoints/Raw_FASTQs/CMO2/CMO2_concatenated_R1.trimmed.fastq.gz \
  | awk 'NR % 4 == 2 {print length($0), $0; if (++n==20) exit}'
```

Expected length:

```text
28 -- data matches expectation
```

---

## 13. Link files to Processed_FASTQs

```bash
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
```

---

## 14. KITE workflow — CMO quantification

Prepare the `cmo_barcodes_for_kite.tsv` file:

```bash
scp /oak/stanford/groups/engreitz/Users/opushkar/broad-nnfc-cmo-multiome/config/10x_15_timepoints_cmo_23_splitcode.config ${output_path}/Metadata/10x_15_timepoints_cmo_23_splitcode.config

awk 'NR>1 && $1 !~ /^>/ {print $1"\t"$2}' \
  ${output_path}/Metadata/10x_15_timepoints_cmo_23_splitcode.config \
  > ${OAK}/Projects/EC_Screen/Data/10x_15_timepoints/Metadata/cmo_barcodes_for_kite.tsv
```

Create index for KITE:
so
```bash
```

Channel 1:

```bash
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
```

Channel 2:

```bash
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
```

---

## Conclusions

Final processed files for downstream analysis:

- channel1:
  - `Processed_FASTQs/channel1_GEX_R1.fastq.gz`
  - `Processed_FASTQs/channel1_GEX_R2.fastq.gz`
  - `Processed_FASTQs/channel1_CMO_R1.fastq.gz` (trimmed)
  - `Processed_FASTQs/channel1_CMO_R2.fastq.gz`
  - `Processed_FASTQs/channel1_ATAC_R1.fastq.gz`
  - `Processed_FASTQs/channel1_ATAC_R2.fastq.gz`
  - `Processed_FASTQs/channel1_ATAC_R3.fastq.gz`
- channel2:
  - `Processed_FASTQs/channel2_GEX_R1.fastq.gz`
  - `Processed_FASTQs/channel2_GEX_R2.fastq.gz`
  - `Processed_FASTQs/channel2_CMO_R1.fastq.gz` (trimmed)
  - `Processed_FASTQs/channel2_CMO_R2.fastq.gz`
  - `Processed_FASTQs/channel2_ATAC_R1.fastq.gz`
  - `Processed_FASTQs/channel2_ATAC_R2.fastq.gz`
  - `Processed_FASTQs/channel2_ATAC_R3.fastq.gz`
- channel1 CMO counts:
  - `CMO_counts/channel1/counts_unfiltered`
- channel2 CMO counts:
  - `CMO_counts/channel2/counts_unfiltered`

All paths above are relative to `${OAK}/Projects/EC_Screen/Data/10x_15_timepoints/`.
