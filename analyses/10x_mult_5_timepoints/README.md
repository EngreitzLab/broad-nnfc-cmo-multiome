# CMO Multiome Processing Pipeline

This README records the commands used to reconstruct per-channel GEX and CMO FASTQs from the mixed 10x multiome / CMO 5 timepoints sequencing data.

Important mapping from the experiment notes:

```text
GEX_1 input:
  CMO-matching reads     -> channel1_CMO
  non-CMO reads          -> channel2_GEX

GEX_2 input:
  CMO-matching reads     -> channel2_CMO
  non-CMO reads          -> channel1_GEX
```

Do not merge channel1 and channel2 cell barcodes before downstream processing.

---

## Environment setup

```bash
conda create -n kb-splitcode -c conda-forge -c bioconda splitcode kb-python
conda activate kb-splitcode
```

---

## 1. Create directories



```bash
mkdir -p ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/GEX1/Raw_FASTQs/H2YL5BGXM
mkdir -p ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/GEX1/Raw_FASTQs/H2YL5BGXM/GEX1/HF7JFBGXM
mkdir -p ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/GEX1/Raw_FASTQs/H2YL5BGXM/GEX2/H2YL5BGXM

mkdir -p ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Metadata
mkdir -p ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode

mkdir -p ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/GEX1/Processed_FASTQs/
```

---

## 2. Download GEX_1 FASTQs from GCS

GEX_1 contains:
- channel1_CMO reads
- channel2_GEX reads

It also has additional sequencing in `HF7JFBGXM`.

```bash
gsutil -m cp \
  'gs:/fc-76565551-28c3-4b6d-a048-e272103bcbd1/data/10x_multiome_5_timepoints/H2YL5BGXM/GEX_1*R*_*.fastq.gz' \
  $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX1/H2YL5BGXM/

gsutil -m cp \
  'gs:/fc-76565551-28c3-4b6d-a048-e272103bcbd1/data/10x_multiome_5_timepoints/HF7JFBGXM/GEX_1*R*_*.fastq.gz' \
  $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX1/HF7JFBGXM/
```

---

## 3. Download GEX_2 FASTQs from GCS

GEX_2 contains:
- channel2_CMO reads
- channel1_GEX reads

```bash
gsutil -m cp \
  'gs:/fc-76565551-28c3-4b6d-a048-e272103bcbd1/data/10x_multiome_5_timepoints/H2YL5BGXM/GEX_2*R*_*.fastq.gz' \
  $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX2/H2YL5BGXM/
```

If there is additional GEX_2 sequencing in another run, download it into a separate subdirectory under `$GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX2/` and include it in the concatenation step below.

---

## 4. Concatenate GEX_1 FASTQs

Concatenate technical lanes / resequencing for the same input.

```bash
cat \
  $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX1/H2YL5BGXM/*_R1_*.fastq.gz \
  $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX1/HF7JFBGXM/*_R1_*.fastq.gz \
  > $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX1/GEX1_concatenated_R1.fastq.gz

cat \
  $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX1/H2YL5BGXM/*_R2_*.fastq.gz \
  $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX1/HF7JFBGXM/*_R2_*.fastq.gz \
  > $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX1/GEX1_concatenated_R2.fastq.gz
```

---

## 5. Concatenate GEX_2 FASTQs

```bash
cat \
  $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX2/H2YL5BGXM/*_R1_*.fastq.gz \
  > $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX2/GEX2_concatenated_R1.fastq.gz

cat \
  $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX2/H2YL5BGXM/*_R2_*.fastq.gz \
  > $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX2/GEX2_concatenated_R2.fastq.gz
```

If there are additional GEX_2 sequencing FASTQs, add them to these `cat` commands.

---

## 6. Split GEX_1 into channel1_CMO and channel2_GEX

```bash
splitcode \
  -c $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/config/cmo_25_splitcode.config \
  --nFastqs=2 \
  --gzip \
  --assign \
  --summary $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/GEX1_split_summary.txt \
  --mapping $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/GEX1_split_mapping.txt \
  -o $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel1_CMO_R1.fastq.gz,$GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel1_CMO_R2.fastq.gz \
  -u $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel2_GEX_R1.fastq.gz,$GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel2_GEX_R2.fastq.gz \
  $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX1/GEX1_concatenated_R1.fastq.gz \
  $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX1/GEX1_concatenated_R2.fastq.gz
```

---

## 7. Split GEX_2 into channel2_CMO and channel1_GEX

```bash
splitcode \
  -c $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/config/cmo_25_splitcode.config \
  --nFastqs=2 \
  --gzip \
  --assign \
  --summary $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/GEX2_split_summary.txt \
  --mapping $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/GEX2_split_mapping.txt \
  -o $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel2_CMO_R1.fastq.gz,$GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel2_CMO_R2.fastq.gz \
  -u $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel1_GEX_R1.fastq.gz,$GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel1_GEX_R2.fastq.gz \
  $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX2/GEX2_concatenated_R1.fastq.gz \
  $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/GEX2/GEX2_concatenated_R2.fastq.gz
```

---

## 8. Count reads in split outputs

```bash
for f in $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/*.fastq.gz; do
  echo "$f"
  gzip -dc "$f" | awk 'END {print NR/4}'
done
```

---

## 9. Inspect CMO barcode distributions

Channel 1 CMO:

```bash
gzip -dc $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel1_CMO_R2.fastq.gz \
  | awk 'NR % 4 == 2 {print substr($0,1,15)}' \
  | sort | uniq -c | sort -nr | head -30
```

Channel 2 CMO:

```bash
gzip -dc $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel2_CMO_R2.fastq.gz \
  | awk 'NR % 4 == 2 {print substr($0,1,15)}' \
  | sort | uniq -c | sort -nr | head -30
```

Expected: the top sequences should match the 25 CMO barcode sequences.

---

## 10. Inspect CMO R1/R2 read structure

Channel 1:

```bash
paste \
  <(gzip -dc $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel1_CMO_R1.fastq.gz | awk 'NR % 4 == 1 || NR % 4 == 2 {print; if (++c==20) exit}') \
  <(gzip -dc $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel1_CMO_R2.fastq.gz | awk 'NR % 4 == 1 || NR % 4 == 2 {print; if (++c==20) exit}')
```

Channel 2:

```bash
paste \
  <(gzip -dc $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel2_CMO_R1.fastq.gz | awk 'NR % 4 == 1 || NR % 4 == 2 {print; if (++c==20) exit}') \
  <(gzip -dc $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel2_CMO_R2.fastq.gz | awk 'NR % 4 == 1 || NR % 4 == 2 {print; if (++c==20) exit}')
```

Observed structure:

```text
CMO R1: variable leading A-run + 10x cell barcode + UMI + ...
CMO R2: CMO barcode + polyA/readthrough
```

---

## 11. Compare candidate CMO barcodes to GEX barcodes using a 1M read subsample

This is only a diagnostic step.

Channel 1:

```bash
gzip -dc $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel1_CMO_R1.fastq.gz \
  | awk 'NR % 4 == 2 {
      seq=$0; sub(/^A+/, "", seq);
      print substr(seq,1,16);
      if (++n==1000000) exit
    }' \
  | sort -u > /tmp/channel1_cmo_cb_1M.txt

gzip -dc $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel1_GEX_R1.fastq.gz \
  | awk 'NR % 4 == 2 {
      print substr($0,1,16);
      if (++n==1000000) exit
    }' \
  | sort -u > /tmp/channel1_gex_cb_1M.txt

comm -12 /tmp/channel1_cmo_cb_1M.txt /tmp/channel1_gex_cb_1M.txt | wc -l
wc -l /tmp/channel1_cmo_cb_1M.txt /tmp/channel1_gex_cb_1M.txt
```

Channel 2:

```bash
gzip -dc $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel2_CMO_R1.fastq.gz \
  | awk 'NR % 4 == 2 {
      seq=$0; sub(/^A+/, "", seq);
      print substr(seq,1,16);
      if (++n==1000000) exit
    }' \
  | sort -u > /tmp/channel2_cmo_cb_1M.txt

gzip -dc $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel2_GEX_R1.fastq.gz \
  | awk 'NR % 4 == 2 {
      print substr($0,1,16);
      if (++n==1000000) exit
    }' \
  | sort -u > /tmp/channel2_gex_cb_1M.txt

comm -12 /tmp/channel2_cmo_cb_1M.txt /tmp/channel2_gex_cb_1M.txt | wc -l
wc -l /tmp/channel2_cmo_cb_1M.txt /tmp/channel2_gex_cb_1M.txt
```

---

## 12. Create trimmed R1 from the last 28 bp of CMO R1

Goal:

```text
trimmed R1 = last 28 bp of CMO R1
```

This step does not perform barcode whitelist matching. It directly keeps the 3' terminal 28 bp from each CMO R1 read.

Channel 1:

```bash
gzip -dc $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel1_CMO_R1.fastq.gz \
  | awk 'NR % 4 == 1 {h=$0}
         NR % 4 == 2 {s=substr($0, length($0)-27, 28)}
         NR % 4 == 3 {p=$0}
         NR % 4 == 0 {
           q=substr($0, length($0)-27, 28)
           print h "\n" s "\n" p "\n" q
         }' \
  | gzip > $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel1_CMO_R1.trimmed.fastq.gz
```

Channel 2:

```bash
gzip -dc $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel2_CMO_R1.fastq.gz \
  | awk 'NR % 4 == 1 {h=$0}
         NR % 4 == 2 {s=substr($0, length($0)-27, 28)}
         NR % 4 == 3 {p=$0}
         NR % 4 == 0 {
           q=substr($0, length($0)-27, 28)
           print h "\n" s "\n" p "\n" q
         }' \
  | gzip > $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel2_CMO_R1.trimmed.fastq.gz
```

---

## 13. Validate channel 1 trimmed R1

```bash
gzip -dc $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel1_CMO_R1.trimmed.fastq.gz \
  | awk 'NR % 4 == 2 {print length($0), $0; if (++n==20) exit}'
```

Expected length:

```text
28
```

---

## 14. Validate channel 2 trimmed R1

```bash
gzip -dc $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel2_CMO_R1.trimmed.fastq.gz \
  | awk 'NR % 4 == 2 {print length($0), $0; if (++n==20) exit}'
```

Expected length:

```text
28
```

---

## 15. KITE workflow

Download onlist file:

```bash
wget https://api.data.igvf.org/tabular-files/IGVFFI8751YQRY/@@download/IGVFFI8751YQRY.tsv.gz \
  -O $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/onlist/IGVFFI8751YQRY.tsv.gz
gunzip $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/onlist/IGVFFI8751YQRY.tsv.gz
```

Prepare the `cmo_barcodes.tsv` file for KITE:

```bash
awk 'NR>1 && $1 !~ /^>/ {print $1"\t"$2}' config/cmo_25_splitcode.config \
  > kite/cmo_barcodes.tsv
```

Create index for KITE:

```bash
kb ref \
  --workflow kite \
  -i kite/cmo.idx \
  -g kite/cmo_t2g.txt \
  -f1 kite/cmo_cdna.fa \
  kite/cmo_barcodes.tsv
```

Channel 1:

```bash
kb count \
  --workflow kite \
  -i kite/cmo.idx \
  -g kite/cmo_t2g.txt \
  -x 0,0,16:0,16,28:1,0,15 \
  -w onlist/IGVFFI8751YQRY.tsv \
  --h5ad \
  -o cmo_counts/channel1 \
  $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel1_CMO_R1.trimmed.fastq.gz \
  $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel1_CMO_R2.fastq.gz
```

Channel 2:

```bash
kb count \
  --workflow kite \
  -i kite/cmo.idx \
  -g kite/cmo_t2g.txt \
  -x 0,0,16:0,16,28:1,0,15 \
  -w onlist/IGVFFI8751YQRY.tsv \
  --h5ad \
  -o cmo_counts/channel2 \
  $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel2_CMO_R1.trimmed.fastq.gz \
  $GROUP_HOME/Data/endothelial-differentiation/10x_multiome_5_timepoints/splitcode/channel2_CMO_R2.fastq.gz
```



gsutil -m cp ch
annel1_GEX_R1.fastq.gz channel1_GEX_R2.fastq.gz channel2_GEX_R1.fastq.gz channel2_GEX_R2.fastq.gz  gs://fc-76565551-28c3-4b6d-a048-e272103bcbd1/data/10x_multiome_5_timepoints/processed/