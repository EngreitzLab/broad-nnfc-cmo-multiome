# CMO Multiome Processing Pipeline

> **Processed data on IGVF portal**
> All processed files for this dataset are available on the IGVF portal as analysis sets:
> - Channel 1: [IGVFDS5477BPOI](https://data.igvf.org/analysis-sets/IGVFDS5477BPOI)
> - Channel 2: [IGVFDS3995WHFT](https://data.igvf.org/analysis-sets/IGVFDS3995WHFT)

This README records the commands used to reconstruct per-channel GEX and CMO FASTQs from the mixed 10x multiome / CMO 5 timepoints sequencing data.

The reason why this is necessary is that the 10x multiome 5 timepoints data was sequenced in a way that mixes CMO and GEX reads in the same FASTQ files. The CMO and GEX reads can be separated based on the presence of CMO barcode sequences, but this requires custom processing steps.

This is how the libraries are structured in the FASTQ files:
```text
GEX_1 R2 file contains two different libraries:
  CMO-matching reads     -> channel1_CMO
  non-CMO reads          -> channel2_GEX

GEX_2 R2 file contains two different libraries:
  CMO-matching reads     -> channel2_CMO
  non-CMO reads          -> channel1_GEX
```

Pay attention to the fact that the GEX-matching reads in GEX_1 correspond to GEX_2 and vice versa.

Do not merge channel1 and channel2 cell barcodes before downstream processing.

---

## Environment setup

The project conda environment contains all required tools (`kb-python`, `splitcode`, `seqspec`, `igvf-utils`). Build it once with:

```bash
sbatch /oak/stanford/groups/engreitz/Users/emattei/git/broad-nnfc-cmo-multiome/config/conda/build_env.sbatch
```

Activate for interactive use:

```bash
source /home/users/emattei/miniforge3/etc/profile.d/conda.sh
conda activate /oak/stanford/groups/engreitz/Users/emattei/git/broad-nnfc-cmo-multiome/env/nnfc-cmo-multiome
```

---

## 1. Create directories



```bash
mkdir -p ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX1/HF7JFBGXM
mkdir -p ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX1/H2YL5BGXM
mkdir -p ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX2/H2YL5BGXM

mkdir -p ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Metadata
mkdir -p ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode
mkdir -p ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Kite
mkdir -p ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/CMO_counts/channel1
mkdir -p ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/CMO_counts/channel2

mkdir -p ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Processed_FASTQs/
```

---

## 2. Moving GEX_1 FASTQs to one central location

There might be multiple copies of the same FASTQs exist in different locations. I will move all GEX_1 FASTQs to one central location and create symlinks back to the original locations. This way, we can be sure that we are using all available GEX_1 sequencing data without creating duplicate files.

GEX_1 R2 file contains:
- channel1_CMO reads
- channel2_GEX reads

It also has additional sequencing in `HF7JFBGXM`.

```bash
NEW_BASE=${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX1

OLD_H2=/oak/stanford/groups/engreitz/Users/munger/Demux/220523_Demultiplexing_EC_Profiling_GEX_CMO/sequence_processing_gex/outs/fastq_path/H2YL5BGXM
OLD_HF=/oak/stanford/groups/engreitz/Users/dulguun/220519_EC_Profiling/sequence_processing_gex/outs/fastq_path/HF7JFBGXM

NEW_H2=${NEW_BASE}/H2YL5BGXM
NEW_HF=${NEW_BASE}/HF7JFBGXM

# Create destination directories
mkdir -p "${NEW_H2}" "${NEW_HF}"

# Move files
mv ${OLD_H2}/GEX_1*R*_*.fastq.gz "${NEW_H2}/"
mv ${OLD_HF}/GEX_1*R*_*.fastq.gz "${NEW_HF}/"

# Create symlinks back in original locations
for f in ${NEW_H2}/GEX_1*R*_*.fastq.gz; do
    ln -s "$f" "${OLD_H2}/$(basename "$f")"
done

for f in ${NEW_HF}/GEX_1*R*_*.fastq.gz; do
    ln -s "$f" "${OLD_HF}/$(basename "$f")"
done
```

---

## 3. Moving GEX_2 FASTQs to one central location

There might be multiple copies of the same FASTQs exist in different locations. I will move all GEX_2 FASTQs to one central location and create symlinks back to the original locations. This way, we can be sure that we are using all available GEX_2 sequencing data without creating duplicate files.

GEX_2 R2 file contains:
- channel2_CMO reads
- channel1_GEX reads

```bash
# Define paths
OLD_GEX2=/oak/stanford/groups/engreitz/Users/munger/Demux/220523_Demultiplexing_EC_Profiling_GEX_CMO/sequence_processing_gex/outs/fastq_path/H2YL5BGXM

NEW_GEX2=${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX2/H2YL5BGXM

# Create destination directory
mkdir -p "${NEW_GEX2}"

# Move GEX_2 R1/R2 FASTQs
mv ${OLD_GEX2}/GEX_2*R*_*.fastq.gz "${NEW_GEX2}/"

# Create symlinks back to original location
for f in ${NEW_GEX2}/GEX_2*R*_*.fastq.gz; do
    ln -s "$f" "${OLD_GEX2}/$(basename "$f")"
done
```

If there is additional GEX_2 sequencing in another run, download it into a separate subdirectory under `${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX2/<new_subdirectory>` and include it in the concatenation step below.

---

## 4. Concatenate GEX_1 FASTQs

Concatenate technical lanes / resequencing for the same input.

```bash
cat \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX1/H2YL5BGXM/*_R1_*.fastq.gz \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX1/HF7JFBGXM/*_R1_*.fastq.gz \
  > ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX1/GEX1_concatenated_R1.fastq.gz

cat \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX1/H2YL5BGXM/*_R2_*.fastq.gz \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX1/HF7JFBGXM/*_R2_*.fastq.gz \
  > ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX1/GEX1_concatenated_R2.fastq.gz
```

---

## 5. Concatenate GEX_2 FASTQs

```bash
cat \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX2/H2YL5BGXM/*_R1_*.fastq.gz \
  > ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX2/GEX2_concatenated_R1.fastq.gz

cat \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX2/H2YL5BGXM/*_R2_*.fastq.gz \
  > ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX2/GEX2_concatenated_R2.fastq.gz
```

If there are additional GEX_2 sequencing FASTQs, add them to these `cat` commands.

---


## Splitcode

[Splitcode](https://github.com/pachterlab/splitcode) is a tool created by Delaney Sullivan and Lior Pachter that enables flexible and efficient parsing, interpreting and editing of sequencing reads according to a user’s specifications

Reference:

Sullivan DK, Pachter L. (2024). Flexible parsing, interpretation, and editing of technical sequences with splitcode. Bioinformatics. https://doi.org/10.1093/bioinformatics/btae331

For our CMO case, splitcode will read a configuration file like this one:

```text
tags    ids     groups  distance        locations       minFindsG       maxFindsG
TTGTCACGGTAATTA CMO01   CMO     2       1:0:40  1       1
ATCGAACCGACAGAG CMO02   CMO     2       1:0:40  1       1
```

splitcode will go through each read in the R2 file, look for sequences that match the `tags` (CMO barcode sequences) within the specified `distance` and `locations`, and assign reads to groups based on the matching results. Then, it will output separate FASTQ files for each group.
All the reads matching a `tag` will be written to the CMO `group` FASTQ output file, and all the reads that do not match any `tag` will be written to the GEX FASTQ output file.

## 6. Split GEX_1 into channel1_CMO and channel2_GEX

**Reminder** GEX_1 R2 file contains:
- channel1_CMO reads
- channel2_GEX reads

```bash
splitcode \
  -c ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Metadata/10x_5_timepoints_cmo_25_splitcode.config \
  --nFastqs=2 \
  --gzip \
  --assign \
  --summary ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/GEX1_splitcode_summary.log.txt \
  --mapping ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/GEX1_splitcode_mapping.log.txt \
  -o ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel1_CMO_R1.fastq.gz,${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel1_CMO_R2.fastq.gz \
  -u ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel2_GEX_R1.fastq.gz,${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel2_GEX_R2.fastq.gz \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX1/GEX1_concatenated_R1.fastq.gz \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX1/GEX1_concatenated_R2.fastq.gz
```

---

## 7. Split GEX_2 into channel2_CMO and channel1_GEX

**Reminder** GEX_2 R2 filecontains:
- channel2_CMO reads
- channel1_GEX reads

```bash
splitcode \
  -c ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Metadata/10x_5_timepoints_cmo_25_splitcode.config \
  --nFastqs=2 \
  --gzip \
  --assign \
  --summary ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/GEX2_splitcode_summary.log.txt \
  --mapping ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/GEX2_splitcode_mapping.log.txt \
  -o ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel2_CMO_R1.fastq.gz,${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel2_CMO_R2.fastq.gz \
  -u ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel1_GEX_R1.fastq.gz,${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel1_GEX_R2.fastq.gz \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX2/GEX2_concatenated_R1.fastq.gz \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Raw_FASTQs/GEX2/GEX2_concatenated_R2.fastq.gz
```

---

## 8. Count reads in split outputs

```bash
for f in ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/*.fastq.gz; do
  echo "$f"
  gzip -dc "$f" | awk 'END {print NR/4}'
done
```

---

## 9. Inspect CMO barcode distributions

Question: Are the top sequences in the newly created CMO FASTQs the expected CMO barcode sequences?
Expectation: The top sequences should match the 25 CMO barcode sequences.

Channel 1 CMO:

```bash
gzip -dc ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel1_CMO_R2.fastq.gz \
  | awk 'NR % 4 == 2 {print substr($0,1,15)}' \
  | sort | uniq -c | sort -nr | head -30
```

Channel 2 CMO:

```bash
gzip -dc ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel2_CMO_R2.fastq.gz \
  | awk 'NR % 4 == 2 {print substr($0,1,15)}' \
  | sort | uniq -c | sort -nr | head -30
```



---

## 10. Inspect CMO R1/R2 read structure

Channel 1:

```bash
paste \
  <(gzip -dc ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel1_CMO_R1.fastq.gz | awk 'NR % 4 == 1 || NR % 4 == 2 {print; if (++c==20) exit}') \
  <(gzip -dc ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel1_CMO_R2.fastq.gz | awk 'NR % 4 == 1 || NR % 4 == 2 {print; if (++c==20) exit}')
```

Channel 2:

```bash
paste \
  <(gzip -dc ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel2_CMO_R1.fastq.gz | awk 'NR % 4 == 1 || NR % 4 == 2 {print; if (++c==20) exit}') \
  <(gzip -dc ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel2_CMO_R2.fastq.gz | awk 'NR % 4 == 1 || NR % 4 == 2 {print; if (++c==20) exit}')
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
gzip -dc ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel1_CMO_R1.fastq.gz \
  | awk 'NR % 4 == 2 {
      seq=$0; sub(/^A+/, "", seq);
      print substr(seq,1,16);
      if (++n==1000000) exit
    }' \
  | sort -u > /tmp/channel1_cmo_cb_1M.txt

gzip -dc ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel1_GEX_R1.fastq.gz \
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
gzip -dc ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel2_CMO_R1.fastq.gz \
  | awk 'NR % 4 == 2 {
      seq=$0; sub(/^A+/, "", seq);
      print substr(seq,1,16);
      if (++n==1000000) exit
    }' \
  | sort -u > /tmp/channel2_cmo_cb_1M.txt

gzip -dc ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel2_GEX_R1.fastq.gz \
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
The `kb kite` workflow used later to quantify CMO counts per barcode will work better with 28bp R1 reads that contain the full cell barcode and UMI sequences without leading A-runs or other technical sequences. This is because the `kb kite` workflow expects the cell barcode and UMI to be in specific positions within the R1 read, and having extra leading sequences can interfere with this.

Channel 1:

```bash
gzip -dc ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel1_CMO_R1.fastq.gz \
  | awk 'NR % 4 == 1 {h=$0}
         NR % 4 == 2 {s=substr($0, length($0)-27, 28)}
         NR % 4 == 3 {p=$0}
         NR % 4 == 0 {
           q=substr($0, length($0)-27, 28)
           print h "\n" s "\n" p "\n" q
         }' \
  | gzip > ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel1_CMO_R1.trimmed.fastq.gz
```

Channel 2:

```bash
gzip -dc ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel2_CMO_R1.fastq.gz \
  | awk 'NR % 4 == 1 {h=$0}
         NR % 4 == 2 {s=substr($0, length($0)-27, 28)}
         NR % 4 == 3 {p=$0}
         NR % 4 == 0 {
           q=substr($0, length($0)-27, 28)
           print h "\n" s "\n" p "\n" q
         }' \
  | gzip > ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel2_CMO_R1.trimmed.fastq.gz
```

---

## 13. Validate channel 1 trimmed R1

```bash
gzip -dc ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel1_CMO_R1.trimmed.fastq.gz \
  | awk 'NR % 4 == 2 {print length($0), $0; if (++n==20) exit}'
```

Expected length:

```text
28
```

---

## 14. Validate channel 2 trimmed R1

```bash
gzip -dc ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel2_CMO_R1.trimmed.fastq.gz \
  | awk 'NR % 4 == 2 {print length($0), $0; if (++n==20) exit}'
```

Expected length:

```text
28
```

---

## Link files to Processed_FASTQs for easier access in downstream steps

```bash
# GEX
ln -s ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel1_GEX_R1.fastq.gz \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Processed_FASTQs/channel1_GEX_R1.fastq.gz

ln -s ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel1_GEX_R2.fastq.gz \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Processed_FASTQs/channel1_GEX_R2.fastq.gz

ln -s ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel2_GEX_R1.fastq.gz \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Processed_FASTQs/channel2_GEX_R1.fastq.gz

ln -s ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel2_GEX_R2.fastq.gz \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Processed_FASTQs/channel2_GEX_R2.fastq.gz

# CMO
ln -s ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel1_CMO_R1.trimmed.fastq.gz \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Processed_FASTQs/channel1_CMO_R1.fastq.gz

ln -s ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel1_CMO_R2.fastq.gz \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Processed_FASTQs/channel1_CMO_R2.fastq.gz

ln -s ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel2_CMO_R1.trimmed.fastq.gz \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Processed_FASTQs/channel2_CMO_R1.fastq.gz

ln -s ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Splitcode/channel2_CMO_R2.fastq.gz \
  ${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Processed_FASTQs/channel2_CMO_R2.fastq.gz
```


## 15. KITE workflow

CMO quantification is automated via `cmo_quantification/quantify_cmo_tags.py`. The script downloads all inputs from IGVF (FASTQs, CMO barcodes, barcode onlist), derives the read format from the seqspec, builds the KITE index, and runs `kb count` for both channels. See [`cmo_quantification/README.md`](cmo_quantification/README.md) for full details.

| Input | IGVF accession |
|---|---|
| Channel 1 FASTQs | `IGVFDS2186TMQE` |
| Channel 2 FASTQs | `IGVFDS0370XXMJ` |
| CMO barcodes | `IGVFFI9308MNHO` |
| Cell barcode onlist | `IGVFFI8751YQRY` |

Submit the Slurm job:

```bash
sbatch analyses/10x_multi_5_timepoints/cmo_quantification/quantify_cmo_tags.sbatch
```

Outputs are written to `${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/CMO_counts/channel{1,2}/`.

## Conclusions

Final processed files for downstream analysis:

**GEX FASTQs** (split by splitcode, steps 6–7):
- `${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Processed_FASTQs/channel1_GEX_R1.fastq.gz`
- `${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Processed_FASTQs/channel1_GEX_R2.fastq.gz`
- `${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Processed_FASTQs/channel2_GEX_R1.fastq.gz`
- `${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/Processed_FASTQs/channel2_GEX_R2.fastq.gz`

**CMO counts** (step 15, automated):
- `${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/CMO_counts/channel1/counts_unfiltered/adata.h5ad`
- `${OAK}/Projects/EC_Screen/Data/10x_5_timepoints/CMO_counts/channel2/counts_unfiltered/adata.h5ad`