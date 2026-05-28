BASE="${HOME}/GitHub/broad-nnfc-cmo-multiome"
LANE1="${BASE}/data/channel2"

R1=$(ls ${LANE1}/*GEX_2*_R1_*.fastq.gz)
R2=$(ls ${LANE1}/*GEX_2*_R2_*.fastq.gz)

splitcode \
  -c "${BASE}/config/cmo_25_splitcode.config" \
  --nFastqs=2 \
  --gzip \
  --assign \
  --mapping lane2_cmo_mapping.txt \
  --summary lane2_cmo_summary.txt \
  -o lane2_CMO_R1.fastq.gz,lane2_CMO_R2.fastq.gz \
  -u lane2_notCMO_R1.fastq.gz,lane2_notCMO_R2.fastq.gz \
  $R1 $R2

