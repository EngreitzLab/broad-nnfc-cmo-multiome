#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 2 || $# -gt 4 ]]; then
  cat <<'EOF'
Usage: benchmark_cb_anchor.bash <cmo_r1.fastq.gz> <barcode_whitelist.txt[.gz]> [n_reads] [window_half_width]

Compares three barcode detection strategies on CMO R1 reads:
1. exact tail anchor: barcode is the first 16 bp of the last 28 bp
2. narrow window: search for the barcode near the tail-anchor position
3. full scan: search the whole read for the first whitelist barcode match

Outputs recall of the first two methods relative to the full scan.
EOF
  exit 1
fi

R1_FASTQ=$1
WHITELIST=$2
N_READS=${3:-100000}
WINDOW_HALF_WIDTH=${4:-4}

if [[ ! -f "$R1_FASTQ" ]]; then
  echo "Missing FASTQ: $R1_FASTQ" >&2
  exit 1
fi

if [[ ! -f "$WHITELIST" ]]; then
  echo "Missing whitelist: $WHITELIST" >&2
  exit 1
fi

if [[ "$WHITELIST" == *.gz ]]; then
  WL_CMD=(gzip -dc "$WHITELIST")
else
  WL_CMD=(cat "$WHITELIST")
fi

gzip -dc "$R1_FASTQ" | awk \
  -v max_reads="$N_READS" \
  -v window_half_width="$WINDOW_HALF_WIDTH" \
  -v wl_cmd="${WL_CMD[*]}" '
BEGIN {
  while ((wl_cmd | getline line) > 0) {
    gsub(/\r/, "", line)
    if (line == "") {
      continue
    }
    sub(/-1$/, "", line)
    whitelist[line] = 1
    whitelist_count++
  }
  close(wl_cmd)
}

NR % 4 != 2 {
  next
}

{
  seq = $0
  len = length(seq)
  max_start = len - 15
  anchor_start = len - 27

  total_reads++
  read_length[len]++

  full_found = 0
  full_pos = 0
  for (i = 1; i <= max_start; i++) {
    candidate = substr(seq, i, 16)
    if (candidate in whitelist) {
      full_found = 1
      full_pos = i
      full_hits++
      full_hit_pos[i]++
      break
    }
  }

  if (!full_found) {
    full_misses++
  }

  tail_found = 0
  if (anchor_start >= 1) {
    tail_candidate = substr(seq, anchor_start, 16)
    if (tail_candidate in whitelist) {
      tail_found = 1
      tail_hits++
    }
  }

  narrow_found = 0
  narrow_start = anchor_start - window_half_width
  narrow_end = anchor_start + window_half_width

  if (narrow_start < 1) {
    narrow_start = 1
  }
  if (narrow_end > max_start) {
    narrow_end = max_start
  }

  for (i = narrow_start; i <= narrow_end; i++) {
    candidate = substr(seq, i, 16)
    if (candidate in whitelist) {
      narrow_found = 1
      narrow_hits++
      narrow_hit_pos[i]++
      break
    }
  }

  if (full_found) {
    if (!tail_found) {
      tail_missed_vs_full++
      tail_miss_pos[full_pos]++
    }
    if (!narrow_found) {
      narrow_missed_vs_full++
      narrow_miss_pos[full_pos]++
    }
    anchor_delta[full_pos - anchor_start]++
  }

  if (total_reads >= max_reads) {
    exit
  }
}

function print_histogram(title, values,    key) {
  print title
  for (key in values) {
    print key "\t" values[key]
  }
  print ""
}

END {
  print "reads_tested\t" total_reads
  print "whitelist_entries\t" whitelist_count
  print "full_scan_hits\t" full_hits
  print "full_scan_misses\t" full_misses
  print "tail28_exact_hits\t" tail_hits
  print "narrow_window_hits\t" narrow_hits

  if (full_hits > 0) {
    print "tail28_misses_vs_full\t" (full_hits - tail_hits)
    print "narrow_window_misses_vs_full\t" (full_hits - narrow_hits)
    print "narrow_window_reads_rescued_vs_tail28\t" (narrow_hits - tail_hits)
    printf("tail28_recall_vs_full\t%.6f\n", tail_hits / full_hits)
    printf("narrow_window_recall_vs_full\t%.6f\n", narrow_hits / full_hits)
  }

  print ""
  print_histogram("read_length_histogram", read_length)
  print_histogram("full_scan_hit_positions", full_hit_pos)
  print_histogram("narrow_window_hit_positions", narrow_hit_pos)
  print_histogram("anchor_delta_full_minus_tail_start", anchor_delta)
  print_histogram("tail28_miss_positions_from_full", tail_miss_pos)
  print_histogram("narrow_window_miss_positions_from_full", narrow_miss_pos)
}'