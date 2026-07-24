#!/usr/bin/env python3
"""
End-to-end CMO quantification for the 10x multiome 5 timepoints McGinnis dataset.

Single channel. CMO library structure (from seqspec IGVFFI6462KGWB):
  R1 (IGVFFI5411OJRK): cell barcode (0-16bp) + UMI (16-28bp)
  R2 (IGVFFI4642YKWI): constant capture/handle sequence — NOT the tag, unused
  R3 (IGVFFI5223HYRP): MULTI-seq barcode / CMO tag (0-8bp) + polyA
  I1 (IGVFFI6088QALC): i7 index — not used for CMO quantification

NOTE: the CMO tag lives in R3, not R2. R2 (IGVFFI4642YKWI) is a near-constant
handle (every read ~ TG.TCTCGGTGGTCGCCGTATCAT); using it as the tag gave only
1.3% pseudoalignment and collapsed all calls onto one tag. R3's first 8bp match
the 23-tag panel ~92% exactly (fwd, pos 0-8). This agrees with seqspec
IGVFFI6462KGWB, whose kb output places the 8bp tag in file index 2 (= R3).

CMO barcode columns: 'barcode' (sequence) and 'sample description' (name).

Usage:
    python quantify_cmo_tags.py --work-dir $L_SCRATCH/cmo_quant

For Slurm, use the accompanying .sbatch file.
"""

import argparse
import gzip
import logging
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "python"))

from utils import (
    build_index,
    get_fastq_hrefs_ordered,
    make_connection,
    prepare_cmo_barcodes,
    quantify_channel,
    stream_download,
)

logging.basicConfig(
    stream=sys.stderr,
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

# IGVF accessions
DEFAULT_AUX_ACCESSION = "IGVFDS4882ESVG"    # CMO auxiliary set
DEFAULT_CMO_ACCESSION = "IGVFFI5955PKRW"    # CMO barcode-to-sample map
DEFAULT_ONLIST_ACCESSION = "IGVFFI8751YQRY" # 10x multiome cell barcode whitelist

# CMO barcode file column names (differ from the 5 timepoints dataset)
CMO_SEQUENCE_COL = "barcode"
CMO_NAME_COL = "sample description"

# kb read format: file0=R1 (barcode 0-16, UMI 16-28), file1=R3 (CMO tag 0-8).
# We pass fastqs as [R1, R3], so file1 is R3 and the tag range is 1,0,8. seqspec
# IGVFFI6462KGWB places the 8bp tag in file index 2 (R3): its kb output segment
# '...,2,0,8:...' was the correct clue; the earlier hardcoding to R2 was the bug.
READ_FORMAT = "0,0,16:0,16,28:1,0,8"

DEFAULT_OUTPUT_DIR = Path("/oak/stanford/groups/engreitz/Projects/EC_Screen/Data/10x_5_timepoints_mcginnis/CMO_counts")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Quantify CMO tags for the McGinnis 10x multiome 5 timepoints dataset."
    )
    parser.add_argument(
        "--aux-accession", default=DEFAULT_AUX_ACCESSION,
        help=f"IGVF auxiliary set with CMO FASTQs. Default: {DEFAULT_AUX_ACCESSION}",
    )
    parser.add_argument(
        "--cmo-accession", default=DEFAULT_CMO_ACCESSION,
        help=f"IGVF tabular file for CMO barcodes. Default: {DEFAULT_CMO_ACCESSION}",
    )
    parser.add_argument(
        "--onlist-accession", default=DEFAULT_ONLIST_ACCESSION,
        help=f"IGVF tabular file for cell barcode whitelist. Default: {DEFAULT_ONLIST_ACCESSION}",
    )
    parser.add_argument(
        "--work-dir", type=Path, required=True,
        help="Working directory for downloads and index (use $L_SCRATCH on Slurm).",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
        help=f"Results directory. Default: {DEFAULT_OUTPUT_DIR}",
    )
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--memory", type=str, default="4G")
    args = parser.parse_args()

    work = args.work_dir
    work.mkdir(parents=True, exist_ok=True)

    conn = make_connection()

    # 1. CMO barcodes -> KITE TSV
    cmo_barcodes = work / "cmo_barcodes_for_kite.tsv.gz"
    prepare_cmo_barcodes(
        args.cmo_accession, cmo_barcodes, conn,
        sequence_col=CMO_SEQUENCE_COL,
        name_col=CMO_NAME_COL,
    )

    # 2. KITE index
    index_dir = work / "kite_index"
    build_index(cmo_barcodes, index_dir, temp_dir=work / "tmp")

    # 3. Barcode onlist (bustools correct requires uncompressed)
    onlist_gz = work / f"{args.onlist_accession}.tsv.gz"
    onlist_meta = conn.get(f"tabular-files/{args.onlist_accession}", frame="object")
    stream_download(onlist_meta["href"], onlist_gz, conn)
    onlist_dest = onlist_gz.with_suffix("")  # .tsv
    if not onlist_dest.exists():
        logging.info("Decompressing onlist -> %s", onlist_dest)
        onlist_dest.write_bytes(gzip.decompress(onlist_gz.read_bytes()))

    # 4. FASTQs: R1 (barcode+UMI) and R3 (CMO tag) — file0 and file1
    r1_href, r3_href = get_fastq_hrefs_ordered(
        args.aux_accession, conn, ["R1", "R3"]
    )
    fastq_dir = work / "fastqs"
    fastq_dir.mkdir(exist_ok=True)
    r1 = fastq_dir / r1_href.strip("/").split("/")[-1]
    r3 = fastq_dir / r3_href.strip("/").split("/")[-1]
    stream_download(r1_href, r1, conn)
    stream_download(r3_href, r3, conn)

    # 5. Quantify
    quantify_channel(
        channel="cmo",
        fastqs=[r1, r3],
        index_dir=index_dir,
        onlist=onlist_dest,
        read_format=READ_FORMAT,
        output_dir=args.output_dir,
        threads=args.threads,
        memory=args.memory,
        temp_dir=work / "tmp",
    )


if __name__ == "__main__":
    main()
