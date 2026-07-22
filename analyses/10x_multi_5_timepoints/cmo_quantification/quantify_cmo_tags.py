#!/usr/bin/env python3
"""
End-to-end CMO quantification from IGVF accessions.

Downloads CMO barcodes, FASTQs, and barcode onlist from IGVF via igvf-utils,
builds a KITE index, then runs kb count for channel1 and channel2.

Usage:
    python quantify_cmo_tags.py --work-dir $L_SCRATCH/cmo_quant

For Slurm, use the accompanying .sbatch file.
"""

import argparse
import logging
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "python"))

from utils import (
    build_index,
    find_seqspec,
    get_fastq_hrefs,
    get_read_format,
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

# IGVF accession defaults
DEFAULT_CH1_ACCESSION = "IGVFDS2186TMQE"
DEFAULT_CH2_ACCESSION = "IGVFDS0370XXMJ"
DEFAULT_CMO_ACCESSION = "IGVFFI9308MNHO"
DEFAULT_ONLIST_ACCESSION = "IGVFFI8751YQRY"

DEFAULT_OUTPUT_DIR = Path("/oak/stanford/groups/engreitz/Projects/EC_Screen/Data/10x_5_timepoints/CMO_counts")
DEFAULT_SEQSPEC_DIR = Path("/oak/stanford/groups/engreitz/Projects/EC_Screen/Data/10x_5_timepoints/seqspecs")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Quantify CMO tags for channel1 and channel2 from IGVF accessions."
    )
    parser.add_argument(
        "--ch1-accession", default=DEFAULT_CH1_ACCESSION,
        help=f"IGVF auxiliary set for channel 1. Default: {DEFAULT_CH1_ACCESSION}",
    )
    parser.add_argument(
        "--ch2-accession", default=DEFAULT_CH2_ACCESSION,
        help=f"IGVF auxiliary set for channel 2. Default: {DEFAULT_CH2_ACCESSION}",
    )
    parser.add_argument(
        "--cmo-accession", default=DEFAULT_CMO_ACCESSION,
        help=f"IGVF tabular file accession for CMO barcodes. Default: {DEFAULT_CMO_ACCESSION}",
    )
    parser.add_argument(
        "--onlist-accession", default=DEFAULT_ONLIST_ACCESSION,
        help=f"IGVF tabular file accession for the cell barcode whitelist. Default: {DEFAULT_ONLIST_ACCESSION}",
    )
    parser.add_argument(
        "--seqspec-dir", type=Path, default=DEFAULT_SEQSPEC_DIR,
        help=f"Directory containing per-channel seqspec YAML.gz files. Default: {DEFAULT_SEQSPEC_DIR}",
    )
    parser.add_argument(
        "--channels", nargs="+", default=["channel1", "channel2"],
        choices=["channel1", "channel2"],
        help="Channels to quantify. Default: both.",
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
    prepare_cmo_barcodes(args.cmo_accession, cmo_barcodes)

    # 2. KITE index
    index_dir = work / "kite_index"
    build_index(cmo_barcodes, index_dir, temp_dir=work / "tmp")

    # 3. Barcode onlist
    onlist_dest = work / f"{args.onlist_accession}.tsv.gz"
    onlist_meta = conn.get(f"tabular-files/{args.onlist_accession}", frame="object")
    stream_download(onlist_meta["href"], onlist_dest, conn)

    # 4. FASTQs + quantification per channel
    channel_accessions = {
        "channel1": args.ch1_accession,
        "channel2": args.ch2_accession,
    }
    for channel in args.channels:
        acc = channel_accessions[channel]

        seqspec_gz = find_seqspec(acc, args.seqspec_dir)
        read_format = get_read_format(seqspec_gz)
        logging.info("%s read format (from seqspec): %s", channel, read_format)

        r1_href, r2_href = get_fastq_hrefs(acc, conn)

        ch_dir = work / channel
        ch_dir.mkdir(exist_ok=True)
        r1 = ch_dir / r1_href.strip("/").split("/")[-1]
        r2 = ch_dir / r2_href.strip("/").split("/")[-1]

        stream_download(r1_href, r1, conn)
        stream_download(r2_href, r2, conn)

        quantify_channel(
            channel=channel,
            r1=r1,
            r2=r2,
            index_dir=index_dir,
            onlist=onlist_dest,
            read_format=read_format,
            output_dir=args.output_dir,
            threads=args.threads,
            memory=args.memory,
            temp_dir=work / "tmp",
        )


if __name__ == "__main__":
    main()
