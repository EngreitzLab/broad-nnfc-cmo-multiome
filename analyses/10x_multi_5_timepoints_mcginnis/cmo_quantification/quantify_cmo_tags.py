#!/usr/bin/env python3
"""
End-to-end CMO quantification for the 10x multiome 5 timepoints McGinnis dataset.

Single channel, three-read CMO library (R1=cell barcode+UMI, R2=CMO tag, R3=cDNA).
Seqspec is fetched from IGVF (IGVFFI6462KGWB) and used to derive the kb read format.
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
    download_seqspec,
    get_fastq_hrefs_ordered,
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

# IGVF accessions
DEFAULT_AUX_ACCESSION = "IGVFDS4882ESVG"       # CMO auxiliary set (R1, R2, R3, I1)
DEFAULT_CMO_ACCESSION = "IGVFFI5955PKRW"        # CMO barcode-to-sample map
DEFAULT_ONLIST_ACCESSION = "IGVFFI8751YQRY"     # 10x multiome cell barcode whitelist
DEFAULT_SEQSPEC_ACCESSION = "IGVFFI6462KGWB"    # seqspec for CMO library

# CMO barcode file uses different column names than the 5 timepoints dataset
CMO_SEQUENCE_COL = "barcode"
CMO_NAME_COL = "sample description"

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
        "--seqspec-accession", default=DEFAULT_SEQSPEC_ACCESSION,
        help=f"IGVF configuration file for the CMO seqspec. Default: {DEFAULT_SEQSPEC_ACCESSION}",
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
        args.cmo_accession, cmo_barcodes,
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

    # 4. Seqspec -> kb read format (seqspec lives on IGVF for this dataset)
    seqspec_gz = work / f"{args.seqspec_accession}.yaml.gz"
    download_seqspec(args.seqspec_accession, seqspec_gz, conn)
    read_format = get_read_format(seqspec_gz)
    logging.info("Read format (from seqspec): %s", read_format)

    # 5. FASTQs: R1 (barcode+UMI), R2 (CMO tag), R3 (cDNA) — order matters for kb
    r1_href, r2_href, r3_href = get_fastq_hrefs_ordered(
        args.aux_accession, conn, ["R1", "R2", "R3"]
    )
    fastq_dir = work / "fastqs"
    fastq_dir.mkdir(exist_ok=True)
    r1 = fastq_dir / r1_href.strip("/").split("/")[-1]
    r2 = fastq_dir / r2_href.strip("/").split("/")[-1]
    r3 = fastq_dir / r3_href.strip("/").split("/")[-1]
    stream_download(r1_href, r1, conn)
    stream_download(r2_href, r2, conn)
    stream_download(r3_href, r3, conn)

    # 6. Quantify — pass R1, R2, R3 in order; the seqspec format string
    #    references file indices 0 (R1) and 2 (R3) for barcode/UMI and
    #    file 0 again for the feature, so all three must be present.
    quantify_channel(
        channel="cmo",
        fastqs=[r1, r2, r3],
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
