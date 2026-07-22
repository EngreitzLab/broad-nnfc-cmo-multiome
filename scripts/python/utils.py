import gzip
import logging
import os
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import requests
from dotenv import find_dotenv, load_dotenv
from igvf_utils.connection import Connection

IGVF_API_BASE = "https://api.data.igvf.org"

# scripts/python/ -> scripts/ -> project root -> parent (git/) -> atomic-workflows
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
RUN_KITE = (
    _PROJECT_ROOT.parent
    / "atomic-workflows"
    / "modules"
    / "igvf-kite-cmo"
    / "run_kite.py"
)


# ---------------------------------------------------------------------------
# IGVF connection and download
# ---------------------------------------------------------------------------

def make_connection() -> Connection:
    """Load .env credentials and return an authenticated igvf-utils Connection."""
    load_dotenv(find_dotenv(usecwd=True))
    return Connection(igvf_mode="prod", no_log_file=True)


def stream_download(href: str, dest: Path, conn: Connection) -> None:
    """
    Stream a file from IGVF to dest using igvf-utils credentials.
    Skips if dest already exists. Uses 8 MB chunks for speed on large FASTQs.
    """
    if dest.exists():
        logging.info("Already downloaded: %s", dest)
        return
    url = f"{IGVF_API_BASE}{href}" if href.startswith("/") else href
    logging.info("Downloading %s -> %s", url, dest)
    tmp = dest.with_suffix(".tmp")
    with requests.get(url, auth=conn.auth, stream=True, allow_redirects=True) as r:
        r.raise_for_status()
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(chunk_size=8 << 20):
                fh.write(chunk)
    tmp.rename(dest)
    logging.info("Saved: %s (%.1f MB)", dest, dest.stat().st_size / 1e6)


def get_fastq_hrefs(aux_accession: str, conn: Connection) -> tuple[str, str]:
    """Return (href_R1, href_R2) for the sequence files in an IGVF auxiliary set."""
    hrefs = get_fastq_hrefs_ordered(aux_accession, conn, ["R1", "R2"])
    return hrefs[0], hrefs[1]


def get_fastq_hrefs_ordered(
    aux_accession: str, conn: Connection, read_types: list[str]
) -> list[str]:
    """Return hrefs from an auxiliary set ordered by read_types.

    Use when kb count needs more than two reads (e.g. R1, R2, R3).
    """
    meta = conn.get(f"auxiliary-sets/{aux_accession}", frame="object")
    if not meta:
        raise ValueError(f"Auxiliary set not found: {aux_accession}")
    seq_accessions = [
        p.strip("/").split("/")[-1]
        for p in meta.get("files", [])
        if "sequence-files" in p
    ]
    hrefs: dict[str, str] = {}
    for acc in seq_accessions:
        fm = conn.get(f"sequence-files/{acc}", frame="object")
        rt = fm.get("illumina_read_type", "")
        if rt in read_types:
            hrefs[rt] = fm["href"]
    missing = [rt for rt in read_types if rt not in hrefs]
    if missing:
        raise ValueError(
            f"Could not find {missing} in auxiliary set {aux_accession}. Found: {list(hrefs)}"
        )
    return [hrefs[rt] for rt in read_types]


# ---------------------------------------------------------------------------
# seqspec
# ---------------------------------------------------------------------------

def download_seqspec(accession: str, dest: Path, conn: Connection) -> None:
    """Download a seqspec YAML.gz from an IGVF configuration-file accession."""
    if dest.exists():
        logging.info("Already downloaded seqspec: %s", dest)
        return
    meta = conn.get(f"configuration-files/{accession}", frame="object")
    stream_download(meta["href"], dest, conn)


def find_seqspec(accession: str, seqspec_dir: Path) -> Path:
    """Find the seqspec YAML.gz for a given auxiliary set accession."""
    matches = list(seqspec_dir.glob(f"*{accession}*.yaml.gz"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected exactly 1 seqspec matching '{accession}' in {seqspec_dir}, "
            f"found: {matches}"
        )
    return matches[0]


def get_read_format(seqspec_gz: Path, modality: str = "tag") -> str:
    """
    Derive the kb read format string (e.g. '0,0,16:0,16,28:1,0,15') from a
    seqspec YAML.gz via: seqspec index -m <modality> -s file -t kb
    """
    with tempfile.NamedTemporaryFile(suffix=".yaml", delete=False) as tmp:
        tmp_path = Path(tmp.name)
        tmp.write(gzip.decompress(seqspec_gz.read_bytes()))
    try:
        result = subprocess.run(
            ["seqspec", "index", "-m", modality, "-s", "file", "-t", "kb", str(tmp_path)],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()
    finally:
        tmp_path.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Shell execution
# ---------------------------------------------------------------------------

def run_shell_cmd(cmd: str) -> str:
    """
    Run a bash command in a child process group.
    Logs PID, PGID, duration, stdout, and stderr.
    Kills the process group and raises on non-zero exit code.
    Returns stdout (newlines stripped).
    Adapted from https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/26eeda81a0540dc793fc69b0c390d232ca7ca50a/src/encode_lib_common.py#L331
    """
    p = subprocess.Popen(
        ["/bin/bash", "-o", "pipefail"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=os.setsid,
    )
    pid = p.pid
    pgid = os.getpgid(pid)
    logging.info("run_shell_cmd: PID=%d, PGID=%d, CMD=%s", pid, pgid, cmd)
    t0 = time.monotonic()
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    dur = time.monotonic() - t0
    err_str = (
        "PID={pid}, PGID={pgid}, RC={rc}, DURATION_SEC={dur:.1f}\n"
        "STDERR={stde}\nSTDOUT={stdo}"
    ).format(pid=pid, pgid=pgid, rc=rc, dur=dur, stde=stderr.strip(), stdo=stdout.strip())
    if rc:
        try:
            os.killpg(pgid, signal.SIGKILL)
        except Exception:
            pass
        finally:
            raise Exception(err_str)
    else:
        logging.info(err_str)
    return stdout.strip("\n")


# ---------------------------------------------------------------------------
# KITE pipeline steps (delegated to the igvf-kite-cmo atomic workflow)
# ---------------------------------------------------------------------------

def prepare_cmo_barcodes(
    accession: str,
    dest: Path,
    sequence_col: str = "multiseq_bc",
    name_col: str = "CMO ID",
) -> None:
    """
    Download a CMO barcode CSV from IGVF and write a KITE-format gzipped TSV.
    Delegates to: run_kite prepare-barcodes --accession <accession> --output <dest>
    Auth is read from .env by run_kite.
    """
    if dest.exists():
        logging.info("CMO barcodes already prepared: %s", dest)
        return
    run_shell_cmd(
        f"{sys.executable} {RUN_KITE} prepare-barcodes "
        f"--accession {accession} "
        f"--sequence-col '{sequence_col}' "
        f"--name-col '{name_col}' "
        f"--output {dest}"
    )


def build_index(cmo_barcodes: Path, index_dir: Path, temp_dir: Path | None) -> None:
    """
    Build the KITE index; skips if already built.
    Delegates to: run_kite index --output_dir <index_dir> --cmo_barcodes <cmo_barcodes>
    """
    if (index_dir / "cmo.idx").exists():
        logging.info("KITE index already exists: %s", index_dir)
        return
    index_dir.mkdir(parents=True, exist_ok=True)
    tmp = f"--temp_dir {temp_dir}" if temp_dir else ""
    run_shell_cmd(
        f"{sys.executable} {RUN_KITE} index "
        f"--output_dir {index_dir} "
        f"--cmo_barcodes {cmo_barcodes} "
        f"{tmp}"
    )


def quantify_channel(
    channel: str,
    fastqs: list[Path],
    index_dir: Path,
    onlist: Path,
    read_format: str,
    output_dir: Path,
    threads: int,
    memory: str,
    temp_dir: Path | None,
) -> None:
    """
    Quantify CMO tags for one channel; skips if h5ad output already exists.
    Delegates to: run_kite quantify --index_dir ... --output_dir ... <fastqs...>
    fastqs are passed in the order expected by the kb read format string.
    """
    out = output_dir / channel
    if (out / "counts_unfiltered" / "adata.h5ad").exists():
        logging.info("Already done: %s", channel)
        return
    out.mkdir(parents=True, exist_ok=True)
    tmp = f"--temp_dir {temp_dir}" if temp_dir else ""
    fastq_args = " ".join(str(f) for f in fastqs)
    run_shell_cmd(
        f"{sys.executable} {RUN_KITE} quantify "
        f"--index_dir {index_dir} "
        f"--output_dir {out} "
        f"--barcode_onlist {onlist} "
        f"--read_format {read_format} "
        f"--threads {threads} "
        f"--memory {memory} "
        f"{tmp} "
        f"{fastq_args}"
    )
    logging.info("Done: %s -> %s", channel, out)
