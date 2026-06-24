# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "igvf-utils==3.1.1",
#     "jinja2==3.1.6",
#     "marimo>=0.23.10",
#     "python-dotenv==1.2.2",
#     "pyyaml==6.0.3",
#     "requests==2.34.2",
#     "seqspec==0.3.1",
# ]
# ///

import marimo

__generated_with = "0.23.10"
app = marimo.App(
    css_file="/usr/local/_marimo/custom.css",
    auto_download=["html"],
)


@app.cell
def _():
    import marimo as mo

    return


@app.cell
def _():
    import csv
    import gzip
    import hashlib
    import io
    import json
    import os
    from pathlib import Path
    from dotenv import load_dotenv
    from igvf_utils.connection import Connection

    return Connection, Path, gzip, hashlib, io, json, load_dotenv, os


@app.cell
def _(Connection, load_dotenv, require_env):
    # This notebook expect a .env file in the repository root folder with the igvf credentials
    #IGVF_API_KEY=
    #IGVF_SECRET_KEY=

    # Connecting to the IGVF portal
    load_dotenv(override=True)

    require_env("IGVF_API_KEY")
    require_env("IGVF_SECRET_KEY")

    domain = "api.data.igvf.org"
    api_base = f"https://{domain}"
    conn = Connection(domain)
    return (conn,)


@app.cell
def _():
    # You can change this to be your IGVF intermediate analyses accessions
    # For 10X multiome 5 timepoints:
    # - IGVFDS5477BPOI
    # - IGVFDS3995WHFT
    intermediate_analysis_sets = [
        "IGVFDS5477BPOI",
        "IGVFDS3995WHFT"
    ]
    return (intermediate_analysis_sets,)


@app.cell
def _(conn, intermediate_analysis_sets):
    # Fetch both analysis sets and collect their input_file_sets
    analysis_sets_raw = {acc: conn.get(f"/analysis-sets/{acc}/") for acc in intermediate_analysis_sets}

    # Map channel label -> list of input file sets (measurement + auxiliary + curated)
    input_file_sets_by_channel = {}
    for _acc, _md in analysis_sets_raw.items():
        _alias = (_md.get("aliases") or [_acc])[0]
        _ch_label = _alias.split("_ch")[-1] if "_ch" in _alias else _acc
        input_file_sets_by_channel[f"ch{_ch_label}_{_acc}"] = _md.get("input_file_sets", [])

    for _ch, _fsets in input_file_sets_by_channel.items():
        print(f"\n=== {_ch} ===")
        for _fs in _fsets:
            print(f"  {_fs['accession']}  type={_fs.get('file_set_type')}  assay={_fs.get('assay_term', {}).get('term_name', 'N/A')}  alias={_fs.get('aliases', ['?'])[0]}")
    return (input_file_sets_by_channel,)


@app.cell
def _(conn, input_file_sets_by_channel):
    # For each measurement/auxiliary set, fetch its sequence files and key metadata
    FIELDS_OF_INTEREST = ["accession", "file_size", "md5sum", "href", "submitted_file_name",
                          "sequencing_run", "read_type", "content_type", "sequencing_platform",
                          "sequencing_kit", "illumina_read_type", "file_format", "aliases"]

    def _fetch_file_set(accession, file_set_type):
        if file_set_type in ("experimental data",):
            return conn.get(f"/measurement-sets/{accession}/")
        elif file_set_type == "lipid-conjugated oligo sequencing":
            return conn.get(f"/auxiliary-sets/{accession}/")
        else:
            return conn.get(f"/curated-sets/{accession}/")

    def _platform_str(platform_obj):
        """Convert portal sequencing_platform object to seqspec enum string."""
        _term = platform_obj.get("term_name", "")
        _raw_id = platform_obj.get("@id", "").strip("/").split("/")[-1]
        _efo_id = _raw_id.replace("_", ":", 1)
        return f"{_term} ({_efo_id})"

    fastq_metadata = {}  # keyed by (channel_key, set_accession)

    for _ch_key, _fsets in input_file_sets_by_channel.items():
        for _fs in _fsets:
            _acc = _fs["accession"]
            _ftype = _fs.get("file_set_type", "")
            if _ftype == "barcodes":
                continue
            _set_md = _fetch_file_set(_acc, _ftype)
            _files = _set_md.get("files", [])
            _file_details = []
            for _f in _files:
                _f_acc = _f if isinstance(_f, str) else _f.get("accession") or _f.get("@id", "")
                if not _f_acc:
                    continue
                if "/" in _f_acc:
                    _f_acc = _f_acc.strip("/").split("/")[-2]
                try:
                    _fmd = conn.get(f"/sequence-files/{_f_acc}/")
                    _file_details.append({k: _fmd.get(k) for k in FIELDS_OF_INTEREST})
                except Exception as _e:
                    _file_details.append({"accession": _f_acc, "error": str(_e)})

            # Derive sequencing protocol/kit from the first file with that info
            _seq_protocol = ""
            _seq_kit = ""
            for _fd in _file_details:
                _plat = _fd.get("sequencing_platform")
                if _plat:
                    _seq_protocol = _platform_str(_plat)
                    _seq_kit = _fd.get("sequencing_kit", "")
                    break

            fastq_metadata[(_ch_key, _acc)] = {
                "alias":             (_fs.get("aliases") or [_acc])[0],
                "assay":             _fs.get("assay_term", {}).get("term_name", _ftype),
                "lab":               _set_md.get("lab", {}).get("@id", ""),
                "award":             _set_md.get("award", {}).get("@id", ""),
                "sequencing_protocol": _seq_protocol,
                "sequencing_kit":    _seq_kit,
                "files":             _file_details,
            }

    # Summary
    for (_ch, _set_acc), _info in fastq_metadata.items():
        print(f"\n{_ch} | {_set_acc} | {_info['assay']}")
        print(f"  lab={_info['lab']}  award={_info['award']}")
        print(f"  seq_protocol={_info['sequencing_protocol']}  kit={_info['sequencing_kit']!r}")
        for _f in _info["files"]:
            print(f"  {_f.get('accession')}  read={_f.get('illumina_read_type')}  run={_f.get('sequencing_run')}")

    return (fastq_metadata,)


@app.cell
def _(conn, fastq_metadata):
    # Summarise available files per channel/modality and inspect download URL
    from collections import defaultdict

    summary = defaultdict(lambda: defaultdict(list))

    for (_ch, _set_acc), _info in fastq_metadata.items():
        _channel = _ch.split("_")[0]   # ch1 or ch2
        _modality = _info["assay"]
        for _f in _info["files"]:
            summary[_channel][_modality].append({
                "accession": _f.get("accession"),
                "read":      _f.get("illumina_read_type"),
                "run":       _f.get("sequencing_run"),
                "size":      _f.get("file_size"),
                "md5":       _f.get("md5sum"),
                "href":      _f.get("href"),
            })

    for _channel in sorted(summary):
        for _modality, _files in summary[_channel].items():
            _reads = sorted(set(_f["read"] for _f in _files))
            print(f"{_channel} | {_modality}: {len(_files)} files, reads={_reads}")

    # Spot-check the first available file's download URL
    _sample_info = next(iter(fastq_metadata.values()))
    _sample_file = _sample_info["files"][0]
    _sample_acc = _sample_file["accession"]
    _sample_md = conn.get(f"/sequence-files/{_sample_acc}/")
    print(f"\nSample file {_sample_acc}:")
    print(f"  href            : {_sample_md.get('href')}")
    print(f"  submitted_file_name: {_sample_md.get('submitted_file_name')}")
    print(f"  cloud_metadata  : {_sample_md.get('cloud_metadata')}")

    return


@app.cell
def _(Path, conn, fastq_metadata, gzip, input_file_sets_by_channel, io):
    from jinja2 import Environment, FileSystemLoader

    _TEMPLATE_DIR = Path("templates")
    _jinja_env = Environment(loader=FileSystemLoader(str(_TEMPLATE_DIR)), keep_trailing_newline=True)

    TEMPLATE_BY_ASSAY = {
        "single-nucleus RNA sequencing assay": "10x_multiome_gex_seqspec.yaml.j2",
        "single-nucleus ATAC-seq":             "10x_multiome_atac_seqspec.yaml.j2",
        "lipid-conjugated oligo sequencing":   "10x_multiome_cmo_seqspec.yaml.j2",
    }

    LIBRARY_PROTOCOL_BY_ASSAY = {
        "single-nucleus RNA sequencing assay": "single-nucleus RNA sequencing assay (OBI:0003109)",
        "single-nucleus ATAC-seq":             "single-nucleus ATAC-seq (OBI:0002762)",
        "lipid-conjugated oligo sequencing":   "Custom",
    }

    LIB_STRUCT = "https://teichlab.github.io/scg_lib_structs/methods_html/10xChromiumMultiome.html"
    LIB_KIT    = "Chromium Next GEM Multiome ATAC + Gene Expression"


    def _tabular_ctx(acc, prefix):
        _md = conn.get(f"/tabular-files/{acc}/")
        _url = f"https://api.data.igvf.org{_md['href']}"
        return {
            f"{prefix}_file_id":      _md["accession"],
            f"{prefix}_file_name":    _url,
            f"{prefix}_file_type":    _md["file_format"],
            f"{prefix}_file_size":    _md["file_size"],
            f"{prefix}_file_url":     _url,
            f"{prefix}_file_urltype": "https",
            f"{prefix}_file_md5":     _md["md5sum"],
        }

    def _seq_file_ctx(f):
        _url = f"https://api.data.igvf.org{f['href']}"
        return dict(
            file_id=f["accession"],
            filename=f["accession"] + ".fastq.gz",
            filetype="fastq",
            filesize=f["file_size"],
            url=_url,
            urltype="https",
            md5=f["md5sum"],
        )

    # GEX and ATAC: standard 10x Multiome reference files registered with stable aliases
    _gex_bc_ctx  = _tabular_ctx("igvf:10X_Multiome_GEX_cell_barcode_inclusion_list", "cell_barcode")
    _atac_bc_ctx = _tabular_ctx("igvf:10X_Multiome_ATAC_cell_barcode_inclusion_list", "cell_barcode")

    # CMO: discovered from the barcodes curated set in input_file_sets_by_channel
    _cmo_bc_acc = None
    for _bc_ch, _bc_fsets in input_file_sets_by_channel.items():
        for _bc_fs in _bc_fsets:
            if _bc_fs.get("file_set_type") == "barcodes":
                _bc_rec = conn.get(f"/curated-sets/{_bc_fs['accession']}/")
                for _bf in _bc_rec.get("files", []):
                    if isinstance(_bf, dict) and _bf.get("content_type") == "barcode to sample mapping":
                        _cmo_bc_acc = _bf["accession"]
                        break

    if not _cmo_bc_acc:
        raise RuntimeError("Could not find CMO barcode file in the barcodes curated set")
    _cmo_bc_ctx = _tabular_ctx(_cmo_bc_acc, "cmo_barcode")
    print(f"GEX barcode: {_gex_bc_ctx['cell_barcode_file_id']}")
    print(f"ATAC barcode: {_atac_bc_ctx['cell_barcode_file_id']}")
    print(f"CMO barcode: {_cmo_bc_ctx['cmo_barcode_file_id']}")

    _seqspec_dir = Path("seqspecs")
    _seqspec_dir.mkdir(exist_ok=True)

    seqspec_paths = {}

    for (_ch_key, _set_acc), _info in fastq_metadata.items():
        _assay = _info["assay"]
        _tmpl_name = TEMPLATE_BY_ASSAY.get(_assay)
        if not _tmpl_name:
            print(f"WARNING: no template for assay '{_assay}', skipping {_set_acc}")
            continue

        _channel = _ch_key.split("_")[0]  # "ch1", "ch2", etc.

        _by_read = {}
        for _f in _info["files"]:
            _rt = _f.get("illumina_read_type")
            if _rt:
                _by_read.setdefault(_rt, []).append(_seq_file_ctx(_f))

        _ctx = {
            "library_structure":   LIB_STRUCT,
            "library_protocol":    LIBRARY_PROTOCOL_BY_ASSAY[_assay],
            "library_kit":         LIB_KIT,
            "sequencing_protocol": _info["sequencing_protocol"],
            "sequencing_kit":      _info["sequencing_kit"],
            "read1_id":   "R1", "read1_name": "Read 1",
            "read2_id":   "R2", "read2_name": "Read 2",
            "read1_files": _by_read.get("R1", []),
            "read2_files": _by_read.get("R2", []),
        }

        if _assay == "single-nucleus RNA sequencing assay":
            _ctx["assay_id"] = f"{_channel}_gex_{_set_acc}"
            _ctx.update(_gex_bc_ctx)

        elif _assay == "single-nucleus ATAC-seq":
            _ctx["assay_id"] = f"{_channel}_atac_{_set_acc}"
            _ctx.update(_atac_bc_ctx)
            _ctx.update({
                "read3_id": "R3", "read3_name": "Read 3",
                "read_i1_id": "I1", "read_i1_name": "Index 1",
                "read3_files":   _by_read.get("R3", []),
                "read_i1_files": _by_read.get("I1", []),
            })

        elif _assay == "lipid-conjugated oligo sequencing":
            _ctx["assay_id"] = f"{_channel}_cmo_{_set_acc}"
            _ctx.update(_gex_bc_ctx)
            _ctx.update(_cmo_bc_ctx)

        _yaml_text = _jinja_env.get_template(_tmpl_name).render(**_ctx)
        _gz_buf = io.BytesIO()
        with gzip.GzipFile(fileobj=_gz_buf, mode="wb", mtime=0) as _gz:
            _gz.write(_yaml_text.encode())
        _gz_path = _seqspec_dir / f"{_set_acc}.yaml.gz"
        _gz_path.write_bytes(_gz_buf.getvalue())
        seqspec_paths[_ctx["assay_id"]] = _gz_path
        print(f"written: {_gz_path}")

    seqspec_paths

    return (seqspec_paths,)


@app.cell
def _(os):
    def require_env(var_name: str) -> str:
        value = os.getenv(var_name)
        if not value:
            raise RuntimeError(
                f"Missing environment variable '{var_name}'. "
                "Set it in your shell before running this script."
            )
        return value

    return (require_env,)


@app.cell
def _(Path, fastq_metadata, hashlib, json, seqspec_paths):
    SEQSPEC_DIR = Path("seqspecs")

    def _config_alias(file_set_alias):
        for prefix in ("measurement_", "auxiliaryset_", "curated_set_"):
            if prefix in file_set_alias:
                lab, rest = file_set_alias.split(":", 1)
                return f"{lab}:configuration_file_{rest.replace(prefix, '', 1)}"
        return file_set_alias

    _rows = []
    for _assay_id, _gz_path in seqspec_paths.items():
        _set_acc = _assay_id.rsplit("_", 1)[-1]
        _info = next(v for (_, acc), v in fastq_metadata.items() if acc == _set_acc)
        _md5 = hashlib.md5(_gz_path.read_bytes()).hexdigest()
        _file_accs = [_f["accession"] for _f in _info["files"] if _f.get("accession")]
        _rows.append({
            "aliases":             json.dumps([_config_alias(_info["alias"])]),
            "award":               _info["award"],
            "lab":                 _info["lab"],
            "md5sum":              _md5,
            "file_format":         "yaml",
            "file_set":            _set_acc,
            "content_type":        "seqspec",
            "submitted_file_name": _gz_path.name,
            "seqspec_of":          json.dumps(_file_accs),
        })

    manifest_tsv = SEQSPEC_DIR / "seqspec_yamls_manifest.tsv"
    _cols = list(_rows[0].keys())
    with open(manifest_tsv, "w") as _f:
        _f.write("\t".join(_cols) + "\n")
        for _row in _rows:
            _f.write("\t".join(_row[c] for c in _cols) + "\n")

    manifest_tsv

    return


if __name__ == "__main__":
    app.run()
