"""Synthetic fixtures and a fake download sandbox for the M10-01 tests (PXD078459 download).

Nothing here reaches the network, the real ``$GIPPA_DATA_DIR``/``$GIPPA_ROOT`` or real MS data:

* ``px_xml`` writes a small ProteomeXchange record (PX XML 1.4 layout) from a list of file specs;
* ``make_tdf``/``make_run`` build a minimal Bruker ``analysis.tdf`` (SQLite) with the tables the
  spec names (``Frames``, ``DiaFrameMsMsWindows``, ``DiaFrameMsMsInfo``, ``PasefFrameMsMsInfo``,
  ``GlobalMetadata``); ``MsMsType`` 0 = MS1, 8 = ddaPASEF MS/MS, 9 = diaPASEF MS/MS, as the
  AlphaTims documentation of the TDF format describes it;
* ``zip_run``/``tar_run`` and the ``unsafe_*`` builders write tiny (possibly hostile) archives;
* ``pxd_sandbox`` extends ``helpers.pilot_kit.pilot_sandbox`` with a fake ``curl`` that serves
  registered URLs from files (GET, HEAD, ``--continue-at -``, ``--fail``, ``-w``), a fake ``df``
  whose free space is ``FAKE_DF_AVAIL_BYTES``, and a ``wget`` that refuses to run.

Conventions settled by the test author (M10-01 spec, docs/specs/M10-01-pxd078459-download.md):

* Record vocabulary (matched by cvParam/userParam **name**, case-insensitive; accessions in the
  fixtures are illustrative): raw files carry ``Associated raw file URI``; the Spectronaut output
  carries ``Search engine output file URI``; ``file size`` holds the byte count; ``MD5`` (or
  ``SHA-1``) holds the checksum; a ``characteristics[genotype]`` userParam is the SDRF-style sample
  annotation. The instrument is timsTOF when an ``Instrument`` cvParam name contains ``timsTOF``.
* The download job uses ``curl`` (the spec names ``curl --continue-at -``); ``df`` gives the free
  space. The committed manifest is ``validation/pilot/inputs/pxd078459_files.tsv`` in the sandbox
  checkout; tests overwrite it there and commit it inside the sandbox.
"""

from __future__ import annotations

import hashlib
import importlib.util
import io
import json
import sqlite3
import subprocess
import sys
import tarfile
import zipfile
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any
from xml.sax.saxutils import quoteattr

from helpers.pilot_kit import PilotSandbox, pilot_sandbox
from helpers.slurm_sandbox import REPO, write_exe

PILOT = REPO / "validation" / "pilot"
INPUTS = PILOT / "inputs"
MANIFEST_REL = Path("validation/pilot/inputs/pxd078459_files.tsv")
ACQUISITION_REL = Path("validation/pilot/inputs/acquisition.json")
PX_RECORD = PILOT / "px_record.py"
ACQUISITION_MODE = PILOT / "acquisition_mode.py"
DATA = REPO / "tests" / "data" / "pilot"

MANIFEST_COLUMNS = (
    "run",
    "file_name",
    "url",
    "url_record",
    "bytes",
    "checksum_record",
    "sha256",
    "sample",
    "selected_because",
)
RECORD_URL = (
    "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD078459&outputMode=XML"
)
IPROX = "download.iprox.cn"
BASE = f"https://{IPROX}/IPX0012345000/IPX0012345001"
GB = 10**9

RAW_URI = ("PRIDE", "PRIDE:0000404", "Associated raw file URI")
SEARCH_URI = ("PRIDE", "PRIDE:0000405", "Search engine output file URI")


# ------------------------------------------------------------------ module loading


def load(path: Path) -> ModuleType:
    """Import a validation script by path (``validation/`` is not a package)."""
    assert path.is_file(), f"{path.relative_to(REPO)} does not exist yet (M10-01 deliverable)"
    spec = importlib.util.spec_from_file_location(f"m10_{path.stem}", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------------ ProteomeXchange record


@dataclass(frozen=True)
class FileSpec:
    name: str
    uri: str
    raw: bool = True
    size: int | None = None
    md5: str | None = None
    genotype: str | None = None


def md5_of(text: str) -> str:
    return hashlib.md5(text.encode()).hexdigest()


def px_xml(
    files: Sequence[FileSpec],
    *,
    instrument: str = "timsTOF HT",
    accession: str = "PXD078459",
) -> bytes:
    """A ProteomeXchange XML record (format 1.4 layout) listing ``files``."""
    out = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<ProteomeXchangeDataset xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'
        f' id="{accession}" formatVersion="1.4.0"'
        ' xsi:noNamespaceSchemaLocation="proteomeXchange-1.4.0.xsd">',
        '  <DatasetSummary announceDate="2026-01-15" hostingRepository="iProX"'
        ' title="Synthetic test record: E. coli Nissle 1917 proteome">',
        "    <Description>Synthetic fixture for the GIPPA M10-01 tests; not real data."
        "</Description>",
        '    <ReviewLevel><cvParam cvRef="PRIDE" accession="PRIDE:0000414"'
        ' name="Peer-reviewed dataset"/></ReviewLevel>',
        "  </DatasetSummary>",
        "  <DatasetIdentifierList><DatasetIdentifier>",
        f'    <cvParam cvRef="MS" accession="MS:1001919" name="ProteomeXchange accession number"'
        f' value="{accession}"/>',
        "  </DatasetIdentifier></DatasetIdentifierList>",
        "  <SpeciesList><Species>",
        '    <cvParam cvRef="MS" accession="MS:1001469" name="taxonomy: scientific name"'
        ' value="Escherichia coli Nissle 1917"/>',
        "  </Species></SpeciesList>",
        "  <InstrumentList>",
        '    <Instrument id="Instrument_1">',
        f'      <cvParam cvRef="MS" accession="MS:1000031" name={quoteattr(instrument)}/>',
        "    </Instrument>",
        "  </InstrumentList>",
        "  <FullDatasetLinkList><FullDatasetLink>",
        '    <cvParam cvRef="PRIDE" accession="PRIDE:0000411" name="Dataset FTP location"'
        f' value="ftp://{IPROX}/IPX0012345000/"/>',
        "  </FullDatasetLink></FullDatasetLinkList>",
        "  <DatasetFileList>",
    ]
    for i, f in enumerate(files, start=1):
        ref, acc, name = RAW_URI if f.raw else SEARCH_URI
        out.append(f'    <DatasetFile id="FILE_{i}" name={quoteattr(f.name)}>')
        out.append(
            f'      <cvParam cvRef="{ref}" accession="{acc}" name="{name}"'
            f" value={quoteattr(f.uri)}/>"
        )
        if f.size is not None:
            out.append(
                f'      <cvParam cvRef="NCIT" accession="NCIT:C25683" name="file size"'
                f' value="{f.size}" unitName="byte"/>'
            )
        if f.md5 is not None:
            out.append(
                f'      <cvParam cvRef="MS" accession="MS:1000568" name="MD5" value="{f.md5}"/>'
            )
        if f.genotype is not None:
            out.append(
                f'      <userParam name="characteristics[genotype]" value={quoteattr(f.genotype)}/>'
            )
        out.append("    </DatasetFile>")
    out += [
        "  </DatasetFileList>",
        "  <RepositoryRecordList><RepositoryRecord repositoryID="
        f'"iProX" uri="https://www.iprox.cn/page/project.html?id=IPX0012345000"'
        f' label="{accession}" name="synthetic"/></RepositoryRecordList>',
        "</ProteomeXchangeDataset>",
        "",
    ]
    return "\n".join(out).encode("utf-8")


def example1_files() -> list[FileSpec]:
    """Spec example 1: three WT and three dLon runs, https URIs, sizes and MD5s, one .sne."""
    files = []
    for i, stem in enumerate(["WT_1", "WT_2", "WT_3", "dLon_1", "dLon_2", "dLon_3"], start=1):
        name = f"{stem}.d.zip"
        files.append(FileSpec(name, f"{BASE}/{name}", size=1_000_000_000 + i, md5=md5_of(name)))
    files.append(FileSpec("Report.sne", f"{BASE}/Report.sne", raw=False, size=5_000))
    return files


def example2_files() -> list[FileSpec]:
    """Spec example 2: SDRF genotype ``wild type`` for S01-S03, ``ΔfliC`` for S04-S06."""
    return [
        FileSpec(
            f"S0{i}.d.zip",
            f"{BASE}/S0{i}.d.zip",
            size=2_000_000_000 + i,
            md5=md5_of(f"S0{i}"),
            genotype="wild type" if i <= 3 else "ΔfliC",
        )
        for i in range(1, 7)
    ]


def example3_files() -> list[FileSpec]:
    """Spec example 3: four wild-type files in record order WT_1, WT_10, WT_2, WT_3."""
    return [
        FileSpec(
            f"WT_{n}.d.zip", f"{BASE}/WT_{n}.d.zip", size=3_000_000_000 + n, md5=md5_of(str(n))
        )
        for n in (1, 10, 2, 3)
    ]


def fixture_xml(n: int) -> Path:
    return DATA / f"px_record_example{n}.xml"


# ------------------------------------------------------------------ analysis.tdf


def make_tdf(
    path: Path,
    *,
    msms: Mapping[int, int] | None = None,
    dia_groups: int = 12,
    windows_per_group: int = 2,
    pasef_rows: int = 0,
    metadata: Mapping[str, str] | None = None,
    frame_dt: float = 0.1,
    tables: Iterable[str] = ("Frames", "DiaFrameMsMsWindows", "PasefFrameMsMsInfo"),
) -> dict[str, Any]:
    """Write a minimal ``analysis.tdf``; returns the numbers a reader should report.

    ``msms`` maps MsMsType -> frame count (default {0: 1200, 9: 10800}: spec example 9). Window
    ``k`` (0-based over all groups) isolates ``400 + 800 * k / n`` with width ``800 / n``, so the
    m/z range is exactly [400, 1200].
    """
    msms = dict({0: 1200, 9: 10800} if msms is None else msms)
    meta = dict(
        {
            "InstrumentName": "timsTOF HT",
            "InstrumentFamily": "9",
            "AcquisitionSoftwareVersion": "timsControl 5.1.8",
            "AcquisitionDateTime": "2026-01-01T10:00:00.000+01:00",
        }
        if metadata is None
        else metadata
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(path)
    tables = set(tables)
    con.execute("CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT)")
    con.executemany("INSERT INTO GlobalMetadata VALUES (?, ?)", list(meta.items()))
    con.execute(
        "CREATE TABLE Frames (Id INTEGER PRIMARY KEY, Time REAL NOT NULL, Polarity CHAR(1),"
        " ScanMode INTEGER, MsMsType INTEGER, TimsId INTEGER, MaxIntensity INTEGER,"
        " SummedIntensities INTEGER, NumScans INTEGER, NumPeaks INTEGER, MzCalibration INTEGER,"
        " T1 REAL, T2 REAL, TimsCalibration INTEGER, PropertyGroup INTEGER,"
        " AccumulationTime REAL, RampTime REAL)"
    )
    # Interleave the frame types the way an acquisition cycles them.
    order: list[int] = []
    remaining = dict(msms)
    while any(remaining.values()):
        for t in sorted(remaining):
            if remaining[t]:
                order.append(t)
                remaining[t] -= 1
    rows = []
    for i, t in enumerate(order, start=1):
        scan_mode = {0: 9 if 9 in msms else 8, 8: 8, 9: 9}.get(t, 0)
        rows.append((i, round(i * frame_dt, 6), "+", scan_mode, t, i, 1000, 10000, 927, 50, 1,
                     25.0, 25.0, 1, 1, 100.0, 100.0))  # fmt: skip
    con.executemany(f"INSERT INTO Frames VALUES ({','.join('?' * 17)})", rows)
    n_windows = dia_groups * windows_per_group
    if "DiaFrameMsMsWindows" in tables:
        con.execute(
            "CREATE TABLE DiaFrameMsMsWindows (WindowGroup INTEGER NOT NULL, ScanNumBegin INTEGER,"
            " ScanNumEnd INTEGER, IsolationMz REAL, IsolationWidth REAL, CollisionEnergy REAL)"
        )
        con.execute(
            "CREATE TABLE DiaFrameMsMsInfo (Frame INTEGER PRIMARY KEY, WindowGroup INTEGER)"
        )
        width = 800.0 / n_windows if n_windows else 0.0
        win = []
        for k in range(n_windows):
            g = k // windows_per_group + 1
            win.append((g, 100 + k * 30, 125 + k * 30, 400.0 + width * (k + 0.5), width, 30.0))
        con.executemany("INSERT INTO DiaFrameMsMsWindows VALUES (?, ?, ?, ?, ?, ?)", win)
        dia_frames = [r[0] for r in rows if r[4] == 9]
        if dia_groups:
            con.executemany(
                "INSERT INTO DiaFrameMsMsInfo VALUES (?, ?)",
                [(f, j % dia_groups + 1) for j, f in enumerate(dia_frames)],
            )
    if "PasefFrameMsMsInfo" in tables:
        con.execute(
            "CREATE TABLE PasefFrameMsMsInfo (Frame INTEGER, ScanNumBegin INTEGER,"
            " ScanNumEnd INTEGER, IsolationMz REAL, IsolationWidth REAL, CollisionEnergy REAL,"
            " Precursor INTEGER)"
        )
        dda_frames = [r[0] for r in rows if r[4] == 8] or [1]
        con.executemany(
            "INSERT INTO PasefFrameMsMsInfo VALUES (?, ?, ?, ?, ?, ?, ?)",
            [
                (dda_frames[j % len(dda_frames)], 100, 140, 500.0 + j, 2.0, 30.0, j + 1)
                for j in range(pasef_rows)
            ],
        )
    con.commit()
    con.close()
    return {
        "msms_type_counts": {str(k): v for k, v in sorted(msms.items()) if v},
        "window_groups": dia_groups if "DiaFrameMsMsWindows" in tables and n_windows else 0,
        "mz_range": [400.0, 1200.0],
        "gradient_length_s": round(len(rows) * frame_dt, 6),
        "metadata": meta,
    }


def make_run(run_dir: Path, *, tdf_bin: bool = True, **tdf: Any) -> dict[str, Any]:
    """A ``.d`` directory: ``analysis.tdf`` plus (by default) a stand-in ``analysis.tdf_bin``."""
    run_dir.mkdir(parents=True, exist_ok=True)
    info = make_tdf(run_dir / "analysis.tdf", **tdf)
    if tdf_bin:
        (run_dir / "analysis.tdf_bin").write_bytes(b"\x00stub tdf_bin\x00" * 8)
    return info


# ------------------------------------------------------------------ archives


def zip_run(archive: Path, folder: str, run_dir: Path) -> bytes:
    """Zip ``run_dir`` under ``folder/`` (e.g. ``WT_1.d`` or ``WT_1/20260101_EcN_WT_1.d``)."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(run_dir.rglob("*")):
            if p.is_file():
                zf.write(p, f"{folder}/{p.relative_to(run_dir).as_posix()}")
    archive.parent.mkdir(parents=True, exist_ok=True)
    archive.write_bytes(buf.getvalue())
    return buf.getvalue()


def tar_run(archive: Path, folder: str, run_dir: Path, *, gz: bool = False) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz" if gz else "w") as tf:
        tf.add(run_dir, arcname=folder)
    archive.parent.mkdir(parents=True, exist_ok=True)
    archive.write_bytes(buf.getvalue())
    return buf.getvalue()


def unsafe_zip(archive: Path, member: str, *, symlink_to: str | None = None) -> bytes:
    """A zip with one hostile member (``..`` path, absolute path, or a symlink)."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        if symlink_to is not None:
            zf.writestr("WT_1.d/analysis.tdf", b"x")
            info = zipfile.ZipInfo(member)
            info.create_system = 3  # unix: unzip honours the mode bits below
            info.external_attr = (0o120777 << 16) | 0  # S_IFLNK
            zf.writestr(info, symlink_to)
            zf.writestr(f"{member}/escaped.txt", b"written through the link\n")
        else:
            zf.writestr("WT_1.d/analysis.tdf", b"x")
            zf.writestr(zipfile.ZipInfo(member), b"#!/bin/sh\necho pwned\n")
    archive.parent.mkdir(parents=True, exist_ok=True)
    archive.write_bytes(buf.getvalue())
    return buf.getvalue()


def unsafe_tar(archive: Path, *, symlink_to: str) -> bytes:
    """A tar whose ``WT_1.d/link`` symlinks outside, then writes ``WT_1.d/link/escaped.txt``."""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as tf:
        d = tarfile.TarInfo("WT_1.d")
        d.type = tarfile.DIRTYPE
        d.mode = 0o755
        tf.addfile(d)
        ln = tarfile.TarInfo("WT_1.d/link")
        ln.type = tarfile.SYMTYPE
        ln.linkname = symlink_to
        tf.addfile(ln)
        payload = b"written through the link\n"
        f = tarfile.TarInfo("WT_1.d/link/escaped.txt")
        f.size = len(payload)
        tf.addfile(f, io.BytesIO(payload))
    archive.parent.mkdir(parents=True, exist_ok=True)
    archive.write_bytes(buf.getvalue())
    return buf.getvalue()


def zip_bomb(archive: Path, *, unpacked_bytes: int) -> bytes:
    """A small zip that inflates to ``unpacked_bytes`` of zeros (plus a valid run layout)."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        zf.writestr("WT_1.d/analysis.tdf", b"x")
        with zf.open("WT_1.d/analysis.tdf_bin", "w", force_zip64=True) as fh:
            chunk = b"\x00" * (1 << 20)
            left = unpacked_bytes
            while left > 0:
                fh.write(chunk[: min(left, len(chunk))])
                left -= len(chunk)
    archive.parent.mkdir(parents=True, exist_ok=True)
    archive.write_bytes(buf.getvalue())
    return buf.getvalue()


# ------------------------------------------------------------------ manifests


def manifest_text(rows: Sequence[Mapping[str, object]]) -> str:
    lines = ["\t".join(MANIFEST_COLUMNS)]
    for r in rows:
        lines.append("\t".join(str(r.get(c, "NA")) for c in MANIFEST_COLUMNS))
    return "\n".join(lines) + "\n"


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    header = lines[0].split("\t")
    return header, [dict(zip(header, ln.split("\t"), strict=True)) for ln in lines[1:]]


def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


# ------------------------------------------------------------------ fakes

_CURL = r"""#!{python}
'''Fake curl for the M10-01 tests: serves $FAKE_SLURM_DIR/server.json, never the network.

server.json: {url: {"responses": [{"status": 200, "body": "/path", "content_type": "..."}],
                     "head": {"status": 200, "content_length": N}}}
The n-th GET of a URL gets responses[min(n, len-1)]; unknown URLs answer 404.
'''
import json, os, sys

state = os.environ["FAKE_SLURM_DIR"]
argv = sys.argv[1:]
VALUE_LONG = {"--output", "--continue-at", "--write-out", "--retry", "--retry-delay",
              "--retry-max-time", "--max-time", "--connect-timeout", "--header", "--request",
              "--dump-header", "--output-dir", "--range", "--user-agent", "--referer", "--user",
              "--proto", "--proto-redir", "--cacert", "--capath", "--limit-rate", "--speed-limit",
              "--speed-time", "--max-filesize", "--stderr", "--trace", "--trace-ascii",
              "--config", "--max-redirs", "--time-cond", "--etag-save", "--etag-compare",
              "--expect100-timeout", "--keepalive-time", "--dns-servers"}
VALUE_SHORT = set("oCwmHXDrAeuTdyYzK")
opts = {"output": None, "continue_at": None, "write_out": None, "head": False, "fail": False,
        "insecure": False, "dump_header": None, "output_dir": None, "remote_name": False,
        "request": None}
urls = []

def setopt(name, val):
    if name in ("-o", "--output"): opts["output"] = val
    elif name in ("-C", "--continue-at"): opts["continue_at"] = val
    elif name in ("-w", "--write-out"): opts["write_out"] = val
    elif name in ("-D", "--dump-header"): opts["dump_header"] = val
    elif name == "--output-dir": opts["output_dir"] = val
    elif name in ("-X", "--request"): opts["request"] = val

def flag(name):
    if name in ("-I", "--head"): opts["head"] = True
    elif name in ("-f", "--fail", "--fail-with-body"): opts["fail"] = True
    elif name in ("-k", "--insecure"): opts["insecure"] = True
    elif name in ("-O", "--remote-name"): opts["remote_name"] = True

i = 0
while i < len(argv):
    a = argv[i]
    if a.startswith("--"):
        if "=" in a:
            k, v = a.split("=", 1); setopt(k, v)
        elif a in VALUE_LONG:
            setopt(a, argv[i + 1] if i + 1 < len(argv) else ""); i += 1
        else:
            flag(a)
    elif a.startswith("-") and len(a) > 1:
        j = 1
        while j < len(a):
            c = a[j]
            if c in VALUE_SHORT:
                val = a[j + 1:] if j + 1 < len(a) else (argv[i + 1] if i + 1 < len(argv) else "")
                if j + 1 >= len(a): i += 1
                setopt("-" + c, val); break
            flag("-" + c); j += 1
    else:
        urls.append(a)
    i += 1
if (opts["request"] or "").upper() == "HEAD":
    opts["head"] = True

spath = os.path.join(state, "server.json")
server = json.load(open(spath)) if os.path.exists(spath) else {}
cpath = os.path.join(state, "served.json")
served = json.load(open(cpath)) if os.path.exists(cpath) else {}

rc = 0
log = []
for url in urls:
    entry = server.get(url)
    rec = {"tool": "curl", "argv": argv, "url": url, "head": opts["head"],
           "continue_at": opts["continue_at"], "insecure": opts["insecure"],
           "output": opts["output"], "resumed_from": None, "created": []}
    if entry is None:
        status, body, ctype = 404, b"<html>404 Not Found</html>\n", "text/html"
        length = len(body)
    elif opts["head"]:
        first = entry["responses"][0]
        h = entry.get("head", {})
        status = h.get("status", first.get("status", 200))
        body = b""
        ctype = h.get("content_type", first.get("content_type", "application/zip"))
        length = h.get("content_length")
        if length is None:
            length = os.path.getsize(first["body"]) if first.get("body") else 0
    else:
        n = served.get(url, 0)
        served[url] = n + 1
        resp = entry["responses"][min(n, len(entry["responses"]) - 1)]
        status = resp.get("status", 200)
        body = open(resp["body"], "rb").read() if resp.get("body") else b""
        ctype = resp.get("content_type", "application/octet-stream")
        length = len(body)
    rec["status"] = status
    headers = (f"HTTP/1.1 {status} X\r\nContent-Length: {length}\r\n"
               f"Content-Type: {ctype}\r\n\r\n").encode()
    if opts["dump_header"]:
        if opts["dump_header"] == "-":
            sys.stdout.buffer.write(headers)
        else:
            open(opts["dump_header"], "wb").write(headers)
    target = opts["output"]
    if target is None and (opts["remote_name"] or opts["output_dir"]):
        target = os.path.join(opts["output_dir"] or os.getcwd(), url.rstrip("/").split("/")[-1])
    failed = status >= 400 and opts["fail"]
    if failed:
        rc = 22
        print(f"curl: (22) The requested URL returned error: {status}", file=sys.stderr)
    elif opts["head"]:
        if target:
            open(target, "wb").write(headers)
        else:
            sys.stdout.buffer.write(headers)
    else:
        if target:
            d = os.path.dirname(target)
            if d:
                os.makedirs(d, exist_ok=True)
            offset = 0
            if opts["continue_at"] == "-" and os.path.exists(target):
                offset = os.path.getsize(target)
                rec["resumed_from"] = offset
            with open(target, "ab" if offset else "wb") as fh:
                fh.write(body[offset:])
            rec["created"].append(os.path.abspath(target))
        else:
            sys.stdout.buffer.write(body)
    if opts["write_out"]:
        w = opts["write_out"]
        for k, v in (("%{http_code}", str(status)), ("%{response_code}", str(status)),
                     ("%{content_type}", ctype), ("%{size_download}", str(len(body))),
                     ("%{url_effective}", url), ("%header{content-length}", str(length)),
                     ("%header{Content-Length}", str(length)),
                     ("%header{content-type}", ctype), ("\\n", "\n")):
            w = w.replace(k, v)
        sys.stdout.write(w)
    log.append(rec)

json.dump(served, open(cpath, "w"))
with open(os.path.join(state, "curl.jsonl"), "a") as fh:
    for rec in log:
        fh.write(json.dumps(rec) + "\n")
sys.exit(rc)
"""

_WGET = r"""#!/bin/bash
printf '%s\n' "{\"tool\": \"wget\", \"argv\": \"$*\"}" >> "$FAKE_SLURM_DIR/curl.jsonl"
echo "fake wget: the M10-01 download uses curl (--continue-at -)" >&2
exit 1
"""

_DF = r"""#!{python}
'''Fake df: every path has FAKE_DF_AVAIL_BYTES free (default 1 PB) on a 2 PB filesystem.'''
import json, os, sys

state = os.environ["FAKE_SLURM_DIR"]
argv = sys.argv[1:]
with open(os.path.join(state, "df.jsonl"), "a") as fh:
    fh.write(json.dumps(argv) + "\n")
avail = int(os.environ.get("FAKE_DF_AVAIL_BYTES", str(10**15)))
size = max(2 * 10**15, avail)
block, output, paths, posix = 1024, None, [], False
i = 0
units = {"K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4, "KB": 1000, "MB": 1000**2,
         "GB": 1000**3, "1": 1}
def parse_block(v):
    v = v.strip()
    if v.isdigit():
        return int(v)
    num = "".join(ch for ch in v if ch.isdigit())
    suf = v[len(num):].upper()
    return int(num or 1) * units.get(suf, 1024)
while i < len(argv):
    a = argv[i]
    if a in ("-k",):
        block = 1024
    elif a in ("-m",):
        block = 1024**2
    elif a == "-P" or a == "--portability":
        posix = True
    elif a in ("-B", "--block-size"):
        block = parse_block(argv[i + 1]); i += 1
    elif a.startswith("--block-size="):
        block = parse_block(a.split("=", 1)[1])
    elif a.startswith("-B") and len(a) > 2:
        block = parse_block(a[2:])
    elif a.startswith("--output="):
        output = a.split("=", 1)[1].split(",")
    elif a == "--output":
        output = ["source", "size", "used", "avail", "pcent", "target"]
    elif a.startswith("-") and len(a) > 1 and not a.startswith("--"):
        for c in a[1:]:
            if c == "k": block = 1024
            elif c == "P": posix = True
            elif c == "m": block = 1024**2
    elif not a.startswith("-"):
        paths.append(a)
    i += 1
paths = paths or ["/"]
def fmt(n):
    return str(-(-n // block))
if output:
    names = {"source": "Filesystem", "size": "1B-blocks", "used": "Used", "avail": "Avail",
             "pcent": "Use%", "target": "Mounted on", "fstype": "Type"}
    print(" ".join(names.get(o, o) for o in output))
    for p in paths:
        vals = {"source": "fakefs", "size": fmt(size), "used": fmt(size - avail),
                "avail": str(avail // block), "pcent": "50%", "target": "/", "fstype": "nfs"}
        print(" ".join(vals.get(o, "-") for o in output))
else:
    print(f"Filesystem {block}-blocks Used Available Capacity Mounted on")
    for p in paths:
        print(f"fakefs {fmt(size)} {fmt(size - avail)} {avail // block} 50% /")
"""


@dataclass
class PxdSandbox:
    pilot: PilotSandbox

    @property
    def sb(self):  # type: ignore[no-untyped-def]
        return self.pilot.sb

    @property
    def dest(self) -> Path:
        return self.sb.root / "data" / "pxd078459"

    @property
    def payloads(self) -> Path:
        p = self.sb.tmp / "payloads"
        p.mkdir(exist_ok=True)
        return p

    def serve(
        self,
        url: str,
        *bodies: bytes | None,
        status: int | Sequence[int] = 200,
        head: Mapping[str, object] | None = None,
        content_type: str = "application/zip",
    ) -> None:
        """Register ``url``; the n-th GET gets ``bodies[n]`` (the last one repeats)."""
        spath = self.sb.state / "server.json"
        server = json.loads(spath.read_text()) if spath.is_file() else {}
        statuses = [status] * len(bodies) if isinstance(status, int) else list(status)
        responses = []
        for k, (body, st) in enumerate(zip(bodies, statuses, strict=True)):
            entry: dict[str, object] = {"status": st, "content_type": content_type}
            if body is not None:
                f = self.payloads / f"{hashlib.sha1(url.encode()).hexdigest()}.{k}"
                f.write_bytes(body)
                entry["body"] = str(f)
            responses.append(entry)
        server[url] = {"responses": responses}
        if head is not None:
            server[url]["head"] = dict(head)
        spath.write_text(json.dumps(server))

    def write_manifest(self, rows: Sequence[Mapping[str, object]]) -> Path:
        path = self.sb.repo / MANIFEST_REL
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(manifest_text(rows))
        self._commit("test manifest")
        return path

    def remove_manifest(self) -> None:
        path = self.sb.repo / MANIFEST_REL
        if path.exists():
            path.unlink()
            self._commit("no manifest")

    def _commit(self, msg: str) -> None:
        git = ["git", "-C", str(self.sb.repo)]
        subprocess.run([*git, "add", "-A"], check=True, capture_output=True, env=self.sb.env)
        subprocess.run([*git, "commit", "-q", "--allow-empty", "-m", msg], check=True,
                       capture_output=True, env=self.sb.env)  # fmt: skip

    def run(self, *args: str, **env: str) -> subprocess.CompletedProcess[str]:
        return self.pilot.run("fetch_pxd078459", *args, **env)

    def calls(self) -> list[dict[str, Any]]:
        return self.sb.jsonl("curl.jsonl")

    def gets(self) -> list[dict[str, Any]]:
        return [c for c in self.calls() if c.get("tool") == "curl" and not c["head"]]

    def heads(self) -> list[dict[str, Any]]:
        return [c for c in self.calls() if c.get("tool") == "curl" and c["head"]]

    def df_calls(self) -> list[list[str]]:
        return self.sb.jsonl("df.jsonl")

    def repo_is_clean(self) -> str:
        return self.pilot.repo_is_clean()


def pxd_sandbox(tmp: Path) -> PxdSandbox:
    pilot = pilot_sandbox(tmp, raws=0)
    py = sys.executable
    write_exe(pilot.sb.bin / "curl", _CURL.replace("{python}", py))
    write_exe(pilot.sb.bin / "wget", _WGET)
    write_exe(pilot.sb.bin / "df", _DF.replace("{python}", py))
    for var in ("PXD078459_FILES", "PXD078459_BASE_URL"):
        pilot.sb.env.pop(var, None)
    return PxdSandbox(pilot)


def out(res: subprocess.CompletedProcess[str]) -> str:
    return f"exit {res.returncode}\n--- stdout\n{res.stdout}\n--- stderr\n{res.stderr}"
