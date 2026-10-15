#!/usr/bin/env python3
"""Import a FarmShare tool recording into tests/data/recorded/ (M2-10; GIAP_SPEC.md §11.4, R11).

A recording directory ($GIAP_ROOT/recordings/<tool>/<sample>/, written by the record_*.sbatch
jobs) holds the raw outputs plus ``command.txt``, ``versions.txt`` (image digest or env-lock
sha256), ``runtime.txt`` (wall time, node, job) and ``.giap-complete``. This script copies the
kept outputs gzipped (no name, mtime 0) and writes ``PROVENANCE.md``: command, image digest or env
lock, DB, input sha256s, date (UTC), host, and a ``sha256sum`` line for every committed file.

  import_recording.py bakta SRC DEST [--strip-json-key sequences.nt] [--tsv]
  import_recording.py terminatornet SRC DEST
  import_recording.py provenance SRC DEST --input FILE [--note TEXT]   (files already in DEST)
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shutil
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

MAX_FILE = 5 * 1024 * 1024  # R11


def gz(data: bytes) -> bytes:
    """``gzip -n -9`` equivalent: identical input gives identical bytes."""
    return gzip.compress(data, compresslevel=9, mtime=0)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strip_json_key(js: dict[str, Any], dotted: str) -> None:
    """Remove ``a.b`` from ``js`` (``a`` may be a list of objects: removed from every item)."""
    head, _, key = dotted.partition(".")
    items = js.get(head)
    for it in items if isinstance(items, list) else [items]:
        if isinstance(it, dict):
            it.pop(key, None)


def write_gz(dest: Path, name: str, data: bytes) -> None:
    out = gz(data)
    if len(out) > MAX_FILE:
        raise SystemExit(
            f"{name}.gz is {len(out) / 2**20:.1f} MB > 5 MB (R11); consider --strip-json-key"
        )
    (dest / f"{name}.gz").write_bytes(out)


def environment(versions: str) -> str:
    """The image digest or env-lock line of ``versions.txt``."""
    m = re.search(r"sha256:[0-9a-f]{64}", versions)
    if m:
        return f"- Image: `{m.group(0)}` (versions.txt: {versions.splitlines()[0].strip()})"
    m = re.search(r"^([0-9a-f]{64})\s+(\S+\.explicit\.txt)$", versions, re.MULTILINE)
    if not m:
        raise SystemExit("versions.txt has neither an image digest nor an env lock sha256")
    commit = re.search(r"commit\s+([0-9a-f]{40})", versions)
    extra = f"; tool commit `{commit.group(1)}`" if commit else ""
    return f"- Env lock: `workflow/envs/locks/{m.group(2)}` sha256 `{m.group(1)}`{extra}"


def provenance(src: Path, dest: Path, inputs: list[Path], db: str, notes: list[str]) -> None:
    versions = (src / "versions.txt").read_text()
    runtime = dict(kv.split("=", 1) for kv in (src / "runtime.txt").read_text().split())
    done = datetime.fromtimestamp((src / ".giap-complete").stat().st_mtime, UTC)
    files = sorted(p for p in dest.rglob("*") if p.is_file() and p.name != "PROVENANCE.md")
    sums = "\n".join(f"{sha(p)}  {p.relative_to(dest).as_posix()}" for p in files)
    lines = (src / "command.txt").read_text().splitlines()
    cmd = " ".join(ln.strip() for ln in lines if ln.strip() and not ln.startswith("#"))
    remarks = "".join(f" ({ln.lstrip('# ').strip()})" for ln in lines if ln.startswith("#"))
    cwd = Path.cwd()
    shown = [p.relative_to(cwd) if p.is_relative_to(cwd) else p for p in inputs]
    ins = "\n".join(f"  - `{s}`: sha256 `{sha(p)}`" for s, p in zip(shown, inputs, strict=True))
    text = (
        f"# {dest.parent.parent.name} {dest.parent.name} on {dest.name}"
        " (M2-10, GIAP_SPEC.md §11.4)\n\n"
        f"- Command: `{cmd}`{remarks}\n"
        f"{environment(versions)}\n"
        + (f"- Database: {db}\n" if db else "")
        + f"- Date: {done:%Y-%m-%dT%H:%M:%SZ} (run finished, UTC)\n"
        f"- Host: {runtime.get('node', 'unknown')} (FarmShare; SLURM_JOB_ID "
        f"{runtime.get('job', 'not recorded')}; {runtime.get('cpus', '16')} CPUs; wall time "
        f"{runtime.get('wall_s', '?')} s)\n"
        f"- Recording: `{src}`\n"
        f"- Inputs:\n{ins}\n"
        + "".join(f"- {n}\n" for n in notes)
        + f"\n## Files (sha256sum)\n\n```\n{sums}\n```\n"
    )
    (dest / "PROVENANCE.md").write_text(text)


def bakta(src: Path, dest: Path, strip: list[str], tsv: bool) -> None:
    (out,) = (src / "bakta").glob("*.json")
    prefix = out.stem
    js = json.loads(out.read_bytes())
    notes = [f"Bakta outputs from `bakta/{prefix}.*`; raw JSON sha256 `{sha(out)}`"]
    data = out.read_bytes()
    if strip:
        for key in strip:
            strip_json_key(js, key)
        data = json.dumps(js, indent=1).encode() + b"\n"
        notes.append(
            f"JSON keys stripped (R11, M2-10 Q1): {', '.join(strip)}; the genome is "
            "tests/data/genomes/<sample>.fna.gz"
        )
    write_gz(dest, f"{prefix}.json", data)
    for ext in ("gff3", "gbff", "faa") + (("tsv",) if tsv else ()):
        write_gz(dest, f"{prefix}.{ext}", (src / "bakta" / f"{prefix}.{ext}").read_bytes())
    pre = src / "preflight" / "work" / "preflight"
    for name in ("replicons.tsv", "contig_map.tsv"):
        shutil.copyfile(pre / name, dest / name)
    v = js["version"]["db"]
    vj = json.loads((src / "versions.txt").read_text().split("\n", 2)[2])
    db = f"Bakta DB {v['type']} v{v['version']} ({vj.get('date')}, doi {vj.get('doi')})"
    genome = Path((src / "input.sha256").read_text().split()[1])
    notes.append(
        f"Preflight (`giap preflight`) wrote replicons.tsv and contig_map.tsv from `{genome.name}`"
    )
    provenance(src, dest, [genome, pre / "genome.norm.fa"], db, notes)


def terminatornet(src: Path, dest: Path) -> None:
    write_gz(dest, "terms.tsv", (src / "terminatornet.tsv").read_bytes())
    notes = [
        f"`terms.tsv.gz` = `terminatornet.tsv` verbatim (sha256 `{sha(src / 'terminatornet.tsv')}`)"
    ]
    provenance(src, dest, [(src / "genome.fa").resolve()], "", notes)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=("bakta", "terminatornet", "provenance"))
    ap.add_argument("src", type=Path)
    ap.add_argument("dest", type=Path)
    ap.add_argument("--strip-json-key", action="append", default=[])
    ap.add_argument("--tsv", action="store_true", help="bakta: keep the TSV too")
    ap.add_argument("--input", action="append", type=Path, default=[])
    ap.add_argument("--db", default="")
    ap.add_argument("--note", action="append", default=[])
    a = ap.parse_args()
    a.dest.mkdir(parents=True, exist_ok=True)
    if a.mode == "bakta":
        bakta(a.src, a.dest, a.strip_json_key, a.tsv)
    elif a.mode == "terminatornet":
        terminatornet(a.src, a.dest)
    else:
        provenance(a.src, a.dest, a.input, a.db, a.note)


if __name__ == "__main__":
    main()
