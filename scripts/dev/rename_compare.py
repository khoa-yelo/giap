#!/usr/bin/env python3
"""Outputs byte-identical apart from names: the ADR-0064 item 8 comparator.

    python scripts/dev/rename_compare.py B_OUT A_OUT [--exclude GLOB ...] [--no-default-excludes]

``B_OUT`` is a run (or any output tree) made at the pre-rename commit B, ``A_OUT`` the same run at
the rename commit A. Exit 0 iff there are zero differences; every difference is printed.

- Files are paired by relative path with ``N`` applied to the path; both sets must be equal.
- ``N`` maps ``GIPPA``/``GIAP`` -> ``@N@``, ``Gippa``/``Giap`` -> ``@Nn@``, ``gippa``/``giap`` ->
  ``@n@`` (an uppercase hit inside a longer uppercase run is left alone: it is sequence, except
  the identifier prefixes ``GIAPORF``/``GIAPNOVEL``).
- Text is compared exactly after ``N``: line endings, the final newline and whitespace count.
- ``.gz`` is decompressed and compared by its inner name; Parquet is read to its Arrow schema
  (with key-value metadata) and a canonical TSV (dtypes, header, every row in file order).
- Sequence guard: FASTA sequence lines, TSV/CSV/Parquet columns and JSON keys whose name holds
  ``seq``/``peptide``/``translation`` (any case; not ``*_id``/``*_count``-like) and
  ``Precursor.Id`` are compared exactly, without ``N`` (GIAP is a peptide).
- CSV records are compared as parsed tuples (quoting and column count count) and as raw text.
- PNG, PDF and every other non-UTF-8 file must be byte-equal.
- JSON is compared as raw text under ``N`` plus its parsed values in document order (key order
  counts). Only a ``run_manifest.json`` drops fields: its top-level ``git_commit``,
  ``git_dirty`` and ``provenance_digest`` and the ``VOLATILE`` pointers the R6 byte-identity
  check ignores (``gippa.manifest.stable_view``; precomputed paths reduced to their names); it
  and any JSON with a recomputed digest are compared parsed only (key order still counts).
- JSON Lines (``logs/events.jsonl``, ``logs/steps/*.events.jsonl``): raw text under ``N`` with
  only each record's top-level ``ts`` value masked.
- Every sha256 a JSON manifest records next to a path (``{"path": p, "sha256": h}``,
  ``{p: h}``, ``"hash": "sha256:h"``) whose file exists in the tree is recomputed on that side
  and must match; the value is then compared as ``@sha256-verified@``. A digest whose file is
  not in the tree (an input outside the run) is compared as written.
- Default excludes (``DEFAULT_EXCLUDES``, each with its reason; ``--no-default-excludes`` drops
  them) are listed with every excluded file when the comparator runs.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import gzip
import hashlib
import io
import json
import re
import sys
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Any

UPPER = re.compile(r"(?<![A-Z])(?:GIPPA|GIAP)(?=ORF|NOVEL|[^A-Z]|$)")
TITLE = re.compile(r"Gippa|Giap")
LOWER = re.compile(r"gippa|giap")
SEQ_KEY = re.compile(r"(?i)seq|peptide|translation|(^|[._ ])aa($|[._ ])")
NOT_SEQ_KEY = re.compile(r"(?i)[._ ](id|ids|n|count|len|length|class|type|status)$")
PRECURSOR_KEY = re.compile(r"(?i)precursor[._ ]?id")
IGNORED_KEYS = frozenset({"git_commit", "git_dirty", "provenance_digest"})  # manifest top level
#: ``gippa.manifest.VOLATILE`` minus ``NAME_ONLY`` (tests/unit/test_rename_compare.py checks it)
VOLATILE: tuple[str, ...] = (
    "/command_line", "/executor", "/profile", "/resume", "/outdir", "/started_at",
    "/finished_at", "/failed_log", "/steps", "/inputs/*/path", "/launch",
)  # fmt: skip
NAME_ONLY = "/config/annotation/precomputed/*"
#: glob -> why it is not compared (a glob matches the relative path or its first component)
EXCLUDE_REASONS: dict[str, str] = {
    "work/*": "Nextflow task directories: hashed names, absolute paths, wall-clock files",
    ".nextflow*": "Nextflow cache, history and .nextflow.log: session IDs and timestamps",
    "logs/trace.tsv": "Nextflow trace: task hashes, native job IDs, times, CPU and memory use",
    "logs/timeline.html": "Nextflow timeline: wall-clock task times",
    "logs/nf_report.html": "Nextflow execution report: wall-clock times and resource use",
}
DEFAULT_EXCLUDES: tuple[str, ...] = tuple(EXCLUDE_REASONS)
BYTE_EXT = frozenset({".png", ".pdf"})
HEX64 = re.compile(r"[0-9a-f]{64}")
PATH_KEYS = ("path", "file", "filename", "relpath", "name")
TS = re.compile(r'("ts"\s*:\s*)"(?:[^"\\]|\\.)*"')
LINE = re.compile(r"[^\n]*\n|[^\n]+\Z")


def N(s: str) -> str:
    """The neutral name token map, applied to both sides."""
    return LOWER.sub("@n@", TITLE.sub("@Nn@", UPPER.sub("@N@", s)))


def is_seq_key(key: str) -> bool:
    if PRECURSOR_KEY.fullmatch(key):
        return True
    return bool(SEQ_KEY.search(key)) and not NOT_SEQ_KEY.search(key)


def lines(text: str) -> list[str]:
    """``text`` split after each ``\\n``, endings kept: ``"".join(lines(t)) == t``."""
    return LINE.findall(text)


def _ending(ln: str) -> tuple[str, str]:
    body = ln.rstrip("\r\n")
    return body, ln[len(body) :]


# ---------------------------------------------------------------- canonical forms


def tsv_lines(text: str) -> list[str]:
    """Comment lines (``#``) under N; the first other line is the header."""
    out: list[str] = []
    seq: set[int] | None = None
    for ln in lines(text):
        body, end = _ending(ln)
        if seq is None and body.startswith("#"):
            out.append(N(ln))
            continue
        cells = body.split("\t")
        if seq is None:
            seq = {i for i, h in enumerate(cells) if is_seq_key(h)}
            out.append(N(ln))
        else:
            out.append("\t".join(c if i in seq else N(c) for i, c in enumerate(cells)) + end)
    return out


def csv_lines(text: str) -> list[str]:
    """Every record as a JSON list (sequence columns exact), then the raw text under N."""
    rows = list(csv.reader(io.StringIO(text, newline="")))
    out: list[str] = []
    if rows:
        seq = {i for i, h in enumerate(rows[0]) if is_seq_key(h)}
        for r in rows:
            cells = [c if i in seq and r is not rows[0] else N(c) for i, c in enumerate(r)]
            out.append(json.dumps(cells, ensure_ascii=False))
    return [*out, "#raw", *(N(ln) for ln in lines(text))]


def fasta_lines(text: str) -> list[str]:
    return [N(ln) if ln.startswith(">") else ln for ln in lines(text)]


def gff_lines(text: str) -> list[str]:
    out, in_fasta = [], False
    for ln in lines(text):
        in_fasta = in_fasta or ln.startswith("##FASTA")
        out.append(ln if in_fasta and not ln.startswith((">", "##")) else N(ln))
    return out


def parquet_lines(data: bytes) -> list[str]:
    import pyarrow.parquet as pq

    table = pq.read_table(io.BytesIO(data))
    schema = [f"#schema {N(ln)}" for ln in str(table.schema).splitlines()]
    df = table.to_pandas()
    dtypes = "\t".join(f"{N(str(c))}:{t}" for c, t in zip(df.columns, df.dtypes, strict=True))
    buf = io.StringIO()
    df.to_csv(buf, sep="\t", index=False, na_rep="<NA>", lineterminator="\n")
    return [*schema, f"#dtypes {dtypes}", *tsv_lines(buf.getvalue())]


def jsonl_lines(text: str) -> list[str]:
    """Each record under N with only its top-level ``ts`` value masked."""
    out = []
    for ln in lines(text):
        masked = TS.sub(r'\1"@ts@"', ln, count=1)
        try:
            ok = masked != ln and json.loads(masked).get("ts") == "@ts@"
        except (json.JSONDecodeError, AttributeError):
            ok = False
        out.append(N(masked if ok else ln))
    return out


@dataclass
class Side:
    """One tree: its root and the problems its manifests' digests showed."""

    root: Path
    label: str
    digest_problems: list[str] = field(default_factory=list)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _resolve(side: Side, base: Path, ref: str) -> Path | None:
    if not ref or "\n" in ref or len(ref) > 4096:
        return None
    p = Path(ref)
    for cand in (p,) if p.is_absolute() else (base / p, side.root / p):
        try:
            if cand.is_file() and side.root.resolve() in cand.resolve().parents:
                return cand
        except OSError:
            continue
    return None


def _verify(side: Side, base: Path, ref: str, digest: str, where: str) -> str:
    target = _resolve(side, base, ref)
    if target is None:
        return digest
    got = _sha(target)
    if got != digest:
        side.digest_problems.append(f"{side.label}: {where}: recorded sha256 {digest} != {got}")
        return f"@sha256-mismatch@{digest}"
    return "@sha256-verified@"


def _digests(node: Any, side: Side, base: Path, where: str) -> Any:
    """Recompute and replace every sha256 recorded next to a path that exists in the tree."""
    if isinstance(node, list):
        return [_digests(v, side, base, f"{where}/{i}") for i, v in enumerate(node)]
    if not isinstance(node, dict):
        return node
    out = {k: _digests(v, side, base, f"{where}/{k}") for k, v in node.items()}
    ref = next((node[k] for k in PATH_KEYS if isinstance(node.get(k), str)), None)
    for key in ("sha256", "hash"):
        val = out.get(key)
        if isinstance(val, str) and ref is not None:
            hexd = val.removeprefix("sha256:")
            if HEX64.fullmatch(hexd):
                out[key] = _verify(side, base, ref, hexd, f"{where}/{key}")
    for k, v in list(out.items()):  # {"tables/genes.tsv": "<hex>"}
        hexv = isinstance(v, str) and HEX64.fullmatch(v) and k not in ("sha256", "hash")
        if hexv and _resolve(side, base, k) is not None:
            out[k] = _verify(side, base, k, v, f"{where}/{k}")
    return out


def _drop(node: Any, parts: list[str], name_only: bool) -> None:
    if not parts:
        return
    head, rest = parts[0], parts[1:]
    if isinstance(node, dict):
        keys = list(node) if head == "*" else [head] if head in node else []
        for k in keys:
            if not rest:
                if name_only:
                    v = node[k]
                    node[k] = PurePosixPath(v).name if isinstance(v, str) else v
                else:
                    del node[k]
            else:
                _drop(node[k], rest, name_only)
    elif isinstance(node, list) and head == "*":
        for v in node:
            _drop(v, rest, name_only)


def _json_flat(node: Any, ptr: str, key: str, out: list[str]) -> None:
    """Every leaf in document order (key order counts); sequence-keyed strings exact."""
    if isinstance(node, dict):
        if not node:
            out.append(f"{N(ptr)}\t{{}}")
        for k, v in node.items():
            _json_flat(v, f"{ptr}/{k}", k, out)
    elif isinstance(node, list):
        if not node:
            out.append(f"{N(ptr)}\t[]")
        for i, v in enumerate(node):
            _json_flat(v, f"{ptr}/{i}", key, out)
    else:
        val = json.dumps(node, ensure_ascii=False)
        out.append(f"{N(ptr)}\t{val if isinstance(node, str) and is_seq_key(key) else N(val)}")


def json_lines(data: str, side: Side, rel: str) -> list[str]:
    doc = json.loads(data)
    checked = _digests(doc, side, (side.root / rel).parent, rel)
    parsed_only = checked != doc
    if PurePosixPath(rel).name == "run_manifest.json" and isinstance(checked, dict):
        parsed_only = True
        for k in IGNORED_KEYS:
            checked.pop(k, None)
        for ptr in VOLATILE:
            _drop(checked, ptr.split("/")[1:], name_only=False)
        _drop(checked, NAME_ONLY.split("/")[1:], name_only=True)
    out: list[str] = []
    _json_flat(checked, "", "", out)
    return out if parsed_only else [*out, "#raw", *(N(ln) for ln in lines(data))]


def canonical(rel: str, data: bytes, side: Side) -> list[str] | bytes:
    """Lines to compare, or bytes when the file must be byte-equal."""
    name = rel
    if name.endswith(".gz"):
        try:
            data = gzip.decompress(data)
        except (OSError, EOFError):
            return data
        name = name[:-3]
    suffix = PurePosixPath(name).suffix.lower()
    if suffix == ".parquet":
        return parquet_lines(data)
    if suffix in BYTE_EXT or b"\0" in data:
        return data
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return data
    if suffix in {".fa", ".faa", ".fasta", ".fna", ".ffn", ".frn"} or text.startswith(">"):
        return fasta_lines(text)
    if suffix in {".tsv", ".tab"}:
        return tsv_lines(text)
    if suffix == ".csv":
        return csv_lines(text)
    if suffix in {".gff", ".gff3"}:
        return gff_lines(text)
    if suffix == ".jsonl":
        return jsonl_lines(text)
    if suffix == ".json":
        try:
            return json_lines(text, side, name)
        except json.JSONDecodeError:
            pass
    return [N(ln) for ln in lines(text)]


# ---------------------------------------------------------------- the trees


def _excluded(rel: str, excludes: Sequence[str]) -> bool:
    return any(fnmatch.fnmatch(rel, g) or fnmatch.fnmatch(rel.split("/")[0], g) for g in excludes)


def files(root: Path, excludes: Sequence[str]) -> dict[str, str]:
    """``N(relative path)`` -> relative path, for every file not excluded."""
    out: dict[str, str] = {}
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(root).as_posix()
        if _excluded(rel, excludes):
            continue
        key = N(rel)
        if key in out:
            raise SystemExit(f"{root}: {out[key]} and {rel} pair to the same name {key}")
        out[key] = rel
    return out


def excluded_files(root: Path, excludes: Sequence[str]) -> list[str]:
    """Relative paths of the files ``excludes`` keeps out of the comparison."""
    rels = (p.relative_to(root).as_posix() for p in sorted(root.rglob("*")) if p.is_file())
    return [r for r in rels if _excluded(r, excludes)]


def _first_diff(a: list[str], b: list[str]) -> Iterator[str]:
    shown = 0
    for i in range(max(len(a), len(b))):
        x = a[i] if i < len(a) else "<no line>"
        y = b[i] if i < len(b) else "<no line>"
        if x != y:
            yield f"    line {i + 1}:\n      B: {x[:300]!r}\n      A: {y[:300]!r}"
            shown += 1
            if shown == 5:
                yield "    …"
                return


def compare(b_root: Path, a_root: Path, excludes: Sequence[str]) -> list[str]:
    """Every difference, readable; empty when the trees are equal apart from names."""
    sb, sa = Side(b_root, "B"), Side(a_root, "A")
    fb, fa = files(b_root, excludes), files(a_root, excludes)
    diffs = [f"only in B: {fb[k]}" for k in sorted(fb.keys() - fa.keys())]
    diffs += [f"only in A: {fa[k]}" for k in sorted(fa.keys() - fb.keys())]
    for key in sorted(fb.keys() & fa.keys()):
        rb, ra = fb[key], fa[key]
        cb = canonical(rb, (b_root / rb).read_bytes(), sb)
        ca = canonical(ra, (a_root / ra).read_bytes(), sa)
        if cb == ca:
            continue
        if isinstance(cb, bytes) or isinstance(ca, bytes):
            diffs.append(f"bytes differ: {rb} | {ra}")
        else:
            diffs.append("\n".join([f"content differs: {rb} | {ra}", *_first_diff(cb, ca)]))
    return diffs + sb.digest_problems + sa.digest_problems


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("b_out", type=Path, help="output tree at the pre-rename commit B")
    ap.add_argument("a_out", type=Path, help="the same output tree at the rename commit A")
    ap.add_argument("--exclude", action="append", default=[], metavar="GLOB")
    ap.add_argument("--no-default-excludes", action="store_true")
    a = ap.parse_args(argv)
    excludes = [*([] if a.no_default_excludes else DEFAULT_EXCLUDES), *a.exclude]
    diffs = compare(a.b_out, a.a_out, excludes)
    for d in diffs:
        print(d)
    for g in excludes:
        print(f"excluded {g}: {EXCLUDE_REASONS.get(g, 'given with --exclude')}")
    for label, root in (("B", a.b_out), ("A", a.a_out)):
        for rel in excluded_files(root, excludes):
            print(f"excluded file {label}: {rel}")
    n = len(files(a.b_out, excludes))
    print(f"rename_compare: {n} files in B, {len(diffs)} differences")
    return 1 if diffs else 0


if __name__ == "__main__":
    sys.exit(main())
