#!/usr/bin/env python3
"""Outputs byte-identical apart from names: the ADR-0064 item 8 comparator.

    python scripts/dev/rename_compare.py B_OUT A_OUT [--exclude GLOB ...] [--no-default-excludes]

``B_OUT`` is a run (or any output tree) made at the pre-rename commit B, ``A_OUT`` the same run at
the rename commit A. Exit 0 iff there are zero differences; every difference is printed.

- Files are paired by relative path with ``N`` applied to the path; both sets must be equal.
- ``N`` maps ``GIPPA``/``GIAP`` -> ``@N@``, ``Gippa``/``Giap`` -> ``@Nn@``, ``gippa``/``giap`` ->
  ``@n@`` (an uppercase hit inside a longer uppercase run is left alone: it is sequence, except
  the identifier prefixes ``GIAPORF``/``GIAPNOVEL``/``GIAPSIM``).
- ``.gz`` is decompressed and compared by its inner name; Parquet is read to a canonical TSV
  (a dtype line, the header, every row in file order).
- Sequence guard: FASTA sequence lines, TSV/CSV/Parquet sequence or peptide columns and JSON
  sequence/peptide/translation values are compared exactly, without ``N`` (GIAP is a peptide).
- PNG, PDF and every other non-UTF-8 file must be byte-equal.
- JSON: only ``git_commit``, ``git_dirty`` and ``provenance_digest`` (at any depth) and, in a
  ``run_manifest.json``, the ``VOLATILE`` pointers that the R6 byte-identity check ignores
  (``gippa.manifest.stable_view``; precomputed paths reduced to their names) are ignored.
- Every sha256 a JSON manifest records next to a path (``{"path": p, "sha256": h}``,
  ``{p: h}``, ``"hash": "sha256:h"``) whose file exists in the tree is recomputed on that side
  and must match; the value is then compared as ``@sha256-verified@``. A digest whose file is
  not in the tree (an input outside the run) is compared as written.
- Default excludes (``--no-default-excludes`` drops them): ``work/``, ``.nextflow*`` and
  ``logs/`` (Nextflow internals and wall-clock logs; the TE7 check excludes them too).
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

UPPER = re.compile(r"(?<![A-Z])(?:GIPPA|GIAP)(?=ORF|NOVEL|SIM|[^A-Z]|$)")
TITLE = re.compile(r"Gippa|Giap")
LOWER = re.compile(r"gippa|giap")
SEQ_KEY = re.compile(r"(?i)(^|_)(seq|sequence|peptide|translation|aa)($|_)")
NOT_SEQ_KEY = re.compile(r"(?i)_(id|ids|n|count|len|length|class|type|status)$")
IGNORED_KEYS = frozenset({"git_commit", "git_dirty", "provenance_digest"})
#: ``gippa.manifest.VOLATILE`` at main 45bc3599 (src/giap/manifest.py)
VOLATILE: tuple[str, ...] = (
    "/command_line", "/executor", "/profile", "/resume", "/outdir", "/started_at",
    "/finished_at", "/failed_log", "/steps", "/inputs/*/path",
)  # fmt: skip
NAME_ONLY = "/config/annotation/precomputed/*"
DEFAULT_EXCLUDES = ("work/*", ".nextflow*", "logs/*")
BYTE_EXT = frozenset({".png", ".pdf"})
HEX64 = re.compile(r"[0-9a-f]{64}")
PATH_KEYS = ("path", "file", "filename", "relpath", "name")


def N(s: str) -> str:
    """The neutral name token map, applied to both sides."""
    return LOWER.sub("@n@", TITLE.sub("@Nn@", UPPER.sub("@N@", s)))


def is_seq_key(key: str) -> bool:
    return bool(SEQ_KEY.search(key)) and not NOT_SEQ_KEY.search(key)


# ---------------------------------------------------------------- canonical forms


def _table_lines(rows: list[list[str]], sep: str) -> list[str]:
    if not rows:
        return []
    header = rows[0]
    seq = {i for i, h in enumerate(header) if is_seq_key(h)}
    out = [sep.join(N(h) for h in header)]
    for r in rows[1:]:
        out.append(sep.join(c if i in seq else N(c) for i, c in enumerate(r)))
    return out


def tsv_lines(text: str, sep: str = "\t") -> list[str]:
    """Comment lines (``#``) under N; the first other line is the header."""
    comments: list[str] = []
    body: list[str] = []
    for ln in text.splitlines():
        (comments if ln.startswith("#") and not body else body).append(ln)
    rows = [ln.split(sep) for ln in body] if sep == "\t" else list(csv.reader(body))
    return [N(c) for c in comments] + _table_lines(rows, sep)


def fasta_lines(text: str) -> list[str]:
    return [N(ln) if ln.startswith(">") else ln for ln in text.splitlines()]


def gff_lines(text: str) -> list[str]:
    out, in_fasta = [], False
    for ln in text.splitlines():
        in_fasta = in_fasta or ln.startswith("##FASTA")
        out.append(ln if in_fasta and not ln.startswith((">", "##")) else N(ln))
    return out


def parquet_lines(data: bytes) -> list[str]:
    import pyarrow.parquet as pq

    table = pq.read_table(io.BytesIO(data))
    df = table.to_pandas()
    dtypes = "\t".join(f"{N(str(c))}:{t}" for c, t in zip(df.columns, df.dtypes, strict=True))
    buf = io.StringIO()
    df.to_csv(buf, sep="\t", index=False, na_rep="<NA>", lineterminator="\n")
    return [f"#dtypes {dtypes}", *tsv_lines(buf.getvalue())]


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


def _strip_ignored(node: Any) -> Any:
    if isinstance(node, dict):
        return {k: _strip_ignored(v) for k, v in node.items() if k not in IGNORED_KEYS}
    if isinstance(node, list):
        return [_strip_ignored(v) for v in node]
    return node


def _json_flat(node: Any, ptr: str, key: str, out: list[str]) -> None:
    if isinstance(node, dict):
        if not node:
            out.append(f"{N(ptr)}\t{{}}")
        for k in sorted(node):
            _json_flat(node[k], f"{ptr}/{k}", k, out)
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
    doc = _digests(doc, side, (side.root / rel).parent, rel)
    if PurePosixPath(rel).name == "run_manifest.json":
        for ptr in VOLATILE:
            _drop(doc, ptr.split("/")[1:], name_only=False)
        _drop(doc, NAME_ONLY.split("/")[1:], name_only=True)
    out: list[str] = []
    _json_flat(_strip_ignored(doc), "", "", out)
    return out


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
        return tsv_lines(text, ",")
    if suffix in {".gff", ".gff3"}:
        return gff_lines(text)
    if suffix == ".json":
        try:
            return json_lines(text, side, name)
        except json.JSONDecodeError:
            pass
    return [N(ln) for ln in text.splitlines()]


# ---------------------------------------------------------------- the trees


def files(root: Path, excludes: Sequence[str]) -> dict[str, str]:
    """``N(relative path)`` -> relative path, for every file not excluded."""
    out: dict[str, str] = {}
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(root).as_posix()
        if any(fnmatch.fnmatch(rel, g) or fnmatch.fnmatch(rel.split("/")[0], g) for g in excludes):
            continue
        key = N(rel)
        if key in out:
            raise SystemExit(f"{root}: {out[key]} and {rel} pair to the same name {key}")
        out[key] = rel
    return out


def _first_diff(a: list[str], b: list[str]) -> Iterator[str]:
    shown = 0
    for i in range(max(len(a), len(b))):
        x = a[i] if i < len(a) else "<no line>"
        y = b[i] if i < len(b) else "<no line>"
        if x != y:
            yield f"    line {i + 1}:\n      B: {x[:300]}\n      A: {y[:300]}"
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
    n = len(files(a.b_out, excludes))
    print(f"rename_compare: {n} files in B, {len(diffs)} differences (excluded: {excludes})")
    return 1 if diffs else 0


if __name__ == "__main__":
    sys.exit(main())
