"""In-pipeline checks of a step's written outputs (GIAP_SPEC.md §6.3, §9.5; M7-15 item 2).

A step calls these after its writers; any problem raises ``ValueError`` naming the output key and
the file, so the step fails (``giap run`` exit 4 naming the step) instead of a downstream reader.

* ``check_tsv``: the TSV against its Frictionless resource (``giap/schemas/tables/<stem>``): the
  schema columns present, the integer/number/boolean types, ``required``, ``enum`` and
  ``minimum``. The dev-only ``frictionless`` package is not needed.
* ``check_vcf``: the VCF 4.x structure GIAP writes (header, ``#CHROM`` line, eight or more columns
  per record, a declared contig, a positive integer POS).
* ``check_gff3``: the M1-05 ``read_gff3`` reader on the file.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from functools import cache
from importlib import resources
from pathlib import Path
from typing import Any

from giap.io._open import open_text
from giap.io.gff3 import read_gff3
from giap.io.tsv import read_tsv
from giap.model.coords import Replicon

__all__ = ["check_gff3", "check_tsv", "check_vcf"]

_BOOLS = frozenset({"true", "false"})


@cache
def _schema(stem: str) -> dict[str, Any]:
    res = resources.files("giap.schemas").joinpath("tables", f"{stem}.schema.json")
    data: dict[str, Any] = json.loads(res.read_text(encoding="utf-8"))
    return data


def _cell_problem(value: str, field: Mapping[str, Any]) -> str | None:
    kind = field.get("type", "string")
    cons = field.get("constraints", {})
    num: float | None = None
    try:
        if kind == "integer":
            num = int(value)
        elif kind == "number":
            num = float(value)
    except ValueError:
        return f"{value!r} is not of type {kind}"
    if kind == "boolean" and value not in _BOOLS:
        return f"{value!r} is not of type boolean"
    if "enum" in cons and value not in cons["enum"]:
        return f"{value!r} is not one of {cons['enum']}"
    if num is not None and "minimum" in cons and num < cons["minimum"]:
        return f"{value!r} is below the minimum {cons['minimum']}"
    return None


def check_tsv(key: str, path: Path, stem: str) -> None:
    """``path`` (gzip allowed) conforms to the ``stem`` table schema, else ``ValueError``."""
    try:
        header, rows = read_tsv(path)
    except ValueError as exc:
        raise ValueError(f"output {key} ({path}): {exc}") from exc
    fields = _schema(stem)["fields"]
    absent = [f["name"] for f in fields if f["name"] not in header]
    if absent:
        raise ValueError(f"output {key} ({path}): schema column(s) missing: {absent}")
    cols = [(header.index(f["name"]), f) for f in fields]
    for lineno, row in enumerate(rows, start=2):
        for j, f in cols:
            v = row[j]
            if v is None:
                if f.get("constraints", {}).get("required") and "missingValues" not in f:
                    raise ValueError(f"output {key} ({path}): line {lineno}: {f['name']} is NA")
                continue
            problem = _cell_problem(v, f)
            if problem is not None:
                raise ValueError(f"output {key} ({path}): line {lineno}: {f['name']} {problem}")


def check_vcf(key: str, path: Path) -> None:
    """``path`` is a structurally valid VCF as GIAP writes it, else ``ValueError``."""
    with open_text(path) as fh:
        lines = fh.read().split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    where = f"output {key} ({path})"
    if not lines or not lines[0].startswith("##fileformat=VCF"):
        raise ValueError(f"{where}: no ##fileformat=VCF header line")
    contigs: set[str] = set()
    n_cols = 0
    for lineno, line in enumerate(lines, start=1):
        if line.startswith("##contig=<ID="):
            contigs.add(line[len("##contig=<ID=") :].split(",")[0].rstrip(">"))
        elif line.startswith("#CHROM"):
            n_cols = len(line.split("\t"))
        elif not line.startswith("#"):
            cells = line.split("\t")
            if not n_cols:
                raise ValueError(f"{where}: line {lineno}: record before the #CHROM line")
            if len(cells) != n_cols or n_cols < 8:
                raise ValueError(f"{where}: line {lineno}: {len(cells)} columns, want {n_cols}")
            if cells[0] not in contigs:
                raise ValueError(f"{where}: line {lineno}: undeclared contig {cells[0]!r}")
            if not cells[1].isdigit() or int(cells[1]) < 1:
                raise ValueError(f"{where}: line {lineno}: POS {cells[1]!r} is not >= 1")
    if not n_cols:
        raise ValueError(f"{where}: no #CHROM line")


def check_gff3(key: str, path: Path, replicons: Mapping[str, Replicon]) -> None:
    """``path`` reads back with the M1-05 GFF3 reader, else ``ValueError``."""
    try:
        read_gff3(path, replicons)
    except ValueError as exc:
        raise ValueError(f"output {key} ({path}): {exc}") from exc
