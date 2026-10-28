"""Precomputed annotation layers (``annotation.precomputed``; GIPPA_SPEC.md §5.2, §6.4; M2-02).

Each layer is a TSV following ``gippa.io.tsv`` (1-based closed coordinates, ``NA``, ``;`` lists)
and equals the output of the parse step it replaces. ``load_precomputed`` returns a validated
frame with 0-based coordinates: ``start``/``end`` become half-open ``[start, end)`` and single
positions (``end3p``, ``tss``) become 0-based positions. Coordinates are checked with
``check_interval``; origin-crossing rows keep ``end > length`` on circular replicons only.
``rbs_novel`` keys rows by a single ``start`` position.
AMR rows of ``regions`` may have NA ``start``/``end`` (nullable ``Int64``). A bad value raises
``ValueError`` naming the key, column and 1-based data row (header not counted).
``bakta_json`` is a Bakta JSON file, not a TSV, and is parsed by the Bakta parser (M2-03).
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

import pandas as pd

from gippa.io.tsv import read_tsv
from gippa.model.coords import Interval, Replicon, Strand, check_interval

# Column kinds: str, float, int (nullable), list (";"-joined), pos (1-based position),
# start/end (1-based closed interval bounds; converted together).
_SCHEMA: dict[str, dict[str, str]] = {
    "operon_pairs": {"gene_a": "str", "gene_b": "str", "uniop_p": "float", "source": "str"},
    "promoters": {"contig": "str", "strand": "str", "end3p": "pos", "score": "float"},
    "tx_rates": {"contig": "str", "strand": "str", "tss": "pos", "tx_rate": "float"},
    "terminators": {
        "contig": "str",
        "start": "start",
        "end": "end",
        "strand": "str",
        "confidence": "float",
        "source": "str",
    },
    "rbs": {
        "locus_tag": "str",
        "tir_ostir": "float",
        "sd_dg": "float",
        "mrna_dg": "float",
        "spacing_dg": "float",
        "standby_dg": "float",
        "spacer_bp": "int",
        "start_codon": "str",
    },
    "regions": {
        "region_id": "str",
        "region_type": "str",
        "contig": "str",
        "start": "start",
        "end": "end",
        "source_tool": "str",
        "score": "float",
        "genes": "list",
        "notes": "str",
    },
    # second-round OSTIR (ADR-0024): ``start`` = the start codon's first base in transcription
    # direction (on ``-`` its highest coordinate), a single position
    "rbs_novel": {
        "contig": "str",
        "strand": "str",
        "start": "pos",
        "tir_ostir": "float",
        "sd_dg": "float",
        "mrna_dg": "float",
        "spacing_dg": "float",
        "standby_dg": "float",
        "spacer_bp": "int",
        "start_codon": "str",
    },
}

PRECOMPUTED_COLUMNS: Mapping[str, tuple[str, ...]] = {k: tuple(v) for k, v in _SCHEMA.items()}


def _to_float(v: str | None) -> float:
    if v is None:
        return float("nan")
    return float({"INF": "inf", "-INF": "-inf"}.get(v, v))


def _to_int(v: str | None) -> int | None:
    return None if v is None else int(v)


def _to_list(v: str | None) -> list[str]:
    return [] if v is None else v.split(";")


_CONVERT: dict[str, Callable[[str | None], Any]] = {
    "str": lambda v: v,
    "float": _to_float,
    "int": _to_int,
    "list": _to_list,
}
_DTYPE = {"float": "float64", "int": "Int64", "pos": "int64", "start": "int64", "end": "int64"}
# regions.tsv: AMR rows may have NA start/end, so its coordinates are nullable integers.
_DTYPE_OVERRIDE: dict[str, dict[str, str]] = {"regions": {"start": "Int64", "end": "Int64"}}


class _CellError(ValueError):
    """A bad value; ``load_precomputed`` adds the key, column and row."""

    def __init__(self, col: str, msg: str) -> None:
        super().__init__(msg)
        self.col = col


def _replicon(row: Mapping[str, Any], replicons: Mapping[str, Replicon]) -> Replicon:
    contig = row["contig"]
    if contig not in replicons:
        raise _CellError("contig", f"unknown contig {contig!r}")
    return replicons[contig]


def _required_int(col: str, v: str | None) -> int:
    if v is None:
        raise _CellError(col, "may not be NA")
    try:
        return int(v)
    except ValueError:
        raise _CellError(col, f"not an integer: {v!r}") from None


def _convert_row(
    key: str, schema: Mapping[str, str], row: dict[str, str | None], reps: Mapping[str, Replicon]
) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for col, kind in schema.items():
        if kind in _CONVERT:
            try:
                out[col] = _CONVERT[kind](row[col])
            except ValueError as exc:
                raise _CellError(col, f"invalid {kind} {row[col]!r} ({exc})") from None
    if "strand" in schema and out["strand"] not in ("+", "-", "."):
        raise _CellError("strand", f"invalid strand {out['strand']!r}")
    if schema.get("start") == "start":
        rep = _replicon(row, reps)
        if (
            key == "regions"
            and row["region_type"] == "amr"
            and row["start"] is None
            and row["end"] is None
        ):
            out["start"] = out["end"] = None  # AMR hits may lack coordinates (M2 decision 6)
        else:
            start1 = _required_int("start", row["start"])
            end1 = _required_int("end", row["end"])
            strand: Strand = row.get("strand") or "."  # type: ignore[assignment]
            if start1 < 1:
                raise _CellError("start", f"start must be >= 1, got {start1}")
            iv = Interval(rep.id, start1 - 1, end1, strand)
            try:
                check_interval(iv, rep)
            except ValueError as exc:
                raise _CellError("start/end", str(exc)) from None
            out["start"], out["end"] = iv.start, iv.end
    for col, kind in schema.items():
        if kind == "pos":
            rep = _replicon(row, reps)
            pos1 = _required_int(col, row[col])
            if not 1 <= pos1 <= rep.length:
                raise _CellError(col, f"{pos1} outside {rep.id} (1..{rep.length})")
            out[col] = pos1 - 1
    return out


def load_precomputed(key: str, path: Path, replicons: Mapping[str, Replicon]) -> pd.DataFrame:
    """Read the ``annotation.precomputed.<key>`` TSV into a validated 0-based frame."""
    if key not in _SCHEMA:
        raise ValueError(f"precomputed {key}: not a TSV layer (known: {', '.join(_SCHEMA)})")
    schema = _SCHEMA[key]
    header, rows = read_tsv(path)
    for col in schema:
        if col not in header:
            raise ValueError(f"precomputed {key}: missing column {col}")
    records = []
    for i, r in enumerate(rows, start=1):  # 1-based data row; the header is not counted
        try:
            records.append(_convert_row(key, schema, dict(zip(header, r, strict=True)), replicons))
        except _CellError as exc:
            raise ValueError(f"precomputed {key}: row {i}, column {exc.col}: {exc}") from None
    df = pd.DataFrame.from_records(records, columns=list(schema))
    dtypes = {c: _DTYPE[k] for c, k in schema.items() if k in _DTYPE}
    dtypes.update(_DTYPE_OVERRIDE.get(key, {}))
    return df.astype(dtypes)
