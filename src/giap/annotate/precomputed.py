"""Precomputed annotation layers (``annotation.precomputed``; GIAP_SPEC.md §5.2, §6.4; M2-02).

Each layer is a TSV following ``giap.io.tsv`` (1-based closed coordinates, ``NA``, ``;`` lists)
and equals the output of the parse step it replaces. ``load_precomputed`` returns a validated
frame with 0-based coordinates: ``start``/``end`` become half-open ``[start, end)`` and single
positions (``end3p``, ``tss``) become 0-based positions. Coordinates are checked with
``check_interval``; origin-crossing rows keep ``end > length`` on circular replicons only.
``bakta_json`` is a Bakta JSON file, not a TSV, and is parsed by the Bakta parser (M2-03).
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

import pandas as pd

from giap.io.tsv import read_tsv
from giap.model.coords import Interval, Replicon, Strand, check_interval

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


def _replicon(key: str, row: Mapping[str, Any], replicons: Mapping[str, Replicon]) -> Replicon:
    contig = row["contig"]
    if contig not in replicons:
        raise ValueError(f"precomputed {key}: unknown contig {contig!r}")
    return replicons[contig]


def _required_int(key: str, col: str, v: str | None) -> int:
    if v is None:
        raise ValueError(f"precomputed {key}: column {col} may not be NA")
    return int(v)


def _convert_row(
    key: str, schema: Mapping[str, str], row: dict[str, str | None], reps: Mapping[str, Replicon]
) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for col, kind in schema.items():
        if kind in _CONVERT:
            out[col] = _CONVERT[kind](row[col])
    if "start" in schema:
        rep = _replicon(key, row, reps)
        start1 = _required_int(key, "start", row["start"])
        end1 = _required_int(key, "end", row["end"])
        strand: Strand = row.get("strand") or "."  # type: ignore[assignment]
        if start1 < 1:
            raise ValueError(f"precomputed {key}: start must be >= 1, got {start1}")
        iv = Interval(rep.id, start1 - 1, end1, strand)
        check_interval(iv, rep)
        out["start"], out["end"] = iv.start, iv.end
    for col, kind in schema.items():
        if kind == "pos":
            rep = _replicon(key, row, reps)
            pos1 = _required_int(key, col, row[col])
            if not 1 <= pos1 <= rep.length:
                raise ValueError(
                    f"precomputed {key}: {col} {pos1} outside {rep.id} (1..{rep.length})"
                )
            out[col] = pos1 - 1
    if "strand" in schema and out["strand"] not in ("+", "-", "."):
        raise ValueError(f"precomputed {key}: invalid strand {out['strand']!r}")
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
    records = [
        _convert_row(key, schema, dict(zip(header, r, strict=True)), replicons) for r in rows
    ]
    df = pd.DataFrame.from_records(records, columns=list(schema))
    return df.astype({c: _DTYPE[k] for c, k in schema.items() if k in _DTYPE})
