"""antiSMASH 8.0.4 JSON parser: BGC regions (GIAP_SPEC.md §7.1.1, §7.1.7; M2-18).

``records[].id`` = contig; ``records[].areas[]`` ``start``/``end`` are 0-based half-open;
``start > end`` on a circular record is origin-crossing and stored unwrapped (``end + length``),
on a linear record it is corrupt input. Products keep antiSMASH order, duplicates collapsed.
Core genes: ``locus_tag`` of CDS features wholly inside the area whose ``gene_kind`` qualifier
contains exactly ``biosynthetic`` (M2-18 lead decision 1), in feature order.
"""

from __future__ import annotations

import json
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from giap.model.coords import Interval, Replicon, check_interval

# Biopython location parts: ``[a:b](±)``, fuzzy ``[<a:b]`` / ``[a:>b]``, inside ``join{…}``.
_PART = re.compile(r"\[<?(\d+):>?(\d+)\]")


@dataclass(frozen=True, slots=True)
class BgcRegion:
    contig: str
    interval: Interval
    products: tuple[str, ...]
    core_genes: tuple[str, ...]


def _area_interval(area: Mapping[str, Any], rep: Replicon) -> Interval:
    start, end = int(area["start"]), int(area["end"])
    if start > end:
        if rep.topology != "circular":
            raise ValueError(f"antiSMASH: area {rep.id}:{start}-{end} wraps on a linear contig")
        end += rep.length
    iv = Interval(rep.id, start, end, ".")
    check_interval(iv, rep)
    return iv


def _inside(location: str, area: Interval, length: int) -> bool:
    parts = [(int(a), int(b)) for a, b in _PART.findall(location)]
    if not parts:
        raise ValueError(f"antiSMASH: cannot parse CDS location {location!r}")
    return all(
        any(area.start <= a + s and b + s <= area.end for s in (0, length)) for a, b in parts
    )


def _core_genes(features: list[Mapping[str, Any]], area: Interval, length: int) -> tuple[str, ...]:
    tags: dict[str, None] = {}
    for f in features:
        q = f.get("qualifiers", {})
        if f.get("type") != "CDS" or "biosynthetic" not in q.get("gene_kind", []):
            continue
        if _inside(f["location"], area, length):
            tags.update(dict.fromkeys(q.get("locus_tag", [])))
    return tuple(tags)


def parse_antismash_json(path: Path, replicons: Mapping[str, Replicon]) -> list[BgcRegion]:
    """BGC regions of an antiSMASH JSON, in record then area order."""
    data = json.loads(path.read_text())
    out: list[BgcRegion] = []
    for rec in data.get("records", []):
        rid = rec["id"]
        if rid not in replicons:
            raise ValueError(f"antiSMASH: record {rid!r} not in the genome")
        rep = replicons[rid]
        for area in rec.get("areas", []):
            iv = _area_interval(area, rep)
            products = tuple(dict.fromkeys(area.get("products", [])))
            out.append(
                BgcRegion(rid, iv, products, _core_genes(rec.get("features", []), iv, rep.length))
            )
    return out
