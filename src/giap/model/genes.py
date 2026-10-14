"""Gene model shared by Module 2 and Module 3 (issue M3-01; GIAP_SPEC.md §6.5, §7.1.2).

``GeneModel`` is the API through which every M3 module consumes the annotation; parquet is only
its serialisation. M2's ``parse_bakta`` MUST construct this same type (the contract test lives
with M2 against the recorded Bakta 1.12.2 JSON). ``gene_model_from_bakta_json`` is the minimal
reader for the Bakta-like JSON of giapsim (§11.3).

Coordinates are 0-based half-open; origin-crossing features on circular replicons are stored
unwrapped (``end > length``) and linear contigs never wrap (see ``giap.model.coords``).
"""

from __future__ import annotations

import re
from bisect import bisect_left
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Literal, cast

from giap.io._open import read_json
from giap.model.coords import (
    Interval,
    Replicon,
    Strand,
    from_one_based_closed,
    interval_sort_key,
    split_at_origin,
)

if TYPE_CHECKING:  # ``giap.model`` does no file I/O itself; reading goes through ``giap.io``.
    from _typeshed import StrPath

GeneType = Literal["cds", "sorf", "pseudogene"]
NcType = Literal["rrna", "trna", "tmrna", "ncrna", "ncrna_region", "crispr", "oric", "orit", "gap"]

_NC_TYPES: dict[str, NcType] = {
    "rrna": "rrna",
    "trna": "trna",
    "tmrna": "tmrna",
    "ncrna": "ncrna",
    "ncrna-region": "ncrna_region",
    "crispr": "crispr",
    "oric": "oric",
    "orit": "orit",
    "gap": "gap",
}
_HYPOTHETICAL = re.compile(r"^hypothetical protein")


@dataclass(frozen=True, slots=True)
class Gene:
    """A CDS, sORF or pseudogene; ``crosses_origin`` is ``interval.crosses_origin(replicon)``."""

    locus_tag: str
    feature_type: GeneType
    interval: Interval
    aa: str
    nt: str | None
    start_codon: str | None
    gene: str | None = None
    product: str = ""
    partial: bool = False
    hypothetical: bool = False
    db_xrefs: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class NcFeature:
    """A non-coding feature (RNA genes, CRISPR arrays, origins, gaps)."""

    id: str
    type: NcType
    interval: Interval


_Item = Gene | NcFeature


@dataclass(frozen=True)
class _ContigIndex:
    starts: list[int]  # sorted segment starts
    segments: list[tuple[int, int, int]]  # (start, end, item index), sorted
    max_len: int


@dataclass(frozen=True)
class GeneModel:
    """Genes (sorted by ``interval_sort_key`` then locus tag) and non-coding features."""

    replicons: Mapping[str, Replicon]
    genes: tuple[Gene, ...]
    features: tuple[NcFeature, ...]
    _by_locus: dict[str, Gene] = field(init=False, repr=False, compare=False)
    _index: dict[str, _ContigIndex] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        by_locus: dict[str, Gene] = {}
        for g in self.genes:
            if g.locus_tag in by_locus:
                raise ValueError(f"duplicate locus tag {g.locus_tag!r}")
            by_locus[g.locus_tag] = g
        segs: dict[str, list[tuple[int, int, int]]] = {}
        for k, item in enumerate(self._items()):
            rep = self._replicon(item.interval.contig)
            for part in split_at_origin(item.interval, rep):
                segs.setdefault(part.contig, []).append((part.start, part.end, k))
        index: dict[str, _ContigIndex] = {}
        for contig, s in segs.items():
            s.sort()
            index[contig] = _ContigIndex(
                [x[0] for x in s], s, max(end - start for start, end, _ in s)
            )
        object.__setattr__(self, "_by_locus", by_locus)
        object.__setattr__(self, "_index", index)

    def _items(self) -> tuple[_Item, ...]:
        return (*self.genes, *self.features)

    def _replicon(self, contig: str) -> Replicon:
        try:
            return self.replicons[contig]
        except KeyError:
            raise ValueError(f"unknown contig {contig!r} (not among the replicons)") from None

    def by_locus(self, locus_tag: str) -> Gene:
        """Return the gene with this locus tag; ``KeyError`` if absent."""
        return self._by_locus[locus_tag]

    def overlapping(self, iv: Interval) -> list[Gene | NcFeature]:
        """Genes and features overlapping ``iv`` on either strand, origin-aware, in model order."""
        idx = self._index.get(iv.contig)
        if idx is None:
            return []
        hits: set[int] = set()
        for q in split_at_origin(iv, self._replicon(iv.contig)):
            lo = bisect_left(idx.starts, q.start - idx.max_len + 1)
            hi = bisect_left(idx.starts, q.end)
            hits.update(k for _, end, k in idx.segments[lo:hi] if end > q.start)
        items = self._items()
        return [items[k] for k in sorted(hits)]


def _strip_stop(aa: str) -> str:
    return aa[:-1] if aa.endswith("*") else aa


def _parse_interval(f: Mapping[str, Any], name: str, replicons: Mapping[str, Replicon]) -> Interval:
    contig = f.get("sequence", f.get("contig"))
    if contig not in replicons:
        raise ValueError(f"feature {name!r}: unknown contig {contig!r}")
    strand = cast(Strand, f.get("strand", "."))
    try:
        return from_one_based_closed(
            contig, int(f["start"]), int(f["stop"]), strand, replicons[contig]
        )
    except (ValueError, TypeError, KeyError) as exc:
        raise ValueError(f"feature {name!r}: {exc}") from exc


def _parse_gene(f: Mapping[str, Any], replicons: Mapping[str, Replicon]) -> Gene:
    locus = str(f["locus"])
    ftype: GeneType = cast(GeneType, str(f["type"]).lower())
    if ftype == "cds" and (f.get("pseudo") or "pseudogene" in f):
        ftype = "pseudogene"
    product = f.get("product") or ""
    hypothetical = (
        bool(f["hypothetical"]) if "hypothetical" in f else bool(_HYPOTHETICAL.match(product))
    )
    return Gene(
        locus_tag=locus,
        feature_type=ftype,
        interval=_parse_interval(f, locus, replicons),
        aa=_strip_stop(f.get("aa") or ""),
        nt=f.get("nt"),
        start_codon=f.get("start_type"),
        gene=f.get("gene"),
        product=product,
        partial="truncated" in f,
        hypothetical=hypothetical,
        db_xrefs=tuple(f.get("db_xrefs") or ()),
    )


def gene_model_from_bakta_json(path: StrPath, replicons: Mapping[str, Replicon]) -> GeneModel:
    """Read a Bakta-like JSON (§11.3); unknown feature types are skipped (no Issue returned)."""
    data = read_json(path)
    order = {c: k for k, c in enumerate(replicons)}
    genes: list[Gene] = []
    features: list[NcFeature] = []
    for f in data.get("features", []):
        t = str(f.get("type", "")).lower()
        if t in ("cds", "sorf"):
            genes.append(_parse_gene(f, replicons))
        elif t in _NC_TYPES:
            fid = str(f.get("id", ""))
            features.append(NcFeature(fid, _NC_TYPES[t], _parse_interval(f, fid, replicons)))
    genes.sort(key=lambda g: interval_sort_key(g.interval, order, g.locus_tag))
    features.sort(key=lambda x: interval_sort_key(x.interval, order, x.id))
    return GeneModel(replicons=replicons, genes=tuple(genes), features=tuple(features))
