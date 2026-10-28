"""Gene model shared by Module 2 and Module 3 (issue M3-01; GIPPA_SPEC.md §6.5, §7.1.2).

``GeneModel`` is the API through which every M3 module consumes the annotation; parquet is only
its serialisation. M2's ``parse_bakta`` MUST construct this same type (the contract test lives
with M2 against the recorded Bakta 1.12.2 JSON). ``gene_model_from_bakta_json`` delegates to
``gippa.annotate.bakta.parse_bakta_json`` (M2-03 review decision 3).

Coordinates are 0-based half-open; origin-crossing features on circular replicons are stored
unwrapped (``end > length``) and linear contigs never wrap (see ``gippa.model.coords``).
"""

from __future__ import annotations

from bisect import bisect_left
from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING, Literal

from gippa.model.coords import Interval, Replicon, split_at_origin

if TYPE_CHECKING:  # ``gippa.model`` does no file I/O itself; reading goes through ``gippa.io``.
    from _typeshed import StrPath

GeneType = Literal["cds", "sorf", "pseudogene"]
NcType = Literal[
    "rrna", "trna", "tmrna", "ncrna", "ncrna_region", "crispr", "oric", "orit", "oriv", "gap"
]

_NC_TYPES: dict[str, NcType] = {
    "rrna": "rrna",
    "trna": "trna",
    "tmrna": "tmrna",
    "ncrna": "ncrna",
    "ncrna-region": "ncrna_region",
    "crispr": "crispr",
    "oric": "oric",
    "orit": "orit",
    "oriv": "oriv",
    "gap": "gap",
}


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
    expert: tuple[str, ...] = ()  # expert-system sources, e.g. ("amrfinder",) (M2-03)


@dataclass(frozen=True, slots=True)
class NcFeature:
    """A non-coding feature (RNA genes, CRISPR arrays, origins, gaps)."""

    id: str
    type: NcType
    interval: Interval
    gene: str | None = None
    product: str = ""
    amino_acid: str | None = None  # tRNA
    anticodon: str | None = None  # tRNA, DNA alphabet, upper case
    nt: str | None = None  # rRNA (16S -> anti-SD)


_Item = Gene | NcFeature


@dataclass(frozen=True, slots=True)
class _Bucket:
    """Segments of similar length (same ``bit_length``), so the bisect window stays tight."""

    starts: list[int]  # sorted segment starts
    segments: list[tuple[int, int, int]]  # (start, end, item index), sorted
    max_len: int


def _buckets(segs: list[tuple[int, int, int]]) -> tuple[_Bucket, ...]:
    by_class: dict[int, list[tuple[int, int, int]]] = {}
    for seg in segs:
        by_class.setdefault((seg[1] - seg[0]).bit_length(), []).append(seg)
    out = []
    for _, group in sorted(by_class.items()):
        group.sort()
        out.append(_Bucket([x[0] for x in group], group, max(e - s for s, e, _ in group)))
    return tuple(out)


@dataclass(frozen=True, eq=False)
class GeneModel:
    """Genes (sorted by ``interval_sort_key`` then locus tag) and non-coding features.

    ``replicons`` is copied into a read-only mapping (the caller's dict is not aliased).
    Equality compares replicons, genes and features (not the derived indexes); instances are
    unhashable (``__hash__ = None``).
    """

    __hash__ = None  # type: ignore[assignment]

    replicons: Mapping[str, Replicon]
    genes: tuple[Gene, ...]
    features: tuple[NcFeature, ...]
    _by_locus: dict[str, Gene] = field(init=False, repr=False, compare=False)
    _index: dict[str, tuple[_Bucket, ...]] = field(init=False, repr=False, compare=False)
    _all: tuple[_Item, ...] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "replicons", MappingProxyType(dict(self.replicons)))
        object.__setattr__(self, "_all", (*self.genes, *self.features))
        by_locus: dict[str, Gene] = {}
        for g in self.genes:
            if g.locus_tag in by_locus:
                raise ValueError(f"duplicate locus tag {g.locus_tag!r}")
            by_locus[g.locus_tag] = g
        segs: dict[str, list[tuple[int, int, int]]] = {}
        for k, item in enumerate(self._all):
            rep = self._replicon(item.interval.contig)
            for part in split_at_origin(item.interval, rep):
                segs.setdefault(part.contig, []).append((part.start, part.end, k))
        object.__setattr__(self, "_by_locus", by_locus)
        object.__setattr__(self, "_index", {c: _buckets(s) for c, s in segs.items()})

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, GeneModel):
            return NotImplemented
        return (
            dict(self.replicons) == dict(other.replicons)
            and self.genes == other.genes
            and self.features == other.features
        )

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
        buckets = self._index.get(iv.contig)
        if buckets is None:
            return []
        hits: set[int] = set()
        for q in split_at_origin(iv, self._replicon(iv.contig)):
            qs, qe = q.start, q.end
            for b in buckets:
                lo = bisect_left(b.starts, qs - b.max_len + 1)
                hi = bisect_left(b.starts, qe, lo)
                hits.update(k for _, end, k in b.segments[lo:hi] if end > qs)
        items = self._all
        return [items[k] for k in sorted(hits)]


def gene_model_from_bakta_json(path: StrPath, replicons: Mapping[str, Replicon]) -> GeneModel:
    """Read a Bakta(-like) JSON (§11.3) via M2-03's parser: one parser, identical results.

    Unknown feature types are skipped; malformed input raises ``BaktaParseError``
    (a ``ValueError``).
    """
    from gippa.annotate.bakta import parse_bakta_json  # lazy: gippa.annotate builds on gippa.model

    return parse_bakta_json(path, replicons).model
