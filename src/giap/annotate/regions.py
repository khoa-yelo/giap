"""Regions: prophage, plasmid, BGC, CRISPR and AMR (GIAP_SPEC.md §5.2, §6.5, §7.1.7, §9.3; M2-19).

A gene belongs to a region when at least ``annotation.region_min_overlap_frac`` of its bases lie
inside it (origin-aware, ``>=``). Regions are numbered per type by contig order, then start
(``<type>_<n>``, §6.5). Prophage members get a function subset by keyword rules on the Bakta
product plus the geNOMAD marker description of the virus gene with the identical interval; BGC
members carry antiSMASH's core-biosynthetic flag. AMR regions are gene sets (one per contig,
``interval = None``) and are left out of the BED file.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from giap.annotate.antismash import BgcRegion
from giap.annotate.genomad import GenomadResult
from giap.config import GiapConfig
from giap.config.model import RegionType
from giap.io.bed import BedRecord
from giap.io.tsv import write_tsv
from giap.model.coords import Interval, Replicon, interval_sort_key, split_at_origin
from giap.model.genes import Gene, NcFeature

__all__ = [
    "TSV_COLUMNS",
    "Region",
    "RegionGene",
    "RegionType",
    "build_regions",
    "members",
    "overlap_frac",
    "prophage_subset",
    "regions_bed",
    "write_regions_tsv",
]

TSV_COLUMNS = (
    "region_id", "region_type", "contig", "start", "end", "source_tool", "score", "genes", "notes",
)  # fmt: skip
CAS_MAX_GAP = 1000
_TYPE_ORDER: tuple[RegionType, ...] = ("prophage", "plasmid", "bgc", "crispr", "amr")


@dataclass(frozen=True, slots=True)
class Region:
    region_id: str
    region_type: RegionType
    contig: str
    interval: Interval | None  # None = gene set (amr)
    source_tool: str
    score: float | None
    genes: tuple[str, ...]
    notes: str


@dataclass(frozen=True, slots=True)
class RegionGene:
    region_id: str
    locus_tag: str
    overlap_frac: float
    subset: str | None  # prophage only: integrase|repressor|structural|lysis|cargo
    core_biosynthetic: bool  # bgc only


# ------------------------------------------------------------------------------- membership


def _pieces(iv: Interval, rep: Replicon) -> tuple[Interval, ...]:
    return split_at_origin(iv, rep) if rep.topology == "circular" else (iv,)


def overlap_frac(gene: Interval, region: Interval, replicon: Replicon) -> float:
    """Fraction of ``gene``'s bases inside ``region`` (both may cross the origin)."""
    if gene.contig != region.contig:
        return 0.0
    inside = sum(
        max(0, min(g.end, r.end) - max(g.start, r.start))
        for g in _pieces(gene, replicon)
        for r in _pieces(region, replicon)
    )
    return inside / gene.length_bp


def members(
    genes: Sequence[Gene], region: Interval, replicon: Replicon, min_frac: float = 0.5
) -> list[tuple[str, float]]:
    """``(locus_tag, overlap_frac)`` of the genes with ``overlap_frac >= min_frac``, input order."""
    out = []
    for g in genes:
        if g.interval.contig == region.contig:
            f = overlap_frac(g.interval, region, replicon)
            if f >= min_frac:
                out.append((g.locus_tag, f))
    return out


# ------------------------------------------------------------------------------- subsets


def _words(kws: str, flags: int = re.IGNORECASE) -> re.Pattern[str]:
    alt = "|".join(re.escape(k) for k in kws.split(","))
    return re.compile(rf"(?<![A-Za-z0-9])(?:{alt})(?![A-Za-z0-9])", flags)


# Checked in this order; ``cI`` and ``Cro`` are case-sensitive (not ``micro…``, ``circ…``).
_SUBSET_RULES: tuple[tuple[str, tuple[re.Pattern[str], ...]], ...] = (
    ("integrase", (_words("integrase,recombinase,excisionase"),)),
    ("repressor", (_words("repressor,antirepressor"), _words("cI,Cro", 0))),
    (
        "structural",
        (_words("capsid,portal,tail,baseplate,head,terminase,scaffold,tape measure,fiber"),),
    ),
    ("lysis", (_words("holin,lysin,lysozyme,endolysin,spanin,rz"),)),
)


def prophage_subset(product: str, marker_description: str | None) -> str:
    """Function subset of a prophage gene; first matching rule wins, else ``cargo``."""
    text = f"{product} {marker_description or ''}"
    for name, pats in _SUBSET_RULES:
        if any(p.search(text) for p in pats):
            return name
    return "cargo"


# ------------------------------------------------------------------------------- building

_Draft = tuple[str, Interval | None, str, float | None, list[tuple[str, float]], str]
# (contig, interval, source_tool, score, members, notes)


def _gene_key(g: Gene) -> tuple[int, int, int, int, str]:
    return interval_sort_key(g.interval, {g.interval.contig: 0}, g.locus_tag)


def _whole(rep: Replicon) -> Interval:
    return Interval(rep.id, 0, rep.length, ".")


def _prophages(
    gm: GenomadResult, genes: Sequence[Gene], reps: Mapping[str, Replicon], frac: float
) -> list[_Draft]:
    out: list[_Draft] = []
    for p in gm.proviruses:
        rep = reps[p.contig]
        iv, notes = p.interval, p.taxonomy
        if iv is None:  # whole contig classified viral
            iv, notes = _whole(rep), ",".join(x for x in ("whole_contig", p.taxonomy) if x)
        out.append((p.contig, iv, "geNOMAD", p.virus_score, members(genes, iv, rep, frac), notes))
    return out


def _plasmids(
    gm: GenomadResult | None, genes: Sequence[Gene], reps: Mapping[str, Replicon], frac: float
) -> list[_Draft]:
    calls = {c.contig: c.plasmid_score for c in (gm.plasmids if gm else ())}
    out: list[_Draft] = []
    for rep in reps.values():
        tools = [
            t for t, ok in (("Bakta", rep.kind == "plasmid"), ("geNOMAD", rep.id in calls)) if ok
        ]
        if tools:
            iv = _whole(rep)
            score = calls.get(rep.id)
            out.append((rep.id, iv, "+".join(tools), score, members(genes, iv, rep, frac), ""))
    return out


def _bgcs(
    bgcs: Sequence[BgcRegion], genes: Sequence[Gene], reps: Mapping[str, Replicon], frac: float
) -> list[_Draft]:
    return [
        (b.contig, b.interval, "antiSMASH", None,
         members(genes, b.interval, reps[b.contig], frac), "products=" + ",".join(b.products))
        for b in bgcs
    ]  # fmt: skip


def _is_cas(g: Gene) -> bool:
    return bool(re.match("cas", g.gene or "", re.IGNORECASE)) or "CRISPR-associated" in g.product


def _crisprs(
    features: Sequence[NcFeature], genes: Sequence[Gene], reps: Mapping[str, Replicon], frac: float
) -> list[_Draft]:
    out: list[_Draft] = []
    for f in features:
        if f.type != "crispr":
            continue
        a = f.interval
        cas = [
            g.interval
            for g in genes
            if g.interval.contig == a.contig
            and _is_cas(g)
            and max(a.start - g.interval.end, g.interval.start - a.end, 0) <= CAS_MAX_GAP
        ]
        start = min([a.start, *(c.start for c in cas)])
        end = max([a.end, *(c.end for c in cas)])
        iv = Interval(a.contig, start, end, a.strand)
        out.append((a.contig, iv, "Bakta", None, members(genes, iv, reps[a.contig], frac), ""))
    return out


def _amr(genes: Sequence[Gene], reps: Mapping[str, Replicon]) -> list[_Draft]:
    by_contig: dict[str, list[tuple[str, float]]] = {}
    for g in genes:
        if any("amrfinder" in e.lower() for e in g.expert):
            by_contig.setdefault(g.interval.contig, []).append((g.locus_tag, 1.0))
    return [(c, None, "AMRFinderPlus", None, by_contig[c], "") for c in reps if c in by_contig]


def build_regions(
    genes: Sequence[Gene],
    features: Sequence[NcFeature],
    replicons: Mapping[str, Replicon],
    genomad: GenomadResult | None,
    antismash: Sequence[BgcRegion] | None,
    cfg: GiapConfig,
) -> tuple[list[Region], list[RegionGene]]:
    """Regions of every type in ``context.region_types`` and their member genes (§7.1.7)."""
    frac = cfg.annotation.region_min_overlap_frac
    wanted = set(cfg.context.region_types)
    ordered = sorted(genes, key=lambda g: (g.interval.contig, _gene_key(g)))
    by_tag = {g.locus_tag: g for g in ordered}
    drafts: dict[RegionType, list[_Draft]] = {
        "prophage": _prophages(genomad, ordered, replicons, frac) if genomad else [],
        "plasmid": _plasmids(genomad, ordered, replicons, frac),
        "bgc": _bgcs(antismash or (), ordered, replicons, frac),
        "crispr": _crisprs(features, ordered, replicons, frac),
        "amr": _amr(ordered, replicons),
    }
    rank = {c: i for i, c in enumerate(replicons)}
    desc: dict[Interval, str | None] = {}
    for v in genomad.virus_genes if genomad else ():
        desc.setdefault(v.interval, v.description)
    core = {(b.contig, b.interval): set(b.core_genes) for b in antismash or ()}
    regions: list[Region] = []
    rgenes: list[RegionGene] = []
    for rtype in _TYPE_ORDER:
        if rtype not in wanted:
            continue
        ds = sorted(
            drafts[rtype],
            key=lambda d: (rank[d[0]], d[1].start if d[1] else 0, d[1].end if d[1] else 0),
        )
        for n, (contig, iv, tool, score, mem, notes) in enumerate(ds, start=1):
            rid = f"{rtype}_{n}"
            regions.append(
                Region(rid, rtype, contig, iv, tool, score, tuple(t for t, _ in mem), notes)
            )
            cores = core.get((contig, iv), set()) if iv is not None else set()
            for tag, f in mem:
                g = by_tag[tag]
                subset = (
                    prophage_subset(g.product, desc.get(g.interval))
                    if rtype == "prophage"
                    else None
                )
                rgenes.append(RegionGene(rid, tag, f, subset, rtype == "bgc" and tag in cores))
    return regions, rgenes


# ------------------------------------------------------------------------------- writers


def write_regions_tsv(path: Path, regions: Sequence[Region]) -> None:
    """``annotation/regions.tsv`` (§9.3): 1-based closed coordinates, ``NA`` for AMR gene sets."""
    rows = [
        (
            r.region_id, r.region_type, r.contig,
            r.interval.start + 1 if r.interval else None, r.interval.end if r.interval else None,
            r.source_tool, r.score, list(r.genes), r.notes,
        )
        for r in regions
    ]  # fmt: skip
    write_tsv(path, TSV_COLUMNS, rows)


def regions_bed(regions: Sequence[Region], replicons: Mapping[str, Replicon]) -> list[BedRecord]:
    """One BED record per region with an interval (AMR gene sets are left out)."""
    return [
        BedRecord(r.interval, r.region_id)
        for r in regions
        if r.interval is not None and r.interval.contig in replicons
    ]
