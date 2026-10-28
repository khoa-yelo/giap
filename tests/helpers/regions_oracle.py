"""Independent oracles and fixture loaders for the M2-19 region tests (GIPPA_SPEC.md §7.1.7).

Nothing here imports ``gippa.annotate.regions``: overlaps are counted base by base, subsets use
a separately written keyword table, and the recorded inputs are read through the merged
producers (M2-03 ``parse_bakta_json``, M2-18 ``parse_genomad`` / ``parse_antismash_json``).
"""

from __future__ import annotations

import json
import os
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import cast

from gippa.annotate.antismash import BgcRegion, parse_antismash_json
from gippa.annotate.bakta import parse_bakta_json
from gippa.annotate.genomad import GenomadResult, VirusGene, parse_genomad
from gippa.model.coords import Interval, Replicon, RepliconKind, Topology
from gippa.model.genes import Gene, NcFeature

TESTS = Path(__file__).resolve().parents[1]
RECORDED = TESTS / "data" / "recorded"
EXCERPTS = RECORDED / "regions"
GENOMAD_REC = RECORDED / "genomad" / "1.12.0"
ANTISMASH_REC = RECORDED / "antismash" / "8.0.4"
FULL = Path(os.environ.get("GIPPA_ROOT", "/nonexistent")) / "recordings"

TSV_HEADER = [
    "region_id", "region_type", "contig", "start", "end", "source_tool", "score", "genes", "notes",
]  # fmt: skip


# ------------------------------------------------------------------------------- base counting


def bases(iv: Interval, rep: Replicon) -> set[int]:
    """Positions covered by ``iv`` folded into ``[0, length)`` (origin-aware)."""
    return {p % rep.length for p in range(iv.start, iv.end)}


def brute_overlap_frac(gene: Interval, region: Interval, rep: Replicon) -> float:
    if gene.contig != region.contig:
        return 0.0
    g = bases(gene, rep)
    return len(g & bases(region, rep)) / len(g)


def brute_members(
    genes: Iterable[Gene], region: Interval, rep: Replicon, min_frac: float = 0.5
) -> dict[str, float]:
    reg = bases(region, rep)  # counted once per region
    out: dict[str, float] = {}
    for g in genes:
        if g.interval.contig != region.contig:
            continue
        gb = bases(g.interval, rep)
        f = len(gb & reg) / len(gb)
        if f >= min_frac:
            out[g.locus_tag] = f
    return out


# ------------------------------------------------------------------------------- subsets

_CI_RULES: list[tuple[str, tuple[str, ...]]] = [
    ("integrase", ("integrase", "recombinase", "excisionase")),
    ("repressor", ("repressor", "antirepressor")),
    ("structural", ("capsid", "portal", "tail", "baseplate", "head", "terminase", "scaffold",
                    "tape measure", "fiber")),
    ("lysis", ("holin", "lysin", "lysozyme", "endolysin", "spanin", "rz")),
]  # fmt: skip
_CS_REPRESSOR = ("cI", "Cro")


def _word(text: str, kw: str, flags: int) -> bool:
    return re.search(rf"(?<![A-Za-z0-9]){re.escape(kw)}(?![A-Za-z0-9])", text, flags) is not None


def oracle_subset(product: str, marker_description: str | None) -> str:
    text = f"{product} {marker_description or ''}"
    for name, kws in _CI_RULES:
        if any(_word(text, k, re.IGNORECASE) for k in kws):
            return name
        if name == "repressor" and any(_word(text, k, 0) for k in _CS_REPRESSOR):
            return name
    return "cargo"


# ------------------------------------------------------------------------------- CRISPR / AMR


def is_cas(g: Gene) -> bool:
    return bool(re.match(r"cas", g.gene or "", re.IGNORECASE)) or ("CRISPR-associated" in g.product)


def cas_within(genes: Iterable[Gene], array: Interval, max_gap: int = 1000) -> list[Gene]:
    out = []
    for g in genes:
        iv = g.interval
        if iv.contig != array.contig or not is_cas(g):
            continue
        gap = max(array.start - iv.end, iv.start - array.end, 0)
        if gap <= max_gap:
            out.append(g)
    return out


def is_amr(g: Gene) -> bool:
    return any("amrfinder" in e.lower() for e in g.expert)


# ------------------------------------------------------------------------------- loaders


@dataclass(frozen=True)
class Sample:
    name: str
    replicons: dict[str, Replicon]
    genes: tuple[Gene, ...]
    features: tuple[NcFeature, ...]
    genomad: GenomadResult
    bgc: list[BgcRegion]


def _kind(raw: str | None) -> RepliconKind:
    k = (raw or "").lower()
    return cast(RepliconKind, k if k in ("chromosome", "plasmid") else "unknown")


def excerpt_replicons(sample: str) -> dict[str, Replicon]:
    data = json.loads((EXCERPTS / f"{sample}.bakta_excerpt.json").read_text())
    return {
        q["id"]: Replicon(q["id"], q["length"], cast(Topology, q["topology"]), _kind(q["type"]))
        for q in data["sequences"]
    }


def load_excerpt(sample: str) -> Sample:
    """Committed Bakta excerpt + complete committed geNOMAD summaries + trimmed antiSMASH JSON."""
    reps = excerpt_replicons(sample)
    model = parse_bakta_json(EXCERPTS / f"{sample}.bakta_excerpt.json", replicons=reps).model
    gm = parse_genomad(GENOMAD_REC / sample, "genome", reps)
    bgc = parse_antismash_json(ANTISMASH_REC / sample / "bakta.json", reps)
    return Sample(sample, reps, model.genes, model.features, gm, bgc)


def full_available(sample: str) -> bool:
    return all(
        p.exists()
        for p in (
            FULL / "bakta_full" / sample / "bakta" / f"{sample}.json",
            FULL / "genomad" / sample / "genomad" / "genome_summary",
            FULL / "antismash" / sample / "antismash" / "bakta.json",
        )
    )


def load_full(sample: str) -> Sample:
    """Full recordings under ``$GIPPA_ROOT/recordings`` (Bakta full DB, geNOMAD, antiSMASH)."""
    ann = parse_bakta_json(FULL / "bakta_full" / sample / "bakta" / f"{sample}.json")
    reps = dict(ann.model.replicons)
    gm = parse_genomad(FULL / "genomad" / sample / "genomad", "genome", reps)
    bgc = parse_antismash_json(FULL / "antismash" / sample / "antismash" / "bakta.json", reps)
    return Sample(sample, reps, ann.model.genes, ann.model.features, gm, bgc)


def exact_virus_gene(g: Gene, vgenes: Sequence[VirusGene]) -> VirusGene | None:
    for v in vgenes:
        if v.interval == g.interval:
            return v
    return None


def overlapping_virus_genes(g: Gene, vgenes: Sequence[VirusGene]) -> list[VirusGene]:
    iv = g.interval
    return [
        v
        for v in vgenes
        if v.contig == iv.contig and v.interval.start < iv.end and v.interval.end > iv.start
    ]


def contig_rank(reps: Mapping[str, Replicon]) -> dict[str, int]:
    return {c: i for i, c in enumerate(reps)}
