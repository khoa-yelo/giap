"""Variant flags (GIAP_SPEC.md §7.5.1, §7.5.3, §9.4, Appendix B.5;
docs/specs/M7-05-variant-flags.md rules 1-7).

Flags annotate a ``sites`` row; they never change its call, confidence or counts. The mass list is
recomputed at the configured tolerance from the B.5 masses (lead decision 2): the committed
``mass_equivalents.tsv`` is never read here.
"""

from __future__ import annotations

from collections.abc import Collection, Mapping
from typing import Literal

from giap.config.model import GiapConfig
from giap.mapping.masses import DELTA_DECIMALS, MOD_DELTA, RESIDUE_MASS
from giap.mapping.sixframe import translate11_text
from giap.model.genes import GeneModel
from giap.variants.vcf_records import minimal_alt_codon

__all__ = [
    "FLAG_ORDER",
    "M705_FLAGS",
    "flags_for",
    "mass_equivalent",
    "nt_class",
    "site_flags",
    "start_codons",
]

M705_FLAGS: tuple[str, ...] = (
    "mass_equivalent", "snv", "multi_nt", "strain_difference", "database_derived",
    "start_codon", "stop_codon",
)  # fmt: skip
FLAG_ORDER: tuple[str, ...] = (*M705_FLAGS, "vcf_unrepresentable")  # §9.4 order (rule 7)

_StartKey = tuple[str, str, int]  # (contig, strand, wrapped first codon)
_START_TYPES = frozenset({"cds", "sorf"})


def _mass(residues: str) -> float | None:
    """Composition mass (``J`` and ``I`` weigh as ``L``); None if empty or not all residues."""
    folded = residues.replace("J", "L").replace("I", "L")
    if not folded or any(c not in RESIDUE_MASS for c in folded):
        return None
    return sum(RESIDUE_MASS[c] for c in folded)


def mass_equivalent(ref: str, alt: str, tol_da: float) -> bool:
    """Rule 1: ``|m(alt) - (m(ref) + Δ)| <= tol_da`` for Δ in {0} or a B.5 modification delta.

    Directional (the observed ALT is the genome REF plus a modification). Float noise below the
    6-decimal rounding of ``giap.mapping.masses`` counts as 0."""
    m_ref, m_alt = _mass(ref), _mass(alt)
    if m_ref is None or m_alt is None:
        return False
    for delta in (0.0, *MOD_DELTA.values()):
        raw = abs(m_alt - m_ref - delta)
        if min(raw, round(raw, DELTA_DECIMALS)) <= tol_da:
            return True
    return False


def nt_class(ref_codon: str, alt_aa: str) -> Literal["snv", "multi_nt"]:
    """Rule 2: ``snv`` when a table-11 codon at most one base from ``ref_codon`` encodes
    ``alt_aa`` (``J``: any I or L codon), else ``multi_nt``."""
    return "snv" if minimal_alt_codon(ref_codon, alt_aa)[1] <= 1 else "multi_nt"


def start_codons(gm: GeneModel) -> frozenset[_StartKey]:
    """The annotated first codon (transcription direction) of every CDS and sORF, wrapped on
    circular replicons; pseudogenes excluded."""
    out: set[_StartKey] = set()
    for g in gm.genes:
        if g.feature_type not in _START_TYPES:
            continue
        iv, rep = g.interval, gm.replicons[g.interval.contig]
        first = iv.start if iv.strand == "+" else iv.end - 3
        if rep.topology == "circular":
            first %= rep.length
        out.add((iv.contig, iv.strand, first))
    return frozenset(out)


def flags_for(
    site: Mapping[str, object],
    starts: Collection[_StartKey],
    cfg: GiapConfig,
    alt_sources: Collection[tuple[str, str]],
) -> list[str]:
    """:func:`site_flags` with the start codons precomputed (``start_codons``)."""
    kind, alt = str(site["alt_type"]), str(site["alt_aa"])
    ref_aa, ref_nt = str(site["ref_aa"]), str(site["ref_codon"])
    tol = cfg.variants.mass_equivalent_tol_da
    found: set[str] = set()
    if kind == "substitution":
        if mass_equivalent(ref_aa, alt, tol):
            found.add("mass_equivalent")
        found.add(nt_class(ref_nt, alt))
        key = (str(site["contig"]), str(site["strand"]), int(site["codon_start"]))  # type: ignore[call-overload]
        if alt == "M" and key in starts:
            found.add("start_codon")
        if ref_aa == "*":
            found.add("stop_codon")
    elif kind == "block" and mass_equivalent(translate11_text(ref_nt), alt, tol):
        found.add("mass_equivalent")
    if cfg.genome_relation == "related":
        found.add("strain_difference")
    if ("database", "other") in set(alt_sources):
        found.add("database_derived")
    return [f for f in M705_FLAGS if f in found]


def site_flags(
    site: Mapping[str, object],
    gm: GeneModel,
    cfg: GiapConfig,
    alt_sources: Collection[tuple[str, str]],
) -> list[str]:
    """Rules 1-6 for one ``giap.variants.pileup.sites`` row (0-based ``codon_start``).

    ``alt_sources`` holds the ``(source, search_db)`` of the peptides supporting the ALT. Returns
    each flag once in §9.4 order, never ``vcf_unrepresentable`` (M7-03). Swaps and frameshifts
    get rules 3-4 only; blocks also rule 1 on the translation of ``ref_codon``. Pure."""
    return flags_for(site, start_codons(gm), cfg, alt_sources)
