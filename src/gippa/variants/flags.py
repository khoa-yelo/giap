"""Variant flags (GIPPA_SPEC.md §7.5.1, §7.5.3, §9.4, Appendix B.5;
docs/specs/M7-05-variant-flags.md rules 1-7).

Flags annotate a ``sites`` row; they never change its call, confidence or counts. The mass list is
recomputed at the configured tolerance from the B.5 masses (lead decision 2): the committed
``mass_equivalents.tsv`` is never read here. ADR-0049: every swap row and every substitution that an
ALT peptide also reads as a left 2->1 block (:func:`left_block_read`) is ``mass_equivalent``.
"""

from __future__ import annotations

from collections.abc import Collection, Iterable, Mapping
from typing import Literal

from gippa.config.model import GippaConfig
from gippa.mapping.masses import DELTA_DECIMALS, MOD_DELTA, RESIDUE_MASS
from gippa.mapping.sixframe import translate_text
from gippa.model.coords import Replicon
from gippa.model.genes import GeneModel
from gippa.variants.align import HitAlignment, translate_locus
from gippa.variants.vcf_records import minimal_alt_codon

__all__ = [
    "FLAG_ORDER",
    "M705_FLAGS",
    "flags_for",
    "identifiable",
    "left_block_read",
    "left_block_site",
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
    6-decimal rounding of ``gippa.mapping.masses`` counts as 0."""
    m_ref, m_alt = _mass(ref), _mass(alt)
    if m_ref is None or m_alt is None:
        return False
    for delta in (0.0, *MOD_DELTA.values()):
        raw = abs(m_alt - m_ref - delta)
        if min(raw, round(raw, DELTA_DECIMALS)) <= tol_da:
            return True
    return False


def identifiable(ref_aa: str, alt_aa: str, tol_da: float) -> bool:
    """ADR-0049 decision 4: a truth row counts in the G1 recall denominator iff its ALT is not
    :func:`mass_equivalent` to its REF. Pure."""
    return not mass_equivalent(ref_aa, alt_aa, tol_da)


def _fold_il(residues: str) -> str:
    return residues.replace("I", "L")


def left_block_read(peptide: str, q: int, ref: str, k: int, tol_da: float) -> bool:
    """ADR-0049 A6a: peptide residue ``q``, read as a substitution of ``ref[k]``, also aligns as the
    2->1 block ``ref[k-1:k+1] -> peptide[q]`` with its prefix shifted by one: ``q == 0`` or
    ``peptide[:q] == ref[k-1-q:k-1]`` (I and L folded). No right-hand or looser form. False when
    the pair or the shifted prefix runs off ``ref``. Pure."""
    if k < 1 or k >= len(ref) or not 0 <= q < len(peptide):
        return False
    if not mass_equivalent(ref[k - 1 : k + 1], peptide[q], tol_da):
        return False
    lo = k - 1 - q
    return q == 0 or (lo >= 0 and _fold_il(peptide[:q]) == _fold_il(ref[lo : k - 1]))


def _left_window(c: int, q: int, strand: str, rep: Replicon) -> tuple[int, int] | None:
    """``[start, end)`` of the ``q + 2`` codons ending (transcription direction) at codon ``c``;
    None when it runs off a linear contig."""
    if strand == "+":
        s, e = c - 3 * (q + 1), c + 3
        if s < 0:
            if rep.topology != "circular":
                return None
            s, e = s + rep.length, e + rep.length
    else:
        s, e = c, c + 3 * (q + 2)
        if e > rep.length and rep.topology != "circular":
            return None
    return (s, e) if e - s <= rep.length else None


def left_block_site(
    site: Mapping[str, object],
    alignments: Iterable[HitAlignment],
    genome: Mapping[str, str],
    replicons: Mapping[str, Replicon],
    tol_da: float,
    *,
    table: int = 11,
) -> bool:
    """True iff ``site`` is a substitution and one of ``alignments`` (those of its ALT peptides)
    has a substitution op at the site codon that :func:`left_block_read` accepts against the
    reference translated (genetic code ``table``) over the codon and the ``q + 1`` codons before
    it. Pure."""
    if site["alt_type"] != "substitution":
        return False
    contig, strand = str(site["contig"]), str(site["strand"])
    c = int(site["codon_start"])  # type: ignore[call-overload]
    rep = replicons[contig]
    for a in alignments:
        if a.contig != contig or a.strand != strand:
            continue
        for op in a.ops:
            if op.kind != "substitution" or not op.codons:
                continue
            oc = op.codons[0] % rep.length if rep.topology == "circular" else op.codons[0]
            q = op.pep[0]
            win = _left_window(c, q, strand, rep) if oc == c else None
            if win is None:
                continue
            ref = translate_locus(genome[contig], rep, win[0], win[1], strand, table=table)
            if left_block_read(a.peptide, q, ref, q + 1, tol_da):
                return True
    return False


def nt_class(ref_codon: str, alt_aa: str, *, table: int = 11) -> Literal["snv", "multi_nt"]:
    """Rule 2: ``snv`` when a codon of genetic code ``table`` at most one base from ``ref_codon``
    encodes ``alt_aa`` (``J``: any I or L codon), else ``multi_nt``."""
    return "snv" if minimal_alt_codon(ref_codon, alt_aa, table=table)[1] <= 1 else "multi_nt"


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
    cfg: GippaConfig,
    alt_sources: Collection[tuple[str, str]],
    *,
    left_block: bool = False,
) -> list[str]:
    """:func:`site_flags` with the start codons precomputed (``start_codons``). ``left_block``
    (:func:`left_block_site`) flags a substitution ``mass_equivalent`` (ADR-0049 A6a)."""
    kind, alt = str(site["alt_type"]), str(site["alt_aa"])
    ref_aa, ref_nt = str(site["ref_aa"]), str(site["ref_codon"])
    tol, table = cfg.variants.mass_equivalent_tol_da, cfg.annotation.translation_table
    found: set[str] = set()
    if kind == "substitution":
        if left_block or mass_equivalent(ref_aa, alt, tol):
            found.add("mass_equivalent")
        if alt not in ("U", "O"):  # U/O over its stop: zero changes (M7-16 rule 11)
            found.add(nt_class(ref_nt, alt, table=table))
        key = (str(site["contig"]), str(site["strand"]), int(site["codon_start"]))  # type: ignore[call-overload]
        if alt == "M" and key in starts:
            found.add("start_codon")
        if ref_aa == "*":
            found.add("stop_codon")
    elif kind == "swap":
        found.add("mass_equivalent")  # ADR-0049 A6b: every swap
    elif kind == "block" and mass_equivalent(translate_text(ref_nt, table), alt, tol):
        found.add("mass_equivalent")
    if cfg.genome_relation == "related":
        found.add("strain_difference")
    if ("database", "other") in set(alt_sources):
        found.add("database_derived")
    return [f for f in M705_FLAGS if f in found]


def site_flags(
    site: Mapping[str, object],
    gm: GeneModel,
    cfg: GippaConfig,
    alt_sources: Collection[tuple[str, str]],
) -> list[str]:
    """Rules 1-6 for one ``gippa.variants.pileup.sites`` row (0-based ``codon_start``).

    ``alt_sources`` holds the ``(source, search_db)`` of the peptides supporting the ALT. Returns
    each flag once in §9.4 order, never ``vcf_unrepresentable`` (M7-03). Frameshifts get rules
    3-4 only; blocks also rule 1 on the translation of ``ref_codon``; swaps are always
    ``mass_equivalent`` (ADR-0049). Pure."""
    return flags_for(site, start_codons(gm), cfg, alt_sources)
