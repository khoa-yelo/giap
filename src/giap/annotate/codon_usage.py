"""Codon usage per gene: CAI, tAI, ENC, Fop, GC3, codon-pair bias, ``at_2_10``,
``rare_codons_2_10`` (GIAP_SPEC.md §7.1.1, §7.1.6; docs/specs/M2-17-codon-usage.md).

Every index is computed on codons 2 .. n-1 (start and stop dropped) of eligible genes
(``scored_sequence``). CAI and Fop use the scored ribosomal-protein genes as reference; tAI the
tRNA gene copy numbers; codon-pair bias and the genome-wide synonymous shares all scored genes.
The indices come from codon-bias (MIT), always with ``genetic_code=11``.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import pandas as pd
from codonbias.scores import (
    CodonAdaptationIndex,
    CodonPairBias,
    EffectiveNumberOfCodons,
    FrequencyOfOptimalCodons,
    TrnaAdaptationIndex,
)

from giap.annotate.bakta_layers import is_ribosomal_protein
from giap.model.genes import Gene

CAI_REFERENCE_MISSING = "cai_reference_missing"
MIN_CODONS = 10
_WINDOW = slice(3, 30)  # codons 2-10 of the full CDS
_BASES = "TCAG"
_AA11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODE11: dict[str, str] = {
    a + b + c: _AA11[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}


@dataclass(frozen=True, slots=True)
class CodonFeatures:
    """Codon-usage features of one gene (``None`` = not defined, written as null)."""

    cai: float | None
    tai: float | None
    enc: float | None
    fop: float | None
    gc3: float | None
    cpb: float | None
    at_2_10: float | None
    rare_codons_2_10: int | None


def _code(table: int) -> dict[str, str]:
    if table not in (1, 11):  # tables 1 and 11 assign identical amino acids and stops
        raise ValueError(f"translation table {table} is not supported (only 11)")
    return _CODE11


def _codons(seq: str) -> list[str]:
    s = seq.upper()
    return [s[i : i + 3] for i in range(0, len(s) - len(s) % 3, 3)]


def scored_sequence(gene: Gene) -> str | None:
    """Codons 2 .. n-1; ``None`` for pseudogenes, partial genes, ``len % 3 != 0`` or < 10 codons."""
    nt = gene.nt
    if nt is None or gene.feature_type == "pseudogene" or gene.partial:
        return None
    if len(nt) % 3 or len(nt) < 3 * MIN_CODONS:
        return None
    return nt[3:-3]


def synonymous_shares(seqs: Iterable[str], table: int = 11) -> dict[str, float]:
    """Codon → share within its amino-acid family over the sense codons of ``seqs``.

    Codons with non-ACGT letters and stops are skipped; families never observed are absent,
    unseen codons of an observed family have share 0.0.
    """
    code = _code(table)
    n: Counter[str] = Counter()
    for s in seqs:
        n.update(c for c in _codons(s) if code.get(c, "*") != "*")
    families: dict[str, list[str]] = {}
    for c in sorted(code):
        if code[c] != "*":
            families.setdefault(code[c], []).append(c)
    out: dict[str, float] = {}
    for fam in families.values():
        total = sum(n[c] for c in fam)
        if total:
            out.update({c: n[c] / total for c in fam})
    return out


def rare_codons(shares: Mapping[str, float], cutoff: float = 0.10) -> frozenset[str]:
    """Codons whose synonymous share is strictly below ``cutoff``."""
    return frozenset(c for c, v in shares.items() if v < cutoff)


def at_2_10(nt: str) -> float | None:
    """A+T fraction of codons 2-10 of the full CDS over its A/C/G/T bases; ``None`` if < 30 nt."""
    if len(nt) < 3 * MIN_CODONS:
        return None
    w = nt[_WINDOW].upper()
    acgt = sum(w.count(b) for b in "ACGT")
    return None if acgt == 0 else (w.count("A") + w.count("T")) / acgt


def rare_2_10(nt: str, rare: frozenset[str]) -> int | None:
    """Number of rare codons among codons 2-10 of the full CDS; ``None`` if < 30 nt."""
    if len(nt) < 3 * MIN_CODONS:
        return None
    return sum(c in rare for c in _codons(nt[_WINDOW]))


def gc3(seq: str) -> float:
    """G+C fraction of the third codon positions (A/C/G/T only)."""
    thirds = [c[2] for c in _codons(seq) if c[2] in "ACGT"]
    return sum(b in "GC" for b in thirds) / len(thirds)


def _finite(x: object) -> float | None:
    v = float(x)  # type: ignore[arg-type]
    return v if math.isfinite(v) else None


def codon_usage(
    genes: Sequence[Gene], trna_counts: Mapping[str, int], table: int = 11
) -> tuple[dict[str, CodonFeatures], tuple[str, ...]]:
    """(locus_tag → features, run flags); see the module docstring for the references."""
    _code(table)
    scored = {g.locus_tag: scored_sequence(g) for g in genes}
    eligible = sorted(s for s in scored.values() if s is not None)
    ref = sorted(
        s for g in genes if (s := scored[g.locus_tag]) is not None and is_ribosomal_protein(g)
    )
    flags = () if ref else (CAI_REFERENCE_MISSING,)
    rare = rare_codons(synonymous_shares(eligible, table))
    cai = CodonAdaptationIndex(ref_seq=ref, genetic_code=table) if ref else None
    fop = FrequencyOfOptimalCodons(ref_seq=ref, genetic_code=table) if ref else None
    tai = None
    if trna_counts:
        tgcn = pd.DataFrame(
            {
                "anti_codon": sorted(trna_counts),
                "GCN": [trna_counts[a] for a in sorted(trna_counts)],
            }
        )
        tai = TrnaAdaptationIndex(tGCN=tgcn, prokaryote=True, genetic_code=table)
    enc = cpb = None
    if eligible:
        enc = EffectiveNumberOfCodons(genetic_code=table)
        cpb = CodonPairBias(ref_seq=eligible, genetic_code=table)

    def score(model: object | None, s: str) -> float | None:
        return None if model is None else _finite(model.get_score(s))  # type: ignore[attr-defined]

    out: dict[str, CodonFeatures] = {}
    for g in genes:
        s, nt = scored[g.locus_tag], g.nt or ""
        at, rc = at_2_10(nt), rare_2_10(nt, rare)
        if s is None:
            out[g.locus_tag] = CodonFeatures(None, None, None, None, None, None, at, rc)
            continue
        out[g.locus_tag] = CodonFeatures(
            score(cai, s), score(tai, s), score(enc, s), score(fop, s), gc3(s), score(cpb, s),
            at, rc,
        )  # fmt: skip
    return out, flags
