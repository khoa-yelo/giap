"""Strain differences (issue M3-18; spec §11.3): a related "reference" genome with SNVs and indels.

The simulated ``genome.fasta`` is the sample; the reference is the sample with the differences
applied, annotation lifted. Truth rows follow VCF orientation: ``ref`` is the reference allele at
``pos_reference``, ``alt`` the sample allele at ``pos_sample`` (``-`` = empty allele).
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from dataclasses import dataclass, replace
from typing import Any

import numpy as np

from giapsim.codons import CODE_11, STOPS, codons_of, revcomp, translate_cds
from giapsim.features import Feat, random_seq

_SPACING = 12  # minimum distance between two differences
_BOUNDARY_GAP = 10  # indels stay this far from feature ends (coordinates lift cleanly)
_COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


@dataclass(frozen=True)
class StrainConfig:
    snv_rate_per_kb: float = 1.0
    indel_rate_per_kb: float = 0.0
    seed_offset: int = 1
    indel_len: tuple[int, int] = (1, 6)


class _Lift:
    """Sample → reference coordinate map for one contig."""

    def __init__(self, shifts: list[tuple[int, int]]) -> None:
        shifts.sort()
        self.thr = [t for t, _ in shifts]
        self.cum = [0, *np.cumsum([d for _, d in shifts]).tolist()]

    def __call__(self, p: int) -> int:
        return p + int(self.cum[bisect_right(self.thr, p)])


def _gene_at(p: int, length: int, genes: list[Feat], starts: list[int], max_len: int) -> Any:
    """First non-pseudo CDS covering sample position ``p`` and ``p`` unwrapped into it."""
    lo = bisect_left(starts, p - max_len)
    for f in genes[lo : bisect_right(starts, p)]:
        if f.start <= p < f.end:
            return f, p
    for f in genes:
        if f.end > length and f.start <= p + length < f.end:
            return f, p + length
    return None, p


def _aa_change(f: Feat, q: int, ref_base: str) -> str:
    o = q - f.start if f.strand == "+" else f.end - 1 - q
    ci, k = divmod(o, 3)
    codon = f.nt[3 * ci : 3 * ci + 3]
    if len(codon) < 3:
        return ""
    b = ref_base if f.strand == "+" else _COMP[ref_base]
    return f"{CODE_11[codon[:k] + b + codon[k + 1 :]]}{ci + 1}{CODE_11[codon]}"


def derive_reference(
    cfg: StrainConfig, seed: int, seqs: dict[str, str], feats: list[Feat]
) -> tuple[dict[str, str], list[Feat], list[list[Any]]]:
    """Reference sequences, lifted features and ``strain_differences.tsv`` rows."""
    rng = np.random.default_rng([seed, 2, cfg.seed_offset])
    ref_seqs, lifts, rows = {}, {}, []
    for cid, s in seqs.items():
        length = len(s)
        on = [f for f in feats if f.contig == cid]
        bounds = sorted({p % length for f in on for p in (f.start, f.end)})
        genes = sorted(
            (f for f in on if f.kind == "cds" and not f.pseudo_reason), key=lambda f: f.start
        )
        starts = [f.start for f in genes]
        max_len = max((f.end - f.start for f in genes), default=0)
        n_snv = int(rng.poisson(cfg.snv_rate_per_kb * length / 1000))
        n_ind = int(rng.poisson(cfg.indel_rate_per_kb * length / 1000))
        n = min(n_snv + n_ind, max(0, (length - 40) // (2 * _SPACING)))
        pos = rng.choice(length - 40, size=n, replace=False) + 20 if n else np.zeros(0, int)
        is_indel = rng.permutation(np.arange(n) < n_ind)
        out, shifts, events, cursor, last = [], [], [], 0, -_SPACING
        for p, indel in sorted(zip(pos.tolist(), is_indel.tolist(), strict=True)):
            if p - last < _SPACING:
                continue
            if indel:
                k = int(rng.integers(cfg.indel_len[0], cfg.indel_len[1] + 1))
                i = bisect_left(bounds, p - _BOUNDARY_GAP)
                if i < len(bounds) and bounds[i] <= p + k + _BOUNDARY_GAP:
                    continue
                if rng.random() < 0.5:  # sample carries k extra bases
                    typ, ref, alt, shift = "ins", "", s[p : p + k], (p + k, -k)
                else:
                    typ, ref, alt, shift = "del", random_seq(rng, k), "", (p, k)
                shifts.append(shift)
            else:
                typ, alt = "snv", s[p]
                ref = "ACGT".replace(alt, "")[int(rng.integers(3))]
            out.append(s[cursor:p] + ref)
            cursor = p + len(alt)
            events.append((p, typ, ref, alt))
            last = p
        out.append(s[cursor:])
        ref_seqs[cid] = "".join(out)
        lifts[cid] = _Lift(shifts)
        for p, typ, ref, alt in events:
            f, q = _gene_at(p, length, genes, starts, max_len)
            change = _aa_change(f, q, ref) if f is not None and typ == "snv" else ""
            pr = lifts[cid](p) - len(ref) if typ == "del" else lifts[cid](p)
            rows.append(
                [cid, p + 1, pr + 1, ref or "-", alt or "-", typ, f.locus if f else "", change]
            )
    lifted = [_lift(f, lifts[f.contig], len(seqs[f.contig]), ref_seqs[f.contig]) for f in feats]
    return ref_seqs, lifted, rows


def _lift(f: Feat, lift: _Lift, length: int, ref: str) -> Feat:
    rs = lift(f.start)
    re = lift((f.end - 1) % length) + 1 + (len(ref) if f.end > length else 0)
    g = replace(f, start=rs, end=re, true_start=None, extra=dict(f.extra))
    if f.kind not in ("cds", "sorf"):
        return g
    seg = ref[rs:re] if re <= len(ref) else ref[rs:] + ref[: re - len(ref)]
    g.nt = seg if f.strand == "+" else revcomp(seg)
    cod = codons_of(g.nt)
    clean = len(g.nt) % 3 == 0 and cod[-1] in STOPS and not STOPS.intersection(cod[1:-1])
    g.extra["start_type"] = g.nt[:3]
    if clean and not f.pseudo_reason:
        g.aa = translate_cds(cod)
    else:
        g.extra["pseudo"] = True
    return g
