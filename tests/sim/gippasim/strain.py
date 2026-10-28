"""Strain differences (issue M3-18; spec §11.3): a related "reference" genome with SNVs and indels.

The simulated ``genome.fasta`` is the sample; the reference is the sample with the differences
applied, annotation lifted. Truth rows follow VCF orientation: ``ref`` is the reference allele at
``pos_reference``, ``alt`` the sample allele at ``pos_sample``. Indels are left-normalized and
anchored as in VCF: ``pos_*`` is the 1-based position of the base before the event and both
alleles include that anchor base. ``aa_change`` reads ``<reference residue><codon><sample
residue>`` for a non-synonymous SNV and ``<reference residues><position><sample residues>``
(``-`` = none) for an in-frame indel, from the translated proteins (codon 1 reads as M).
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from dataclasses import dataclass, replace
from typing import Any

import numpy as np

from gippasim.codons import code, codons_of, revcomp, stops
from gippasim.features import STARTS, Feat, random_seq
from gippasim.proteomics import digest_spans

_SPACING = 12  # minimum distance between two differences
_BOUNDARY_GAP = 10  # indels stay this far from feature ends (coordinates lift cleanly)
_COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


@dataclass(frozen=True)
class StrainConfig:
    snv_rate_per_kb: float = 1.0
    indel_rate_per_kb: float = 0.0
    seed_offset: int = 1
    indel_len: tuple[int, int] = (1, 6)
    n_stop_loss: int = 0  # M7-16 rule 20: SNVs turning sample terminal stops into sense codons


Digest = tuple[tuple[int, int], int]  # (digest_len, missed_cleavages) of the proteomics layer
STOP_LOSS_COLUMNS = [
    "locus_tag", "contig", "strand", "codon_start", "sample_codon", "reference_codon",
    "reference_residue", "cterm_peptide", "planted",
]  # fmt: skip


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


def _offset(f: Feat, q: int) -> int:
    """Offset of (unwrapped) position ``q`` in ``f`` in transcription direction."""
    return q - f.start if f.strand == "+" else f.end - 1 - q


def _residue(codon: str, index: int, table: int = 11) -> str:
    """Residue of codon ``index`` (0-based) under ``table``; a valid start codon 1 reads as M."""
    return "M" if index == 0 and codon in STARTS else code(table)[codon]


def protein(nt: str, table: int = 11) -> str:
    """Translation of ``nt`` up to its first stop (codon 1 as M when it is a valid start)."""
    out, st = [], stops(table)
    for i, c in enumerate(codons_of(nt)):
        if c in st:
            break
        out.append(_residue(c, i, table))
    return "".join(out)


def is_clean(nt: str, table: int = 11) -> bool:
    """Valid start, length a multiple of 3, a final stop and no internal stop."""
    cod, st = codons_of(nt), stops(table)
    return (
        len(nt) % 3 == 0
        and len(cod) >= 2
        and cod[0] in STARTS
        and cod[-1] in st
        and not st.intersection(cod[1:-1])
    )


def _snv_change(f: Feat, q: int, g: Feat, qr: int, table: int = 11) -> str:
    """``<reference residue><codon number><sample residue>`` for a SNV, "" when synonymous."""
    cs, cr = divmod(_offset(f, q), 3)[0], divmod(_offset(g, qr), 3)[0]
    a, b = f.nt[3 * cs : 3 * cs + 3], g.nt[3 * cr : 3 * cr + 3]
    if len(a) < 3 or len(b) < 3:
        return ""
    ra, rb = _residue(b, cr, table), _residue(a, cs, table)
    return "" if ra == rb else f"{ra}{cs + 1}{rb}"


def _indel_change(f: Feat, g: Feat, table: int = 11) -> str:
    """Protein difference ``<reference residues><first position><sample residues>`` (``-`` =
    none) for an in-frame indel, from the two translated proteins."""
    a, b = protein(g.nt, table), protein(f.nt, table)
    i = 0
    while i < min(len(a), len(b)) and a[i] == b[i]:
        i += 1
    j = 0
    while j < min(len(a), len(b)) - i and a[-1 - j] == b[-1 - j]:
        j += 1
    ra, rb = a[i : len(a) - j], b[i : len(b) - j]
    return f"{ra or '-'}{i + 1}{rb or '-'}"


def _normalize(s: str, p: int, seq: str) -> tuple[int, str]:
    """Shift an event inserting ``seq`` before sample position ``p`` as far left as possible."""
    while p > 2 and s[p - 1] == seq[-1]:
        p, seq = p - 1, s[p - 1] + seq[:-1]
    return p, seq


def derive_reference(
    cfg: StrainConfig, seed: int, seqs: dict[str, str], feats: list[Feat], table: int = 11
) -> tuple[dict[str, str], list[Feat], list[list[Any]]]:
    """Reference sequences, lifted features and ``strain_differences.tsv`` rows.

    Indels are left-normalized and written VCF-style: ``pos`` is the base before the event and
    both alleles start with it.
    """
    ref_seqs, lifted, rows, _ = derive_strain(cfg, seed, seqs, feats, table=table)
    return ref_seqs, lifted, rows


def _terminal(f: Feat, n: int) -> int:
    """Lowest 0-based coordinate (wrapped) of a feature's last codon."""
    return (f.end - 3) % n if f.strand == "+" else f.start % n


def _cterm_peptide(aa: str, dg: Digest) -> str:
    """The shortest C-terminal Trypsin/P product of ``aa`` with >= 10 aa ("" if none)."""
    (lo, hi), mc = dg
    ends = [a for a, b in digest_spans(aa, (max(10, lo), hi), mc) if b == len(aa)]
    return aa[max(ends) :] if ends else ""


def _near(q: int, events: list[tuple[int, str, str, str]], n: int, gap: int) -> bool:
    """Some event position (an indel's anchor and end included) lies within ``gap`` of ``q``."""
    for p, typ, _ref, alt in events:
        for x in (p,) if typ == "snv" else (p - 1, p, p + len(alt)):
            d = abs(q - x) % n
            if min(d, n - d) < gap:
                return True
    return False


def _in_span(events: list[tuple[int, str, str, str]], lo: int, hi: int, n: int) -> bool:
    return any((p - lo) % n < hi - lo for p, *_ in events)


def _plant_stop_losses(
    cfg: StrainConfig, seed: int, seqs: dict[str, str], feats: list[Feat],
    ref_seqs: dict[str, str], lifts: dict[str, _Lift], all_events: dict[str, list[Any]],
    dg: Digest, table: int = 11,
) -> dict[str, str]:  # fmt: skip
    """Rule 20's knob: ``n_stop_loss`` SNVs (own generator) in the terminal stop codon of distinct
    non-pseudo CDS whose reference C-terminal product (>= 10 aa) holds no other event. Edits
    ``ref_seqs`` and ``all_events`` in place; returns locus -> ``cterm_peptide``.

    Genes whose reference protein then differs from the sample's only by the read-through residue
    are tried first, so the wrong peptide is a pure read-through where the build allows it.
    """
    rng = np.random.default_rng([seed, 2, cfg.seed_offset, 1])
    st, tr = stops(table), code(table)
    cds = [f for f in feats if f.kind == "cds" and not f.pseudo_reason]
    on = {c: sorted((f for f in cds if f.contig == c), key=lambda f: f.start) for c in seqs}
    cands: list[tuple[bool, Feat, int, str, str, int]] = []
    for k in rng.permutation(len(cds)).tolist():
        f = cds[k]
        n, ev = len(seqs[f.contig]), all_events[f.contig]
        a = _terminal(f, n)
        g = _lift(f, lifts[f.contig], n, ref_seqs[f.contig], table)
        stop = g.nt[-3:]
        opts = [(i, b) for i in range(3) for b in "ACGT" if b != stop[i]
                and stop[:i] + b + stop[i + 1:] not in st]  # fmt: skip
        i, b = opts[int(rng.integers(len(opts)))]
        q = (a + i if f.strand == "+" else a + 2 - i) % n
        sense = stop[:i] + b + stop[i + 1 :]
        aa = protein(g.nt[:-3] + sense, table)
        pep = _cterm_peptide(aa, dg) if len(aa) == len(g.nt) // 3 else ""
        if not pep or stop not in st or _near(q, ev, n, _SPACING):
            continue
        span = 3 * len(pep)
        lo = a + 3 - span if f.strand == "+" else a
        if _in_span(ev, lo, lo + span, n):
            continue
        genes = on[f.contig]
        hit, _ = _gene_at(
            q, n, genes, [h.start for h in genes], max(h.end - h.start for h in genes)
        )
        if hit is not f:  # the strain row must name this gene
            continue
        base = b if f.strand == "+" else _COMP[b]
        cands.append((aa != f.aa + tr[sense], f, q, base, pep, lo))
    cands.sort(key=lambda c: c[0])  # stable: pure read-through genes first, then the draw order
    out: dict[str, str] = {}
    spans: list[tuple[str, int, int]] = []
    for _impure, f, q, base, pep, lo in cands:
        if len(out) == cfg.n_stop_loss:
            break
        ev, n, hi = all_events[f.contig], len(seqs[f.contig]), lo + 3 * len(pep)
        if (
            _near(q, ev, n, _SPACING)
            or _in_span(ev, lo, hi, n)
            or any(c == f.contig and (q - a) % n < b - a for c, a, b in spans)
        ):  # spacing from, and peptides clean of, the SNVs planted so far
            continue
        spans.append((f.contig, lo, hi))
        r = lifts[f.contig](q)
        ref_seqs[f.contig] = ref_seqs[f.contig][:r] + base + ref_seqs[f.contig][r + 1 :]
        ev.append((q, "snv", base, seqs[f.contig][q]))
        ev.sort()
        out[f.locus] = pep
    if len(out) < cfg.n_stop_loss:
        raise ValueError(f"n_stop_loss {cfg.n_stop_loss}: only {len(out)} eligible genes")
    return out


def stop_loss_rows(
    seqs: dict[str, str], feats: list[Feat], lifted: list[Feat], planted: dict[str, str],
    dg: Digest, table: int = 11,
) -> list[list[Any]]:  # fmt: skip
    """``truth/stop_loss.tsv``: every non-pseudo sample CDS whose terminal stop is a sense codon
    in the reference (planted or from the uniform draws)."""
    by_locus = {g.locus: g for g in lifted}
    rows, st = [], stops(table)
    for f in feats:
        if f.kind != "cds" or f.pseudo_reason:
            continue
        g, n = by_locus[f.locus], len(seqs[f.contig])
        sc, rc = f.nt[-3:], g.nt[-3:]
        if sc not in st or rc in st or len(rc) < 3:
            continue
        pep = planted.get(f.locus) or (
            _cterm_peptide(g.aa, dg) if len(g.aa) == len(g.nt) // 3 else ""
        )
        rows.append([f.locus, f.contig, f.strand, _terminal(f, n), sc, rc, code(table)[rc], pep,
                     str(f.locus in planted).lower()])  # fmt: skip
    return rows


def derive_strain(
    cfg: StrainConfig, seed: int, seqs: dict[str, str], feats: list[Feat],
    dg: Digest = ((7, 30), 1), table: int = 11,
) -> tuple[dict[str, str], list[Feat], list[list[Any]], dict[str, str]]:  # fmt: skip
    """:func:`derive_reference` plus the ``n_stop_loss`` planting (locus -> ``cterm_peptide``;
    empty when ``n_stop_loss`` is 0, which leaves every draw and output unchanged)."""
    rng = np.random.default_rng([seed, 2, cfg.seed_offset])
    ref_seqs, lifts, all_events = {}, {}, {}
    for cid, s in seqs.items():
        length = len(s)
        on = [f for f in feats if f.contig == cid]
        bounds = sorted({p % length for f in on for p in (f.start, f.end)})
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
                if rng.random() < 0.5:  # sample carries k extra bases s[p:p+k]
                    typ, ref = "ins", ""
                    q, alt = _normalize(s, p, s[p : p + k])
                else:  # the reference carries k extra bases before sample position p
                    typ, alt = "del", ""
                    q, ref = _normalize(s, p, random_seq(rng, k))
                if q - 1 - last < _SPACING:
                    continue
                i = bisect_left(bounds, q - 1 - _BOUNDARY_GAP)
                if i < len(bounds) and bounds[i] <= q + k + _BOUNDARY_GAP:
                    continue
                p = q
                shifts.append((p + k, -k) if typ == "ins" else (p, k))
            else:
                typ, alt = "snv", s[p]
                ref = "ACGT".replace(alt, "")[int(rng.integers(3))]
            out.append(s[cursor:p] + ref)
            cursor = p + len(alt)
            events.append((p, typ, ref, alt))
            last = p + len(alt)
        out.append(s[cursor:])
        ref_seqs[cid] = "".join(out)
        lifts[cid] = _Lift(shifts)
        all_events[cid] = events
    planted = (
        _plant_stop_losses(cfg, seed, seqs, feats, ref_seqs, lifts, all_events, dg, table)
        if cfg.n_stop_loss
        else {}
    )
    lifted = [_lift(f, lifts[f.contig], len(seqs[f.contig]), ref_seqs[f.contig], table)
              for f in feats]  # fmt: skip
    by_locus = {g.locus: g for g in lifted}
    rows = []
    for cid, events in all_events.items():
        s, length, lift, rlen = seqs[cid], len(seqs[cid]), lifts[cid], len(ref_seqs[cid])
        genes = sorted(
            (f for f in feats if f.contig == cid and f.kind == "cds" and not f.pseudo_reason),
            key=lambda f: f.start,
        )
        starts = [f.start for f in genes]
        max_len = max((f.end - f.start for f in genes), default=0)
        for p, typ, ref, alt in events:
            a = p if typ == "snv" else p - 1  # 0-based sample position of the row (anchor)
            if typ != "snv":
                ref, alt = s[a] + ref, s[a] + alt
            f, q = _gene_at(a, length, genes, starts, max_len)
            pr = lift(a)
            change = ""
            if f is not None:
                g = by_locus[f.locus]
                qr = pr + (rlen if q != a else 0)
                if typ == "snv":
                    change = _snv_change(f, q, g, qr, table)
                elif (len(ref) - len(alt)) % 3 == 0:
                    change = _indel_change(f, g, table)
            rows.append([cid, a + 1, pr + 1, ref, alt, typ, f.locus if f else "", change])
    return ref_seqs, lifted, rows, planted


def _lift(f: Feat, lift: _Lift, length: int, ref: str, table: int = 11) -> Feat:
    rs = lift(f.start)
    re = lift((f.end - 1) % length) + 1 + (len(ref) if f.end > length else 0)
    g = replace(f, start=rs, end=re, true_start=None, extra=dict(f.extra))
    if f.kind not in ("cds", "sorf"):
        return g
    seg = ref[rs:re] if re <= len(ref) else ref[rs:] + ref[: re - len(ref)]
    g.nt = seg if f.strand == "+" else revcomp(seg)
    g.extra["start_type"] = g.nt[:3]
    g.aa = protein(g.nt, table)  # never a protein the reference ``nt`` does not encode
    if f.pseudo_reason or not is_clean(g.nt, table):
        g.extra["pseudo"] = True
    return g
