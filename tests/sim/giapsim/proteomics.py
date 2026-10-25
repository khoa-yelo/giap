"""Simulated proteomics: proteins, tryptic peptides and database-search errors (issue M3-19).

GIAP_SPEC.md §11.3 "Expression and proteomics". Expressed proteins are the *true* proteins of the
sample genome (non-pseudo CDS and sORFs from their true start, plus the unannotated novel ORFs);
each distinct Trypsin/P product is detected once by a logistic draw in its length and the
protein's log abundance and then identified in every run. Wrong identifications (a share
``pi_wrong`` of the distinct target peptides) come from a residue-shuffled proteome or from the
related strain's proteome and never occur, up to I/L, in a protein of the sample. Scores follow a
two-component normal mixture from which PEP and q-values are computed consistently; decoy rows
are reversed tryptic peptides with null scores.

``wrong_source = "sixframe_noncoding"`` (M7-12 decision 12, opt-in) draws the wrong peptides from
the Trypsin/P products of the six-frame translation of the genome that are in frame with no
annotated CDS/sORF, no N-terminal extension and no novel ORF and occur in no true protein: the
false novel candidates a six-frame search produces. They are drawn in turn from the intergenic,
out-of-frame (same-strand CDS overlap) and antisense (opposite-strand overlap only) territories,
then (decision 22) the products in frame upstream of a CDS/sORF's true start with no stop in
between (``nterm_extension``), so every novel class gets false candidates; their truth rows carry
the source locus, and each lies in an entry of ``giap makedb --six-frame-orfs``. These draws use
their own stream (``[seed, 7]``). Decision 21 (stream ``[seed, 8]``): each null PSM row of a wrong
peptide becomes, with probability 1/2, the decoy row of its makedb entry (reversed peptide,
``rev_GIAPORF_…``, the same null score and mixture q / PEP); gene-level decoys are reversed
products of annotated non-pseudo proteins (``rev_<locus_tag>``, ``n_decoy_rows`` of them).
M7-11b decision 9: in this mode, and in ``decoys = "scored_null"`` mode whatever
``wrong_source`` (MVP-04 step 5), the search database is the makedb FASTA, so target rows name
makedb entries too: a wrong peptide's row its entry, a novel ORF's peptide the entry it lies in
(``GIAPORF_…``), a peptide of a true protein upstream of its annotated start the gene's
``nterm_extension`` entry and any other peptide of an annotated protein its locus tag (makedb's
``db_id``). Truth rows keep the planted IDs. A target occurrence in no makedb entry (e.g. upstream
of a gene whose upstream part is < 7 aa, which makedb drops) is not in the search database, so a
peptide with no other occurrence is not detected.

``decoys = "scored_null"`` (MVP-04 steps 5-6, opt-in) replaces the wrong peptides and decoys: the
null peptides split 50/50 between protein entries and six-frame entries (lead decision 10: half,
rounded down, go to protein entries, capped by the pool of protein nulls; the six-frame side takes
the rest). ``null_split`` in ``truth/proteomics.json`` records ``mode = "fifty_fifty"``, the
realised counts and, for information only, the protein share of the distinct Trypsin/P products of
the simulated search space (the annotated proteins against the six-frame products of the non-coding
pool not also produced by a protein). A protein-entry null (lead decision 9) is a real Trypsin/P
product of a non-pseudo recorded CDS / sORF none of whose products was detected (including those the
sample strain made pseudo, translated on the sample, no product across a stop), the locus drawn
uniformly (without replacement) over those with an undetected product and the product uniformly
among its undetected ones (stream ``[seed, 9]``); its truth rows (``origin = undetected_product``)
are placed like correct rows and its target rows carry the locus tag. Six-frame nulls are drawn as
in ``sixframe_noncoding`` and carry their makedb entry ID. Each null PSM row then goes, on its own
1/2 coin (stream ``[seed, 8]``), to the target or to the decoy of its entry (reversed peptide,
``rev_<entry>``, the same score, q and PEP, intensity NA); a null whose reversal is a target peptide
(I = L) stays a target and the PSMs the coin would have sent over are counted
(``decoy_collisions``). No other decoys are written.

All draws come from ``numpy.random.default_rng([seed, 3])``, so the genome and its truth are
unchanged by the proteomics options. Gene-level abundance is a plain log-normal per protein
(``gene_log_abundance``; M4 replaces it with the TU/region/AR(1) model).
"""

from __future__ import annotations

import json
import logging
import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
from scipy import stats

from giapsim.codons import CODE_11, STOPS, codons_of, revcomp, translate_cds
from giapsim.features import Feat

_LOG = logging.getLogger(__name__)

__all__ = [
    "Protein",
    "ProteinSimConfig",
    "ProteomicsConfig",
    "digest",
    "digest_spans",
    "expressed_proteins",
    "gene_log_abundance",
    "mixture_pep",
    "mixture_q",
    "simulate_proteomics",
]

AA20 = "ACDEFGHIKLMNPQRSTVWY"
_NA = "NA"
#: Mean and SD of the per-protein natural-log abundance (placeholder until M4).
LOG_ABUNDANCE = (1.0, 1.5)
INTENSITY_SD = 0.5
P_CHARGE_2 = 0.6
_STREAM = 3  # numpy SeedSequence key of the proteomics stream (1: features, 2: signals)
_SIXFRAME_STREAM = 7  # numpy SeedSequence key of the six-frame wrong peptides (M7-12)
WrongSource = Literal["shuffled", "other_strain", "sixframe_noncoding"]
_WRONG_SOURCES = ("shuffled", "other_strain", "sixframe_noncoding")
_DECOY_STREAM = 8  # numpy SeedSequence key of the six-frame decoy coin and protein decoys
SIXFRAME_CLASSES = ("intergenic", "cds_outframe", "antisense", "nterm_extension")
DecoyMode = Literal["unscored", "scored_null"]
_DECOY_MODES = ("unscored", "scored_null")
_NULL_STREAM = 9  # numpy SeedSequence key of the scored-null protein nulls (MVP-04)
#: ``truth/peptides.tsv`` origin of a scored-null protein-entry null (MVP-04 lead decision 9).
PROTEIN_NULL_ORIGIN = "undetected_product"


@dataclass(frozen=True)
class ProteinSimConfig:
    """Cross-strain protein copies in ``proteins.faa`` (for the miniprot/tblastn routes)."""

    identity: tuple[float, float] = (0.80, 0.99)
    fraction: float = 0.2  # share of annotated CDS proteins (>= 50 aa) that get a copy
    min_length_aa: int = 50

    def __post_init__(self) -> None:
        object.__setattr__(self, "identity", tuple(float(x) for x in self.identity))
        lo, hi = self.identity
        if not 0.0 < lo <= hi < 1.0:
            raise ValueError("identity must satisfy 0 < lo <= hi < 1")
        if not 0.0 <= self.fraction <= 1.0:
            raise ValueError("fraction must be in [0, 1]")


@dataclass(frozen=True)
class ProteomicsConfig:
    """Peptide simulation knobs (spec §11.3; docs/specs/M3-19-sim-peptides-db-errors.md)."""

    digest_len: tuple[int, int] = (7, 30)
    missed_cleavages: int = 1  # trypsin: after K/R, also before P (Trypsin/P)
    detect_logit: tuple[float, float, float] = (-2.0, 0.08, 0.9)  # intercept, per aa, log abund.
    pi_wrong: float = 0.01
    wrong_source: WrongSource = "shuffled"
    score_mixture: tuple[float, float, float, float] = (3.0, 1.0, 0.0, 1.0)  # mu1, sd1, mu0, sd0
    n_decoy_rows: Literal["match_wrong"] | int = "match_wrong"
    search_db: Literal["same", "other"] = "same"
    n_runs: int = 3
    proteins: ProteinSimConfig | None = None
    decoys: DecoyMode = "unscored"

    def __post_init__(self) -> None:
        for name in ("digest_len", "detect_logit", "score_mixture"):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        lo, hi = self.digest_len
        if not 1 <= lo <= hi:
            raise ValueError("digest_len must satisfy 1 <= lo <= hi")
        if self.missed_cleavages < 0:
            raise ValueError("missed_cleavages must be >= 0")
        if not 0.0 <= self.pi_wrong < 1.0:
            raise ValueError("pi_wrong must be in [0, 1)")
        if self.wrong_source not in _WRONG_SOURCES:
            raise ValueError(f"wrong_source must be one of {_WRONG_SOURCES}")
        if self.search_db not in ("same", "other"):
            raise ValueError("search_db must be 'same' or 'other'")
        if self.score_mixture[1] <= 0 or self.score_mixture[3] <= 0:
            raise ValueError("score_mixture SDs must be > 0")
        if self.n_decoy_rows != "match_wrong" and (
            not isinstance(self.n_decoy_rows, int) or self.n_decoy_rows < 0
        ):
            raise ValueError("n_decoy_rows must be 'match_wrong' or an int >= 0")
        if self.n_runs < 1:
            raise ValueError("n_runs must be >= 1")
        if self.decoys not in _DECOY_MODES:
            raise ValueError(f"decoys must be one of {_DECOY_MODES}, got {self.decoys!r}")


# ---------------------------------------------------------------- digest


def digest_spans(
    protein: str, digest_len: tuple[int, int] = (7, 30), missed_cleavages: int = 1
) -> list[tuple[int, int]]:
    """Spans ``[a, b)`` of the Trypsin/P products of ``protein`` (cut after every K/R), in order
    of ``a`` then ``b``, with ``digest_len[0] <= b - a <= digest_len[1]``."""
    if not protein:
        return []
    cuts = [0, *(i + 1 for i, c in enumerate(protein[:-1]) if c in "KR"), len(protein)]
    lo, hi = digest_len
    out = []
    for i in range(len(cuts) - 1):
        for j in range(i + 1, min(i + 2 + missed_cleavages, len(cuts))):
            if lo <= cuts[j] - cuts[i] <= hi:
                out.append((cuts[i], cuts[j]))
    return out


def digest(
    protein: str, digest_len: tuple[int, int] = (7, 30), missed_cleavages: int = 1
) -> list[str]:
    """Distinct Trypsin/P products of ``protein`` in order of first occurrence."""
    spans = digest_spans(protein, digest_len, missed_cleavages)
    return list(dict.fromkeys(protein[a:b] for a, b in spans))


# ---------------------------------------------------------------- statistics


def mixture_pep(
    scores: np.ndarray, pi: float, mix: tuple[float, float, float, float]
) -> np.ndarray:
    """Posterior error probability pi f0 / (pi f0 + (1 - pi) f1) of each score."""
    mu1, sd1, mu0, sd0 = mix
    s = np.asarray(scores, dtype=float)
    if pi <= 0.0:
        return np.zeros_like(s)
    l0 = math.log(pi) + stats.norm.logpdf(s, mu0, sd0)
    l1 = math.log1p(-pi) + stats.norm.logpdf(s, mu1, sd1)
    out: np.ndarray = np.exp(l0 - np.logaddexp(l0, l1))
    return out


def mixture_q(scores: np.ndarray, pi: float, mix: tuple[float, float, float, float]) -> np.ndarray:
    """q-values: the mixture FDR pi S0 / (pi S0 + (1 - pi) S1) of the score threshold, made
    monotone (minimum over every threshold at or below the score)."""
    mu1, sd1, mu0, sd0 = mix
    s = np.asarray(scores, dtype=float)
    if pi <= 0.0 or s.size == 0:
        return np.zeros_like(s)
    l0 = math.log(pi) + stats.norm.logsf(s, mu0, sd0)
    l1 = math.log1p(-pi) + stats.norm.logsf(s, mu1, sd1)
    fdr: np.ndarray = np.exp(l0 - np.logaddexp(l0, l1))
    order = np.argsort(s, kind="stable")
    q: np.ndarray = np.empty_like(fdr)
    q[order] = np.minimum.accumulate(fdr[order])
    return q


# ---------------------------------------------------------------- proteins


@dataclass(frozen=True)
class Protein:
    """An expressed true protein; 0-based ``[start, end)`` unwrapped (end may exceed the contig)."""

    locus: str
    contig: str
    start: int
    end: int
    strand: str
    aa: str


def _slice(seq: str, s0: int, e0: int, strand: str) -> str:
    n = len(seq)
    nt = seq[s0:e0] if e0 <= n else seq[s0:] + seq[: e0 - n]
    return revcomp(nt) if strand == "-" else nt


def expressed_proteins(
    feats: Sequence[Feat],
    novel: Sequence[Sequence[Any]],
    seqs: dict[str, str],
    readthrough: Mapping[str, tuple[str, str]] | None = None,
) -> list[Protein]:
    """Non-pseudo CDS and sORFs translated from the true start, then the novel ORFs.

    ``readthrough`` (M7-16 rule 14) maps a locus tag to (residue read at the stop, extension
    residues): that protein continues through its annotated stop to the next in-frame stop, and
    its interval grows by the stop codon and the extension codons."""
    out = []
    rt = readthrough or {}
    for f in feats:
        if f.kind not in ("cds", "sorf") or f.pseudo_reason or f.extra.get("pseudo"):
            continue
        s0, e0 = f.span()
        aa = translate_cds(codons_of(_slice(seqs[f.contig], s0, e0, f.strand)))
        if f.sec_residues:  # planted Sec (M7-16c): the TGA is decoded as U, after any extension
            up = ((f.start - s0) if f.strand == "+" else (e0 - f.end)) // 3
            for i in f.sec_residues:
                aa = aa[: up + i] + "U" + aa[up + i + 1 :]
        if f.locus in rt:
            residue, ext = rt[f.locus]
            aa += residue + ext
            grow = 3 * len(ext)
            if f.strand == "+":
                e0 += grow
            else:
                s0 -= grow
                if s0 < 0:  # the extension crosses the origin of a circular replicon
                    s0, e0 = s0 + len(seqs[f.contig]), e0 + len(seqs[f.contig])
        out.append(Protein(f.locus, f.contig, s0, e0, f.strand, aa))
    for orf_id, _cls, contig, s1, e0, strand, aa in novel:
        out.append(Protein(str(orf_id), str(contig), int(s1) - 1, int(e0), str(strand), str(aa)))
    return out


def gene_log_abundance(rng: np.random.Generator, n: int) -> np.ndarray:
    """Natural-log abundance of ``n`` proteins: i.i.d. normal (M4 replaces this model)."""
    out: np.ndarray = np.round(rng.normal(*LOG_ABUNDANCE, size=n), 6)
    return out


def _il(s: str) -> str:
    return s.replace("I", "L")


# ---------------------------------------------------------------- six-frame locus counts

_NT_CODE = np.full(256, 0, dtype=np.int64)
for _i, _b in enumerate(b"ACGT"):
    _NT_CODE[_b] = _i
_CODON_AA = np.frombuffer(
    "".join(CODE_11[a + b + c] for a in "ACGT" for b in "ACGT" for c in "ACGT").encode(),
    dtype=np.uint8,
)
_K = 7  # prefix length of the six-frame index


def _translate_np(nt: bytes) -> bytes:
    m = len(nt) // 3
    c = _NT_CODE[np.frombuffer(nt[: 3 * m], dtype=np.uint8)].reshape(m, 3)
    out: bytes = _CODON_AA[16 * c[:, 0] + 4 * c[:, 1] + c[:, 2]].tobytes()
    return out


def _kmer_codes(aa: bytes, k: int) -> np.ndarray:
    a = np.frombuffer(aa, dtype=np.uint8).astype(np.int64) & 31
    m = len(a) - k + 1
    if m <= 0:
        return np.zeros(0, dtype=np.int64)
    code = np.zeros(m, dtype=np.int64)
    for j in range(k):
        code = (code << 5) | a[j : j + m]
    return code


def sixframe_counts(
    peptides: Sequence[str], seqs: dict[str, str], topology: dict[str, str]
) -> dict[str, int]:
    """Number of distinct loci (contig, plus-strand start, strand) whose six-frame translation is
    exactly each peptide with I = L folded (GIAP's exact matcher); circular contigs are searched
    across the origin. Peptides without a hit count 0."""
    peps = sorted({_il(p) for p in peptides if len(p) >= _K})
    by_prefix: dict[int, list[str]] = {}
    for p in peps:
        by_prefix.setdefault(int(_kmer_codes(p[:_K].encode(), _K)[0]), []).append(p)
    keys = np.array(sorted(by_prefix), dtype=np.int64)
    loci: dict[str, set[tuple[str, int, str]]] = {p: set() for p in peps}
    max_len = max((len(p) for p in peps), default=0)
    for cid, seq in seqs.items():
        n = len(seq)
        circular = topology.get(cid, "circular") == "circular"
        ext = seq + seq[: 3 * max_len] if circular else seq
        for strand, src in (("+", ext), ("-", revcomp(ext))):
            raw = src.encode()
            for off in range(3):
                aa = _translate_np(raw[off:]).replace(b"I", b"L")
                hits = np.flatnonzero(np.isin(_kmer_codes(aa, _K), keys))
                text = aa.decode()
                for i in hits.tolist():
                    for p in by_prefix[(int(_kmer_codes(aa[i : i + _K], _K)[0]))]:
                        if text.startswith(p, i):
                            s = off + 3 * i if strand == "+" else len(ext) - off - 3 * (i + len(p))
                            loci[p].add((cid, s % n if circular else s, strand))
    return {p: len(loci.get(_il(p), ())) for p in peptides}


# ---------------------------------------------------------------- six-frame non-coding peptides

#: A six-frame product: (peptide, contig, 1-based start, 1-based end (closed), strand).
Locus = tuple[str, str, int, int, str]


def _codon_mask(n: int, spans: Sequence[tuple[int, int, str]]) -> dict[str, np.ndarray]:
    """Per strand, ``True`` at the strand-local first base of every codon in frame with a span
    ``[s0, e0)`` (``e0`` may exceed ``n``: crossing the origin) and overlapping it."""
    mask = {"+": np.zeros(n, dtype=bool), "-": np.zeros(n, dtype=bool)}
    for s0, e0, strand in spans:
        if strand == "+":
            mask["+"][np.arange(s0, e0, 3) % n] = True
        else:
            c = np.arange(e0 - 3, s0 - 3, -3) % n  # plus-strand low base of each codon
            loc = n - c - 3
            mask["-"][loc[loc >= 0]] = True
    return mask


def _cover(n: int, spans: Sequence[tuple[int, int, str]]) -> dict[str, np.ndarray]:
    """Per strand, the cumulative count of plus-strand bases covered by a span (length n + 1)."""
    out = {}
    for st in "+-":
        d = np.zeros(n + 1, dtype=np.int64)
        for s0, e0, strand in spans:
            if strand != st:
                continue
            for a, b in ((s0, min(e0, n)), (0, e0 - n)) if e0 > n else ((s0, e0),):
                if b > a:
                    d[a] += 1
                    d[b] -= 1
        out[st] = np.concatenate([[0], np.cumsum((np.cumsum(d[:n]) > 0).astype(np.int64))])
    return out


@dataclass
class SixframePool:
    """Trypsin/P products of six frames by class: ``idx[cls]`` rows (frame, aa start, aa end)."""

    frames: list[tuple[str, str, int, bytes, int]]  # contig, strand, offset, translation, length
    idx: dict[str, np.ndarray]

    def locus(self, cls: str, i: int) -> Locus:
        f, a0, b0 = (int(x) for x in self.idx[cls][i])
        cid, strand, off, aa, n = self.frames[f]
        t0, t1 = off + 3 * a0, off + 3 * b0  # strand-local nt
        s0, e0 = (t0, t1) if strand == "+" else (n - t1, n - t0)
        return aa[a0:b0].decode(), cid, s0 + 1, e0, strand


def sixframe_noncoding_pool(
    seqs: dict[str, str],
    coding: Sequence[tuple[str, int, int, str]],
    annotated: Sequence[tuple[str, int, int, str]],
    digest_len: tuple[int, int],
    missed_cleavages: int,
    true_starts: Sequence[tuple[str, int, str]] = (),
) -> SixframePool:
    """Trypsin/P products (linear, no origin crossing) of the six-frame translation of ``seqs``
    with no codon in frame with a ``coding`` span (contig, s0, e0, strand), by class of their
    overlap with the ``annotated`` CDS/sORF spans (intergenic / cds_outframe / antisense).
    Products in frame upstream of a ``true_starts`` codon (contig, plus-strand position: the
    codon's first base on ``+``, its end on ``-``) with no stop in between are
    ``nterm_extension`` instead (M7-12 decision 22)."""
    lo, hi = digest_len
    frames: list[tuple[str, str, int, bytes, int]] = []
    parts: dict[str, list[np.ndarray]] = {c: [] for c in SIXFRAME_CLASSES}
    cut = np.frombuffer(b"KR*", dtype=np.uint8)
    for cid, seq in seqs.items():
        n = len(seq)
        mask = _codon_mask(n, [(a, b, st) for c, a, b, st in coding if c == cid])
        cov = _cover(n, [(a, b, st) for c, a, b, st in annotated if c == cid])
        for strand, src in (("+", seq), ("-", revcomp(seq))):
            other = "-" if strand == "+" else "+"
            raw = src.encode()
            for off in range(3):
                aa = _translate_np(raw[off:])
                m = len(aa)
                arr = np.frombuffer(aa, dtype=np.uint8)
                stop = arr == ord("*")
                b = np.unique(np.concatenate([[0, m], np.flatnonzero(np.isin(arr, cut)) + 1,
                                              np.flatnonzero(stop)]))  # fmt: skip
                bad = np.concatenate([[0], np.cumsum(mask[strand][off::3][:m] | stop)])
                ups = [p for c, p, st in true_starts if c == cid and st == strand]
                up = _upstream_mask(m, off, n, stop, strand, ups)
                fi = len(frames)
                frames.append((cid, strand, off, aa, n))
                for k in range(1, 2 + missed_cleavages):
                    i, j = b[:-k], b[k:]
                    keep = (j - i >= lo) & (j - i <= hi) & (bad[j] == bad[i])
                    i, j = i[keep], j[keep]
                    t0, t1 = off + 3 * i, off + 3 * j
                    s0, e0 = (t0, t1) if strand == "+" else (n - t1, n - t0)
                    nterm = up[j] - up[i] == j - i
                    same = ~nterm & (cov[strand][e0] > cov[strand][s0])
                    anti = ~nterm & ~same & (cov[other][e0] > cov[other][s0])
                    rows = np.stack([np.full(len(i), fi), i, j], axis=1)
                    for cls, sel in (("cds_outframe", same), ("antisense", anti),
                                     ("intergenic", ~nterm & ~same & ~anti),
                                     ("nterm_extension", nterm)):  # fmt: skip
                        parts[cls].append(rows[sel])
    idx = {c: np.concatenate(v) if v else np.zeros((0, 3), dtype=np.int64)
           for c, v in parts.items()}  # fmt: skip
    return SixframePool(frames, idx)


def _upstream_mask(
    m: int, off: int, n: int, stop: np.ndarray, strand: str, starts: Sequence[int]
) -> np.ndarray:
    """Cumulative count (length ``m + 1``) of the codons of frame ``off`` (strand-local) that lie
    upstream of a start in ``starts`` with no stop in between."""
    mark = np.zeros(m, dtype=bool)
    stops = np.flatnonzero(stop)
    for p in starts:
        loc = p if strand == "+" else n - p
        if not 0 <= loc <= n or (loc - off) % 3:
            continue
        t = min((loc - off) // 3, m)
        k = int(np.searchsorted(stops, t)) - 1
        mark[(int(stops[k]) + 1 if k >= 0 else 0) : t] = True
    return np.concatenate([[0], np.cumsum(mark)])


def _sixframe_draw(
    rng: np.random.Generator, pool: SixframePool, n: int, ok: Callable[[Locus], bool]
) -> list[Locus]:
    """``n`` products taken in turn from each class of ``pool`` (random order within a class)
    that pass ``ok`` and whose peptide was not taken before; a class that runs out is skipped."""
    order = {c: iter(rng.permutation(len(pool.idx[c])).tolist()) for c in SIXFRAME_CLASSES}
    out: list[Locus] = []
    seen: set[str] = set()
    live = [c for c in SIXFRAME_CLASSES if len(pool.idx[c])]
    k = 0
    while len(out) < n and live:
        cls = live[k % len(live)]
        for i in order[cls]:
            loc = pool.locus(cls, i)
            if loc[0] not in seen and ok(loc):
                seen.add(loc[0])
                out.append(loc)
                k += 1
                break
        else:
            live.remove(cls)
    return out


# ---------------------------------------------------------------- simulation


@dataclass
class _Target:
    peptide: str
    wrong: bool
    origin: str  # gene | shuffled | other_strain
    loc: int  # index of the protein whose abundance sets the intensity
    ids: str  # protein_ids


def _detect(
    cfg: ProteomicsConfig, rng: np.random.Generator, prots: list[Protein], la: np.ndarray
) -> tuple[list[str], dict[str, list[tuple[int, int]]], dict[str, int]]:
    """Detected distinct peptides, every (protein, offset) occurrence and the first protein."""
    occ: dict[str, list[tuple[int, int]]] = {}
    first: dict[str, int] = {}
    for k, p in enumerate(prots):
        for a, b in digest_spans(p.aa, cfg.digest_len, cfg.missed_cleavages):
            pep = p.aa[a:b]
            if pep not in occ:
                occ[pep] = []
                first[pep] = k
            occ[pep].append((k, a))
    order = list(occ)
    b0, b1, b2 = cfg.detect_logit
    lengths = np.array([len(x) for x in order], dtype=float)
    abund = la[[first[x] for x in order]] if order else np.zeros(0)
    eta = b0 + b1 * lengths + b2 * abund
    prob = 1.0 / (1.0 + np.exp(-np.clip(eta, -700, 700)))
    hit = rng.random(len(order)) < prob
    detected = [x for x, h in zip(order, hit.tolist(), strict=True) if h]
    return detected, occ, first


def _wrong_peptides(
    cfg: ProteomicsConfig,
    rng: np.random.Generator,
    n: int,
    prots: list[Protein],
    reference: Sequence[Feat],
    hay: str,
    first: Sequence[tuple[str, str]] = (),
) -> list[tuple[str, str]]:
    """``n`` distinct wrong peptides (peptide, source locus) absent from ``hay`` (I = L); with
    ``other_strain``, those of ``first`` (M7-16 rule 20 C-terminal peptides) come first."""
    out: dict[str, str] = {}
    if n <= 0:
        return []

    def ok(pep: str) -> bool:
        return pep not in out and _il(pep) not in hay

    if cfg.wrong_source == "other_strain":
        for pep, loc in first:
            if ok(pep) and len(out) < n:
                out[pep] = loc
        pool = sorted(
            {
                (pep, f.locus)
                for f in reference
                if f.kind in ("cds", "sorf") and f.aa
                for pep in digest(f.aa, cfg.digest_len, cfg.missed_cleavages)
            }
        )
        for i in rng.permutation(len(pool)).tolist():
            if len(out) == n:
                break
            pep, loc = pool[i]
            if ok(pep):
                out[pep] = loc
        return list(out.items())
    usable = [p for p in prots if len(p.aa) >= cfg.digest_len[0]]
    tries = 0
    while usable and len(out) < n and tries < 100 * n + 1000:
        tries += 1
        p = usable[int(rng.integers(len(usable)))]
        shuffled = "".join(rng.permutation(np.array(list(p.aa))).tolist())
        prods = digest(shuffled, cfg.digest_len, cfg.missed_cleavages)
        if not prods:
            continue
        pep = prods[int(rng.integers(len(prods)))]
        if ok(pep):
            out[pep] = p.locus
    return list(out.items())


def _noncoding_pool(
    cfg: ProteomicsConfig, seqs: dict[str, str], feats: Sequence[Feat],
    novel: Sequence[Sequence[Any]],
) -> SixframePool:  # fmt: skip
    """The non-coding six-frame pool: out of frame of every CDS/sORF (true start to stop, pseudo
    genes included) and every novel ORF; classes by the annotated CDS/sORF spans, and
    ``nterm_extension`` upstream of the true start of a CDS/sORF not crossing the origin."""
    cds = [f for f in feats if f.kind in ("cds", "sorf")]
    annotated = [(f.contig, f.start, f.end, f.strand) for f in cds]
    coding = [(f.contig, *f.span(), f.strand) for f in cds] + annotated
    coding += [(str(r[2]), int(r[3]) - 1, int(r[4]), str(r[5])) for r in novel]
    starts = [(f.contig, f.span()[0] if f.strand == "+" else f.span()[1], f.strand)
              for f in cds if f.end <= len(seqs[f.contig])]  # fmt: skip
    return sixframe_noncoding_pool(seqs, coding, annotated, cfg.digest_len, cfg.missed_cleavages,
                                   starts)  # fmt: skip


#: ``giap.makedb.ORF_ID_PREFIX``, ``NTERM_PROTEIN_AA`` and ``NTERM_MIN_UPSTREAM_AA`` (restated:
#: giapsim imports no giap).
ORF_ID_PREFIX = "GIAPORF_"
NTERM_PROTEIN_AA = 30
NTERM_MIN_UPSTREAM_AA = 7  # giap.makedb.build: shorter upstream parts get no entry
_RNA_KINDS = frozenset({"rrna", "trna", "tmrna", "ncrna"})


def _json_pseudo(f: Feat) -> bool:
    """Pseudogene as GIAP's Bakta parser reads ``bakta_like.json`` (a CDS flagged ``pseudo``)."""
    return f.kind == "cds" and (bool(f.extra.get("pseudo")) or "pseudogene" in f.extra)


class OrfEntries:
    """Names the ``giap makedb --six-frame-orfs`` entry a six-frame product lies in
    (giap.makedb.build, restated). An entry is the stop-to-stop ORF of the product's frame,
    classed ``cds_outframe`` / ``antisense`` / ``intergenic`` by its overlap with the annotated
    CDS/sORFs (none when it overlaps an RNA or a same-strand pseudogene); when that ORF is in frame
    with an annotated CDS/sORF, the gene's ``nterm_extension`` entry (upstream part + its first
    30 codons; none when the upstream walk runs off a linear contig or the upstream part is
    shorter than ``NTERM_MIN_UPSTREAM_AA``)."""

    def __init__(self, feats: Sequence[Feat], seqs: dict[str, str],
                 topology: dict[str, str]) -> None:  # fmt: skip
        self.seqs = seqs
        self.rc = {c: revcomp(x) for c, x in seqs.items()}
        self.circ = {c: topology.get(c, "circular") == "circular" for c in seqs}
        genes = [f for f in feats if f.kind in ("cds", "sorf") and not _json_pseudo(f)]
        self.mask, self.cov, self.pcov, self.rcov = {}, {}, {}, {}
        self.by_stop: dict[tuple[str, str, int], Feat] = {}
        pseudo = [f for f in feats if _json_pseudo(f)]
        rna = [f for f in feats if f.kind.lower() in _RNA_KINDS]
        for c, x in seqs.items():
            n = len(x)
            spans = [(f.start, f.end, f.strand) for f in genes if f.contig == c]
            self.mask[c] = _codon_mask(n, spans)
            self.cov[c] = _cover(n, spans)
            self.pcov[c] = _cover(n, [(f.start, f.end, f.strand) for f in pseudo if f.contig == c])
            self.rcov[c] = _cover(n, [(f.start, f.end, "+") for f in rna if f.contig == c])
        for f in genes:
            n = len(seqs[f.contig])
            loc = (f.end - 3) % n if f.strand == "+" else (n - f.start - 3) % n
            self.by_stop[(f.contig, f.strand, loc)] = f

    def _orf(self, cid: str, a: int, b: int, strand: str) -> tuple[int, int, bool] | None:
        """Strand-local ``[lo, hi)`` of the stop-to-stop ORF holding ``[a, b)`` and whether its
        5' end is the end of a linear contig; None without a stop within one turn."""
        t, circ = (self.seqs if strand == "+" else self.rc)[cid], self.circ[cid]
        n = len(t)

        def stop(q: int) -> bool:
            q %= n
            return (t[q : q + 3] if q + 3 <= n else t[q:] + t[: q + 3 - n]) in STOPS

        lo, open5 = a, False
        while True:
            if not circ and lo < 3:
                open5 = True
                break
            if stop(lo - 3):
                break
            lo -= 3
            if b - lo > n:
                return None
        hi = b
        while (circ or hi + 3 <= n) and not stop(hi):
            hi += 3
            if hi - lo > n:
                return None
        return lo, hi, open5

    def _hits(self, cov: np.ndarray, s: int, e: int) -> bool:
        n = len(cov) - 1
        return any(cov[min(y, n)] > cov[x] for x, y in ((s, e), (0, e - n)) if y > x)

    def of(self, loc: Locus) -> str | None:
        """``db_id`` of the entry holding the product ``loc`` in frame, else None."""
        _, cid, s1, e1, strand = loc
        n = len(self.seqs[cid])
        a, b = (s1 - 1, e1) if strand == "+" else (n - e1, n - s1 + 1)
        orf = self._orf(cid, a, b, strand)
        if orf is None:
            return None
        lo, hi, open5 = orf
        cls = ""
        if self.mask[cid][strand][np.arange(lo, hi, 3) % n].any():  # in frame with a gene
            g = self.by_stop.get((cid, strand, hi % n))
            if g is None or open5:
                return None
            gs = lo + (((g.start if strand == "+" else n - g.end) - lo) % n)
            if gs - lo < 3 * NTERM_MIN_UPSTREAM_AA:
                return None
            hi, cls = gs + 3 * min(NTERM_PROTEIN_AA, (g.end - g.start) // 3 - 1), "nterm_extension"
        s, e = (lo, hi) if strand == "+" else (n - hi, n - lo)
        s, e = s % n, e - (s - s % n)
        if not cls:
            other = "-" if strand == "+" else "+"
            if self._hits(self.rcov[cid]["+"], s, e) or self._hits(self.pcov[cid][strand], s, e):
                return None
            same, anti = (self._hits(self.cov[cid][st], s, e) for st in (strand, other))
            cls = "cds_outframe" if same else "antisense" if anti else "intergenic"
        return f"{ORF_ID_PREFIX}{cls}_{cid}_{s + 1}_{e}_{'p' if strand == '+' else 'm'}"


def _protein_decoys(
    cfg: ProteomicsConfig,
    rng: np.random.Generator,
    n: int,
    annotated: Sequence[tuple[str, str, str]],
    target_set: set[str],
    decoys: list[tuple[str, str, int, int]],
) -> None:
    """Append ``n`` reversed Trypsin/P products of annotated proteins (``rev_<locus_tag>``) to
    ``decoys``, so gene-level decoy pairing works in six-frame mode (M7-12 decision 21)."""
    src = [(loc, aa, sp) for loc, _, aa in annotated
           if (sp := digest_spans(aa, cfg.digest_len, cfg.missed_cleavages))]  # fmt: skip
    tries = 0
    while src and len(decoys) < n and tries < 100 * n + 1000:
        tries += 1
        loc, aa, spans = src[int(rng.integers(len(src)))]
        a, b = spans[int(rng.integers(len(spans)))]
        rev = aa[a:b][::-1]
        if rev not in target_set:
            run = int(rng.integers(cfg.n_runs))
            decoys.append((rev, f"rev_{loc}", run, 2 if rng.random() < P_CHARGE_2 else 3))


def _occurrence(p: Protein, a: int, length: int, n: int) -> tuple[int, int]:
    """0-based ``[s0, e0)`` of the peptide at residue ``a`` (``length`` aa) of ``p`` on a contig
    of length ``n`` (``s0 < n``; ``e0`` may exceed ``n`` across the origin)."""
    if p.strand == "+":
        s0 = p.start + 3 * a
        e0 = s0 + 3 * length
    else:
        e0 = p.end - 3 * a
        s0 = e0 - 3 * length
    if s0 >= n:
        s0, e0 = s0 - n, e0 - n
    return s0, e0


def _place(p: Protein, a: int, pep: str, seqs: dict[str, str]) -> tuple[int, int, bool]:
    """0-based ``[s0, e0)`` of ``pep`` at residue ``a`` of ``p`` (folded into the contig) and
    whether it starts at the initiator with a non-ATG start codon read as M."""
    s0, e0 = _occurrence(p, a, len(pep), len(seqs[p.contig]))
    codon_aa = _start_codon_aa(seqs[p.contig], s0, e0, p.strand)
    return s0, e0, a == 0 and pep[0] == "M" and codon_aa != "M"


def _search_space_share(
    pool: SixframePool, annotated: Sequence[tuple[str, str, str]], cfg: ProteomicsConfig
) -> tuple[float, int, int]:
    """(protein share, distinct protein-entry products, distinct six-frame-only products) of the
    simulated makedb search space (MVP-04 Q2; information only, lead decision 10)."""
    prot = {x.encode() for _, _, aa in annotated
            for x in digest(aa, cfg.digest_len, cfg.missed_cleavages)}  # fmt: skip
    six: set[bytes] = set()
    for arr in pool.idx.values():
        for f, a0, b0 in arr.tolist():
            six.add(pool.frames[f][3][a0:b0])
    six -= prot
    n_p, n_s = len(prot), len(six)
    return (n_p / (n_p + n_s) if n_p + n_s else 0.5), n_p, n_s


def _protein_nulls(
    cfg: ProteomicsConfig,
    rng: np.random.Generator,
    n: int,
    prots: list[Protein],
    genes: set[str],
    detected: Sequence[str],
    occ: dict[str, list[tuple[int, int]]],
    seqs: dict[str, str],
) -> tuple[list[tuple[str, str]], dict[str, list[tuple[str, str, int, int, str, bool]]]]:
    """``n`` protein-entry nulls (lead decision 9): ``[(peptide, locus)]`` and their placements.
    ``prots`` may end with recorded loci the sample strain made pseudo (no ``occ`` entries);
    products across their premature stop are skipped."""
    seen = {_il(p) for p in detected}
    hit = {prots[k].locus for p in detected for k, _ in occ[p]}
    pool = []
    for k, p in enumerate(prots):
        if p.locus not in genes or p.locus in hit:
            continue
        spans = [(a, b) for a, b in digest_spans(p.aa, cfg.digest_len, cfg.missed_cleavages)
                 if "*" not in p.aa[a:b]]  # fmt: skip
        prods = sorted({p.aa[a:b] for a, b in spans if _il(p.aa[a:b]) not in seen})
        if prods:
            pool.append((k, prods, spans))
    out: list[tuple[str, str]] = []
    places: dict[str, list[tuple[str, str, int, int, str, bool]]] = {}
    used: set[str] = set()
    for i in rng.permutation(len(pool)).tolist() if n else []:
        k, prods, spans = pool[i]
        free = [x for x in prods if _il(x) not in used]
        if not free:
            continue
        pep = free[int(rng.integers(len(free)))]
        used.add(_il(pep))
        p = prots[k]
        out.append((pep, p.locus))
        rows = []
        for a, b in spans:
            if p.aa[a:b] == pep:
                s0, e0, alt = _place(p, a, pep, seqs)
                rows.append((p.locus, p.contig, s0, e0, p.strand, alt))
        places[pep] = rows
        if len(out) == n:
            break
    if len(out) < n:
        _LOG.warning("proteomics: only %d protein-entry nulls for %d requested", len(out), n)
    return out, places


def _scored_nulls(
    cfg: ProteomicsConfig,
    seed: int,
    xrng: np.random.Generator,
    prots: list[Protein],
    detected: Sequence[str],
    occ: dict[str, list[tuple[int, int]]],
    feats: Sequence[Feat],
    novel: Sequence[Sequence[Any]],
    seqs: dict[str, str],
    topology: dict[str, str],
    hay: str,
    annotated: Sequence[tuple[str, str, str]],
    n_wrong: int,
) -> tuple[
    list[tuple[str, str]],
    dict[str, str],
    dict[str, Locus],
    dict[str, list[tuple[str, str, int, int, str, bool]]],
    dict[str, Any],
]:
    """The ``scored_null`` null peptides: ``[(peptide, entry ID)]``, origins, six-frame loci,
    protein-null placements and the ``null_split`` record (module doc)."""
    pool = _noncoding_pool(cfg, seqs, feats, novel)
    share, n_prod_p, n_prod_s = _search_space_share(pool, annotated, cfg)
    n_p = n_wrong // 2  # lead decision 10: fixed 50/50, capped by the protein-null pool
    genes = {f.locus for f in feats if f.kind in ("cds", "sorf")}
    # recorded loci the sample strain made pseudo are undetected recorded loci too (decision 9)
    strain_pseudo = [
        Protein(f.locus, f.contig, *f.span(), f.strand,
                translate_cds(codons_of(_slice(seqs[f.contig], *f.span(), f.strand))))
        for f in feats
        if f.kind in ("cds", "sorf") and not f.pseudo_reason and f.extra.get("pseudo")
    ]  # fmt: skip
    nrng = np.random.default_rng([seed, _NULL_STREAM])
    prot_nulls, places = _protein_nulls(
        cfg, nrng, n_p, [*prots, *strain_pseudo], genes, detected, occ, seqs
    )
    orfs = OrfEntries(feats, seqs, topology)
    entry: dict[str, str] = {}

    def ok(loc: Locus) -> bool:
        db_id = orfs.of(loc) if _il(loc[0]) not in hay else None
        if db_id is not None:
            entry[loc[0]] = db_id
        return db_id is not None

    drawn = _sixframe_draw(xrng, pool, n_wrong - len(prot_nulls), ok)
    loci = {loc[0]: loc for loc in drawn}
    wrong = prot_nulls + [(w, entry[w]) for w in loci]
    origins = {w: PROTEIN_NULL_ORIGIN for w, _ in prot_nulls}
    origins |= {w: "sixframe_noncoding" for w in loci}
    split = {"mode": "fifty_fifty", "search_space_protein_share": round(share, 6),
             "n_protein_products": n_prod_p, "n_sixframe_products": n_prod_s,
             "n_protein": len(prot_nulls), "n_sixframe": len(loci)}  # fmt: skip
    return wrong, origins, loci, places, split


def _start_codon_aa(seq: str, s0: int, e0: int, strand: str) -> str:
    """Amino acid that code 11 assigns to the first codon of ``[s0, e0)`` on ``strand`` (the span
    may run past the end of a circular contig)."""
    n = len(seq)
    a = s0 if strand == "+" else e0 - 3
    a %= n
    return CODE_11[_slice(seq, a, a + 3, strand)]


def _fmt(x: float) -> str:
    return f"{x:.10g}"


def _copies(
    pc: ProteinSimConfig, rng: np.random.Generator, annotated: list[tuple[str, str, str]]
) -> list[tuple[str, str]]:
    """Cross-strain copies (UniProt-like ID, sequence) at identity drawn from ``pc.identity``."""
    src = [
        (loc, aa) for loc, kind, aa in annotated if kind == "cds" and len(aa) >= pc.min_length_aa
    ]
    k = round(pc.fraction * len(src))
    if src and pc.fraction > 0:
        k = max(k, 1)
    picks = sorted(rng.choice(len(src), size=k, replace=False).tolist()) if k else []
    out, used = [], set()
    alnum = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
    for i in picks:
        aa = src[i][1]
        ident = float(rng.uniform(*pc.identity))
        n_sub = min(max(1, round((1 - ident) * len(aa))), int(0.2 * len(aa)))
        seq = list(aa)
        for j in rng.choice(np.arange(1, len(aa)), size=n_sub, replace=False).tolist():
            alt = [c for c in AA20 if c != seq[j]]
            seq[j] = alt[int(rng.integers(len(alt)))]
        while True:
            acc = (
                "OPQ"[int(rng.integers(3))]
                + str(int(rng.integers(10)))
                + "".join(alnum[int(x)] for x in rng.integers(36, size=3))
                + str(int(rng.integers(10)))
            )
            if acc not in used:
                used.add(acc)
                break
        out.append((f"tr|{acc}|{acc}_SIMUL", "".join(seq)))
    return out


def _write_fasta(path: Path, records: Sequence[tuple[str, str]]) -> None:
    parts = []
    for name, seq in records:
        parts.append(f">{name}\n" + "".join(seq[i : i + 60] + "\n" for i in range(0, len(seq), 60)))
    path.write_text("".join(parts))


def _write_tsv(path: Path, header: Sequence[str], rows: Sequence[Sequence[str]]) -> None:
    lines = ["\t".join(header)] + ["\t".join(r) for r in rows]
    path.write_text("\n".join(lines) + "\n")


PEPTIDE_HEADER = (
    "peptide", "protein_ids", "q_value", "pep", "score", "is_decoy", "intensity", "run", "charge",
    "spectrum_id", "source", "search_db",
)  # fmt: skip
TRUTH_HEADER = (
    "peptide", "is_wrong", "origin", "locus_tag", "contig", "start", "end", "strand",
    "n_loci_genome", "alt_start_nterm",
)  # fmt: skip


def simulate_proteomics(
    cfg: ProteomicsConfig,
    seed: int,
    out: Path,
    feats: Sequence[Feat],
    novel: Sequence[Sequence[Any]],
    seqs: dict[str, str],
    topology: dict[str, str],
    reference: Sequence[Feat] | None = None,
    log_abundance: np.ndarray | None = None,
    first_wrong: Sequence[tuple[str, str]] = (),
    readthrough: Mapping[str, tuple[str, str]] | None = None,
) -> dict[str, Path]:
    """Write ``peptides.tsv``, ``proteins.faa``, ``abundance.tsv`` and ``truth/peptides.tsv``.

    ``log_abundance`` (MVP-04, opt-in) gives the natural-log abundance of every expressed protein
    in ``expressed_proteins`` order instead of drawing ``gene_log_abundance``; ``None`` keeps the
    draw (and so every existing output byte for byte)."""
    if cfg.wrong_source == "other_strain" and reference is None:
        raise ValueError("wrong_source='other_strain' needs SimConfig.strain")
    rng = np.random.default_rng([seed, _STREAM])
    prots = expressed_proteins(feats, novel, seqs, readthrough)
    if log_abundance is None:
        la = gene_log_abundance(rng, len(prots))
    else:
        la = np.asarray(log_abundance, dtype=float)
        if la.shape != (len(prots),):
            raise ValueError(f"log_abundance has shape {la.shape}, expected ({len(prots)},)")
    detected, occ, first = _detect(cfg, rng, prots, la)
    # makedb target IDs (decision 9): six-frame wrong peptides, or the scored_null search
    six = cfg.wrong_source == "sixframe_noncoding" or cfg.decoys == "scored_null"
    orfs = OrfEntries(feats, seqs, topology) if six else None
    novel_ids = {str(r[0]) for r in novel}
    # codons of a true protein upstream of its annotated start (N-terminal extensions)
    upstream = {f.locus: ((f.start - f.span()[0]) if f.strand == "+" else (f.span()[1] - f.end))
                // 3 for f in feats if f.kind in ("cds", "sorf")}  # fmt: skip

    def db_id(pep: str, k: int, a: int) -> str | None:
        """The protein ID a search reports for occurrence ``(k, a)``. In six-frame and
        scored_null mode (decision 9) a novel ORF's peptide names its makedb entry, and an
        annotated protein's peptide reaching upstream of the annotated start names the
        gene's ``nterm_extension`` entry (either None when makedb has no such entry: the
        occurrence is not in the search database); any other peptide of an annotated protein
        names its locus tag, makedb's ``db_id``."""
        p = prots[k]
        if orfs is None or (p.locus not in novel_ids and a >= upstream.get(p.locus, 0)):
            return p.locus
        s0, e0 = _occurrence(p, a, len(pep), len(seqs[p.contig]))
        return orfs.of((pep, p.contig, s0 + 1, e0, p.strand))

    found = {pep: [d for k, a in occ[pep] if (d := db_id(pep, k, a)) is not None]
           for pep in detected}  # fmt: skip
    detected = [pep for pep in detected if found[pep]]  # a search cannot report the others
    annotated = [(f.locus, f.kind, f.aa) for f in feats if f.kind in ("cds", "sorf") and f.aa]
    annotated_loci = {loc for loc, _, _ in annotated}
    hay = "|".join(_il(p.aa) for p in prots) + "|" + "|".join(_il(a) for _, _, a in annotated)
    pi = cfg.pi_wrong
    n_wrong = round(len(detected) * pi / (1 - pi))
    xrng = np.random.default_rng([seed, _SIXFRAME_STREAM])  # six-frame draws (M7-12)
    scored = cfg.decoys == "scored_null"
    six = six and not scored  # the six-frame wrong-peptide draw follows wrong_source alone
    pool: SixframePool | None = None
    loci: dict[str, Locus] = {}
    entry: dict[str, str] = {}  # six-frame wrong peptide -> makedb entry db_id (target + decoy)
    origins: dict[str, str] = {}  # wrong peptide -> truth origin (scored_null)
    places: dict[str, list[tuple[str, str, int, int, str, bool]]] = {}
    null_split: dict[str, Any] | None = None
    if scored:
        nulls = _scored_nulls(cfg, seed, xrng, prots, detected, occ, feats, novel, seqs,
                              topology, hay, annotated, n_wrong)  # fmt: skip
        wrong, origins, loci, places, null_split = nulls
        entry = dict(wrong)
    elif six and orfs is not None:
        pool = _noncoding_pool(cfg, seqs, feats, novel)

        def ok(loc: Locus) -> bool:
            db_id = orfs.of(loc) if _il(loc[0]) not in hay else None
            if db_id is not None:
                entry[loc[0]] = db_id
            return db_id is not None

        drawn = _sixframe_draw(xrng, pool, n_wrong, ok)
        loci = {loc[0]: loc for loc in drawn}
        wrong = [(w, f"sixframe:{c}:{a}-{b}:{st}") for w, c, a, b, st in loci.values()]
    else:
        if cfg.wrong_source == "other_strain" and n_wrong < len(first_wrong):
            raise ValueError(f"{n_wrong} wrong-peptide slots for {len(first_wrong)} stop-losses")
        wrong = _wrong_peptides(cfg, rng, n_wrong, prots, reference or (), hay, first_wrong)
    # PEP and q follow the wrong share actually produced (spec M3-19, review decision 1)
    n_all = len(detected) + len(wrong)
    pi_realised = len(wrong) / n_all if n_all else 0.0
    if len(wrong) < n_wrong:
        _LOG.warning(
            "proteomics: wrong-peptide shortfall: pi_wrong requested %.3f, realised %.3f "
            "(%d of %d wrong peptides)",
            pi,
            pi_realised,
            len(wrong),
            n_wrong,
        )
    origin = cfg.wrong_source
    targets = [
        _Target(
            pep,
            False,
            "gene",
            first[pep],
            ";".join(dict.fromkeys(found[pep])),
        )
        for pep in detected
    ]
    targets += [
        _Target(
            pep, True, origins.get(pep, origin), int(rng.integers(len(prots))), entry.get(pep, loc)
        )
        for pep, loc in wrong
    ]

    # target rows: every target peptide is identified once in every run
    n_t, n_runs = len(targets), cfg.n_runs
    mu1, sd1, mu0, sd0 = cfg.score_mixture
    is_wrong = np.array([t.wrong for t in targets], dtype=bool)
    z = rng.standard_normal((n_runs, n_t))
    scores = np.round(np.where(is_wrong, mu0 + sd0 * z, mu1 + sd1 * z), 6)
    abund = la[[t.loc for t in targets]] if n_t else np.zeros(0)
    intensity = np.exp(abund + INTENSITY_SD * rng.standard_normal((n_runs, n_t)))
    charge = np.where(rng.random((n_runs, n_t)) < P_CHARGE_2, 2, 3)
    peps = mixture_pep(scores.ravel(), pi_realised, cfg.score_mixture).reshape(scores.shape)
    qs = mixture_q(scores.ravel(), pi_realised, cfg.score_mixture).reshape(scores.shape)

    # decoy rows: reversed tryptic peptides of expressed proteins, null scores
    n_dec = int(is_wrong.sum()) * n_runs if cfg.n_decoy_rows == "match_wrong" else cfg.n_decoy_rows
    target_set = {t.peptide for t in targets}
    decoys: list[tuple[str, str, int, int]] = []  # peptide, protein_ids, run, charge
    usable = [
        k for k, p in enumerate(prots) if digest_spans(p.aa, cfg.digest_len, cfg.missed_cleavages)
    ]
    tries = 0
    to_decoy = np.zeros((n_runs, n_t), dtype=bool)
    collisions = 0
    if scored:  # MVP-04 rule 6: a 1/2 coin per null PSM row, no other decoys
        drng = np.random.default_rng([seed, _DECOY_STREAM])
        coin = is_wrong & (drng.random((n_runs, n_t)) < 0.5)
        target_il = {_il(t.peptide) for t in targets}
        collide = np.array([t.wrong and _il(t.peptide[::-1]) in target_il for t in targets],
                           dtype=bool)  # fmt: skip
        to_decoy = coin & ~collide
        collisions = int((coin & collide).sum())
        usable = []
    elif six:  # decision 21: a 1/2 coin per null PSM, protein decoys from annotated proteins only
        drng = np.random.default_rng([seed, _DECOY_STREAM])
        to_decoy = is_wrong & (drng.random((n_runs, n_t)) < 0.5)
        pseudo_loci = {f.locus for f in feats if f.pseudo_reason or f.extra.get("pseudo")}
        real = [x for x in annotated if x[0] not in pseudo_loci]
        _protein_decoys(cfg, drng, n_dec, real, target_set, decoys)
        usable = []
    while usable and len(decoys) < n_dec and tries < 100 * n_dec + 1000:
        tries += 1
        p = prots[usable[int(rng.integers(len(usable)))]]
        spans = digest_spans(p.aa, cfg.digest_len, cfg.missed_cleavages)
        a, b = spans[int(rng.integers(len(spans)))]
        rev = p.aa[a:b][::-1]
        if rev in target_set:
            continue
        run = int(rng.integers(n_runs))
        decoys.append((rev, f"rev_{p.locus}", run, 2 if rng.random() < P_CHARGE_2 else 3))

    rows: list[list[str]] = []
    for r in range(n_runs):
        run_rows: list[list[str]] = []
        for i, t in enumerate(targets):
            if to_decoy[r, i]:  # the null PSM went to the decoy of its makedb entry
                run_rows.append([
                    t.peptide[::-1], f"rev_{entry[t.peptide]}", _fmt(qs[r, i]),
                    _fmt(peps[r, i]), _fmt(scores[r, i]), "true", _NA, f"run{r + 1}",
                    str(charge[r, i]),
                ])  # fmt: skip
                continue
            run_rows.append([
                t.peptide, t.ids, _fmt(qs[r, i]), _fmt(peps[r, i]), _fmt(scores[r, i]), "false",
                f"{intensity[r, i]:.6g}", f"run{r + 1}", str(charge[r, i]),
            ])  # fmt: skip
        for pep, ids, run, ch in decoys:
            if run == r:
                run_rows.append([pep, ids, _NA, _NA, _NA, "true", _NA, f"run{r + 1}", str(ch)])
        for j, i in enumerate(rng.permutation(len(run_rows)).tolist()):
            rows.append([*run_rows[i], f"run{r + 1}.{j + 1}", "database", cfg.search_db])

    # truth
    counts = sixframe_counts([*detected, *(w for w, _ in wrong)], seqs, topology)
    truth: list[list[str]] = []
    for pep in detected:
        n_loci = str(counts[pep])
        for k, a in occ[pep]:
            p = prots[k]
            s0, e0, alt = _place(p, a, pep, seqs)
            truth.append([
                pep, "false", "gene", p.locus, p.contig, str(s0 + 1), str(e0), p.strand, n_loci,
                "true" if alt else "false",
            ])  # fmt: skip
    for pep, _loc in wrong:
        if pep in places:  # scored_null protein nulls: placed like correct rows
            for loc_tag, c, s0, e0, st, alt in places[pep]:
                truth.append([pep, "true", PROTEIN_NULL_ORIGIN, loc_tag, c, str(s0 + 1), str(e0),
                              st, str(counts[pep]), "true" if alt else "false"])  # fmt: skip
            continue
        if pep in loci:
            _, c, a, b, st = loci[pep]
            truth.append([pep, "true", origins.get(pep, origin), _NA, c, str(a), str(b), st,
                          str(counts[pep]), "false"])  # fmt: skip
            continue
        truth.append([pep, "true", origin, _NA, _NA, _NA, _NA, _NA, str(counts[pep]), "false"])
    six_decoys = [t.peptide[::-1] for i, t in enumerate(targets) if to_decoy[:, i].any()]
    for pep in dict.fromkeys([*(d[0] for d in decoys), *six_decoys]):
        truth.append([pep, "true", "decoy", _NA, _NA, _NA, _NA, _NA, _NA, _NA])

    pseudo = {f.locus for f in feats if f.pseudo_reason or f.extra.get("pseudo")}
    faa = [(loc, aa) for loc, _, aa in annotated if loc not in pseudo]
    if cfg.proteins is not None:
        faa += [
            c
            for c in _copies(cfg.proteins, rng, [x for x in annotated if x[0] not in pseudo])
            if c[0] not in annotated_loci
        ]

    paths = {
        "peptides": out / "peptides.tsv",
        "proteins": out / "proteins.faa",
        "abundance": out / "abundance.tsv",
        "truth": out / "truth" / "peptides.tsv",
        "truth_meta": out / "truth" / "proteomics.json",
    }
    _write_tsv(paths["peptides"], PEPTIDE_HEADER, rows)
    _write_fasta(paths["proteins"], faa)
    _write_tsv(
        paths["abundance"],
        ("locus_tag", "log_abundance"),
        [[p.locus, f"{x:.6f}"] for p, x in zip(prots, la.tolist(), strict=True)],
    )
    _write_tsv(paths["truth"], TRUTH_HEADER, truth)
    meta: dict[str, Any] = {
        "pi_wrong": pi,
        "pi_wrong_realised": round(pi_realised, 10),
        "n_correct": len(detected),
        "n_wrong": len(wrong),
        "n_wrong_requested": n_wrong,
    }
    if scored:
        meta |= {"decoys": cfg.decoys, "decoy_collisions": collisions, "null_split": null_split}
    paths["truth_meta"].write_text(json.dumps(meta, indent=2) + "\n")
    return paths
