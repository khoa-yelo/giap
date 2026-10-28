"""Genome features for mapping (issue M3-18; spec §11.3): RNA genes, planted genes, novel ORFs.

Planted features are built as self-contained *blocks* (a sequence on the ``+`` strand plus the
features in it, block-local coordinates) that the layout places between transcription units.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING, Any

import numpy as np

from gippasim.codons import STOPS, code, codons_of, revcomp, stops, translate_cds

if TYPE_CHECKING:
    from gippasim.simulate import _Sim

_LOG = logging.getLogger(__name__)

ANTI_SD_TAIL = "GATCACCTCCTTA"  # 3' end of 16S rRNA (anti-Shine-Dalgarno)
NOVEL_CLASSES = ("intergenic", "cds_outframe", "antisense")
STARTS = frozenset({"ATG", "GTG", "TTG"})
_RRNA = (("16S", "rrs", 1542), ("23S", "rrl", 2904), ("5S", "rrf", 120))
_ITS = ((300, 500), (80, 150))  # 16S-23S and 23S-5S spacers
_AA3 = dict(
    zip(
        "ACDEFGHIKLMNPQRSTVWY",
        re.findall("...", "AlaCysAspGluPheGlyHisIleLysLeuMetAsnProGlnArgSerThrValTrpTyr"),
        strict=True,
    )
)
_MIN_NOVEL_AA = 10


@dataclass
class Feat:
    """An annotated feature; 0-based half-open, ``end`` > contig length when crossing the origin."""

    kind: str  # JSON type: cds | sorf | rRNA | tRNA
    contig: str = ""
    start: int = 0
    end: int = 0
    strand: str = "+"
    nt: str = ""
    aa: str = ""
    extra: dict[str, Any] = field(default_factory=dict)  # further JSON keys
    true_start: int | None = None  # 1-based first base of the true start codon
    pseudo_reason: str = ""
    locus: str = ""
    tu: str = ""  # transcription unit of a laid-out gene (truth/genes.tsv); "" for planted ones
    sec_residues: tuple[int, ...] = ()  # annotated residues expressed as U (M7-16c rule 23)

    def placed(self, contig: str, offset: int) -> Feat:
        ts = None if self.true_start is None else self.true_start + offset
        return replace(
            self, contig=contig, start=self.start + offset, end=self.end + offset, true_start=ts
        )

    def span(self) -> tuple[int, int]:
        """Interval including an N-terminal extension."""
        if self.true_start is None:
            return self.start, self.end
        return (
            (self.true_start - 1, self.end) if self.strand == "+" else (self.start, self.true_start)
        )


@dataclass
class Block:
    seq: str  # + strand
    feats: list[Feat]
    kind: str = ""  # v2 blocks: crispr | oric | bgc | prophage


def _block(parts: list[tuple[str, Feat | None]], strand: str) -> Block:
    """Assemble parts given in transcription order into a block on ``strand``."""
    seq = "".join(p for p, _ in parts)
    n, cur, feats = len(seq), 0, []
    for p, f in parts:
        if f is not None:
            f.start, f.end = (cur, cur + len(p)) if strand == "+" else (n - cur - len(p), n - cur)
            f.strand = strand
            feats.append(f)
        cur += len(p)
    return Block(seq if strand == "+" else revcomp(seq), feats)


def random_seq(rng: np.random.Generator, n: int, gc: float = 0.5) -> str:
    idx = rng.choice(4, size=n, p=[(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2])
    return np.frombuffer(b"ACGT", dtype=np.uint8)[idx].tobytes().decode()


def coding(kind: str, codons: list[str], aa: str | None = None, table: int = 11, **kw: Any) -> Feat:
    extra = {"start_type": codons[0], **kw.pop("extra", {})}
    nt = kw.pop("nt", "".join(codons))
    return Feat(kind, nt=nt, aa=aa or translate_cds(codons, table), extra=extra, **kw)


def rrna_operon(sim: _Sim) -> list[tuple[str, Feat | None]]:
    rng = sim.rng
    parts: list[tuple[str, Feat | None]] = []
    for i, (name, gene, length) in enumerate(_RRNA):
        seq = random_seq(rng, length, 0.53)
        if name == "16S":
            seq = seq[: -len(ANTI_SD_TAIL)] + ANTI_SD_TAIL
        parts.append((seq, Feat("rRNA", extra={"gene": gene, "product": f"{name} ribosomal RNA"})))
        if i < len(_ITS):
            parts.append((random_seq(rng, int(rng.integers(*_ITS[i])), 0.45), None))
    return parts


def trna(sim: _Sim) -> list[tuple[str, Feat | None]]:
    rng = sim.rng
    codon = sorted(sim.usage)[int(rng.integers(len(sim.usage)))]
    ac, aa3 = revcomp(codon), _AA3[code(sim.table)[codon]]
    seq = random_seq(rng, 33, 0.55) + ac + random_seq(rng, 37, 0.55) + "CCA"
    extra = {"anti_codon": ac, "amino_acid": aa3, "product": f"tRNA-{aa3}"}
    return [(seq, Feat("tRNA", extra=extra))]


def sorf(sim: _Sim) -> list[tuple[str, Feat | None]]:
    n_aa = int(sim.rng.integers(10, 30))  # < 30 aa
    codons = [sim.start_codon(), *sim.sense(n_aa - 1), sim.stop_codon()]
    return [("".join(codons), coding("sorf", codons, table=sim.table))]


def pseudogene(sim: _Sim, reason: str) -> list[tuple[str, Feat | None]]:
    """A gene with an internal stop or a 1-nt deletion; ``aa`` is the intact ancestor's."""
    cod = sim.gene().codons
    aa = translate_cds(cod, sim.table)
    i = int(sim.rng.integers(len(cod) // 4 + 1, 3 * len(cod) // 4))
    if reason == "internal_stop":
        cod[i] = sim.stop_codon()
        nt = "".join(cod)
    else:
        nt = "".join(cod)
        nt = nt[: 3 * i + 1] + nt[3 * i + 2 :]
    f = coding("cds", cod, aa, nt=nt, extra={"pseudo": True}, pseudo_reason=reason)
    return [(nt, f)]


def nterm_extension(sim: _Sim) -> list[tuple[str, Feat | None]]:
    """A gene whose annotated start (ATG) lies 3·k nt downstream of the true start."""
    cod = sim.gene().codons
    k = int(sim.rng.integers(5, min(30, len(cod) - 10) + 1))
    cod[k] = "ATG"
    return [("".join(cod[:k]), None), ("".join(cod[k:]), coding("cds", cod[k:], table=sim.table))]


def make_blocks(sim: _Sim) -> list[Block]:
    from gippasim import blocks as v2  # simulator v2 (M2-05)

    cfg, rng = sim.cfg, sim.rng
    rrna: Callable[[], list[tuple[str, Feat | None]]] = lambda: rrna_operon(sim)  # noqa: E731
    trnas: list[Callable[[], list[tuple[str, Feat | None]]]] = [lambda: trna(sim)] * cfg.n_trna
    if cfg.real_rrna_trna:
        rrna = lambda: v2.real_rrna_operon(sim)  # noqa: E731
        trnas = [lambda i=i: v2.real_trna(i) for i in v2.trna_picks(rng, cfg.n_trna)]  # type: ignore[misc]
    makers: list[Callable[[], list[tuple[str, Feat | None]]]] = [
        *[rrna] * cfg.n_rrna_operons,
        *trnas,
        *[lambda: sorf(sim)] * cfg.n_sorfs,
        *[lambda: nterm_extension(sim)] * cfg.n_nterm_extensions,
    ]
    r0 = int(rng.integers(2))
    reasons = [("internal_stop", "frameshift")[(i + r0) % 2] for i in range(cfg.n_pseudogenes)]
    makers += [lambda r=r: pseudogene(sim, r) for r in reasons]  # type: ignore[misc]
    kinds = [""] * len(makers)
    for kind, maker in v2.v2_makers(sim):
        makers.append(maker)
        kinds.append(kind)
    blocks = []
    for j in rng.permutation(len(makers)):
        strand = "+" if rng.random() < 0.5 else "-"
        b = _block(makers[int(j)](), strand)
        b.kind = kinds[int(j)]
        for f in b.feats:
            if f.kind in ("crispr", "oric"):
                f.strand = "."
        for f in b.feats:  # N-terminal extension: true start at the block's 5' end
            if f.kind == "cds" and not f.pseudo_reason and not b.kind and len(f.nt) < len(b.seq):
                f.true_start = 1 if strand == "+" else len(b.seq)
        blocks.append(b)
    return blocks


# ---------------------------------------------------------------- novel ORFs


def _scan(
    seq: str, frames: tuple[int, ...], min_aa: int = _MIN_NOVEL_AA, stop_set: frozenset[str] = STOPS
) -> tuple[int, int] | None:
    """First start..stop ORF (>= ``min_aa`` aa) in ``frames`` of ``seq``: (offset, length)."""
    for fr in frames:
        open_at = None
        for j, c in enumerate(codons_of(seq[fr:])):
            if c in stop_set:
                if open_at is not None and j - open_at >= min_aa:
                    return fr + 3 * open_at, 3 * (j - open_at + 1)
                open_at = None
            elif open_at is None and c in STARTS:
                open_at = j
    return None


def _orf_row(
    cls: str, contig: str, s0: int, e0: int, strand: str, nt: str, table: int = 11
) -> list[Any]:
    return [cls, contig, s0 + 1, e0, strand, translate_cds(codons_of(nt), table)]


def novel_orfs(
    sim: _Sim,
    counts: Mapping[str, int],
    seqs: dict[str, bytearray],
    feats: list[Feat],
    len_aa: tuple[int, int] = (_MIN_NOVEL_AA, 30),
    margin: int = 20,
) -> list[list[Any]]:
    """Plant intergenic ORFs and find out-of-frame/antisense ORFs inside CDSs (none annotated).

    Planted ORFs have ``len_aa`` residues (closed range) and keep ``margin`` bp from the ends of
    their gap; found ORFs have at least ``len_aa[0]`` residues (M7-12; defaults as M3-18). With
    fewer usable gaps than requested intergenic ORFs, one is planted per gap and a WARNING is
    logged (M7-12 decision 15; the truth table records the realised count)."""
    rng, rows = sim.rng, []
    lo, hi = len_aa
    # +1: rng.integers(margin, gap - len - margin) needs a non-empty range (decision 15)
    min_gap = max(160, 3 * (hi + 1) + 2 * margin + 1)
    lengths = {c: len(s) for c, s in seqs.items()}
    known = {f.aa for f in feats if f.kind in ("cds", "sorf")}
    chrom = next(iter(seqs))
    spans = sorted(p for f in feats if f.contig == chrom for p in _split(*f.span(), lengths[chrom]))
    gaps, last = [], 0
    for s, e in spans:
        if s - last >= min_gap:
            gaps.append((last, s))
        last = max(last, e)
    n = counts.get("intergenic", 0)
    if n > len(gaps):
        _LOG.warning("novel ORFs: only %d intergenic gaps for %d intergenic ORFs; planting %d",
                     len(gaps), n, len(gaps))  # fmt: skip
        n = len(gaps)
    for gi in sorted(rng.choice(len(gaps), size=n, replace=False).tolist()) if n else []:
        a, b = gaps[gi]
        n_aa = int(rng.integers(lo, hi + 1))
        nt = "".join([sim.start_codon(), *sim.sense(n_aa - 1), sim.stop_codon()])
        strand = "+" if rng.random() < 0.5 else "-"
        s0 = a + int(rng.integers(margin, b - a - len(nt) - margin))
        seqs[chrom][s0 : s0 + len(nt)] = (nt if strand == "+" else revcomp(nt)).encode()
        rows.append(_orf_row("intergenic", chrom, s0, s0 + len(nt), strand, nt, sim.table))
    cds = [
        f for f in feats if f.kind == "cds" and not f.pseudo_reason and f.end <= lengths[f.contig]
    ]
    for cls, frames, anti in (("cds_outframe", (1, 2), False), ("antisense", (0, 1, 2), True)):
        want = counts.get(cls, 0)
        for k in rng.permutation(len(cds)).tolist() if want else []:
            f = cds[k]
            seq = revcomp(f.nt) if anti else f.nt
            hit = _scan(seq, frames, lo, stops(sim.table))
            if hit is None:
                continue
            o, ln = hit
            nt = seq[o : o + ln]
            if translate_cds(codons_of(nt), sim.table) in known:
                continue
            on_plus = (f.strand == "+") != anti  # orientation of ``seq`` on the genome
            s0, e0 = (f.start + o, f.start + o + ln) if on_plus else (f.end - o - ln, f.end - o)
            rows.append(_orf_row(cls, f.contig, s0, e0, "+" if on_plus else "-", nt, sim.table))
            want -= 1
            if not want:
                break
        if want:
            raise ValueError(f"could not find enough {cls} ORFs")
    return [[f"NOVEL_{i + 1:04d}", *r] for i, r in enumerate(rows)]


def _split(s: int, e: int, length: int) -> list[tuple[int, int]]:
    return [(s, e)] if e <= length else [(s, length), (0, e - length)]
