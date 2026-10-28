"""Planted blocks of simulator v2 (issue M2-05; GIPPA_SPEC.md §11.3).

Real E. coli K-12 MG1655 rRNA/tRNA genes (``data/``, see ``data/PROVENANCE.md``), a CRISPR array,
an oriC-like region next to ``mnmG``, a BGC-like gene block and a prophage gene block whose GC
is shifted relative to the rest of the chromosome. The prophage block is laid out as a
placeholder and filled (``fill_prophage``) once the rest of the chromosome is known.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from functools import cache
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from gippasim.codons import SENSE_USAGE
from gippasim.features import Feat, coding, random_seq
from gippasim.signals import plant_gene_rbs

if TYPE_CHECKING:
    from gippasim.simulate import _Sim

DATA = Path(__file__).parent / "data"
Parts = list[tuple[str, Feat | None]]
_ITS = ((300, 500), (80, 150))
_RRNA = {"rrsB": ("rrs", "16S"), "rrlB": ("rrl", "23S"), "rrfB": ("rrf", "5S")}
_MNMG = "tRNA uridine 5-carboxymethylaminomethyl modification enzyme MnmG"
_BLOCK_SPACER = (100, 200)  # between block genes: far apart, like genes of different TUs
_GC_TOLERANCE = 0.004
_GC_TRIES = 40


def read_fasta(path: Path) -> list[tuple[str, str]]:
    recs: list[tuple[str, list[str]]] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            recs.append((line[1:], []))
        elif line.strip():
            recs[-1][1].append(line.strip())
    return [(h, "".join(s)) for h, s in recs]


@cache
def bundled() -> tuple[tuple[tuple[str, str], ...], tuple[tuple[str, str, str, str], ...]]:
    """(rRNA (name, seq) in operon order, tRNA (name, seq, amino acid, anticodon))."""
    rrna = tuple((h.split()[0], s) for h, s in read_fasta(DATA / "ecoli_rrnB.fa"))
    trna = []
    for h, s in read_fasta(DATA / "ecoli_trna.fa"):
        tags = dict(t.split("=", 1) for t in h.split()[2:])
        trna.append((h.split()[0], s, tags["amino_acid"], tags["anticodon"]))
    return rrna, tuple(trna)


def real_rrna_operon(sim: _Sim) -> Parts:
    parts: Parts = []
    rrna, _ = bundled()
    for i, (name, seq) in enumerate(rrna):
        gene, kind = _RRNA[name]
        parts.append((seq, Feat("rRNA", extra={"gene": gene, "product": f"{kind} ribosomal RNA"})))
        if i < len(_ITS):
            parts.append((random_seq(sim.rng, int(sim.rng.integers(*_ITS[i])), 0.45), None))
    return parts


def real_trna(i: int) -> Parts:
    name, seq, aa3, ac = bundled()[1][i]
    extra = {"gene": name, "anti_codon": ac, "amino_acid": aa3, "product": f"tRNA-{aa3}"}
    return [(seq, Feat("tRNA", extra=extra))]


def trna_picks(rng: np.random.Generator, n: int) -> list[int]:
    if n <= 0:
        return []
    k = len(bundled()[1])
    return [int(x) for x in np.concatenate([rng.permutation(k) for _ in range(-(-n // k))])[:n]]


def crispr(sim: _Sim) -> Parts:
    rng = sim.rng
    n = int(rng.integers(5, 11))
    repeat = random_seq(rng, 29, 0.5)
    arr = repeat + "".join(
        random_seq(rng, int(rng.integers(32, 37))) + repeat for _ in range(n - 1)
    )
    product = f"CRISPR array with {n} repeats of length 29, consensus sequence {repeat}"
    return [(arr, Feat("crispr", extra={"product": product}))]


def oric(sim: _Sim) -> Parts:
    """A 232-bp oriC-like region (DnaA boxes, AT-rich 13-mers) and ``mnmG`` 25 bp downstream."""
    rng = sim.rng
    core = random_seq(rng, 232, 0.4)
    for k, motif in enumerate(["GATCTATTTATTT"] * 3 + ["TTATCCACA"] * 5):
        p = 10 + 25 * k
        core = core[:p] + motif + core[p + len(motif) :]
    cod = ["ATG", *sim.sense(628), sim.stop_codon()]
    return [
        (random_seq(rng, 40, 0.35), None),
        (core, Feat("oric", extra={"product": "oriC"})),
        (random_seq(rng, 25, 0.35), None),
        (
            "".join(cod),
            coding("cds", cod, table=sim.table, extra={"gene": "mnmG", "product": _MNMG}),
        ),
    ]


def gene_parts(sim: _Sim, length: int, gc: float) -> Parts:
    """Genes of ``sim``'s codon model separated by spacers (GC ``gc``), exactly ``length`` bp."""
    rng, used = sim.rng, 0
    parts: Parts = []
    lo, hi = _BLOCK_SPACER
    for _ in range(200):
        sp = int(rng.integers(lo, hi + 1))
        cod = sim.gene().codons
        if used + sp + 3 * len(cod) + lo > length:
            if parts:
                break
            continue
        parts += [
            (random_seq(rng, sp, gc), None),
            ("".join(cod), coding("cds", cod, table=sim.table)),
        ]
        used += sp + 3 * len(cod)
    parts.append((random_seq(rng, length - used, gc), None))
    return parts


def bgc(sim: _Sim) -> Parts:
    return gene_parts(sim, sim.cfg.bgc_bp, 1 - sim.cfg.intergenic_at)


def placeholder(n: int) -> Callable[[], Parts]:
    return lambda: [("N" * n, None)]


def v2_makers(sim: _Sim) -> list[tuple[str, Callable[[], Parts]]]:
    cfg = sim.cfg
    out: list[tuple[str, Callable[[], Parts]]] = []
    if cfg.crispr:
        out.append(("crispr", lambda: crispr(sim)))
    if cfg.oric:
        out.append(("oric", lambda: oric(sim)))
    if cfg.bgc_bp:
        out.append(("bgc", lambda: bgc(sim)))
    if cfg.prophage_bp:
        out.append(("prophage", placeholder(cfg.prophage_bp)))
    return out


def gc_of(s: bytes | bytearray | str) -> float:
    b = s.encode() if isinstance(s, str) else bytes(s)
    return (b.count(b"G") + b.count(b"C")) / max(1, len(b))


def tilted_cum(target: float, usage: Mapping[str, float] = SENSE_USAGE) -> np.ndarray:
    """Cumulative codon probabilities exp-tilted towards a per-base GC of ``target``."""
    w = np.array(list(usage.values()), dtype=float)
    g = np.array([sum(b in "GC" for b in c) / 3 for c in usage])
    lo, hi = -40.0, 40.0
    for _ in range(60):
        lam = (lo + hi) / 2
        p = w * np.exp(lam * (g - 0.5))
        lo, hi = (lam, hi) if (p @ g) / p.sum() < target else (lo, lam)
    return np.cumsum(p / p.sum())


def fill_prophage(
    sim: _Sim, seq: bytearray, pos: int, length: int, strand: str, target: float
) -> tuple[list[Feat], list[tuple[Feat, int]]]:
    """Fill ``seq[pos:pos+length]`` with genes at GC ``target`` (RBS planted at ``sim``'s
    strength); returns the placed features and the planted RBS spacers."""
    from gippasim.features import _block

    sim.sense_cum = tilted_cum(target, sim.usage)
    best: Any = None
    for _ in range(_GC_TRIES):
        trial = bytearray(seq)
        b = _block(gene_parts(sim, length, target), strand)
        trial[pos : pos + length] = b.seq.encode()
        feats = [f.placed("", pos) for f in b.feats]
        rbs = []
        for f in feats:
            sp = plant_gene_rbs(sim.rng, sim.cfg.rbs_strength, trial, f.start, f.end, f.strand)
            if sp is not None:
                rbs.append((f, sp))
        err = abs(gc_of(trial[pos : pos + length]) - target)
        if best is None or err < best[0]:
            best = (err, trial, feats, rbs)
        if err <= _GC_TOLERANCE:
            break
    seq[:] = best[1]
    return best[2], best[3]
