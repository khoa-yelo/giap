"""M0 simulator: one circular chromosome with protein-coding genes grouped in TUs (spec §11.3)."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from giapsim.codons import CODE_11, SENSE_USAGE, STOPS

_COMP = str.maketrans("ACGT", "TGCA")
_SENSE = np.array(list(SENSE_USAGE), dtype="<U3")
_MAX_OVERLAP_TRIES = 200
#: Overlaps (bp) for which a stop codon of the upstream gene can share bases with the
#: start codon (and following codons) of the downstream gene. 2, 5 and multiples of 3 cannot.
_FEASIBLE_OVERLAPS = frozenset(k for k in range(1, 21) if k % 3 and k not in (2, 5))


def _cum(weights: Any) -> np.ndarray:
    w = np.asarray(weights, dtype=float)
    return np.cumsum(w / w.sum())


def _draw(rng: np.random.Generator, cum: np.ndarray, size: int) -> np.ndarray:
    """Indices drawn from the categorical distribution with cumulative probabilities ``cum``."""
    return np.minimum(np.searchsorted(cum, rng.random(size), side="right"), len(cum) - 1)


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


@dataclass(frozen=True)
class SimConfig:
    """Simulator knobs (M0 subset of spec §11.3)."""

    genome_length: int = 50_000
    contig_id: str = "chr1"
    locus_prefix: str = "GSIM"
    mean_tu_size: float = 1.8  # geometric TU sizes
    p_same_strand: float = 0.6  # probability that a TU keeps the previous TU's strand
    gene_len_median_aa: float = 280.0  # log-normal gene lengths
    gene_len_sigma: float = 0.45
    gene_len_range_aa: tuple[int, int] = (50, 1500)
    start_codons: dict[str, float] = field(
        default_factory=lambda: {"ATG": 0.90, "GTG": 0.08, "TTG": 0.02}
    )
    stop_codons: dict[str, float] = field(
        default_factory=lambda: {"TAA": 0.61, "TGA": 0.30, "TAG": 0.09}
    )
    p_intra_tu_overlap: float = 0.3  # otherwise a gap drawn from intra_tu_gap_bp
    overlap_bp: dict[int, float] = field(
        default_factory=lambda: {1: 0.3, 4: 0.5, 8: 0.1, 11: 0.04, 14: 0.03, 17: 0.02, 20: 0.01}
    )
    intra_tu_gap_bp: tuple[int, int] = (0, 40)
    inter_tu_gap_bp: tuple[int, int] = (50, 300)
    intergenic_at: float = 0.65  # AT fraction of non-coding sequence

    def __post_init__(self) -> None:
        bad = set(self.overlap_bp) - _FEASIBLE_OVERLAPS
        if bad:
            raise ValueError(
                f"infeasible overlaps {sorted(bad)}; allowed {sorted(_FEASIBLE_OVERLAPS)}"
            )


@dataclass(frozen=True)
class SimResult:
    outdir: Path
    genome_fasta: Path
    bakta_json: Path
    genes_truth: Path
    n_genes: int


@dataclass
class _Gene:
    codons: list[str]
    tu: int = 0
    strand: str = "+"
    start: int = 0  # 0-based half-open genome coordinates
    end: int = 0

    @property
    def nt(self) -> str:
        return "".join(self.codons)


class _Sim:
    def __init__(self, cfg: SimConfig, rng: np.random.Generator) -> None:
        self.cfg = cfg
        self.rng = rng
        self.sense_cum = _cum(list(SENSE_USAGE.values()))
        self.starts = list(cfg.start_codons)
        self.start_cum = _cum(list(cfg.start_codons.values()))
        self.stops = list(cfg.stop_codons)
        self.stop_cum = _cum(list(cfg.stop_codons.values()))
        self.overlaps = list(cfg.overlap_bp)
        self.overlap_cum = _cum(list(cfg.overlap_bp.values()))

    def sense(self, n: int) -> list[str]:
        out: list[str] = _SENSE[_draw(self.rng, self.sense_cum, n)].tolist()
        return out

    def gene(self) -> _Gene:
        cfg, rng = self.cfg, self.rng
        lo, hi = cfg.gene_len_range_aa
        n_aa = int(
            np.clip(
                round(rng.lognormal(np.log(cfg.gene_len_median_aa), cfg.gene_len_sigma)), lo, hi
            )
        )
        start = self.starts[int(_draw(rng, self.start_cum, 1)[0])]
        stop = self.stops[int(_draw(rng, self.stop_cum, 1)[0])]
        return _Gene([start, *self.sense(n_aa - 1), stop])

    def overlap(self, a: _Gene, b: _Gene, k: int) -> bool:
        """Make the last ``k`` bp of ``a`` equal the first ``k`` bp of ``b`` (same strand).

        ``a`` keeps a stop codon and ``b`` its start codon; codons touched by the overlap are
        redrawn until neither gene has an internal stop. Returns False if no solution was found.
        """
        q = -(-k // 3)  # codons of either gene touched by the overlap
        free = a.codons[-q][: 3 * q - k]  # bases of a's codon that lie before the overlap
        tries = _MAX_OVERLAP_TRIES if q > 1 else 1
        for _ in range(tries):
            head = [b.codons[0], *self.sense(q - 1)]
            tail = free + "".join(head)[:k]
            a_tail = [tail[i : i + 3] for i in range(0, 3 * q, 3)]
            if a_tail[-1] in STOPS and not STOPS.intersection(a_tail[:-1]):
                a.codons[-q:] = a_tail
                b.codons[:q] = head
                return True
        return False

    def tu(self, strand: str, tu_id: int) -> tuple[list[_Gene], int]:
        """Genes of one TU with TU-local coordinates in transcription orientation; TU length."""
        cfg, rng = self.cfg, self.rng
        genes = [self.gene() for _ in range(int(rng.geometric(1.0 / cfg.mean_tu_size)))]
        cursor = 0
        for i, g in enumerate(genes):
            if i:
                k = 0
                if rng.random() < cfg.p_intra_tu_overlap:
                    k = self.overlaps[int(_draw(rng, self.overlap_cum, 1)[0])]
                    if not self.overlap(genes[i - 1], g, k):
                        k = 0
                cursor -= k
                if not k:
                    cursor += int(rng.integers(cfg.intra_tu_gap_bp[0], cfg.intra_tu_gap_bp[1] + 1))
            g.tu, g.strand, g.start = tu_id, strand, cursor
            cursor += 3 * len(g.codons)
            g.end = cursor
        return genes, cursor

    def layout(self) -> list[_Gene]:
        cfg, rng = self.cfg, self.rng
        lo, hi = cfg.inter_tu_gap_bp
        pos = int(rng.integers(lo, hi + 1))
        strand = "+" if rng.random() < 0.5 else "-"
        placed: list[_Gene] = []
        tu_id = 0
        while True:
            if tu_id and rng.random() >= cfg.p_same_strand:
                strand = "-" if strand == "+" else "+"
            genes, length = self.tu(strand, tu_id + 1)
            if pos + length + lo > cfg.genome_length:  # no gene crosses the origin
                return placed
            tu_id += 1
            for g in genes:
                g.start, g.end = (
                    (pos + g.start, pos + g.end)
                    if strand == "+"
                    else (pos + length - g.end, pos + length - g.start)
                )
            placed.extend(genes)
            pos += length + int(rng.integers(lo, hi + 1))

    def background(self) -> bytearray:
        at, gc = self.cfg.intergenic_at / 2, (1 - self.cfg.intergenic_at) / 2
        idx = _draw(self.rng, _cum([at, gc, gc, at]), self.cfg.genome_length)
        return bytearray(np.frombuffer(b"ACGT", dtype=np.uint8)[idx].tobytes())


def _translate(codons: list[str]) -> str:
    return "M" + "".join(CODE_11[c] for c in codons[1:-1])


def simulate(config: SimConfig, seed: int, outdir: Path | str) -> SimResult:
    """Simulate a genome with ``config`` and ``seed``; write inputs and truth under ``outdir``."""
    rng = np.random.default_rng(seed)
    sim = _Sim(config, rng)
    genes = sorted(sim.layout(), key=lambda g: (g.start, g.end))
    genome = sim.background()
    for g in genes:
        nt = g.nt if g.strand == "+" else revcomp(g.nt)
        genome[g.start : g.end] = nt.encode()
    sequence = genome.decode()

    out = Path(outdir)
    (out / "truth").mkdir(parents=True, exist_ok=True)
    loci = [f"{config.locus_prefix}_{5 * (i + 1):05d}" for i in range(len(genes))]

    fasta = out / "genome.fasta"
    lines = [sequence[i : i + 60] for i in range(0, len(sequence), 60)]
    fasta.write_text(f">{config.contig_id} [topology=circular]\n" + "\n".join(lines) + "\n")

    features = [
        {
            "type": "cds",
            "contig": config.contig_id,
            "start": g.start + 1,
            "stop": g.end,
            "strand": g.strand,
            "locus": locus,
            "start_type": g.codons[0],
            "nt": g.nt,
            "aa": _translate(g.codons),
        }
        for g, locus in zip(genes, loci, strict=True)
    ]
    data = {
        "sequences": [
            {
                "id": config.contig_id,
                "type": "chromosome",
                "topology": "circular",
                "complete": True,
                "length": len(sequence),
                "sequence": sequence,
            }
        ],
        "features": features,
    }
    bakta = out / "bakta_like.json"
    bakta.write_text(json.dumps(data, sort_keys=True, indent=1) + "\n")

    truth = out / "truth" / "genes.tsv"
    with truth.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["locus_tag", "tu_id", "start", "end", "strand", "start_codon", "length_aa"])
        for g, locus in zip(genes, loci, strict=True):
            w.writerow(
                [locus, f"TU_{g.tu}", g.start + 1, g.end, g.strand, g.codons[0], len(g.codons) - 1]
            )
    return SimResult(out, fasta, bakta, truth, len(genes))
