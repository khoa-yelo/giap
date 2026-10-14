"""M0 simulator: one circular chromosome with protein-coding genes grouped in TUs (spec §11.3)."""

from __future__ import annotations

import csv
import json
from collections.abc import Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

import numpy as np

from giap.io.gff3 import GffFeature, write_gff3
from giap.model.coords import Interval, Replicon
from giapsim.codons import SENSE_USAGE, STOPS, revcomp
from giapsim.features import NOVEL_CLASSES, Block, Feat, coding, make_blocks, novel_orfs
from giapsim.strain import StrainConfig, derive_reference

__all__ = ["SimConfig", "SimResult", "StrainConfig", "revcomp", "simulate"]
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
    # M3-18 features; all off by default (output then byte-identical to M0)
    n_plasmids: int = 0
    plasmid_length: tuple[int, int] = (5_000, 60_000)
    linear_contig_length: int = 0  # > 0 adds one linear contig (no wrap)
    origin_crossing_gene: bool = False  # place one gene across the chromosome origin
    n_rrna_operons: int = 0
    n_trna: int = 0
    n_pseudogenes: int = 0
    n_sorfs: int = 0
    n_nterm_extensions: int = 0  # annotated start moved 3·k nt downstream (k in [5, 30])
    n_novel_orfs: dict[str, int] = field(default_factory=dict)
    strain: StrainConfig | None = None

    def __post_init__(self) -> None:
        if set(self.n_novel_orfs) - set(NOVEL_CLASSES):
            raise ValueError(f"n_novel_orfs keys must be among {NOVEL_CLASSES}")
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
    annotation_gff3: Path | None = None
    features_truth: Path | None = None
    novel_orfs_truth: Path | None = None
    reference_fasta: Path | None = None
    reference_json: Path | None = None
    strain_truth: Path | None = None


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
        start, stop = self.start_codon(), self.stop_codon()
        return _Gene([start, *self.sense(n_aa - 1), stop])

    def start_codon(self) -> str:
        return self.starts[int(_draw(self.rng, self.start_cum, 1)[0])]

    def stop_codon(self) -> str:
        return self.stops[int(_draw(self.rng, self.stop_cum, 1)[0])]

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

    def layout(
        self,
        blocks: Sequence[Block] = (),
        targets: Sequence[int] = (),
        lead: int = 0,
        tail: int = 0,
    ) -> tuple[list[_Gene], list[tuple[int, Block]]]:
        """TUs from ``lead`` to ``genome_length - tail``; each block goes in once ``pos`` passes
        its target. Returns the genes and the block positions."""
        cfg, rng = self.cfg, self.rng
        lo, hi = cfg.inter_tu_gap_bp
        limit = cfg.genome_length - tail
        pending = list(zip(targets, blocks, strict=True))
        reserve = sum(len(b.seq) + hi for b in blocks) + (hi if blocks else 0)
        put: list[tuple[int, Block]] = []
        pos = lead + int(rng.integers(lo, hi + 1))

        def place(b: Block) -> None:
            nonlocal pos, reserve
            if pos + len(b.seq) + lo > limit:
                raise ValueError(f"genome_length {cfg.genome_length} too short for the features")
            put.append((pos, b))
            pos += len(b.seq) + int(rng.integers(lo, hi + 1))
            reserve -= len(b.seq) + hi

        strand = "+" if rng.random() < 0.5 else "-"
        placed: list[_Gene] = []
        tu_id = 0
        while True:
            if tu_id and rng.random() >= cfg.p_same_strand:
                strand = "-" if strand == "+" else "+"
            while pending and pos >= pending[0][0]:
                place(pending.pop(0)[1])
            genes, length = self.tu(strand, tu_id + 1)
            if pos + length + lo + reserve > limit:  # no gene crosses the origin
                break
            tu_id += 1
            for g in genes:
                g.start, g.end = (
                    (pos + g.start, pos + g.end)
                    if strand == "+"
                    else (pos + length - g.end, pos + length - g.start)
                )
            placed.extend(genes)
            pos += length + int(rng.integers(lo, hi + 1))
        for _, b in pending:
            place(b)
        return placed, put

    def background(self) -> bytearray:
        at, gc = self.cfg.intergenic_at / 2, (1 - self.cfg.intergenic_at) / 2
        idx = _draw(self.rng, _cum([at, gc, gc, at]), self.cfg.genome_length)
        return bytearray(np.frombuffer(b"ACGT", dtype=np.uint8)[idx].tobytes())


def _gene_feat(g: _Gene) -> Feat:
    return coding("cds", g.codons, start=g.start, end=g.end, strand=g.strand)


def _write_genes(genome: bytearray, genes: list[_Gene]) -> list[Feat]:
    for g in genes:
        nt = g.nt if g.strand == "+" else revcomp(g.nt)
        genome[g.start : g.end] = nt.encode()
    return [_gene_feat(g) for g in genes]


def _extra_contig(
    cfg: SimConfig, rng: np.random.Generator, length: int
) -> tuple[bytearray, list[Feat]]:
    sim = _Sim(replace(cfg, genome_length=length), rng)
    genes, _ = sim.layout()
    genome = sim.background()
    return genome, _write_genes(genome, genes)


def _tsv(path: Path, header: list[str], rows: Sequence[Sequence[Any]]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def _write_fasta(path: Path, contigs: list[tuple[str, str, str]], seqs: dict[str, str]) -> None:
    parts = []
    for cid, _, topology in contigs:
        s = seqs[cid]
        lines = [s[i : i + 60] for i in range(0, len(s), 60)]
        parts.append(f">{cid} [topology={topology}]\n" + "\n".join(lines) + "\n")
    path.write_text("".join(parts))


def _json_feature(f: Feat, length: int) -> dict[str, Any]:
    d: dict[str, Any] = {
        "type": f.kind,
        "contig": f.contig,
        "start": f.start + 1,
        "stop": f.end if f.end <= length else f.end - length,  # Bakta-style start > stop
        "strand": f.strand,
        "locus": f.locus,
        **f.extra,
    }
    if f.kind in ("cds", "sorf"):
        d |= {"nt": f.nt, "aa": f.aa}
    else:
        d["id"] = f.locus
    return d


def _write_json(
    path: Path, contigs: list[tuple[str, str, str]], seqs: dict[str, str], feats: list[Feat]
) -> None:
    data = {
        "sequences": [
            {
                "id": cid,
                "type": kind,
                "topology": topology,
                "complete": topology == "circular",
                "length": len(seqs[cid]),
                "sequence": seqs[cid],
            }
            for cid, kind, topology in contigs
        ],
        "features": [_json_feature(f, len(seqs[f.contig])) for f in feats],
    }
    path.write_text(json.dumps(data, sort_keys=True, indent=1) + "\n")


def _write_gff3(
    path: Path, contigs: list[tuple[str, str, str]], seqs: dict[str, str], feats: list[Feat]
) -> None:
    reps = [
        Replicon(cid, len(seqs[cid]), topology, kind if kind != "contig" else "unknown")  # type: ignore[arg-type]
        for cid, kind, topology in contigs
    ]
    out = []
    for f in feats:
        gtype = "pseudogene" if f.pseudo_reason else "CDS" if f.kind in ("cds", "sorf") else f.kind
        attrs = {"ID": f.locus, "locus_tag": f.locus}
        if f.extra.get("product"):
            attrs["product"] = f.extra["product"]
        iv = Interval(f.contig, f.start, f.end, f.strand)  # type: ignore[arg-type]
        phase = 0 if gtype == "CDS" else None
        out.append(GffFeature(iv, "giapsim", gtype, phase=phase, attributes=attrs))
    write_gff3(path, out, reps)


def simulate(config: SimConfig, seed: int, outdir: Path | str) -> SimResult:
    """Simulate a genome with ``config`` and ``seed``; write inputs and truth under ``outdir``."""
    rng = np.random.default_rng(seed)
    frng = np.random.default_rng([seed, 1])  # M3-18 features; M0 defaults never draw from it
    sim, fsim = _Sim(config, rng), _Sim(config, frng)
    length, chrom = config.genome_length, config.contig_id
    origin, lead, tail = None, 0, 0
    if config.origin_crossing_gene:
        origin = fsim.gene()
        n = 3 * len(origin.codons)
        tail = int(frng.integers(n // 4, 3 * n // 4 + 1))
        lead = n - tail
        origin.strand = "+" if frng.random() < 0.5 else "-"
        origin.start, origin.end = length - tail, length - tail + n
    blocks = make_blocks(fsim)
    room = length - tail - sum(len(b.seq) + 300 for b in blocks)
    targets = sorted(frng.integers(lead, max(lead + 1, room), size=len(blocks)).tolist())
    genes, put = sim.layout(blocks, targets, lead, tail)
    genes.sort(key=lambda g: (g.start, g.end))
    genome = sim.background()
    gene_feats = _write_genes(genome, genes)
    feats = list(gene_feats)
    for pos, b in put:
        genome[pos : pos + len(b.seq)] = b.seq.encode()
        feats += [f.placed(chrom, pos) for f in b.feats]
    if origin is not None:
        nt = (origin.nt if origin.strand == "+" else revcomp(origin.nt)).encode()
        genome[origin.start :], genome[:lead] = nt[:tail], nt[tail:]
        feats.append(_gene_feat(origin))
    for f in feats:
        f.contig = chrom
    contigs = [(chrom, "chromosome", "circular")]
    seqs = {chrom: genome}
    extra = [
        ("plasmid", "circular", int(frng.integers(*config.plasmid_length, endpoint=True)))
        for _ in range(config.n_plasmids)
    ]
    if config.linear_contig_length:
        extra.append(("contig", "linear", config.linear_contig_length))
    for i, (kind, topology, n_bp) in enumerate(extra):
        cid = f"plasmid_{i + 1}" if kind == "plasmid" else "linear_1"
        seqs[cid], fs = _extra_contig(config, frng, n_bp)
        contigs.append((cid, kind, topology))
        feats += [replace(f, contig=cid) for f in fs]
    order = {cid: k for k, (cid, _, _) in enumerate(contigs)}
    feats.sort(key=lambda f: (order[f.contig], f.start, f.end))
    for i, f in enumerate(feats):
        f.locus = f"{config.locus_prefix}_{5 * (i + 1):05d}"
    novel = novel_orfs(fsim, config.n_novel_orfs, seqs, feats) if config.n_novel_orfs else []
    text = {cid: s.decode() for cid, s in seqs.items()}

    out = Path(outdir)
    (out / "truth").mkdir(parents=True, exist_ok=True)
    fasta, bakta, gff = out / "genome.fasta", out / "bakta_like.json", out / "annotation.gff3"
    _write_fasta(fasta, contigs, text)
    _write_json(bakta, contigs, text, feats)
    _write_gff3(gff, contigs, text, feats)
    truth = out / "truth" / "genes.tsv"
    _tsv(
        truth,
        ["locus_tag", "tu_id", "start", "end", "strand", "start_codon", "length_aa"],
        [
            [f.locus, f"TU_{g.tu}", g.start + 1, g.end, g.strand, g.codons[0], len(g.codons) - 1]
            for g, f in zip(genes, gene_feats, strict=True)
        ],
    )
    feat_truth = out / "truth" / "features.tsv"
    _tsv(
        feat_truth,
        ["feature_id", "type", "contig", "start", "end", "strand", "true_start", "pseudo_reason"],
        [
            [
                f.locus,
                f.kind,
                f.contig,
                d["start"],
                d["stop"],
                f.strand,
                f.true_start or "",
                f.pseudo_reason,
            ]
            for f in feats
            for d in [_json_feature(f, len(text[f.contig]))]
        ],
    )
    novel_truth = out / "truth" / "novel_orfs.tsv"
    _tsv(novel_truth, ["orf_id", "class", "contig", "start", "end", "strand", "aa"], novel)
    res = SimResult(out, fasta, bakta, truth, len(genes), gff, feat_truth, novel_truth)
    if config.strain is None:
        return res
    ref_seqs, lifted, rows = derive_reference(config.strain, seed, text, feats)
    ref_fasta, ref_json = out / "reference.fasta", out / "reference_bakta_like.json"
    strain_truth = out / "truth" / "strain_differences.tsv"
    _write_fasta(ref_fasta, contigs, ref_seqs)
    _write_json(ref_json, contigs, ref_seqs, lifted)
    header = ["contig", "pos_sample", "pos_reference", "ref", "alt", "type", "gene", "aa_change"]
    _tsv(strain_truth, header, rows)
    return replace(
        res, reference_fasta=ref_fasta, reference_json=ref_json, strain_truth=strain_truth
    )
