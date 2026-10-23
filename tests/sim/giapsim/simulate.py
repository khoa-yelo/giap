"""giapsim: genomes with genes in TUs, planted features, signals and precomputed layers (§11.3).

M0 skeleton, M3-18 features (plasmids, RNA genes, pseudogenes, sORFs, strains) and the M2-05
simulator v2 (signals, CRISPR/oriC/prophage/BGC blocks, precomputed layers and their truth). All
v2 fields default to off; the output is then byte-identical to M0/M3-18 (``truth/features.tsv``
gains its ``anticodon`` column only with ``write_precomputed`` or ``real_rrna_trna``).

Accepted choices (M2-05 review): plasmids are ``plasmid_<n>`` (tests also accept ``p<n>``);
with v2 signals, intra-TU gaps of 8-10 bp become 11-13 bp (no RBS spacer fits them); genes with
an N-terminal extension get no planted RBS; signals are planted at their strengths even when
``write_precomputed`` is off (no truth tables are then written).
"""

from __future__ import annotations

import csv
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

import numpy as np

from giap.io.gff3 import GffFeature, write_gff3
from giap.io.tsv import write_tsv
from giap.model.coords import Interval, Replicon
from giapsim import blocks as v2blocks
from giapsim import layers
from giapsim.codons import SENSE_USAGE, STOPS, revcomp
from giapsim.denovo import DenovoConfig, simulate_denovo
from giapsim.features import NOVEL_CLASSES, Block, Feat, coding, make_blocks, novel_orfs
from giapsim.proteomics import ProteomicsConfig, simulate_proteomics
from giapsim.scenario_truth import plant_novel_rbs, strong_tirs, write_m7_truth
from giapsim.signals import TuPlan, apply_plan, plan_tu, plant_gene_rbs
from giapsim.strain import StrainConfig, derive_reference

__all__ = [
    "SIM100K",
    "SIM100K_SEED",
    "SimConfig",
    "SimResult",
    "StrainConfig",
    "revcomp",
    "simulate",
]
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
    # weight tables are tuples of (key, weight) pairs so the config stays hashable; a mapping
    # passed in is converted
    start_codons: tuple[tuple[str, float], ...] = (("ATG", 0.90), ("GTG", 0.08), ("TTG", 0.02))
    stop_codons: tuple[tuple[str, float], ...] = (("TAA", 0.61), ("TGA", 0.30), ("TAG", 0.09))
    p_intra_tu_overlap: float = 0.3  # otherwise a gap drawn from intra_tu_gap_bp
    overlap_bp: tuple[tuple[int, float], ...] = (
        (1, 0.3), (4, 0.5), (8, 0.1), (11, 0.04), (14, 0.03), (17, 0.02), (20, 0.01)
    )  # fmt: skip
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
    n_novel_orfs: tuple[tuple[str, int], ...] = ()  # (class, count) pairs
    strain: StrainConfig | None = None
    # M2-05 (simulator v2); all off by default
    real_rrna_trna: bool = False  # bundled E. coli MG1655 rrnB / tRNA genes (data/PROVENANCE.md)
    crispr: bool = False
    oric: bool = False  # oriC-like region next to mnmG
    prophage_bp: int = 0
    bgc_bp: int = 0
    prophage_gc_shift: float = -0.08  # prophage GC minus the rest of the chromosome
    rbs_strength: float = 0.0  # per-gene probability of a planted AGGAGG
    promoter_strength: float = 0.0  # per-TU probability of a planted sigma-70 promoter
    terminator_strength: float = 0.0  # per-TU probability of a planted terminator
    layer_false_positive_rate: float = 0.05  # expected false rows per true row of a layer
    write_precomputed: bool = False
    # M3-19 simulated proteomics (peptides, database errors); None = off
    proteomics: ProteomicsConfig | None = None
    # M3-20 de novo reads and reference peptides (needs proteomics); None = off
    denovo: DenovoConfig | None = None
    # M7-12 scenarios (giapsim.scenario_truth); all off by default
    novel_orf_len_aa: tuple[int, int] = (10, 30)  # planted intergenic ORFs; min of found ones
    plant_novel_rbs: bool = False  # AGGAGG before planted intergenic / N-terminal true starts
    m7_truth: bool = False  # alt_sites, indel_sites and the extended novel_orfs truth

    @property
    def v2(self) -> bool:
        """Any simulator-v2 feature on (TUs then carry leaders/trailers for signals)."""
        return bool(
            self.real_rrna_trna or self.crispr or self.oric or self.prophage_bp or self.bgc_bp
            or self.rbs_strength or self.promoter_strength or self.terminator_strength
            or self.write_precomputed
        )  # fmt: skip

    def __post_init__(self) -> None:
        object.__setattr__(self, "novel_orf_len_aa", tuple(self.novel_orf_len_aa))
        lo, hi = self.novel_orf_len_aa
        if not 2 <= lo <= hi:
            raise ValueError("novel_orf_len_aa must satisfy 2 <= lo <= hi")
        for name in ("start_codons", "stop_codons", "overlap_bp", "n_novel_orfs"):
            v = getattr(self, name)
            object.__setattr__(
                self, name, tuple(v.items() if isinstance(v, Mapping) else map(tuple, v))
            )
        if self.denovo is not None and self.proteomics is None:
            raise ValueError("denovo needs SimConfig.proteomics")
        if set(dict(self.n_novel_orfs)) - set(NOVEL_CLASSES):
            raise ValueError(f"n_novel_orfs keys must be among {NOVEL_CLASSES}")
        bad = set(dict(self.overlap_bp)) - _FEASIBLE_OVERLAPS
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
    gff3: Path | None = None
    config_yaml: Path | None = None
    precomputed: Mapping[str, Path] = field(default_factory=dict)
    truth: Mapping[str, Path] = field(default_factory=dict)
    peptides: Path | None = None
    proteins_faa: Path | None = None
    abundance: Path | None = None


@dataclass
class _Gene:
    codons: list[str]
    tu: int = 0
    strand: str = "+"
    start: int = 0  # 0-based half-open genome coordinates
    end: int = 0
    feat: Feat | None = None

    @property
    def nt(self) -> str:
        return "".join(self.codons)


class _Sim:
    def __init__(
        self, cfg: SimConfig, rng: np.random.Generator, srng: np.random.Generator | None = None
    ) -> None:
        self.cfg = cfg
        self.rng = rng
        self.srng = srng  # v2 signals; TUs get leaders/trailers when set
        self.plans: list[TuPlan] = []
        self.sense_cum = _cum(list(SENSE_USAGE.values()))
        self.starts = [k for k, _ in cfg.start_codons]
        self.start_cum = _cum([w for _, w in cfg.start_codons])
        self.stops = [k for k, _ in cfg.stop_codons]
        self.stop_cum = _cum([w for _, w in cfg.stop_codons])
        self.overlaps = [k for k, _ in cfg.overlap_bp]
        self.overlap_cum = _cum([w for _, w in cfg.overlap_bp])

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
                    gap = int(rng.integers(cfg.intra_tu_gap_bp[0], cfg.intra_tu_gap_bp[1] + 1))
                    if self.srng is not None and 8 <= gap <= 10:
                        gap += 3  # no RBS spacer fits gaps of 8-10 bp (see signals.py)
                    cursor += gap
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
        max_genes: int | None = None,
    ) -> tuple[list[_Gene], list[tuple[int, Block]]]:
        """TUs from ``lead`` to ``genome_length - tail``; each block goes in once ``pos`` passes
        its target. Returns the genes and the block positions. With ``srng`` set, each TU's
        signals are planned (``self.plans``)."""
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
            if max_genes is not None:
                if len(placed) >= max_genes:
                    break
                genes = genes[: max_genes - len(placed)]
                length = genes[-1].end
            plan = None
            if self.srng is not None:
                plan = plan_tu(self.srng, cfg, genes, length)
                length = plan.length
            if pos + length + lo + reserve > limit:  # no gene crosses the origin
                break
            tu_id += 1
            if plan is not None:
                plan.pos = pos
                self.plans.append(plan)
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
    tu = f"TU_{g.tu}" if g.tu else ""
    return coding("cds", g.codons, start=g.start, end=g.end, strand=g.strand, tu=tu)


def _write_genes(genome: bytearray, genes: list[_Gene]) -> list[Feat]:
    for g in genes:
        nt = g.nt if g.strand == "+" else revcomp(g.nt)
        genome[g.start : g.end] = nt.encode()
        g.feat = _gene_feat(g)
    return [g.feat for g in genes if g.feat is not None]


def _extra_contig(
    cfg: SimConfig,
    rng: np.random.Generator,
    length: int,
    srng: np.random.Generator | None = None,
    max_genes: int | None = None,
) -> tuple[bytearray, list[Feat], list[TuPlan]]:
    sim = _Sim(replace(cfg, genome_length=length), rng, srng)
    genes, _ = sim.layout(max_genes=max_genes)
    genome = sim.background()
    return genome, _write_genes(genome, genes), sim.plans


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


_JSON_TYPE = {"oric": "oriC"}
_GFF_TYPE = {"crispr": "repeat_region", "oric": "origin_of_replication"}
_VERSION = {"bakta": "sim", "db": {"version": "sim", "type": "sim"}}


def _slice(seq: str, f: Feat) -> str:
    n = len(seq)
    s = seq[f.start : f.end] if f.end <= n else seq[f.start :] + seq[: f.end - n]
    return revcomp(s) if f.strand == "-" else s


def _json_feature(f: Feat, length: int, seq: str | None = None) -> dict[str, Any]:
    """One feature; ``seq`` (the contig) given selects Bakta 1.12 field names (simulator v2)."""
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
    if seq is None:
        return d
    d["sequence"] = d.pop("contig")
    d["type"] = _JSON_TYPE.get(f.kind, f.kind)
    if f.kind in ("cds", "sorf"):
        d = {"gene": None, "product": "hypothetical protein", "db_xrefs": []} | d
    else:
        d["nt"] = _slice(seq, f)
        if f.strand == ".":
            d["strand"] = "?"
    return d


def _write_json(
    path: Path,
    contigs: list[tuple[str, str, str]],
    seqs: dict[str, str],
    feats: list[Feat],
    bakta: bool = False,
) -> None:
    key = "nt" if bakta else "sequence"
    data: dict[str, Any] = {
        "sequences": [
            {
                "id": cid,
                "type": kind,
                "topology": topology,
                "complete": topology == "circular",
                "length": len(seqs[cid]),
                key: seqs[cid],
            }
            for cid, kind, topology in contigs
        ],
        "features": [
            _json_feature(f, len(seqs[f.contig]), seqs[f.contig] if bakta else None) for f in feats
        ],
    }
    if bakta:
        data["version"] = _VERSION
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
        gtype = _GFF_TYPE.get(gtype, gtype)
        attrs = {"ID": f.locus, "locus_tag": f.locus}
        if f.extra.get("product"):
            attrs["product"] = f.extra["product"]
        if f.kind == "sorf":  # Bakta writes sORFs as CDS
            attrs["sorf"] = "true"
        iv = Interval(f.contig, f.start, f.end, f.strand)  # type: ignore[arg-type]
        phase = 0 if gtype == "CDS" else None
        out.append(GffFeature(iv, "giapsim", gtype, phase=phase, attributes=attrs))
    write_gff3(path, out, reps)


@dataclass
class _V2:
    """Simulator-v2 state gathered before loci are assigned."""

    plans: dict[str, list[TuPlan]]
    rbs: list[tuple[Feat, int]] = field(default_factory=list)  # genes outside TUs
    blocks: list[tuple[str, int, int]] = field(default_factory=list)  # kind, start, end on chr


def _plant_signals(
    cfg: SimConfig,
    srng: np.random.Generator,
    seqs: dict[str, bytearray],
    feats: list[Feat],
    st: _V2,
    chrom: str,
    put: list[tuple[int, Block]],
) -> list[Feat]:
    """Apply TU signals, plant RBSs of block genes and fill the prophage; returns new feats."""
    for cid, plans in st.plans.items():
        for p in plans:
            apply_plan(seqs[cid], p)
    for f in feats:
        # genes with an N-terminal extension (``true_start`` set) get no RBS
        if f.kind in ("cds", "sorf") and not f.tu and f.true_start is None:
            sp = plant_gene_rbs(srng, cfg.rbs_strength, seqs[f.contig], f.start, f.end, f.strand)
            if sp is not None:
                st.rbs.append((f, sp))
    new: list[Feat] = []
    for pos, b in put:
        if b.kind in ("prophage", "bgc"):
            st.blocks.append((b.kind, pos, pos + len(b.seq)))
        if b.kind != "prophage":
            continue
        g = seqs[chrom]
        rest = g[:pos] + g[pos + len(b.seq) :]
        target = v2blocks.gc_of(rest) + cfg.prophage_gc_shift
        strand = "+" if srng.random() < 0.5 else "-"
        pfeats, rbs = v2blocks.fill_prophage(_Sim(cfg, srng), g, pos, len(b.seq), strand, target)
        for f in pfeats:
            f.contig = chrom
        st.rbs += rbs
        new += pfeats
    return new


def _signal_truth(st: _V2, chrom: str, feats: list[Feat], seqs: dict[str, str]) -> layers.Truth:
    t = layers.Truth()
    for cid, plans in st.plans.items():
        for p in plans:
            gf = [g.feat for g in p.genes if g.feat is not None]
            tss = None if p.tss is None else p.span(p.tss, p.tss + 1)[0] + 1
            term = None if p.term is None else p.span(*p.term)
            t.tus.append([
                gf[0].tu, cid, min(f.start for f in gf) + 1, max(f.end for f in gf), p.strand,
                [f.locus for f in gf], tss is not None, tss, term is not None,
            ])  # fmt: skip
            if tss is not None:
                t.promoters.append((cid, p.strand, tss))
            if term is not None:
                t.terminators.append((cid, term[0] + 1, term[1], p.strand))
            t.rbs += [(g.feat.locus, sp, g.feat.extra["start_type"]) for g, sp in p.rbs if g.feat]
    t.rbs += [(f.locus, sp, f.extra["start_type"]) for f, sp in st.rbs]
    for kind, s, e in st.blocks:
        t.regions.append((kind, chrom, s + 1, e))
    for f in feats:
        if f.kind == "crispr":
            t.regions.append(("crispr", f.contig, f.start + 1, f.end))
    for cid, n in ((c, len(s)) for c, s in seqs.items()):
        if cid.startswith("plasmid_"):
            t.regions.append(("plasmid", cid, 1, n))
    return t


def simulate(config: SimConfig, seed: int, outdir: Path | str) -> SimResult:
    """Simulate a genome with ``config`` and ``seed``; write inputs and truth under ``outdir``."""
    prot = config.proteomics
    if prot is not None and prot.wrong_source == "other_strain" and config.strain is None:
        raise ValueError("ProteomicsConfig(wrong_source='other_strain') needs SimConfig.strain")
    rng = np.random.default_rng(seed)
    frng = np.random.default_rng([seed, 1])  # M3-18 features; M0 defaults never draw from it
    srng = np.random.default_rng([seed, 2]) if config.v2 else None  # M2-05 signals and layers
    sim, fsim = _Sim(config, rng, srng), _Sim(config, frng)
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
        feats[-1].extra["edge"] = True  # Bakta's flag for a feature crossing the origin
    for f in feats:
        f.contig = chrom
    contigs = [(chrom, "chromosome", "circular")]
    seqs = {chrom: genome}
    st = _V2({chrom: sim.plans})
    extra = [
        ("plasmid", "circular", int(frng.integers(*config.plasmid_length, endpoint=True)))
        for _ in range(config.n_plasmids)
    ]
    if config.linear_contig_length:
        extra.append(("contig", "linear", config.linear_contig_length))
    for i, (kind, topology, n_bp) in enumerate(extra):
        cid = f"plasmid_{i + 1}" if kind == "plasmid" else "linear_1"  # tests also accept p1
        max_genes = int(srng.integers(1, 4)) if srng is not None and kind == "plasmid" else None
        seqs[cid], fs, st.plans[cid] = _extra_contig(config, frng, n_bp, srng, max_genes)
        contigs.append((cid, kind, topology))
        for f in fs:
            f.contig, f.tu = cid, f"{cid}:{f.tu}"
        feats += fs
    if srng is not None:
        feats += _plant_signals(config, srng, seqs, feats, st, chrom, put)
    order = {cid: k for k, (cid, _, _) in enumerate(contigs)}
    feats.sort(key=lambda f: (order[f.contig], f.start, f.end))
    for i, f in enumerate(feats):
        f.locus = f"{config.locus_prefix}_{5 * (i + 1):05d}"
    novel = (
        novel_orfs(fsim, dict(config.n_novel_orfs), seqs, feats, config.novel_orf_len_aa,
                   40 if config.plant_novel_rbs else 20)
        if config.n_novel_orfs
        else []
    )  # fmt: skip
    mrng = np.random.default_rng([seed, _M7_STREAM])  # M7-12 draws; unused by older configs
    planted = plant_novel_rbs(mrng, seqs, novel, feats) if config.plant_novel_rbs else set()
    text = {cid: s.decode() for cid, s in seqs.items()}

    out = Path(outdir)
    (out / "truth").mkdir(parents=True, exist_ok=True)
    fasta, bakta, gff = out / "genome.fasta", out / "bakta_like.json", out / "annotation.gff3"
    _write_fasta(fasta, contigs, text)
    _write_json(bakta, contigs, text, feats, bakta=config.v2)
    _write_gff3(gff, contigs, text, feats)
    novel_truth = out / "truth" / "novel_orfs.tsv"
    _tsv(novel_truth, ["orf_id", "class", "contig", "start", "end", "strand", "aa"], novel)
    truth, feat_truth = out / "truth" / "genes.tsv", out / "truth" / "features.tsv"
    # every CDS (pseudogenes included); ``end`` unwrapped (> contig length across the origin)
    _tsv(
        truth,
        ["locus_tag", "tu_id", "start", "end", "strand", "start_codon", "length_aa"],
        [
            [f.locus, f.tu, f.start + 1, f.end, f.strand, f.extra["start_type"], len(f.aa)]
            for f in feats
            if f.kind == "cds"
        ],
    )
    ref_pseudo: dict[str, bool] = {}
    lifted: list[Feat] | None = None
    strain_out: tuple[list[list[Any]], list[Feat], dict[str, str]] | None = None
    res = SimResult(
        out, fasta, bakta, truth, len(genes), gff, feat_truth, novel_truth, gff3=gff,
        truth={"genes": truth, "features": feat_truth},
    )  # fmt: skip
    if config.strain is not None:
        ref_seqs, lifted, rows = derive_reference(config.strain, seed, text, feats)
        ref_fasta, ref_json = out / "reference.fasta", out / "reference_bakta_like.json"
        strain_truth = out / "truth" / "strain_differences.tsv"
        _write_fasta(ref_fasta, contigs, ref_seqs)
        _write_json(ref_json, contigs, ref_seqs, lifted, bakta=config.v2)
        header = [
            "contig",
            "pos_sample",
            "pos_reference",
            "ref",
            "alt",
            "type",
            "gene",
            "aa_change",
        ]
        _tsv(strain_truth, header, rows)
        ref_pseudo = {g.locus: bool(g.extra.get("pseudo")) for g in lifted}
        strain_out = (rows, lifted, ref_seqs)
        res = replace(
            res, reference_fasta=ref_fasta, reference_json=ref_json, strain_truth=strain_truth
        )
    # ``anticodon`` only with v2 output, so M0/M3-18 truth stays byte-identical
    anticodon = config.write_precomputed or config.real_rrna_trna
    header = ["feature_id", "type", "contig", "start", "end", "strand", "true_start"]
    header += ["pseudo_reason", "ref_pseudo", *(["anticodon"] if anticodon else [])]
    _tsv(
        feat_truth,
        header,
        [
            [
                f.locus,
                f.kind,
                f.contig,
                f.start + 1,
                f.end,
                f.strand,
                f.true_start or "",
                f.pseudo_reason,
                "" if f.locus not in ref_pseudo else str(ref_pseudo[f.locus]).lower(),
                *([f.extra.get("anti_codon", "")] if anticodon else []),
            ]
            for f in feats
        ],
    )
    if srng is not None and config.write_precomputed:
        res = _write_v2(out, config, srng, st, chrom, feats, text, bakta, res)
        _write_rbs_novel(out, seed, text, novel, feats, bakta, res, planted, mrng)
    if prot is not None:
        topo_map = {cid: topo for cid, _, topo in contigs}
        paths = simulate_proteomics(prot, seed, out, feats, novel, text, topo_map, lifted)
        res = replace(
            res,
            peptides=paths["peptides"],
            proteins_faa=paths["proteins"],
            abundance=paths["abundance"],
            truth={**res.truth, "peptides": paths["truth"], "proteomics": paths["truth_meta"]},
        )
        if config.denovo is not None:
            dp = simulate_denovo(config.denovo, seed, out)
            truth = {"denovo": dp["denovo"], "denovo_precision": dp["denovo_precision"]}
            res = replace(res, truth={**res.truth, **truth})
    if config.m7_truth:
        m7 = write_m7_truth(out, novel, feats, text, planted, strain_out)
        res = replace(res, truth={**res.truth, **m7})
    return res


def _write_v2(
    out: Path,
    cfg: SimConfig,
    srng: np.random.Generator,
    st: _V2,
    chrom: str,
    feats: list[Feat],
    text: dict[str, str],
    bakta: Path,
    res: SimResult,
) -> SimResult:
    t = _signal_truth(st, chrom, feats, text)
    tus = out / "truth" / "tus.tsv"
    write_tsv(tus, layers.TUS_COLUMNS, t.tus)
    lengths = {c: len(s) for c, s in text.items()}
    coding = [f for f in feats if f.kind in ("cds", "sorf")]
    rows = layers.build(srng, cfg.layer_false_positive_rate, t, lengths, coding)
    pre, tru = layers.write(out, rows)
    cfg_yaml = out / "config.yaml"
    layers.write_config(cfg_yaml, bakta, pre)
    return replace(
        res, config_yaml=cfg_yaml, precomputed=pre, truth={**res.truth, "tus": tus, **tru}
    )


_RBS_NOVEL_STREAM = 5  # numpy SeedSequence key of the rbs_novel layer (M7-08)
_M7_STREAM = 6  # numpy SeedSequence key of the M7-12 draws (planted RBSs and their TIRs)


def _write_rbs_novel(
    out: Path,
    seed: int,
    text: dict[str, str],
    novel: Sequence[Sequence[Any]],
    feats: Sequence[Feat],
    bakta: Path,
    res: SimResult,
    planted: set[tuple[str, str, int]] | None = None,
    mrng: np.random.Generator | None = None,
) -> None:
    """``precomputed/rbs_novel.tsv`` from its own RNG stream (the other layers are unchanged)
    and its ``annotation.precomputed.rbs_novel`` key in ``config.yaml`` (M7-08 rule 6). Planted
    novel starts (M7-12) get TIRs from ``scenario_truth.PLANTED_TIR`` drawn with ``mrng``."""
    truth = {(str(r[2]), str(r[5]), int(r[3]) - 1 if r[5] == "+" else int(r[4]) - 1)
             for r in novel}  # fmt: skip
    truth |= {(f.contig, f.strand, f.true_start - 1) for f in feats
              if f.kind == "cds" and f.true_start}  # fmt: skip
    rows = layers.rbs_novel_rows(np.random.default_rng([seed, _RBS_NOVEL_STREAM]), text, truth)
    if planted and mrng is not None:
        strong_tirs(mrng, rows, planted)
    path = layers.write_rbs_novel(out, rows)
    assert res.config_yaml is not None
    layers.write_config(res.config_yaml, bakta, {**res.precomputed, "rbs_novel": path})


#: The 100-kb genome of M2-05 (recorded with real tools by M2-10).
SIM100K = SimConfig(
    genome_length=100_000,
    n_plasmids=1,
    plasmid_length=(8_000, 8_000),
    n_rrna_operons=1,
    n_trna=10,
    real_rrna_trna=True,
    crispr=True,
    oric=True,
    prophage_bp=12_000,
    bgc_bp=10_000,
    n_pseudogenes=2,
    n_sorfs=3,
    rbs_strength=0.8,
    promoter_strength=0.8,
    terminator_strength=0.8,
    write_precomputed=True,
)
SIM100K_SEED = 20261014
