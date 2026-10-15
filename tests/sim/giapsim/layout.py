"""Sequence-free genome layouts for M4 calibration (issue M4-07; GIAP_SPEC.md §11.3, §6.5).

``simulate_layout`` draws replicons (circular chromosome ``chr1``, circular plasmids ``p<i>``,
linear draft contigs ``ctg<i>``, in that order), genes in planted TUs (geometric sizes, strand kept
with ``p_same_strand``), log-normal protein lengths with sORFs and pseudogenes, back-translated
``nt``, UniOP-like pair probabilities, region blocks and xrefs. No nucleotide genome is built:
gene coordinates are synthetic but consistent (span 3·(len_aa + 1) bp; -20...+40 bp between genes
of a TU, 50-300 bp between TUs). ``layout_from_simulation`` turns the M0 full simulator output into
the same frames. Randomness only from the generator passed in; nothing is imported from ``giap``.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from giapsim.codons import CODE_11, SENSE_USAGE
from giapsim.contract_tables import REQUIRED_COLUMNS as CONTRACT_COLUMNS

__all__ = [
    "DEFAULT_REGIONS",
    "ECOLI_AA",
    "PROPHAGE_SUBSETS",
    "Layout",
    "LayoutConfig",
    "RegionSpec",
    "layout_from_simulation",
    "simulate_layout",
]

#: E. coli K-12 proteome amino-acid composition (percent; initiator M not included).
ECOLI_AA: Mapping[str, float] = {
    "A": 9.5, "R": 5.5, "N": 3.9, "D": 5.1, "C": 1.2, "Q": 4.4, "E": 5.8, "G": 7.4, "H": 2.3,
    "I": 6.0, "L": 10.7, "K": 4.4, "M": 2.8, "F": 3.9, "P": 4.4, "S": 5.8, "T": 5.4, "W": 1.5,
    "Y": 2.8, "V": 7.1,
}  # fmt: skip
#: Prophage function subsets (§7.1.7) and their member probabilities.
PROPHAGE_SUBSETS: tuple[tuple[str, float], ...] = (
    ("integrase", 0.05), ("repressor", 0.05), ("structural", 0.4), ("lysis", 0.1),
    ("cargo", 0.4),
)  # fmt: skip
_STARTS = (("ATG", 0.90), ("GTG", 0.08), ("TTG", 0.02))
_STOPS = (("TAA", 0.61), ("TGA", 0.30), ("TAG", 0.09))
_SORF_AA = (10, 29)  # sORF lengths (aa, inclusive): < 30 aa
_MIN_AA = 30  # CDS and pseudogene lengths are clipped below at 30 aa (never sORF-sized)
_SHORT_FACTOR = 0.3  # attack null 2: median length factor inside ``short_genes_in`` regions
_INTRA_GAP = (-20, 40)
_INTER_GAP = (50, 300)
_P_XREF = 0.7
_SOURCE_TOOL = "giapsim"


@dataclass(frozen=True)
class RegionSpec:
    """A planted region on the chromosome: a block of consecutive genes, or a scattered set."""

    region_type: str  # prophage, bgc, crispr, amr
    n_genes: int
    scattered: bool = False


DEFAULT_REGIONS: tuple[RegionSpec, ...] = (
    RegionSpec("prophage", 40),
    RegionSpec("bgc", 25),
    RegionSpec("crispr", 6),
    RegionSpec("amr", 10, scattered=True),
)


@dataclass(frozen=True)
class LayoutConfig:
    """Knobs of the sequence-free layout (spec M4-07 "Interface")."""

    n_genes_chromosome: int = 1000
    plasmids: tuple[int, ...] = ()  # genes per plasmid (circular)
    linear_contigs: tuple[int, ...] = ()  # genes per linear contig (draft assemblies)
    mean_tu_size: float = 1.8
    p_same_strand: float = 0.6
    gene_len_median_aa: float = 280.0
    gene_len_sigma: float = 0.45
    frac_sorf: float = 0.02
    frac_pseudogene: float = 0.01
    frac_hypothetical: float = 0.25
    regions: tuple[RegionSpec, ...] = DEFAULT_REGIONS
    uniop_beta_in: tuple[float, float] = (8, 2)  # p for true intra-TU pairs
    uniop_beta_out: tuple[float, float] = (2, 8)  # p for adjacent same-strand pairs across TUs
    # amino-acid composition as (residue, weight) pairs (a mapping passed in is converted)
    aa_composition: tuple[tuple[str, float], ...] = tuple(ECOLI_AA.items())
    short_genes_in: str | None = None  # attack null 2: region type whose genes are drawn short

    def __post_init__(self) -> None:
        for name in ("plasmids", "linear_contigs", "regions", "uniop_beta_in", "uniop_beta_out"):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        comp = self.aa_composition
        items = comp.items() if isinstance(comp, Mapping) else comp
        object.__setattr__(self, "aa_composition", tuple((str(a), float(w)) for a, w in items))


@dataclass(frozen=True)
class Layout:
    """Frames of a layout (internal representation: 0-based half-open, Python list cells).

    ``replicons``, ``genes``, ``tus``, ``regions``, ``region_genes`` and ``xrefs`` carry the
    contract columns (C1, C2, C5, C7, C8, C3); ``operon_pairs`` has ``gene_a`` (earlier in genome
    order), ``gene_b``, ``uniop_p``; ``truth`` has one row per gene (``locus_tag``, ``tu_id``,
    ``region_ids`` and the ``region_effect`` slot filled by the expression simulator).
    """

    replicons: pd.DataFrame
    genes: pd.DataFrame
    tus: pd.DataFrame
    operon_pairs: pd.DataFrame
    regions: pd.DataFrame
    region_genes: pd.DataFrame
    truth: pd.DataFrame
    xrefs: pd.DataFrame


# --------------------------------------------------------------------------- codons


def _cum(weights: Any) -> np.ndarray:
    w = np.asarray(weights, dtype=float)
    return np.cumsum(w / w.sum())


def _draw(rng: np.random.Generator, cum: np.ndarray, size: int) -> np.ndarray:
    return np.minimum(np.searchsorted(cum, rng.random(size), side="right"), len(cum) - 1)


def _codon_tables(residues: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Per residue: cumulative synonymous-codon weights (padded with 1) and codon bytes."""
    syn = [[c for c in SENSE_USAGE if CODE_11[c] == a] for a in residues]
    width = max(len(s) for s in syn)
    cum = np.ones((len(residues), width))
    codes = np.zeros((len(residues), width, 3), dtype=np.uint8)
    for i, s in enumerate(syn):
        if not s:
            raise ValueError(f"no sense codon for residue {residues[i]!r}")
        cum[i, : len(s)] = _cum([SENSE_USAGE[c] for c in s])
        for j, c in enumerate(s):
            codes[i, j] = np.frombuffer(c.encode(), dtype=np.uint8)
        codes[i, len(s) :] = codes[i, len(s) - 1]
    return cum, codes


def _codon_bytes(codons: tuple[tuple[str, float], ...]) -> tuple[np.ndarray, np.ndarray]:
    arr = np.array([np.frombuffer(c.encode(), dtype=np.uint8) for c, _ in codons])
    return _cum([w for _, w in codons]), arr


# --------------------------------------------------------------------------- validation


def _check(cfg: LayoutConfig) -> None:
    if cfg.n_genes_chromosome < 2:
        raise ValueError(f"n_genes_chromosome must be >= 2, got {cfg.n_genes_chromosome}")
    if any(n < 1 for n in (*cfg.plasmids, *cfg.linear_contigs)):
        raise ValueError("plasmids and linear contigs need at least one gene each")
    if cfg.mean_tu_size < 1:
        raise ValueError(f"mean_tu_size must be >= 1, got {cfg.mean_tu_size}")
    if not 0 <= cfg.p_same_strand <= 1:
        raise ValueError(f"p_same_strand must be in [0, 1], got {cfg.p_same_strand}")
    if cfg.frac_sorf < 0 or cfg.frac_pseudogene < 0 or cfg.frac_sorf + cfg.frac_pseudogene > 1:
        raise ValueError("frac_sorf and frac_pseudogene must be >= 0 with sum <= 1")
    if not 0 <= cfg.frac_hypothetical <= 1:
        raise ValueError(f"frac_hypothetical must be in [0, 1], got {cfg.frac_hypothetical}")
    n = cfg.n_genes_chromosome
    block = sum(r.n_genes for r in cfg.regions if not r.scattered)
    scattered = sum(r.n_genes for r in cfg.regions if r.scattered)
    if any(r.n_genes < 1 for r in cfg.regions):
        raise ValueError("every region needs at least one gene")
    if block > n or scattered > n - block:
        raise ValueError(
            f"regions need {block} block and {scattered} scattered genes; the chromosome has {n}"
        )


# --------------------------------------------------------------------------- pieces


def _tus_of_contig(
    rng: np.random.Generator, n: int, cfg: LayoutConfig
) -> tuple[list[int], list[str]]:
    """TU sizes (geometric, mean ``mean_tu_size``, last one truncated) and TU strands."""
    sizes_all = rng.geometric(1.0 / cfg.mean_tu_size, size=n)
    k = int(np.searchsorted(np.cumsum(sizes_all), n)) + 1
    sizes = sizes_all[:k].astype(int).tolist()
    sizes[-1] -= sum(sizes) - n
    first = int(rng.random() < 0.5)
    flips = (rng.random(k - 1) >= cfg.p_same_strand).astype(int)
    idx = (first + np.concatenate([[0], np.cumsum(flips)])) % 2
    return sizes, ["+" if i == 0 else "-" for i in idx.tolist()]


def _place_regions(
    rng: np.random.Generator, n: int, specs: tuple[RegionSpec, ...]
) -> list[tuple[RegionSpec, list[int]]]:
    """Non-overlapping blocks at random positions, then scattered sets outside the blocks."""
    blocks = [r for r in specs if not r.scattered]
    order = rng.permutation(len(blocks)).tolist()
    free = n - sum(r.n_genes for r in blocks)
    cuts = np.sort(rng.integers(0, free + 1, size=len(blocks))).tolist()
    placed: dict[int, list[int]] = {}
    used = 0
    for cut, b in zip(cuts, order, strict=True):
        start = cut + used
        placed[b] = list(range(start, start + blocks[b].n_genes))
        used += blocks[b].n_genes
    taken = {i for m in placed.values() for i in m}
    out: list[tuple[RegionSpec, list[int]]] = []
    bi = 0
    for r in specs:
        if r.scattered:
            pool = np.array([i for i in range(n) if i not in taken])
            members = np.sort(rng.choice(pool, size=r.n_genes, replace=False)).tolist()
            taken.update(members)
            out.append((r, members))
        else:
            out.append((r, placed[bi]))
            bi += 1
    return out


def _region_frames(
    rows: list[tuple[str, str, list[str], int | None, int | None]],
    contig_rank: dict[str, int],
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """C7 and C8 from (type, contig, members, start, end); IDs ``<type>_<n>`` by contig, start."""
    keyed = sorted(
        rows,
        key=lambda r: (contig_rank[r[1]], -1 if r[3] is None else r[3], r[0]),
    )
    counter: dict[str, int] = {}
    reg_rows, rg_rows = [], []
    subset_names = [s for s, _ in PROPHAGE_SUBSETS]
    subset_cum = _cum([p for _, p in PROPHAGE_SUBSETS])
    for typ, contig, members, start, end in keyed:
        counter[typ] = counter.get(typ, 0) + 1
        rid = f"{typ}_{counter[typ]}"
        reg_rows.append((rid, typ, contig, start, end, _SOURCE_TOOL, None))
        subsets: list[str | None] = [None] * len(members)
        core = [False] * len(members)
        if typ == "prophage":
            subsets = [subset_names[i] for i in _draw(rng, subset_cum, len(members)).tolist()]
        elif typ == "bgc":
            core = (rng.random(len(members)) < 0.3).tolist()
        for g, s, c in zip(members, subsets, core, strict=True):
            rg_rows.append((rid, g, 1.0, s, bool(c)))
    regions = pd.DataFrame(reg_rows, columns=list(CONTRACT_COLUMNS["regions"]))
    regions = regions.astype({"start": "Int64", "end": "Int64", "score": "float64"})
    region_genes = pd.DataFrame(rg_rows, columns=list(CONTRACT_COLUMNS["region_genes"]))
    region_genes = region_genes.astype({"overlap_frac": "float64", "core_biosynthetic": bool})
    return regions, region_genes


def _tus_frame(
    tags: list[str], contig: list[str], start: np.ndarray, end: np.ndarray,
    strand: list[str], tu_index: list[int],
) -> tuple[pd.DataFrame, list[str]]:  # fmt: skip
    """C5 from genes in genome order and their TU index (TUs consecutive, numbered in order)."""
    members: dict[int, list[int]] = {}
    for i, t in enumerate(tu_index):
        members.setdefault(t, []).append(i)
    rows = []
    gene_tu = [""] * len(tags)
    for k, t in enumerate(sorted(members, key=lambda t: members[t][0])):
        idx = members[t]
        tid = f"TU_{k + 1}"
        for i in idx:
            gene_tu[i] = tid
        s = strand[idx[0]]
        order = idx if s == "+" else idx[::-1]
        rows.append((tid, contig[idx[0]], int(start[idx].min()), int(end[idx].max()), s,
                     [tags[i] for i in order]))  # fmt: skip
    tus = pd.DataFrame(rows, columns=list(CONTRACT_COLUMNS["tus"]))
    return tus.astype({"start": "int64", "end": "int64"}), gene_tu


def _pairs_frame(
    tags: list[str], contig: list[str], strand: list[str], gene_tu: list[str],
    p_of: Any,
) -> pd.DataFrame:  # fmt: skip
    """Every adjacent same-strand pair in genome order (no wrap pair); ``p_of(same_tu)`` arrays."""
    a, b, same = [], [], []
    for i in range(len(tags) - 1):
        if contig[i] == contig[i + 1] and strand[i] == strand[i + 1]:
            a.append(tags[i])
            b.append(tags[i + 1])
            same.append(gene_tu[i] == gene_tu[i + 1])
    p = p_of(np.array(same, dtype=bool))
    return pd.DataFrame({"gene_a": a, "gene_b": b, "uniop_p": np.asarray(p, dtype=float)})


def _xrefs(rng: np.random.Generator, tags: list[str]) -> pd.DataFrame:
    """70 % of genes get ``COG####``; of those, about 60 % also GO terms and KEGG orthologs."""
    n = len(tags)
    has = rng.random(n) < _P_XREF
    cog = rng.integers(1, 5000, size=n)
    n_go = np.where(rng.random(n) < 0.6, rng.integers(1, 4, size=n), 0)
    go = rng.integers(1, 100_000, size=(n, 3))
    kegg = np.where(rng.random(n) < 0.6, rng.integers(1, 30_000, size=n), 0)
    rows = []
    for i in np.flatnonzero(has).tolist():
        rows.append((tags[i], "COG", f"COG{cog[i]:04d}"))
        rows += [(tags[i], "GO", f"GO:{int(t):07d}") for t in sorted(set(go[i, : n_go[i]]))]
        if kegg[i]:
            rows.append((tags[i], "KEGG", f"K{int(kegg[i]):05d}"))
    return pd.DataFrame(rows, columns=list(CONTRACT_COLUMNS["xrefs"]))


def _truth(tags: list[str], gene_tu: list[str], region_genes: pd.DataFrame) -> pd.DataFrame:
    of: dict[str, list[str]] = {}
    for rid, g in zip(region_genes["region_id"], region_genes["locus_tag"], strict=True):
        of.setdefault(g, []).append(rid)
    return pd.DataFrame(
        {
            "locus_tag": tags,
            "tu_id": gene_tu,
            "region_ids": [sorted(of.get(g, [])) for g in tags],
            "region_effect": np.zeros(len(tags)),
        }
    )


# --------------------------------------------------------------------------- simulate_layout


def simulate_layout(cfg: LayoutConfig, rng: np.random.Generator) -> Layout:
    """A sequence-free layout for ``cfg`` drawn from ``rng`` (pure function of both)."""
    _check(cfg)
    reps = [("chr1", cfg.n_genes_chromosome, "circular", "chromosome")]
    reps += [(f"p{i + 1}", n, "circular", "plasmid") for i, n in enumerate(cfg.plasmids)]
    reps += [(f"ctg{i + 1}", n, "linear", "unknown") for i, n in enumerate(cfg.linear_contigs)]
    contig_rank = {c: i for i, (c, *_) in enumerate(reps)}

    # TUs and strands (genes in genome order, contig by contig)
    contig: list[str] = []
    strand: list[str] = []
    tu_index: list[int] = []
    t = 0
    for cid, n, _, _ in reps:
        sizes, strands = _tus_of_contig(rng, n, cfg)
        for size, s in zip(sizes, strands, strict=True):
            contig += [cid] * size
            strand += [s] * size
            tu_index += [t] * size
            t += 1
    total = len(contig)

    # regions on the chromosome (gene indices 0..n-1 are the chromosome genes)
    placed = _place_regions(rng, cfg.n_genes_chromosome, cfg.regions)

    # gene classes and lengths
    u = rng.random(total)
    is_pseudo = u < cfg.frac_sorf + cfg.frac_pseudogene
    ftype = np.where(u < cfg.frac_sorf, "sorf", np.where(is_pseudo, "pseudogene", "cds"))
    hypothetical = rng.random(total) < cfg.frac_hypothetical
    median = np.full(total, float(cfg.gene_len_median_aa))
    if cfg.short_genes_in is not None:
        for spec, members in placed:
            if spec.region_type == cfg.short_genes_in:
                median[members] *= _SHORT_FACTOR
    draws = rng.lognormal(np.log(median), cfg.gene_len_sigma)
    len_aa = np.maximum(np.rint(draws).astype(np.int64), _MIN_AA)
    sorf_len = rng.integers(_SORF_AA[0], _SORF_AA[1] + 1, size=total)
    len_aa = np.where(ftype == "sorf", sorf_len, len_aa)

    # sequences: residues after the initiator, back-translated with the bundled codon usage
    residues = [a for a, _ in cfg.aa_composition]
    res_cum = _cum([w for _, w in cfg.aa_composition])
    syn_cum, syn_codes = _codon_tables(residues)
    body_n = len_aa - 1
    n_body = int(body_n.sum())
    res_idx = _draw(rng, res_cum, n_body)
    pick = (rng.random(n_body)[:, None] > syn_cum[res_idx]).sum(axis=1)
    body = syn_codes[res_idx, np.minimum(pick, syn_codes.shape[1] - 1)]  # (n_body, 3)
    letters = np.frombuffer("".join(residues).encode(), dtype=np.uint8)[res_idx]
    start_cum, start_codes = _codon_bytes(_STARTS)
    stop_cum, stop_codes = _codon_bytes(_STOPS)
    first = _draw(rng, start_cum, total)
    last = _draw(rng, stop_cum, total)
    pseudo = np.flatnonzero(ftype == "pseudogene")
    internal_stop = _draw(rng, stop_cum, len(pseudo))
    stop_at = (rng.random(len(pseudo)) * (body_n[pseudo])).astype(np.int64)  # body codon index
    aa_na = rng.random(len(pseudo)) < 0.5
    offsets = np.concatenate([[0], np.cumsum(body_n)])
    body = body.copy()
    for k, g in enumerate(pseudo.tolist()):
        body[offsets[g] + stop_at[k]] = stop_codes[internal_stop[k]]
    body_bytes = body.reshape(-1).tobytes()
    letter_bytes = letters.tobytes()
    nts: list[str] = []
    aas: list[str | None] = []
    start_codon: list[str] = []
    na_of = dict(zip(pseudo.tolist(), aa_na.tolist(), strict=True))
    for g in range(total):
        sc = start_codes[first[g]].tobytes().decode()
        o0, o1 = int(offsets[g]), int(offsets[g + 1])
        nts.append(
            sc + body_bytes[3 * o0 : 3 * o1].decode() + stop_codes[last[g]].tobytes().decode()
        )
        aas.append(None if na_of.get(g, False) else "M" + letter_bytes[o0:o1].decode())
        start_codon.append(sc)

    # coordinates
    span = 3 * (len_aa + 1)
    gap_intra = rng.integers(_INTRA_GAP[0], _INTRA_GAP[1] + 1, size=total)
    gap_inter = rng.integers(_INTER_GAP[0], _INTER_GAP[1] + 1, size=total + len(reps))
    start = np.zeros(total, dtype=np.int64)
    end = np.zeros(total, dtype=np.int64)
    lengths: list[int] = []
    partial = np.zeros(total, dtype=bool)
    g = 0
    for r, (_cid, n, topology, _kind) in enumerate(reps):
        cursor = int(gap_inter[total + r])
        for j in range(g, g + n):
            if j > g:
                same = tu_index[j] == tu_index[j - 1]
                cursor = int(end[j - 1]) + int(gap_intra[j] if same else gap_inter[j])
            start[j] = cursor
            end[j] = cursor + span[j]
        lengths.append(int(end[g + n - 1]) + int(gap_inter[g]))
        if topology == "linear":
            partial[g] = partial[g + n - 1] = True
        g += n

    tags = [f"GSIM_{5 * (i + 1):05d}" for i in range(total)]
    product = ["hypothetical protein" if h else f"simulated protein {i + 1}"
               for i, h in enumerate(hypothetical.tolist())]  # fmt: skip
    genes = pd.DataFrame(
        {
            "locus_tag": tags,
            "gene": pd.Series([None] * total, dtype=object),
            "product": product,
            "feature_type": ftype.tolist(),
            "contig": contig,
            "start": start,
            "end": end,
            "strand": strand,
            "partial": partial,
            "crosses_origin": np.zeros(total, dtype=bool),
            "hypothetical": hypothetical,
            "aa": pd.Series(aas, dtype=object),
            "nt": nts,
            "start_codon": start_codon,
        },
        columns=list(CONTRACT_COLUMNS["genes"]),
    )
    replicons = pd.DataFrame(
        [(cid, lengths[i], topo, kind) for i, (cid, _, topo, kind) in enumerate(reps)],
        columns=list(CONTRACT_COLUMNS["replicons"]),
    ).astype({"length": "int64"})

    tus, gene_tu = _tus_frame(tags, contig, start, end, strand, tu_index)
    bi, bo = cfg.uniop_beta_in, cfg.uniop_beta_out

    def p_of(same: np.ndarray) -> np.ndarray:
        p_in = rng.beta(bi[0], bi[1], size=len(same))
        p_out = rng.beta(bo[0], bo[1], size=len(same))
        return np.where(same, p_in, p_out)

    pairs = _pairs_frame(tags, contig, strand, gene_tu, p_of)

    region_rows: list[tuple[str, str, list[str], int | None, int | None]] = []
    for spec, members in placed:
        names = [tags[i] for i in members]
        if spec.scattered:
            region_rows.append((spec.region_type, "chr1", names, None, None))
        else:
            region_rows.append((spec.region_type, "chr1", names, int(start[members].min()),
                                int(end[members].max())))  # fmt: skip
    offset = 0
    for cid, n, _, kind in reps:
        if kind == "plasmid":
            region_rows.append(("plasmid", cid, tags[offset : offset + n], 0,
                                lengths[contig_rank[cid]]))  # fmt: skip
        offset += n
    regions, region_genes = _region_frames(region_rows, contig_rank, rng)
    xrefs = _xrefs(rng, tags)
    return Layout(replicons, genes, tus, pairs, regions, region_genes,
                  _truth(tags, gene_tu, region_genes), xrefs)  # fmt: skip


# --------------------------------------------------------------------------- M0 adapter


def _empty_regions() -> tuple[pd.DataFrame, pd.DataFrame]:
    regions = pd.DataFrame({c: pd.Series([], dtype=object) for c in CONTRACT_COLUMNS["regions"]})
    rg = pd.DataFrame({c: pd.Series([], dtype=object) for c in CONTRACT_COLUMNS["region_genes"]})
    return regions, rg


def _xref_rows(tag: str, xrefs: list[str]) -> list[tuple[str, str, str]]:
    rows = []
    for x in xrefs:
        db, _, ident = str(x).partition(":")
        if ident:
            rows.append((tag, db, x if db == "GO" else ident))
    return rows


def layout_from_simulation(sim: Any) -> Layout:
    """Layout frames from an M0 ``SimResult`` (``bakta_like.json`` and ``truth/genes.tsv``).

    Genes are the CDS and sORF features (pseudogenes flagged ``pseudo``) with 0-based coordinates
    (``end`` unwrapped across the origin); TUs follow ``tu_id`` of the truth, genes without one are
    singleton TUs, and all TUs are renumbered ``TU_<n>`` by contig order then start. Pair
    probabilities are the truth (1 within a TU, 0 across); regions are empty until M2 extends the
    full simulator. Reads only; the simulator output is not modified.
    """
    data = json.loads(Path(sim.bakta_json).read_text(encoding="utf-8"))
    with Path(sim.genes_truth).open(newline="", encoding="utf-8") as fh:
        truth_tu = {r["locus_tag"]: r["tu_id"] for r in csv.DictReader(fh, delimiter="\t")}
    kinds = {"chromosome": "chromosome", "plasmid": "plasmid"}
    reps = [
        (s["id"], int(s["length"]), s["topology"], kinds.get(s.get("type", ""), "unknown"))
        for s in data["sequences"]
    ]
    length = {c: n for c, n, _, _ in reps}
    rank = {c: i for i, (c, *_) in enumerate(reps)}
    rows = []
    xref_rows: list[tuple[str, str, str]] = []
    for f in data["features"]:
        if f.get("type") not in ("cds", "sorf"):
            continue
        cid = f.get("contig", f.get("sequence"))
        s0, stop = int(f["start"]) - 1, int(f["stop"])
        cross = stop <= s0
        e = stop + length[cid] if cross else stop
        pseudo = bool(f.get("pseudo") or f.get("pseudogene"))
        ftype = "pseudogene" if pseudo else str(f["type"])
        aa = f.get("aa")
        aa = aa.rstrip("*") if isinstance(aa, str) and aa else None
        product = f.get("product")
        nt = str(f["nt"])
        rows.append((
            f["locus"], f.get("gene"), product, ftype, cid, s0, e, f["strand"],
            bool(f.get("partial", False)), cross,
            product is None or product == "hypothetical protein", aa, nt,
            f.get("start_type") or nt[:3],
        ))  # fmt: skip
        xref_rows += _xref_rows(f["locus"], f.get("db_xrefs") or [])
    genes = pd.DataFrame(rows, columns=list(CONTRACT_COLUMNS["genes"]))
    genes = (
        genes.assign(_r=genes["contig"].map(rank))
        .sort_values(["_r", "start", "end", "strand", "locus_tag"], kind="mergesort")
        .drop(columns="_r")
        .reset_index(drop=True)
        .astype({"start": "int64", "end": "int64"})
    )
    tags = genes["locus_tag"].tolist()
    contig = genes["contig"].tolist()
    strand = genes["strand"].tolist()
    key: dict[str, int] = {}
    tu_index = []
    for g, c in zip(tags, contig, strict=True):
        t = truth_tu.get(g) or f"\0{g}"  # genes without a TU: singletons
        tu_index.append(key.setdefault(f"{c}\0{t}", len(key)))
    start, end = genes["start"].to_numpy(), genes["end"].to_numpy()
    tus, gene_tu = _tus_frame(tags, contig, start, end, strand, tu_index)
    pairs = _pairs_frame(tags, contig, strand, gene_tu, lambda same: same.astype(float))
    replicons = pd.DataFrame(reps, columns=list(CONTRACT_COLUMNS["replicons"]))
    regions, region_genes = _empty_regions()
    xrefs = pd.DataFrame(xref_rows, columns=list(CONTRACT_COLUMNS["xrefs"]))
    return Layout(replicons, genes, tus, pairs, regions, region_genes,
                  _truth(tags, gene_tu, region_genes), xrefs)  # fmt: skip
