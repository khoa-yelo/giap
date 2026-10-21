"""Independent checks for giapsim's sequence-free layout (issue M4-07; GIAP_SPEC.md §11.3).

Written from docs/specs/M4-07-sim-layout-tables.md and M2-20's column definitions, without reading
the simulator: a table-11 translator, the expected ``uniop_p_prev/next`` in transcription order,
gap and TU structure checks. Frames are the internal representation (0-based half-open, Python
lists). Every ``check_*`` function returns a list of human-readable problems (empty = fine).
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

BASES = "TCAG"
AA_11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE_11 = {
    a + b + c: AA_11[16 * i + 4 * j + k]
    for i, a in enumerate(BASES)
    for j, b in enumerate(BASES)
    for k, c in enumerate(BASES)
}
STOPS = frozenset(c for c, a in CODE_11.items() if a == "*")
STARTS = frozenset({"ATG", "GTG", "TTG"})
SUBSETS = ("integrase", "repressor", "structural", "lysis", "cargo")
SUBSET_P = (0.05, 0.05, 0.4, 0.1, 0.4)
LOCUS_RE = re.compile(r"^GSIM_(\d{5,})$")
TU_RE = re.compile(r"^TU_(\d+)$")
XREF_RE = {
    "COG": re.compile(r"^COG\d{4}$"),
    "GO": re.compile(r"^GO:\d{7}$"),
    "KEGG": re.compile(r"^K\d{5}$"),
    "EC": re.compile(r"^\d+\.(\d+|-)\.(\d+|-)\.(\d+|-)$"),
}
#: The M2-side files ``write_contract_tables`` writes (C1-C8 and C15; C4 is JSON).
WRITTEN_TABLES = ("replicons", "genes", "xrefs", "tus", "layers", "regions", "region_genes")
M2_FILES = (
    "replicons.parquet", "genes.parquet", "xrefs.parquet", "annotation_meta.json",
    "tus.parquet", "layers.parquet", "regions.parquet", "region_genes.parquet",
    "features.parquet",  # C15 (M5-04): written always, header-only without an oriC
)  # fmt: skip
M3_FILES = ("peptide_hits.tsv", "protein_hits.tsv", "unmapped.tsv")


def translate(nt: str) -> str:
    """Table-11 translation of whole codons (stops shown as ``*``, initiator kept as coded)."""
    return "".join(CODE_11.get(nt[i : i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def codons(nt: str) -> list[str]:
    return [nt[i : i + 3] for i in range(0, len(nt), 3)]


def is_na(v: object) -> bool:
    if isinstance(v, list | tuple | np.ndarray):
        return False
    try:
        return bool(pd.isna(v))
    except (TypeError, ValueError):
        return False


def as_list(v: object) -> list[Any]:
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, list | tuple):
        return list(v)
    raise AssertionError(f"expected a list cell, got {type(v).__name__}: {v!r}")


def contig_order(replicons: pd.DataFrame) -> dict[str, int]:
    return {c: i for i, c in enumerate(replicons["contig"].tolist())}


def genes_by_contig(genes: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Genes per contig in genome order (start, end)."""
    return {
        str(c): g.sort_values(["start", "end"], kind="mergesort").reset_index(drop=True)
        for c, g in genes.groupby("contig", sort=False)
    }


def tu_of(tus: pd.DataFrame) -> dict[str, str]:
    return {
        g: str(t) for t, cell in zip(tus["tu_id"], tus["genes"], strict=True) for g in as_list(cell)
    }


def pair_prob_column(pairs: pd.DataFrame) -> str:
    for c in ("uniop_p", "p"):
        if c in pairs.columns:
            return c
    raise AssertionError(f"operon_pairs has no probability column: {list(pairs.columns)}")


def pair_map(pairs: pd.DataFrame) -> dict[frozenset[str], float]:
    col = pair_prob_column(pairs)
    return {
        frozenset((a, b)): float(p)
        for a, b, p in zip(pairs["gene_a"], pairs["gene_b"], pairs[col], strict=True)
    }


def adjacent_same_strand(
    genes: pd.DataFrame, replicons: pd.DataFrame, wrap: bool = False
) -> list[tuple[str, str]]:
    """Adjacent same-strand gene pairs in genome order (``wrap``: add circular last→first)."""
    topo = dict(zip(replicons["contig"], replicons["topology"], strict=True))
    out = []
    for c, g in genes_by_contig(genes).items():
        tags, strands = g["locus_tag"].tolist(), g["strand"].tolist()
        n = len(tags)
        for i in range(n - 1):
            if strands[i] == strands[i + 1]:
                out.append((tags[i], tags[i + 1]))
        if wrap and n > 2 and topo[c] == "circular" and strands[-1] == strands[0]:
            out.append((tags[-1], tags[0]))
    return out


def expected_uniop(
    genes: pd.DataFrame, replicons: pd.DataFrame, pairs: pd.DataFrame
) -> dict[str, tuple[float | None, float | None]]:
    """M2-20: p to the previous / next gene in transcription order (NA: other strand or absent).

    The neighbour is the genome-order neighbour on the same contig (on a circular contig the wrap
    neighbour counts only when the simulator lists that pair).
    """
    pm = pair_map(pairs)
    topo = dict(zip(replicons["contig"], replicons["topology"], strict=True))
    out: dict[str, tuple[float | None, float | None]] = {}
    for c, g in genes_by_contig(genes).items():
        tags, strands = g["locus_tag"].tolist(), g["strand"].tolist()
        n = len(tags)
        circ = topo[c] == "circular" and n > 2

        def p(i: int, j: int, tags: list[str] = tags, strands: list[str] = strands,
              n: int = n, circ: bool = circ) -> float | None:  # fmt: skip
            if not (0 <= j < n) and not circ:
                return None
            j %= n
            if j == i or strands[i] != strands[j]:
                return None
            return pm.get(frozenset((tags[i], tags[j])))

        for i in range(n):
            up, down = p(i, i - 1), p(i, i + 1)
            out[tags[i]] = (up, down) if strands[i] == "+" else (down, up)
    return out


def check_gene_coordinates(genes: pd.DataFrame, replicons: pd.DataFrame) -> list[str]:
    """Spec: gene span 3·(len_aa+1) bp, nt consistent with aa, sorted, inside the replicon."""
    errs: list[str] = []
    length = dict(zip(replicons["contig"], replicons["length"], strict=True))
    for r in genes.itertuples(index=False):
        nt, aa, ft = r.nt, r.aa, r.feature_type
        tag = r.locus_tag
        if r.end - r.start != len(nt):
            errs.append(f"{tag}: span {r.end - r.start} != len(nt) {len(nt)}")
        if len(nt) % 3:
            errs.append(f"{tag}: len(nt) {len(nt)} not a multiple of 3")
        cs = codons(nt)
        if not cs or cs[-1] not in STOPS:
            errs.append(f"{tag}: nt does not end with a stop codon")
        if cs and not is_na(r.start_codon) and cs[0] != r.start_codon:
            errs.append(f"{tag}: start_codon {r.start_codon} != first codon {cs[0]}")
        if cs and cs[0] not in STARTS:
            errs.append(f"{tag}: first codon {cs[0]} is not a start codon")
        internal = [i for i, c in enumerate(cs[:-1]) if c in STOPS]
        if ft == "pseudogene":
            if not internal:
                errs.append(f"{tag}: pseudogene without an internal in-frame stop")
        else:
            if internal:
                errs.append(f"{tag}: {ft} has an internal stop at codon {internal[0]}")
            if is_na(aa):
                errs.append(f"{tag}: {ft} with aa = NA")
            else:
                if len(nt) != 3 * (len(aa) + 1):
                    errs.append(f"{tag}: len(nt) {len(nt)} != 3*(len(aa)+1) = {3 * (len(aa) + 1)}")
                tr = translate(nt[3:-3])
                if aa[1:] != tr:
                    errs.append(f"{tag}: aa does not match the translation of nt")
                if aa[:1] != "M":
                    errs.append(f"{tag}: aa does not start with M")
            if ft == "sorf" and not is_na(aa) and len(aa) >= 30:
                errs.append(f"{tag}: sORF of {len(aa)} aa (needs < 30)")
        if not r.crosses_origin and r.end > length[r.contig]:
            errs.append(f"{tag}: end {r.end} beyond {r.contig} length {length[r.contig]}")
    for c, g in genes_by_contig(genes).items():
        starts = g["start"].to_numpy()
        if len(starts) > 1 and not (np.diff(starts) > 0).all():
            errs.append(f"{c}: gene starts not strictly increasing")
    return errs


def check_gaps(genes: pd.DataFrame, tus: pd.DataFrame) -> list[str]:
    """Spec: -20...+40 bp between genes of one TU, 50-300 bp between TUs (same contig)."""
    errs: list[str] = []
    tmap = tu_of(tus)
    for c, g in genes_by_contig(genes).items():
        tags, s, e = g["locus_tag"].tolist(), g["start"].tolist(), g["end"].tolist()
        for i in range(len(tags) - 1):
            gap = s[i + 1] - e[i]
            same = tmap.get(tags[i]) == tmap.get(tags[i + 1])
            lo, hi = (-20, 40) if same else (50, 300)
            if not lo <= gap <= hi:
                kind = "intra" if same else "inter"
                errs.append(f"{c}: {kind}-TU gap {gap} between {tags[i]} and {tags[i + 1]}")
    return errs


def check_tus(genes: pd.DataFrame, tus: pd.DataFrame, replicons: pd.DataFrame) -> list[str]:
    """TUs: partition, single strand and contig, consecutive genes, transcription order, bounds,
    ``TU_<n>`` numbered by contig order then start."""
    errs: list[str] = []
    g = genes.set_index("locus_tag")
    order = contig_order(replicons)
    seen: dict[str, int] = {}
    pos = {
        t: (c, i)
        for c, gg in genes_by_contig(genes).items()
        for i, t in enumerate(gg["locus_tag"].tolist())
    }
    keys = []
    for r in tus.itertuples(index=False):
        members = as_list(r.genes)
        if not members:
            errs.append(f"{r.tu_id}: empty TU")
            continue
        for m in members:
            seen[m] = seen.get(m, 0) + 1
        if any(m not in g.index for m in members):
            errs.append(f"{r.tu_id}: unknown genes")
            continue
        sub = g.loc[members]
        if set(sub["strand"]) != {r.strand}:
            errs.append(f"{r.tu_id}: strands {set(sub['strand'])} vs TU strand {r.strand}")
        if set(sub["contig"]) != {r.contig}:
            errs.append(f"{r.tu_id}: contigs {set(sub['contig'])} vs TU contig {r.contig}")
        starts = sub["start"].tolist()
        want = sorted(starts) if r.strand == "+" else sorted(starts, reverse=True)
        if starts != want:
            errs.append(f"{r.tu_id}: genes not in transcription order")
        idx = sorted(pos[m][1] for m in members)
        if idx != list(range(idx[0], idx[0] + len(idx))):
            errs.append(f"{r.tu_id}: genes not consecutive in genome order")
        if r.start != sub["start"].min() or r.end != sub["end"].max():
            errs.append(f"{r.tu_id}: bounds [{r.start},{r.end}) vs genes")
        if not TU_RE.match(str(r.tu_id)):
            errs.append(f"{r.tu_id}: id not TU_<n>")
        keys.append((order.get(r.contig, 10**9), int(sub["start"].min()), str(r.tu_id)))
    missing = sorted(set(g.index) - set(seen))
    multi = sorted(t for t, n in seen.items() if n > 1)
    if missing:
        errs.append(f"genes in no TU: {missing[:5]}")
    if multi:
        errs.append(f"genes in several TUs: {multi[:5]}")
    ids = [k[2] for k in sorted(keys)]
    if ids != [f"TU_{i + 1}" for i in range(len(ids))]:
        errs.append(f"TU ids not TU_1..TU_n by contig order then start: {ids[:6]} …")
    return errs


def check_locus_tags(genes: pd.DataFrame, replicons: pd.DataFrame) -> list[str]:
    """``GSIM_<5·i:05d>``, i = 1, 2, … across replicons in contig order then start (as M0)."""
    order = contig_order(replicons)
    g = genes.assign(_o=genes["contig"].map(order)).sort_values(["_o", "start"], kind="mergesort")
    want = [f"GSIM_{5 * (i + 1):05d}" for i in range(len(g))]
    got = g["locus_tag"].tolist()
    return [] if got == want else [f"locus tags {got[:4]} … != {want[:4]} …"]


def circular_consecutive(idx: Iterable[int], n: int) -> bool:
    """Indices form one run of consecutive genes on a circle of ``n`` genes."""
    s = sorted(set(idx))
    k = len(s)
    if k == 0 or k == n:
        return k == n
    breaks = sum(1 for i in range(k) if (s[(i + 1) % k] - s[i]) % n != 1)
    return breaks == 1


def members(region_genes: pd.DataFrame, region_id: str) -> list[str]:
    return region_genes.loc[region_genes["region_id"] == region_id, "locus_tag"].tolist()


def overlap_fraction(gs: int, ge: int, rs: int, re_: int, length: int, circular: bool) -> float:
    """Fraction of gene [gs, ge) inside region [rs, re_) (unwrapped if the region crosses)."""
    best = 0
    shifts = (0, length, -length) if circular else (0,)
    for sh in shifts:
        best = max(best, min(ge + sh, re_) - max(gs + sh, rs))
    return max(best, 0) / (ge - gs)


def chi2_gof_p(counts: Sequence[int], probs: Sequence[float]) -> float:
    from scipy.stats import chisquare

    n = sum(counts)
    return float(chisquare(counts, [n * p for p in probs]).pvalue)


def within_sd(k: float, n: int, p: float, z: float = 4.0) -> bool:
    sd = math.sqrt(n * p * (1 - p))
    return abs(k - n * p) <= z * sd + 1


def write_m3_headers(work: Path) -> None:
    """Header-only M3 tables (nothing placed) so the M4-06 loader sees a complete work dir."""
    from giap.genetable.contract import TABLES

    for key, name in (("peptide_hits", M3_FILES[0]), ("protein_hits", M3_FILES[1]),
                      ("unmapped", M3_FILES[2])):  # fmt: skip
        path = work / name
        if not path.exists():
            path.write_text("\t".join(c.name for c in TABLES[key]) + "\n", encoding="utf-8")


def load_work(work: Path) -> tuple[Any, list[Any]]:
    from giap.io.work_tables import default_paths, load_m4_inputs

    write_m3_headers(work)
    return load_m4_inputs(default_paths(work))


def errors(issues: Iterable[Any]) -> list[Any]:
    return [i for i in issues if i.severity == "error"]


def read_parquets(work: Path) -> Mapping[str, pd.DataFrame]:
    from giap.io.work_tables import FILE_NAMES

    return {k: pd.read_parquet(work / FILE_NAMES[k]) for k in WRITTEN_TABLES}
