"""The pair table: classes, orientation, distances and between-gene evidence (GIAP_SPEC.md
§7.4.1, §7.1.3-§7.1.5, §9.3; docs/specs/M6-04-pair-table.md).

Adjacent pairs (OP, NO, CV, DV) are the M2-12 chain with mixed-strand pairs kept; SW decoys are
genes i and i + 2 on one strand with gene i + 1 on the other. OP, NO and SW are oriented in the
direction of transcription (a upstream); CV and DV have a = the left gene. Evidence comes from
the M2-13/14/15 ``*_between`` helpers. Pure: no I/O, coordinates 0-based half-open.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import pandas as pd

from giap.annotate.operons import adjacent_pairs, intergenic_distance
from giap.annotate.promoter_calculator import tx_rate_between_arrays
from giap.annotate.promoters import promoter_between_tracks, promoter_tracks
from giap.annotate.terminators import Terminator, terminator_between
from giap.annotate.uniop_adapter import genome_order
from giap.config import GiapConfig
from giap.model.coords import Replicon
from giap.model.genes import Gene

PairClass = Literal["OP", "NO", "CV", "DV", "SW"]

#: §9.2 test columns, filled by M6-05/06; NA here with ``test = "none"``
TEST_COLUMNS: tuple[str, ...] = (
    "stat", "effect", "ci_low", "ci_high", "p_value", "q_value", "method", "null_model",
    "n_perm", "lfdr",
)  # fmt: skip
#: §9.3 ``pairs.tsv`` columns, in table order
PAIR_COLUMNS: tuple[str, ...] = (
    "pair_id", "class", "gene_a", "gene_b", "strand", "intergenic_bp", "uniop_p",
    "promoter_between", "tx_rate_between", "terminator_between", "detected_a", "detected_b",
    "p_detect_a", "p_detect_b", "abundance_log2_a", "abundance_log2_b", "delta_log2", "test",
    *TEST_COLUMNS, "flags",
)  # fmt: skip
RNA_BETWEEN = "rna_between"
RNA_TYPES: tuple[str, ...] = ("rrna", "trna", "ncrna", "tmrna")
_NA = math.nan
_CLASS_RANK = {"OP": 0, "NO": 1, "CV": 2, "DV": 3, "SW": 4}


@dataclass(frozen=True, slots=True)
class Evidence:
    """Between-gene evidence layers, built once per run; ``None`` = layer missing.

    ``windows`` is ``promoter_windows.parquet``, ``tx`` the Promoter Calculator TSS rates and
    ``rna`` the C15 rows of an RNA type (the caller filters by type)."""

    windows: pd.DataFrame | None
    terminators: tuple[Terminator, ...] | None
    tx: pd.DataFrame | None
    rna: pd.DataFrame


def adjacent_chain(
    genes: Sequence[Gene], replicons: Mapping[str, Replicon]
) -> list[tuple[Gene, Gene]]:
    """All adjacent pairs in genome order (left, right), wrap pair included on circular
    replicons; a 2-gene ring has one pair, not two."""
    seen: set[frozenset[str]] = set()
    out: list[tuple[Gene, Gene]] = []
    for a, b in adjacent_pairs(genes, replicons, same_strand=False):
        key = frozenset((a.locus_tag, b.locus_tag))
        if key not in seen:
            seen.add(key)
            out.append((a, b))
    return out


def sandwich_pairs(
    genes: Sequence[Gene], replicons: Mapping[str, Replicon]
) -> list[tuple[Gene, Gene, Gene]]:
    """Genome-order triples (i, i + 1, i + 2) with i, i + 2 on one strand and i + 1 on the other;
    wrapping on circular replicons with at least 3 genes."""
    out: list[tuple[Gene, Gene, Gene]] = []
    for contig, gs in genome_order(genes, list(replicons)).items():
        n = len(gs)
        ring = replicons[contig].topology == "circular" and n >= 3
        for i in range(n if ring else max(n - 2, 0)):
            x, y, z = gs[i], gs[(i + 1) % n], gs[(i + 2) % n]
            if x.interval.strand == z.interval.strand != y.interval.strand:
                out.append((x, y, z))
    return out


def classify(a: Gene, b: Gene, uniop_p: float | None, threshold: float) -> PairClass:
    """Class of the adjacent pair ``a`` (left) - ``b`` (right); OP iff ``uniop_p ≥ threshold``."""
    sa, sb = a.interval.strand, b.interval.strand
    if sa == sb:
        if uniop_p is None or math.isnan(uniop_p):
            raise ValueError(
                f"same-strand pair {a.locus_tag}-{b.locus_tag} has no uniop_p in pairs_raw"
            )
        return "OP" if uniop_p >= threshold else "NO"
    return "CV" if sa == "+" else "DV"


#: candidate margin around a gap for the terminator prefilter: ``terminator_between``'s window is
#: the gap widened by ``slack_bp`` (50) on the upstream side
_TERM_MARGIN_BP = 64


def _arrays(df: pd.DataFrame | None) -> dict[tuple[str, str], tuple[np.ndarray, np.ndarray]] | None:
    if df is None:
        return None
    return {
        (str(c), str(s)): (sub["tss"].to_numpy(), sub["tx_rate"].to_numpy())
        for (c, s), sub in df.groupby(["contig", "strand"], sort=False)
    }


def _by_strand(
    terms: tuple[Terminator, ...] | None,
) -> dict[tuple[str, str], tuple[np.ndarray, list[Terminator]]] | None:
    """Terminators per (contig, strand), sorted by start, with the start array for searchsorted."""
    if terms is None:
        return None
    out: dict[tuple[str, str], list[Terminator]] = {}
    for t in terms:
        out.setdefault((t.interval.contig, t.interval.strand), []).append(t)
    for v in out.values():
        v.sort(key=lambda t: t.interval.start)
    return {k: (np.array([t.interval.start for t in v], dtype=np.int64), v) for k, v in out.items()}


def _opt(v: float | None) -> float:
    return _NA if v is None else float(v)


_EMPTY = np.empty(0)


class _Lookup:
    """Evidence layers indexed once per run (per contig and strand), so each ``*_between`` call
    sees only its own strand and, for terminators, only the calls near the pair."""

    def __init__(self, ev: Evidence) -> None:
        self.tracks = None if ev.windows is None else promoter_tracks(ev.windows)
        self.tx = _arrays(ev.tx)
        self.terms = _by_strand(ev.terminators)
        self.rna = {
            str(c): (sub["start"].to_numpy(np.int64), sub["end"].to_numpy(np.int64))
            for c, sub in ev.rna.groupby("contig", sort=False)
        }

    def _near(self, key: tuple[str, str], lo: int, hi: int, rep: Replicon) -> list[Terminator]:
        """Terminators of ``key`` whose start lies in ``[lo, hi]`` (modulo the ring): a superset of
        what ``terminator_between`` can select for a gap inside ``[lo, hi]``."""
        assert self.terms is not None
        hit = self.terms.get(key)
        if hit is None:
            return []
        starts, terms = hit
        n = rep.length
        if rep.topology != "circular":
            ranges = [(lo, hi)]
        elif hi - lo + 1 >= n:
            return terms
        else:
            a = lo % n
            b = a + (hi - lo)
            ranges = [(a, b)] if b < n else [(a, n - 1), (0, b - n)]
        out: list[Terminator] = []
        for x, y in ranges:
            out += terms[np.searchsorted(starts, x, "left") : np.searchsorted(starts, y, "right")]
        return out

    def evidence(
        self, up: Gene, down: Gene, left: Gene, right: Gene, rep: Replicon
    ) -> tuple[float, float, float]:
        """Promoter, Tx-rate and terminator evidence of ``up`` -> ``down`` (``left``/``right`` =
        the same two genes in genome order)."""
        key = (up.interval.contig, up.interval.strand)
        prom = tx = term = _NA
        if self.tracks is not None:
            prom = _opt(promoter_between_tracks(up, down, self.tracks, rep))
        if self.tx is not None:
            tss, rate = self.tx.get(key, (_EMPTY, _EMPTY))
            tx = _opt(tx_rate_between_arrays(up, down, tss, rate, rep))
        if self.terms is not None:
            lo = left.interval.end
            hi = lo + intergenic_distance(left, right, rep)
            near = self._near(
                key, min(lo, hi) - _TERM_MARGIN_BP, max(lo, hi) + _TERM_MARGIN_BP, rep
            )
            # yes = 0.0: the raw maximum; the cut-off is applied in M6-07 (decision 5)
            term = _opt(terminator_between(up, down, near, rep, yes=0.0))
        return prom, tx, term

    def rna_between(self, gaps: Sequence[tuple[Gene, Gene]], rep: Replicon) -> bool:
        """Any RNA feature overlapping the open gap between consecutive genes (any strand)."""
        hit = self.rna.get(rep.id)
        if hit is None:
            return False
        s, e = hit
        shifts = (-rep.length, 0, rep.length) if rep.topology == "circular" else (0,)
        for left, right in gaps:
            lo = left.interval.end
            hi = lo + intergenic_distance(left, right, rep)
            if hi <= lo:
                continue
            if any(bool(((s + k < hi) & (e + k > lo)).any()) for k in shifts):
                return True
        return False


def _gene_columns(genes_tbl: pd.DataFrame, genes: Sequence[Gene]) -> pd.DataFrame:
    tbl = genes_tbl.drop_duplicates("locus_tag").set_index("locus_tag")
    missing = [g.locus_tag for g in genes if g.locus_tag not in tbl.index]
    if missing:
        raise ValueError(f"genes absent from the gene table: {', '.join(missing[:10])}")
    return tbl[["detected", "p_detect", "abundance_log2"]]


def _empty_frame() -> pd.DataFrame:
    return pd.DataFrame({c: pd.Series(dtype=object) for c in PAIR_COLUMNS})


def pair_frame(
    genes_tbl: pd.DataFrame,
    genes: Sequence[Gene],
    replicons: Mapping[str, Replicon],
    pairs_raw: pd.DataFrame,
    ev: Evidence,
    cfg: GiapConfig,
) -> pd.DataFrame:
    """The ``tables/pairs.tsv`` frame: one row per adjacent pair plus the SW decoys, ``pair_id``
    ``PAIR_<n>`` over contig (C1 order), ``start_a``, ``start_b``, class (never ``uniop_p``)."""
    per_gene = _gene_columns(genes_tbl, genes)
    threshold = cfg.annotation.uniop.pair_threshold
    raw = {
        (str(a), str(b)): float(p)
        for a, b, p in pairs_raw[["gene_a", "gene_b", "uniop_p"]].itertuples(index=False)
    }
    look = _Lookup(ev)
    order = {c: k for k, c in enumerate(replicons)}
    rows: list[dict[str, Any]] = []

    def emit(
        cls: PairClass, a: Gene, b: Gene, d: int, p: float, gaps: Sequence[tuple[Gene, Gene]]
    ) -> None:
        rep = replicons[a.interval.contig]
        oriented = cls not in ("CV", "DV")
        ends = (gaps[0][0], gaps[-1][1])  # the pair in genome order
        prom, tx, term = look.evidence(a, b, *ends, rep) if oriented else (_NA, _NA, _NA)
        rows.append({
            "_contig": order[a.interval.contig], "_start_a": a.interval.start,
            "_start_b": b.interval.start, "_rank": _CLASS_RANK[cls],
            "class": cls, "gene_a": a.locus_tag, "gene_b": b.locus_tag,
            "strand": a.interval.strand if oriented else ".", "intergenic_bp": d,
            "uniop_p": p, "promoter_between": prom, "tx_rate_between": tx,
            "terminator_between": term,
            "flags": [RNA_BETWEEN] if look.rna_between(gaps, rep) else [],
        })  # fmt: skip

    for left, right in adjacent_chain(genes, replicons):
        rep = replicons[left.interval.contig]
        same = left.interval.strand == right.interval.strand
        p = raw.get((left.locus_tag, right.locus_tag)) if same else None
        cls = classify(left, right, p, threshold)
        a, b = (right, left) if same and left.interval.strand == "-" else (left, right)
        emit(cls, a, b, intergenic_distance(left, right, rep), _opt(p), [(left, right)])
    for x, y, z in sandwich_pairs(genes, replicons):
        rep = replicons[x.interval.contig]
        d = intergenic_distance(x, y, rep) + intergenic_distance(y, z, rep)
        a, b = (z, x) if x.interval.strand == "-" else (x, z)
        emit("SW", a, b, d, _NA, [(x, y), (y, z)])

    if not rows:
        return _empty_frame()
    df = pd.DataFrame(rows).sort_values(
        ["_contig", "_start_a", "_start_b", "_rank", "gene_a", "gene_b"], kind="stable"
    )
    df = df.drop(columns=["_contig", "_start_a", "_start_b", "_rank"]).reset_index(drop=True)
    df.insert(0, "pair_id", [f"PAIR_{n}" for n in range(1, len(df) + 1)])
    for side in ("a", "b"):
        got = per_gene.loc[df[f"gene_{side}"]]
        for col in ("detected", "p_detect", "abundance_log2"):
            df[f"{col}_{side}"] = got[col].to_numpy()
    df["delta_log2"] = df["abundance_log2_b"].astype(float) - df["abundance_log2_a"].astype(float)
    df["test"] = "none"
    for col in TEST_COLUMNS:
        df[col] = _NA if col not in ("method", "null_model", "n_perm") else None
    df["n_perm"] = pd.array([pd.NA] * len(df), dtype="Int64")
    df["intergenic_bp"] = df["intergenic_bp"].astype("int64")
    return df[list(PAIR_COLUMNS)]
