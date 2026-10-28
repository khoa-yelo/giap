"""Operon scenarios: within-TU t increments, planted co-abundance and discordance, the pair table
and its truth (issue M6-03; GIPPA_SPEC.md §11.3, §8.4 ``discordance`` and ``missed_operons``, §8.5
attack 1, §7.4.5, §7.4.7).

``simulate_operons((LayoutConfig, ExpressionConfig), seed)`` draws a layout and its expression
truth (``simulate_layout``, ``simulate_expression``) and builds the §9.3 pair table the simulator
can know.
With an ``OperonConfig`` carrying ``operon_t = (loc, scale, df)``, every TU's first gene keeps its
merged-model abundance (TU level + stoichiometry + noise) and each later gene is the previous one
plus ``loc + scale * t_df`` (lead decision 8), so Delta of a true intra-TU pair (``b`` downstream)
is exactly that t. Planted pairs change that increment:

* ``n_missed`` true-NO pairs (adjacent, same strand, UniOP p < 0.5, different TUs) get
  ``Delta = missed_effect_sd * SD * z`` (``z ~ N(0, 1)``, SD the Delta SD of the both-detected true
  NO pairs before planting; lead decision 5): the downstream TU then follows its upstream partner;
* ``n_discordant`` true OP pairs get ``Delta = +-effect + loc + scale * t``: ceil(n/2) with a
  planted internal promoter (``+``, ``downstream_higher``), floor(n/2) with an internal terminator
  (``-``, ``upstream_higher``; lead decision 9).

Planted pairs are drawn with both genes detected first, then from the rest. With ``detection =
"abundance"`` (lead decision 10) the detection drawn by ``simulate_expression`` is provisional (it
picks the planted pairs and the background SD); after the walk and the planting the logit is
recomputed on the final ``abundance_log2_true`` and ``detected_true`` is redrawn from the operon
stream, so detection selects on the abundance the pairs compare. True elements (§11.3): a promoter
upstream of every TU start and a terminator after every TU end, and nothing inside a TU unless
planted (lead decision 1). Randomness: ``SeedSequence(seed).spawn(3)`` (layout, expression, operon
planting), so a plain ``ExpressionConfig`` gives the merged expression model unchanged. Nothing is
imported from ``gippa``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from gippasim.expression import (
    Expression,
    ExpressionConfig,
    detection_terms,
    simulate_expression,
)
from gippasim.layout import Layout, LayoutConfig, simulate_layout

__all__ = ["OperonConfig", "OperonSim", "simulate_operons"]

_RANK = {"OP": 0, "NO": 1, "CV": 2, "DV": 3, "SW": 4}
PAIR_TRUTH = ("true_operon", "promoter_between_true", "terminator_between_true")


@dataclass(frozen=True)
class OperonConfig(ExpressionConfig):
    """``ExpressionConfig`` plus the operon plan read by ``simulate_operons`` (only there)."""

    operon_t: tuple[float, float, float] | None = None  # (loc, scale, df) of the TU increment
    n_missed: int = 0
    missed_effect_sd: float = 0.0
    n_discordant: int = 0
    discordant_effect_log2: float = 0.0

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.operon_t is not None:
            object.__setattr__(self, "operon_t", tuple(float(v) for v in self.operon_t))
            if self.operon_t[1] <= 0 or self.operon_t[2] <= 0:
                raise ValueError(f"operon_t scale and df must be > 0, got {self.operon_t}")
        if min(self.n_missed, self.n_discordant) < 0 or self.missed_effect_sd < 0:
            raise ValueError("planted counts and missed_effect_sd must be >= 0")


@dataclass(frozen=True)
class OperonSim:
    """One replicate: layout, expression truth, the pair table and the truth frames."""

    layout: Layout
    expression: Expression
    pairs: pd.DataFrame
    truth: Mapping[str, pd.DataFrame]


def _skeleton(layout: Layout) -> pd.DataFrame:
    """Adjacent pairs (wrap pair on circular replicons of >= 3 genes) and SW triples, classed and
    oriented as ``gippa.operons.pairs`` does; ``pair_id`` in its order."""
    g = layout.genes
    reps = layout.replicons
    length = dict(zip(reps["contig"].astype(str), reps["length"].astype(int), strict=True))
    circ = set(reps.loc[reps["topology"] == "circular", "contig"].astype(str))
    rank = {c: k for k, c in enumerate(length)}
    tags = g["locus_tag"].astype(str).to_numpy()
    contig = g["contig"].astype(str).to_numpy()
    start = g["start"].to_numpy(np.int64)
    end = g["end"].to_numpy(np.int64)
    strand = g["strand"].astype(str).to_numpy()
    op = layout.operon_pairs
    raw = {(str(a), str(b)): float(p)
           for a, b, p in zip(op["gene_a"], op["gene_b"], op["uniop_p"], strict=True)}  # fmt: skip
    tu = dict(zip(layout.truth["locus_tag"].astype(str), layout.truth["tu_id"].astype(str),
                  strict=True))  # fmt: skip
    rows: list[tuple[Any, ...]] = []

    def dist(x: int, y: int, wrap: bool) -> int:
        return int(start[y] + (length[contig[x]] if wrap else 0) - end[x])

    for c in length:
        idx = np.flatnonzero(contig == c)
        idx = idx[np.argsort(start[idx], kind="stable")].tolist()
        n = len(idx)
        ring = c in circ and n >= 3
        chain = [(idx[i], idx[i + 1], False) for i in range(n - 1)]
        chain += [(idx[-1], idx[0], True)] if ring else []
        for left, right, wrap in chain:
            d = dist(left, right, wrap)
            if strand[left] == strand[right]:
                p = raw[(tags[left], tags[right])]
                cls = "OP" if p >= 0.5 else "NO"
                a, b = (right, left) if strand[left] == "-" else (left, right)
                same = tu[tags[a]] == tu[tags[b]]
                rows.append((cls, a, b, strand[a], d, p, same, not same, not same))
            else:
                cv = strand[left] == "+"
                rows.append(("CV" if cv else "DV", left, right, ".", d, math.nan, False,
                             not cv, cv))  # fmt: skip
        for i in range(n if ring else max(n - 2, 0)):
            x, y, z = idx[i], idx[(i + 1) % n], idx[(i + 2) % n]
            if strand[x] == strand[z] != strand[y]:
                d = dist(x, y, i == n - 1) + dist(y, z, i == n - 2)
                a, b = (z, x) if strand[x] == "-" else (x, z)
                rows.append(("SW", a, b, strand[a], d, math.nan, False, True, True))
    cols = ["class", "ia", "ib", "strand", "intergenic_bp", "uniop_p", *PAIR_TRUTH]
    df = pd.DataFrame(rows, columns=cols)
    df["_c"] = [rank[contig[a]] for a in df["ia"]]
    df["_sa"], df["_sb"] = start[df["ia"].to_numpy(int)], start[df["ib"].to_numpy(int)]
    df["_r"] = df["class"].map(_RANK)
    df["gene_a"], df["gene_b"] = tags[df["ia"].to_numpy(int)], tags[df["ib"].to_numpy(int)]
    df = df.sort_values(["_c", "_sa", "_sb", "_r", "gene_a", "gene_b"], kind="stable")
    df = df.drop(columns=["_c", "_sa", "_sb", "_r"]).reset_index(drop=True)
    df.insert(0, "pair_id", [f"PAIR_{k}" for k in range(1, len(df) + 1)])
    df["intergenic_bp"] = df["intergenic_bp"].astype("int64")
    for col in PAIR_TRUTH:
        df[col] = df[col].astype(bool)
    return df


def _pick(rng: np.random.Generator, cand: np.ndarray, both: np.ndarray, k: int) -> np.ndarray:
    """``k`` of the row indices ``cand``: both-detected rows first, random order within each."""
    key = rng.random(len(cand))
    order = np.lexsort((key, ~both[cand]))
    return np.sort(cand[order[:k]])


def _resolve(base: np.ndarray, parent: np.ndarray, inc: np.ndarray) -> np.ndarray:
    """``x[g] = x[parent[g]] + inc[g]`` (``base[g]`` where ``parent[g] < 0``)."""
    val = base.copy()
    done = parent < 0
    for g in range(len(base)):
        stack = []
        while not done[g]:
            stack.append(g)
            g = int(parent[g])
            if len(stack) > len(base):
                raise ValueError("planted pairs form a cycle")
        for h in reversed(stack):
            val[h] = val[parent[h]] + inc[h]
            done[h] = True
    return val


def _empty(cols: tuple[str, ...]) -> pd.DataFrame:
    return pd.DataFrame({c: pd.Series([], dtype=object) for c in cols})


def simulate_operons(cfgs: tuple[LayoutConfig, ExpressionConfig], seed: int) -> OperonSim:
    """One operon-scenario replicate of ``cfgs`` drawn from ``seed`` (see the module docstring)."""
    lcfg, ecfg = cfgs
    r_layout, r_expr, r_op = (np.random.default_rng(s)
                              for s in np.random.SeedSequence(seed).spawn(3))  # fmt: skip
    layout = simulate_layout(lcfg, r_layout)
    expr = simulate_expression(layout, ecfg, r_expr)
    pairs = _skeleton(layout)
    if not (pairs["class"] == "CV").any():
        raise ValueError("no convergent (CV) pair: the CV/SW decoys would not exist")
    ocfg = ecfg if isinstance(ecfg, OperonConfig) else OperonConfig()
    genes = expr.genes
    n = len(genes)
    det = genes["detected_true"].to_numpy(bool)
    ia, ib = pairs["ia"].to_numpy(int), pairs["ib"].to_numpy(int)
    both = det[ia] & det[ib]
    cls = pairs["class"].to_numpy()
    true_op = pairs["true_operon"].to_numpy(bool)
    t_draw = r_op.standard_t(ocfg.operon_t[2], n) if ocfg.operon_t else np.zeros(n)
    z_draw = r_op.standard_normal(n)
    base = genes["abundance_log2_true"].to_numpy(float).copy()
    parent = np.full(n, -1, dtype=np.int64)
    inc = np.zeros(n)
    if ocfg.operon_t is not None:
        loc, scale, _ = ocfg.operon_t
        pos = {t: i for i, t in enumerate(genes["locus_tag"].astype(str))}
        for members in layout.tus["genes"]:
            m = [pos[str(x)] for x in members]
            parent[m[1:]] = m[:-1]
            inc[m[1:]] = loc + scale * t_draw[m[1:]]
    disc = pd.DataFrame()
    if ocfg.n_discordant:
        cand = np.flatnonzero((cls == "OP") & true_op)
        if ocfg.n_discordant > len(cand):
            raise ValueError(f"n_planted={ocfg.n_discordant} discordant pairs exceed the "
                             f"{len(cand)} available true OP pairs")  # fmt: skip
        rows = _pick(r_op, cand, both, ocfg.n_discordant)
        n_prom = (ocfg.n_discordant + 1) // 2
        is_prom = np.zeros(len(rows), dtype=bool)
        is_prom[r_op.permutation(len(rows))[:n_prom]] = True
        eff = np.where(is_prom, 1.0, -1.0) * ocfg.discordant_effect_log2
        inc[ib[rows]] += eff
        pairs.loc[rows[is_prom], "promoter_between_true"] = True
        pairs.loc[rows[~is_prom], "terminator_between_true"] = True
        disc = pd.DataFrame({
            "pair_id": pairs["pair_id"].to_numpy()[rows],
            "gene_a": pairs["gene_a"].to_numpy()[rows],
            "gene_b": pairs["gene_b"].to_numpy()[rows],
            "direction": np.where(is_prom, "downstream_higher", "upstream_higher"),
            "element": np.where(is_prom, "promoter", "terminator"),
            "effect_log2": eff.astype(float),
        })  # fmt: skip
    ab = _resolve(base, parent, inc)
    missed = pd.DataFrame()
    if ocfg.n_missed:
        cand = np.flatnonzero((cls == "NO") & ~true_op & (pairs["uniop_p"].to_numpy() < 0.5))
        if ocfg.n_missed > len(cand):
            raise ValueError(f"n_planted={ocfg.n_missed} co-abundant pairs exceed the "
                             f"{len(cand)} available NO pairs")  # fmt: skip
        bg = cand[both[cand]]
        sd = float(np.std(ab[ib[bg]] - ab[ia[bg]], ddof=1)) if len(bg) > 1 else 1.0
        rows = _pick(r_op, cand, both, ocfg.n_missed)
        parent[ib[rows]] = ia[rows]
        inc[ib[rows]] = ocfg.missed_effect_sd * sd * z_draw[ib[rows]]
        ab = _resolve(base, parent, inc)
        missed = pairs.loc[rows, ["pair_id", "gene_a", "gene_b"]].reset_index(drop=True)
    if ocfg.operon_t is not None or ocfg.n_missed or ocfg.n_discordant:
        genes = genes.copy()
        genes["abundance_log2_true"] = ab
        if ocfg.detection == "abundance":  # lead decision 10: detect on the final abundance
            logit, shift, _, _ = detection_terms(layout, ocfg)
            logit += ocfg.beta["log2_abundance"] * ab
            logit += shift
            genes["p_detect_true"] = 1.0 / (1.0 + np.exp(-logit))
            det = r_op.random(n) < genes["p_detect_true"].to_numpy(float)
            genes["detected_true"] = det
        expr = Expression(genes, expr.islands)
    p_det = genes["p_detect_true"].to_numpy(float)
    obs = np.where(det, ab, np.nan)
    for side, ix in (("a", ia), ("b", ib)):
        pairs[f"detected_{side}"] = det[ix]
        pairs[f"p_detect_{side}"] = p_det[ix]
        pairs[f"abundance_log2_{side}"] = obs[ix]
    pairs["delta_log2"] = pairs["abundance_log2_b"] - pairs["abundance_log2_a"]
    pairs = pairs.drop(columns=["ia", "ib"])
    head = ("pair_id", "class", "gene_a", "gene_b", "strand", "intergenic_bp", "uniop_p")
    rest = [c for c in pairs.columns if c not in head and c not in PAIR_TRUTH]
    pairs = pairs[[*head, *rest, *PAIR_TRUTH]]
    t = ocfg.operon_t
    truth = {
        "pairs": pairs[["pair_id", "class", "gene_a", "gene_b", *PAIR_TRUTH]].copy(),
        "missed_operons": missed if len(missed) else _empty(("pair_id", "gene_a", "gene_b")),
        "discordant_pairs": disc if len(disc) else _empty(
            ("pair_id", "gene_a", "gene_b", "direction", "element", "effect_log2")),
        "operon_t": (pd.DataFrame({"loc": [t[0]], "scale": [t[1]], "df": [t[2]]}) if t
                     else _empty(("loc", "scale", "df"))),
    }  # fmt: skip
    return OperonSim(layout, expr, pairs, truth)
