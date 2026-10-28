"""Function annotations, ORA/GSEA nulls and attacks, planted determinants (issue M8-14;
GIPPA_SPEC.md §7.7.5, §8.4, §8.5 attacks 2 and 3, §11.3, §11.4).

Every synthetic ID lies in a reserved range (lead decision 1): GO terms ``GO:9000000`` up to
``GO:9999998`` (``GO:9999999`` is the fixtures' absent-term marker and is never emitted, decision
3 gap 10), KEGG orthologs ``K9xxxx`` and KEGG maps ``map9xxxx``; COG keeps its letters on
synthetic ``COG<digits>`` IDs (gap 9). Nothing here is derived from real GO, KEGG or COG content.

``synthetic_ontology`` draws a GO-shaped DAG (three namespaces, one root each, ``is_a`` and
``part_of``) and annotates 30-60 % of the genes directly (M8-03 propagates). ``draw(pair, rng)``
turns a ``gippasim.scenarios`` pair into a :class:`World`: the layout (its C3 ``xrefs`` replaced by
the synthetic GO/KEGG/COG xrefs), the expression truth, the OBO text, the KEGG link rows, the
gene sets M8-03 would build (``sets``), the determinants features and the truth frames.

Sets (``sets``, one row per member) follow M8-03: GO annotations propagated over ``is_a`` and
``part_of`` and split by namespace, KEGG maps through the link rows, COG letters through the
``cog_def`` table; kept when they hold 10-500 TUs (the synthetic roots, 30-60 % of the genes, are
above that at the calibration sizes and drop out, decision 3 gap 3). ``gsea_attack_operon_packed``
lists its packed KEGG sets only (unfiltered): 1-3 packed operons of 10-30 genes plus
``PACKED_SINGLES`` single-gene TUs each, whole TUs; packed TUs draw their level like every other
TU (the shared TU level is the attack; lead decision 4, B2).

Detection is ``Bernoulli(p_true)`` with ``p_true`` the M4-09 covariate logistic; it is drawn
from its own stream, so the attack at shift 0 reproduces ``ora_null``, and the determinants null
and planted scenarios share layout and detection for one seed. ``ora_attack_confounded`` carries
``member_shift`` through ``n_obs`` (``CONFOUND_COVARIATE``; lead decision 6, B1): each confounded
member with an ``aa`` and ``n_obs > 0`` has its K/R residues moved (with their codons in ``nt``;
composition, hence GRAVY and length, unchanged) so that its recomputed observable-peptide count
gives the logit shift nearest ``member_shift``; ``n_obs_true`` and ``p_detect_true`` follow the
edited sequence, and the logit stays the covariate model's. The edits are deterministic (no draw),
so every other scenario is unchanged. Determinants: ``abundance_log2_true
= mu + TU effect (sd 1) + noise (noise_sd) + sum_j beta_j z_j`` with ``z_j`` the population
(ddof 0) z-score of feature ``j``; the true R² is ``sum beta² / (sum beta² + 1 + noise_sd²)``.
Planted features: a sigma-70 box score per TU (``TTGACA``/``TATAAT`` matches, spacer penalty), a
Shine-Dalgarno score per gene (``AGGAGG`` matches, spacer penalty) and the CAI of the gene's
back-translated ``nt`` against the bundled codon usage. Nothing is imported from ``gippa``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from types import MappingProxyType
from typing import Any, Literal

import numpy as np
import pandas as pd

from gippasim.codons import CODE_11, SENSE_USAGE
from gippasim.expression import (
    Expression,
    ExpressionConfig,
    FrozenMap,
    simple_covariates,
    simulate_expression,
)
from gippasim.layout import Layout, LayoutConfig, simulate_layout

__all__ = [
    "BLOCKS",
    "CONFOUND_COVARIATE",
    "FEATURES",
    "PACKED_TU_SIZE",
    "RESERVED_PREFIXES",
    "SET_TU_RANGE",
    "FunctionExpressionConfig",
    "FunctionSpec",
    "World",
    "draw",
    "synthetic_ontology",
    "write_truth",
]

#: M4-13 covariate carrying ``ora_attack_confounded``'s shift (lead decision 6, B1).
CONFOUND_COVARIATE: Literal["n_obs"] = "n_obs"
_PEPTIDE_LEN = (7, 30)  # observable peptide lengths (gippasim.expression / GIPPA defaults)
#: Reserved ID prefixes of the synthetic annotations (lead decision 1).
RESERVED_PREFIXES: tuple[str, ...] = ("GO:9", "K9", "map9")  # GO term, KEGG KO, KEGG map
_GO_BASE, _GO_LAST = 9_000_000, 9_999_998  # GO:9999999 is the fixtures' absent marker
_KO_BASE = _MAP_BASE = 90_000
NAMESPACES = ("biological_process", "molecular_function", "cellular_component")
_COLLECTION = dict(zip(NAMESPACES, ("go_bp", "go_mf", "go_cc"), strict=True))
#: Determinants features and their blocks (spec §11.3).
FEATURES: tuple[str, ...] = ("promoter", "rbs", "cai")
BLOCKS: Mapping[str, str] = MappingProxyType(
    {"promoter": "transcription", "rbs": "translation", "cai": "codon"}
)
SET_TU_RANGE = (10, 500)  # sets kept when they hold this many TUs (M8-03 gsea_size)
PACKED_TU_SIZE = (10, 30)  # genes per packed operon
PACKED_SINGLES = 40  # single-gene chr1 TUs per packed set (fewer on small layouts, >= 9)
_DET_TU_SD = 1.0  # determinants: TU random effect SD
_COG_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
_P_COG, _P_KO, _P_NS = 0.65, 0.4, 0.75
_LEVEL_GROWTH = 1.6  # expected terms per level grow geometrically with depth
_P_EXTRA_PARENT = 0.25
_SIGMA70 = "TTGACA" + "TATAAT"
_SIGMA70_W = np.array([0.9, 1.0, 0.8, 0.7, 0.6, 0.5, 0.6, 1.0, 0.8, 0.7, 0.6, 0.9])
_SD_W = np.array([0.6, 1.0, 0.9, 1.0, 0.8, 0.5])
_XREF_COLUMNS = ["locus_tag", "db", "id"]
SET_COLUMNS = ["set_id", "collection", "locus_tag", "confounded"]


@dataclass(frozen=True)
class FunctionSpec:
    """What ``draw`` adds to the M4 layout and expression for one M8-14 scenario."""

    kind: Literal["ora", "gsea", "packed", "determinants"]
    attack: bool = False  # ora: confounded categories present
    member_shift: float = 0.0  # ora attack: logit shift of the confounded members
    n_sets: int = 0  # packed: number of packed sets
    coefs: Mapping[str, float] = FrozenMap()  # determinants: beta per feature
    features: Literal["planted", "gaussian"] = "planted"

    def __post_init__(self) -> None:
        object.__setattr__(self, "coefs", FrozenMap(self.coefs))


@dataclass(frozen=True)
class FunctionExpressionConfig(ExpressionConfig):
    """An ``ExpressionConfig`` carrying the scenario's :class:`FunctionSpec`."""

    function: FunctionSpec | None = None


@dataclass(frozen=True)
class World:
    """One draw of a function scenario (see the module docstring)."""

    layout: Layout
    expression: Expression
    obo: str
    kegg_links: pd.DataFrame  # ko (``ko:K9xxxx``), pathway (``path:map9xxxx``)
    kegg_names: pd.DataFrame  # pathway, name
    cog_def: pd.DataFrame  # cog_id, cog_category
    cog_names: pd.DataFrame  # category, name
    sets: pd.DataFrame  # SET_COLUMNS, one row per member
    features: pd.DataFrame | None  # determinants: locus_tag, promoter, rbs, cai
    truth: Mapping[str, pd.DataFrame]


# --------------------------------------------------------------------------- ontology


@dataclass(frozen=True)
class _Dag:
    ids: list[str]
    ns: np.ndarray  # namespace index per term
    level: np.ndarray  # 0 = root
    is_a: list[list[int]]
    part_of: list[list[int]]


def _pick(u: float, pool: Sequence[int]) -> int:
    return pool[min(int(u * len(pool)), len(pool) - 1)]


def _dag(n_terms: int, depth: int, rng: np.random.Generator) -> _Dag:
    if not 3 <= n_terms <= _GO_LAST - _GO_BASE + 1:
        raise ValueError(f"n_terms must be in [3, {_GO_LAST - _GO_BASE + 1}], got {n_terms}")
    if depth < 1:
        raise ValueError(f"depth must be >= 1, got {depth}")
    m = n_terms - 3
    ns = np.concatenate([[0, 1, 2], rng.integers(0, 3, size=m)])
    w = _LEVEL_GROWTH ** np.arange(1, depth + 1)
    raw = 1 + np.searchsorted(np.cumsum(w / w.sum()), rng.random(m), side="right")
    raw = np.minimum(raw, depth)
    level = np.zeros(n_terms, dtype=np.int64)
    for k in range(3):  # sorted levels per namespace, made contiguous (no empty level)
        idx = np.flatnonzero(ns[3:] == k) + 3
        top = 0
        for i, r in zip(idx.tolist(), np.sort(raw[idx - 3]).tolist(), strict=True):
            top = min(int(r), top + 1)
            level[i] = top
    u = rng.random((4, n_terms))
    pools: dict[tuple[int, int], list[int]] = {}
    for i in range(n_terms):
        pools.setdefault((int(ns[i]), int(level[i])), []).append(i)
    is_a: list[list[int]] = [[] for _ in range(n_terms)]
    part_of: list[list[int]] = [[] for _ in range(n_terms)]
    for i in range(3, n_terms):
        k, lv = int(ns[i]), int(level[i])
        primary = _pick(u[0, i], pools[(k, lv - 1)])
        is_a[i].append(primary)
        if u[1, i] < _P_EXTRA_PARENT and lv >= 2:
            above = [j for d in range(lv) for j in pools[(k, d)] if j != primary]
            if above:
                extra = _pick(u[2, i], above)
                (part_of if u[3, i] < 0.5 else is_a)[i].append(extra)
    ids = [f"GO:{_GO_BASE + i:07d}" for i in range(n_terms)]
    return _Dag(ids, ns, level, is_a, part_of)


def _obo(dag: _Dag) -> str:
    out = ["format-version: 1.2", "data-version: gippasim/synthetic", "ontology: go", ""]
    for i, t in enumerate(dag.ids):
        ns = NAMESPACES[int(dag.ns[i])]
        name = f"synthetic {ns} root" if i < 3 else f"synthetic term {i}"
        out += ["[Term]", f"id: {t}", f"name: {name}", f"namespace: {ns}"]
        out += [f"is_a: {dag.ids[p]}" for p in dag.is_a[i]]
        out += [f"relationship: part_of {dag.ids[p]}" for p in dag.part_of[i]]
        out.append("")
    return "\n".join(out)


def _annotate(dag: _Dag, genes: Sequence[str], rng: np.random.Generator) -> pd.DataFrame:
    """Direct GO annotations: 30-60 % of ``genes``, one term in each namespace w.p. 0.75."""
    tags = list(dict.fromkeys(str(g) for g in genes))
    n = len(tags)
    share = rng.uniform(0.35, 0.55)
    on = rng.random(n) < share
    take = (rng.random((3, n)) < _P_NS) & on
    take[0] |= on & ~take.any(axis=0)  # every annotated gene gets a term
    u = rng.random((3, n))
    rows: list[tuple[int, str]] = []
    for k in range(3):
        pool = [i for i in range(3, len(dag.ids)) if dag.ns[i] == k] or [k]
        pick = np.minimum((u[k] * len(pool)).astype(np.int64), len(pool) - 1)
        rows += [(g, dag.ids[pool[int(pick[g])]]) for g in np.flatnonzero(take[k]).tolist()]
    rows = sorted(set(rows))
    return pd.DataFrame(
        {"locus_tag": [tags[g] for g, _ in rows], "db": ["GO"] * len(rows),
         "id": [t for _, t in rows]},
        columns=_XREF_COLUMNS,
    )  # fmt: skip


def synthetic_ontology(
    n_terms: int, depth: int, rng: np.random.Generator, *, genes: Sequence[str] = ()
) -> tuple[str, pd.DataFrame]:
    """(OBO text of ``n_terms`` live ``GO:9xxxxxx`` terms of height <= ``depth``, the C3-shaped
    ``locus_tag, db, id`` table of direct annotations of ``genes``)."""
    dag = _dag(n_terms, depth, rng)
    return _obo(dag), _annotate(dag, genes, rng)


# --------------------------------------------------------------------------- annotations


def _n_terms(n_genes: int) -> int:
    return int(min(600, max(60, n_genes // 6)))


def _kegg(n_genes: int, rng: np.random.Generator) -> tuple[pd.DataFrame, list[str]]:
    """Link rows of ``n_ko`` synthetic KOs to 2-20 synthetic maps, and the KO list."""
    n_path = int(rng.integers(2, 21))
    n_ko = max(20, n_genes // 20)
    paths = [f"map{_MAP_BASE + i:05d}" for i in range(n_path)]
    kos = [f"K{_KO_BASE + i:05d}" for i in range(n_ko)]
    first = np.arange(n_ko) % n_path
    first = np.where(np.arange(n_ko) < n_path, first, rng.integers(0, n_path, size=n_ko))
    second = rng.integers(0, n_path, size=n_ko)
    extra = rng.random(n_ko) < 0.3
    rows = []
    for i, ko in enumerate(kos):
        for j in sorted({int(first[i])} | ({int(second[i])} if extra[i] else set())):
            rows.append((f"ko:{ko}", f"path:{paths[j]}"))
    return pd.DataFrame(rows, columns=["ko", "pathway"]), kos


def _names(paths: Sequence[str]) -> pd.DataFrame:
    return pd.DataFrame({"pathway": [f"path:{p}" for p in paths],
                         "name": [f"synthetic pathway {p[3:]}" for p in paths]})  # fmt: skip


def _cog(tags: list[str], rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame]:
    n_cog = max(26, len(tags) // 12)
    ids = [f"COG{i + 1:04d}" for i in range(n_cog)]
    letters = [_COG_LETTERS[i] for i in rng.integers(0, len(_COG_LETTERS), size=n_cog)]
    on = rng.random(len(tags)) < _P_COG
    pick = rng.integers(0, n_cog, size=len(tags))
    xr = pd.DataFrame({"locus_tag": [tags[i] for i in np.flatnonzero(on)],
                       "db": "COG", "id": [ids[int(pick[i])] for i in np.flatnonzero(on)]},
                      columns=_XREF_COLUMNS)  # fmt: skip
    return xr, pd.DataFrame({"cog_id": ids, "cog_category": letters})


# --------------------------------------------------------------------------- sets


def _closure(dag: _Dag) -> list[frozenset[int]]:
    up: list[frozenset[int]] = [frozenset()] * len(dag.ids)
    for i in np.argsort(dag.level, kind="stable").tolist():  # parents have lower levels
        s = {i}
        for p in (*dag.is_a[i], *dag.part_of[i]):
            s |= up[p]
        up[i] = frozenset(s)
    return up


def _collect(
    tags: list[str], tu_of: Mapping[str, str], xrefs: pd.DataFrame, dag: _Dag,
    links: pd.DataFrame, cog_def: pd.DataFrame,
) -> pd.DataFrame:  # fmt: skip
    """M8-03-shaped sets (propagated GO by namespace, KEGG maps, COG letters), 10-500 TUs."""
    index = {t: i for i, t in enumerate(dag.ids)}
    up = _closure(dag)
    path: dict[str, set[str]] = {}
    for ko, pw in zip(links["ko"], links["pathway"], strict=True):
        path.setdefault(str(ko).removeprefix("ko:"), set()).add(str(pw).removeprefix("path:"))
    letter = dict(zip(cog_def["cog_id"], cog_def["cog_category"], strict=True))
    members: dict[tuple[str, str], set[str]] = {}
    for g, db, xid in zip(xrefs["locus_tag"], xrefs["db"], xrefs["id"], strict=True):
        if db == "GO":
            for t in up[index[xid]]:
                key = (_COLLECTION[NAMESPACES[int(dag.ns[t])]], dag.ids[t])
                members.setdefault(key, set()).add(g)
        elif db == "KEGG":
            for mp in path.get(xid, ()):
                members.setdefault(("kegg", mp), set()).add(g)
        elif db == "COG":
            members.setdefault(("cog", letter[xid]), set()).add(g)
    rank = {g: i for i, g in enumerate(tags)}
    order = {"cog": 0, "kegg": 1, "go_bp": 2, "go_mf": 3, "go_cc": 4}
    rows: list[tuple[str, str, str, bool]] = []
    lo, hi = SET_TU_RANGE
    for (c, sid), m in sorted(members.items(), key=lambda kv: (order[kv[0][0]], kv[0][1])):
        if lo <= len({tu_of[g] for g in m}) <= hi:
            rows += [(sid, c, g, False) for g in sorted(m, key=rank.__getitem__)]
    return pd.DataFrame(rows, columns=SET_COLUMNS)


# --------------------------------------------------------------------------- features


def _cai(nts: Sequence[Any]) -> np.ndarray:
    """CAI of each ``nt`` against ``SENSE_USAGE`` (single-codon residues and stops skipped)."""
    best: dict[str, float] = {}
    for c, f in SENSE_USAGE.items():
        best[CODE_11[c]] = max(best.get(CODE_11[c], 0.0), f)
    syn = {a: sum(CODE_11[c] == a for c in SENSE_USAGE) for a in best}
    logw = np.full(64, np.nan)
    code = {ord(b): i for i, b in enumerate("TCAG")}
    for c, f in SENSE_USAGE.items():
        if syn[CODE_11[c]] > 1:
            logw[16 * code[ord(c[0])] + 4 * code[ord(c[1])] + code[ord(c[2])]] = math.log(
                f / best[CODE_11[c]]
            )
    lut = np.zeros(256, dtype=np.int64)
    for byte, i in code.items():
        lut[byte] = i
    strs = [str(s)[: len(str(s)) - len(str(s)) % 3] for s in nts]
    buf = np.frombuffer("".join(strs).encode(), dtype=np.uint8)
    b = lut[buf].reshape(-1, 3)
    lw = logw[16 * b[:, 0] + 4 * b[:, 1] + b[:, 2]]
    ok = ~np.isnan(lw)
    lens = np.array([len(s) // 3 for s in strs], dtype=np.int64)
    seg = np.repeat(np.arange(len(strs)), lens)
    tot = np.bincount(seg[ok], weights=lw[ok], minlength=len(strs))
    cnt = np.bincount(seg[ok], minlength=len(strs))
    out: np.ndarray = np.exp(tot / np.maximum(cnt, 1))
    return out


def _box_score(rng: np.random.Generator, n: int, weights: np.ndarray, spacers: tuple[int, int],
               best: int, penalty: float) -> np.ndarray:  # fmt: skip
    """Weighted consensus matches of ``n`` drawn boxes minus a spacer penalty: each position
    matches w.p. q + (1 - q) / 4 with a per-box strength q ~ Beta(2, 2)."""
    q = rng.beta(2.0, 2.0, size=n)
    match = rng.random((n, len(weights))) < (q + (1.0 - q) / 4.0)[:, None]
    spacer = rng.integers(spacers[0], spacers[1] + 1, size=n)
    out: np.ndarray = match.astype(float) @ weights - penalty * np.abs(spacer - best)
    return out


def _features(layout: Layout, mode: str, rng: np.random.Generator) -> pd.DataFrame:
    tags = layout.genes["locus_tag"].astype(str).tolist()
    n = len(tags)
    if mode == "gaussian":
        z = rng.standard_normal((n, 3))
        return pd.DataFrame({"locus_tag": tags, "promoter": z[:, 0], "rbs": z[:, 1],
                             "cai": z[:, 2]})  # fmt: skip
    tus = layout.tus
    tu_score = _box_score(rng, len(tus), _SIGMA70_W, (15, 19), 17, 0.5)
    by_tu = {str(g): float(tu_score[i]) for i, gs in enumerate(tus["genes"]) for g in gs}
    rbs = _box_score(rng, n, _SD_W, (5, 10), 7, 0.3)
    cai = _cai(layout.genes["nt"].tolist())
    return pd.DataFrame({"locus_tag": tags, "promoter": [by_tu[g] for g in tags], "rbs": rbs,
                         "cai": cai})  # fmt: skip


def _z(x: Any) -> np.ndarray:
    a = np.asarray(x, dtype=float)
    sd = a.std()
    out: np.ndarray = (a - a.mean()) / sd if sd > 0 else np.zeros_like(a)
    return out


# --------------------------------------------------------------------------- draw


def _sigmoid(x: np.ndarray) -> np.ndarray:
    out: np.ndarray = 1.0 / (1.0 + np.exp(-x))
    return out


def _tu_of(layout: Layout) -> dict[str, str]:
    return {str(g): str(t) for t, gs in zip(layout.tus["tu_id"], layout.tus["genes"],
                                            strict=True) for g in gs}  # fmt: skip


def _redetect(genes: pd.DataFrame, delta: np.ndarray, rng: np.random.Generator) -> pd.DataFrame:
    p = genes["p_detect_true"].to_numpy(float)
    logit = np.log(p) - np.log1p(-p)
    q = _sigmoid(logit + delta)
    return genes.assign(p_detect_true=q, detected_true=rng.random(len(q)) < q)


def _chunks(xs: list[int], kr: list[int], sizes: list[int]) -> tuple[list[int], int]:
    """(position 0, then ``sizes[j]`` non-K/R positions and K/R ``j`` per chunk; positions used)."""
    order, at = [0], 0
    for j, size in enumerate(sizes):
        order += [*xs[at : at + size], kr[j]]
        at += size
    return order, at


def _kr_orders(aa: str, k: int) -> list[list[int]]:
    """Residue orders of ``aa`` (position 0 fixed; non-K/R residues keep their order) with ``k``
    leading K/R-terminated fragments of 7-30 aa, then either one fragment longer than 30 aa and
    the remaining K/R as single residues at the C terminus, or (short proteins) the rest cut by
    the remaining K/R into fragments shorter than 7 aa; candidates in preference order."""
    lo, hi = _PEPTIDE_LEN
    kr = [i for i in range(1, len(aa)) if aa[i] in "KR"]
    xs = [i for i in range(1, len(aa)) if aa[i] not in "KR"]
    m, nx = len(kr), len(xs)
    if k > m:
        return []
    # chunk j holds x_j non-K/R residues (+ M for j = 0) and one K/R
    mins = [lo - 2] + [lo - 1] * (k - 1) if k else []
    maxs = [hi - 2] + [hi - 1] * (k - 1) if k else []
    smin, smax = sum(mins), sum(maxs)
    lead, trail = (1 if k == 0 else 0), (1 if m > k else 0)

    def tail_len(total: int) -> int:
        return lead + nx - total + trail

    totals = [t for t in range(smin, min(smax, nx) + 1)
              if not lo <= tail_len(t) <= hi]  # fmt: skip
    totals.sort(key=lambda t: (tail_len(t) < lo, t))  # long last fragment first
    out: list[list[int]] = []
    for total in totals[:4]:
        sizes, extra = list(mins), total - smin
        for j in range(k):  # spread the extra residues evenly, capped at the maxima
            add = min(maxs[j] - sizes[j], -(-extra // (k - j)))
            sizes[j] += add
            extra -= add
        order, at = _chunks(xs, kr, sizes)
        out.append(order + xs[at:] + kr[k:])
    # short proteins: fill the k chunks, cut the rest into pieces of < lo aa by the remaining K/R
    total = min(smax, nx)
    rest, r = nx - total, m - k
    if total >= smin and rest <= r * (lo - 2) + (lo - 1) - lead:
        sizes, extra = list(mins), total - smin
        for j in range(k):
            add = min(maxs[j] - sizes[j], extra)
            sizes[j] += add
            extra -= add
        order, at = _chunks(xs, kr, sizes)
        for j in range(r + 1):  # piece j: up to lo - 2 residues (+ M when k = 0) and one K/R
            cap = (lo - 2 if j < r else lo - 1) - (lead if j == 0 else 0)
            take = min(cap, -(-(nx - at) // (r + 1 - j)))
            order += [*xs[at : at + take], *kr[k + j : k + j + 1]]  # no K/R after the last
            at += take
        if at == nx:
            out.append(order)
    return out


#: Largest |realised mean member shift - member_shift| accepted (M8-14 decision 6 band).
SHIFT_TOL = 0.15


def _move_kr(aa: str, nt: Any, n0: int, shift: float, b: float) -> tuple[str, Any, int] | None:
    """(edited ``aa``, ``nt`` with codons moved alongside, new ``n_obs``) whose logit shift
    ``b * (log1p(n) - log1p(n0))`` is nearest ``shift`` (same sign), or None when none is
    reachable by moving K/R."""
    ks = [k for k in range(0, 2 * n0 + 2) if (k - n0) * shift > 0]
    ks.sort(key=lambda k: (abs(b * (math.log1p(k) - math.log1p(n0)) - shift), k))
    for k in ks:
        for order in _kr_orders(aa, k):
            new = "".join(aa[i] for i in order)
            if simple_covariates(new)[0] != k:
                continue
            if isinstance(nt, str) and len(nt) >= 3 * len(aa):
                nt = "".join(nt[3 * i : 3 * i + 3] for i in order) + nt[3 * len(aa) :]
            return new, nt, k
    return None


def _confound(
    layout: Layout, genes: pd.DataFrame, inside: np.ndarray, shift: float, b: float
) -> tuple[Layout, pd.DataFrame, np.ndarray]:
    """Carry ``shift`` through the members' ``n_obs`` (B1): (edited layout, genes with the new
    ``n_obs_true``, logit change per gene)."""
    delta = np.zeros(len(genes))
    if shift == 0.0 or not inside.any():
        return layout, genes, delta
    lg = layout.genes
    row = {str(t): i for i, t in enumerate(lg["locus_tag"])}
    aa, nt = lg["aa"].tolist(), lg["nt"].tolist()
    n_obs = genes["n_obs_true"].to_numpy(np.int64).copy()
    movable: list[int] = []
    for i in np.flatnonzero(inside).tolist():
        r = row[str(genes["locus_tag"].iat[i])]
        a, n0 = aa[r], int(n_obs[i])
        if not (isinstance(a, str) and a) or n0 <= 0:
            continue  # immovable: stays a member, unchanged
        movable.append(i)
        moved = _move_kr(a, nt[r], n0, shift, b)
        if moved is None:
            continue  # a short sORF with too few K/R: shift 0, still counted (decision 8)
        aa[r], nt[r], n_obs[i] = moved
        delta[i] = b * (math.log1p(n_obs[i]) - math.log1p(n0))
    if movable and abs(float(delta[movable].mean()) - shift) > SHIFT_TOL:
        raise ValueError(
            f"member_shift {shift} is not reachable by K/R moves: the mean over "
            f"{len(movable)} movable members is {float(delta[movable].mean()):.3f}"
        )
    lg = lg.assign(aa=pd.Series(aa, index=lg.index, dtype=object),
                   nt=pd.Series(nt, index=lg.index, dtype=lg["nt"].dtype))  # fmt: skip
    return replace(layout, genes=lg), genes.assign(n_obs_true=n_obs), delta


def _confounded_sets(sets: pd.DataFrame, n_genes: int, rng: np.random.Generator) -> list[str]:
    """One GO set of 8-20 % of the genes (else the one nearest 12 %), drawn at random."""
    go = sets[sets["collection"].str.startswith("go_")]
    size = go.groupby("set_id")["locus_tag"].size()
    ok = sorted(size.index[(size >= 0.08 * n_genes) & (size <= 0.2 * n_genes)])
    if not ok:
        if size.empty:
            return []
        ok = [str((size - 0.12 * n_genes).abs().idxmin())]
    return [ok[int(rng.integers(0, len(ok)))]]


def _packed(layout: Layout, n_sets: int, rng: np.random.Generator) -> list[list[str]]:
    """TU ids of each packed set: 1, 2, 3, 1, ... chr1 TUs of 10-30 genes, then
    up to ``PACKED_SINGLES`` single-gene chr1 TUs (at least 9, disjoint across sets)."""
    tus = layout.tus
    lo, hi = PACKED_TU_SIZE
    cand = [str(t) for t, c, gs in zip(tus["tu_id"], tus["contig"], tus["genes"], strict=True)
            if c == "chr1" and lo <= len(gs) <= hi]  # fmt: skip
    need = [1 + k % 3 for k in range(n_sets)]
    if len(cand) < sum(need):
        raise ValueError(f"layout has {len(cand)} packed TUs, {sum(need)} needed")
    single = [str(t) for t, c, gs in zip(tus["tu_id"], tus["contig"], tus["genes"], strict=True)
              if c == "chr1" and len(gs) == 1]  # fmt: skip
    per = min(PACKED_SINGLES, len(single) // n_sets)
    if per < 9:
        raise ValueError(f"layout has {len(single)} single-gene TUs, {9 * n_sets} needed")
    order = [cand[i] for i in rng.permutation(len(cand))]
    singles = [single[i] for i in rng.permutation(len(single))]
    out, at = [], 0
    for j, k in enumerate(need):
        out.append(order[at : at + k] + singles[j * per : (j + 1) * per])
        at += k
    return out


def _seeds(rng: np.random.Generator) -> list[np.random.Generator]:
    return [np.random.default_rng(int(s)) for s in rng.integers(0, 2**63 - 1, size=7)]


def draw(pair: tuple[LayoutConfig, ExpressionConfig], rng: np.random.Generator) -> World:
    """One world of a ``gippasim.scenarios`` function pair (a pure function of both)."""
    lcfg, ecfg = pair
    spec = getattr(ecfg, "function", None)
    if not isinstance(spec, FunctionSpec):
        raise ValueError("draw needs a gippasim.scenarios function pair (FunctionSpec missing)")
    r_lay, r_ont, r_pick, r_expr, r_det, r_feat, _ = _seeds(rng)
    layout = simulate_layout(lcfg, r_lay)
    tags = layout.genes["locus_tag"].astype(str).tolist()
    tu_of = _tu_of(layout)
    dag = _dag(_n_terms(len(tags)), 5, r_ont)
    go = _annotate(dag, tags, r_ont)
    cog_x, cog_def = _cog(tags, r_ont)
    links, kos = _kegg(len(tags), r_ont)
    packed: list[list[str]] = []
    if spec.kind == "packed":
        packed = _packed(layout, spec.n_sets, r_pick)
        kos = [f"K{_KO_BASE + k:05d}" for k in range(len(packed))]
        paths = [f"map{_MAP_BASE + k:05d}" for k in range(len(packed))]
        links = pd.DataFrame({"ko": [f"ko:{k}" for k in kos],
                              "pathway": [f"path:{p}" for p in paths]})  # fmt: skip
        tu_genes = {str(t): [str(g) for g in gs] for t, gs in
                    zip(layout.tus["tu_id"], layout.tus["genes"], strict=True)}  # fmt: skip
        kx = [(g, kos[k]) for k, ts in enumerate(packed) for t in ts for g in tu_genes[t]]
    else:
        on = np.flatnonzero(r_ont.random(len(tags)) < _P_KO)
        pick = r_ont.integers(0, len(kos), size=len(tags))
        kx = [(tags[i], kos[int(pick[i])]) for i in on.tolist()]
    kegg_x = pd.DataFrame({"locus_tag": [g for g, _ in kx], "db": "KEGG",
                           "id": [k for _, k in kx]}, columns=_XREF_COLUMNS)  # fmt: skip
    rank = {g: i for i, g in enumerate(tags)}
    dbr = {"COG": 0, "GO": 1, "KEGG": 2}
    xrefs = pd.concat([cog_x, go, kegg_x], ignore_index=True)
    xrefs = (xrefs.assign(_g=xrefs["locus_tag"].map(rank), _d=xrefs["db"].map(dbr))
             .sort_values(["_g", "_d", "id"], kind="mergesort").drop(columns=["_g", "_d"])
             .drop_duplicates().reset_index(drop=True))  # fmt: skip
    layout = replace(layout, xrefs=xrefs)
    paths_all = sorted({str(p).removeprefix("path:") for p in links["pathway"]})
    if spec.kind == "packed":  # the packed sets only, unfiltered
        rank_of = rank.__getitem__
        sets = pd.DataFrame(
            [(paths[k], "kegg", g, False) for k, ts in enumerate(packed)
             for g in sorted((g for t in ts for g in tu_genes[t]), key=rank_of)],
            columns=SET_COLUMNS,
        )  # fmt: skip
    else:
        sets = _collect(tags, tu_of, xrefs, dag, links, cog_def)

    expr = simulate_expression(layout, ecfg, r_expr)
    genes = expr.genes
    inside = np.zeros(len(tags), dtype=bool)
    if spec.kind == "ora" and spec.attack:
        conf = _confounded_sets(sets, len(tags), r_pick)
        flag = sets["set_id"].isin(conf).to_numpy()
        sets = sets.assign(confounded=flag)
        inside = genes["locus_tag"].isin(set(sets.loc[flag, "locus_tag"])).to_numpy()
    shift = spec.member_shift if spec.attack else 0.0
    layout, genes, delta = _confound(layout, genes, inside, shift, ecfg.beta["log1p_n_obs"])
    genes = _redetect(genes, delta, r_det)
    features: pd.DataFrame | None = None
    truth: dict[str, pd.DataFrame] = {"sets": sets}
    if spec.kind == "determinants":
        features = _features(layout, spec.features, r_feat)
        beta = np.array([spec.coefs.get(f, 0.0) for f in FEATURES])
        z = np.column_stack([_z(features[f]) for f in FEATURES])
        genes = genes.assign(abundance_log2_true=genes["abundance_log2_true"] + z @ beta)
        resid = _DET_TU_SD**2 + ecfg.gene_sd**2
        s2 = float((beta**2).sum())
        r2 = s2 / (s2 + resid)
        truth = {
            "determinants": pd.DataFrame(
                {
                    "feature": list(FEATURES),
                    "block": [BLOCKS[f] for f in FEATURES],
                    "beta": beta,
                    "true_r2": r2,
                    "residual_var": resid,
                }
            ),
            "determinants_genes": pd.DataFrame(
                {
                    "locus_tag": features["locus_tag"],
                    "dominant_block": np.array([BLOCKS[f] for f in FEATURES])[
                        np.abs(z * beta).argmax(axis=1)
                    ],
                }
            ),
        }
    letters = sorted(set(cog_def["cog_category"]))
    return World(
        layout=layout,
        expression=replace(expr, genes=genes),
        obo=_obo(dag),
        kegg_links=links,
        kegg_names=_names(paths_all),
        cog_def=cog_def,
        cog_names=pd.DataFrame({"category": letters,
                                "name": [f"synthetic COG category {c}" for c in letters]}),
        sets=sets,
        features=features,
        truth=MappingProxyType(truth),
    )  # fmt: skip


def write_truth(world: World, outdir: Path | str) -> dict[str, Path]:
    """Write every truth frame to ``outdir/truth/<name>.tsv``; returns ``{name: path}``."""
    d = Path(outdir) / "truth"
    d.mkdir(parents=True, exist_ok=True)
    out: dict[str, Path] = {}
    for name in sorted(world.truth):
        path = d / f"{name}.tsv"
        world.truth[name].to_csv(path, sep="\t", index=False, lineterminator="\n")
        out[name] = path
    return out
