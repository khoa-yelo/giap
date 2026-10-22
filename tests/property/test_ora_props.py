"""M8-04b properties of the TU-counted ORA (docs/specs/M8-04 acceptance criteria 3-6).

- TU membership never depends on ``detected``: flipping any detected value leaves every
  ``n_in_set`` (and ``n_in_set_genes``) unchanged;
- relabelling genes outside C, or permuting the labels of TUs outside C that have equal q, leaves
  C's row unchanged; splitting a C-free TU into one TU with the same q plus singletons of
  negligible p (1e-300, undetected) changes, among the TU columns, only ``n_universe``;
- BH families are exactly (collection, direction); a NaN p stays out (q NaN); ``p_classic`` and
  ``p_noelim`` never enter BH;
- a GO term none of whose (tested) descendants is elim-significant has ``p_value == p_noelim``.

Universes are small (<= 12 TUs) and use the raw q (no folds), so each example is cheap. Rows are
compared to 1e-10 relative: each side is a tilted evaluation within the decision-11 contract.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numpy as np
import pandas as pd
from helpers.ora_data import Gene, frame, gene_sets, ora
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.stats.multitest import bh

SETTINGS = settings(max_examples=40, deadline=None, suppress_health_check=list(HealthCheck))
ALPHA = 0.05
SIZE = (1, 500)
#: two tilted evaluations, each within the decision-11 contract (1e-11) of the exact law
RTOL = 1e-10
TU_COLS = ("n_in_set", "n_detected", "k_detected_in_set", "p_noelim", "p_value", "effect", "stat")

_p = st.floats(min_value=0.05, max_value=0.95, allow_nan=False)


@st.composite
def universes(draw: Any) -> tuple[dict[str, list[Gene]], list[str]]:
    n_tu = draw(st.integers(min_value=3, max_value=12))
    tus: dict[str, list[Gene]] = {}
    tags: list[str] = []
    for i in range(n_tu):
        size = draw(st.integers(min_value=1, max_value=3))
        genes = []
        for _ in range(size):
            tag = f"g{len(tags) + 1}"
            det = draw(st.booleans())
            genes.append(Gene(tag, 10.0 if det else None, det, p_detect=draw(_p)))
            tags.append(tag)
        tus[f"TU_{i + 1}"] = genes
    return tus, tags


@st.composite
def gene_subsets(draw: Any, tags: list[str], n: int) -> dict[str, set[str]]:
    out = {}
    for j in range(n):
        out[f"T{j}"] = set(draw(st.lists(st.sampled_from(tags), min_size=1, max_size=6)))
    return out


def _run(genes: pd.DataFrame, colls: list[Any]) -> pd.DataFrame:
    return ora().detection_enrichment(genes, colls, SIZE, ALPHA, folds=None)


def _rows(out: pd.DataFrame, term: str) -> pd.DataFrame:
    return out[out["term"] == term].sort_values("direction").reset_index(drop=True)


def _same(a: pd.DataFrame, b: pd.DataFrame, cols: tuple[str, ...]) -> None:
    for c in cols:
        x, y = a[c].to_numpy(float), b[c].to_numpy(float)
        atol = 1e-10 if c == "effect" else 0.0  # k - E[K | N]: a difference of O(1) numbers
        np.testing.assert_allclose(x, y, rtol=RTOL, atol=atol, equal_nan=True, err_msg=c)


@SETTINGS
@given(data=st.data())
def test_membership_never_depends_on_detection(data: Any) -> None:
    tus, tags = data.draw(universes())
    sets = data.draw(gene_subsets(tags, 3))
    colls = [gene_sets("kegg", sets)]
    out = _run(frame(tus), colls)
    flip = data.draw(st.lists(st.sampled_from(tags), min_size=1, unique=True))
    g2 = frame(tus)
    sel = g2["locus_tag"].isin(flip)
    g2.loc[sel, "detected"] = ~g2.loc[sel, "detected"].astype(bool)
    out2 = _run(g2, colls)
    key = ["collection", "term", "direction"]
    a = out.sort_values(key).reset_index(drop=True)
    b = out2.sort_values(key).reset_index(drop=True)
    assert list(a["term"]) == list(b["term"])
    np.testing.assert_array_equal(a["n_in_set"], b["n_in_set"])
    np.testing.assert_array_equal(a["n_in_set_genes"], b["n_in_set_genes"])
    np.testing.assert_array_equal(a["n_universe"], b["n_universe"])


@SETTINGS
@given(data=st.data())
def test_relabelling_genes_outside_c_leaves_its_row(data: Any) -> None:
    tus, tags = data.draw(universes())
    c = set(data.draw(st.lists(st.sampled_from(tags), min_size=1, max_size=5)))
    other = data.draw(gene_subsets(tags, 2))
    out = _run(frame(tus), [gene_sets("kegg", {"C": c, **other})])
    ren = {t: (t if t in c else f"z{t}") for t in tags}
    tus2 = {tu: [Gene(ren[x.tag], x.value, x.detected, p_detect=x.p_detect) for x in gs]
            for tu, gs in tus.items()}  # fmt: skip
    other2 = {k: {ren[t] for t in v} for k, v in other.items()}
    out2 = _run(frame(tus2), [gene_sets("kegg", {"C": c, **other2})])
    a, b = _rows(out, "C"), _rows(out2, "C")
    _same(a, b, (*TU_COLS, "odds_weight", "p_classic", "q_value", "n_universe"))


@SETTINGS
@given(
    n_out=st.integers(min_value=3, max_value=10),
    p_out=_p,
    det_out=st.lists(st.booleans(), min_size=10, max_size=10),
    c_tus=st.lists(st.tuples(st.integers(1, 3), _p, st.booleans()), min_size=1, max_size=4),
    perm_seed=st.integers(0, 2**31 - 1),
)
def test_permuting_equal_q_tus_outside_c_leaves_its_row(
    n_out: int, p_out: float, det_out: list[bool], c_tus: list[tuple[int, float, bool]],
    perm_seed: int,
) -> None:  # fmt: skip
    tus: dict[str, list[Gene]] = {}
    c: set[str] = set()
    for i, (size, p, det) in enumerate(c_tus):
        genes = [Gene(f"c{i}_{j}", None, det and j == 0, p_detect=p) for j in range(size)]
        tus[f"TU_{i + 1}"] = genes
        c |= {x.tag for x in genes}
    base = len(c_tus)
    out_ids = [f"TU_{base + j + 1}" for j in range(n_out)]
    for j, tu in enumerate(out_ids):
        tus[tu] = [Gene(f"o{j}", None, det_out[j], p_detect=p_out)]
    out = _run(frame(tus), [gene_sets("kegg", {"C": c})])
    perm = np.random.default_rng(perm_seed).permutation(n_out)
    tus2 = {k: v for k, v in tus.items() if k not in out_ids}
    for j, tu in enumerate(out_ids):
        tus2[out_ids[int(perm[j])]] = tus[tu]
    out2 = _run(frame(tus2), [gene_sets("kegg", {"C": c})])
    _same(_rows(out, "C"), _rows(out2, "C"), (*TU_COLS, "odds_weight", "p_classic", "n_universe"))


@SETTINGS
@given(data=st.data())
def test_splitting_a_c_free_tu_changes_only_n_universe(data: Any) -> None:
    tus, _ = data.draw(universes())
    c_tu = data.draw(st.sampled_from(sorted(tus)))
    c = {x.tag for x in tus[c_tu]}
    free = [t for t in sorted(tus) if t != c_tu and len(tus[t]) >= 2]
    if not free:
        tus["TU_99"] = [Gene("s1", None, False, p_detect=0.3), Gene("s2", 10.0, True,
                                                                    p_detect=0.4)]  # fmt: skip
        free = ["TU_99"]
    t = data.draw(st.sampled_from(free))
    out = _run(frame(tus), [gene_sets("kegg", {"C": c})])
    genes = tus[t]
    q = -math.expm1(sum(math.log1p(-x.p_detect) for x in genes))
    det = any(x.detected for x in genes)
    tus2 = {k: v for k, v in tus.items() if k != t}
    tus2[t] = [Gene(genes[0].tag, 10.0 if det else None, det, p_detect=q)]
    for j, x in enumerate(genes[1:]):
        tus2[f"TU_{1000 + j}"] = [Gene(x.tag, None, False, p_detect=1e-300)]
    out2 = _run(frame(tus2), [gene_sets("kegg", {"C": c})])
    a, b = _rows(out, "C"), _rows(out2, "C")
    _same(a, b, TU_COLS)
    np.testing.assert_array_equal(b["n_universe"], a["n_universe"] + len(genes) - 1)


@st.composite
def go_dags(draw: Any, tags: list[str]) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """A random DAG over 2-6 terms; each term's genes include its children's (propagated)."""
    n = draw(st.integers(min_value=2, max_value=6))
    names = [f"GO:{i}" for i in range(n)]
    parents: dict[str, set[str]] = {}
    for i, t in enumerate(names):
        parents[t] = set(draw(st.lists(st.sampled_from(names[:i]), max_size=2))) if i else set()
    own = {t: set(draw(st.lists(st.sampled_from(tags), min_size=1, max_size=4))) for t in names}
    sets = {t: set(own[t]) for t in names}
    for t in reversed(names):  # children have larger indices: push genes up
        for p in parents[t]:
            sets[p] |= sets[t]
    return sets, parents


def _descendants(parents: Mapping[str, set[str]]) -> dict[str, set[str]]:
    kids: dict[str, set[str]] = {t: set() for t in parents}
    for t, ps in parents.items():
        for p in ps:
            kids[p].add(t)

    def walk(t: str) -> set[str]:
        out: set[str] = set()
        for k in kids[t]:
            out |= {k} | walk(k)
        return out

    return {t: walk(t) for t in parents}


@SETTINGS
@given(data=st.data())
def test_bh_families_and_elim_identity(data: Any) -> None:
    tus, tags = data.draw(universes())
    sets, parents = data.draw(go_dags(tags))
    kegg_sets = data.draw(gene_subsets(tags, 3))
    colls = [gene_sets("go_bp", sets, parents=parents), gene_sets("kegg", kegg_sets)]
    out = _run(frame(tus), colls)
    assert set(out["direction"]) <= {"over", "under"}
    for (_, _), fam in out.groupby(["collection", "direction"]):
        p = fam["p_value"].to_numpy(float)
        q = fam["q_value"].to_numpy(float)
        np.testing.assert_array_equal(np.isnan(q), np.isnan(p))
        np.testing.assert_allclose(q, bh(p), rtol=1e-15, equal_nan=True)
    go_over = out[(out["collection"] == "go_bp") & (out["direction"] == "over")]
    elim_sig = set(go_over.loc[go_over["p_value"] < 0.01, "term"])
    desc = _descendants(parents)
    for _, r in go_over.iterrows():
        if not desc[r["term"]] & elim_sig:
            assert r["p_value"] == r["p_noelim"] or (
                math.isnan(r["p_value"]) and math.isnan(r["p_noelim"])
            ), r["term"]
    rest = out[~((out["collection"] == "go_bp") & (out["direction"] == "over"))]
    np.testing.assert_array_equal(rest["p_value"].to_numpy(float), rest["p_noelim"].to_numpy(float))
