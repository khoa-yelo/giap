"""M8-05 properties of the TU-collapsed GSEA (docs/specs/M8-05-gsea-tu-collapsed.md).

Acceptance criteria: permuting gene rows changes nothing; relabelling TUs changes nothing;
adding an undetected gene to a set changes nothing; every leading-edge gene is a detected member
of its set; set sizes are counted in TUs; the redundancy re-test p replaces the parent's p in its
collection's BH vector (lead decision 3).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.gsea_data import Gene, as_list, frame, gene_sets
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.config.model import GseaConfig
from giap.rng import for_step
from giap.stats.fgsea import SimpleResult
from giap.stats.multitest import bh

SETTINGS = settings(max_examples=25, deadline=None, suppress_health_check=list(HealthCheck))
ALPHA = 0.05
N_PERM = 60


def _g() -> Any:
    from giap.function import gsea

    return gsea


def _run(genes: pd.DataFrame, colls: list[Any], seed: int = 3, size: tuple[int, int] = (2, 500)
         ) -> pd.DataFrame:  # fmt: skip
    cfg = GseaConfig(n_perm_simple=N_PERM)
    return _g().gsea_table(genes, colls, size, cfg, ALPHA, for_step("function.gsea", seed))


def _canon(out: pd.DataFrame) -> pd.DataFrame:
    df = out.sort_values(["collection", "term"]).reset_index(drop=True).copy()
    for c in ("leading_edge", "flags"):
        df[c] = df[c].map(lambda v: tuple(as_list(v)))
    return df


@st.composite
def worlds(draw: st.DrawFn, ties: bool = True) -> tuple[pd.DataFrame, dict[str, set[str]]]:
    """A gene frame of 8-30 TUs (1-3 genes, some undetected) and 1-4 random gene sets."""
    n_tu = draw(st.integers(8, 30))
    seed = draw(st.integers(0, 2**31 - 1))
    rng = np.random.default_rng(seed)
    tus: dict[str, list[Gene]] = {}
    k = 0
    for i in range(n_tu):
        genes = []
        for _ in range(int(rng.integers(1, 4))):
            k += 1
            det = bool(rng.random() < 0.8) or not genes
            v = float(np.round(rng.normal(20, 2), 0 if ties else 6)) + (0 if ties else 1e-6 * k)
            genes.append(
                Gene(f"g{k}", v if det else None, detected=det, imputed=bool(rng.random() < 0.2))
            )
        tus[f"TU_{i + 1}"] = genes
    df = frame(tus)
    tags = list(df["locus_tag"])
    sets = {
        f"S{j}": {
            tags[int(x)]
            for x in rng.choice(len(tags), int(rng.integers(2, len(tags) // 2 + 2)), replace=False)
        }
        for j in range(int(rng.integers(1, 5)))
    }
    return df, sets


def _colls(sets: Mapping[str, set[str]]) -> list[Any]:
    return [gene_sets("cog", sets), gene_sets("kegg", {f"K{k}": v for k, v in sets.items()})]


@SETTINGS
@given(worlds(), st.integers(0, 2**31 - 1))
def test_permuting_gene_rows_changes_nothing(world: Any, perm_seed: int) -> None:
    genes, sets = world
    shuffled = genes.sample(frac=1.0, random_state=perm_seed % (2**32)).reset_index(drop=True)
    a, b = _canon(_run(genes, _colls(sets))), _canon(_run(shuffled, _colls(sets)))
    pd.testing.assert_frame_equal(a, b)


@SETTINGS
@given(worlds(), st.integers(0, 2**31 - 1))
def test_tu_ranking_is_independent_of_row_order(world: Any, perm_seed: int) -> None:
    genes, _ = world
    g = _g()
    shuffled = genes.sample(frac=1.0, random_state=perm_seed % (2**32)).reset_index(drop=True)
    a, b = g.tu_ranking(genes), g.tu_ranking(shuffled)
    assert tuple(a.tu_ids) == tuple(b.tu_ids)
    np.testing.assert_array_equal(np.asarray(a.scores), np.asarray(b.scores))
    assert {k: tuple(v) for k, v in a.tu_genes.items()} == {
        k: tuple(v) for k, v in b.tu_genes.items()
    }


@SETTINGS
@given(worlds(), st.integers(0, 2**31 - 1))
def test_scores_are_centred_tu_means(world: Any, seed: int) -> None:
    genes, _ = world
    rank = _g().tu_ranking(genes)
    det = genes[genes["detected"] & genes["abundance_log2"].notna()]
    means = det.groupby("tu_id")["abundance_log2"].mean()
    assert set(rank.tu_ids) == set(means.index)
    expect = means - float(np.median(means.to_numpy()))
    got = dict(zip(rank.tu_ids, np.asarray(rank.scores, dtype=float), strict=True))
    for t, v in expect.items():
        assert got[t] == pytest.approx(float(v), abs=1e-9)
    assert float(np.median(np.asarray(rank.scores, dtype=float))) == pytest.approx(0.0, abs=1e-9)


@SETTINGS
@given(worlds(ties=False), st.integers(0, 2**31 - 1))
def test_relabelling_tus_changes_nothing(world: Any, perm_seed: int) -> None:
    """Without tied scores, renaming the TUs by a bijection leaves every statistic unchanged."""
    genes, sets = world
    ids = sorted(set(genes["tu_id"]))
    new = np.random.default_rng(perm_seed).permutation(len(ids))
    rename = {t: f"TU_{int(j) + 1}" for t, j in zip(ids, new, strict=True)}
    renamed = genes.assign(tu_id=genes["tu_id"].map(rename))
    a, b = _canon(_run(genes, _colls(sets))), _canon(_run(renamed, _colls(sets)))
    pd.testing.assert_frame_equal(a, b)


@SETTINGS
@given(worlds(), st.integers(0, 3))
def test_adding_an_undetected_gene_to_a_set_changes_nothing(world: Any, which: int) -> None:
    genes, sets = world
    tu = str(genes["tu_id"].iloc[which % len(genes)])
    extra = frame({tu: [Gene("u_new", None, detected=False)]}).assign(tu_position=99)
    more = pd.concat([genes, extra], ignore_index=True)
    key = sorted(sets)[which % len(sets)]
    bigger = {k: (v | {"u_new"} if k == key else v) for k, v in sets.items()}
    a, b = _canon(_run(genes, _colls(sets))), _canon(_run(more, _colls(bigger)))
    pd.testing.assert_frame_equal(a, b)


@SETTINGS
@given(worlds())
def test_leading_edge_genes_are_detected_members_and_sizes_count_tus(world: Any) -> None:
    genes, sets = world
    g = _g()
    rank = g.tu_ranking(genes)
    coll = gene_sets("cog", sets)
    tus = g.tu_sets(rank, coll)
    detected = set(genes.loc[genes["detected"] & genes["abundance_log2"].notna(), "locus_tag"])
    tu_of = dict(zip(genes["locus_tag"], genes["tu_id"], strict=True))
    out = _run(genes, [coll])
    for _, r in out.iterrows():
        term = r["term"]
        assert r["size_tus"] == len(tus[term])
        assert tus[term] == {tu_of[x] for x in sets[term] if x in detected}
        edge = as_list(r["leading_edge"])
        assert len(edge) == len(set(edge))
        for x in edge:
            assert x in sets[term] and x in detected, (term, x)
        assert {tu_of[x] for x in edge} <= tus[term]
        if r["es"] == 0:
            assert edge == []


@SETTINGS
@given(worlds())
def test_q_values_are_bh_of_each_collection(world: Any) -> None:
    genes, sets = world
    out = _run(genes, _colls(sets))
    for _, part in out.groupby("collection"):
        p = part["p_value"].to_numpy(float)
        np.testing.assert_allclose(part["q_value"].to_numpy(float), bh(p), atol=1e-15)


# --------------------------------------------------------------------------- BH substitution


class _Fake:
    """``fgsea_simple`` with p keyed by set size in TUs; records the sizes of each call."""

    def __init__(self, p_by_size: Mapping[int, float]) -> None:
        self.p = dict(p_by_size)
        self.calls: list[list[int]] = []

    def __call__(self, prep: Any, n_perm: int, sample_size: int, rng: Any) -> SimpleResult:
        sizes = [int(s) for s in np.asarray(prep.sizes)]
        self.calls.append(sizes)
        k = len(sizes)
        p = np.array([self.p[s] for s in sizes], dtype=float)
        edges = tuple(
            np.asarray(prep.order)[np.asarray(m)[:1] - 1].astype(np.int64) for m in prep.members
        )
        return SimpleResult(
            np.full(k, 0.5), np.full(k, 1.5), p, np.full(k, 0.1), np.zeros(k, dtype=np.int64),
            np.full(k, 500, dtype=np.int64), np.zeros(k, dtype=bool), edges,
        )  # fmt: skip


p_values = st.floats(1e-4, 1.0, allow_nan=False)


@settings(max_examples=150, deadline=None, suppress_health_check=list(HealthCheck))
@given(
    p_child=p_values,
    p_parent=p_values,
    p_retest=p_values,
    p_null=st.lists(p_values, min_size=0, max_size=8),
)
def test_redundant_iff_bh_with_the_retest_p_substituted_exceeds_alpha(
    monkeypatch: pytest.MonkeyPatch,
    p_child: float,
    p_parent: float,
    p_retest: float,
    p_null: list[float],
) -> None:
    """C = 12 TUs, P = C + 15 TUs (re-test size 15), nulls of sizes 40, 41, ...: the parent is
    re-tested only when both are significant, and is redundant exactly when the BH q of its
    collection's main p-vector with P's p replaced by the re-test p exceeds alpha."""
    import giap.stats.fgsea as fg

    n = 400
    genes = frame({f"TU_{i + 1}": [Gene(f"g{i + 1}", 20.0 + 0.01 * i)] for i in range(n)})
    tags = list(genes["locus_tag"])
    sets = {"GO:C": tags[:12], "GO:P": tags[:27]}
    parents: dict[str, set[str]] = {"GO:C": {"GO:P"}, "GO:P": set()}
    table = {12: p_child, 27: p_parent, 15: p_retest}
    for i, p in enumerate(p_null):
        sets[f"GO:N{i}"] = tags[n - 40 - i :]
        parents[f"GO:N{i}"] = set()
        table[40 + i] = p
    fake = _Fake(table)
    monkeypatch.setattr(fg, "fgsea_simple", fake)
    mod = _g()
    if hasattr(mod, "fgsea_simple"):
        monkeypatch.setattr(mod, "fgsea_simple", fake)
    out = _run(genes, [gene_sets("go_bp", sets, parents=parents)], size=(10, 500))

    main = np.array([p_child, p_parent, *p_null])
    q = bh(main)
    by = {r["term"]: r for _, r in out.iterrows()}
    assert by["GO:P"]["q_value"] == pytest.approx(q[1], abs=1e-15)
    retested = q[0] <= ALPHA and q[1] <= ALPHA
    assert (fake.calls[-1] == [15]) == retested
    sub = main.copy()
    sub[1] = p_retest
    expect = retested and bool(bh(sub)[1] > ALPHA)
    assert bool(by["GO:P"]["redundant"]) == expect
    assert not bool(by["GO:C"]["redundant"])
    assert not any(bool(by[f"GO:N{i}"]["redundant"]) for i in range(len(p_null)))
    assert all(not math.isnan(float(r["q_value"])) for r in by.values())
