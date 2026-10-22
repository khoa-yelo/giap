"""M8-05: abundance GSEA with a TU-collapsed ranking (docs/specs/M8-05-gsea-tu-collapsed.md).

Examples 1-8, every degenerate row and the column contract of ``giap.function.gsea``:
``tu_ranking(genes) -> TuRanking``, ``tu_sets(rank, coll) -> {term: frozenset[TU]}`` and
``gsea_table(genes, colls, size, cfg, alpha, rng) -> DataFrame``.

Resolutions pinned here (the spec leaves them open):

* Tie order. Lead decision 4 (M8-02 lead decision 7) overrides rule 2 and example 4: TUs enter
  fgsea in a seeded random order drawn from ``rng`` (stream ``function.gsea``), never by TU ID.
  ``TuRanking.tu_ids`` itself stays in natural TU-ID order (TU_9 before TU_10), so the ranking does
  not depend on the row order of the gene table.
* ``TuRanking.tu_genes`` holds the genes that score the TU (detected, eligible, non-NA
  ``abundance_log2``) in ``tu_position`` order; only ranked TUs carry genes.
* Run flags (``layer_missing:abundance``, ``too_few_ranked_tus``) are returned in
  ``table.attrs["run_flags"]``, a list of strings; an empty table keeps every column.
* List cells (``leading_edge``, ``flags``) are sequences of strings; the writer joins them.
  Every row's ``flags`` carries ``ci_unestimated`` (lead decision 4).
* ``q_value`` is the BH q of the main-run p within the collection. The redundancy re-test changes
  only ``redundant`` (decision 10); its own q is not a column.
* Tests that need exact control of p-values replace ``fgsea_simple`` (rule 4 says it is the
  engine) by a fake keyed on set size in TUs.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.gsea_data import (
    GSEA_COLUMNS,
    Gene,
    as_list,
    frame,
    gene_sets,
    singleton_frame,
)

from giap.config.model import GseaConfig
from giap.rng import for_step
from giap.stats.fgsea import SimpleResult, gsea_es, prepare
from giap.stats.multitest import bh

STREAM = "function.gsea"
ALPHA = 0.05


def _g() -> Any:
    from giap.function import gsea

    return gsea


def _run(
    genes: pd.DataFrame,
    colls: Sequence[Any],
    size: tuple[int, int] = (10, 500),
    n_perm: int = 1000,
    seed: int = 1,
    alpha: float = ALPHA,
) -> pd.DataFrame:
    cfg = GseaConfig(n_perm_simple=n_perm)
    return _g().gsea_table(genes, colls, size, cfg, alpha, for_step(STREAM, seed))


def _row(out: pd.DataFrame, term: str, collection: str | None = None) -> pd.Series:
    sel = out["term"] == term
    if collection is not None:
        sel &= out["collection"] == collection
    rows = out[sel]
    assert len(rows) == 1, f"expected one row for {term}, got {len(rows)}"
    return rows.iloc[0]


def _na(v: Any) -> bool:
    return v is None or (isinstance(v, float) and math.isnan(v)) or bool(pd.isna(v))


def _flags(v: Any) -> list[str]:
    return as_list(v)


def _ranked_tus(genes: pd.DataFrame) -> list[str]:
    """Singleton frame helper: TU IDs from highest to lowest abundance."""
    return list(genes.sort_values("abundance_log2", ascending=False)["tu_id"])


def _tag(genes: pd.DataFrame, tus: Sequence[str]) -> list[str]:
    by = dict(zip(genes["tu_id"], genes["locus_tag"], strict=True))
    return [by[t] for t in tus]


def _distinct_frame(n: int, seed: int = 3) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    return singleton_frame(20.0 + 0.01 * rng.permutation(n))


# --------------------------------------------------------------------------- tu_ranking


def example_1() -> pd.DataFrame:
    return frame(
        {
            "TU_1": [Gene("g1", 18.0), Gene("g2", 16.0)],
            "TU_2": [Gene("g3", 15.0)],
            "TU_3": [Gene("g4", None, detected=False)],
        }
    )


def test_example_1_tu_scores_are_median_centred_means() -> None:
    rank = _g().tu_ranking(example_1())
    assert tuple(rank.tu_ids) == ("TU_1", "TU_2")
    np.testing.assert_allclose(np.asarray(rank.scores, dtype=float), [1.0, -1.0], atol=1e-12)
    assert np.asarray(rank.scores).dtype == np.float64
    assert tuple(rank.tu_genes["TU_1"]) == ("g1", "g2")
    assert tuple(rank.tu_genes["TU_2"]) == ("g3",)
    assert "TU_3" not in rank.tu_ids
    assert tuple(rank.tu_genes.get("TU_3", ())) == ()


def test_odd_number_of_tus_centres_on_the_middle_score() -> None:
    genes = singleton_frame([10.0, 13.0, 11.0])
    rank = _g().tu_ranking(genes)
    by = dict(zip(rank.tu_ids, np.asarray(rank.scores, dtype=float), strict=True))
    assert by == pytest.approx({"TU_1": -1.0, "TU_2": 2.0, "TU_3": 0.0}, abs=1e-12)


def test_imputed_abundances_count_in_the_tu_mean() -> None:
    """Lead decision 2 (spec decision 9): imputed values are included."""
    genes = frame(
        {
            "TU_1": [Gene("a", 20.0), Gene("b", 10.0, imputed=True)],
            "TU_2": [Gene("c", 14.0)],
            "TU_3": [Gene("d", 16.0, imputed=True)],
        }
    )
    rank = _g().tu_ranking(genes)
    by = dict(zip(rank.tu_ids, np.asarray(rank.scores, dtype=float), strict=True))
    # means 15, 14, 16 -> median 15
    assert by == pytest.approx({"TU_1": 0.0, "TU_2": -1.0, "TU_3": 1.0}, abs=1e-12)
    assert tuple(rank.tu_genes["TU_1"]) == ("a", "b")


def test_undetected_ineligible_and_na_genes_do_not_score() -> None:
    genes = frame(
        {
            "TU_1": [
                Gene("a", 12.0),
                Gene("u", 40.0, detected=False),  # undetected: never scored, whatever its value
                Gene("p", 40.0, partial=True),  # partial: not eligible
                Gene("n", 40.0, p_detect=math.nan),  # no p_detect: not eligible
                Gene("x", None),  # detected without a value
            ],
            "TU_2": [Gene("b", 10.0)],
        }
    )
    rank = _g().tu_ranking(genes)
    by = dict(zip(rank.tu_ids, np.asarray(rank.scores, dtype=float), strict=True))
    assert by == pytest.approx({"TU_1": 1.0, "TU_2": -1.0}, abs=1e-12)
    assert tuple(rank.tu_genes["TU_1"]) == ("a",)


def test_example_4_tu_ids_in_natural_order() -> None:
    """Natural sort of TU_<n>: TU_2 < TU_9 < TU_10, whatever the row order."""
    genes = frame(
        {
            "TU_10": [Gene("a", 5.0, imputed=True)],
            "TU_9": [Gene("b", 5.0, imputed=True)],
            "TU_2": [Gene("c", 7.0)],
        }
    )
    assert tuple(_g().tu_ranking(genes).tu_ids) == ("TU_2", "TU_9", "TU_10")


def test_genes_within_a_tu_follow_tu_position_not_locus_tag() -> None:
    genes = frame({"TU_1": [Gene("z1", 3.0), Gene("a1", 5.0)], "TU_2": [Gene("b", 1.0)]})
    genes = genes.iloc[::-1].reset_index(drop=True)
    assert tuple(_g().tu_ranking(genes).tu_genes["TU_1"]) == ("z1", "a1")


def test_tu_ranking_does_not_mutate_its_input() -> None:
    genes = example_1()
    before = genes.copy(deep=True)
    _g().tu_ranking(genes)
    pd.testing.assert_frame_equal(genes, before)


# --------------------------------------------------------------------------- tu_sets


def test_example_2_a_tu_is_in_a_set_through_any_detected_gene() -> None:
    g = _g()
    rank = g.tu_ranking(example_1())
    out = g.tu_sets(rank, gene_sets("cog", {"S": {"g2"}}))
    assert out["S"] == frozenset({"TU_1"})


def test_example_3_an_undetected_member_puts_no_tu_in_the_set() -> None:
    g = _g()
    rank = g.tu_ranking(example_1())
    out = g.tu_sets(rank, gene_sets("cog", {"S": {"g4"}}))
    assert out.get("S", frozenset()) == frozenset()


def test_set_size_counts_tus_not_genes() -> None:
    g = _g()
    rank = g.tu_ranking(example_1())
    out = g.tu_sets(rank, gene_sets("cog", {"S": {"g1", "g2", "g3", "g4"}}))
    assert out["S"] == frozenset({"TU_1", "TU_2"})


# --------------------------------------------------------------------------- gsea_table


def _strong(n: int = 200, k: int = 15) -> tuple[pd.DataFrame, list[str]]:
    genes = _distinct_frame(n)
    return genes, _ranked_tus(genes)


def test_columns_and_constant_fields() -> None:
    genes, order = _strong()
    top = _tag(genes, order[:15])
    mid = _tag(genes, order[90:105])
    coll = gene_sets("kegg", {"map00010": top, "map00020": mid}, coverage=0.4)
    out = _run(genes, [coll])
    assert set(GSEA_COLUMNS) <= set(out.columns), set(GSEA_COLUMNS) - set(out.columns)
    assert list(out.columns)[-1] == "flags"  # ADR-A: flags is the last column
    assert len(out) == 2
    for _, r in out.iterrows():
        assert r["collection"] == "kegg"
        assert r["name"] == f"name of {r['term']}"
        assert r["size_tus"] == 15
        assert r["method"] == "fgsea_simple"
        assert r["null_model"] == "tu_set_permutation"
        assert r["n_perm"] == 1000
        assert r["stat"] == r["es"]
        assert r["effect"] == r["nes"] or (_na(r["effect"]) and _na(r["nes"]))
        assert _na(r["ci_low"]) and _na(r["ci_high"])
        assert r["coverage"] == pytest.approx(0.4)
        assert r["redundant"] is False or r["redundant"] == False  # noqa: E712
        assert _na(r["parent_group"])  # non-GO
        assert "ci_unestimated" in _flags(r["flags"])
        assert 1 / 1001 <= r["p_value"] <= 1
    t = _row(out, "map00010")
    assert t["es"] > 0 and t["nes"] > 0
    assert t["p_value"] <= 0.01


def test_n_perm_follows_the_config() -> None:
    genes, order = _strong()
    out = _run(genes, [gene_sets("cog", {"K": _tag(genes, order[:15])})], n_perm=200)
    r = _row(out, "K")
    assert r["n_perm"] == 200
    assert r["p_value"] >= 1 / 201


def test_es_matches_the_fgsea_port_on_the_centred_tu_scores() -> None:
    """With distinct scores the tie order is irrelevant, so ES is that of M8-02's ``gsea_es`` on
    the median-centred TU means."""
    genes, order = _strong()
    s = genes["abundance_log2"].to_numpy(float)
    centred = s - np.median(s)
    labels = list(genes["tu_id"])
    rng = np.random.default_rng(5)
    sets = {f"S{i}": [order[j] for j in rng.choice(200, 12 + i, replace=False)] for i in range(5)}
    prep = prepare(centred, labels, sets, 10, 500)
    expect = {n: gsea_es(prep.stats, m)[0] for n, m in zip(prep.names, prep.members, strict=True)}
    coll = gene_sets("go_mf", {k: _tag(genes, v) for k, v in sets.items()}, parents={})
    out = _run(genes, [coll])
    for name, es in expect.items():
        r = _row(out, name)
        assert r["es"] == pytest.approx(es, abs=1e-12)
        assert r["size_tus"] == len(sets[name])


def test_q_value_is_bh_within_each_collection() -> None:
    genes, order = _strong(300)
    rng = np.random.default_rng(9)
    colls = []
    for c, top in (("cog", True), ("kegg", False), ("go_slim", True)):
        sets = {
            f"{c}_{i}": _tag(genes, [order[j] for j in rng.choice(300, 12, replace=False)])
            for i in range(6)
        }
        if top:
            sets[f"{c}_top"] = _tag(genes, order[:12])
        colls.append(gene_sets(c, sets))
    out = _run(genes, colls)
    assert set(out["collection"]) == {"cog", "kegg", "go_slim"}
    for _c, part in out.groupby("collection"):
        np.testing.assert_allclose(
            part["q_value"].to_numpy(float), bh(part["p_value"].to_numpy(float)), atol=1e-15
        )
    # one family per collection, not one over the run
    pooled = bh(out["p_value"].to_numpy(float))
    assert not np.allclose(pooled, out["q_value"].to_numpy(float))


def test_example_5_size_filter_in_tus() -> None:
    genes, order = _strong(700)
    coll = gene_sets(
        "cog",
        {
            "big": _tag(genes, order[:600]),
            "max": _tag(genes, order[100:600]),
            "ok": _tag(genes, order[:50]),
            "ten": _tag(genes, order[200:210]),
            "nine": _tag(genes, order[300:309]),
        },
    )
    out = _run(genes, [coll], n_perm=100)
    assert set(out["term"]) == {"max", "ok", "ten"}
    assert _row(out, "max")["size_tus"] == 500
    assert _row(out, "ten")["size_tus"] == 10


def _operon_frame(n_single: int = 200, n_op: int = 30) -> pd.DataFrame:
    rng = np.random.default_rng(17)
    tus = {
        f"TU_{i + 1}": [Gene(f"s{i + 1}", 20.0 + 0.01 * float(v))]
        for i, v in enumerate(rng.permutation(n_single))
    }
    tus[f"TU_{n_single + 1}"] = [Gene(f"op{j + 1}", 25.0 + 0.001 * j) for j in range(n_op)]
    return frame(tus)


def test_example_7_an_operon_packed_set_is_one_tu_and_dropped() -> None:
    genes = _operon_frame()
    packed = [f"op{j + 1}" for j in range(30)]
    spread = [f"s{i + 1}" for i in range(30)]
    out = _run(genes, [gene_sets("cog", {"packed": packed, "spread": spread})])
    assert set(out["term"]) == {"spread"}
    assert _row(out, "spread")["size_tus"] == 30
    assert _g().tu_sets(_g().tu_ranking(genes), gene_sets("cog", {"packed": packed}))[
        "packed"
    ] == frozenset({"TU_201"})
    # the gene-level ranking would test the packed set at 30 genes (why GIAP ranks TUs)
    vals = genes["abundance_log2"].to_numpy(float)
    gene_prep = prepare(
        vals - np.median(vals), list(genes["locus_tag"]), {"packed": packed}, 10, 500
    )
    assert list(gene_prep.names) == ["packed"] and list(gene_prep.sizes) == [30]


def test_example_8_same_seed_same_table() -> None:
    genes, order = _strong()
    rng = np.random.default_rng(2)
    coll = gene_sets(
        "go_bp",
        {
            f"T{i}": _tag(genes, [order[j] for j in rng.choice(200, 15, replace=False)])
            for i in range(8)
        },
        parents={},
    )
    a = _run(genes, [coll], seed=7)
    b = _run(genes, [coll], seed=7)
    pd.testing.assert_frame_equal(a, b)


def test_example_8_identical_across_numba_thread_counts() -> None:
    numba = pytest.importorskip("numba")
    max_threads = int(numba.config.NUMBA_NUM_THREADS)
    if max_threads < 2:
        pytest.skip("needs NUMBA_NUM_THREADS >= 2")
    genes, order = _strong()
    rng = np.random.default_rng(2)
    coll = gene_sets(
        "kegg",
        {
            f"T{i}": _tag(genes, [order[j] for j in rng.choice(200, 15, replace=False)])
            for i in range(8)
        },
    )
    before = numba.get_num_threads()
    try:
        numba.set_num_threads(1)
        one = _run(genes, [coll], seed=7)
        numba.set_num_threads(min(8, max_threads))
        many = _run(genes, [coll], seed=7)
    finally:
        numba.set_num_threads(before)
    pd.testing.assert_frame_equal(one, many)


def test_tied_tus_enter_fgsea_in_a_seeded_random_order() -> None:
    """Lead decision 4 / M8-02 decision 7: ties are broken by a seeded random order, never by TU
    ID. A set of the ten lowest-ID tied TUs would have one fixed ES under TU-ID order; under a
    random order its ES varies with the seed."""
    tus = {f"TU_{i + 1}": [Gene(f"d{i + 1}", 21.0 + 0.01 * i)] for i in range(30)}
    tus |= {f"TU_{i + 1}": [Gene(f"t{i + 1}", 20.0, imputed=True)] for i in range(30, 60)}
    genes = frame(tus)
    tied = [f"t{i + 1}" for i in range(30, 40)]
    coll = gene_sets("cog", {"tied": tied})
    es = {float(_row(_run(genes, [coll], n_perm=50, seed=s), "tied")["es"]) for s in range(8)}
    assert len(es) > 1, f"ES of a tied set does not depend on the seed: {es}"


# --------------------------------------------------------------------------- leading edge


def test_leading_edge_expands_tus_to_detected_member_genes() -> None:
    """Rule 6: leading-edge TUs in fgsea order, each expanded to its detected member genes in
    tu_position order; non-member and undetected genes of those TUs are left out."""
    tus: dict[str, list[Gene]] = {
        "TU_1": [Gene("z1", 30.0), Gene("a1", 30.2), Gene("u1", None, detected=False)],
    }
    for i in range(2, 11):
        tus[f"TU_{i}"] = [Gene(f"a{i}", 30.0 - 0.1 * i), Gene(f"b{i}", 30.0 - 0.1 * i)]
    rng = np.random.default_rng(4)
    for i, v in enumerate(rng.permutation(140)):
        tus[f"TU_{i + 11}"] = [Gene(f"s{i + 11}", 20.0 + 0.01 * float(v))]
    genes = frame(tus)
    members = (
        {"z1", "a1", "u1"} | {f"a{i}" for i in range(2, 11)} | {f"s{i}" for i in range(11, 16)}
    )
    out = _run(genes, [gene_sets("cog", {"S": members})])
    r = _row(out, "S")
    assert r["size_tus"] == 15 and r["es"] > 0
    # independent expectation: TU means, centred, ranked by the port
    means = genes[genes["detected"]].groupby("tu_id")["abundance_log2"].mean()
    labels = list(means.index)
    prep = prepare(
        means.to_numpy(float) - float(np.median(means)),
        labels,
        {"S": [t for t in labels if t in {f"TU_{i}" for i in range(1, 16)}]},
        10,
        500,
    )
    es, edge = gsea_es(prep.stats, prep.members[0])
    assert r["es"] == pytest.approx(es, abs=1e-12)
    edge_tus = [labels[int(prep.order[k - 1])] for k in edge]
    member_genes = {
        "TU_1": ["z1", "a1"],
        **{f"TU_{i}": [f"a{i}"] for i in range(2, 11)},
        **{f"TU_{i}": [f"s{i}"] for i in range(11, 16)},
    }
    expect = [gn for t in edge_tus for gn in member_genes[t]]
    assert as_list(r["leading_edge"]) == expect
    assert as_list(r["leading_edge"])[:2] == ["z1", "a1"]


def test_es_zero_tie_gives_an_empty_leading_edge_and_a_normal_row() -> None:
    """Scores 3, 2, 1 centre to +1, 0, -1; the single middle TU walks to tops 0.5, bottoms -0.5
    (fgsea's NR = 0 branch): ES = 0 exactly, empty edge, the row is otherwise estimated."""
    genes = singleton_frame([3.0, 2.0, 1.0])
    out = _run(genes, [gene_sets("cog", {"mid": {"g2"}})], size=(1, 500))
    r = _row(out, "mid")
    assert r["es"] == 0.0
    assert as_list(r["leading_edge"]) == []
    assert np.isfinite(r["p_value"]) and np.isfinite(r["q_value"])
    assert "unestimated" not in _flags(r["flags"])


# --------------------------------------------------------------------------- degenerate rows


def _assert_empty_with_flag(out: pd.DataFrame, flag: str) -> None:
    assert len(out) == 0
    assert set(GSEA_COLUMNS) <= set(out.columns)
    flags = list(out.attrs.get("run_flags", []))
    assert flag in flags, flags


def test_abundance_scale_column_absent_gives_no_rows_and_layer_missing() -> None:
    genes, order = _strong()
    coll = gene_sets("cog", {"K": _tag(genes, order[:15])})
    out = _run(genes.drop(columns=["abundance_scale"]), [coll])
    _assert_empty_with_flag(out, "layer_missing:abundance")


def test_no_abundance_input_gives_no_rows_and_layer_missing() -> None:
    """The M4 no-abundance table: abundance_log2 NA and abundance_scale null on every row."""
    genes, order = _strong()
    genes = genes.assign(abundance_log2=math.nan, abundance_scale=None)
    coll = gene_sets("cog", {"K": _tag(genes, order[:15])})
    out = _run(genes, [coll])
    _assert_empty_with_flag(out, "layer_missing:abundance")


@pytest.mark.parametrize("n_ranked", [0, 1])
def test_fewer_than_two_ranked_tus(n_ranked: int) -> None:
    tus = {f"TU_{i + 1}": [Gene(f"g{i + 1}", 10.0)] for i in range(n_ranked)}
    tus["TU_9"] = [Gene("u", None, detected=False)]
    out = _run(frame(tus), [gene_sets("cog", {"K": {"g1", "u"}})], size=(1, 500))
    _assert_empty_with_flag(out, "too_few_ranked_tus")


def test_relative_abundance_scale_still_runs() -> None:
    genes, order = _strong()
    genes = genes.assign(abundance_scale="relative")
    out = _run(genes, [gene_sets("cog", {"K": _tag(genes, order[:15])})])
    assert len(out) == 1
    assert "layer_missing:abundance" not in list(out.attrs.get("run_flags", []))


def test_empty_collection_gives_no_rows_for_it() -> None:
    genes, order = _strong()
    out = _run(genes, [gene_sets("cog", {}), gene_sets("kegg", {"m": _tag(genes, order[:15])})])
    assert set(out["collection"]) == {"kegg"}


def test_gsea_table_does_not_mutate_its_inputs() -> None:
    genes, order = _strong()
    before = genes.copy(deep=True)
    coll = gene_sets("cog", {"K": _tag(genes, order[:15])})
    _run(genes, [coll])
    pd.testing.assert_frame_equal(genes, before)


# --------------------------------------------------------------------------- fake engine

FakeP = Callable[[int], float]


class FakeEngine:
    """Stands in for ``fgsea_simple``: p by set size in TUs; records every call."""

    def __init__(self, p_by_size: Mapping[int, float], *, nes_nan: Sequence[int] = ()) -> None:
        self.p_by_size = dict(p_by_size)
        self.nes_nan = set(nes_nan)
        self.calls: list[tuple[list[int], int]] = []

    def __call__(self, prep: Any, n_perm: int, sample_size: int, rng: Any) -> SimpleResult:
        sizes = [int(s) for s in np.asarray(prep.sizes)]
        self.calls.append((sizes, n_perm))
        k = len(sizes)
        p = np.array([self.p_by_size.get(s, 0.5) for s in sizes], dtype=float)
        bad = np.array([s in self.nes_nan for s in sizes], dtype=bool) | np.isnan(p)
        es = np.full(k, 0.5)
        nes = np.where(bad, np.nan, 1.5)
        p = np.where(bad, np.nan, p)
        log2err = np.where(bad, np.nan, 0.1)
        edges = tuple(
            np.asarray(prep.order)[np.asarray(m)[:1] - 1].astype(np.int64) for m in prep.members
        )
        return SimpleResult(
            es,
            nes,
            p,
            log2err,
            np.zeros(k, dtype=np.int64),
            np.where(bad, 3, 500).astype(np.int64),
            np.zeros(k, dtype=bool),
            edges,
        )


@pytest.fixture
def fake(monkeypatch: pytest.MonkeyPatch) -> Callable[..., FakeEngine]:
    def install(p_by_size: Mapping[int, float], **kw: Any) -> FakeEngine:
        import giap.stats.fgsea as fg

        mod = _g()  # import gsea before patching fg, or its import binds the fake
        eng = FakeEngine(p_by_size, **kw)
        monkeypatch.setattr(fg, "fgsea_simple", eng)
        if hasattr(mod, "fgsea_simple"):
            monkeypatch.setattr(mod, "fgsea_simple", eng)
        return eng

    return install


def _go_frame(n: int = 400) -> tuple[pd.DataFrame, list[str]]:
    genes = _distinct_frame(n, seed=8)
    return genes, list(genes["locus_tag"])


def test_unestimated_rows_are_na_flagged_and_out_of_bh(fake: Callable[..., FakeEngine]) -> None:
    """modeFraction < 10 or NES undefined: p, q, NES, log2err NA; flags unestimated; the other
    q-values are BH without it."""
    genes, tags = _go_frame()
    eng = fake({11: 0.001, 12: 0.04, 13: 0.2, 14: math.nan}, nes_nan=[15])
    coll = gene_sets("cog", {f"s{k}": tags[:k] for k in (11, 12, 13, 14, 15)})
    out = _run(genes, [coll])
    assert eng.calls, "gsea_table did not call fgsea_simple"
    for t in ("s14", "s15"):
        r = _row(out, t)
        assert all(_na(r[c]) for c in ("p_value", "q_value", "nes", "log2err")), r
        assert "unestimated" in _flags(r["flags"])
        assert np.isfinite(r["es"])
    ok = [_row(out, t) for t in ("s11", "s12", "s13")]
    np.testing.assert_allclose([r["q_value"] for r in ok], bh([0.001, 0.04, 0.2]), atol=1e-15)
    assert all("unestimated" not in _flags(r["flags"]) for r in ok)


def _family(
    tags: Sequence[str], child: int, extra: int, nulls: Sequence[int]
) -> tuple[dict[str, list[str]], dict[str, set[str]]]:
    """GO sets: C = tags[:child], P = C + tags[child:child + extra], null sets of the given sizes
    taken from the far end. Distinct sizes let the fake engine key p by size."""
    sets = {"GO:C": list(tags[:child]), "GO:P": list(tags[: child + extra])}
    for i, k in enumerate(nulls):
        sets[f"GO:N{i}"] = list(tags[len(tags) - k - 3 * i : len(tags) - 3 * i])
    parents = {"GO:C": {"GO:P"}, "GO:P": set()} | {f"GO:N{i}": set() for i in range(len(nulls))}
    return sets, parents


def test_example_6_parent_retested_without_child_genes(fake: Callable[..., FakeEngine]) -> None:
    """P = C + 3 TUs; both significant; P re-tested on the 3 others (size 3 < 10: still tested,
    flagged size_filtered); its re-test p replaces its p in the BH vector; q > 0.05 -> redundant."""
    genes, tags = _go_frame()
    sets, parents = _family(tags, 15, 3, [40, 41, 42])
    eng = fake({15: 0.001, 18: 0.002, 3: 0.3, 40: 0.5, 41: 0.6, 42: 0.7})
    out = _run(genes, [gene_sets("go_bp", sets, parents=parents)])
    p, c = _row(out, "GO:P"), _row(out, "GO:C")
    assert p["q_value"] <= ALPHA and c["q_value"] <= ALPHA
    assert p["redundant"] == True  # noqa: E712
    assert c["redundant"] == False  # noqa: E712
    assert "size_filtered" in _flags(p["flags"])
    assert "size_filtered" not in _flags(c["flags"])
    # q_value keeps the main-run BH (the re-test decides only `redundant`)
    np.testing.assert_allclose(p["q_value"], bh([0.001, 0.002, 0.5, 0.6, 0.7])[1], atol=1e-15)
    # re-test: one set of 3 TUs, same n_perm, after the main run
    assert eng.calls[-1] == ([3], 1000)
    assert any(18 in s for s, _ in eng.calls[:-1])


def test_parent_kept_when_its_own_genes_carry_the_signal(fake: Callable[..., FakeEngine]) -> None:
    genes, tags = _go_frame()
    sets, parents = _family(tags, 15, 12, [40, 41, 42])
    fake({15: 0.001, 27: 0.001, 12: 0.002, 40: 0.5, 41: 0.6, 42: 0.7})
    out = _run(genes, [gene_sets("go_bp", sets, parents=parents)])
    p = _row(out, "GO:P")
    assert p["redundant"] == False  # noqa: E712
    assert "size_filtered" not in _flags(p["flags"])


def test_no_retest_without_a_significant_descendant(fake: Callable[..., FakeEngine]) -> None:
    genes, tags = _go_frame()
    sets, parents = _family(tags, 15, 3, [40, 41, 42])
    eng = fake({15: 0.6, 18: 0.001, 40: 0.5, 41: 0.6, 42: 0.7})
    out = _run(genes, [gene_sets("go_bp", sets, parents=parents)])
    assert _row(out, "GO:P")["redundant"] == False  # noqa: E712
    assert all(3 not in s for s, _ in eng.calls)


def test_retest_removes_genes_of_all_significant_descendants(
    fake: Callable[..., FakeEngine],
) -> None:
    """P <- M <- G with G significant and M not: P loses G's genes (transitive descendants)."""
    genes, tags = _go_frame()
    sets = {"GO:G": tags[:11], "GO:M": tags[:24], "GO:P": tags[:30]}
    parents = {"GO:G": {"GO:M"}, "GO:M": {"GO:P"}, "GO:P": set()}
    eng = fake({11: 0.001, 24: 0.4, 30: 0.002, 19: 0.003})
    out = _run(genes, [gene_sets("go_bp", sets, parents=parents)])
    assert eng.calls[-1][0] == [19]  # 30 - 11
    assert _row(out, "GO:P")["redundant"] == False  # noqa: E712
    assert _row(out, "GO:M")["redundant"] == False  # noqa: E712


def test_a_tu_stays_when_another_of_its_genes_remains(fake: Callable[..., FakeEngine]) -> None:
    """Rule 7: removing the child's genes keeps a TU whose other detected gene is in P."""
    tus = {
        f"TU_{i}": [Gene(f"x{i}", 20.0 + 0.01 * i), Gene(f"y{i}", 20.0 + 0.01 * i)]
        for i in range(1, 13)
    }
    tus |= {f"TU_{i}": [Gene(f"s{i}", 10.0 + 0.01 * i)] for i in range(13, 300)}
    genes = frame(tus)
    child = [f"x{i}" for i in range(1, 13)]  # 12 TUs
    parent = child + [f"y{i}" for i in range(1, 5)] + [f"s{i}" for i in range(13, 23)]  # 22 TUs
    sets = {"GO:C": child, "GO:P": parent}
    eng = fake({12: 0.001, 22: 0.001, 14: 0.001})
    _run(genes, [gene_sets("go_bp", sets, parents={"GO:C": {"GO:P"}, "GO:P": set()})])
    assert eng.calls[-1][0] == [14]  # TU_1..TU_4 kept through y1..y4, plus 10 singletons


def test_parent_group_is_the_most_general_significant_ancestor(
    fake: Callable[..., FakeEngine],
) -> None:
    genes, tags = _go_frame()
    sets = {"GO:G": tags[:11], "GO:M": tags[:24], "GO:P": tags[:30], "GO:N": tags[-45:]}
    parents = {"GO:G": {"GO:M"}, "GO:M": {"GO:P"}, "GO:P": set(), "GO:N": set()}
    fake({11: 0.001, 24: 0.001, 30: 0.001, 45: 0.9, 13: 0.001, 6: 0.001})
    out = _run(genes, [gene_sets("go_bp", sets, parents=parents)])
    assert _row(out, "GO:G")["parent_group"] == "GO:P"
    assert _row(out, "GO:M")["parent_group"] == "GO:P"
    assert _row(out, "GO:P")["parent_group"] == "GO:P"
    assert _na(_row(out, "GO:N")["parent_group"])


def test_non_go_collections_are_never_redundant(fake: Callable[..., FakeEngine]) -> None:
    genes, tags = _go_frame()
    sets = {"C": tags[:15], "P": tags[:18]}
    fake({15: 0.001, 18: 0.001, 3: 0.9})
    for c in ("cog", "kegg", "go_slim"):
        out = _run(genes, [gene_sets(c, sets)])
        assert not out["redundant"].astype(bool).any()
        assert out["parent_group"].map(_na).all()
