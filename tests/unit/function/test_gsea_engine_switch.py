"""M8-17b T6 and T8: ``gsea_table`` with an injected engine (docs/specs/M8-17-r-fgsea-engine.md
§5, §7, §8, examples 2-8 and 10, A5, A7; ADR-0038; lead decision 3).

Interface pinned here (test-author reading of §7 and §8):

- ``gsea_table(genes, colls, size, cfg, alpha, rng, replicons=None, *, engine=SIMPLE_ENGINE)``;
  ``SIMPLE_ENGINE`` (``method = "fgsea_simple"``) is exported by ``gippa.function.gsea``.
- Main call: one ``GseaJob("main:<collection>", sets, size[0], size[1])`` per collection, in the
  order of ``colls``, holding only the terms with ``size[0] <= size_tus <= size[1]`` and
  ``size_tus < n``, in the collection's term order, each set's member TUs in natural TU order
  (the order of ``tu_ranking``). No job at all: the engine is not called.
- ``EngineCall.tu_ids`` is the ranking in ``rng.permutation(n)`` order, ``scores`` the same order,
  ``params`` is ``cfg`` and ``seed`` the next ``int(rng.integers(1, 2**31 - 1))`` draw (§5).
- Re-tests: one second call with every ``retest:<collection>:<term>`` job (single set, ``[1, n]``)
  of every GO collection, in collection then term order; same ranking and seed. No re-test job:
  no second call. ``k = 0``: no job, ``redundant`` True.
- Row mapping (§8 table), flags ``unestimated`` (p NA) and the new ``log2err_unestimated`` (p
  finite, log2err NA), the leading edge back to genes, and a whole-ranking set's NA row.
- An engine row whose ``size`` differs from Python's count fails the call: ``ValueError`` (the
  pure module cannot raise the io error class; an assumption flagged to the lead).
"""

from __future__ import annotations

import ast
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.fgsea_switch_data import ExplodingEngine, RecordingEngine, flag_set
from helpers.gsea_data import Gene, as_list, frame, gene_sets, singleton_frame, spread
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.config.model import GseaConfig
from gippa.rng import for_step
from gippa.stats.multitest import bh

GSEA_PY = Path(__file__).resolve().parents[3] / "src" / "gippa" / "function" / "gsea.py"
CFG = GseaConfig()
SIZE = (10, 500)
ALPHA = 0.05
MAX_SEED = 2**31 - 2


def _g() -> Any:
    from gippa.function import gsea

    return gsea


def _run(
    genes: pd.DataFrame,
    colls: Sequence[Any],
    engine: Any,
    *,
    seed: int = 7,
    size: tuple[int, int] = SIZE,
    cfg: GseaConfig = CFG,
) -> pd.DataFrame:
    return _g().gsea_table(
        genes, list(colls), size, cfg, ALPHA, for_step("function.gsea", seed), engine=engine
    )


def _values(n: int, seed: int = 11) -> np.ndarray:
    return spread(n, np.random.default_rng(seed))


def _genes(n: int) -> list[str]:
    return [f"g{i}" for i in range(1, n + 1)]


def _row(df: pd.DataFrame, term: str) -> pd.Series:
    sub = df[df["term"] == term]
    assert len(sub) == 1, (term, df[["collection", "term"]])
    return sub.iloc[0]


def _expected_ranking(n: int, seed: int) -> tuple[np.ndarray, int]:
    rng = for_step("function.gsea", seed)
    perm = rng.permutation(n)
    return perm, int(rng.integers(1, 2**31 - 1))


def _natural_tus(genes: Sequence[str]) -> tuple[str, ...]:
    """TUs of singleton genes g<i> (TU_<i>) in natural order."""
    return tuple(f"TU_{i}" for i in sorted(int(g[1:]) for g in genes))


# ---------------------------------------------------------------- the interface (§7)


def test_simple_engine_is_exported_and_named() -> None:
    eng = _g().SIMPLE_ENGINE
    assert eng.method == "fgsea_simple"
    assert callable(eng)


def test_gsea_table_takes_engine_keyword_only() -> None:
    import inspect

    sig = inspect.signature(_g().gsea_table)
    p = sig.parameters["engine"]
    assert p.kind is inspect.Parameter.KEYWORD_ONLY
    assert p.default is _g().SIMPLE_ENGINE


def test_gsea_module_imports_nothing_from_io() -> None:
    tree = ast.parse(GSEA_PY.read_text(encoding="utf-8"))
    mods = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            mods.append(node.module)
        elif isinstance(node, ast.Import):
            mods += [a.name for a in node.names]
    assert not [m for m in mods if m == "gippa.io" or m.startswith("gippa.io.")], mods
    assert not [m for m in mods if m.startswith("gippa.steps")], mods


# ---------------------------------------------------------------- example 2: one call


def _example2() -> tuple[pd.DataFrame, list[Any]]:
    """n = 300 singleton TUs; cog 14 in-range terms + 2 too small, kegg 9 in-range + 1 small."""
    n = 300
    genes = singleton_frame(_values(n))
    r = np.random.default_rng(2)
    tags = _genes(n)

    def pick(k: int) -> list[str]:
        return [str(t) for t in r.choice(tags, k, replace=False)]

    cog = {f"C{i:02d}": pick(10 + i) for i in range(1, 15)}
    cog |= {"S1": pick(3), "S2": pick(8)}
    kegg = {f"K{i:05d}": pick(12 + i) for i in range(1, 10)}
    kegg |= {"KS": pick(5)}
    return genes, [gene_sets("cog", cog), gene_sets("kegg", kegg)]


def test_example_2_exactly_one_call_with_two_jobs() -> None:
    genes, colls = _example2()
    eng = RecordingEngine()
    out = _run(genes, colls, eng, seed=7)
    assert len(eng.calls) == 1
    call = eng.calls[0]
    assert [(j.job_id, j.min_size, j.max_size) for j in call.jobs] == [
        ("main:cog", 10, 500),
        ("main:kegg", 10, 500),
    ]
    assert list(call.jobs[0].sets) == [f"C{i:02d}" for i in range(1, 15)]
    assert list(call.jobs[1].sets) == [f"K{i:05d}" for i in range(1, 10)]
    assert call.params == CFG
    assert len(out) == 23
    assert set(out["method"]) == {"fgsea_multilevel"}
    assert set(out["n_perm"]) == {1000}


def test_example_2_sets_are_natural_order_member_tus() -> None:
    genes, colls = _example2()
    eng = RecordingEngine()
    _run(genes, colls, eng, seed=7)
    job_of = {j.job_id: j for j in eng.calls[0].jobs}
    for coll in colls:
        job = job_of[f"main:{coll.collection}"]
        for term, members in job.sets.items():
            assert isinstance(members, tuple)
            assert members == _natural_tus(coll.sets[term]), term


def test_example_2_stats_in_perm_order_and_seed_is_second_draw() -> None:
    genes, colls = _example2()
    eng = RecordingEngine()
    _run(genes, colls, eng, seed=7)
    call = eng.calls[0]
    rank = _g().tu_ranking(genes)
    perm, seed = _expected_ranking(300, 7)
    assert call.tu_ids == tuple(rank.tu_ids[int(i)] for i in perm)
    np.testing.assert_array_equal(np.asarray(call.scores), np.asarray(rank.scores)[perm])
    assert call.seed == seed
    assert 1 <= call.seed <= MAX_SEED


def test_example_2_gsea_draws_only_perm_and_seed() -> None:
    """§5: the fakes do not consume ``rng``; ``gsea_table`` draws the permutation, then the seed."""
    genes, colls = _example2()
    rng = for_step("function.gsea", 7)
    _g().gsea_table(genes, colls, SIZE, CFG, ALPHA, rng, engine=RecordingEngine())
    ref = for_step("function.gsea", 7)
    ref.permutation(300)
    ref.integers(1, 2**31 - 1)
    assert rng.random() == ref.random()


def test_seed_changes_with_cfg_seed() -> None:
    genes, colls = _example2()
    a, b = RecordingEngine(), RecordingEngine()
    _run(genes, colls, a, seed=7)
    _run(genes, colls, b, seed=8)
    assert a.calls[0].seed == _expected_ranking(300, 7)[1]
    assert b.calls[0].seed == _expected_ranking(300, 8)[1]
    assert a.calls[0].seed != b.calls[0].seed


def test_engine_rows_map_to_columns() -> None:
    """§8 row-mapping table, one scripted row."""
    genes, colls = _example2()
    eng = RecordingEngine(
        rows={("main:kegg", "K00003"): {"es": -0.42, "nes": -1.7, "p": 0.004, "log2err": 0.31}}
    )
    out = _run(genes, colls, eng, seed=7)
    r = _row(out, "K00003")
    assert r["collection"] == "kegg"
    assert r["es"] == -0.42 and r["stat"] == -0.42
    assert r["nes"] == -1.7 and r["effect"] == -1.7
    assert r["p_value"] == 0.004
    assert r["log2err"] == 0.31
    assert r["size_tus"] == len(colls[1].sets["K00003"])
    assert r["method"] == "fgsea_multilevel"
    assert r["null_model"] == "tu_set_permutation"
    assert r["n_perm"] == CFG.n_perm_simple
    assert math.isnan(r["ci_low"]) and math.isnan(r["ci_high"])
    assert flag_set(r["flags"]) >= {"ci_unestimated"}
    kegg = out[out["collection"] == "kegg"].set_index("term")
    p = kegg["p_value"].to_numpy(float)
    np.testing.assert_allclose(kegg["q_value"].to_numpy(float), bh(p), rtol=0, atol=0)


def test_method_column_is_the_engines() -> None:
    genes, colls = _example2()
    out = _run(genes, colls, RecordingEngine(method="fgsea_simple"), seed=7)
    assert set(out["method"]) == {"fgsea_simple"}


def test_n_perm_is_cfg_n_perm_simple() -> None:
    genes, colls = _example2()
    cfg = GseaConfig(n_perm_simple=200)
    eng = RecordingEngine()
    out = _run(genes, colls, eng, seed=7, cfg=cfg)
    assert set(out["n_perm"]) == {200}
    assert eng.calls[0].params == cfg


def test_size_range_is_passed_to_the_main_jobs() -> None:
    genes, colls = _example2()
    eng = RecordingEngine()
    _run(genes, colls, eng, seed=7, size=(15, 20))
    for job in eng.calls[0].jobs:
        assert (job.min_size, job.max_size) == (15, 20)
        assert all(15 <= len(m) <= 20 for m in job.sets.values())


# ---------------------------------------------------------------- examples 7 and 8: zero calls


def test_example_7_no_abundance_zero_calls() -> None:
    genes, colls = _example2()
    genes = genes.assign(abundance_log2=math.nan, abundance_scale=None)
    out = _run(genes, colls, ExplodingEngine())
    assert len(out) == 0
    assert list(out.columns) == list(_g().COLUMNS)
    assert out.attrs["run_flags"] == ["layer_missing:abundance"]


def test_example_8_every_set_out_of_range_zero_calls() -> None:
    genes = singleton_frame(_values(60))
    colls = [
        gene_sets("cog", {f"C{i}": _genes(60)[i : i + 3 + i % 6] for i in range(10)}),
        gene_sets("kegg", {}),
    ]
    out = _run(genes, colls, ExplodingEngine())
    assert len(out) == 0
    assert out.attrs["run_flags"] == []


def test_example_8_no_collection_zero_calls() -> None:
    out = _run(singleton_frame(_values(30)), [], ExplodingEngine())
    assert len(out) == 0 and out.attrs["run_flags"] == []


def test_too_few_ranked_tus_zero_calls() -> None:
    genes = singleton_frame([20.0])
    out = _run(genes, [gene_sets("cog", {"J": ["g1"]})], ExplodingEngine(), size=(1, 500))
    assert len(out) == 0
    assert out.attrs["run_flags"] == ["too_few_ranked_tus"]


def test_collection_without_in_range_term_gets_no_job() -> None:
    genes = singleton_frame(_values(60))
    colls = [
        gene_sets("cog", {"tiny": _genes(60)[:3]}),
        gene_sets("kegg", {"K00001": _genes(60)[:12]}),
    ]
    eng = RecordingEngine()
    _run(genes, colls, eng)
    assert len(eng.calls) == 1
    assert [j.job_id for j in eng.calls[0].jobs] == ["main:kegg"]


def test_non_finite_score_raises_before_any_engine_call() -> None:
    """§8 edge case "NA scores": a non-finite score reaches ``EngineCall``, which raises."""
    vals = list(_values(30))
    vals[4] = math.inf
    genes = singleton_frame(vals)
    with pytest.raises(ValueError):
        _run(genes, [gene_sets("kegg", {"K1": _genes(30)[:12]})], ExplodingEngine())


# ---------------------------------------------------------------- example 6: batched redundancy


def _go_world(q_minus_d: int = 0) -> tuple[pd.DataFrame, list[Any]]:
    """n = 60. go_bp: GO:P (TU 1-22) with significant child GO:C (TU 13-22): P - C = 12 TUs.
    go_mf: GO:Q (TU 30-41, plus ``q_minus_d`` more) with significant child GO:D (TU 30-41)."""
    genes = singleton_frame(_values(60))
    tags = _genes(60)
    bp = gene_sets("go_bp", {"GO:P": tags[0:22], "GO:C": tags[12:22]}, parents={"GO:C": ["GO:P"]})
    mf = gene_sets(
        "go_mf",
        {"GO:Q": tags[29 : 41 + q_minus_d], "GO:D": tags[29:41]},
        parents={"GO:D": ["GO:Q"]},
    )
    return genes, [bp, mf]


def _sig(retest_p: Mapping[tuple[str, str], float] | None = None) -> RecordingEngine:
    rows: dict[tuple[str, str], dict[str, Any]] = {
        ("main:go_bp", "GO:P"): {"p": 0.001},
        ("main:go_bp", "GO:C"): {"p": 0.002},
        ("main:go_mf", "GO:Q"): {"p": 0.001},
        ("main:go_mf", "GO:D"): {"p": 0.002},
    }
    for k, p in (retest_p or {}).items():
        rows[k] = {"p": p}
    return RecordingEngine(rows=rows)


def test_example_6_two_calls_second_is_the_batched_retest() -> None:
    genes, colls = _go_world()
    eng = _sig({("retest:go_bp:GO:P", "GO:P"): 0.2})
    out = _run(genes, colls, eng, seed=7)
    assert len(eng.calls) == 2
    main, re = eng.calls
    assert [j.job_id for j in main.jobs] == ["main:go_bp", "main:go_mf"]
    assert len(re.jobs) == 1
    job = re.jobs[0]
    assert job.job_id == "retest:go_bp:GO:P"
    assert (job.min_size, job.max_size) == (1, 60)
    assert dict(job.sets) == {"GO:P": tuple(f"TU_{i}" for i in range(1, 13))}
    assert re.tu_ids == main.tu_ids
    np.testing.assert_array_equal(np.asarray(re.scores), np.asarray(main.scores))
    assert re.seed == main.seed == _expected_ranking(60, 7)[1]
    assert re.params == main.params
    # q > alpha once p_re = 0.2 is substituted: redundant; Q's set minus D is empty: redundant
    assert bool(_row(out, "GO:P")["redundant"]) is True
    q = _row(out, "GO:Q")
    assert bool(q["redundant"]) is True
    assert "size_filtered" in flag_set(q["flags"])
    assert bool(_row(out, "GO:C")["redundant"]) is False
    assert bool(_row(out, "GO:D")["redundant"]) is False


def test_example_6_retest_p_below_alpha_keeps_the_parent() -> None:
    genes, colls = _go_world()
    out = _run(genes, colls, _sig({("retest:go_bp:GO:P", "GO:P"): 0.001}))
    assert bool(_row(out, "GO:P")["redundant"]) is False


def test_example_6_nan_retest_p_leaves_redundant_na() -> None:
    genes, colls = _go_world()
    out = _run(genes, colls, _sig({("retest:go_bp:GO:P", "GO:P"): math.nan}))
    assert pd.isna(_row(out, "GO:P")["redundant"])


def test_retests_of_all_go_collections_share_one_call() -> None:
    genes, colls = _go_world(q_minus_d=4)
    eng = _sig()
    out = _run(genes, colls, eng)
    assert len(eng.calls) == 2
    re = eng.calls[1]
    assert [(j.job_id, j.min_size, j.max_size) for j in re.jobs] == [
        ("retest:go_bp:GO:P", 1, 60),
        ("retest:go_mf:GO:Q", 1, 60),
    ]
    assert dict(re.jobs[1].sets) == {"GO:Q": tuple(f"TU_{i}" for i in range(42, 46))}
    assert "size_filtered" in flag_set(_row(out, "GO:Q")["flags"])  # k = 4 < 10


def test_no_significant_go_term_no_second_call() -> None:
    genes, colls = _go_world()
    eng = RecordingEngine(
        rows={
            (f"main:{c}", t): {"p": 0.5}
            for c, t in [("go_bp", "GO:P"), ("go_bp", "GO:C"), ("go_mf", "GO:Q"), ("go_mf", "GO:D")]
        }
    )
    _run(genes, colls, eng)
    assert len(eng.calls) == 1


def test_non_go_collections_are_never_retested() -> None:
    genes, _ = _go_world()
    tags = _genes(60)
    kegg = gene_sets("kegg", {"P": tags[0:22], "C": tags[12:22]}, parents={"C": ["P"]})
    eng = RecordingEngine(rows={("main:kegg", "P"): {"p": 0.001}, ("main:kegg", "C"): {"p": 0.002}})
    _run(genes, [kegg], eng)
    assert len(eng.calls) == 1


# ---------------------------------------------------------------- examples 3 and 4: flags


def _kegg_world() -> tuple[pd.DataFrame, list[Any]]:
    """n = 60; four kegg sets of 10 TUs spaced 6 apart (no genomic cluster, all detected)."""
    genes = singleton_frame(_values(60))
    tags = _genes(60)
    sets = {f"K0000{k}": [tags[(k - 1) + 6 * i] for i in range(10)] for k in range(1, 5)}
    return genes, [gene_sets("kegg", sets)]


def test_example_3_unbalanced_set_is_unestimated() -> None:
    genes, colls = _kegg_world()
    eng = RecordingEngine(
        rows={
            ("main:kegg", "K00001"): {
                "p": math.nan,
                "nes": math.nan,
                "log2err": math.nan,
                "es": 0.81,
            },
            ("main:kegg", "K00002"): {"p": 0.01},
            ("main:kegg", "K00003"): {"p": 0.02},
            ("main:kegg", "K00004"): {"p": 0.04},
        }
    )
    out = _run(genes, colls, eng)
    r = _row(out, "K00001")
    for c in ("p_value", "q_value", "nes", "effect", "log2err"):
        assert math.isnan(r[c]), c
    assert r["es"] == 0.81 and r["stat"] == 0.81
    assert as_list(r["flags"]) == ["ci_unestimated", "unestimated"]
    rest = out[out["term"] != "K00001"].set_index("term").loc[["K00002", "K00003", "K00004"]]
    np.testing.assert_array_equal(rest["q_value"].to_numpy(float), bh(np.array([0.01, 0.02, 0.04])))
    for t in ("K00002", "K00003", "K00004"):
        assert "unestimated" not in flag_set(_row(out, t)["flags"])


def test_example_4_log2err_not_estimated() -> None:
    genes, colls = _kegg_world()
    eng = RecordingEngine(rows={("main:kegg", "K00002"): {"p": 3.2e-15, "log2err": math.nan}})
    out = _run(genes, colls, eng)
    r = _row(out, "K00002")
    assert r["p_value"] == 3.2e-15
    assert math.isnan(r["log2err"])
    assert as_list(r["flags"]) == ["ci_unestimated", "log2err_unestimated"]
    assert math.isfinite(r["q_value"])
    p = out.set_index("term").loc[["K00001", "K00002", "K00003", "K00004"], "p_value"]
    np.testing.assert_array_equal(
        out.set_index("term").loc[p.index, "q_value"].to_numpy(float), bh(p.to_numpy(float))
    )
    for t in ("K00001", "K00003", "K00004"):
        assert as_list(_row(out, t)["flags"]) == ["ci_unestimated"], t


def test_finite_log2err_has_no_new_flag() -> None:
    genes, colls = _kegg_world()
    out = _run(genes, colls, RecordingEngine())
    assert all("log2err_unestimated" not in flag_set(f) for f in out["flags"])


def test_size_mismatch_fails_the_call() -> None:
    genes, colls = _kegg_world()
    eng = RecordingEngine(rows={("main:kegg", "K00003"): {"size": 11}})
    with pytest.raises(ValueError):
        _run(genes, colls, eng)


# ---------------------------------------------------------------- example 5: leading edge


def _edge_world() -> tuple[pd.DataFrame, Any]:
    tus: dict[str, list[Gene]] = {}
    vals = _values(13, seed=5)
    for i in range(1, 13):
        if i == 7:
            tus["TU_7"] = [Gene("b0007a", float(vals[i])), Gene("b0007b", float(vals[i]) + 0.003)]
        else:
            tus[f"TU_{i}"] = [Gene(f"b{i:04d}", float(vals[i]))]
    tus["TU_9"] = [Gene("b0099", float(vals[9]))]
    genes = frame(tus)
    coll = gene_sets(
        "kegg", {"T": ["b0007b", "b0002", "b0099"], "T2": ["b0007a", "b0007b", "b0002"]}
    )
    return genes, coll


def test_example_5_leading_edge_tus_back_to_set_genes() -> None:
    genes, coll = _edge_world()
    eng = RecordingEngine(
        rows={
            ("main:kegg", "T"): {"leading_edge": ("TU_7", "TU_2")},
            ("main:kegg", "T2"): {"leading_edge": ("TU_2", "TU_7")},
        }
    )
    out = _run(genes, [coll], eng, size=(2, 500))
    assert dict(eng.calls[0].jobs[0].sets) == {
        "T": ("TU_2", "TU_7", "TU_9"),
        "T2": ("TU_2", "TU_7"),
    }
    assert as_list(_row(out, "T")["leading_edge"]) == ["b0007b", "b0002"]
    # fgsea's TU order first, then tu_position order within a TU
    assert as_list(_row(out, "T2")["leading_edge"]) == ["b0002", "b0007a", "b0007b"]


def test_example_5_zero_es_or_empty_edge_gives_empty_cell() -> None:
    genes, coll = _edge_world()
    eng = RecordingEngine(
        rows={
            ("main:kegg", "T"): {"es": 0.0, "leading_edge": ("TU_7",)},
            ("main:kegg", "T2"): {"leading_edge": ()},
        }
    )
    out = _run(genes, [coll], eng, size=(2, 500))
    assert as_list(_row(out, "T")["leading_edge"]) == []
    assert as_list(_row(out, "T2")["leading_edge"]) == []


# ---------------------------------------------------------------- example 10: whole ranking


def test_example_10_set_covering_every_ranked_tu() -> None:
    genes = singleton_frame(_values(40))
    tags = _genes(40)
    coll = gene_sets("cog", {"J": tags, "K": tags[::4]})
    eng = RecordingEngine()
    out = _run(genes, [coll], eng)
    assert len(eng.calls) == 1
    assert list(eng.calls[0].jobs[0].sets) == ["K"]
    r = _row(out, "J")
    assert r["size_tus"] == 40
    for c in ("es", "nes", "p_value", "q_value"):
        assert math.isnan(r[c]), c
    assert as_list(r["flags"]) == ["ci_unestimated", "unestimated"]
    k = _row(out, "K")
    assert k["q_value"] == bh(np.array([k["p_value"]]))[0]  # J is not in K's BH count


# ---------------------------------------------------------------- SIMPLE_ENGINE


def _simple_call(seed: int) -> Any:
    g = _g()
    ids = tuple(f"TU_{i}" for i in range(1, 61))
    scores = np.linspace(2.0, -2.0, 60)
    job = g.GseaJob("main:kegg", {"A": ids[:12], "B": ids[20:40:2]}, 10, 500)
    return g.EngineCall(ids, scores, (job,), CFG, seed)


def test_simple_engine_seeds_itself_from_the_call() -> None:
    eng = _g().SIMPLE_ENGINE
    a, b = eng(_simple_call(101)), eng(_simple_call(101))
    assert a == b
    c = eng(_simple_call(202))
    assert [r.p for r in a["main:kegg"]] != [r.p for r in c["main:kegg"]]
    assert [r.es for r in a["main:kegg"]] == [r.es for r in c["main:kegg"]]


def test_simple_engine_returns_rows_per_job() -> None:
    res = _g().SIMPLE_ENGINE(_simple_call(5))
    assert set(res) == {"main:kegg"}
    rows = res["main:kegg"]
    assert [r.term for r in rows] == ["A", "B"]
    assert [r.size for r in rows] == [12, 10]
    for r in rows:
        assert 0 < r.p <= 1
        assert set(r.leading_edge) <= set(_simple_call(5).jobs[0].sets[r.term])


def test_default_engine_is_simple() -> None:
    genes, colls = _kegg_world()
    rng_a, rng_b = for_step("function.gsea", 3), for_step("function.gsea", 3)
    a = _g().gsea_table(genes, colls, SIZE, CFG, ALPHA, rng_a)
    b = _g().gsea_table(genes, colls, SIZE, CFG, ALPHA, rng_b, engine=_g().SIMPLE_ENGINE)
    pd.testing.assert_frame_equal(a, b)
    assert set(a["method"]) == {"fgsea_simple"}


# ---------------------------------------------------------------- properties (A5)


@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    n=st.integers(min_value=12, max_value=80),
    seed=st.integers(min_value=0, max_value=2**32 - 1),
    data=st.data(),
)
def test_property_call_invariants(n: int, seed: int, data: st.DataObject) -> None:
    genes = singleton_frame(_values(n, seed=seed % 1000))
    tags = _genes(n)
    sizes = data.draw(st.lists(st.integers(1, n), min_size=1, max_size=6), label="sizes")
    sets = {f"S{i}": tags[i % 3 : i % 3 + k] for i, k in enumerate(sizes)}
    lo = data.draw(st.integers(1, 12), label="lo")
    eng = RecordingEngine()
    _run(genes, [gene_sets("kegg", sets)], eng, seed=seed, size=(lo, 500))
    in_range = [t for t in sorted(sets) if lo <= len(set(sets[t])) < n]
    if not in_range:
        assert eng.calls == []
        return
    assert len(eng.calls) == 1
    call = eng.calls[0]
    perm, r_seed = _expected_ranking(n, seed)
    assert call.seed == r_seed and 1 <= call.seed <= MAX_SEED
    assert sorted(call.tu_ids) == sorted(f"TU_{i}" for i in range(1, n + 1))
    rank = _g().tu_ranking(genes)
    assert call.tu_ids == tuple(rank.tu_ids[int(i)] for i in perm)
    assert np.isfinite(np.asarray(call.scores)).all()
    (job,) = call.jobs
    assert list(job.sets) == in_range
    for t, m in job.sets.items():
        assert m == _natural_tus(sets[t])
        assert lo <= len(m) < n


@settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(seed=st.integers(min_value=0, max_value=2**32 - 1))
def test_property_stats_follow_the_seeded_permutation(seed: int) -> None:
    genes, colls = _kegg_world()
    eng = RecordingEngine()
    _run(genes, colls, eng, seed=seed)
    rank = _g().tu_ranking(genes)
    perm, _ = _expected_ranking(60, seed)
    assert eng.calls[0].tu_ids == tuple(rank.tu_ids[int(i)] for i in perm)
