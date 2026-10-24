"""M8-13 §1.2: GIAP's R fgsea call against the direct R fgsea run on the golden data (ADR-0038).

The reference is ``make_golden_r.R``'s one ``fgseaMultilevel`` over the golden stats in file
order (``outputs/fgsea.tsv.gz``, 119 rows: SET_027 drops below minSize). ``gsea_table`` runs on
the fixture of ``helpers.golden_gsea_fixture`` with ``rng = for_step("function.gsea",
GOLDEN_SEED)``.

Written now (they do not need M8-12b):

- **G1 call construction (CI, pure):** one job ``main:golden`` [10, 500]; its sets are the golden
  sets restricted to ranked TUs in golden set order, SET_027 absent; params 101 / 0 / 1000; the
  seed is the second draw after the tie permutation; tie blocks follow the seeded permutation, not
  TU-ID order; scores are the golden scores after GIAP's centring (``tu_ranking`` subtracts the
  median); ``write_call`` round-trips them.
- **G2 result mapping (CI, replay):** the golden rows map to the table exactly (es, nes, p,
  log2err, size, leading edges in fgsea's order, q = BH, constants, flags); the written table
  validates against its schema.
- The ES oracle reproduces ``calc_gsea_stat.tsv.gz`` (CI) and the §3 p rule's NA classes (CI).
- **G3b, the M8-12b-independent part (Slurm, slow + requires_r):** on GIAP's own call the live
  ES equals the oracle under GIAP's order, and leading edges equal the oracle's as sets.

**Not written, because the golden GSEA scores are not median-exact (F1, M8-12b missing):** G1's
exact ``{tu_id: score}`` equality and sorted-score equality, example 9 (the centring mutation),
G3a (exact strings against ``fgsea.tsv.gz``) and G3b's comparisons with the golden ES / NES / p.
"""

from __future__ import annotations

import math
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from helpers.fgsea_r_data import real_prefix
from helpers.golden_gsea_fixture import (
    GENE,
    CapturingEngine,
    GoldenReplayEngine,
    golden_fgsea,
    golden_world,
    log10p_check,
    na_log2err_over_cap,
    oracle_es,
    run_gsea,
)
from helpers.golden_r_gate import GOLDEN, golden_r_tol, golden_seed, num, read_gz, update_report

N_TESTED = 119
DROPPED = "SET_027"


def _natural(tus: list[str]) -> list[str]:
    from giap.function.gsea import _natural_key

    return sorted(tus, key=_natural_key)


@pytest.fixture(scope="module")
def world():
    return golden_world()


@pytest.fixture(scope="module")
def call(world):
    eng = CapturingEngine()
    run_gsea(world, eng, golden_seed())
    assert len(eng.calls) == 1, "a non-GO collection makes exactly one engine call"
    return eng.calls[0]


# ---------------------------------------------------------------- G1 call construction


def test_g1_one_main_job_with_the_size_range(call) -> None:
    assert [(j.job_id, j.min_size, j.max_size) for j in call.jobs] == [("main:golden", 10, 500)]


def test_g1_sets_are_the_golden_sets_on_ranked_tus(world, call) -> None:
    sets = call.jobs[0].sets
    gold = golden_fgsea()
    ranked = set(world.stats)
    assert DROPPED not in sets and DROPPED not in gold
    assert list(sets) == list(gold) == [s for s in world.sets if s != DROPPED]
    assert len(sets) == N_TESTED
    for term, members in sets.items():
        want = [m for m in world.sets[term] if m in ranked]
        assert set(members) == set(want) and len(members) == len(set(want)), term
        assert list(members) == _natural(list(members)), term
        assert len(members) == int(gold[term]["size"]), term


def test_g1_params(call) -> None:
    p = call.params
    assert (p.sample_size, p.eps, p.n_perm_simple) == (101, 0.0, 1000)


def test_g1_seed_is_the_second_draw(world, call) -> None:
    from giap.rng import for_step

    rng = for_step("function.gsea", golden_seed())
    rng.permutation(len(world.stats))
    want = int(rng.integers(1, 2**31 - 1))
    assert type(call.seed) is int and 1 <= call.seed <= 2**31 - 2
    assert call.seed == want


def test_g1_tie_blocks_follow_the_seeded_permutation(world, call) -> None:
    from giap.rng import for_step

    natural = _natural(list(world.stats))
    perm = for_step("function.gsea", golden_seed()).permutation(len(natural))
    assert call.tu_ids == tuple(natural[int(i)] for i in perm)
    blocks = [v for v, c in Counter(world.stats.values()).items() if c > 1]
    assert len(blocks) == 3  # two values x 50 TUs and 6 zeros (F2)
    for v in blocks:
        order = [t for t in call.tu_ids if world.stats[t] == v]
        assert order != _natural(order), f"tie block {v} is in TU-ID order"


def test_g1_scores_are_the_centred_golden_scores(world, call) -> None:
    """TU collapse and centring: each TU's score is golden - median(golden), exactly. Once M8-12b
    makes the median exactly 0 this is G1's exact equality."""
    gold = np.array(list(world.stats.values()), dtype=np.float64)
    med = float(np.median(gold))
    got = dict(zip(call.tu_ids, np.asarray(call.scores).tolist(), strict=True))
    assert set(got) == set(world.stats)
    bad = {t: (got[t], s) for t, s in world.stats.items() if got[t] != s - med}
    assert bad == {}


def test_g1_write_call_round_trips_the_scores(call) -> None:
    from giap.io.fgsea_r import write_call

    with tempfile.TemporaryDirectory() as d:
        write_call(call, Path(d))
        lines = (Path(d) / "stats.tsv").read_text(encoding="utf-8").splitlines()[1:]
    rows = [line.split("\t") for line in lines]
    assert tuple(r[0] for r in rows) == call.tu_ids
    assert [float(r[1]) for r in rows] == np.asarray(call.scores).tolist()


# ---------------------------------------------------------------- G2 result mapping (replay)


@pytest.fixture(scope="module")
def replayed(world):
    return run_gsea(world, GoldenReplayEngine(), golden_seed())


def test_g2_numbers_equal_golden(replayed) -> None:
    gold = golden_fgsea()
    assert replayed["term"].tolist() == list(gold)
    for r in replayed.itertuples():
        g = gold[r.term]
        for col, key in (("es", "ES"), ("nes", "NES"), ("p_value", "pval"), ("log2err", "log2err")):
            assert getattr(r, col) == num(g[key]), (r.term, col)
        assert r.size_tus == int(g["size"]), r.term


def test_g2_leading_edges_in_fgsea_order(replayed) -> None:
    gold = golden_fgsea()
    for r in replayed.itertuples():
        want = [GENE + t for t in gold[r.term]["leadingEdge"].split(";")]
        assert list(r.leading_edge) == want, r.term


def test_g2_q_is_bh_of_golden_p(replayed) -> None:
    from giap.stats.multitest import bh

    gold = golden_fgsea()
    want = bh(np.array([num(gold[t]["pval"]) for t in replayed["term"]], dtype=np.float64))
    assert replayed["q_value"].to_numpy().tolist() == want.tolist()


def test_g2_constants_and_flags(replayed) -> None:
    assert set(replayed["method"]) == {"fgsea_multilevel"}
    assert set(replayed["null_model"]) == {"tu_set_permutation"}
    assert set(replayed["n_perm"]) == {1000}
    assert all("ci_unestimated" in f for f in replayed["flags"])
    assert not any({"log2err_unestimated", "unestimated"} & set(f) for f in replayed["flags"])


def test_g2_na_log2err_and_pval_rows_are_flagged(world) -> None:
    """log2err reporting: an NA log2err keeps p and gets ``log2err_unestimated``; NA p gets
    ``unestimated`` (golden rows altered by the replay engine)."""
    eng = GoldenReplayEngine(override={"SET_001": {"log2err": math.nan},
                                       "SET_002": {"p": math.nan, "nes": math.nan}})  # fmt: skip
    df = run_gsea(world, eng, golden_seed()).set_index("term")
    assert math.isnan(df.loc["SET_001", "log2err"])
    assert df.loc["SET_001", "p_value"] == num(golden_fgsea()["SET_001"]["pval"])
    assert "log2err_unestimated" in df.loc["SET_001", "flags"]
    assert "unestimated" in df.loc["SET_002", "flags"]
    assert math.isnan(df.loc["SET_002", "q_value"])


def test_g2_written_table_validates(tmp_path: Path) -> None:
    """The schema's ``collection`` is an enum, so this table uses ``kegg`` (not ``golden``)."""
    from helpers.fgsea_switch_data import validate_gsea

    from giap.io.function_tables import write_enrichment_gsea_tsv

    df = run_gsea(golden_world("kegg"), GoldenReplayEngine(), golden_seed())
    assert len(df) == N_TESTED
    out = tmp_path / "enrichment_gsea.tsv"
    write_enrichment_gsea_tsv(out, df)
    validate_gsea(out)


# ---------------------------------------------------------------- the ES oracle and the p rule


def test_es_oracle_reproduces_calc_gsea_stat(world) -> None:
    """Golden (file) order; every set with >= 1 ranked member (120), at fgsea_es_abs_tol."""
    tol = golden_r_tol()["fgsea_es_abs_tol"]
    tus, scores = list(world.stats), list(world.stats.values())
    calc = read_gz(GOLDEN / "outputs" / "calc_gsea_stat.tsv.gz")
    assert len(calc) == N_TESTED + 1
    for r in calc:
        members = [m for m in world.sets[r["pathway"]] if m in world.stats]
        es, _ = oracle_es(tus, scores, members)
        assert abs(es - float(r["ES"])) <= tol, r["pathway"]


def test_es_oracle_leading_edges_match_golden(world) -> None:
    tus, scores = list(world.stats), list(world.stats.values())
    for term, g in golden_fgsea().items():
        _, edge = oracle_es(tus, scores, world.sets[term])
        assert edge == frozenset(g["leadingEdge"].split(";")), term


@pytest.mark.parametrize(
    ("args", "want"),
    [
        ((1e-3, 0.1, 1e-3, 0.1), "ok"),
        ((1e-3, 0.1, 1e-4, 0.1), "exceed"),  # Δlog10 p = 1 > 3·0.1414/3.32 = 0.128
        ((1e-3, math.nan, 1e-3, 0.1), "na_log2err_one"),
        ((1e-3, math.nan, 1e-3, math.nan), "na_log2err_both"),
        ((math.nan, 0.1, 1e-3, 0.1), "na_fail"),
        ((math.nan, math.nan, math.nan, math.nan), "na_pval_both"),
    ],
)
def test_log10p_rule_classes(args, want) -> None:
    k = golden_r_tol()["fgsea_log10p_k"]
    assert log10p_check(*args, k)[0] == want


def test_log10p_rule_bounds() -> None:
    k = 3.0
    _, r = log10p_check(1e-3, 0.3, 1e-3 * 10 ** (k * 0.5 / math.log2(10)), 0.4, k)
    assert r == pytest.approx(1.0, rel=1e-12)  # sqrt(0.3² + 0.4²) = 0.5
    _, r1 = log10p_check(1e-3, math.nan, 1e-3 * 10 ** (k * 0.4 / math.log2(10)), 0.4, k)
    assert r1 == pytest.approx(1.0, rel=1e-12)  # the finite side alone
    assert not na_log2err_over_cap(["na_log2err_one"] * 5 + ["ok"] * 114, N_TESTED)
    assert na_log2err_over_cap(["na_log2err_one"] * 3 + ["na_log2err_both"] * 3, N_TESTED)


# ---------------------------------------------------------------- G3b on real R (Slurm)

PREFIX = real_prefix()
NO_R = (
    "the giap-fgsea env is not built (sbatch scripts/slurm/build_envs.sbatch fgsea); "
    "run on Slurm with scripts/slurm/test_branch.sbatch -m requires_r"
)


@pytest.mark.slow
@pytest.mark.requires_r
@pytest.mark.skipif(PREFIX is None, reason=NO_R)
def test_g3b_live_es_and_leading_edges_follow_giap_order(world, tmp_path: Path) -> None:
    """GIAP's own seed and tie order: live ES equals the oracle on the sent ranking at
    fgsea_es_abs_tol for all 119 sets, leading edges equal it as sets, and no row has NA pval
    or NES (the golden side has none: a one-sided NA fails)."""
    from giap.io.fgsea_r import r_fgsea_engine, resolve_launcher

    tol = golden_r_tol()
    live = r_fgsea_engine(resolve_launcher({"GIAP_FGSEA_PREFIX": str(PREFIX)}), tmp_root=tmp_path)
    eng = CapturingEngine(inner=live)
    df = run_gsea(world, eng, golden_seed()).set_index("term")
    (call,) = eng.calls
    assert len(df) == N_TESTED
    gold = golden_fgsea()
    classes, ratios, d_nes = [], [], []
    for term, members in call.jobs[0].sets.items():
        es, edge = oracle_es(call.tu_ids, call.scores, members)
        row = df.loc[term]
        assert abs(row["es"] - es) <= tol["fgsea_es_abs_tol"], term
        assert {g.removeprefix(GENE) for g in row["leading_edge"]} == edge, term
        assert not math.isnan(row["p_value"]) and not math.isnan(row["nes"]), term
        g = gold[term]
        c, r = log10p_check(row["p_value"], row["log2err"], num(g["pval"]), num(g["log2err"]),
                            tol["fgsea_log10p_k"])  # fmt: skip
        classes.append(c)
        ratios.append(r)
        d_nes.append(abs(row["nes"] - num(g["NES"])))
    prov = live.provenance
    update_report(
        "gsea_live_pre_m8_12b",
        {
            "note": "golden scores not median-exact (F1): NES/p deltas are not the gate's",
            "seed": call.seed,
            "r": prov.r if prov else None,
            "fgsea": prov.fgsea if prov else None,
            "es_oracle_pass": True,
            "max_abs_delta_nes": float(np.max(d_nes)),
            "p95_abs_delta_nes": float(np.percentile(d_nes, 95)),
            "p_exceedances": int(sum(r > 1.0 for r in ratios if not math.isnan(r))),
            "max_ratio": float(np.nanmax(ratios)),
            "na_log2err": dict(Counter(c for c in classes if c.startswith("na_"))),
            "na_log2err_over_cap": na_log2err_over_cap(classes, N_TESTED),
        },
    )
    assert "na_fail" not in classes
