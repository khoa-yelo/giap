"""M8-13 §1.2: GIPPA's R fgsea call against the direct R fgsea run on the golden data (ADR-0038).

The reference is ``make_golden_r.R``'s one ``fgseaMultilevel`` over the golden stats in file
order (``outputs/fgsea.tsv.gz``, 119 rows: SET_027 drops below minSize). ``gsea_table`` runs on
the fixture of ``helpers.golden_gsea_fixture`` with ``rng = for_step("function.gsea",
GOLDEN_SEED)``.

The golden GSEA scores are median-exact (M8-12b, Slurm job 1879857), so every check is exact.

- **G1 call construction (CI, pure):** one job ``main:golden`` [10, 500]; its sets are the golden
  sets restricted to ranked TUs in golden set order, SET_027 absent; params 101 / 0 / 1000; the
  seed is the second draw after the tie permutation; tie blocks follow the seeded permutation, not
  TU-ID order; ``{tu_id: score}`` equals the golden stats exactly and ``write_call`` round-trips
  them.
- **G2 result mapping (CI, replay):** the golden rows map to the table exactly (es, nes, p,
  log2err, size, leading edges in fgsea's order, q = BH, constants, flags); the written table
  validates against its schema.
- **Example 9 (mutation, A2/A6):** a ``tu_ranking`` that skips the median, on a fixture whose
  abundances are shifted by F1's old median: G1's score equality fails on every TU (CI) and G3a
  fails on every ES (Slurm). Injected by monkeypatch; ``src/`` is untouched.
- The ES oracle reproduces ``calc_gsea_stat.tsv.gz`` (CI), the §3 p rule's NA classes (CI) and
  the A5 vacuity floor (CI).
- **G3a exact wiring (slow + requires_r, A3/example 6):** the live engine behind
  :class:`NormalisingEngine` (seed -> GOLDEN_SEED, tie blocks -> file order): R's
  ``results.tsv`` equals ``fgsea.tsv.gz`` as strings, row for row, and G2's table rules hold.
- **G3b GIPPA's own run (slow + requires_r, A4/examples 7-8):** gated on the deterministic rules
  (live ES = oracle on GIPPA's order, tie-insensitive ES and leading edges = golden, no one-sided
  NA pval/NES); |ΔNES| and the p rule are recorded in the ``gsea_live`` report section.
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
    NormalisingEngine,
    golden_fgsea,
    golden_file_order,
    golden_world,
    log10p_check,
    na_log2err_over_cap,
    oracle_es,
    run_gsea,
)
from helpers.golden_r_gate import (
    GOLDEN,
    decade,
    golden_r_tol,
    golden_seed,
    num,
    read_gz,
    update_report,
)

N_TESTED = 119
DROPPED = "SET_027"
F1_MEDIAN = 0.005705  # the old golden median (F1): the shift of the example-9 mutation
RESULT_COLUMNS = ("pathway", "pval", "log2err", "ES", "NES", "size", "leadingEdge")


def _natural(tus: list[str]) -> list[str]:
    from gippa.function.gsea import _natural_key

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
    from gippa.rng import for_step

    rng = for_step("function.gsea", golden_seed())
    rng.permutation(len(world.stats))
    want = int(rng.integers(1, 2**31 - 1))
    assert type(call.seed) is int and 1 <= call.seed <= 2**31 - 2
    assert call.seed == want


def test_g1_tie_blocks_follow_the_seeded_permutation(world, call) -> None:
    from gippa.rng import for_step

    natural = _natural(list(world.stats))
    perm = for_step("function.gsea", golden_seed()).permutation(len(natural))
    assert call.tu_ids == tuple(natural[int(i)] for i in perm)
    blocks = [v for v, c in Counter(world.stats.values()).items() if c > 1]
    assert len(blocks) == 3  # two values x 50 TUs and 6 zeros (F2)
    for v in blocks:
        order = [t for t in call.tu_ids if world.stats[t] == v]
        assert order != _natural(order), f"tie block {v} is in TU-ID order"


def _score_mismatches(stats: dict[str, float], call) -> dict[str, tuple[float, float]]:
    """G1's exact rule: TU -> (sent score, golden score) wherever they differ (float ``==``)."""
    got = dict(zip(call.tu_ids, np.asarray(call.scores).tolist(), strict=True))
    assert set(got) == set(stats)
    return {t: (got[t], s) for t, s in stats.items() if got[t] != s}


def test_g1_scores_equal_the_golden_scores(world, call) -> None:
    """TU collapse and centring on median-exact data (M8-12b): {tu_id: score} == golden."""
    assert float(np.median(list(world.stats.values()))) == 0.0
    assert _score_mismatches(world.stats, call) == {}


def test_g1_sorted_scores_equal_the_golden_sorted_scores(world, call) -> None:
    got = sorted(np.asarray(call.scores).tolist(), reverse=True)
    assert got == sorted(world.stats.values(), reverse=True)


def _skip_median(monkeypatch: pytest.MonkeyPatch) -> None:
    """Example 9's fault: ``tu_ranking`` returns the TU means without subtracting the median."""
    import gippa.function.gsea as gsea

    real = gsea.tu_ranking

    def uncentred(genes):
        rank = real(genes)
        means = genes.groupby("tu_id")["abundance_log2"].mean()
        scores = np.array([means[t] for t in rank.tu_ids], dtype=np.float64)
        return gsea.TuRanking(rank.tu_ids, scores, rank.tu_genes)

    monkeypatch.setattr(gsea, "tu_ranking", uncentred)


def test_example9_shift_alone_is_centred_away() -> None:
    """Control: with the real ``tu_ranking`` the shifted fixture still sends golden scores (to
    float rounding of the shift), so the mutation below fails because of the fault alone."""
    world = golden_world(shift=F1_MEDIAN)
    eng = CapturingEngine()
    run_gsea(world, eng, golden_seed())
    (call,) = eng.calls
    bad = _score_mismatches(world.stats, call)
    assert max((abs(a - b) for a, b in bad.values()), default=0.0) <= 1e-15


def test_example9_uncentred_scores_fail_g1(monkeypatch: pytest.MonkeyPatch) -> None:
    _skip_median(monkeypatch)
    world = golden_world(shift=F1_MEDIAN)
    eng = CapturingEngine()
    run_gsea(world, eng, golden_seed())
    (call,) = eng.calls
    bad = _score_mismatches(world.stats, call)
    assert len(bad) == len(world.stats), "G1 must catch the fault on every TU"
    assert all(abs(a - b - F1_MEDIAN) <= 1e-12 for a, b in bad.values())


def test_g1_write_call_round_trips_the_scores(call) -> None:
    from gippa.io.fgsea_r import write_call

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


def _numbers_equal_golden(df) -> None:
    gold = golden_fgsea()
    assert df["term"].tolist() == list(gold)
    for r in df.itertuples():
        g = gold[r.term]
        for col, key in (("es", "ES"), ("nes", "NES"), ("p_value", "pval"), ("log2err", "log2err")):
            assert getattr(r, col) == num(g[key]), (r.term, col)
        assert r.size_tus == int(g["size"]), r.term


def _leading_edges_in_fgsea_order(df) -> None:
    gold = golden_fgsea()
    for r in df.itertuples():
        want = [GENE + t for t in gold[r.term]["leadingEdge"].split(";")]
        assert list(r.leading_edge) == want, r.term


def _q_is_bh_of_golden_p(df) -> None:
    from gippa.stats.multitest import bh

    gold = golden_fgsea()
    want = bh(np.array([num(gold[t]["pval"]) for t in df["term"]], dtype=np.float64))
    assert df["q_value"].to_numpy().tolist() == want.tolist()


def _constants_and_flags(df) -> None:
    assert set(df["method"]) == {"fgsea_multilevel"}
    assert set(df["null_model"]) == {"tu_set_permutation"}
    assert set(df["n_perm"]) == {1000}
    assert all("ci_unestimated" in f for f in df["flags"])
    assert not any({"log2err_unestimated", "unestimated"} & set(f) for f in df["flags"])


G2_RULES = (_numbers_equal_golden, _leading_edges_in_fgsea_order, _q_is_bh_of_golden_p,
            _constants_and_flags)  # fmt: skip


@pytest.mark.parametrize("rule", G2_RULES, ids=lambda f: f.__name__.lstrip("_"))
def test_g2_table_rules(replayed, rule) -> None:
    rule(replayed)


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

    from gippa.io.function_tables import write_enrichment_gsea_tsv

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


def _oracles(
    world, call
) -> dict[str, tuple[tuple[float, frozenset[str]], tuple[float, frozenset[str]]]]:
    """Per sent set: (oracle on GIPPA's sent order, oracle on golden file order)."""
    tus, scores = list(world.stats), list(world.stats.values())
    return {
        term: (oracle_es(call.tu_ids, call.scores, members), oracle_es(tus, scores, members))
        for term, members in call.jobs[0].sets.items()
    }


def _tie_insensitive(oracles) -> set[str]:
    """Sets whose oracle ES agrees under GIPPA's and golden order at fgsea_es_abs_tol (§1.2)."""
    tol = golden_r_tol()["fgsea_es_abs_tol"]
    return {t for t, ((es, _), (es_g, _)) in oracles.items() if abs(es - es_g) <= tol}


def _edge_tie_insensitive(oracles) -> set[str]:
    """Sets whose oracle leading edge is the same set under both orders. Not implied by ES
    agreement: SET_016's ES is identical under both orders, but its edge (a negative ES, edge
    after the minimum) swaps TU_0795 for TU_1719, two TUs tied at -0.734567. §4's leading-edge
    rule is therefore applied on these sets (a tie-sensitive edge may legitimately differ)."""
    return {t for t, ((_, edge), (_, edge_g)) in oracles.items() if edge == edge_g}


def test_a5_gsea_vacuity_floor(world, call) -> None:
    """A5: >= 5 golden sets with p_R < 1e-10 and >= 1 of them tie-insensitive under GIPPA's
    order (computable without R: the oracle on the captured call)."""
    gold = golden_fgsea()
    deep = {t for t, g in gold.items() if num(g["pval"]) < 1e-10}
    insensitive = _tie_insensitive(_oracles(world, call))
    update_report(
        "gsea_wiring",
        {
            "n_sets": len(gold),
            "n_p_below_1e-10": len(deep),
            "n_tie_insensitive": len(insensitive),
            "n_tie_insensitive_p_below_1e-10": len(deep & insensitive),
        },
        merge=True,
    )
    assert len(gold) == N_TESTED
    assert len(deep) >= 5
    assert len(deep & insensitive) >= 1


G1_CHECKS = (test_g1_one_main_job_with_the_size_range, test_g1_params,
             test_g1_write_call_round_trips_the_scores)  # fmt: skip
G1_WORLD_CHECKS = (test_g1_sets_are_the_golden_sets_on_ranked_tus, test_g1_seed_is_the_second_draw,
                   test_g1_tie_blocks_follow_the_seeded_permutation,
                   test_g1_scores_equal_the_golden_scores,
                   test_g1_sorted_scores_equal_the_golden_sorted_scores)  # fmt: skip


def _passes(check, *args) -> bool:
    try:
        check(*args)
    except AssertionError:
        return False
    return True


def test_gsea_wiring_report_records_g1_and_g2(world, call, replayed) -> None:
    """§5 lists "G1/G2 pass" (R3): every G1 check on the captured call and every G2 table rule
    on the replayed table, recorded under ``gsea_wiring``."""
    g1 = [_passes(c, call) for c in G1_CHECKS] + [_passes(c, world, call) for c in G1_WORLD_CHECKS]
    g2 = [_passes(rule, replayed) for rule in G2_RULES]
    update_report("gsea_wiring", {"g1_pass": all(g1), "g2_pass": all(g2)}, merge=True)
    assert all(g1) and all(g2), (g1, g2)


# ---------------------------------------------------------------- G3 on real R (Slurm)

PREFIX = real_prefix()
NO_R = (
    "the gippa-fgsea env is not built (sbatch scripts/slurm/build_envs.sbatch fgsea); "
    "run on Slurm with scripts/slurm/test_branch.sbatch -m requires_r"
)


def _live(tmp_path: Path):
    from gippa.io.fgsea_r import r_fgsea_engine, resolve_launcher

    return r_fgsea_engine(resolve_launcher({"GIPPA_FGSEA_PREFIX": str(PREFIX)}), tmp_root=tmp_path)


def _capture_results(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """The raw ``results.tsv`` text of every live call, before the wrapper parses and deletes
    it (G3a compares strings)."""
    import gippa.io.fgsea_r as fr

    texts: list[str] = []
    real = fr.read_results

    def capturing(path: Path, call):
        texts.append(path.read_text(encoding="utf-8"))
        return real(path, call)

    monkeypatch.setattr(fr, "read_results", capturing)
    return texts


def _string_rows(text: str) -> list[dict[str, str]]:
    import csv
    import io

    return list(csv.DictReader(io.StringIO(text), delimiter="\t", quoting=csv.QUOTE_NONE))


def _string_mismatches(text: str) -> list[tuple[int, str, str, str]]:
    """(row, column, live, golden) for every cell of RESULT_COLUMNS that differs."""
    got = _string_rows(text)
    want = read_gz(GOLDEN / "outputs" / "fgsea.tsv.gz")
    assert len(got) == len(want) == N_TESTED, (len(got), len(want))
    assert {r["job_id"] for r in got} == {"main:golden"}
    return [
        (i, c, g[c], w[c])
        for i, (g, w) in enumerate(zip(got, want, strict=True))
        for c in RESULT_COLUMNS
        if g[c] != w[c]
    ]


@pytest.mark.slow
@pytest.mark.requires_r
@pytest.mark.skipif(PREFIX is None, reason=NO_R)
def test_g3a_normalised_live_run_equals_golden_strings(
    world, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Example 6: live R, seed -> 50000, ties -> file order: the 119 rows of R's results equal
    ``fgsea.tsv.gz`` as strings, and G2's table rules hold on the live table."""
    texts = _capture_results(monkeypatch)
    live = _live(tmp_path)
    eng = NormalisingEngine(live, golden_seed(), golden_file_order())
    df = run_gsea(world, eng, golden_seed())
    (text,) = texts
    (sent,) = eng.calls
    assert sent.seed == golden_seed()
    bad = _string_mismatches(text)
    prov = live.provenance
    update_report(
        "gsea_live",
        {
            "g3a_pass": not bad,
            "g3a_mismatches": len(bad),
            "g3a_first_mismatches": bad[:20],
            "r": prov.r if prov else None,
            "fgsea": prov.fgsea if prov else None,
        },
        merge=True,
    )
    assert bad == [], f"{len(bad)} cells differ, first: {bad[:10]}"
    for rule in G2_RULES:
        rule(df)


@pytest.mark.slow
@pytest.mark.requires_r
@pytest.mark.skipif(PREFIX is None, reason=NO_R)
def test_example9_uncentred_scores_fail_g3a_on_every_es(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Example 9 on real R: with the median skipped (and F1's shift in the fixture) G3a's
    string comparison fails on every ES."""
    _skip_median(monkeypatch)
    texts = _capture_results(monkeypatch)
    eng = NormalisingEngine(_live(tmp_path), golden_seed(), golden_file_order())
    run_gsea(golden_world(shift=F1_MEDIAN), eng, golden_seed())
    (text,) = texts
    es_bad = {i for i, c, _, _ in _string_mismatches(text) if c == "ES"}
    assert es_bad == set(range(N_TESTED))


def _deep() -> set[str]:
    """Golden sets with p_R < 1e-10 (A5's deep tail)."""
    return {t for t, g in golden_fgsea().items() if num(g["pval"]) < 1e-10}


def _g3b_floor_failures(insensitive: set[str], edge_insensitive: set[str]) -> list[str]:
    """G3b's own vacuity floor (R1, PR #297 review). Under a fault that moves every ES (example
    9) every set is tie-sensitive, so the rules against golden would check nothing; the weekly
    ``gsea`` task runs only G3, so G3b holds A5's floor itself. Half the sets is robust to the
    seed: the tie blocks cover 106 of 2,000 TUs."""
    half = math.ceil(N_TESTED / 2)
    bad = []
    if len(insensitive) < half:
        bad.append(f"floor: {len(insensitive)} tie-insensitive sets < {half}")
    if len(edge_insensitive) < half:
        bad.append(f"floor: {len(edge_insensitive)} edge-tie-insensitive sets < {half}")
    if not _deep() & insensitive:
        bad.append("floor: no tie-insensitive set with golden p < 1e-10")
    return bad


def test_example9_uncentred_scores_fail_the_g3b_floor(monkeypatch: pytest.MonkeyPatch) -> None:
    """Example 9 against G3b (CI, no R: the classification uses only the oracle on the captured
    call). With the median skipped every ES moves, so no set is tie-insensitive and the floor
    fails; without the floor G3b passed under this fault (PR #297 review)."""
    _skip_median(monkeypatch)
    world = golden_world(shift=F1_MEDIAN)
    eng = CapturingEngine()
    run_gsea(world, eng, golden_seed())
    (call,) = eng.calls
    oracles = _oracles(world, call)
    insensitive = _tie_insensitive(oracles)
    bad = _g3b_floor_failures(insensitive, _edge_tie_insensitive(oracles))
    assert insensitive == set()
    assert f"floor: 0 tie-insensitive sets < {math.ceil(N_TESTED / 2)}" in bad
    assert "floor: no tie-insensitive set with golden p < 1e-10" in bad


def test_g3b_floor_holds_on_the_golden_call(world, call) -> None:
    """Control: on GIPPA's real call (its own seed and tie order) the floor holds."""
    oracles = _oracles(world, call)
    assert _g3b_floor_failures(_tie_insensitive(oracles), _edge_tie_insensitive(oracles)) == []


def _nes_p_summary(terms: list[str], d_nes: dict[str, float], ratios: dict[str, float],
                   nes_tol: float) -> dict[str, object]:  # fmt: skip
    """|ΔNES| and the p rule over one family of sets (R2: reported per tie class)."""
    dn = [d_nes[t] for t in terms if not math.isnan(d_nes[t])]
    rs = [ratios[t] for t in terms if not math.isnan(ratios[t])]
    return {
        "n_sets": len(terms),
        "max_abs_delta_nes": float(max(dn, default=0.0)),
        "p95_abs_delta_nes": float(np.percentile(dn, 95)) if dn else 0.0,
        "n_nes_over_tol": int(sum(d > nes_tol for d in dn)),
        "p_exceedances": int(sum(r > 1.0 for r in rs)),
        "p_exceedance_sets": sorted(t for t in terms if ratios[t] > 1.0),
        "max_ratio": float(max(rs, default=0.0)),
    }


def _g3b(world, call, df) -> tuple[list[str], dict[str, object]]:
    """G3b's rules on one live run: (failures, the ``gsea_live`` report fields)."""
    tol = golden_r_tol()
    gold = golden_fgsea()
    oracles = _oracles(world, call)
    insensitive = _tie_insensitive(oracles)
    edge_insensitive = _edge_tie_insensitive(oracles)
    bad = _g3b_floor_failures(insensitive, edge_insensitive)
    classes: list[str] = []
    ratios: dict[str, float] = {}
    d_nes: dict[str, float] = {}
    d_es: list[float] = []
    for term, ((es, edge), _) in oracles.items():
        row = df.loc[term]
        g = gold[term]
        live_edge = {x.removeprefix(GENE) for x in row["leading_edge"]}
        if not abs(row["es"] - es) <= tol["fgsea_es_abs_tol"]:
            bad.append(f"{term}: live ES != oracle on the sent order")
        if live_edge != edge:
            bad.append(f"{term}: live edge != oracle edge on the sent order")
        if term in insensitive:
            d_es.append(abs(row["es"] - num(g["ES"])))
            if not d_es[-1] <= tol["fgsea_es_abs_tol"]:
                bad.append(f"{term}: tie-insensitive live ES != golden")
        if term in edge_insensitive and live_edge != set(g["leadingEdge"].split(";")):
            bad.append(f"{term}: edge-tie-insensitive live edge != golden")
        if math.isnan(row["p_value"]) != math.isnan(num(g["pval"])):
            bad.append(f"{term}: one-sided NA pval")
        if math.isnan(row["nes"]) != math.isnan(num(g["NES"])):
            bad.append(f"{term}: one-sided NA NES")
        c, r = log10p_check(row["p_value"], row["log2err"], num(g["pval"]), num(g["log2err"]),
                            tol["fgsea_log10p_k"])  # fmt: skip
        classes.append(c)
        ratios[term] = r
        d_nes[term] = abs(row["nes"] - num(g["NES"]))
    if "na_fail" in classes:
        bad.append("one-sided NA p (na_fail)")
    terms = list(oracles)
    nes_tol = tol["fgsea_nes_abs_tol"]
    overall = _nes_p_summary(terms, d_nes, ratios, nes_tol)
    ins = _nes_p_summary([t for t in terms if t in insensitive], d_nes, ratios, nes_tol)
    sens = _nes_p_summary([t for t in terms if t not in insensitive], d_nes, ratios, nes_tol)
    report: dict[str, object] = {
        "seed": call.seed,
        "es_oracle_pass": not any("oracle" in b for b in bad),
        "g3b_pass": not bad,
        "g3b_failures": bad[:20],
        "n_tie_insensitive": len(insensitive),
        "n_tie_sensitive": N_TESTED - len(insensitive),
        "n_edge_tie_insensitive": len(edge_insensitive),
        "n_tie_insensitive_p_below_1e-10": len(_deep() & insensitive),
        "es_tie_sensitive_sets": sorted(set(terms) - insensitive),
        "edge_tie_sensitive_sets": sorted(set(terms) - edge_insensitive),
        "max_abs_delta_es_tie_insensitive": float(max(d_es, default=0.0)),
        "max_abs_delta_nes": overall["max_abs_delta_nes"],
        "p95_abs_delta_nes": overall["p95_abs_delta_nes"],
        "nes_abs_tol": nes_tol,
        "n_nes_over_tol": overall["n_nes_over_tol"],
        "log10p_k": tol["fgsea_log10p_k"],
        "p_exceedances": overall["p_exceedances"],
        "max_ratio": overall["max_ratio"],
        "tie_insensitive": ins,
        "tie_sensitive": sens,
        # §3's P(>= 3 of n) = 0.4 % holds only where both runs estimate the same ES (R2)
        "p_exceedances_alert": ins["p_exceedances"] >= 3,
        "na_log2err": dict(Counter(c for c in classes if c.startswith("na_"))),
        "na_log2err_over_cap": na_log2err_over_cap(classes, N_TESTED),
        "p_decades": dict(Counter(decade(float(p)) for p in df["p_value"])),
    }
    return bad, report


@pytest.mark.slow
@pytest.mark.requires_r
@pytest.mark.skipif(PREFIX is None, reason=NO_R)
def test_g3b_gippa_own_run(world, tmp_path: Path) -> None:
    """Examples 7-8, GIPPA's own seed and tie order. Gated: the floor (>= ceil(119/2) tie- and
    edge-tie-insensitive sets, >= 1 tie-insensitive set with golden p < 1e-10); live ES equals
    the oracle on the sent ranking (all 119 sets) and its leading edge as a set; on
    tie-insensitive sets live ES equals golden at fgsea_es_abs_tol; on edge-tie-insensitive sets
    the leading edge equals golden as a set; no one-sided NA pval or NES. Recorded, not gated
    (question 4): |ΔNES| and the p rule at k = fgsea_log10p_k, overall and per tie class (the
    exceedance alert on tie-insensitive sets only), NA log2err classes, the live p decades."""
    live = _live(tmp_path)
    eng = CapturingEngine(inner=live)
    df = run_gsea(world, eng, golden_seed()).set_index("term")
    (call,) = eng.calls
    assert len(df) == N_TESTED
    bad, report = _g3b(world, call, df)
    prov = live.provenance
    report |= {"r": prov.r if prov else None, "fgsea": prov.fgsea if prov else None}
    update_report("gsea_live", report, merge=True)
    assert bad == [], bad[:10]


@pytest.mark.slow
@pytest.mark.requires_r
@pytest.mark.skipif(PREFIX is None, reason=NO_R)
def test_example9_uncentred_scores_fail_g3b(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Example 9 on real R against G3b's whole rule set: every ES moves, so no set is
    tie-insensitive and only the floor can fail (live ES still equals the oracle on the sent,
    uncentred ranking). Before the floor (PR #297 review, R1) G3b passed under this fault."""
    _skip_median(monkeypatch)
    world = golden_world(shift=F1_MEDIAN)
    eng = CapturingEngine(inner=_live(tmp_path))
    df = run_gsea(world, eng, golden_seed()).set_index("term")
    (call,) = eng.calls
    bad, report = _g3b(world, call, df)
    assert report["n_tie_insensitive"] == 0
    assert bad and all(b.startswith("floor: ") for b in bad), bad[:10]
    assert "floor: no tie-insensitive set with golden p < 1e-10" in bad
