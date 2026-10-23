"""giapsim operon scenarios (issue M6-03; GIAP_SPEC.md §11.3, §8.4 ``discordance`` and
``missed_operons``, §8.5 attack 1, §7.4.5, §7.4.7).

Examples 1-7, the degenerate-case rows and lead decision 1 of
docs/specs/M6-03-sim-operon-scenarios.md, the truth-file round trip through
``giapsim.contract_tables``, determinism and Hypothesis invariants. The generator-level
distributional claims (no excess of small |Delta| under the null, an excess under the power
scenarios) are in tests/calibration/test_sim_operon_calibration.py.

Interface assumed beyond the spec's five signatures (test-author choices, reported to the lead):

* the five scenarios live in ``giapsim.scenarios`` and return ``Pair = (LayoutConfig,
  ExpressionConfig)`` like the M5-04 scenarios (pure functions of their arguments, so equal
  arguments give equal pairs); whatever they plant is carried by the configs;
* ``giapsim.scenarios.simulate_operons(cfgs: Pair, seed: int) -> OperonSim`` draws one replicate
  (all randomness from ``seed``) and returns an object with
  - ``layout`` (``giapsim.layout.Layout``) and ``expression`` (``giapsim.expression.Expression``),
  - ``pairs``: one row per pair, at least ``PAIR_KEYS`` below: the §9.3 ``pairs.tsv`` columns the
    sim can know (``pair_id`` unique; ``class`` OP/NO/CV/DV/SW as giap.operons.pairs would call it:
    same-strand adjacent pairs are OP iff ``uniop_p >= 0.5``; OP, NO and SW oriented in the
    direction of transcription with ``a`` upstream, CV and DV with ``a`` the left gene;
    ``abundance_log2_*`` equal to ``abundance_log2_true`` for detected genes and ``delta_log2 =
    abundance_log2_b - abundance_log2_a`` when both are detected) plus three truth columns:
    ``true_operon`` (both genes in one planted TU), ``promoter_between_true`` and
    ``terminator_between_true`` (a true element planted in the gap between the two genes),
  - ``truth``: a mapping of truth frames, keys ``TRUTH_KEYS``: ``pairs`` (``pair_id, class,
    gene_a, gene_b, true_operon, promoter_between_true, terminator_between_true``, every row of
    ``pairs``), ``missed_operons`` (``pair_id, gene_a, gene_b``: the planted co-abundant NO
    pairs), ``discordant_pairs`` (``pair_id, gene_a, gene_b, direction, element, effect_log2``;
    ``direction`` is ``downstream_higher`` or ``upstream_higher`` (§7.4.5 wording), ``element``
    ``promoter`` or ``terminator``, ``effect_log2`` the signed planted Delta) and ``operon_t``
    (one row ``loc, scale, df``: the generating Student-t of the within-operon Delta); scenarios
    without a planted effect have empty ``missed_operons`` / ``discordant_pairs`` frames;
* ``giapsim.contract_tables.OPERON_TRUTH_TYPES`` declares the truth files as
  ``CONTRACT_TYPES`` declares C1-C8 and C15 (key -> ((column, logical dtype), ...)), with at least
  the columns above in that order; ``write_operon_truth(sim, out_dir) -> dict[str, Path]`` writes
  ``out_dir/truth/<key>.tsv`` (giap TSV conventions: ``true``/``false``, ``NA``) and
  ``read_operon_truth(out_dir) -> dict[str, DataFrame]`` reads them back;
* example 2's "Delta of OP pairs" and example 3's background are the *true* classes (class OP and
  ``true_operon``; class NO and not ``true_operon``), both genes detected, so a UniOP error the
  scenario plants on purpose (item 6) does not contaminate the generating distribution;
* example 3 "the 50 planted pairs have the smallest |Delta|" is read as: the planted pairs' median
  |Delta| is at most a quarter of the unplanted NO pairs' median and below the OP, CV and SW
  medians (both genes detected); with ~600 independent NO pairs a literal "the 50 smallest" cannot
  hold at |Delta| ~ 0.1 SD, since ~8 % of independent pairs fall below 0.1 SD by chance;
* M6-09 example 1 ("≥ 60 % of planted pairs recovered") needs the planted pairs measurable, so at
  least 60 % of planted pairs have both genes detected (examples 3 and 4);
* the "no convergent pair" degenerate row is reached with ``p_same_strand = 1`` (the only strand
  knob of ``LayoutConfig``): ``simulate_operons`` raises ``ValueError`` naming convergent pairs.
"""

from __future__ import annotations

import dataclasses
import sys
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import expression_oracle as eo
import oric_oracle as oo

N = 1800
SCENARIOS = ("operon_null", "missed_operon_power", "discordance_power", "operon_attack_ar1",
             "operon_attack_detection")  # fmt: skip
PAIR_KEYS = (
    "pair_id", "class", "gene_a", "gene_b", "strand", "intergenic_bp", "uniop_p", "detected_a",
    "detected_b", "abundance_log2_a", "abundance_log2_b", "delta_log2", "true_operon",
    "promoter_between_true", "terminator_between_true",
)  # fmt: skip
TRUTH_KEYS = ("pairs", "missed_operons", "discordant_pairs", "operon_t")
TRUTH_COLUMNS = {
    "pairs": ("pair_id", "class", "gene_a", "gene_b", "true_operon", "promoter_between_true",
              "terminator_between_true"),
    "missed_operons": ("pair_id", "gene_a", "gene_b"),
    "discordant_pairs": ("pair_id", "gene_a", "gene_b", "direction", "element", "effect_log2"),
    "operon_t": ("loc", "scale", "df"),
}  # fmt: skip


def _sc() -> Any:
    import giapsim.scenarios

    return giapsim.scenarios


def _call(name: str, *args: Any) -> tuple[Any, Any]:
    out: tuple[Any, Any] = getattr(_sc(), name)(*args)
    return out


@cache
def sim(name: str, args: tuple[Any, ...], seed: int) -> Any:
    return _sc().simulate_operons(_call(name, *args), seed)


def _b(s: pd.Series) -> np.ndarray:
    return np.asarray(s.astype(bool).to_numpy(), dtype=bool)


def both(p: pd.DataFrame) -> pd.DataFrame:
    return p[_b(p["detected_a"]) & _b(p["detected_b"])]


def true_ab(s: Any) -> dict[str, float]:
    g = s.expression.genes
    return {str(t): float(v) for t, v in zip(g["locus_tag"], g["abundance_log2_true"],
                                             strict=True)}  # fmt: skip


def true_delta(s: Any, p: pd.DataFrame) -> np.ndarray:
    ab = true_ab(s)
    return np.array([ab[str(b)] - ab[str(a)] for a, b in zip(p["gene_a"], p["gene_b"],
                                                             strict=True)])  # fmt: skip


def op_true(p: pd.DataFrame) -> pd.DataFrame:
    return p[(p["class"] == "OP").to_numpy() & _b(p["true_operon"])]


def no_true(p: pd.DataFrame) -> pd.DataFrame:
    return p[(p["class"] == "NO").to_numpy() & ~_b(p["true_operon"])]


def t_params(s: Any) -> tuple[float, float, float]:
    t = s.truth["operon_t"]
    assert len(t) == 1
    return float(t["loc"].iloc[0]), float(t["scale"].iloc[0]), float(t["df"].iloc[0])


# ----------------------------------------------------------------------------- interface


def test_every_scenario_and_the_driver_exist() -> None:
    m = _sc()
    for name in (*SCENARIOS, "simulate_operons"):
        assert callable(getattr(m, name, None)), name


@pytest.mark.parametrize(
    "name,args",
    [
        ("operon_null", (600,)),
        ("missed_operon_power", (600, 10, 0.1)),
        ("discordance_power", (600, 10, 2.0)),
        ("operon_attack_ar1", (600, 0.5)),
        ("operon_attack_detection", (600, "prophage")),
    ],
)
def test_scenarios_return_pure_config_pairs(name: str, args: tuple[Any, ...]) -> None:
    from giapsim.expression import ExpressionConfig
    from giapsim.layout import LayoutConfig

    cfgs = _call(name, *args)
    assert isinstance(cfgs, tuple) and len(cfgs) == 2
    assert isinstance(cfgs[0], LayoutConfig) and isinstance(cfgs[1], ExpressionConfig)
    assert _call(name, *args) == cfgs  # a pure function of its arguments
    assert cfgs[0].n_genes_chromosome == args[0]
    assert tuple(cfgs[0].plasmids) == () and tuple(cfgs[0].linear_contigs) == ()


def test_operon_null_has_a_rho_keyword() -> None:
    a, b = _call("operon_null", 600), _sc().operon_null(600, rho=0.5)
    assert a[1].ar1_rho == 0.0 and b[1].ar1_rho == 0.5


def test_pairs_and_truth_frames_have_the_declared_columns() -> None:
    s = sim("operon_null", (600,), 1)
    assert set(PAIR_KEYS) <= set(s.pairs.columns), set(PAIR_KEYS) - set(s.pairs.columns)
    assert set(s.truth) >= set(TRUTH_KEYS)
    for key, cols in TRUTH_COLUMNS.items():
        assert tuple(s.truth[key].columns[: len(cols)]) == cols, key
    assert s.pairs["pair_id"].is_unique
    assert set(s.pairs["class"]) <= {"OP", "NO", "CV", "DV", "SW"}


def test_pairs_agree_with_the_layout() -> None:
    """Classes, orientation and ``true_operon`` follow the layout (TUs, strands, UniOP p)."""
    s = sim("operon_null", (600,), 2)
    g = s.layout.genes.set_index("locus_tag")
    tu = {str(m): str(t) for t, gs in zip(s.layout.tus["tu_id"], s.layout.tus["genes"],
                                          strict=True) for m in gs}  # fmt: skip
    op = s.layout.operon_pairs
    raw = {frozenset((str(a), str(b))): float(p)
           for a, b, p in zip(op["gene_a"], op["gene_b"], op["uniop_p"], strict=True)}  # fmt: skip
    half = float(s.layout.replicons["length"].iloc[0]) / 2.0
    for r in s.pairs.to_dict("records"):
        a, b, cls = str(r["gene_a"]), str(r["gene_b"]), str(r["class"])
        sa, sb = g.at[a, "strand"], g.at[b, "strand"]
        if cls in ("OP", "NO"):
            assert sa == sb
            p = raw[frozenset((a, b))]
            assert float(r["uniop_p"]) == pytest.approx(p)
            assert (cls == "OP") == (p >= 0.5)
            assert bool(r["true_operon"]) == (tu[a] == tu[b])
        else:
            assert not bool(r["true_operon"])
        if cls in ("OP", "NO", "SW"):  # transcription order: a upstream (wrap pairs reversed)
            assert sa == sb
            sa0, sb0 = int(g.at[a, "start"]), int(g.at[b, "start"])
            wrap = abs(sa0 - sb0) > half
            assert ((sa0 < sb0) == (sa == "+")) != wrap, r
        if cls == "CV":
            assert (sa, sb) == ("+", "-")
        if cls == "DV":
            assert (sa, sb) == ("-", "+")


def test_delta_is_b_minus_a_of_the_true_abundance_when_both_detected() -> None:
    s = sim("operon_null", (600,), 3)
    p = both(s.pairs)
    assert len(p) > 0
    np.testing.assert_allclose(p["delta_log2"].astype(float).to_numpy(), true_delta(s, p),
                               atol=1e-9)  # fmt: skip
    det = dict(zip(s.expression.genes["locus_tag"], s.expression.genes["detected_true"],
                   strict=True))  # fmt: skip
    for side in ("a", "b"):
        assert (s.pairs[f"detected_{side}"].astype(bool)
                == s.pairs[f"gene_{side}"].map(det).astype(bool)).all()  # fmt: skip


def test_truth_pairs_cover_every_pair_row() -> None:
    s = sim("discordance_power", (600, 10, 2.0), 4)
    tp = s.truth["pairs"].set_index("pair_id")
    pp = s.pairs.set_index("pair_id")
    assert set(tp.index) == set(pp.index)
    for col in TRUTH_COLUMNS["pairs"][1:]:
        assert (tp[col].astype(str) == pp.loc[tp.index, col].astype(str)).all(), col
    ids = set(pp.index)
    assert set(s.truth["discordant_pairs"]["pair_id"]) <= ids
    assert set(s.truth["missed_operons"]["pair_id"]) <= ids


# ----------------------------------------------------------------------------- examples


def test_example1_operon_null_pair_counts() -> None:
    """``operon_null(1800)``: >= 300 OP, >= 200 NO, >= 100 CV, >= 30 SW pairs with both genes
    detected (every seed of 3)."""
    for seed in range(3):
        counts = both(sim("operon_null", (N,), seed).pairs)["class"].value_counts()
        assert counts.get("OP", 0) >= 300, (seed, dict(counts))
        assert counts.get("NO", 0) >= 200, (seed, dict(counts))
        assert counts.get("CV", 0) >= 100, (seed, dict(counts))
        assert counts.get("SW", 0) >= 30, (seed, dict(counts))


def test_example2_op_delta_follows_the_generating_t() -> None:
    """``fit_student_t`` on the both-detected true OP Delta recovers ``truth["operon_t"]``
    (loc within 0.15 scale, scale within 20 %, df within a factor 2.5 or both >= 30); KS against
    the generating t at p > 0.001."""
    from giap.stats.tdist import fit_student_t

    s = sim("operon_null", (N,), 0)
    loc, scale, df = t_params(s)
    assert scale > 0 and df > 0
    d = both(op_true(s.pairs))["delta_log2"].astype(float).to_numpy()
    assert len(d) >= 300
    fit = fit_student_t(d, df_bounds=(1.0, 200.0))
    assert fit.converged
    assert abs(fit.loc - loc) <= 0.15 * scale, (fit, loc)
    assert abs(fit.scale / scale - 1.0) <= 0.2, (fit, scale)
    assert (df / 2.5 <= fit.df <= df * 2.5) or (df >= 30 and fit.df >= 30), (fit, df)
    ks = stats.kstest(d, stats.t(df, loc=loc, scale=scale).cdf)
    assert ks.pvalue > 0.001, ks


def test_example3_planted_pairs_have_the_smallest_delta() -> None:
    """``missed_operon_power(1800, 50, 0.1)``: in >= 80 % of 5 seeds the planted pairs' median
    |Delta| is <= 1/4 of the unplanted NO median and below the OP, CV and SW medians."""
    ok = 0
    for seed in range(5):
        s = sim("missed_operon_power", (N, 50, 0.1), seed)
        planted = set(s.truth["missed_operons"]["pair_id"])
        assert len(planted) == 50
        p = both(s.pairs)
        is_p = p["pair_id"].isin(planted).to_numpy()
        assert is_p.sum() >= 30, is_p.sum()  # >= 60 % measurable (M6-09 example 1)
        m = {
            "planted": np.median(np.abs(true_delta(s, p[is_p]))),
            "NO": np.median(np.abs(true_delta(s, no_true(p[~is_p])))),
            "OP": np.median(np.abs(true_delta(s, op_true(p)))),
            "CV": np.median(np.abs(true_delta(s, p[p["class"] == "CV"]))),
            "SW": np.median(np.abs(true_delta(s, p[p["class"] == "SW"]))),
        }
        ok += m["planted"] <= 0.25 * m["NO"] and all(m["planted"] < m[k] for k in
                                                       ("OP", "CV", "SW"))  # fmt: skip
    assert ok >= 4, ok


def test_example3_planted_pairs_are_uniop_false_negatives() -> None:
    """Item 2/6: planted pairs are adjacent same-strand NO pairs with UniOP p < 0.5 and are not
    in one TU (a planted co-abundant NO pair is exactly a UniOP false negative)."""
    s = sim("missed_operon_power", (N, 50, 0.1), 0)
    pp = s.pairs.set_index("pair_id").loc[list(s.truth["missed_operons"]["pair_id"])]
    assert (pp["class"] == "NO").all()
    assert (pp["uniop_p"].astype(float) < 0.5).all()
    assert not _b(pp["true_operon"]).any()
    tr = s.truth["missed_operons"].set_index("pair_id")
    assert (tr["gene_a"] == pp["gene_a"]).all() and (tr["gene_b"] == pp["gene_b"]).all()
    assert s.truth["discordant_pairs"].empty


def test_example4_discordance_power_truth() -> None:
    """``discordance_power(1800, 40, 2.0)``: 20 planted promoters with Delta > 0, 20 planted
    terminators with Delta < 0, all in the truth file; >= 38 of 40 true Delta have the planted
    sign (the within-operon t noise may flip a pair)."""
    s = sim("discordance_power", (N, 40, 2.0), 0)
    tr = s.truth["discordant_pairs"]
    assert len(tr) == 40 and tr["pair_id"].is_unique
    prom = tr[tr["element"] == "promoter"]
    term = tr[tr["element"] == "terminator"]
    assert len(prom) == 20 and len(term) == 20
    assert set(prom["direction"]) == {"downstream_higher"}
    assert set(term["direction"]) == {"upstream_higher"}
    assert np.allclose(prom["effect_log2"].astype(float), 2.0)
    assert np.allclose(term["effect_log2"].astype(float), -2.0)
    pp = s.pairs.set_index("pair_id")
    for part, element in ((prom, "promoter"), (term, "terminator")):
        rows = pp.loc[list(part["pair_id"])]
        assert (rows["class"] == "OP").all() and _b(rows["true_operon"]).all()
        assert _b(rows[f"{element}_between_true"]).all()
        other = "terminator" if element == "promoter" else "promoter"
        assert not _b(rows[f"{other}_between_true"]).any()
    d = true_delta(s, pp.loc[list(tr["pair_id"])].reset_index())
    sign = np.sign(tr["effect_log2"].astype(float).to_numpy())
    assert int((np.sign(d) == sign).sum()) >= 38
    measurable = both(pp.loc[list(tr["pair_id"])])
    assert len(measurable) >= 24  # >= 60 % with both genes detected
    assert s.truth["missed_operons"].empty


def test_example4_only_planted_op_pairs_carry_internal_elements() -> None:
    s = sim("discordance_power", (N, 40, 2.0), 1)
    planted = set(s.truth["discordant_pairs"]["pair_id"])
    op = s.pairs[_b(s.pairs["true_operon"])]
    inner = _b(op["promoter_between_true"]) | _b(op["terminator_between_true"])
    assert set(op.loc[inner, "pair_id"]) == planted


def test_example5_ar1_attack_morans_i_and_no_planted_effect() -> None:
    """``operon_attack_ar1(1800, 0.9)``: lag-1 Moran's I of chr1 TU abundance >= 0.5 (5 seeds);
    the truth has no planted co-abundant or discordant pair."""
    for seed in range(5):
        s = sim("operon_attack_ar1", (N, 0.9), seed)
        x = eo.tu_series(s.layout, s.expression.genes)["chr1"]
        assert oo.morans_i_lag1(x) >= 0.5, seed
        assert s.truth["missed_operons"].empty
        assert s.truth["discordant_pairs"].empty
    assert _call("operon_attack_ar1", N, 0.9)[1].ar1_rho == 0.9


def _block_genes(s: Any, kind: str) -> set[str]:
    ids = set(s.layout.regions.loc[s.layout.regions["region_type"] == kind, "region_id"])
    rg = s.layout.region_genes
    return set(rg.loc[rg["region_id"].isin(ids), "locus_tag"].astype(str))


def test_example6_detection_attack_confounds_detection_and_distance() -> None:
    """``operon_attack_detection(1800, "prophage")``: the detected share inside the prophage block
    is <= half the share outside; the median gap between consecutive block genes is below the
    chr1 median gap (3 seeds)."""
    lcfg, _ = _call("operon_attack_detection", N, "prophage")
    assert lcfg.short_genes_in == "prophage"
    for seed in range(3):
        s = sim("operon_attack_detection", (N, "prophage"), seed)
        block = _block_genes(s, "prophage")
        assert len(block) >= 20
        ex = s.expression.genes
        inside = ex["locus_tag"].astype(str).isin(block).to_numpy()
        det = ex["detected_true"].astype(bool).to_numpy()
        assert det[inside].mean() <= 0.5 * det[~inside].mean(), seed
        g = s.layout.genes[s.layout.genes["contig"] == "chr1"].sort_values("start")
        tags = g["locus_tag"].astype(str).to_numpy()
        gaps = g["start"].to_numpy()[1:] - g["end"].to_numpy()[:-1]
        in_b = np.array([t in block for t in tags])
        both_in = in_b[1:] & in_b[:-1]
        assert both_in.sum() >= 10
        assert np.median(gaps[both_in]) < np.median(gaps), seed


def _frames(s: Any) -> list[pd.DataFrame]:
    return [s.layout.genes, s.layout.tus, s.expression.genes, s.pairs,
            *(s.truth[k] for k in TRUTH_KEYS)]  # fmt: skip


@pytest.mark.parametrize(
    "name,args",
    [("operon_null", (600,)), ("missed_operon_power", (600, 10, 0.1)),
     ("discordance_power", (600, 10, 2.0)), ("operon_attack_ar1", (600, 0.9)),
     ("operon_attack_detection", (600, "prophage"))],
)  # fmt: skip
def test_example7_same_seed_identical_outputs(
    name: str, args: tuple[Any, ...], tmp_path: Path
) -> None:
    from giapsim.contract_tables import write_operon_truth

    m = _sc()
    one = m.simulate_operons(_call(name, *args), 11)
    two = m.simulate_operons(_call(name, *args), 11)
    for x, y in zip(_frames(one), _frames(two), strict=True):
        pd.testing.assert_frame_equal(x, y)
    p1 = write_operon_truth(one, tmp_path / "a")
    p2 = write_operon_truth(two, tmp_path / "b")
    assert set(p1) == set(p2)
    for k in p1:
        assert Path(p1[k]).read_bytes() == Path(p2[k]).read_bytes(), k
    other = m.simulate_operons(_call(name, *args), 12)
    assert not other.expression.genes.equals(one.expression.genes)


# ----------------------------------------------------------------------- items 6, 7, decision 1


def test_decision1_operon_null_plants_no_internal_elements() -> None:
    """Lead decision 1: no true promoter or terminator between members of one TU; the elements
    of §11.3 sit at TU boundaries (most true NO pairs have a promoter or a terminator between)."""
    for seed in range(3):
        p = sim("operon_null", (N,), seed).pairs
        op = p[_b(p["true_operon"])]
        assert len(op) > 0
        assert not _b(op["promoter_between_true"]).any(), seed
        assert not _b(op["terminator_between_true"]).any(), seed
        no = p[(p["class"].isin(["OP", "NO"])).to_numpy() & ~_b(p["true_operon"])]
        boundary = _b(no["promoter_between_true"]) | _b(no["terminator_between_true"])
        assert boundary.mean() >= 0.5, seed


def test_item6_uniop_errors_are_recorded_separately() -> None:
    """The planted truth (``true_operon``) is independent of the UniOP class: the default layout
    gives true-operon pairs on both sides of 0.5 (Beta(8, 2) has ~2 % mass below 0.5)."""
    p = sim("missed_operon_power", (N, 50, 0.1), 0).pairs
    same = p[p["class"].isin(["OP", "NO"])]
    tru = _b(same["true_operon"])
    assert tru.any() and (~tru).any()
    assert ((same["class"] == "OP").to_numpy() & tru).sum() > 0
    assert ((same["class"] == "NO").to_numpy() & ~tru).sum() > 0
    assert ((same["class"] == "NO").to_numpy() & tru).sum() + (
        (same["class"] == "OP").to_numpy() & ~tru).sum() > 0  # fmt: skip


def test_item7_one_detected_op_pairs_exist() -> None:
    """Peptide-level truth for M6-02: OP pairs with one gene undetected and its partner detected
    (>= 20 in ``operon_null(1800)``)."""
    p = op_true(sim("operon_null", (N,), 0).pairs)
    one = _b(p["detected_a"]) ^ _b(p["detected_b"])
    assert one.sum() >= 20


# ----------------------------------------------------------------------------- contract


def test_truth_types_are_declared() -> None:
    from giapsim.contract_tables import CONTRACT_TYPES, OPERON_TRUTH_TYPES

    assert set(OPERON_TRUTH_TYPES) >= set(TRUTH_KEYS)
    assert not set(OPERON_TRUTH_TYPES) & set(CONTRACT_TYPES)  # truth, not giap contract tables
    for key, cols in TRUTH_COLUMNS.items():
        names = tuple(c for c, _ in OPERON_TRUTH_TYPES[key])
        assert names[: len(cols)] == cols, key
        for _, dtype in OPERON_TRUTH_TYPES[key]:
            assert dtype.rstrip("?") in {"str", "int", "float", "bool", "list"}


@pytest.mark.parametrize(
    "name,args",
    [("operon_null", (600,)), ("missed_operon_power", (600, 10, 0.1)),
     ("discordance_power", (600, 10, 2.0))],
)  # fmt: skip
def test_truth_files_round_trip(name: str, args: tuple[Any, ...], tmp_path: Path) -> None:
    from giapsim.contract_tables import OPERON_TRUTH_TYPES, read_operon_truth, write_operon_truth

    s = sim(name, args, 5)
    paths = write_operon_truth(s, tmp_path)
    assert set(paths) == set(OPERON_TRUTH_TYPES)
    for key, path in paths.items():
        assert Path(path) == tmp_path / "truth" / f"{key}.tsv"
        header = Path(path).read_text(encoding="utf-8").splitlines()[0].split("\t")
        assert header == [c for c, _ in OPERON_TRUTH_TYPES[key]], key
    back = read_operon_truth(tmp_path)
    assert set(back) == set(OPERON_TRUTH_TYPES)
    for key in OPERON_TRUTH_TYPES:
        a = s.truth[key].reset_index(drop=True)
        b = back[key].reset_index(drop=True)
        assert list(b.columns) == list(a.columns), key
        pd.testing.assert_frame_equal(a, b, check_dtype=False, check_exact=False, rtol=1e-5)
    text = (tmp_path / "truth" / "pairs.tsv").read_text(encoding="utf-8")
    assert "\ttrue" in text or "\tfalse" in text  # giap TSV booleans


def test_missed_operons_file_lists_the_planted_ids(tmp_path: Path) -> None:
    from giapsim.contract_tables import write_operon_truth

    s = sim("missed_operon_power", (N, 50, 0.1), 2)
    path = write_operon_truth(s, tmp_path)["missed_operons"]
    df = pd.read_csv(path, sep="\t")
    assert len(df) == 50
    assert set(df["pair_id"]) == set(s.truth["missed_operons"]["pair_id"])


# ----------------------------------------------------------------------------- degenerate cases


def test_too_many_planted_no_pairs_names_both_counts() -> None:
    with pytest.raises(ValueError) as e:
        _sc().simulate_operons(_call("missed_operon_power", 60, 500, 0.1), 0)
    msg = str(e.value)
    assert "500" in msg and "NO" in msg
    others = [int(t) for t in "".join(c if c.isdigit() else " " for c in msg).split()]
    assert any(n < 500 for n in others), msg  # the available count


def test_too_many_planted_op_pairs_names_both_counts() -> None:
    """The OP analogue for ``discordance_power`` (not a spec row; same behaviour assumed)."""
    with pytest.raises(ValueError) as e:
        _sc().simulate_operons(_call("discordance_power", 60, 500, 2.0), 0)
    msg = str(e.value)
    assert "500" in msg and "OP" in msg


@pytest.mark.parametrize(
    "name,args",
    [("operon_null", (60,)), ("missed_operon_power", (60, 1, 0.1)),
     ("discordance_power", (60, 2, 2.0)), ("operon_attack_ar1", (60, 0.5)),
     ("operon_attack_detection", (200, "prophage"))],
)  # fmt: skip
def test_small_genome_allowed_and_documented(name: str, args: tuple[Any, ...]) -> None:
    """Too small for 20 OP pairs: allowed; each docstring says the suite needs >= 1,000 genes."""
    doc = getattr(_sc(), name).__doc__ or ""
    assert "1,000" in doc or "1000" in doc, name
    s = sim(name, args, 0)
    assert len(s.layout.genes) == args[0]


@pytest.mark.parametrize("seed", [0, 1])
def test_rho0_ar1_attack_is_operon_null(seed: int) -> None:
    assert _call("operon_attack_ar1", 600, 0.0) == _call("operon_null", 600)
    a = sim("operon_attack_ar1", (600, 0.0), seed)
    b = sim("operon_null", (600,), seed)
    for x, y in zip(_frames(a), _frames(b), strict=True):
        pd.testing.assert_frame_equal(x, y)


def test_no_convergent_pair_is_a_value_error() -> None:
    lcfg, ecfg = _call("operon_null", 600)
    one_strand = dataclasses.replace(lcfg, p_same_strand=1.0)
    with pytest.raises(ValueError, match=r"(?i)convergent"):
        _sc().simulate_operons((one_strand, ecfg), 0)


# ----------------------------------------------------------------------------- properties


@settings(max_examples=8, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(seed=st.integers(0, 2**31 - 1), n=st.integers(200, 500),
       name=st.sampled_from(["operon_null", "operon_attack_ar1"]))  # fmt: skip
def test_property_pair_invariants(seed: int, n: int, name: str) -> None:
    """Any seed: OP iff UniOP p >= 0.5 on same-strand pairs, ``true_operon`` only on OP/NO,
    unique pair ids, no planted internal element in a null, truth ``pairs`` == ``pairs``."""
    args: tuple[Any, ...] = (n,) if name == "operon_null" else (n, 0.5)
    s = _sc().simulate_operons(_call(name, *args), seed)
    p = s.pairs
    assert p["pair_id"].is_unique
    same = p["class"].isin(["OP", "NO"]).to_numpy()
    u = p["uniop_p"].astype(float).to_numpy()
    assert ((p["class"] == "OP").to_numpy()[same] == (u[same] >= 0.5)).all()
    assert not _b(p["true_operon"])[~same].any()
    tru = _b(p["true_operon"])
    assert not (_b(p["promoter_between_true"]) & tru).any()
    assert not (_b(p["terminator_between_true"]) & tru).any()
    assert set(s.truth["pairs"]["pair_id"]) == set(p["pair_id"])
    assert s.truth["missed_operons"].empty and s.truth["discordant_pairs"].empty
