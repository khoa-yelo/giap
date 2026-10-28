"""M4-16 basic detection summary (docs/specs/M4-16-detection-summary.md, examples 1-5).

Spec: GIPPA_SPEC.md §7.3.5(a) (n_genes, n_detected, % detected, E = sum p_detect, O - E with a 95 %
TU-block bootstrap CI, O/E per replicon, feature class and COG category), §8.2 ("resample TUs, not
genes"), §9.3 (long format), §6.6 (``gippa.rng.for_step``).

Resolved ambiguities (pinned here):
- Every basic group also reports a metric row ``n_tu`` (distinct TUs of the group's eligible
  genes) and the ``o_minus_e`` row carries the flag ``ci_few_units`` when ``n_tu < 10``
  (docs/methods/bootstrap.md; M4-04 forwarded item, lead brief) or min(E, n_genes - E) < 10
  (PR #87 lead decision 4; tests in test_detection_summary_review.py), and not otherwise.
- Rows other than ``o_minus_e`` / ``o_over_e`` carry NA ``ci_low`` / ``ci_high``.
- ``flags`` may be a ``;``-joined string or a list in memory; compared as a set.
- Groups are emitted only for categories with >= 1 eligible gene; a COG group label is the letter,
  ``none`` for genes without a category.
- Single-peptide genes (``single_peptide_unestimated``, ``detected = false``) are eligible and
  undetected (primary analysis).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detection_summary_fixture as fx
from detection_summary_fixture import gene, metric

B = 1000
SEED = 20261015


def _summary() -> Any:
    from gippa.context import summary

    return summary


def _rng(seed: int = SEED) -> np.random.Generator:
    from gippa.rng import for_step

    return for_step("context_stats", seed)


def _basic(genes: pd.DataFrame, b: int = B, seed: int = SEED) -> pd.DataFrame:
    return _summary().basic_rows(genes, b, _rng(seed))  # type: ignore[no-any-return]


def _ex1(**kw: Any) -> pd.DataFrame:
    d = [True, False, True, True]
    p = [0.9, 0.2, 0.5, 0.4]
    tus = ["A", "A", "B", "C"]
    return fx.genes_frame(gene(f"G{i + 1}", tus[i], d[i], p[i], **kw) for i in range(4))


# --------------------------------------------------------------------------- interface


def test_summary_columns_constant() -> None:
    assert tuple(_summary().SUMMARY_COLUMNS) == fx.SUMMARY_COLUMNS


def test_basic_rows_columns_and_dtypes() -> None:
    out = _basic(_ex1())
    assert tuple(out.columns) == fx.SUMMARY_COLUMNS
    assert set(out["section"].astype(str)) == {"basic"}
    for c in ("value", "ci_low", "ci_high"):
        assert pd.api.types.is_float_dtype(out[c]), c
    assert pd.api.types.is_integer_dtype(out["n"]), out["n"].dtype
    for c in ("section", "group_type", "group", "metric"):
        assert all(isinstance(v, str) for v in out[c]), c
    assert set(out["group_type"].astype(str)) <= {"replicon", "feature_class", "cog_category"}


def test_metric_set_per_group() -> None:
    out = _basic(_ex1())
    for (gt, g), sub in out.groupby(["group_type", "group"]):
        assert sorted(sub["metric"].astype(str)) == sorted([*fx.BASIC_METRICS, "n_tu"]), (gt, g)


# --------------------------------------------------------------------------- example 1


@pytest.mark.parametrize(
    ("group_type", "group"),
    [("replicon", "chr1"), ("feature_class", "cds"), ("cog_category", "none")],
)
def test_example_1_point_estimates(group_type: str, group: str) -> None:
    out = _basic(_ex1())
    want = {"n_genes": 4, "n_detected": 3, "pct_detected": 75.0, "expected": 2.0,
            "o_minus_e": 1.0, "o_over_e": 1.5, "n_tu": 3}  # fmt: skip
    for name, v in want.items():
        assert metric(out, group_type, group, name)["value"] == pytest.approx(v, abs=1e-12), name
    for name in [*fx.BASIC_METRICS, "n_tu"]:
        assert int(metric(out, group_type, group, name)["n"]) == 4, name


def test_example_1_groups_present() -> None:
    out = _basic(_ex1())
    keys = set(zip(out["group_type"].astype(str), out["group"].astype(str), strict=True))
    assert keys == {("replicon", "chr1"), ("feature_class", "cds"), ("cog_category", "none")}


def test_count_rows_have_no_ci() -> None:
    out = _basic(_ex1())
    for name in ("n_genes", "n_detected", "pct_detected", "expected", "n_tu"):
        r = metric(out, "replicon", "chr1", name)
        assert math.isnan(r["ci_low"]) and math.isnan(r["ci_high"]), name


# --------------------------------------------------------------------------- example 2


def test_example_2_ci_within_attainable_range() -> None:
    r = metric(_basic(_ex1()), "replicon", "chr1", "o_minus_e")
    assert r["value"] == pytest.approx(1.0)
    lo, hi = float(r["ci_low"]), float(r["ci_high"])
    assert -0.3 - 1e-9 <= lo <= hi <= 1.8 + 1e-9, (lo, hi)
    assert lo < hi  # three distinct TU sums: the bootstrap is not degenerate


def test_example_2_matches_tu_block_bootstrap() -> None:
    """The basic CI is the M4-04 TU-block bootstrap of the group sum of d - p (one group)."""
    from gippa.stats.bootstrap import tu_block_bootstrap

    genes = _ex1()
    out = _basic(genes)
    r = metric(out, "replicon", "chr1", "o_minus_e")
    # same sum statistic, B and level; the CI must lie in the same attainable support
    ci = tu_block_bootstrap(
        genes.assign(v=genes["detected"].astype(float) - genes["p_detect"]),
        "tu_id", lambda d: float(d["v"].sum()), B, np.random.default_rng(1),
    )  # fmt: skip
    assert ci.estimate == pytest.approx(r["value"])
    # 1000 replicates of 3 draws: both percentile CIs are close (MC error only)
    assert abs(ci.ci_low - r["ci_low"]) <= 0.35 and abs(ci.ci_high - r["ci_high"]) <= 0.35


def test_example_2_few_units_flag_and_n_tu() -> None:
    out = _basic(_ex1())
    assert metric(out, "replicon", "chr1", "n_tu")["value"] == 3
    assert fx.FEW_UNITS_FLAG in fx.flag_set(metric(out, "replicon", "chr1", "o_minus_e")["flags"])


def test_no_few_units_flag_with_ten_tus() -> None:
    # 10 TUs of 4 genes, p = 0.5: E = n - E = 20, so only the TU count can set the flag
    # (PR #87 lead decision 4 added min(E, n - E) < 10; aligned from 10 single-gene TUs).
    rows = [gene(f"G{i}", f"TU{i // 4}", i % 2 == 0, 0.5) for i in range(40)]
    out = _basic(fx.genes_frame(rows), b=200)
    r = metric(out, "replicon", "chr1", "o_minus_e")
    assert metric(out, "replicon", "chr1", "n_tu")["value"] == 10
    assert fx.FEW_UNITS_FLAG not in fx.flag_set(r["flags"])
    rows9 = rows[:36]
    r9 = metric(_basic(fx.genes_frame(rows9), b=200), "replicon", "chr1", "o_minus_e")
    assert fx.FEW_UNITS_FLAG in fx.flag_set(r9["flags"])


def test_single_tu_group_ci_is_degenerate() -> None:
    rows = [gene("G1", "A", True, 0.3), gene("G2", "A", False, 0.6), gene("G3", "A", True, 0.8)]
    r = metric(_basic(fx.genes_frame(rows), b=200), "replicon", "chr1", "o_minus_e")
    est = 2 - 1.7
    assert r["value"] == pytest.approx(est)
    assert r["ci_low"] == pytest.approx(est) and r["ci_high"] == pytest.approx(est)


def test_only_group_genes_inside_a_drawn_tu_count() -> None:
    """TU A holds a K gene and a T gene: K's single-TU CI is the K gene's d - p only."""
    rows = [
        gene("G1", "A", True, 0.25, cog_category="K"),
        gene("G2", "A", False, 0.9, cog_category="T"),
        gene("G3", "B", True, 0.5, cog_category="T"),
    ]
    out = _basic(fx.genes_frame(rows), b=300)
    k = metric(out, "cog_category", "K", "o_minus_e")
    assert k["value"] == pytest.approx(0.75)
    assert k["ci_low"] == pytest.approx(0.75) and k["ci_high"] == pytest.approx(0.75)
    t = metric(out, "cog_category", "T", "o_minus_e")
    # T TU sums: A -> -0.9, B -> 0.5; replicates are sums of two draws from {-0.9, 0.5}
    assert t["value"] == pytest.approx(-0.4)
    assert -1.8 - 1e-9 <= t["ci_low"] <= t["ci_high"] <= 1.0 + 1e-9


def test_ci_is_tu_block_not_gene_level() -> None:
    """One TU of 6 identical genes plus 6 singleton TUs: TU resampling gives the attainable
    support of sums of 7 TU draws; check the CI width against the TU-block bootstrap."""
    rows = [gene(f"A{i}", "A", True, 0.0) for i in range(6)]
    rows += [gene(f"S{i}", f"S{i}", False, 0.0) for i in range(6)]
    out = _basic(fx.genes_frame(rows), b=1000)
    r = metric(out, "replicon", "chr1", "o_minus_e")
    # TU-block: O - E = 6 * Binomial(7, 1/7) per replicate; 95 % CI about [0, 18].
    # Gene-level resampling would give 12 draws with success 1/2: about [3, 9].
    assert r["ci_low"] <= 0.0 + 1e-9
    assert r["ci_high"] >= 12.0


# --------------------------------------------------------------------------- example 3 (COG)


def test_example_3_multi_letter_cog_counts_in_each_letter() -> None:
    rows = [
        gene("G1", "A", True, 0.5, cog_category="KT"),
        gene("G2", "B", False, 0.5, cog_category="K"),
        gene("G3", "C", True, 0.5),
    ]
    out = _basic(fx.genes_frame(rows), b=100)
    assert metric(out, "cog_category", "K", "n_genes")["value"] == 2
    assert metric(out, "cog_category", "T", "n_genes")["value"] == 1
    assert metric(out, "cog_category", "T", "n_detected")["value"] == 1
    assert metric(out, "cog_category", "none", "n_genes")["value"] == 1
    assert fx.rows_of(out, "cog_category", "KT").empty


def test_no_cog_annotation_only_none_rows() -> None:
    out = _basic(_ex1())
    cog = out[out["group_type"].astype(str) == "cog_category"]
    assert set(cog["group"].astype(str)) == {"none"}


# --------------------------------------------------------------------------- example 4 (E = 0)


def test_example_4_zero_expected_gives_na_ratio() -> None:
    rows = [gene("G1", "A", True, 0.0), gene("G2", "B", False, 0.0)]
    out = _basic(fx.genes_frame(rows), b=100)
    assert metric(out, "replicon", "chr1", "expected")["value"] == 0.0
    assert math.isnan(metric(out, "replicon", "chr1", "o_over_e")["value"])
    assert metric(out, "replicon", "chr1", "o_minus_e")["value"] == pytest.approx(1.0)


# --------------------------------------------------------------------------- example 5 (partial)


def test_example_5_partial_genes_excluded_everywhere() -> None:
    rows = list(_ex1().to_dict("records"))
    rows.append(gene("P1", "D", False, None, status="partial_excluded", cog_category="K",
                     hypothetical=True))  # fmt: skip
    rows.append(gene("P2", "E", False, None, contig="ctg9", status="partial_excluded"))
    out = _basic(fx.genes_frame(rows))
    assert metric(out, "replicon", "chr1", "n_genes")["value"] == 4
    assert metric(out, "replicon", "chr1", "expected")["value"] == pytest.approx(2.0)
    assert metric(out, "feature_class", "cds", "n_genes")["value"] == 4
    assert fx.rows_of(out, "cog_category", "K").empty
    assert fx.rows_of(out, "feature_class", "hypothetical").empty
    assert fx.rows_of(out, "replicon", "ctg9").empty
    assert not out["value"].isna().all()
    assert np.isfinite(metric(out, "replicon", "chr1", "o_minus_e")["ci_low"])


def test_single_peptide_genes_are_eligible_and_undetected() -> None:
    rows = list(_ex1().to_dict("records"))
    rows.append(gene("S1", "D", False, 0.3, status="single_peptide_unestimated"))
    out = _basic(fx.genes_frame(rows))
    assert metric(out, "replicon", "chr1", "n_genes")["value"] == 5
    assert metric(out, "replicon", "chr1", "n_detected")["value"] == 3
    assert metric(out, "replicon", "chr1", "expected")["value"] == pytest.approx(2.3)


# --------------------------------------------------------------------------- feature classes


def test_feature_classes_and_hypothetical_any_type() -> None:
    rows = [
        gene("G1", "A", True, 0.5, feature_type="cds", hypothetical=True),
        gene("G2", "B", False, 0.5, feature_type="sorf", hypothetical=True),
        gene("G3", "C", True, 0.5, feature_type="pseudogene"),
        gene("G4", "D", True, 0.5, feature_type="cds", contig="p1"),
    ]
    out = _basic(fx.genes_frame(rows), b=100)
    n = {g: metric(out, "feature_class", g, "n_genes")["value"] for g in fx.FEATURE_CLASSES}
    assert n == {"cds": 2, "sorf": 1, "pseudogene": 1, "hypothetical": 2}
    assert metric(out, "replicon", "chr1", "n_genes")["value"] == 3
    assert metric(out, "replicon", "p1", "n_genes")["value"] == 1


def test_pct_detected_is_percentage() -> None:
    rows = [gene(f"G{i}", f"T{i}", i < 1, 0.5) for i in range(8)]
    out = _basic(fx.genes_frame(rows), b=50)
    assert metric(out, "replicon", "chr1", "pct_detected")["value"] == pytest.approx(12.5)


# --------------------------------------------------------------------------- determinism


def test_deterministic_with_same_seed() -> None:
    a = _basic(_ex1(), seed=3)
    b = _basic(_ex1(), seed=3)
    pd.testing.assert_frame_equal(a, b)


def test_different_seed_changes_only_cis() -> None:
    rows = [gene(f"G{i}", f"T{i // 2}", i % 3 == 0, 0.1 * (i % 9)) for i in range(40)]
    g = fx.genes_frame(rows)
    a, b = _basic(g, b=200, seed=1), _basic(g, b=200, seed=2)
    key = ["group_type", "group", "metric"]
    a, b = a.sort_values(key).reset_index(drop=True), b.sort_values(key).reset_index(drop=True)
    pd.testing.assert_series_equal(a["value"], b["value"])
    oe = a["metric"] == "o_minus_e"
    assert not np.allclose(a.loc[oe, "ci_low"], b.loc[oe, "ci_low"])


def test_point_estimates_independent_of_row_order() -> None:
    rows = [gene(f"G{i}", f"T{i // 3}", i % 2 == 0, 0.05 * (i % 19),
                 cog_category=["K", "T", None, "KT"][i % 4]) for i in range(30)]  # fmt: skip
    g = fx.genes_frame(rows)
    shuffled = g.sample(frac=1.0, random_state=5).reset_index(drop=True)
    key = ["group_type", "group", "metric"]
    a = _basic(g, b=50).sort_values(key).reset_index(drop=True)
    b = _basic(shuffled, b=50).sort_values(key).reset_index(drop=True)
    pd.testing.assert_frame_equal(a[[*key, "value", "n"]], b[[*key, "value", "n"]])


def test_input_not_mutated() -> None:
    g = _ex1()
    copy = g.copy(deep=True)
    _basic(g, b=50)
    pd.testing.assert_frame_equal(g, copy)


def test_matches_oracle_on_mixed_frame() -> None:
    rows = [gene(f"G{i}", f"T{i // 2}", i % 3 != 0, 0.03 * (i % 30),
                 contig="chr1" if i < 25 else "p1",
                 feature_type=["cds", "cds", "sorf", "pseudogene"][i % 4],
                 hypothetical=i % 5 == 0,
                 cog_category=["K", "T", None, "KT", "J"][i % 5]) for i in range(36)]  # fmt: skip
    g = fx.genes_frame(rows)
    out = _basic(g, b=50)
    groups = fx.group_members(g)
    keys = set(zip(out["group_type"].astype(str), out["group"].astype(str), strict=True))
    assert keys == set(groups)
    for (gt, grp), mask in groups.items():
        want = fx.point_estimates(g, mask)
        for name, v in want.items():
            got = metric(out, gt, grp, name)["value"]
            if math.isnan(v):
                assert math.isnan(got), (gt, grp, name)
            else:
                assert got == pytest.approx(v, abs=1e-9), (gt, grp, name)
