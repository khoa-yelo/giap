"""Properties of the M4-18 ``regions_stats`` rows (docs/specs/M4-18-regions-stats-table.md; lead
decisions 2-8 of 2026-10-20).

- P1 plasmid rows equal the oracles on random genomes: own O - E counts; HL = raw log2 HL of the
  detected, non-imputed genes against the ``kind = chromosome`` replicons; every row untested
  (p, q NA, ``null_model = none``, ``n_perm = 0``, ``method = replicon_contrast``);
- P2 adding c to every plasmid log2 abundance shifts HL and its CI by exactly c (same seed) and
  leaves the O - E row unchanged: a per-replicon residualised HL (statistician's mutant) fails;
- P3 adding c to every gene's log2 abundance leaves HL and its CI unchanged;
- P4 the depth-adjusted row is HL - log2(plasmid depth / length-weighted reference depth), CI
  shifted by the same constant;
- P5 prophage subset rows partition the members that carry a subset (n_genes, observed, expected
  add up) and are descriptive;
- P6 in the assembled table ``null_model = none`` <=> untested for every M4-18 row; a p-value
  implies ``null_model != none``, a q-value implies a p-value; M4-17 q-values are copied.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import region_tests_data as rd
import regions_stats_data as sd
from region_tests_data import Contig

SETTINGS = settings(max_examples=30, deadline=None, suppress_health_check=[HealthCheck.too_slow])
P_VALUES = [0.0, 0.1, 0.3, 0.5, 0.8, 0.95]


def _cfg() -> Any:
    return rd.config(bootstrap=20)


@st.composite
def contig(draw: st.DrawFn, name: str, kind: str, lo: int, hi: int) -> Contig:
    n = draw(st.integers(lo, hi))
    d = draw(st.lists(st.integers(0, 1), min_size=n, max_size=n))
    p = draw(st.lists(st.sampled_from(P_VALUES), min_size=n, max_size=n))
    a = draw(st.lists(st.integers(0, 400).map(lambda x: 5.0 + x / 16), min_size=n, max_size=n))
    cuts = draw(st.lists(st.booleans(), min_size=n, max_size=n))
    tu, k = [], 0
    for i in range(n):
        k += int(cuts[i] and i > 0)
        tu.append(f"TU_{name}_{k}")
    return Contig(name, d, p, a=a, tu=tu, kind=kind)


@st.composite
def genome_case(draw: st.DrawFn) -> dict[str, Any]:
    chrom = draw(contig("chr1", "chromosome", 2, 30))
    pl = draw(contig("p1", "plasmid", 3, 20))  # >= M4-17 SMALL_REGION: statistics computed
    det = [f"p1_g{i:05d}" for i, x in enumerate(pl.d) if x] + [
        f"chr1_g{i:05d}" for i, x in enumerate(chrom.d) if x
    ]
    imputed = draw(st.lists(st.sampled_from(det), unique=True, max_size=6)) if det else []
    return {"contigs": [chrom, pl], "imputed": imputed, "seed": draw(st.integers(0, 2**31 - 1))}


def _plasmids(case: dict[str, Any], contigs: list[Contig] | None = None,
              depth: pd.DataFrame | None = None) -> tuple[sd.Genome, pd.DataFrame]:  # fmt: skip
    g = sd.genome(contigs or case["contigs"], plasmids=["p1"], imputed=case["imputed"])
    T = sd.module()
    return g, T.plasmid_rows(g.genes, g.replicons, depth, _cfg(), sd.rng(case["seed"]))


def _same(x: Any, y: Any, tol: float = 1e-9) -> bool:
    if sd.na(x) or sd.na(y):
        return sd.na(x) and sd.na(y)
    return abs(float(x) - float(y)) <= tol


def _shifted(c: Contig, by: float) -> Contig:
    return Contig(c.name, c.d, c.p, a=[x + by for x in (c.a or [])], tu=c.tu, kind=c.kind)


@SETTINGS
@given(case=genome_case())
def test_p1_plasmid_rows_match_oracles(case: dict[str, Any]) -> None:
    g, res = _plasmids(case)
    assert sorted(res["stat"]) == ["hl_shift", "o_minus_e"]
    want = sd.detection_counts(g.genes, ["p1"])
    hl_want = sd.plasmid_hl(g.genes, "p1", ["chr1"])
    for _, row in res.iterrows():
        assert sd.na(row["p_value"]) and sd.na(row["q_value"])
        assert row["null_model"] == "none" and int(row["n_perm"]) == 0
        assert row["method"] == sd.PLASMID_METHOD
        if want["n_eligible"] > 0:
            for c in ("observed", "expected", "o_minus_e"):
                assert _same(row[c], want[c]), (row["stat"], c)
            assert _same(row["o_over_e"], want["o_over_e"])
            assert _same(row["hl_shift_log2"], hl_want)
    hl = sd.unadjusted(res, "plasmid_1", "hl_shift")
    assert _same(hl["effect"], hl_want, 1e-12)
    oe = sd.unadjusted(res, "plasmid_1", "o_minus_e")
    assert _same(oe["effect"], want["o_minus_e"])


@SETTINGS
@given(case=genome_case(), c=st.sampled_from([-3.0, 0.5, 4.0, math.log2(20.0)]))
def test_p2_plasmid_dosage_shift_is_carried_by_hl(case: dict[str, Any], c: float) -> None:
    chrom, pl = case["contigs"]
    _, a = _plasmids(case)
    _, b = _plasmids(case, [chrom, _shifted(pl, c)])
    ha, hb = sd.unadjusted(a, "plasmid_1", "hl_shift"), sd.unadjusted(b, "plasmid_1", "hl_shift")
    for col in ("effect", "ci_low", "ci_high", "hl_shift_log2"):
        x = ha[col]
        assert _same(hb[col], math.nan if sd.na(x) else float(x) + c), col
    oa, ob = sd.unadjusted(a, "plasmid_1", "o_minus_e"), sd.unadjusted(b, "plasmid_1", "o_minus_e")
    for col in ("effect", "ci_low", "ci_high", "observed", "expected"):
        assert _same(oa[col], ob[col], 0.0), col


@SETTINGS
@given(case=genome_case(), c=st.sampled_from([-7.0, 2.5, 11.0]))
def test_p3_common_shift_leaves_hl_unchanged(case: dict[str, Any], c: float) -> None:
    chrom, pl = case["contigs"]
    _, a = _plasmids(case)
    _, b = _plasmids(case, [_shifted(chrom, c), _shifted(pl, c)])
    ha, hb = sd.unadjusted(a, "plasmid_1", "hl_shift"), sd.unadjusted(b, "plasmid_1", "hl_shift")
    for col in ("effect", "ci_low", "ci_high"):
        assert _same(ha[col], hb[col]), col


@SETTINGS
@given(case=genome_case(), dc=st.floats(1.0, 500.0), dp=st.floats(1.0, 5000.0))
def test_p4_depth_offset(case: dict[str, Any], dc: float, dp: float) -> None:
    _, res = _plasmids(case, depth=sd.depth_frame({"chr1": dc, "p1": dp}))
    hl = sd.unadjusted(res, "plasmid_1", "hl_shift")
    adj = res[res["label"].astype(str) == sd.ADJUSTED_LABEL]
    if sd.na(hl["effect"]):  # lead decision 13: no adjusted row when HL is NA
        assert adj.empty
        return
    assert len(adj) == 1
    assert sd.DEPTH_FLAG in sd.flags(adj.iloc[0]["flags"])
    row = adj.iloc[0]
    off = math.log2(dp / dc)
    for col in ("effect", "ci_low", "ci_high"):
        x = hl[col]
        assert _same(row[col], math.nan if sd.na(x) else float(x) - off), col
    assert sd.na(row["p_value"]) and row["null_model"] == "none" and int(row["n_perm"]) == 0


@st.composite
def prophage_case(draw: st.DrawFn) -> dict[str, Any]:
    n = draw(st.integers(8, 40))
    d = draw(st.lists(st.integers(0, 1), min_size=n, max_size=n))
    p = draw(st.lists(st.sampled_from(P_VALUES), min_size=n, max_size=n))
    lo = draw(st.integers(0, n - 3))
    hi = draw(st.integers(lo + 2, n))
    subs = draw(st.lists(st.sampled_from([*sd.SUBSETS, None]), min_size=hi - lo, max_size=hi - lo))
    return {"contig": Contig("chr1", d, p, tu=sd.pairs("chr1", n)), "lo": lo, "hi": hi,
            "subsets": subs}  # fmt: skip


@SETTINGS
@given(case=prophage_case())
def test_p5_subset_rows_partition_the_prophage(case: dict[str, Any]) -> None:
    c = case["contig"]
    tags = [f"chr1_g{i:05d}" for i in range(case["lo"], case["hi"])]
    sub = {t: s for t, s in zip(tags, case["subsets"], strict=True) if s is not None}
    g = sd.genome([c], {"prophage_1": ("prophage", tags)}, subsets=sub)
    res = sd.module().subset_rows(g.genes, g.regions, g.region_genes)
    assert set(res["subset"]) == set(sub.values())
    assert int(res["n_genes"].astype(int).sum()) == len(sub)
    genes = g.genes.set_index("locus_tag")
    assert float(res["observed"].astype(float).sum()) == float(
        genes.loc[list(sub), "detected"].sum()
    )
    e = float(genes.loc[list(sub), "p_detect"].astype(float).sum())
    assert float(res["expected"].astype(float).sum()) == pytest.approx(e, abs=1e-9)
    for _, row in res.iterrows():
        assert "descriptive" in sd.flags(row["flags"])
        assert row["method"] == sd.DESCRIPTIVE and row["label"] == sd.DESCRIPTIVE
        assert sd.na(row["ci_low"]) and sd.na(row["ci_high"])
        assert row["stat"] == "o_minus_e" and row["null_model"] == "none"
        assert sd.na(row["p_value"]) and sd.na(row["q_value"]) and int(row["n_perm"]) == 0


@settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(case=genome_case(), lo=st.integers(0, 10), width=st.integers(1, 12))
def test_p6_table_null_model_rule(case: dict[str, Any], lo: int, width: int) -> None:
    chrom, pl = case["contigs"]
    n = len(chrom.d)
    idx = list(range(min(lo, n - 1), min(lo + width, n)))
    tags = [f"chr1_g{i:05d}" for i in idx]
    g = sd.genome([chrom, pl], {"bgc_1": ("bgc", tags)}, plasmids=["p1"],
                  imputed=case["imputed"])  # fmt: skip
    T = sd.module()
    tests = sd.tests(g, case["seed"])
    res = T.regions_stats(
        tests,
        T.plasmid_rows(g.genes, g.replicons, sd.depth_frame({"chr1": 30.0, "p1": 90.0}), _cfg(),
                       sd.rng(case["seed"])),
        T.subset_rows(g.genes, g.regions, g.region_genes),
        g.regions,
    )  # fmt: skip
    assert tuple(res.columns) == sd.REGIONS_STATS_COLUMNS
    sd.check_row_order(res.reset_index(drop=True))
    for _, row in res.iterrows():
        if row["null_model"] == "none":
            assert sd.na(row["p_value"]) and sd.na(row["q_value"]) and int(row["n_perm"]) == 0
        if not sd.na(row["p_value"]):
            assert row["null_model"] != "none"
        if not sd.na(row["q_value"]):
            assert not sd.na(row["p_value"])
    m17 = res[res["region_id"] == "bgc_1"].sort_values("stat")
    src = tests.sort_values("stat")
    np.testing.assert_array_equal(m17["q_value"].to_numpy(float), src["q_value"].to_numpy(float))
    assert (res.loc[res["region_type"] == "plasmid", "null_model"] == "none").all()
