"""Hypothesis properties of the M4-16 basic detection summary.

- §7.3.6 property "O - E is unchanged by relabelling genes outside the region", applied to groups:
  changing ``detected`` / ``p_detect`` of genes outside a group (membership kept) leaves every row
  of that group unchanged, CIs included (same seed; the draws depend only on TU counts).
- Point estimates equal the oracle; o_minus_e = n_detected - expected; o_over_e = O / E.
- The CI lies in the attainable support: n_tu * min(TU sum) <= ci_low <= ci_high <= n_tu * max.
- Same seed, same output (determinism).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detection_summary_fixture as fx
from detection_summary_fixture import gene

B = 40
COGS = ["K", "T", None, "KT", "J"]
FTYPES = ["cds", "sorf", "pseudogene"]
STATUSES = ["detected", "not_detected", "single_peptide_unestimated", "partial_excluded"]


def _basic(genes: pd.DataFrame, seed: int = 11) -> pd.DataFrame:
    from giap.context.summary import basic_rows
    from giap.rng import for_step

    return basic_rows(genes, B, for_step("context_stats", seed))  # type: ignore[no-any-return]


@st.composite
def gene_frames(draw: Any) -> pd.DataFrame:
    n = draw(st.integers(1, 14))
    rows = []
    tu = 0
    for i in range(n):
        if i and draw(st.booleans()):
            tu += 1
        status = draw(st.sampled_from(STATUSES))
        det = status == "detected"
        p = (
            None
            if status == "partial_excluded"
            else draw(st.sampled_from([0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]))
        )
        rows.append(gene(
            f"G{i:02d}", f"TU{tu}", det, p, status=status,
            contig=draw(st.sampled_from(["chr1", "p1"])),
            feature_type=draw(st.sampled_from(FTYPES)),
            hypothetical=draw(st.booleans()),
            cog_category=draw(st.sampled_from(COGS)),
        ))  # fmt: skip
    return fx.genes_frame(rows)


def _same(a: Any, b: Any) -> bool:
    if fx.is_na(a) and fx.is_na(b):
        return True
    if isinstance(a, float) or isinstance(b, float):
        return math.isclose(float(a), float(b), rel_tol=0, abs_tol=1e-12)
    return bool(a == b)


SETTINGS = settings(max_examples=40, deadline=None,
                    suppress_health_check=[HealthCheck.too_slow])  # fmt: skip


@SETTINGS
@given(gene_frames(), st.data())
def test_relabelling_outside_group_leaves_group_rows_unchanged(g: pd.DataFrame, data: Any) -> None:
    groups = fx.group_members(g)
    if not groups:
        return
    key = data.draw(st.sampled_from(sorted(groups)))
    mask = groups[key]
    elig = (g["detection_status"] != "partial_excluded").to_numpy()
    h = g.copy()
    for i in np.flatnonzero(~mask & elig):
        det = data.draw(st.booleans())
        h.loc[h.index[i], "detected"] = det
        h.loc[h.index[i], "detection_status"] = "detected" if det else "not_detected"
        h.loc[h.index[i], "p_detect"] = data.draw(st.sampled_from([0.0, 0.3, 0.6, 1.0]))
    a = fx.rows_of(_basic(g), *key).sort_values("metric").reset_index(drop=True)
    b = fx.rows_of(_basic(h), *key).sort_values("metric").reset_index(drop=True)
    assert len(a) == len(b) > 0
    for col in fx.SUMMARY_COLUMNS:
        for x, y in zip(a[col], b[col], strict=True):
            if col == "flags":
                assert fx.flag_set(x) == fx.flag_set(y)
            else:
                assert _same(x, y), (key, col, x, y)


@SETTINGS
@given(gene_frames())
def test_point_estimates_and_ci_support(g: pd.DataFrame) -> None:
    out = _basic(g)
    groups = fx.group_members(g)
    keys = set(zip(out["group_type"].astype(str), out["group"].astype(str), strict=True))
    assert keys == set(groups)
    for key, mask in groups.items():
        want = fx.point_estimates(g, mask)
        for name, v in want.items():
            got = float(fx.metric(out, *key, name)["value"])
            assert (math.isnan(v) and math.isnan(got)) or math.isclose(got, v, abs_tol=1e-9), (
                key, name, got, v)  # fmt: skip
        r = fx.metric(out, *key, "o_minus_e")
        sums = list(fx.tu_sums(g, mask).values())
        k = len(sums)
        lo, hi = k * min(sums) - 1e-9, k * max(sums) + 1e-9
        assert lo <= float(r["ci_low"]) <= float(r["ci_high"]) <= hi, (key, r, sums)
        assert int(r["n"]) == int(mask.sum())
        e = float(g["p_detect"].to_numpy()[mask].sum())
        want_flag = fx.few_units(k, e, int(mask.sum()))  # PR #87 decision 4
        assert (fx.FEW_UNITS_FLAG in fx.flag_set(r["flags"])) == want_flag


@settings(max_examples=15, deadline=None)
@given(gene_frames(), st.integers(0, 2**31 - 1))
def test_same_seed_same_output(g: pd.DataFrame, seed: int) -> None:
    pd.testing.assert_frame_equal(_basic(g, seed), _basic(g, seed))
