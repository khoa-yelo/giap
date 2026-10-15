"""Hypothesis properties of M4-14 gene-level abundance (docs/specs/M4-14-abundance.md).

Acceptance criterion "Property": scaling every intensity by c > 0 shifts every non-imputed
``abundance_log2`` by log2(c); iBAQ <= sum of intensities; imputed genes are a subset of detected
genes; ``exclude`` never sets ``abundance_imputed``. Plus: agreement with the independent oracle
(tests/helpers/abundance_oracle.py) for ``peptide_abundance`` and ``impute``, Top3 bounds,
invariance to the row order of ``accepted``, and the same scale shift on the linear protein route.
Conventions as in test_abundance_gene.py.
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
import abundance_oracle as ao
from abundance_oracle import acc, accepted_frame, c10, pmap_frame
from test_abundance_gene import TAGS, cfg, covariates_frame, detection, genes_frame, na

SETTINGS = settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
POLICY = st.sampled_from(ao.POLICIES)
INTENSITY = st.one_of(
    st.floats(min_value=1e-3, max_value=1e12, allow_nan=False, allow_infinity=False),
    st.just(0.0),
    st.none(),
)


def _ab() -> Any:
    from giap.genetable import abundance

    return abundance


@st.composite
def scenarios(draw: st.DrawFn) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for t in TAGS:
        n = draw(st.integers(min_value=0, max_value=6))
        for j in range(n):
            uniq = draw(st.sampled_from(["unique", "unique", "unique", "multi_identical"]))
            rows.append(acc(t, f"P{t}X{j}K", draw(INTENSITY), uniqueness=uniq, start=30 * j))
    n_obs = {t: draw(st.integers(min_value=0, max_value=8)) for t in TAGS}
    detected = tuple(t for t in TAGS if draw(st.booleans()))
    partial = tuple(t for t in TAGS if t not in detected and draw(st.booleans()))
    return {
        "rows": rows,
        "n_obs": n_obs,
        "detected": detected,
        "partial": partial,
        "policy": draw(POLICY),
        "method": draw(st.sampled_from(["ibaq", "top3"])),
    }


def _run(sc: dict[str, Any], rows: list[dict[str, Any]] | None = None) -> Any:
    rows = sc["rows"] if rows is None else rows
    return _ab().gene_abundance(
        genes_frame(TAGS, sc["partial"]),
        accepted_frame(rows),
        covariates_frame(sc["n_obs"], TAGS),
        detection(TAGS, sc["detected"], sc["partial"]),
        None,
        pmap_frame([]),
        cfg(method=sc["method"]),
        policy=sc["policy"],
    )


def _scaled(rows: list[dict[str, Any]], c: float) -> list[dict[str, Any]]:
    out = []
    for r in rows:
        r2 = dict(r)
        if not ao.is_na(r2["intensity"]):
            r2["intensity"] = r2["intensity"] * c
        out.append(r2)
    return out


@SETTINGS
@given(scenarios(), st.floats(min_value=1e-3, max_value=1e3))
def test_scaling_intensities_shifts_log2_by_log2_c(sc: dict[str, Any], c: float) -> None:
    a = _run(sc).table.set_index("locus_tag")
    b = _run(sc, _scaled(sc["rows"], c)).table.set_index("locus_tag")
    shift = math.log2(c)
    for t in TAGS:
        if na(a.at[t, "abundance_log2"]) or bool(a.at[t, "abundance_imputed"]):
            continue
        assert not bool(b.at[t, "abundance_imputed"])
        assert math.isclose(
            float(b.at[t, "abundance_log2"]), float(a.at[t, "abundance_log2"]) + shift,
            rel_tol=0, abs_tol=1e-9,
        )  # fmt: skip
        if not na(a.at[t, "top3_log2"]):
            assert math.isclose(
                float(b.at[t, "top3_log2"]), float(a.at[t, "top3_log2"]) + shift, abs_tol=1e-9
            )


@SETTINGS
@given(scenarios())
def test_imputed_subset_of_detected_and_exclude_never_imputes(sc: dict[str, Any]) -> None:
    res = _run(sc)
    t = res.table
    imputed = t["abundance_imputed"].fillna(False).astype(bool)
    imputed_tags = set(t.loc[imputed, "locus_tag"])
    assert imputed_tags <= set(sc["detected"])
    if sc["policy"] == "exclude":
        assert not imputed_tags
        assert res.fill_value is None
    for _, r in t.iterrows():
        f = set() if ao.is_na(r["abundance_flags"]) else {str(x) for x in r["abundance_flags"]}
        assert ("imputed_abundance" in f) == (r["locus_tag"] in imputed_tags)
        if r["locus_tag"] in imputed_tags:
            assert float(r["abundance_log2"]) == res.fill_value


@SETTINGS
@given(scenarios())
def test_gene_abundance_matches_oracle(sc: dict[str, Any]) -> None:
    res = _run(sc)
    t = res.table.set_index("locus_tag")
    vals = ao.peptide_values(accepted_frame(sc["rows"]))
    if not vals:
        return  # source "none" or nothing quantified: covered by the example tests
    raw = {
        g: (ao.ibaq(v, sc["n_obs"][g]) if sc["method"] == "ibaq" else ao.top3(v))
        for g, v in vals.items()
    }
    log2 = {g: (math.log2(raw[g]) if g in raw else math.nan) for g in TAGS}
    det = {g: g in sc["detected"] for g in TAGS}
    filled, mask, fill = ao.impute(log2, det, sc["policy"])
    assert (res.fill_value is None) == (fill is None)
    if fill is not None:
        assert math.isclose(res.fill_value, fill, rel_tol=1e-12, abs_tol=1e-9)
    for g in TAGS:
        got = t.at[g, "abundance_log2"]
        exp = filled[g]
        assert (na(got) and na(exp)) or math.isclose(float(got), exp, abs_tol=1e-9), (g, got, exp)
        assert bool(t.at[g, "abundance_imputed"]) == mask[g]
        assert t.at[g, "imputation_policy"] == sc["policy"]
        lt3 = g in vals and len(vals[g]) < 3
        f = t.at[g, "abundance_flags"]
        fl = set() if ao.is_na(f) else {str(x) for x in f}
        assert ("top3_lt3" in fl) == lt3


@SETTINGS
@given(scenarios())
def test_peptide_abundance_matches_oracle_and_bounds(sc: dict[str, Any]) -> None:
    accepted = accepted_frame(sc["rows"])
    out = _ab().peptide_abundance(accepted, pd.Series(sc["n_obs"]))
    vals = ao.peptide_values(accepted)
    by = out.set_index("locus_tag")
    for g, v in vals.items():
        ibaq, top3 = float(by.at[g, "ibaq"]), float(by.at[g, "top3"])
        assert math.isclose(ibaq, ao.ibaq(v, sc["n_obs"][g]), rel_tol=1e-12)
        assert math.isclose(top3, ao.top3(v), rel_tol=1e-12)
        assert int(by.at[g, "n_quant_peptides"]) == len(v)
        assert ibaq <= math.fsum(v) * (1 + 1e-12)
        top = sorted(v, reverse=True)[:3]
        assert min(top) * (1 - 1e-12) <= top3 <= max(top) * (1 + 1e-12)
    for g in set(by.index) - set(vals):
        n = by.at[g, "n_quant_peptides"]
        assert na(n) or int(n) == 0
        assert na(by.at[g, "ibaq"]) and na(by.at[g, "top3"])


@SETTINGS
@given(scenarios(), st.randoms(use_true_random=False))
def test_row_order_invariance(sc: dict[str, Any], rnd: Any) -> None:
    rows = list(sc["rows"])
    rnd.shuffle(rows)
    pd.testing.assert_frame_equal(_run(sc).table, _run(sc, rows).table)


@SETTINGS
@given(
    st.lists(
        st.one_of(st.none(), st.floats(min_value=-50, max_value=50, allow_nan=False)),
        min_size=1,
        max_size=30,
    ),
    st.lists(st.booleans(), min_size=30, max_size=30),
    POLICY,
)
def test_impute_matches_oracle(values: list[float | None], det: list[bool], policy: str) -> None:
    idx = [f"G{i:03d}" for i in range(len(values))]
    log2 = pd.Series([math.nan if v is None else v for v in values], index=idx, dtype="float64")
    detected = pd.Series(det[: len(values)], index=idx)
    before = log2.copy()
    filled, mask, fill = _ab().impute(log2, detected, policy)
    o_filled, o_mask, o_fill = ao.impute(log2.to_dict(), detected.to_dict(), policy)
    assert (fill is None) == (o_fill is None)
    if fill is not None:
        assert math.isclose(fill, o_fill, abs_tol=1e-9)  # type: ignore[arg-type]
    filled = pd.Series(filled, index=idx)
    mask = pd.Series(mask, index=idx).astype(bool)
    assert mask.to_dict() == o_mask
    for g in idx:
        assert (na(filled[g]) and na(o_filled[g])) or math.isclose(
            float(filled[g]), float(o_filled[g]), abs_tol=1e-9
        )
    assert not (mask & ~detected).any()
    pd.testing.assert_series_equal(log2, before)  # input not mutated


@SETTINGS
@given(
    st.dictionaries(
        st.sampled_from(["P1", "P2", "P3", "P4", "P5"]),
        st.floats(min_value=1e-3, max_value=1e12),
        min_size=1,
    ),
    st.floats(min_value=1e-3, max_value=1e3),
)
def test_protein_route_linear_scale_shift(values: dict[str, float], c: float) -> None:
    owners = {"P1": TAGS[0], "P2": TAGS[0], "P3": TAGS[1], "P4": TAGS[2], "P5": TAGS[3]}
    pm = pmap_frame([(p, owners[p], "protein_hit") for p in values])
    ab = _ab()
    a = ab.protein_abundance(c10(values), pm, "maxlfq", "linear").set_index("locus_tag")
    b = ab.protein_abundance(c10({k: v * c for k, v in values.items()}), pm, "maxlfq", "linear")
    b = b.set_index("locus_tag")
    for g in {owners[p] for p in values}:
        exp = max(math.log2(values[p]) for p in values if owners[p] == g)
        assert math.isclose(float(a.at[g, "value_log2"]), exp, abs_tol=1e-9)
        assert math.isclose(float(b.at[g, "value_log2"]), exp + math.log2(c), abs_tol=1e-9)
        n = sum(owners[p] == g for p in values)
        assert int(a.at[g, "n_proteins"]) == n
        assert ("two_proteins_one_gene" in [str(f) for f in a.at[g, "flags"]]) == (n >= 2)


def test_no_nan_leak_into_numpy_types() -> None:
    # a quick non-Hypothesis anchor: abundance_log2 is a float column
    sc = {"rows": [acc(TAGS[0], "PEPK", 8.0)], "n_obs": dict.fromkeys(TAGS, 1),
          "detected": (TAGS[0],), "partial": (), "policy": "median", "method": "ibaq"}  # fmt: skip
    t = _run(sc).table
    assert np.issubdtype(t["abundance_log2"].dtype, np.floating)
    assert np.issubdtype(t["top3_log2"].dtype, np.floating)
