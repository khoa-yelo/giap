"""Hypothesis properties of the gippasim expression model (issue M4-08; GIPPA_SPEC.md §11.3).

For any small layout and any valid ``ExpressionConfig`` (rho in [0, 1), plasmid copy numbers,
region abundance and detection shifts, islands, both detection modes): one row per layout gene;
``p_detect_true`` is exactly the spec logistic (oracle, tests/helpers/expression_oracle.py);
``tu_log2`` is one value per TU; ``in_island`` marks exactly the planted TUs; planted shifts are
additive (same seed with and without them); the frames are a pure function of (layout, config,
seed); ``simple_covariates`` equals the M4-10 definitions on any protein string.
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import expression_oracle as eo

SETTINGS = settings(max_examples=25, deadline=None,
                    suppress_health_check=[HealthCheck.too_slow])  # fmt: skip


#: IslandSpec object id -> (first TU index, n_tus), so tests need not know the field names.
ISLAND_TUS: dict[int, tuple[int, int]] = {}


def _mods() -> tuple[Any, Any]:
    import gippasim.expression as ex
    import gippasim.layout as la

    return ex, la


@st.composite
def scenarios(draw: st.DrawFn) -> tuple[Any, dict[str, Any], int]:
    ex, la = _mods()
    n = draw(st.integers(120, 400))
    plasmids = tuple(draw(st.lists(st.integers(1, 40), max_size=2)))
    linear = tuple(draw(st.lists(st.integers(1, 15), max_size=1)))
    lay = la.simulate_layout(
        la.LayoutConfig(n_genes_chromosome=n, plasmids=plasmids, linear_contigs=linear,
                        origin_crossing=draw(st.booleans())),
        np.random.default_rng(draw(st.integers(0, 2**31))),
    )  # fmt: skip
    shift = st.floats(-4, 4, allow_nan=False)
    kw: dict[str, Any] = {
        "ar1_rho": draw(st.floats(0.0, 0.95)),
        "plasmid_copy_log2": {f"p{i + 1}": draw(shift) for i in range(len(plasmids))},
        "region_effect_log2": draw(st.dictionaries(st.sampled_from(["prophage", "bgc", "crispr",
                                                                    "amr"]), shift, max_size=2)),
        "region_detect_shift": draw(st.dictionaries(st.sampled_from(["prophage", "bgc", "amr"]),
                                                    st.floats(-2, 2), max_size=2)),
        "detection": draw(st.sampled_from(["covariates_only", "abundance"])),
    }  # fmt: skip
    n_tus = len(eo.tu_order(lay)["chr1"])
    if draw(st.booleans()) and n_tus >= 5:
        k = draw(st.integers(5, min(20, n_tus)))
        first = draw(st.integers(0, n_tus - k))
        kw["islands"] = (ex.IslandSpec("chr1", first, k, draw(st.floats(-2, 2))),)
        ISLAND_TUS[id(kw["islands"][0])] = (first, k)
    if kw["detection"] == "abundance":
        beta = dict(ex.DEFAULT_BETA)
        beta["log2_abundance"] = draw(st.floats(0.0, 1.0))
        beta["intercept"] -= beta["log2_abundance"] * 20.0
        kw["beta"] = beta
    return lay, kw, draw(st.integers(0, 2**31))


def _run(lay: Any, kw: dict[str, Any], seed: int) -> Any:
    ex, _ = _mods()
    return ex.simulate_expression(lay, ex.ExpressionConfig(**kw), np.random.default_rng(seed))


@SETTINGS
@given(scenarios())
def test_rows_probabilities_and_logit(sc: tuple[Any, dict[str, Any], int]) -> None:
    lay, kw, seed = sc
    ex, _ = _mods()
    e = _run(lay, kw, seed)
    m = eo.merged(lay, e.genes)
    p = m["p_detect_true"].to_numpy(float)
    assert np.isfinite(p).all() and ((p > 0) & (p < 1)).all()
    assert m["detected_true"].isin([True, False]).all()
    has = ~m["aa"].map(eo.is_na)
    beta = kw.get("beta", ex.DEFAULT_BETA)
    x = eo.expected_logit(m[has], beta, kw["detection"], kw["region_detect_shift"],
                          eo.region_types_of_genes(lay))  # fmt: skip
    np.testing.assert_allclose(p[has.to_numpy()], eo.sigmoid(x), atol=1e-9)


@SETTINGS
@given(scenarios())
def test_tu_level_and_islands(sc: tuple[Any, dict[str, Any], int]) -> None:
    lay, kw, seed = sc
    e = _run(lay, kw, seed)
    g = e.genes.set_index("locus_tag")
    mem = eo.tu_members(lay)
    for members in mem.values():
        assert g.loc[members, "tu_log2"].nunique() == 1
    want: set[str] = set()
    for spec in kw.get("islands", ()):
        first, k = ISLAND_TUS[id(spec)]
        for t in eo.tu_order(lay)["chr1"][first : first + k]:
            want.update(mem[t])
    assert set(g.index[g["in_island"].astype(bool)]) == want
    assert len(e.islands) == len(kw.get("islands", ()))


@SETTINGS
@given(scenarios())
def test_planted_shifts_are_additive(sc: tuple[Any, dict[str, Any], int]) -> None:
    """Plasmid copy and region effects: same seed with vs without = the planted shift on genes of
    TUs entirely inside (or outside) the shifted set."""
    lay, kw, seed = sc
    plain = {k: v for k, v in kw.items() if k not in ("plasmid_copy_log2", "region_effect_log2")}
    a = _run(lay, plain, seed).genes.set_index("locus_tag")
    b = _run(lay, kw, seed).genes.set_index("locus_tag")
    rt = eo.region_types_of_genes(lay)
    contig = dict(zip(lay.genes["locus_tag"], lay.genes["contig"], strict=True))
    for members in eo.tu_members(lay).values():
        types = [rt.get(g, set()) & set(kw["region_effect_log2"]) for g in members]
        if any(t != types[0] for t in types):
            continue  # TU straddles a region boundary: not pinned by the spec
        want = kw["plasmid_copy_log2"].get(contig[members[0]], 0.0)
        want += sum(kw["region_effect_log2"][t] for t in types[0])
        for gname in members:
            assert math.isclose(b.at[gname, "tu_log2"] - a.at[gname, "tu_log2"], want,
                                abs_tol=1e-9)  # fmt: skip


@SETTINGS
@given(scenarios())
def test_pure_function_of_layout_config_seed(sc: tuple[Any, dict[str, Any], int]) -> None:
    lay, kw, seed = sc
    a, b = _run(lay, kw, seed), _run(lay, kw, seed)
    pd.testing.assert_frame_equal(a.genes, b.genes)
    pd.testing.assert_frame_equal(a.islands, b.islands)


@settings(max_examples=200, deadline=None)
@given(st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=1, max_size=400))
def test_simple_covariates_match_m4_10_definitions(aa: str) -> None:
    ex, _ = _mods()
    n, gr = ex.simple_covariates(aa)
    assert n == eo.n_obs(aa)
    assert math.isclose(gr, eo.gravy(aa), abs_tol=1e-9)
    assert 0 <= n <= len(aa) // 7


@settings(max_examples=30, deadline=None)
@given(st.floats(0.0, 1.0, exclude_max=True))
def test_any_rho_in_unit_interval_accepted(rho: float) -> None:
    ex, la = _mods()
    lay = la.simulate_layout(la.LayoutConfig(n_genes_chromosome=60, regions=()),
                             np.random.default_rng(1))  # fmt: skip
    e = ex.simulate_expression(lay, ex.ExpressionConfig(ar1_rho=rho), np.random.default_rng(1))
    assert np.isfinite(e.genes["tu_log2"].to_numpy(float)).all()
