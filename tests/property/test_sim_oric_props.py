"""Hypothesis properties of giapsim's oriC, replichore truth and ori-ter gradient (issue M5-04).

For any chromosome (100-400 genes, optional plasmid / linear contig, origin crossing on or off),
``oric_frac`` in [0, 1) and seed: one 232-bp ``oric`` row on chr1, centred in the nearest chr1 gap
it fits in, overlapping no gene, with 0 <= start < L and end <= 2 L. For any circle, oriC and gene
midpoints: the replichore truth equals the spec formula (tests/helpers/oric_oracle.py), positions in
[0, 1]. For any gradient size and shape: a paired run differs from ``oriter_log2 = 0`` by exactly
g(p) of the gene's TU (no extra draws).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import expression_oracle as eo
import oric_oracle as oo

SETTINGS = settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow])


@st.composite
def layout_kw(draw: st.DrawFn) -> dict[str, Any]:
    return {
        "n_genes_chromosome": draw(st.integers(100, 400)),
        "plasmids": tuple(draw(st.lists(st.integers(1, 20), max_size=1))),
        "linear_contigs": tuple(draw(st.lists(st.integers(1, 20), max_size=1))),
        "origin_crossing": draw(st.booleans()),
        "oric_frac": draw(st.floats(0.0, 1.0, exclude_max=True, allow_nan=False)),
    }


def _layout(kw: dict[str, Any], seed: int) -> Any:
    from giapsim.layout import LayoutConfig, simulate_layout

    return simulate_layout(LayoutConfig(**kw), np.random.default_rng(seed))


@SETTINGS
@given(kw=layout_kw(), seed=st.integers(0, 2**32 - 1))
def test_oric_placement_invariants(kw: dict[str, Any], seed: int) -> None:
    lay = _layout(kw, seed)
    f = lay.features
    assert tuple(f.columns) == oo.FEATURE_COLUMNS
    assert f["type"].tolist() == ["oric"] and f["contig"].tolist() == ["chr1"]
    assert f["strand"].tolist() == ["."]
    assert f["feature_id"].tolist() == ["chr1_oric1"]
    n = int(lay.replicons.set_index("contig").loc["chr1", "length"])
    s, e = int(f["start"].iloc[0]), int(f["end"].iloc[0])
    assert e - s == oo.ORIC_LEN and 0 <= s < n and e <= 2 * n
    g = lay.genes[lay.genes["contig"] == "chr1"]
    for gs, ge in zip(g["start"], g["end"], strict=True):
        assert not oo.overlaps_circular((s, e), (int(gs), int(ge)), n)
    gaps = oo.chr1_gaps(lay.genes, n)
    home = [(a, b) for a, b in gaps for k in (-1, 0, 1) if a <= s + k * n and e + k * n <= b]
    assert home, (s, e)
    a, b = home[0]
    assert oo.circ_dist((s + e) / 2, (a + b) / 2, n) <= 1.0
    target = kw["oric_frac"] * n
    best = min(oo.circ_dist((x + y) / 2, target, n) for x, y in gaps if y - x >= oo.ORIC_LEN)
    assert oo.circ_dist((a + b) / 2, target, n) <= best + 1.0


Circle = tuple[int, tuple[int, int], list[tuple[int, int]], int]


@st.composite
def circles(draw: st.DrawFn) -> Circle:
    length = draw(st.integers(2_000, 200_000))
    o_start = draw(st.integers(0, length - 1))
    mids = draw(st.lists(st.integers(10, length - 11), min_size=1, max_size=12))
    cross = draw(st.integers(0, 9))  # half-length of an origin-crossing gene (0 = none)
    return length, (o_start, o_start + oo.ORIC_LEN), [(m - 10, m + 10) for m in mids], cross


@SETTINGS
@given(c=circles())
def test_replichore_truth_formula(c: Circle) -> None:
    from giapsim.expression import ExpressionConfig, simulate_expression

    length, oric, genes, cross = c
    if cross:
        genes = [*genes, (length - cross, length + cross)]  # midpoint exactly at the origin
    lay = oo.hand_layout(length, genes, oric)
    g = simulate_expression(lay, ExpressionConfig(), np.random.default_rng(0)).genes
    o = oo.midpoint(*oric, length)
    for (s, e), arm, p in zip(genes, g["replichore_true"], g["replichore_position_true"],
                              strict=True):  # fmt: skip
        want_arm, want_p = oo.replichore(oo.midpoint(s, e, length), o, length)
        assert (arm, round(float(p), 9)) == (want_arm, round(want_p, 9))
        assert 0.0 <= float(p) <= 1.0


@SETTINGS
@given(
    seed=st.integers(0, 2**32 - 1),
    frac=st.floats(0.0, 1.0, exclude_max=True, allow_nan=False),
    oriter=st.floats(-4.0, 4.0, allow_nan=False).filter(lambda v: v != 0.0),
    shape=st.sampled_from(["cosine", "linear"]),
    rho=st.sampled_from([0.0, 0.5, 0.9]),
)
def test_gradient_paired_difference(
    seed: int, frac: float, oriter: float, shape: str, rho: float
) -> None:
    from giapsim.expression import ExpressionConfig, simulate_expression

    lay = _layout({"n_genes_chromosome": 150, "plasmids": (10,), "oric_frac": frac}, seed)
    a = simulate_expression(lay, ExpressionConfig(ar1_rho=rho), np.random.default_rng(seed))
    b = simulate_expression(
        lay, ExpressionConfig(ar1_rho=rho, oriter_log2=oriter, oriter_shape=shape),
        np.random.default_rng(seed),
    )  # fmt: skip
    _, _, o = oo.oric_of(lay) or (0, 0, 0.0)
    n = int(lay.replicons.set_index("contig").loc["chr1", "length"])
    tus = lay.tus.set_index("tu_id")
    want = {}
    for tid, members in eo.tu_members(lay).items():
        row = tus.loc[tid]
        shift = 0.0
        if row["contig"] == "chr1":
            p = oo.replichore(oo.midpoint(int(row["start"]), int(row["end"]), n), o, n)[1]
            shift = oo.gradient(p, oriter, shape)
        want.update(dict.fromkeys(members, shift))
    diff = b.genes.set_index("locus_tag")["tu_log2"] - a.genes.set_index("locus_tag")["tu_log2"]
    np.testing.assert_allclose(diff.to_numpy(float), pd.Series(want)[diff.index].to_numpy(float),
                               atol=1e-9)  # fmt: skip
    pd.testing.assert_series_equal(a.genes["detected_true"], b.genes["detected_true"])
