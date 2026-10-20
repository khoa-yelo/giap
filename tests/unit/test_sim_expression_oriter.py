"""giapsim replichore truth and the ori-ter gradient (issue M5-04; GIAP_SPEC.md §7.3.5(e), §11.3).

Examples 3-6 and the expression half of the determinism criterion of
docs/specs/M5-04-sim-oric-gradient.md (the spec's ``tests/sim/test_expression_oriter.py``; lives in
tests/unit like the M4-08 tests). Expected values come from tests/helpers/oric_oracle.py.

Interpretations pinned here (reported to the lead):
- ``replichore_true`` / ``replichore_position_true`` are columns of ``Expression.genes`` (the spec
  lists them under expression.py), from each gene's midpoint (mod L) and the oriC midpoint;
- the gradient is a TU-level shift from the TU midpoint (``tus`` start/end, unwrapped, mod L); it
  is additive, so a paired run (same layout and seed) differs from ``oriter_log2 = 0`` by exactly
  g(p) in ``tu_log2`` and ``abundance_log2_true`` and in nothing else under ``covariates_only``;
- replicons other than chr1 get no shift; an unknown ``oriter_shape`` is a ``ValueError``.
"""

from __future__ import annotations

import hashlib
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import expression_oracle as eo
import oric_oracle as oo

MU = 20.0
NEW_COLUMNS = ("replichore_true", "replichore_position_true")
# ``simulate_expression`` gene frames before M5-04 (main at ecd3f6f), as CSV sha256:
# default = LayoutConfig() seed 1, ExpressionConfig() seed 1; mixed = the M4-08 mixed layout
# (400 + plasmids 25, 1 + linear 12, 1, origin crossing) seed 1 with MIXED_EXPR, seed 2.
PRE_M5_GENES_SHA256 = {
    "default": "32963ce5efe2d24f2a22f3488fd40e9762956824f43de793ba7c558bb1d547a3",
    "mixed": "fb1c2426d96b141ab71c5476e95b19ac2b921dbbe4f5fb8a0222fa39324579b2",
}
PRE_M5_MIXED_ISLANDS_SHA256 = "f63cdebe931aa7da6d05cadb80d91a3f811c4fefca3b25a007171120944895f1"
MIXED_KW: dict[str, Any] = {
    "n_genes_chromosome": 400, "plasmids": (25, 1), "linear_contigs": (12, 1),
    "origin_crossing": True,
}  # fmt: skip
FLAT: dict[str, Any] = {"tu_sd": 0.0, "gene_sd": 0.0, "stoich_sd": 0.0, "position_decay_log2": 0.0}


def _ex() -> Any:
    import giapsim.expression

    return giapsim.expression


def layout(seed: int = 1, **kw: Any) -> Any:
    from giapsim.layout import LayoutConfig, simulate_layout

    return simulate_layout(LayoutConfig(**kw), np.random.default_rng(seed))


def express(lay: Any, seed: int = 1, **kw: Any) -> Any:
    m = _ex()
    return m.simulate_expression(lay, m.ExpressionConfig(**kw), np.random.default_rng(seed))


def mixed_expr() -> dict[str, Any]:
    return {"ar1_rho": 0.9, "detection": "abundance",
            "islands": (_ex().IslandSpec("chr1", 10, 5, 1.5),),
            "region_effect_log2": {"prophage": -2.0}}  # fmt: skip


def sha(df: pd.DataFrame) -> str:
    return hashlib.sha256(df.to_csv(index=False).encode()).hexdigest()


def by_tag(expr: Any) -> pd.DataFrame:
    return expr.genes.set_index("locus_tag")


def tu_positions(lay: Any) -> dict[str, float]:
    """TU id -> replichore position of the TU midpoint (chr1 only; needs an oriC)."""
    oc = oo.oric_of(lay)
    assert oc is not None
    o = oc[2]
    n = int(lay.replicons.set_index("contig").loc["chr1", "length"])
    out = {}
    for tid, c, s, e in zip(lay.tus["tu_id"], lay.tus["contig"], lay.tus["start"], lay.tus["end"],
                            strict=True):  # fmt: skip
        if c == "chr1":
            out[str(tid)] = oo.replichore(oo.midpoint(int(s), int(e), n), o, n)[1]
    return out


def expected_shift(lay: Any, oriter: float, shape: str) -> pd.Series:
    """Per gene (locus tag): g(p) of its TU; 0 off chr1."""
    pos = tu_positions(lay)
    out = {}
    for tid, gs in eo.tu_members(lay).items():
        g = oo.gradient(pos[tid], oriter, shape) if tid in pos else 0.0
        for t in gs:
            out[t] = g
    return pd.Series(out)


# ----------------------------------------------------------------------------- interface


def test_config_defaults() -> None:
    c = _ex().ExpressionConfig()
    assert c.oriter_log2 == 0.0
    assert c.oriter_shape == "cosine"
    hash(_ex().ExpressionConfig(oriter_log2=2.0, oriter_shape="linear"))


def test_unknown_shape_is_value_error() -> None:
    with pytest.raises(ValueError, match="oriter_shape"):
        _ex().ExpressionConfig(oriter_shape="quadratic")


def test_gene_frame_gains_replichore_columns() -> None:
    lay = layout(1, n_genes_chromosome=200, oric_frac=0.4)
    g = express(lay, 1).genes
    for c in NEW_COLUMNS:
        assert c in g.columns, c
    assert set(g["replichore_true"].dropna()) <= {"right", "left"}
    assert pd.api.types.is_float_dtype(g["replichore_position_true"])


# ----------------------------------------------------------------------------- examples 3, 4


def test_example3_replichore_truth() -> None:
    """o = 900 on L = 1,000; midpoints 100, 400, 500, 650 -> (right, 0.4), (right, 1.0),
    (left, 0.8), (left, 0.5)."""
    lay = oo.hand_layout(1000, [(90, 110), (390, 410), (490, 510), (640, 660)], (784, 1016))
    g = express(lay, 1).genes
    assert g["replichore_true"].tolist() == ["right", "right", "left", "left"]
    np.testing.assert_allclose(g["replichore_position_true"].to_numpy(float), [0.4, 1.0, 0.8, 0.5],
                               atol=1e-12)  # fmt: skip


def test_example4_gene_crossing_the_origin() -> None:
    """Start 990, end 1,010 (unwrapped) on L = 1,000 with o = 0: midpoint 0 -> (right, 0.0)."""
    lay = oo.hand_layout(1000, [(240, 260), (740, 760), (990, 1010)], (884, 1116))
    g = by_tag(express(lay, 1))
    assert bool(lay.genes["crosses_origin"].iloc[2])
    tags = lay.genes["locus_tag"].tolist()
    assert (g.loc[tags[2], "replichore_true"], g.loc[tags[2], "replichore_position_true"]) == (
        "right", 0.0,
    )  # fmt: skip
    assert g.loc[tags[0], "replichore_true"] == "right"
    assert g.loc[tags[0], "replichore_position_true"] == pytest.approx(0.5)
    assert g.loc[tags[1], "replichore_true"] == "left"
    assert g.loc[tags[1], "replichore_position_true"] == pytest.approx(0.5)


def test_replichore_na_on_plasmids_and_without_oric() -> None:
    lay = oo.hand_layout(1000, [(90, 110), (390, 410)], (784, 1016),
                         plasmid=(500, [(10, 40), (200, 230)]))  # fmt: skip
    g = express(lay, 1).genes
    on_p = (lay.genes["contig"] == "p1").to_numpy()
    assert g.loc[on_p, "replichore_true"].isna().all()
    assert g.loc[on_p, "replichore_position_true"].isna().all()
    assert g.loc[~on_p, "replichore_true"].notna().all()
    none = express(layout(1, n_genes_chromosome=200), 1).genes
    assert none["replichore_true"].isna().all()
    assert none["replichore_position_true"].isna().all()


@pytest.mark.parametrize("frac", [0.0, 0.37, 0.8])
def test_replichore_truth_matches_oracle_on_simulated_layouts(frac: float) -> None:
    lay = layout(4, n_genes_chromosome=500, plasmids=(30,), origin_crossing=True, oric_frac=frac)
    g = by_tag(express(lay, 4))
    _, _, o = oo.oric_of(lay) or (0, 0, math.nan)
    n = int(lay.replicons.set_index("contig").loc["chr1", "length"])
    for t, c, s, e in zip(lay.genes["locus_tag"], lay.genes["contig"], lay.genes["start"],
                          lay.genes["end"], strict=True):  # fmt: skip
        if c != "chr1":
            assert pd.isna(g.loc[t, "replichore_true"]), t
            continue
        arm, p = oo.replichore(oo.midpoint(int(s), int(e), n), o, n)
        assert g.loc[t, "replichore_true"] == arm, t
        assert g.loc[t, "replichore_position_true"] == pytest.approx(p, abs=1e-12), t
    pos = g["replichore_position_true"].dropna()
    assert ((pos >= 0) & (pos <= 1)).all()
    arms = g["replichore_true"].dropna()
    assert 0.3 < (arms == "right").mean() < 0.7  # the oriC splits chr1 into two halves


# ----------------------------------------------------------------------------- example 5


def _example5_layout() -> Any:
    """L = 2,000, o = 1,000; one-gene TUs at p = 0, 0.25, 0.5 (right), 0.5 (left) and 1."""
    genes = [(990, 1010), (1240, 1260), (1490, 1510), (490, 510), (1990, 2010)]
    return oo.hand_layout(2000, genes, (884, 1116))


@pytest.mark.parametrize(
    ("shape", "want"),
    [("cosine", [1.0, math.sqrt(0.5), 0.0, 0.0, -1.0]), ("linear", [1.0, 0.5, 0.0, 0.0, -1.0])],
)
def test_example5_tu_shift_at_known_positions(shape: str, want: list[float]) -> None:
    """``oriter_log2 = 2``, no noise, no decay: p = 0 -> mu + 1, p = 1 -> mu - 1, p = 0.5 -> mu."""
    lay = _example5_layout()
    g = express(lay, 1, oriter_log2=2.0, oriter_shape=shape, **FLAT).genes
    np.testing.assert_allclose(g["tu_log2"].to_numpy(float), MU + np.array(want), atol=1e-9)
    np.testing.assert_allclose(g["abundance_log2_true"].to_numpy(float), MU + np.array(want),
                               atol=1e-9)  # fmt: skip
    np.testing.assert_allclose(g["replichore_position_true"].to_numpy(float),
                               [0.0, 0.25, 0.5, 0.5, 1.0], atol=1e-12)  # fmt: skip


@pytest.mark.parametrize("shape", ["cosine", "linear"])
def test_example5_on_a_simulated_layout(shape: str) -> None:
    lay = layout(2, n_genes_chromosome=600, plasmids=(40,), linear_contigs=(20,), oric_frac=0.3)
    g = by_tag(express(lay, 2, oriter_log2=2.0, oriter_shape=shape, **FLAT))
    want = MU + expected_shift(lay, 2.0, shape)
    np.testing.assert_allclose(g.loc[want.index, "tu_log2"].to_numpy(float), want.to_numpy(),
                               atol=1e-9)  # fmt: skip
    off = lay.genes.loc[lay.genes["contig"] != "chr1", "locus_tag"]
    assert (g.loc[off, "tu_log2"] == MU).all()  # no gradient off chr1
    on = g.loc[lay.genes.loc[lay.genes["contig"] == "chr1", "locus_tag"]]
    assert on["tu_log2"].max() == pytest.approx(MU + 1.0, abs=0.05)
    assert on["tu_log2"].min() == pytest.approx(MU - 1.0, abs=0.05)


def test_example5_sign_genes_near_ori_higher() -> None:
    lay = layout(3, n_genes_chromosome=800, oric_frac=0.8)
    g = express(lay, 3, oriter_log2=4.0).genes
    r = oo.pearson(g["replichore_position_true"], g["abundance_log2_true"])
    assert r < -0.2


@pytest.mark.parametrize("shape", ["cosine", "linear"])
def test_gradient_is_additive_and_draws_nothing(shape: str) -> None:
    """Paired runs (same layout, same seed, AR(1), noise, islands and region effects on) differ by
    exactly g(p) in ``tu_log2`` and ``abundance_log2_true`` and nowhere else."""
    lay = layout(5, n_genes_chromosome=500, plasmids=(30,), oric_frac=0.55)
    kw: dict[str, Any] = {"ar1_rho": 0.7, "islands": (_ex().IslandSpec("chr1", 40, 6, 2.0),),
                          "region_effect_log2": {"prophage": -1.5}}  # fmt: skip
    a = by_tag(express(lay, 6, **kw))
    b = by_tag(express(lay, 6, oriter_log2=1.5, oriter_shape=shape, **kw))
    shift = expected_shift(lay, 1.5, shape).reindex(a.index)
    for c in ("tu_log2", "abundance_log2_true"):
        np.testing.assert_allclose((b[c] - a[c]).to_numpy(float), shift.to_numpy(float),
                                   atol=1e-9, err_msg=c)  # fmt: skip
    rest = [c for c in a.columns if c not in ("tu_log2", "abundance_log2_true")]
    pd.testing.assert_frame_equal(a[rest], b[rest])


def test_gradient_under_abundance_detection_changes_only_p_detect() -> None:
    lay = layout(5, n_genes_chromosome=300, oric_frac=0.55)
    beta = dict(_ex().DEFAULT_BETA)
    beta.update(log2_abundance=0.3, intercept=beta["intercept"] - 0.3 * MU)
    a = express(lay, 6, detection="abundance", beta=beta).genes
    b = express(lay, 6, detection="abundance", beta=beta, oriter_log2=2.0).genes
    keep = [c for c in a.columns if c not in ("tu_log2", "abundance_log2_true", "p_detect_true",
                                              "detected_true")]  # fmt: skip
    pd.testing.assert_frame_equal(a[keep], b[keep])
    assert not np.allclose(a["p_detect_true"], b["p_detect_true"])


# ----------------------------------------------------------------------------- example 6


@pytest.mark.parametrize("oriter", [1.0, -1.0, 1e-9])
def test_example6_gradient_without_oric_is_value_error(oriter: float) -> None:
    lay = layout(1, n_genes_chromosome=200)
    with pytest.raises(ValueError, match=r"(?i)oric"):
        express(lay, 1, oriter_log2=oriter)


def test_zero_gradient_without_oric_is_fine() -> None:
    express(layout(1, n_genes_chromosome=200), 1, oriter_log2=0.0)


# ----------------------------------------------------------------------------- byte identity


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="pinned float formatting")
@pytest.mark.parametrize("name", ["default", "mixed"])
def test_oriter_zero_keeps_pre_m5_expression(name: str) -> None:
    """``oriter_log2 = 0`` (default): every pre-M5 column is byte-identical to before M5-04."""
    if name == "default":
        ex = express(layout(1), 1)
    else:
        ex = express(layout(1, **MIXED_KW), 2, **mixed_expr())
        assert sha(ex.islands) == PRE_M5_MIXED_ISLANDS_SHA256
    old = [c for c in ex.genes.columns if c not in NEW_COLUMNS]
    assert sha(ex.genes[old]) == PRE_M5_GENES_SHA256[name]


def test_oric_with_zero_gradient_changes_only_replichore_columns() -> None:
    a = express(layout(7, n_genes_chromosome=400), 7)
    b = express(layout(7, n_genes_chromosome=400, oric_frac=0.25), 7)
    old = [c for c in a.genes.columns if c not in NEW_COLUMNS]
    pd.testing.assert_frame_equal(a.genes[old], b.genes[old])
    assert b.genes["replichore_true"].notna().all()


# ----------------------------------------------------------------------------- determinism


def test_same_seed_identical_expression_frame() -> None:
    lay = layout(8, n_genes_chromosome=300, plasmids=(20,), oric_frac=0.8)
    np.random.seed(1)  # noqa: NPY002 (global numpy state must not matter)
    a = express(lay, 8, oriter_log2=1.0)
    np.random.seed(2)  # noqa: NPY002
    b = express(lay, 8, oriter_log2=1.0)
    assert sha(a.genes) == sha(b.genes)
    pd.testing.assert_frame_equal(a.genes, b.genes)


def test_layout_not_mutated_by_gradient() -> None:
    lay = layout(8, n_genes_chromosome=300, oric_frac=0.8)
    before = {k: getattr(lay, k).copy() for k in ("genes", "tus", "features")}
    express(lay, 8, oriter_log2=2.0)
    for k, v in before.items():
        pd.testing.assert_frame_equal(getattr(lay, k), v, obj=k)
