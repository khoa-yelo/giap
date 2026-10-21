"""Hypothesis properties of the M5-08 ori-ter analysis (acceptance list of docs/specs/M5-08).

Invariants:

* **geometry** -- ``replichore`` positions lie in [0, 1]; 0 at the oriC and 1 at the terminus;
  a point and its mirror image across the oriC share a position on opposite arms; the label is
  ``right`` exactly when the circular distance oriC -> point (increasing coordinate) is at most
  L / 2 (spec "Replichore");
* **rotation** -- rotating every gene coordinate *and* the oriC by the same offset leaves
  (``replichore``, ``replichore_position``) unchanged for every gene, and leaves the rows
  unchanged (the geometry is intrinsic to the circle);
* **strand** -- swapping strands never matters (the midpoint does not depend on it);
* **coverage of the gene columns** -- one row per gene of the gene table, NA exactly off the
  analysed replicons (circular chromosomes with an oriC call);
* **row contract (example 8)** -- three rows per analysed replicon with an oriC (one with
  ``no_oric`` otherwise), ``p_value`` / ``q_value`` / ``lag`` / the null columns NA, flag
  ``descriptive`` always, ``not_estimable`` exactly when the value or the CI is NA, ``n_perm``
  in [0, ``stats.bootstrap``], and never the dropped ``uncalibrated`` (ADR-0021 Amendment 1);
* **estimator** -- the slope is equivariant under ``y -> c y + b`` with c > 0 (OLS, and the
  block length is a function of Moran's I, which is affine-invariant, so the draws do not
  change); with every TU on one arm the pooled and that arm's rows agree;
* **oriC** -- one call gives its midpoint; rotating the calls rotates the oriC (circular mean);
  calls more than 0.05 L apart give no oriC (lead decision 31).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import oric_oracle as oo
import oriter_data as od
from oriter_data import flag_set, is_na

SETTINGS = settings(
    max_examples=25,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)
CFG = {"bootstrap": 20, "permutations": 9}


def _oriter() -> Any:
    from giap.context import oriter

    return oriter


@st.composite
def circle(draw: Any) -> tuple[int, float, list[float]]:
    """(length, oriC midpoint, gene midpoints) on one circular replicon."""
    length = draw(st.integers(2_000, 200_000))
    oric = float(draw(st.integers(0, length - 1)))
    n = draw(st.integers(1, 12))
    mids = draw(st.lists(st.integers(0, length - 1), min_size=n, max_size=n, unique=True))
    return length, oric, [float(m) for m in mids]


# --------------------------------------------------------------------------- geometry


@SETTINGS
@given(circle())
def test_positions_are_in_the_unit_interval(spec: tuple[int, float, list[float]]) -> None:
    length, oric, mids = spec
    arm, pos = _oriter().replichore(mids, oric, length)
    p = np.asarray(pos, dtype=float)
    assert np.all((p >= 0.0) & (p <= 1.0))
    assert set(str(a) for a in arm) <= {"right", "left"}
    for a, x, m in zip(arm, p, mids, strict=True):
        assert (str(a), float(x)) == pytest.approx(oo.replichore(m, oric, float(length)),
                                                   abs=1e-12)  # fmt: skip


@SETTINGS
@given(circle())
def test_origin_is_zero_and_terminus_is_one(spec: tuple[int, float, list[float]]) -> None:
    length, oric, _ = spec
    arm, pos = _oriter().replichore([oric, (oric + length / 2.0) % length], oric, length)
    assert float(np.asarray(pos, dtype=float)[0]) == pytest.approx(0.0)
    assert str(arm[0]) == "right"
    assert float(np.asarray(pos, dtype=float)[1]) == pytest.approx(1.0)


@SETTINGS
@given(circle())
def test_mirror_images_share_a_position_on_opposite_arms(
    spec: tuple[int, float, list[float]],
) -> None:
    length, oric, mids = spec
    mirrored = [(2.0 * oric - m) % length for m in mids]
    arm_a, pos_a = _oriter().replichore(mids, oric, length)
    arm_b, pos_b = _oriter().replichore(mirrored, oric, length)
    assert np.allclose(np.asarray(pos_a, float), np.asarray(pos_b, float))
    for a, b, p in zip(arm_a, arm_b, np.asarray(pos_a, float), strict=True):
        if 0.0 < p < 1.0:  # the oriC and the terminus are their own mirror image
            assert {str(a), str(b)} == {"right", "left"}


@SETTINGS
@given(circle(), st.integers(0, 10**6))
def test_rotation_leaves_the_gene_columns_unchanged(
    spec: tuple[int, float, list[float]], shift: int
) -> None:
    """Rotating every coordinate and the oriC by the same offset changes nothing."""
    length, oric, mids = spec
    genes, _, replicons = od.placed([od.Placed("chr1", length, mids, width=2)])
    rot_mids = [(m + shift) % length for m in mids]
    rot_genes, _, rot_replicons = od.placed([od.Placed("chr1", length, rot_mids, width=2)])
    base = od.gene_columns(genes, replicons, od.oric_at("chr1", oric, width=2))
    rot = od.gene_columns(
        rot_genes, rot_replicons, od.oric_at("chr1", (oric + shift) % length, width=2)
    )
    assert [str(x) for x in base["replichore"]] == [str(x) for x in rot["replichore"]]
    assert base["replichore_position"].to_numpy(float) == pytest.approx(
        rot["replichore_position"].to_numpy(float), abs=1e-9
    )


@SETTINGS
@given(circle(), st.lists(st.sampled_from(["+", "-"]), min_size=1, max_size=12))
def test_strand_never_matters(spec: tuple[int, float, list[float]], strands: list[str]) -> None:
    length, oric, mids = spec
    genes, _tus, replicons = od.placed([od.Placed("chr1", length, mids, width=2)])
    feats = od.oric_at("chr1", oric, width=2)
    flipped = genes.copy()
    flipped["strand"] = [strands[i % len(strands)] for i in range(len(flipped))]
    a = od.gene_columns(genes, replicons, feats)
    b = od.gene_columns(flipped, replicons, feats)
    pd.testing.assert_frame_equal(a, b)


# --------------------------------------------------------------------------- the gene columns


@SETTINGS
@given(circle(), st.booleans(), st.sampled_from(["chromosome", "plasmid", "unknown"]),
       st.sampled_from(["circular", "linear"]))  # fmt: skip
def test_gene_columns_cover_every_gene_and_are_na_off_the_analysis(
    spec: tuple[int, float, list[float]], with_oric: bool, kind: str, topology: str
) -> None:
    length, oric, mids = spec
    genes, _, replicons = od.placed(
        [od.Placed("chr1", length, mids, width=2, kind=kind, topology=topology)]
    )
    feats = od.oric_at("chr1", oric, width=2) if with_oric else od.features()
    out = od.gene_columns(genes, replicons, feats)
    assert list(out.columns) == ["locus_tag", "replichore", "replichore_position"]
    assert sorted(out["locus_tag"].astype(str)) == sorted(genes["locus_tag"].astype(str))
    analysed = with_oric and kind == "chromosome" and topology == "circular"
    na = out["replichore"].map(is_na)
    assert (~na).all() if analysed else na.all()
    assert all(
        is_na(a) == is_na(p)
        for a, p in zip(out["replichore"], out["replichore_position"], strict=True)
    )


# --------------------------------------------------------------------------- the rows


@st.composite
def genome(draw: Any) -> tuple[Any, ...]:
    """A circular chromosome with TUs on both arms, optionally an oriC, plus a plasmid."""
    length = draw(st.integers(20_000, 200_000))
    oric = float(draw(st.integers(0, length - 1)))
    n = draw(st.integers(1, 14))
    pos = draw(st.lists(st.floats(0.02, 0.98), min_size=n, max_size=n, unique=True))
    side = draw(st.lists(st.booleans(), min_size=n, max_size=n))
    mids = [
        (oric + (1.0 if r else -1.0) * p * length / 2.0) % length
        for p, r in zip(pos, side, strict=True)
    ]
    rng = np.random.default_rng(draw(st.integers(0, 2**31 - 1)))
    values = list(20.0 - 3.0 * np.asarray(pos) + rng.standard_normal(n))
    detected = draw(st.lists(st.booleans(), min_size=n, max_size=n))
    specs = [
        od.Placed(
            "chr1", length, mids, values=values, detected=[int(d) for d in detected], width=2
        ),
        od.Placed("p1", 6_000, [500, 2_500, 4_500], kind="plasmid", width=2),
    ]
    genes, tus, replicons = od.placed(specs)
    feats = od.oric_at("chr1", oric, width=2) if draw(st.booleans()) else od.features()
    return genes, tus, replicons, feats


@SETTINGS
@given(genome())
def test_row_contract_holds_on_every_genome(frames: tuple[Any, ...]) -> None:
    """Example 8 and the degenerate table as a property of any genome."""
    from giap.context.autocorrelation import AUTOCORR_COLUMNS

    genes, tus, replicons, feats = frames
    cfg = od.config(**CFG)
    out = od.rows(genes, tus, replicons, feats, cfg)
    assert list(out.columns) == list(AUTOCORR_COLUMNS)
    assert set(out["contig"].astype(str)) <= {"chr1"}
    has_oric = bool(len(feats))
    assert len(out) == (3 if has_oric else 1)
    stats = sorted(out["stat"].astype(str))
    assert stats == sorted(f"{od.STATISTIC}_{a}" for a in (od.ARMS if has_oric else ("pooled",)))
    for _, row in out.iterrows():
        flags = flag_set(row["flags"])
        assert od.DESCRIPTIVE in flags
        assert "uncalibrated" not in flags
        assert is_na(row["p_value"]) and is_na(row["q_value"]) and is_na(row["lag"])
        assert is_na(row["null_mean"]) and is_na(row["null_lo"]) and is_na(row["null_hi"])
        assert str(row["method"]) == od.METHOD and str(row["null_model"]) == od.NULL_MODEL
        assert 0 <= int(row["n_perm"]) <= int(cfg.stats.bootstrap)
        if is_na(row["value"]) or is_na(row["ci_low"]):
            assert flags & {od.NOT_ESTIMABLE, od.NO_ORIC}
        else:
            assert flags == {od.DESCRIPTIVE}
            assert float(row["effect"]) == pytest.approx(float(row["value"]))
            assert float(row["ci_low"]) <= float(row["ci_high"])


@SETTINGS
@given(genome())
def test_rows_match_the_oracle(frames: tuple[Any, ...]) -> None:
    genes, tus, replicons, feats = frames
    cfg = od.config(**CFG)
    got = od.rows(genes, tus, replicons, feats, cfg)
    want = od.oracle_rows(genes, tus, replicons, feats, cfg)
    pd.testing.assert_frame_equal(
        got.reset_index(drop=True), want.reset_index(drop=True), check_dtype=False
    )


@SETTINGS
@given(genome(), st.floats(0.25, 4.0), st.floats(-5.0, 5.0))
def test_slope_is_equivariant_under_affine_maps_of_abundance(
    frames: tuple[Any, ...], scale: float, shift: float
) -> None:
    """y -> c y + b (c > 0) scales every slope and CI by c and changes no flag."""
    genes, tus, replicons, feats = frames
    cfg = od.config(**CFG)
    scaled = genes.copy()
    scaled["abundance_log2"] = scale * genes["abundance_log2"].to_numpy(float) + shift
    a = od.rows(genes, tus, replicons, feats, cfg)
    b = od.rows(scaled, tus, replicons, feats, cfg)
    assert list(a["flags"]) == list(b["flags"])
    assert list(a["n_perm"]) == list(b["n_perm"])
    for col in ("value", "effect", "ci_low", "ci_high"):
        assert b[col].to_numpy(float) == pytest.approx(
            scale * a[col].to_numpy(float), nan_ok=True, rel=1e-9, abs=1e-9
        )


@SETTINGS
@given(circle())
def test_single_arm_pooled_equals_that_arm(spec: tuple[int, float, list[float]]) -> None:
    """With every TU on the right arm the pooled fit is the right fit."""
    length, oric, _ = spec
    n = 6
    mids = od.right_arm_mids(oric, length, [0.1, 0.25, 0.4, 0.55, 0.7, 0.85])
    rng = np.random.default_rng(1)
    genes, tus, replicons = od.placed(
        [od.Placed("chr1", length, mids, values=list(20.0 + rng.standard_normal(n)), width=2)]
    )
    out = od.rows(genes, tus, replicons, od.oric_at("chr1", oric, width=2), od.config(**CFG))
    pooled, right = od.row_of(out, "pooled"), od.row_of(out, "right")
    assert float(pooled["value"]) == pytest.approx(float(right["value"]))
    assert is_na(od.row_of(out, "left")["value"])


# --------------------------------------------------------------------------- oriC


@SETTINGS
@given(circle())
def test_single_oric_call_is_its_midpoint(spec: tuple[int, float, list[float]]) -> None:
    length, oric, _ = spec
    feats = od.oric("chr1", (int(oric) - 10, int(oric) + 10))
    assert _oriter().oric_position(feats, "chr1", length) == pytest.approx(float(int(oric)))


@SETTINGS
@given(circle(), st.integers(0, 10**6), st.floats(0.0, 0.04))
def test_oric_rotates_with_its_calls(
    spec: tuple[int, float, list[float]], shift: int, gap: float
) -> None:
    """Rotating two close calls rotates their circular mean by the same offset."""
    length, oric, _ = spec
    d = gap * length
    base = od.oric_at("chr1", oric, (oric + d) % length, width=2)
    rot = od.oric_at("chr1", (oric + shift) % length, (oric + d + shift) % length, width=2)
    a = _oriter().oric_position(base, "chr1", length)
    b = _oriter().oric_position(rot, "chr1", length)
    assert a is not None and b is not None
    assert oo.circ_dist(float(b), float(a) + shift, float(length)) == pytest.approx(0.0, abs=1e-6)


@SETTINGS
@given(circle(), st.floats(0.051, 0.49))
def test_calls_spread_wider_than_five_percent_give_no_oric(
    spec: tuple[int, float, list[float]], gap: float
) -> None:
    length, oric, _ = spec
    feats = od.oric_at("chr1", oric, (oric + gap * length) % length, width=2)
    assert _oriter().oric_position(feats, "chr1", length) is None


@SETTINGS
@given(circle())
def test_other_feature_types_are_never_an_oric(spec: tuple[int, float, list[float]]) -> None:
    length, oric, _ = spec
    feats = od.features(
        [(t, "chr1", int(oric) - 10, int(oric) + 10) for t in ("oriv", "orit", "crispr", "rrna")]
    )
    assert _oriter().oric_position(feats, "chr1", length) is None
