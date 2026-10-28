"""PR #179 code-review blockers on the TU-counted ORA (docs/specs/M8-04 lead decision 13).

Pinned here (the statistician's p-value domain, convergence and attack-arm files cover the rest):

- **Skip path.** When the TU recalibration is skipped (``folds=None``, a fold whose training TUs
  have fewer than 10 of a class, or a fit without an MLE), every row carries
  ``tu_recal_skipped`` and ``ora_uncalibrated``. ``folds`` is a required keyword of
  ``tu_universe`` and ``detection_enrichment``: omitting it is a ``TypeError``; ``folds=None``
  is the caller's explicit choice.
- **``ora_uncalibrated`` is on** (ADR-0035 decision 3): ``CALIBRATION_FLAGS ==
  ("ora_uncalibrated",)`` and every row of every run carries it.
- **Degenerate run** (degenerate-case table; review B2). An undetected TU of 700 genes at
  p̂ 0.9 has log comp = 700 log 0.1 ≈ -1611.8, so comp underflows to 0 and q_raw rounds to 1:
  the observation has probability 0. Every row has ``p_value``, ``p_noelim``, ``q_value`` and
  ``effect`` NA, flags ``unestimated`` and ``null_degenerate``; ``null_degenerate`` is a run flag
  and a run-level ``UserWarning`` is emitted. The same TU *detected* is not degenerate (the
  ``& ~detected`` condition): the run is tested normally and does not warn. A detected TU with
  q == 0 is not reachable from a gene frame (p̂ is clipped to >= 1e-12, so q_raw >= 1e-12; the
  fitted q is a logistic value kept open by ``tilt_pair``), so it stays pinned on the primitive
  (M8-04a).
- **p_detect outside [0, 1]** is a ``ValueError`` naming the locus tags (fail closed, M4-16).
  Exact 0 and 1 are still clipped and flagged ``p_detect_clipped``.
"""

from __future__ import annotations

import math
import warnings
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.gsea_data import as_list
from helpers.ora_data import T6, Gene, frame, kegg, ora, singletons, t6_frame

from gippa.genetable.detectability import tu_folds

ALPHA = 0.05
ALL = (1, 500)
SETS = {"x": {"g1", "g4", "g7"}, "y": {"g3", "g6"}}


def _run(genes: pd.DataFrame, sets: Any = None, **kw: Any) -> pd.DataFrame:
    return ora().detection_enrichment(genes, [kegg(sets or SETS)], ALL, ALPHA, **kw)


def _flags(v: Any) -> set[str]:
    return set(as_list(v))


def _all_rows_have(out: pd.DataFrame, *flags: str) -> None:
    assert len(out) > 0
    for i, f in enumerate(out["flags"]):
        missing = set(flags) - _flags(f)
        assert not missing, (i, out.iloc[i]["term"], out.iloc[i]["direction"], sorted(_flags(f)))


def _big_tu_frame(detected: bool) -> pd.DataFrame:
    """T6 plus TU_7: 700 genes h1..h700 at p̂ 0.9, all detected or all undetected."""
    tus = {
        tu: [Gene(g, 10.0 if g in {"g2", "g4", "g6", "g9"} else None,
                  g in {"g2", "g4", "g6", "g9"}, p_detect=p) for g, p in gs]
        for tu, gs in T6.items()
    }  # fmt: skip
    tus["TU_7"] = [
        Gene(f"h{i}", 10.0 if detected else None, detected, p_detect=0.9) for i in range(1, 701)
    ]
    return frame(tus)


BIG_SETS = {"x": {"g1", "g4", "g7"}, "y": {"g3", "g6", "h1"}, "z": {"g10", "h5"}}


def _few_positive_singletons() -> tuple[pd.DataFrame, dict[str, int]]:
    """60 single-gene TUs, 8 detected: every fold's training TUs have < 10 detected."""
    rng = np.random.default_rng(17901)
    p = rng.uniform(0.05, 0.4, 60)
    det = sorted(rng.choice(60, 8, replace=False).tolist())
    g = singletons(p, det)
    return g, tu_folds([f"TU_{i}" for i in range(1, 61)], 5, rng)


def _separated_singletons() -> tuple[pd.DataFrame, dict[str, int]]:
    """200 single-gene TUs detected exactly when p̂ > 0.5 (complete separation on logit q_raw):
    each class has ~80 training TUs per fold, but the logistic MLE does not exist (the
    likelihood increases without bound along the slope), so the fit cannot converge."""
    rng = np.random.default_rng(17902)
    p = rng.uniform(0.05, 0.95, 200)
    g = singletons(p, np.flatnonzero(p > 0.5).tolist())
    return g, tu_folds([f"TU_{i}" for i in range(1, 201)], 5, rng)


# ------------------------------------------------------------------- 1. skip path row flag


def test_folds_none_flags_every_row_skipped() -> None:
    out = _run(t6_frame(), folds=None)
    assert "tu_recal_skipped" in out.attrs["run_flags"]
    _all_rows_have(out, "tu_recal_skipped", "ora_uncalibrated")


def test_a_fold_short_of_a_class_flags_every_row_skipped() -> None:
    g, folds = _few_positive_singletons()
    sets = {"a": {f"g{i}" for i in range(1, 21)}, "b": {f"g{i}" for i in range(30, 45)}}
    u = ora().tu_universe(g, folds=folds)
    assert "tu_recal_skipped" in u.flags
    out = _run(g, sets, folds=folds)
    assert "tu_recal_skipped" in out.attrs["run_flags"]
    _all_rows_have(out, "tu_recal_skipped", "ora_uncalibrated")


def test_a_fit_without_an_mle_is_skipped_and_flags_every_row() -> None:
    g, folds = _separated_singletons()
    u = ora().tu_universe(g, folds=folds)
    assert "tu_recal_skipped" in u.flags
    np.testing.assert_array_equal(u.q, u.q_raw)
    sets = {"a": {f"g{i}" for i in range(1, 41)}, "b": {f"g{i}" for i in range(100, 130)}}
    out = _run(g, sets, folds=folds)
    assert "tu_recal_skipped" in out.attrs["run_flags"]
    _all_rows_have(out, "tu_recal_skipped", "ora_uncalibrated")


def test_a_fitted_run_does_not_flag_rows_skipped() -> None:
    from helpers.ora_calsim import collection, simulate, whole_tu_sets

    rng = np.random.default_rng(17903)
    sim = simulate(rng, c=1.0)
    out = ora().detection_enrichment(
        sim.genes, [collection(whole_tu_sets(sim, rng, "random", 6))], (5, 500), ALPHA,
        folds=sim.folds,
    )  # fmt: skip
    assert "tu_recal_skipped" not in out.attrs["run_flags"]
    assert len(out) > 0
    assert not any("tu_recal_skipped" in _flags(f) for f in out["flags"])
    _all_rows_have(out, "ora_uncalibrated")


def test_tu_universe_folds_is_a_required_keyword() -> None:
    with pytest.raises(TypeError):
        ora().tu_universe(t6_frame())  # type: ignore[call-arg]
    u = ora().tu_universe(t6_frame(), folds=None)
    assert "tu_recal_skipped" in u.flags


def test_detection_enrichment_folds_is_a_required_keyword() -> None:
    with pytest.raises(TypeError):
        ora().detection_enrichment(t6_frame(), [kegg(SETS)], ALL, ALPHA)  # type: ignore[call-arg]
    out = ora().detection_enrichment(t6_frame(), [kegg(SETS)], ALL, ALPHA, folds=None)
    assert len(out) == 4


# ------------------------------------------------------------------ 2. ora_uncalibrated on


def test_calibration_flags_is_ora_uncalibrated() -> None:
    assert ora().CALIBRATION_FLAGS == ("ora_uncalibrated",)


def test_every_row_of_every_run_carries_ora_uncalibrated() -> None:
    from helpers.ora_calsim import collection, simulate, whole_tu_sets

    rng = np.random.default_rng(17904)
    sim = simulate(rng, c=2.5)
    clipped = t6_frame()
    clipped.loc[clipped["locus_tag"] == "g10", "p_detect"] = 1.0
    runs = {
        "t6_raw": _run(t6_frame(), folds=None),
        "no_detection": _run(t6_frame(detected=[]), folds=None),
        "clipped": _run(clipped, folds=None),
        "fitted": ora().detection_enrichment(
            sim.genes,
            [collection(whole_tu_sets(sim, rng, "large", 5))],
            (5, 500),
            ALPHA,
            folds=sim.folds,
        ),
    }
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        runs["degenerate"] = _run(_big_tu_frame(detected=False), BIG_SETS, folds=None)
    for name, out in runs.items():
        assert len(out) > 0, name
        for f in out["flags"]:
            assert "ora_uncalibrated" in _flags(f), (name, sorted(_flags(f)))


# ---------------------------------------------------------------------- 3. degenerate run


def test_big_undetected_tu_underflows_its_complement() -> None:
    u = ora().tu_universe(_big_tu_frame(detected=False), folds=None)
    i = u.gene_tu["h1"]
    assert u.log_comp[i] == pytest.approx(700 * math.log(0.1), rel=1e-12)
    assert u.log_comp[i] < -1600
    assert u.comp[i] == 0.0 and u.q[i] == 1.0 and not u.detected[i]


def test_degenerate_run_is_all_na_flagged_and_warns() -> None:
    with pytest.warns(UserWarning):
        out = _run(_big_tu_frame(detected=False), BIG_SETS, folds=None)
    assert len(out) == 6
    for c in ("p_value", "p_noelim", "q_value", "effect"):
        assert out[c].isna().all(), c
    _all_rows_have(out, "unestimated", "null_degenerate")
    assert "null_degenerate" in out.attrs["run_flags"]


def test_degenerate_run_warns_once_per_run_not_per_row() -> None:
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        _run(_big_tu_frame(detected=False), BIG_SETS, folds=None)
    user = [w for w in rec if issubclass(w.category, UserWarning)]
    assert len(user) == 1, [str(w.message) for w in user]


def test_a_detected_certain_tu_is_not_degenerate() -> None:
    """comp == 0 on a *detected* TU is consistent with the data: the ``& ~detected`` condition."""
    g = _big_tu_frame(detected=True)
    u = ora().tu_universe(g, folds=None)
    i = u.gene_tu["h1"]
    assert u.comp[i] == 0.0 and u.detected[i]
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        out = _run(g, BIG_SETS, folds=None)
    assert not [w for w in rec if "degenerate" in str(w.message).lower()]
    assert "null_degenerate" not in out.attrs["run_flags"]
    assert np.isfinite(out["p_value"]).all() and np.isfinite(out["p_noelim"]).all()
    for f in out["flags"]:
        assert not _flags(f) & {"null_degenerate", "unestimated"}, sorted(_flags(f))


def test_a_regular_run_does_not_warn() -> None:
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        out = _run(t6_frame(), folds=None)
    assert not [w for w in rec if issubclass(w.category, UserWarning)]
    assert "null_degenerate" not in out.attrs["run_flags"]


# -------------------------------------------------------------- 5. p_detect outside [0, 1]


@pytest.mark.parametrize("bad", [1.5, -0.2, 1.0 + 1e-9, -1e-300])
def test_p_detect_outside_the_unit_interval_fails_closed(bad: float) -> None:
    g = t6_frame()
    g.loc[g["locus_tag"] == "g7", "p_detect"] = bad
    with pytest.raises(ValueError, match=r"\bg7\b"):
        ora().tu_universe(g, folds=None)
    with pytest.raises(ValueError, match=r"\bg7\b"):
        _run(g, folds=None)


def test_every_out_of_range_locus_tag_is_named() -> None:
    g = t6_frame()
    g.loc[g["locus_tag"] == "g3", "p_detect"] = 1.5
    g.loc[g["locus_tag"] == "g9", "p_detect"] = -0.2
    with pytest.raises(ValueError) as e:
        _run(g, folds=None)
    msg = str(e.value)
    assert "g3" in msg and "g9" in msg
    assert "g1," not in msg and "g10" not in msg  # only the offending genes


def test_an_out_of_range_partial_gene_is_ignored() -> None:
    """Fail closed applies to eligible genes only, as for NaN p_detect (spec §2)."""
    g = t6_frame()
    g.loc[g["locus_tag"] == "g7", ["p_detect", "partial"]] = [1.5, True]
    out = _run(g, folds=None)
    assert len(out) == 4


@pytest.mark.parametrize(("gene", "p", "detected"), [("g10", 1.0, False), ("g3", 0.0, True)])
def test_exact_0_and_1_are_still_clipped_and_flagged(gene: str, p: float, detected: bool) -> None:
    g = t6_frame()
    g.loc[g["locus_tag"] == gene, ["p_detect", "detected"]] = [p, detected]
    out = _run(g, {"has": {gene, "g1", "g4"}, "not": {"g6", "g7"}}, folds=None)
    assert "p_detect_clipped" in out.attrs["run_flags"]
    has = out[out["term"] == "has"]
    assert all("p_detect_clipped" in _flags(f) for f in has["flags"])
    assert np.isfinite(out["p_value"]).all()
