"""M4-24: unscored and partly scored decoys in gene detection (docs/specs/M4-24-decoy-unscored.md).

Spec: GIAP_SPEC.md §7.3.2 picked steps 1-3 (targets and decoys compared on a scale both carry; a
tie counts as a decoy), the tier table ("Otherwise, q-values or PEPs exist" -> ``mayu``), §8
preamble / R4 (no silent claim), §5.3 (``q_value``, ``pep``, ``score`` optional on C9), ADR-0009
(1/f unchanged), ADR-0047 (decoys above ``q_value_max`` leave the picked competition; the filter
applies only when some C9 target row carries q). Lead decisions 2026-10-22: 2 (a partly NA decoy
column is dropped for everyone, with ``decoy_unscored``), 5 (one warning for both cases), 6 (a
sparse target column stays in S).

Pinned here:
- ``KEY_COLUMNS`` = ("q_value", "pep", "score"); ``shared_key_columns(targets, paired)`` returns the
  columns carried by the targets (>= 1 non-NA) and complete on the paired decoys (no NA), in
  ``KEY_COLUMNS`` order; an empty paired frame leaves every carried column (vacuous).
- Rule 1: the picked ranking uses only S; rule 2: ``decoy_unscored`` iff guard passes, target side
  non-empty and (S != C_T or S empty); rule 3: S empty -> the no-decoy chain, the tier rows still
  report the decoy counts, f and 1/f, ``q_floor`` NA; rule 4: the check order of ``choose_tier``;
  rule 5: one NA paired decoy drops the column; unpaired decoys never count; rule 6: the decoy q
  filter needs target q (symmetry).
- Examples 1-10 with the spec's numbers. Examples 1, 4 and 6 use target q = 1e-5 i instead of the
  table's 1e-4 i: at the spec's own ``q_value_max`` = 0.01, 1e-4 i rejects G101-G150 as targets
  (§7.3.1), while the numbers in the table assume all 150 accepted. 1e-5 i keeps every target
  accepted and the stated numbers exact (the decoy of example 4 moves with it: q = 1e-5).
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
import detection_tiers_fixtures as fx
from detection_tiers_fixtures import NAN, Builder, cfg, tag

UNSCORED = "decoy_unscored"
TIER_ROWS = ("n_decoys_paired", "n_decoys_unpaired", "decoy_paired_share", "decoy_win_weight")


def _det() -> Any:
    from giap.genetable import detection

    return detection


def _alpha(a: float, **kw: Any) -> Any:
    return cfg(stats={"alpha": {"identification": a}}, **kw)


def _frame(**cols: list[Any]) -> pd.DataFrame:
    """A key-column frame (``q_value``, ``pep``, ``score`` float, ``gene`` object)."""
    n = max(len(v) for v in cols.values())
    out = pd.DataFrame(index=range(n))
    for c in ("q_value", "pep", "score"):
        out[c] = pd.Series(cols.get(c, [NAN] * n), dtype="float64")
    out["gene"] = pd.Series(cols.get("gene", ["G"] * n), dtype=object)
    return out


def _qs(r: fx.Run, tags: list[str]) -> np.ndarray:
    return np.array([r.q(g) for g in tags])


def _pairing(r: fx.Run, *, q_value_max: float | None = None) -> pd.DataFrame:
    return _det().pair_decoys(r.peptides, r.det.protein_gene_map, cfg().input.decoy_prefixes,
                              q_value_max=q_value_max,
                              contaminant_prefixes=cfg().input.contaminant_prefixes)  # fmt: skip


def _choose(r: fx.Run, pairing: pd.DataFrame, accepted: pd.DataFrame | None) -> tuple[str, list]:
    tier, warns = _det().choose_tier(r.peptides, r.peptide_hits, r.protein_hits, pairing,
                                     accepted=accepted)  # fmt: skip
    return tier, list(warns)


def assert_same_as_decoy_free(with_decoys: fx.Run, without: fx.Run) -> None:
    """Rule 3: ``table`` and the ``q_threshold`` rows equal the run without decoy rows."""
    assert with_decoys.det.tier == without.det.tier
    pd.testing.assert_frame_equal(with_decoys.det.table, without.det.table)
    a = fx.curve_rows(with_decoys.det).reset_index(drop=True)
    b = fx.curve_rows(without.det).reset_index(drop=True)
    assert len(a), "the fallback tier writes q_threshold rows"
    pd.testing.assert_frame_equal(a, b)


# ---------------------------------------------------------------------------- public interface


def test_key_columns_constant() -> None:
    assert _det().KEY_COLUMNS == ("q_value", "pep", "score")


def test_shared_key_columns_follow_key_order() -> None:
    t = _frame(score=[1.0, 2.0], q_value=[0.001, 0.002], pep=[0.1, 0.2])
    d = _frame(pep=[0.3], score=[5.0], q_value=[0.004])
    assert _det().shared_key_columns(t, d) == ("q_value", "pep", "score")
    assert _det().shared_key_columns(t, _frame(score=[1.0], pep=[0.1])) == ("pep", "score")


def test_one_na_decoy_value_drops_the_column() -> None:
    t = _frame(q_value=[0.001, 0.002], score=[10.0, 20.0])
    d = _frame(q_value=[0.001, NAN, 0.003], score=[1.0, 2.0, 3.0])
    assert _det().shared_key_columns(t, d) == ("score",)


def test_one_non_na_target_value_carries_the_column() -> None:
    t = _frame(q_value=[NAN, NAN, 0.001], score=[10.0, 20.0, 30.0])
    d = _frame(q_value=[0.002, 0.003], score=[1.0, 2.0])
    assert _det().shared_key_columns(t, d) == ("q_value", "score")


def test_column_carried_only_by_decoys_is_not_shared() -> None:
    t = _frame(q_value=[0.001, 0.002])
    d = _frame(q_value=[0.002], pep=[0.1], score=[3.0])
    assert _det().shared_key_columns(t, d) == ("q_value",)


def test_empty_paired_frame_keeps_every_carried_column() -> None:
    t = _frame(q_value=[0.001], score=[10.0])
    assert _det().shared_key_columns(t, _frame(q_value=[]).iloc[:0]) == ("q_value", "score")


def test_targets_carrying_nothing_share_nothing() -> None:
    t = _frame(q_value=[NAN, NAN])
    d = _frame(q_value=[0.001], pep=[0.1], score=[3.0])
    assert _det().shared_key_columns(t, d) == ()


# ---------------------------------------------------------------------------- example 1 (all NA)


def _ex1(with_decoys: bool) -> Builder:
    b = Builder(150)
    t = b.tags or []
    for i, g in enumerate(t, start=1):
        b.target(g, 1e-5 * i, None, None)
    if with_decoys:
        for g in t:
            b.decoy([f"rev_{g}"], None, None, None)
    return b


def test_example1_all_unscored_decoys_fall_back_to_mayu() -> None:
    r = _ex1(True).run()
    assert r.det.tier == "mayu"
    assert UNSCORED in r.det.warnings
    assert set(r.det.table["detection_tier"].astype(str)) == {"mayu"}
    assert_same_as_decoy_free(r, _ex1(False).run())


def test_example1_tier_rows_report_the_decoys() -> None:
    r = _ex1(True).run()
    assert fx.tier_metric(r.det, "n_decoys_paired") == 150
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 0
    assert fx.tier_metric(r.det, "decoy_paired_share") == 1.0
    assert fx.tier_metric(r.det, "decoy_win_weight") == 1.0
    assert math.isnan(fx.tier_metric(r.det, "q_floor")), "q_floor is NA outside decoy_picked"


def test_example1_decoy_free_run_has_no_decoy_unscored() -> None:
    r = _ex1(False).run()
    assert r.det.tier == "mayu" and UNSCORED not in r.det.warnings


# ---------------------------------------------------------------------------- example 2 (score)


def test_example2_score_only_decoys_rank_on_score() -> None:
    b = Builder(4)
    t = b.tags or []
    for g, q, s in zip(t, (0.001, 0.002, 0.003, 0.004), (50.0, 40.0, 30.0, 20.0), strict=True):
        b.target(g, q, None, s)
    for g, s in ((t[0], 10.0), (t[2], 35.0), (t[3], 25.0)):
        b.decoy([f"rev_{g}"], None, None, s)
    r = b.run()
    assert r.det.tier == "decoy_picked"
    assert UNSCORED in r.det.warnings
    np.testing.assert_allclose(_qs(r, t[:2]), 0.5, rtol=1e-12)
    for g in t[2:]:
        assert r.status(g) == "not_detected" and math.isnan(r.q(g)), g
    assert not any(r.detected(g) for g in t)


# ---------------------------------------------------------------------------- example 3


def test_example3_score_only_decoys_beat_their_targets() -> None:
    b = Builder(300)
    t = b.tags or []
    for i, g in enumerate(t, start=1):
        b.target(g, 1e-5 * i, None, 301.0 - i)
    for i in range(201, 301):
        b.decoy([f"rev_{t[i - 1]}"], None, None, 301.0 - i + 0.5)
    r = b.run()
    assert r.det.tier == "decoy_picked"
    assert UNSCORED in r.det.warnings
    np.testing.assert_allclose(_qs(r, t[:200]), 1 / 200, rtol=1e-12)
    assert all(r.detected(g) for g in t[:200])
    for g in t[200:]:
        assert r.status(g) == "not_detected" and math.isnan(r.q(g)), g
    assert fx.tier_metric(r.det, "n_detected") == 200
    assert fx.tier_metric(r.det, "q_floor") == pytest.approx(1 / 200, rel=1e-12)


# ---------------------------------------------------------------------------- example 4 (PEP NA)


def test_example4_decoy_without_pep_ties_on_q_and_wins() -> None:
    b = Builder(150)
    t = b.tags or []
    for i, g in enumerate(t, start=1):
        b.target(g, 1e-5 * i, 1e-5 * i, None)
    b.decoy([f"rev_{t[0]}"], 1e-5, None, None)
    r = b.run()
    assert r.det.tier == "decoy_picked"
    assert UNSCORED in r.det.warnings
    assert r.status(t[0]) == "not_detected" and math.isnan(r.q(t[0]))
    np.testing.assert_allclose(_qs(r, t[1:]), 2 / 149, rtol=1e-12)
    assert fx.tier_metric(r.det, "n_detected") == 0


# ---------------------------------------------------------------------------- example 5 (q partial)


def test_example5_q_on_some_decoys_only_drops_q() -> None:
    b = Builder(400)
    t = b.tags or []
    for i, g in enumerate(t, start=1):
        b.target(g, 1e-5 * i, None, 1000.0 - i)
    b.decoy([f"rev_{t[100]}"], 5e-4, None, 949.5)
    b.decoy([f"rev_{t[101]}"], None, None, 949.4)
    r = b.run()
    assert r.det.tier == "decoy_picked"
    assert UNSCORED in r.det.warnings
    for g in (t[100], t[101]):
        assert r.status(g) == "not_detected" and math.isnan(r.q(g)), g
    rest = [g for g in t if g not in (t[100], t[101])]
    np.testing.assert_allclose(_qs(r, rest), 3 / 398, rtol=1e-12)
    assert all(r.detected(g) for g in rest)
    assert fx.tier_metric(r.det, "n_detected") == 398


# ---------------------------------------------------------------------------- example 6


def _ex6(with_decoys: bool) -> Builder:
    b = Builder(150)
    t = b.tags or []
    for i, g in enumerate(t, start=1):
        b.target(g, 1e-5 * i, None, None)
    if with_decoys:
        for g, q in zip(t[:4], (5e-4, 6e-4, 7e-4, None), strict=True):
            b.decoy([f"rev_{g}"], q, None, None)
    return b


def test_example6_one_na_q_decoy_empties_s_and_falls_back() -> None:
    r = _ex6(True).run()
    assert r.det.tier == "mayu"
    assert UNSCORED in r.det.warnings
    assert_same_as_decoy_free(r, _ex6(False).run())
    assert fx.tier_metric(r.det, "n_decoys_paired") == 4
    assert math.isnan(fx.tier_metric(r.det, "q_floor"))


# ---------------------------------------------------------------------------- example 7 (target NA)


def test_example7_sparse_target_column_is_unchanged_from_main() -> None:
    b = Builder(151)
    t = b.tags or []
    for i, g in enumerate(t[:150], start=1):
        b.target(g, 1e-5 * i, None, None)
    b.target(t[150], None, None, None)
    b.decoy([f"rev_{t[150]}"], 0.005, None, None)
    r = b.run()
    assert r.det.tier == "decoy_picked"
    assert UNSCORED not in r.det.warnings
    np.testing.assert_allclose(_qs(r, t[:150]), 1 / 150, rtol=1e-12)
    assert r.status(t[150]) == "not_detected" and math.isnan(r.q(t[150]))
    assert fx.tier_metric(r.det, "n_detected") == 150


# ---------------------------------------------------------------------------- example 8 (rule 4a)


def _ex8() -> Builder:
    b = Builder(12)
    t = b.tags or []
    for i, g in enumerate(t, start=1):
        b.target(g, 1e-4 * i, None, None)
    for g in t[:7]:
        b.decoy([f"rev_{g}"], None, None, None)
    for k in range(3):
        b.decoy([f"rev_UNKNOWN{k}"], None, None, None)
    return b


def test_example8_low_pairing_short_circuits_the_key_check() -> None:
    r = _ex8().run()
    assert r.det.tier == "mayu"
    assert "decoy_pairing_low" in r.det.warnings
    assert UNSCORED not in r.det.warnings
    pairing = _pairing(r, q_value_max=0.01)
    assert len(pairing) == 10 and int(pairing["gene"].notna().sum()) == 7
    assert _choose(r, pairing, r.accepted) == ("mayu", ["decoy_pairing_low"])


# ---------------------------------------------------------------------------- example 9 (rule 6)


def test_example9_no_target_q_means_no_decoy_q_filter() -> None:
    b = Builder(150)
    t = b.tags or []
    for i, g in enumerate(t, start=1):
        b.target(g, None, None, 151.0 - i)
    b.decoy([f"rev_{t[149]}"], 0.05, None, 0.5)
    r = b.run()
    assert fx.tier_metric(r.det, "n_decoys_paired") == 1
    assert r.det.tier == "decoy_picked"
    assert UNSCORED not in r.det.warnings
    np.testing.assert_allclose(_qs(r, t), 1 / 150, rtol=1e-12)
    assert all(r.detected(g) for g in t)


def test_decoy_q_filter_still_applies_when_targets_carry_q() -> None:
    """ADR-0047: with target q, a decoy above ``q_value_max`` is neither paired nor unpaired."""
    b = Builder(150)
    t = b.tags or []
    for i, g in enumerate(t, start=1):
        b.target(g, 1e-5 * i, None, 151.0 - i)
    b.decoy([f"rev_{t[149]}"], 0.05, None, 0.5)
    b.decoy([f"rev_{t[148]}"], 1e-3, None, 0.4)
    b.decoy([f"rev_{t[147]}"], None, None, 0.3)  # NA q passes, row by row (§7.3.1)
    r = b.run()
    assert fx.tier_metric(r.det, "n_decoys_paired") == 2
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 0
    assert r.det.tier == "decoy_picked"


def test_one_target_q_value_switches_the_decoy_filter_on() -> None:
    """Rule 6 reads *any* C9 non-decoy q, also on a target that is not accepted."""
    b = Builder(150)
    t = b.tags or []
    for i, g in enumerate(t, start=1):
        b.target(g, None, None, 151.0 - i)
    b.target(None, 0.5, None, 1.0, protein_ids=["OTHER"], placed=False)  # unplaced, has q
    b.decoy([f"rev_{t[149]}"], 0.05, None, 0.5)
    r = b.run()
    assert fx.tier_metric(r.det, "n_decoys_paired") == 0
    assert r.det.tier == "two_peptide"


def test_decoy_q_filter_symmetry_on_unpaired_counts() -> None:
    """Without target q the filtered-out decoy would have changed f: it now counts."""
    b = Builder(10)
    t = b.tags or []
    for i, g in enumerate(t, start=1):
        b.target(g, None, None, 100.0 - i)
    for g in t[:4]:
        b.decoy([f"rev_{g}"], 0.001, None, 1.0)
    b.decoy(["rev_UNKNOWN"], 0.5, None, 1.0)
    r = b.run()
    assert fx.tier_metric(r.det, "n_decoys_paired") == 4
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 1
    assert r.det.tier == "decoy_picked"  # share 0.8, guard passes


# ---------------------------------------------------------------------------- example 10 (DIA-NN)


def test_example10_diann_shaped_input_is_unchanged() -> None:
    """q and PEP on every row of both sides, no score: S = C_T = (q, PEP); M4-12 example 1 shape
    (q 0.5, 0.5, NA, 2/3, NA, 0.75) with a q/PEP tie going to the decoy."""
    b = Builder(6)
    t = b.tags or []
    qs = (0.001, 0.002, 0.003, 0.003, 0.004, 0.005)
    peps = (0.01, 0.02, 0.03, 0.035, 0.04, 0.05)
    for g, q, p in zip(t, qs, peps, strict=True):
        b.target(g, q, p, None)
    b.decoy([f"rev_{t[0]}"], 0.006, 0.06, None)
    b.decoy([f"rev_{t[2]}"], 0.0025, 0.02, None)
    b.decoy([f"rev_{t[4]}"], 0.004, 0.04, None)  # tie -> decoy
    r = b.run(_alpha(0.7))
    assert r.det.tier == "decoy_picked"
    assert UNSCORED not in r.det.warnings
    assert _det().shared_key_columns(r.accepted, _pairing(r, q_value_max=0.01)) == (
        "q_value",
        "pep",
    )
    np.testing.assert_allclose(_qs(r, t), [0.5, 0.5, NAN, 2 / 3, NAN, 0.75], rtol=1e-12)
    assert [r.detected(g) for g in t] == [True, True, False, True, False, False]


def test_sim50k_fixture_shares_every_target_key_column(tmp_path: Path, golden_seed: int) -> None:
    """Example 10: the sim50k M4 golden (M4-09 decoys carry q and PEP like the targets) keeps
    ``decoy_picked`` with S = C_T and no ``decoy_unscored``; its bytes are pinned by
    tests/golden/test_m4_golden.py (case ``sim50k``, unchanged)."""
    import m4_chain as mc

    from giap.genetable.covariates import covariate_table
    from giap.genetable.detection import detect
    from giap.genetable.evidence import accepted_peptides, gene_evidence
    from giap.io.work_tables import default_paths, load_m4_inputs

    work = mc.build(tmp_path / "work", seed=golden_seed, size="sim50k")
    m4, _ = load_m4_inputs(default_paths(work))
    c = mc.config(golden_seed, stats={"alpha": {"identification": 0.05}})
    acc = accepted_peptides(m4.peptide_hits, m4.peptides, c)
    ev = gene_evidence(m4.genes, m4.peptide_hits, acc, m4.replicons)
    cov = covariate_table(m4.genes, c)
    d = detect(m4.genes, ev, acc, m4.peptides, m4.peptide_hits, m4.protein_hits, cov, c)
    assert d.tier == "decoy_picked"
    assert UNSCORED not in d.warnings
    pairing = _det().pair_decoys(m4.peptides, d.protein_gene_map, c.input.decoy_prefixes,
                                 q_value_max=c.input.q_value_max,
                                 contaminant_prefixes=c.input.contaminant_prefixes)  # fmt: skip
    paired = pairing[pairing["gene"].notna()]
    assert len(paired)
    carried = tuple(k for k in _det().KEY_COLUMNS if acc[k].astype(float).notna().any())
    assert _det().shared_key_columns(acc, paired) == carried
    assert {"q_value", "pep"} <= set(carried)


# ---------------------------------------------------------------------------- rule 4 (check order)


def _scored_pairs(decoy: tuple[Any, Any, Any], target: tuple[Any, Any, Any]) -> fx.Run:
    b = Builder(4)
    t = b.tags or []
    for i, g in enumerate(t, start=1):
        b.target(g, 1e-3 * i if target[0] else None, 1e-2 * i if target[1] else None,
                 100.0 - i if target[2] else None)  # fmt: skip
    for g in t:
        b.decoy([f"rev_{g}"], *decoy)
    return b.run()


def test_rule4b_zero_target_rows_stays_picked_without_warning() -> None:
    r = _scored_pairs((None, None, None), (True, False, False))
    pairing = _pairing(r, q_value_max=0.01)
    empty = r.accepted.iloc[:0]
    assert _choose(r, pairing, empty) == ("decoy_picked", [])


@pytest.mark.parametrize(
    ("decoy", "target", "tier", "warn"),
    [
        ((0.5e-3, 0.005, 1.0), (True, True, True), "decoy_picked", False),  # S = C_T
        ((None, 0.005, 1.0), (True, True, True), "decoy_picked", True),  # S = (pep, score)
        ((None, None, 1.0), (True, True, True), "decoy_picked", True),  # S = (score,)
        ((None, None, None), (True, True, True), "mayu", True),  # S empty
        ((None, None, None), (False, True, False), "mayu", True),  # PEP-only targets
        ((None, 0.005, None), (False, True, False), "decoy_picked", False),
        ((0.5e-3, None, None), (False, False, True), "two_peptide", True),  # S empty, score only
    ],
)
def test_rule4c_choose_tier_on_the_shared_key(
    decoy: tuple[Any, Any, Any], target: tuple[bool, bool, bool], tier: str, warn: bool
) -> None:
    r = _scored_pairs(decoy, target)
    assert r.det.tier == tier
    assert (UNSCORED in r.det.warnings) == warn
    pairing = _pairing(r, q_value_max=0.01 if target[0] else None)
    got_tier, got_warns = _choose(r, pairing, r.accepted)
    assert got_tier == tier
    assert got_warns == ([UNSCORED] if warn else [])


def test_rule4c_without_accepted_reads_c9_targets() -> None:
    """``choose_tier`` without ``accepted``: the target side is the C9 non-decoy rows."""
    r = _scored_pairs((None, None, 1.0), (True, True, True))
    got_tier, got_warns = _choose(r, _pairing(r, q_value_max=0.01), None)
    assert (got_tier, got_warns) == ("decoy_picked", [UNSCORED])


def test_rule4d_no_pairing_rows_is_unchanged() -> None:
    b = Builder(3)
    for i, g in enumerate(b.tags or [], start=1):
        b.target(g, 1e-3 * i, None, None)
    b.decoy([], None, None, None)  # no protein_ids: no pairing row
    r = b.run()
    assert r.det.tier == "mayu" and UNSCORED not in r.det.warnings


# ---------------------------------------------------------------------------- degenerate cases


def test_zero_accepted_targets_with_paired_decoys() -> None:
    """Rule 4(b): every target fails q <= q_value_max; tier stays ``decoy_picked``."""
    b = Builder(3)
    t = b.tags or []
    for g in t:
        b.target(g, 0.5, None, None)
    for g in t:
        b.decoy([f"rev_{g}"], 0.001, None, None)
    r = b.run()
    assert len(r.accepted) == 0
    assert r.det.tier == "decoy_picked"
    assert UNSCORED not in r.det.warnings
    assert all(r.status(g) == "not_detected" and math.isnan(r.q(g)) for g in t)


def test_targets_without_key_column_give_two_peptide_and_warn() -> None:
    b = Builder(3)
    t = b.tags or []
    for g in t:
        b.target(g, None, None, None)
        b.target(g, None, None, None)
    for g in t:
        b.decoy([f"rev_{g}"], 0.001, 0.01, 5.0)
    r = b.run()
    assert r.det.tier == "two_peptide"
    assert UNSCORED in r.det.warnings
    assert all(r.status(g) == "detected" for g in t)
    assert fx.tier_metric(r.det, "n_decoys_paired") == 3


def test_targets_without_key_column_and_no_c11_give_protein_list() -> None:
    b = Builder(3)
    t = b.tags or []
    for g in t:
        b.target(g, None, None, None)
        b.decoy([f"rev_{g}"], 0.001, 0.01, 5.0)
    r = b.run()
    pairing = _pairing(r)
    no_hits = r.peptide_hits.iloc[:0]
    tier, warns = _det().choose_tier(r.peptides, no_hits, r.protein_hits, pairing,
                                     accepted=r.accepted)  # fmt: skip
    assert (tier, list(warns)) == ("protein_list", [UNSCORED])


def test_decoy_only_column_does_not_warn_and_q_tie_goes_to_decoy() -> None:
    b = Builder(3)
    t = b.tags or []
    for g, q in zip(t, (0.001, 0.002, 0.003), strict=True):
        b.target(g, q, None, None)
    b.decoy([f"rev_{t[0]}"], 0.001, None, 1.0)  # score only on the decoy side
    r = b.run(_alpha(1.0))
    assert r.det.tier == "decoy_picked"
    assert UNSCORED not in r.det.warnings
    assert r.status(t[0]) == "not_detected" and math.isnan(r.q(t[0]))


def test_decoy_only_pep_is_ignored_so_target_score_decides() -> None:
    """S = (q, score): a decoy-only PEP does not break the q tie; the higher target score wins."""
    b = Builder(1)
    g = (b.tags or [])[0]
    b.target(g, 0.001, None, 50.0)
    b.decoy([f"rev_{g}"], 0.001, 0.5, 1.0)
    r = b.run(_alpha(1.0))
    assert UNSCORED not in r.det.warnings
    assert r.detected(g)


def test_decoy_on_a_partial_gene_counts_for_completeness() -> None:
    b = Builder(4, partial={tag(3)})
    t = b.tags or []
    for i, g in enumerate(t[:3], start=1):
        b.target(g, 1e-3 * i, None, 100.0 - i)
    b.decoy([f"rev_{t[0]}"], 1e-4, None, 10.0)
    b.decoy([f"rev_{t[3]}"], None, None, 5.0)  # partial gene, NA q: drops q for everyone
    r = b.run(_alpha(1.0))
    assert r.det.tier == "decoy_picked"
    assert UNSCORED in r.det.warnings
    assert r.detected(t[0]), "ranked on score: the target (99) beats its decoy (10)"


def test_unpaired_decoys_never_enter_s() -> None:
    b = Builder(10)
    t = b.tags or []
    for i, g in enumerate(t, start=1):
        b.target(g, 1e-4 * i, None, None)
    for g in t[:9]:
        b.decoy([f"rev_{g}"], 0.002, None, None)
    b.decoy(["rev_UNKNOWN"], None, None, None)  # unpaired and unscored
    r = b.run()
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 1
    assert r.det.tier == "decoy_picked"
    assert UNSCORED not in r.det.warnings


def test_ties_on_integer_scores_go_to_the_decoy() -> None:
    b = Builder(2)
    t = b.tags or []
    b.target(t[0], 0.001, None, 7.0)
    b.target(t[1], 0.002, None, 9.0)
    b.decoy([f"rev_{t[0]}"], None, None, 7.0)
    r = b.run(_alpha(1.0))
    assert UNSCORED in r.det.warnings
    assert r.status(t[0]) == "not_detected" and math.isnan(r.q(t[0]))
    assert r.detected(t[1])


def test_fallback_tier_rows_report_partial_pairing() -> None:
    """Rule 3 with f < 1: the ``mayu`` tier rows still carry the decoy counts, f and 1/f."""
    b = Builder(12)
    t = b.tags or []
    for i, g in enumerate(t, start=1):
        b.target(g, 1e-4 * i, None, None)
    for g in t[:9]:
        b.decoy([f"rev_{g}"], None, None, None)
    b.decoy(["rev_UNKNOWN"], None, None, None)
    r = b.run()
    assert r.det.tier == "mayu" and UNSCORED in r.det.warnings
    vals = {m: fx.tier_metric(r.det, m) for m in TIER_ROWS}
    assert vals["n_decoys_paired"] == 9 and vals["n_decoys_unpaired"] == 1
    assert vals["decoy_paired_share"] == pytest.approx(0.9, rel=1e-12)
    assert vals["decoy_win_weight"] == pytest.approx(1 / 0.9, rel=1e-12)
    assert math.isnan(fx.tier_metric(r.det, "q_floor"))


def test_row_order_does_not_change_s_or_the_tier() -> None:
    def build() -> Builder:
        b = Builder(40)
        t = b.tags or []
        for i, g in enumerate(t, start=1):
            b.target(g, 1e-5 * i, None, 100.0 - i)
        for k, g in enumerate(t[20:]):
            b.decoy([f"rev_{g}"], None if k == 7 else 1e-4, None, 85.0 - k)
        return b

    a = build().run()
    for seed in range(3):
        s = build().run(shuffle=np.random.default_rng(seed))
        assert (s.det.tier, s.det.warnings) == (a.det.tier, a.det.warnings)
        pd.testing.assert_frame_equal(s.det.table, a.det.table)
    assert UNSCORED in a.det.warnings
