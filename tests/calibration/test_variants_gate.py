"""M7 gate, variants half (M7-13a): G1 variant F1, G2 same-strain, attack 6 (A6a-c) and the
reference-proteome arm (docs/specs/M7-13-calibration-gate.md § Checks, rules 1-4, 7, 9-11, risks
1-2, lead decisions 1-3; GIAP_SPEC.md §7.5.6, §8.4, §8.5 attack 6, §13 row M7; ADR-0015, ADR-0042).

R is ``replicates_fast/full.variants`` (20 / 100, ADR-0042); every check is fixed-seed through
``helpers.m7_gate_records`` (seed bases there). Fast checks are nightly (> 60 s); the per-PR
subset is the one-replicate smoke test. Full tier: array task 13 (``variants``).

Lead decision 4 (ADR-0049, proposed, pending Khoa): G1, A6a, A6b and A6c are strict xfails
(``raises=AssertionError``) until the method ADR is decided; each writes its diagnostics record
(estimate, CI, status) before it asserts, so the nightly shows the number. Vacuity raises
``RuntimeError`` (never an xfail). No threshold is widened (R3).

* G1: pooled F1 of ``high`` substitutions (rules 2-4) >= ``f1_min.variants_high`` through
  ``calibrate_pooled`` (point estimate, ADR-0015 blocks), guards sum n_truth >= 50 and sum n_pred
  >= 50 (rule 7). Fast at 100 kb, full at 1 Mb. Recorded: the non-``mass_equivalent`` subset's F1,
  the high block/swap rows, the masked count, the jackknife SE (rule 10).
* G2: zero ``high`` rows lacking both ``start_codon`` and ``stop_codon``, in every replicate.
* A6a: every isobaric-only substitution/block row carries ``mass_equivalent``, in every replicate
  (and at least one such row over the block: non-vacuity).
* A6b: FWER of "a folded high swap token" <= 0.05 + 2 SE (``calibrate``, criterion ``fwer``).
* A6c: as A6b after the systematic transform; expected to fail (risk 1): strict xfail naming the
  pending call-rule ADR (lead decision 2).
* Reference-proteome arm (Q8; a diagnostic, lead decision 4): per class, the mean rate over the
  replicates with ``class_error_n >= 10`` is >= 0.95 (pooled, not per replicate; CI recorded).
"""

from __future__ import annotations

import json
import os
from collections.abc import Callable
from typing import Any

import _criteria as C
import numpy as np
import pytest
import yaml
from helpers import m7_gate_records as M

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("variants")]
NIGHTLY = pytest.mark.calibration_nightly  # fast tests > 60 s (rule 11)

SUITE = "variants"
#: statuses that are not a failure; ``vacuous`` is deliberately absent (not a pass, rule 7)
NOT_FAILED = frozenset({"pass", "unconfirmed", C.SHARD_STATUS, C.SKIPPED_STATUS})
G1_GUARD = {"n_truth": 50, "n_pred": 50}
REF_MIN_RATE, REF_MIN_N = 0.95, 10
ADR_0049 = "ADR-0049 (proposed, pending Khoa)"
#: lead decision 6: a ``--shard`` or ``GIAP_CALIBRATION_CONFIRM`` run has no verdict (``skipped``:
#: a check not listed for block 2), so a strict xfail must not XPASS there
NO_VERDICT = frozenset({C.SHARD_STATUS, C.PENDING_STATUS, C.SKIPPED_STATUS})
SHARD_SKIP = "shard run: the verdict comes from --merge"


def skip_without_verdict(status: str) -> None:
    """Called after the record is written, before the assertion (lead decision 6)."""
    if status in NO_VERDICT:
        pytest.skip(SHARD_SKIP)


def xfail_0049(why: str) -> pytest.MarkDecorator:
    return pytest.mark.xfail(strict=True, raises=AssertionError, reason=f"{ADR_0049}: {why}")


G1_XFAIL = xfail_0049("G1 F1 below f1_min (risk 2)")
A6A_XFAIL = xfail_0049("isobaric blocks aligned as substitutions")
A6B_XFAIL = xfail_0049("coincident random swaps reach high")
A6C_XFAIL = xfail_0049("systematic swaps reach high (risk 1)")


def k_se() -> float:
    return C.load_tolerances().se_multiplier


def share_ci(p: float | None, n: int) -> list[float] | None:
    """``p +- k sqrt(p (1 - p) / n)`` clipped to [0, 1] (recorded only)."""
    if p is None or n < 1:
        return None
    h = k_se() * (p * (1 - p) / n) ** 0.5
    return [max(0.0, p - h), min(1.0, p + h)]


def f1_min() -> float:
    return float(yaml.safe_load(C.TOLERANCES_PATH.read_text())["f1_min"]["variants_high"])


def diagnostics(record: dict[str, Any]) -> None:
    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / "variants_gate_diagnostics.jsonl").open("a") as fh:
        fh.write(json.dumps(record, default=float) + "\n")


def by_index(stream: str, R: int, fn: Callable[[int], float]) -> Callable[[Any], float]:
    """``calibrate``'s ``fn(rng)`` for a function of the replicate index (rule 1: the index keys
    the shared simulation, the stream only names the check)."""
    index = {C.suite_seed(stream, i): i for i in range(2 * R)}
    return lambda rng: fn(index[int(rng.bit_generator.seed_seq.entropy)])


def reps(R: int) -> range:
    return C.replicate_range(R)


def ran_range(rec: dict[str, Any]) -> range:
    """Replicates a ``calibrate`` / ``calibrate_pooled`` record ran (none when skipped)."""
    return range(0) if rec["status"] == C.SKIPPED_STATUS else range(rec["start"], rec["stop"])


# --------------------------------------------------------------------------- per-PR smoke


def test_smoke_one_replicate_record_is_complete() -> None:
    rec = M.g1_record(0, "100k")
    for key in ("g1", "g1_non_mass_equivalent", "mutant_any_confidence"):
        assert set(rec[key]) == set(M.G1_COUNT_KEYS), (key, rec[key])
        assert all(isinstance(v, int) and v >= 0 for v in rec[key].values())
        assert rec[key]["tp_p"] <= rec[key]["n_pred"] and rec[key]["tp_r"] <= rec[key]["n_truth"]
    assert rec["g1"]["n_truth"] > 0 and rec["masked"] == rec["g1"]["masked"]


# --------------------------------------------------------------------------- G1


def _g1(full: bool) -> None:
    R, size = C.replicates(SUITE, full), "1M" if full else "100k"
    rec = C.calibrate_pooled(SUITE, M.G1_BASE, "g1", f1_min(), lambda i: M.g1_record(i, size), R,
                             guard=G1_GUARD)  # fmt: skip
    ran = [M.g1_record(i, size) for i in ran_range(rec)]
    if ran:
        sub = C.pooled_f1([r["g1_non_mass_equivalent"] for r in ran])
        other: dict[str, int] = {}
        for r in ran:
            for k, v in r["high_other"].items():
                other[k] = other.get(k, 0) + v
        diagnostics({"check": "G1", "full": full, "size": size, "status": rec["status"],
                     "f1": rec["observed"], "precision": rec["precision"],
                     "recall": rec["recall"], "se": rec["se"],
                     "ci": None if rec["observed"] is None or rec["se"] is None else
                     [rec["observed"] - k_se() * rec["se"], rec["observed"] + k_se() * rec["se"]],
                     "non_mass_equivalent": dict(zip(("precision", "recall", "f1"), sub,
                                                     strict=True)),
                     "high_other": other, "masked": sum(r["masked"] for r in ran),
                     "seconds_per_replicate": np.mean([r["seconds"] for r in ran])})  # fmt: skip
    skip_without_verdict(rec["status"])
    assert rec["status"] in NOT_FAILED, rec


@NIGHTLY
@G1_XFAIL
def test_g1_variant_f1() -> None:
    _g1(False)


@pytest.mark.calibration_full
@G1_XFAIL
def test_g1_variant_f1_full() -> None:
    _g1(True)


# --------------------------------------------------------------------------- G2


def _g2(full: bool) -> None:
    got = [M.g2_record(i) for i in reps(C.replicates(SUITE, full))]
    diagnostics({"check": "G2", "full": full, "R_run": len(got),
                 "high": sum(r["high"] for r in got),
                 "start_codon_high": sum(r["start_codon_high"] for r in got)})  # fmt: skip
    bad = {r["replicate"]: r["violations"] for r in got if r["violations"]}
    assert not bad, f"G2: high ALT rows without start/stop flag in replicates {bad}"


@NIGHTLY
def test_g2_same_strain_no_internal_high() -> None:
    _g2(False)


@pytest.mark.calibration_full
def test_g2_same_strain_no_internal_high_full() -> None:
    _g2(True)


# --------------------------------------------------------------------------- A6a / A6b


def _a6a(full: bool) -> None:
    got = [M.a6_record(i) for i in reps(C.replicates(SUITE, full))]
    n_rows = sum(r["a6a"]["n_rows"] for r in got)
    n_unflagged = sum(r["a6a"]["n_unflagged"] for r in got)
    bad = {r["replicate"]: r["a6a"] for r in got if r["a6a"]["n_unflagged"]}
    share = len(bad) / len(got) if got else None
    diagnostics({"check": "A6a", "full": full, "R_run": len(got), "n_rows": n_rows,
                 "n_unflagged": n_unflagged, "share_replicates_unflagged": share,
                 "ci": share_ci(share, len(got)), "status": "fail" if bad else "pass"})  # fmt: skip
    if got and n_rows < 1:  # not an AssertionError: never an xfail
        raise RuntimeError("A6a vacuous: no isobaric-only substitution/block row in the block")
    if not bad and os.environ.get(C.SHARD_ENV):  # a clean partial range is no verdict
        pytest.skip(SHARD_SKIP)
    assert not bad, f"A6a: isobaric-only ALT rows without mass_equivalent: {bad}"


@NIGHTLY
@A6A_XFAIL
def test_a6a_isobaric_rows_flagged() -> None:
    _a6a(False)


@pytest.mark.calibration_full
@A6A_XFAIL
def test_a6a_isobaric_rows_flagged_full() -> None:
    _a6a(True)


def _a6_fwer(full: bool, stream: str, record: Callable[[int], dict[str, Any]], key: str) -> None:
    R = C.replicates(SUITE, full)
    rec = C.calibrate(SUITE, "fwer", 0.05,
                      by_index(stream, R, lambda i: float(record(i)[key]["folded"] > 0)), R,
                      stream=stream)  # fmt: skip
    ran = [record(i) for i in ran_range(rec)]
    diagnostics({"check": key, "full": full, "status": rec["status"], "share": rec["observed"],
                 "bound": rec["bound"], "ci": share_ci(rec["observed"], len(ran)),
                 "unfolded_any": float(np.mean([r[key]["unfolded"] > 0 for r in ran]))
                 if ran else None})  # fmt: skip
    skip_without_verdict(rec["status"])
    assert rec["status"] in NOT_FAILED, rec


@NIGHTLY
@A6B_XFAIL
def test_a6b_random_swaps_fwer() -> None:
    _a6_fwer(False, "variants.a6b", M.a6_record, "a6b")


@pytest.mark.calibration_full
@A6B_XFAIL
def test_a6b_random_swaps_fwer_full() -> None:
    _a6_fwer(True, "variants.a6b", M.a6_record, "a6b")


# --------------------------------------------------------------------------- A6c (risk 1)


def _a6c(full: bool) -> None:
    R = C.replicates(SUITE, full)
    first = M.a6c_record(next(iter(reps(R))))
    if first["transform"].get("n_sites", 0) < 1:  # not an AssertionError: never an xfail
        raise RuntimeError(f"A6c vacuous: the systematic transform drew no site {first}")
    _a6_fwer(full, M.A6C_BASE, M.a6c_record, "a6c")


@NIGHTLY
@A6C_XFAIL
def test_a6c_systematic_swaps_fwer() -> None:
    _a6c(False)


@pytest.mark.calibration_full
@A6C_XFAIL
def test_a6c_systematic_swaps_fwer_full() -> None:
    _a6c(True)


# --------------------------------------------------------------------------- reference arm


@NIGHTLY
def test_reference_proteome_arm_is_degenerate() -> None:
    got = [M.ref_arm_record(i) for i in reps(C.replicates(SUITE, False))]
    per_class: dict[str, list[float]] = {}
    for r in got:
        for c, (rate, n) in r["rates"].items():
            if n >= REF_MIN_N:
                per_class.setdefault(c, []).append(float(rate))
    pooled = {}
    for c, v in sorted(per_class.items()):
        m = float(np.mean(v))
        h = k_se() * float(np.std(v, ddof=1)) / len(v) ** 0.5 if len(v) > 1 else None
        pooled[c] = {"mean_rate": m, "R": len(v), "ci": None if h is None else [m - h, m + h]}
    diagnostics({"check": "reference_proteome_arm", "R_run": len(got),
                 "n_dropped": [r["n_dropped"] for r in got], "pooled": pooled,
                 "rates": [r["rates"] for r in got]})  # fmt: skip
    assert all(r["n_dropped"] > 0 for r in got), "transform removed nothing: no ALT pairing left"
    assert pooled, f"no class with class_error_n >= {REF_MIN_N}: {got}"
    low = {c: p for c, p in pooled.items() if p["mean_rate"] < REF_MIN_RATE}
    assert not low, f"reference-proteome arm: pooled mean rate < {REF_MIN_RATE}: {low}"
