"""M7 gate, variants half (M7-13a): G1 variant F1, G2 same-strain, attack 6 (A6a-c) and the
reference-proteome arm (docs/specs/M7-13-calibration-gate.md § Checks, rules 1-4, 7, 9-11, risks
1-2, lead decisions 1-3; GIAP_SPEC.md §7.5.6, §8.4, §8.5 attack 6, §13 row M7; ADR-0015, ADR-0042).

R is ``replicates_fast/full.variants`` (20 / 100, ADR-0042); every check is fixed-seed through
``helpers.m7_gate_records`` (seed bases there). Fast checks are nightly (> 60 s); the per-PR
subset is the one-replicate smoke test. Full tier: array task 13 (``variants``).

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
* Reference-proteome arm (Q8): every class with ``class_error_n >= 10`` has rate >= 0.95.
"""

from __future__ import annotations

import json
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
A6C_ADR = ("M7-13 risk 1: swap-ALT call-rule ADR (pending Khoa; options (i) a swap ALT is never"
           " high, (ii) high needs class_error_rate <= a bound, (iii) restate attack 6 as random"
           " swaps only)")  # fmt: skip


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
                     "non_mass_equivalent": dict(zip(("precision", "recall", "f1"), sub,
                                                     strict=True)),
                     "high_other": other, "masked": sum(r["masked"] for r in ran),
                     "seconds_per_replicate": np.mean([r["seconds"] for r in ran])})  # fmt: skip
    assert rec["status"] in NOT_FAILED, rec


@NIGHTLY
def test_g1_variant_f1() -> None:
    _g1(False)


@pytest.mark.calibration_full
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
    diagnostics({"check": "A6a", "full": full, "R_run": len(got), "n_rows": n_rows})
    bad = {r["replicate"]: r["a6a"] for r in got if r["a6a"]["n_unflagged"]}
    assert not bad, f"A6a: isobaric-only ALT rows without mass_equivalent: {bad}"
    assert n_rows >= 1, "A6a vacuous: no isobaric-only substitution/block row in the block"


@NIGHTLY
def test_a6a_isobaric_rows_flagged() -> None:
    _a6a(False)


@pytest.mark.calibration_full
def test_a6a_isobaric_rows_flagged_full() -> None:
    _a6a(True)


def _a6_fwer(full: bool, stream: str, record: Callable[[int], dict[str, Any]], key: str) -> None:
    R = C.replicates(SUITE, full)
    rec = C.calibrate(SUITE, "fwer", 0.05,
                      by_index(stream, R, lambda i: float(record(i)[key]["folded"] > 0)), R,
                      stream=stream)  # fmt: skip
    ran = [record(i) for i in ran_range(rec)]
    diagnostics({"check": key, "full": full, "status": rec["status"], "share": rec["observed"],
                 "unfolded_any": float(np.mean([r[key]["unfolded"] > 0 for r in ran]))
                 if ran else None})  # fmt: skip
    assert rec["status"] in NOT_FAILED, rec


@NIGHTLY
def test_a6b_random_swaps_fwer() -> None:
    _a6_fwer(False, "variants.a6b", M.a6_record, "a6b")


@pytest.mark.calibration_full
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
@pytest.mark.xfail(strict=True, raises=AssertionError, reason=A6C_ADR)
def test_a6c_systematic_swaps_fwer() -> None:
    _a6c(False)


@pytest.mark.calibration_full
@pytest.mark.xfail(strict=True, raises=AssertionError, reason=A6C_ADR)
def test_a6c_systematic_swaps_fwer_full() -> None:
    _a6c(True)


# --------------------------------------------------------------------------- reference arm


@NIGHTLY
def test_reference_proteome_arm_is_degenerate() -> None:
    got = [M.ref_arm_record(i) for i in reps(C.replicates(SUITE, False))]
    diagnostics({"check": "reference_proteome_arm", "R_run": len(got),
                 "n_dropped": [r["n_dropped"] for r in got],
                 "rates": [r["rates"] for r in got]})  # fmt: skip
    assert all(r["n_dropped"] > 0 for r in got), "transform removed nothing: no ALT pairing left"
    informative = [(r["replicate"], c, rn) for r in got for c, rn in r["rates"].items()
                   if rn[1] >= REF_MIN_N]  # fmt: skip
    assert informative, f"no class with class_error_n >= {REF_MIN_N}: {got}"
    low = [x for x in informative if x[2][0] < REF_MIN_RATE]
    assert not low, f"reference-proteome arm: class rate < {REF_MIN_RATE}: {low}"
