"""The MVP-04 K-12 evidence run: every acceptance gate (docs/specs/MVP-04-ecoli-evidence-run.md,
"Acceptance criteria and expected numbers", Q9, Q10, lead decisions 1 and 3).

Reads ``$GIAP_MVP04_OUT/evidence_summary.json``, written by ``scripts/dev/mvp04_summary.py`` from
the run's published outputs, the giapsim inputs and the second giapsim call (gate 9). Skips when
``$GIAP_MVP04_OUT`` is unset; fails when it is set but the summary is missing. Runs weekly (Q9)
and once for the acceptance run, never in CI or in the session allocation.

The gates are constants **here**, checked against the raw values in the summary, independently of
the script's own ``status`` (which is checked too). Q10 / decision 1: after the first measured run
the lead re-sets the expected ranges from it and this file goes on ``ci/protected.txt``; a gate is
never widened (R3). A value outside its expected range but inside its gate is reported as a
pytest warning, not a failure.
"""

from __future__ import annotations

import json
import math
import os
import warnings
from pathlib import Path
from typing import Any

import pytest

pytestmark = [pytest.mark.e2e, pytest.mark.slow]

ENV = "GIAP_MVP04_OUT"
GATES = [str(n) for n in range(1, 17)]
FORBIDDEN = ("decoy_unscored", "decoy_pairing_low", "q_floor_above_alpha", "no_eligible_genes")

# ---- gates (spec table, column "Gate"); never widened (R3)
MIN_PAIRED_SHARE = 0.8  # 3
MAX_Q_FLOOR = 1 / 110  # 4: T >= 110
MIN_DECOY_WINS = 1  # 5
MIN_DETECTED, MAX_DETECTED_SHARE = 1000, 0.90  # 6
MAX_GENE_FDP = 0.02  # 7
TRUTH_SHARE_BAND = (0.55, 0.70)  # 8
MIN_HIGH, MIN_PRECISION, MIN_RECALL = 30, 0.95, 0.90  # 11
MIN_NOVEL_CLASSES, MIN_NOVEL_PRECISION = 2, 0.95  # 12
MAX_AUTOCORR_P, MIN_ISLANDS, MAX_SPURIOUS = 0.01, 4, 2  # 13
MIN_TUS = 1000  # 14

# ---- expected ranges (spec column "Expected"; re-set after the first run, decision 1)
EXPECTED: dict[str, tuple[float | None, float | None]] = {
    "3": (0.98, 1.0), "4": (2300, 3100), "5": (20, 300), "7": (0.0, 0.01), "10": (0.95, 1.0),
    "11": (60, 300), "13": (5, 6), "14": (1500, None), "16": (3, 8),
}  # fmt: skip


@pytest.fixture(scope="module")
def summary() -> dict[str, Any]:
    out = os.environ.get(ENV)
    if not out:
        pytest.skip(f"${ENV} is unset (the MVP-04 run outdir)")
    path = Path(out) / "evidence_summary.json"
    assert path.is_file(), f"{path} missing: run scripts/dev/mvp04_summary.py on the run first"
    data: dict[str, Any] = json.loads(path.read_text())
    return data


def g(summary: dict[str, Any], n: str) -> dict[str, Any]:
    entry: dict[str, Any] = summary["gates"][n]
    return entry


def _num(x: Any) -> float:
    assert x is not None, "value missing (NA)"
    v = float(x)
    assert math.isfinite(v), v
    return v


def test_every_gate_is_in_the_summary(summary: dict[str, Any]) -> None:
    assert sorted(summary["gates"], key=int) == GATES


@pytest.mark.parametrize("n", GATES)
def test_gate_status_is_not_fail(summary: dict[str, Any], n: str) -> None:
    e = g(summary, n)
    assert e["status"] in ("pass", "outside_expected"), (n, e)
    if e["status"] == "outside_expected":
        warnings.warn(f"MVP-04 gate {n} outside its expected range: {e['value']} vs "
                      f"{e['expected']}", stacklevel=1)  # fmt: skip
    lo, hi = EXPECTED.get(n, (None, None))
    if (lo, hi) != (None, None) and isinstance(e["value"], int | float):
        v = float(e["value"])
        if (lo is not None and v < lo) or (hi is not None and v > hi):
            warnings.warn(f"MVP-04 gate {n}: {v} outside [{lo}, {hi}]", stacklevel=1)


def test_gate1_run_completes(summary: dict[str, Any]) -> None:
    e = g(summary, "1")
    assert e["value"] == 0 and e["run_status"] == "completed" and not e["failed_steps"], e


def test_gate2_tier_and_warnings(summary: dict[str, Any]) -> None:
    e = g(summary, "2")
    assert e["value"] == "decoy_picked", e
    assert not set(e["warnings_found"]) & set(FORBIDDEN), e


def test_gate3_paired_share(summary: dict[str, Any]) -> None:
    assert _num(g(summary, "3")["value"]) >= MIN_PAIRED_SHARE


def test_gate4_q_floor(summary: dict[str, Any]) -> None:
    e = g(summary, "4")
    assert _num(e["q_floor"]) <= MAX_Q_FLOOR + 1e-12, e


def test_gate5_decoys_compete(summary: dict[str, Any]) -> None:
    assert _num(g(summary, "5")["value"]) >= MIN_DECOY_WINS


def test_gate6_detected_not_saturated(summary: dict[str, Any]) -> None:
    e = g(summary, "6")
    n = _num(e["value"])
    assert n >= MIN_DETECTED, e
    assert n <= MAX_DETECTED_SHARE * _num(e["n_eligible"]), e


def test_gate7_realised_gene_fdp(summary: dict[str, Any]) -> None:
    e = g(summary, "7")
    assert _num(e["value"]) <= MAX_GENE_FDP, e
    assert e["value"] == pytest.approx(e["n_false_detected"] / e["n_detected"])


def test_gate8_truth_detected_share(summary: dict[str, Any]) -> None:
    lo, hi = TRUTH_SHARE_BAND
    assert lo <= _num(g(summary, "8")["value"]) <= hi


def test_gate9_te7_inputs_byte_identical(summary: dict[str, Any]) -> None:
    assert g(summary, "9")["value"] == []


def test_gate10_tolerant_route_used(summary: dict[str, Any]) -> None:
    sh = g(summary, "10")["shares"]
    assert _num(sh["tolerant"]) > 0, sh
    assert sum(_num(sh[k]) for k in ("exact", "tolerant", "unmapped")) == pytest.approx(1.0)


def test_gate11_variants(summary: dict[str, Any]) -> None:
    e = g(summary, "11")
    assert _num(e["value"]) >= MIN_HIGH, e
    assert _num(e["precision"]) >= MIN_PRECISION, e
    assert _num(e["recall"]) >= MIN_RECALL, e
    assert e["n_rows_missing_flags"] == 0 and e["vcf_valid"] is True, e


def test_gate12_novel_orfs(summary: dict[str, Any]) -> None:
    e = g(summary, "12")
    live = {c: k for c, k in e["classes"].items() if k["n_calls"] >= 1}
    assert len(live) >= MIN_NOVEL_CLASSES, e
    for c, k in live.items():
        assert _num(k["precision"]) >= MIN_NOVEL_PRECISION, (c, k)


def test_gate13_context(summary: dict[str, Any]) -> None:
    e = g(summary, "13")
    ac = e["autocorrelation"]
    assert _num(ac["morans_i"]) > 0 and _num(ac["p_value"]) <= MAX_AUTOCORR_P, ac
    assert _num(e["value"]) >= MIN_ISLANDS and e["n_spurious"] <= MAX_SPURIOUS, e
    assert e["regions_ok"] is True, e


def test_gate14_operons(summary: dict[str, Any]) -> None:
    e = g(summary, "14")
    assert e["tests_run"] is True and e["too_few_op_pairs"] is False, e
    assert _num(e["n_tus"]) >= MIN_TUS, e


def test_gate15_report(summary: dict[str, Any]) -> None:
    e = g(summary, "15")
    assert e["empty_sections"] == [] and e["q_floor_caveat"] is False, e


def test_gate16_runtime_is_recorded(summary: dict[str, Any]) -> None:
    """Budget overruns are recorded (backlog), never a failure."""
    e = g(summary, "16")
    assert e["rows"], e
    if e["over_budget"]:
        warnings.warn(f"MVP-04 runtime over budget: {e['value']} h wall, {e['cpu_h']} CPU-h, "
                      f"{e['peak_rss_gb']} GB peak", stacklevel=1)  # fmt: skip
