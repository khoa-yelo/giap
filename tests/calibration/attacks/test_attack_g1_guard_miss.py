"""Attack null of G1's guard (docs/specs/M7-13-calibration-gate.md rule 7, lead decision 7;
statistician review of #255, B3): variant calling broken, so no ``high`` substitution is predicted
in any replicate (``n_pred = 0``). Rule 7 calls this a broken scenario, not a gate verdict, so
``test_g1_variant_f1`` must *error* (``RuntimeError``), never report the ADR-0049 strict XFAIL
(``raises=AssertionError``) that the real F1 0.838 gives. Before the fix it reported 1 xfailed.

The gate test runs in a pytest subprocess with the statistician's stub plugin: ``g1_record`` is
replaced by empty-prediction counts (20 truth rows per replicate, so only the ``n_pred`` guard
misses); the real ``calibrate_pooled`` runs on them (no GIAP run). Unmarked: every commit.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import _criteria as C

GATE = C.REPO / "tests" / "calibration" / "test_variants_gate.py"

#: the statistician's stub (review comment 6365982711, B3), verbatim in substance
STUB = textwrap.dedent(
    """
    # attack null: variant calling broken (no high substitution at all) -> G1 guard miss
    from helpers import m7_gate_records as M

    def _rec(i, size):
        z = {"n_pred": 0, "tp_p": 0, "n_truth": 20, "tp_r": 0, "masked": 0}
        return {"replicate": i, "g1": dict(z), "g1_non_mass_equivalent": dict(z),
                "mutant_any_confidence": dict(z), "high_other": {}, "masked": 0,
                "seconds": 0.0}

    M.g1_record = _rec
    """
)

PLUGIN = """
import json
out = {"outcomes": [], "text": ""}
def pytest_runtest_logreport(report):
    if report.when == "call" or report.outcome != "passed":
        xf = getattr(report, "wasxfail", None) is not None
        out["outcomes"].append("xfailed" if xf and report.skipped else
                               "xpassed" if xf else report.outcome)
        out["text"] += str(report.longrepr or "")
def pytest_sessionfinish(session, exitstatus):
    out["exit"] = int(exitstatus)
    open(REPORT, "w").write(json.dumps(out))
"""


def test_g1_with_no_predictions_errors_not_xfails(tmp_path: Path) -> None:
    (tmp_path / "stub_g1.py").write_text(STUB)
    report = tmp_path / "report.json"
    (tmp_path / "count_g1.py").write_text(f"REPORT = {str(report)!r}\n" + PLUGIN)
    keep = {k: v for k, v in os.environ.items() if k not in (C.SHARD_ENV, C.BLOCK_ENV,
                                                             C.CONFIRM_ENV)}  # fmt: skip
    path = os.pathsep.join([str(tmp_path), str(GATE.parent), keep.get("PYTHONPATH", "")])
    cmd = [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "-p", "stub_g1",
           "-p", "count_g1", "-o", "addopts=", f"{GATE}::test_g1_variant_f1"]  # fmt: skip
    env = {**keep, "PYTHONPATH": path, C.LOG_DIR_ENV: str(tmp_path / "log")}
    proc = subprocess.run(cmd, cwd=C.REPO, env=env, capture_output=True, text=True, timeout=300,
                          check=False)  # fmt: skip
    assert report.exists(), proc.stdout[-2000:] + proc.stderr[-2000:]
    out = json.loads(report.read_text())
    assert "xfailed" not in out["outcomes"], out
    assert out["outcomes"] == ["failed"] and out["exit"] == 1, out
    assert "RuntimeError" in out["text"] and "guard" in out["text"], out["text"][-2000:]
