"""Shard and confirm runs of the variants gate give SKIPPED, never XPASS (docs/specs/
M7-13-calibration-gate.md lead decision 6). Under ``--shard`` or ``GIAP_CALIBRATION_CONFIRM`` the
harness returns ``shard`` / ``pending`` / ``skipped`` (no verdict); the strict-xfail gate tests
(G1, A6a, A6b, A6c) must write their record and then skip, or every task-13 shard exits 1.

Each case runs one gate test in a pytest subprocess with a stub plugin: the harness call and the
replicate records are replaced (no GIAP run), so only the test's own status handling is exercised.
Unmarked: runs on every commit.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import _criteria as C
import pytest

GATE = C.REPO / "tests" / "calibration" / "test_variants_gate.py"

STUB = textwrap.dedent(
    """
    import os
    import _criteria as C
    from helpers import m7_gate_records as M

    STATUS = os.environ["STUB_STATUS"]

    def _rec(*a, **k):
        return {"status": STATUS, "start": 0, "stop": 0, "observed": None, "bound": None,
                "precision": None, "recall": None, "se": None}

    C.calibrate_pooled = _rec
    C.calibrate = _rec
    clean = {"n_rows": 5, "n_unflagged": 0}
    M.a6_record = lambda i: {"replicate": i, "a6a": dict(clean),
                             "a6b": {"folded": 0, "unfolded": 0}}
    M.a6c_record = lambda i: {"replicate": i, "transform": {"n_sites": 3},
                              "a6c": {"folded": 0, "unfolded": 0}}
    """
)


def run_gate(tmp: Path, test: str, status: str, env: dict[str, str]) -> dict[str, int]:
    (tmp / "stub_gate.py").write_text(STUB)
    report = tmp / "report.json"
    keep = {k: v for k, v in os.environ.items() if k not in (C.SHARD_ENV, C.BLOCK_ENV,
                                                             C.CONFIRM_ENV)}  # fmt: skip
    path = os.pathsep.join([str(tmp), str(GATE.parent), keep.get("PYTHONPATH", "")])
    plugin = (
        "import json, sys\n"
        "out = {}\n"
        "def pytest_runtest_logreport(report):\n"
        "    if report.when == 'call' or report.outcome == 'skipped':\n"
        "        k = 'xpassed' if getattr(report, 'wasxfail', None) is not None and report.passed"
        " else 'xfailed' if getattr(report, 'wasxfail', None) is not None else report.outcome\n"
        "        out[k] = out.get(k, 0) + 1\n"
        "def pytest_sessionfinish(session, exitstatus):\n"
        f"    out['exit'] = int(exitstatus); open({str(report)!r}, 'w').write(json.dumps(out))\n"
    )
    (tmp / "count_outcomes.py").write_text(plugin)
    cmd = [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "-p", "stub_gate",
           "-p", "count_outcomes", "-o", "addopts=", f"{GATE}::{test}"]  # fmt: skip
    full = {**keep, "PYTHONPATH": path, "STUB_STATUS": status, C.LOG_DIR_ENV: str(tmp / "log"),
            **env}  # fmt: skip
    proc = subprocess.run(cmd, cwd=C.REPO, env=full, capture_output=True, text=True, timeout=300,
                          check=False)  # fmt: skip
    assert report.exists(), proc.stdout[-2000:] + proc.stderr[-2000:]
    out: dict[str, int] = json.loads(report.read_text())
    return out


SHARD = {C.SHARD_ENV: "0:10", C.BLOCK_ENV: "1"}
CONFIRM = {C.CONFIRM_ENV: json.dumps([["variants.other", "fwer", 0.05, 20]])}


@pytest.mark.parametrize(
    ("test", "status", "env"),
    [
        ("test_g1_variant_f1", C.SHARD_STATUS, SHARD),
        ("test_g1_variant_f1", C.SKIPPED_STATUS, CONFIRM),
        ("test_g1_variant_f1", C.PENDING_STATUS, {}),
        ("test_a6b_random_swaps_fwer", C.SHARD_STATUS, SHARD),
        ("test_a6b_random_swaps_fwer", C.SKIPPED_STATUS, CONFIRM),
        ("test_a6c_systematic_swaps_fwer", C.SHARD_STATUS, SHARD),
        ("test_a6a_isobaric_rows_flagged", C.SHARD_STATUS, SHARD),
    ],
)
def test_no_verdict_run_is_skipped_not_xpass(
    tmp_path: Path, test: str, status: str, env: dict[str, str]
) -> None:
    out = run_gate(tmp_path, test, status, env)
    assert out.get("skipped") == 1 and "xpassed" not in out, out
    assert out["exit"] == 0, out


def test_a_verdict_still_xpasses_strictly(tmp_path: Path) -> None:
    """Teeth: with a passing verdict (``pass``), the strict xfail still turns XPASS into exit 1."""
    out = run_gate(tmp_path, "test_a6b_random_swaps_fwer", "pass", {})
    assert out.get("skipped") is None and out["exit"] == 1, out
