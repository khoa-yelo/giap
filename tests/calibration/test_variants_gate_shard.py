"""Shard and confirm runs of the variants gate give SKIPPED, never XPASS (docs/specs/
M7-13-calibration-gate.md lead decision 6). Under ``--shard`` or ``GIAP_CALIBRATION_CONFIRM`` the
harness returns ``shard`` / ``pending`` / ``skipped`` (no verdict); the strict-xfail gate tests
(G1, A6a, A6b, A6c) must write their record and then skip, or every task-13 shard exits 1.

Lead decision 7 adds: a dirty A6a shard is a verdict (XFAIL, not a skip); under a confirm run
the tests that bypass the harness (G2, A6a, the A6a / A6c vacuity paths, M7-10's three full tests)
skip, while G1, A6b and A6c still reach the harness (their block 2 is logged there).

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
    A6A = os.environ.get("STUB_A6A", "clean")  # clean | unflagged (replicate 3) | vacuous
    A6C_SITES = int(os.environ.get("STUB_A6C_SITES", "3"))
    G2_BAD = int(os.environ.get("STUB_G2_VIOLATIONS", "0"))
    CALLS = os.environ.get("STUB_CALLS")

    def _rec(*a, **k):
        if CALLS:
            open(CALLS, "a").write("harness\\n")
        return {"status": STATUS, "start": 0, "stop": 0, "observed": None, "bound": None,
                "precision": None, "recall": None, "se": None}

    def _a6a(i):
        n_rows = 0 if A6A == "vacuous" else 5
        return {"n_rows": n_rows, "n_unflagged": int(A6A == "unflagged" and i == 3)}

    C.calibrate_pooled = _rec
    C.calibrate = _rec
    M.a6_record = lambda i: {"replicate": i, "a6a": _a6a(i),
                             "a6b": {"folded": 0, "unfolded": 0}}
    M.a6c_record = lambda i: {"replicate": i, "transform": {"n_sites": A6C_SITES},
                              "a6c": {"folded": 0, "unfolded": 0}}
    M.g2_record = lambda i: {"replicate": i, "violations": G2_BAD, "high": G2_BAD,
                             "start_codon_high": 0}
    """
)

#: M7-10's full tests bypass the harness; under a confirm run their body must not run. The stub
#: makes any body that gets past ``replicate_range`` fail at once (no GIAP run).
M710_STUB = textwrap.dedent(
    """
    import pytest

    def _boom(*a, **k):
        raise RuntimeError("M7-10 body ran under GIAP_CALIBRATION_CONFIRM")

    @pytest.hookimpl(tryfirst=True)
    def pytest_runtest_call(item):
        for name in ("_replicate", "_known_rate_replicate"):
            if hasattr(item.module, name):
                setattr(item.module, name, _boom)
    """
)


def run_gate(
    tmp: Path, test: str, status: str, env: dict[str, str], target: Path = GATE
) -> dict[str, int]:
    (tmp / "stub_gate.py").write_text(STUB)
    (tmp / "stub_m710.py").write_text(M710_STUB)
    report = tmp / "report.json"
    keep = {k: v for k, v in os.environ.items() if k not in (C.SHARD_ENV, C.BLOCK_ENV,
                                                             C.CONFIRM_ENV)}  # fmt: skip
    # tests/ for the -p plugins' `helpers` imports (test_branch.sbatch sets only src:tests/sim)
    tests_dir = GATE.parent.parent
    parts = [str(tmp), str(GATE.parent), str(tests_dir), keep.get("PYTHONPATH", "")]
    path = os.pathsep.join(parts)
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
           "-p", "stub_m710", "-p", "count_outcomes", "-o", "addopts=",
           f"{target}::{test}"]  # fmt: skip
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


# --------------------------------------------------------------------------- lead decision 7

#: a dirty shard is a verdict for A6a: one unflagged row anywhere fails the zero-tolerance check,
#: so the shard reports the ADR-0049 XFAIL, never a skip (code review B3, surviving mutant)


def test_a6a_shard_with_one_unflagged_row_is_xfailed_not_skipped(tmp_path: Path) -> None:
    out = run_gate(tmp_path, "test_a6a_isobaric_rows_flagged", C.SHARD_STATUS,
                   {**SHARD, "STUB_A6A": "unflagged"})  # fmt: skip
    assert out.get("xfailed") == 1 and "skipped" not in out, out
    assert out["exit"] == 0, out


#: confirm runs rerun block 2 of the listed harness checks only; a test that bypasses the harness
#: (G2, A6a, the A6a / A6c vacuity paths) has no block 2, so it skips there (lead decision 7)
CONFIRM_SKIPS = [
    ("test_g2_same_strain_no_internal_high", {"STUB_G2_VIOLATIONS": "1"}),
    ("test_g2_same_strain_no_internal_high_full", {"STUB_G2_VIOLATIONS": "1"}),
    ("test_g2_same_strain_no_internal_high_full", {}),
    ("test_a6a_isobaric_rows_flagged", {"STUB_A6A": "unflagged"}),
    ("test_a6a_isobaric_rows_flagged_full", {"STUB_A6A": "unflagged"}),
    ("test_a6a_isobaric_rows_flagged_full", {}),
    ("test_a6a_isobaric_rows_flagged_full", {"STUB_A6A": "vacuous"}),
    ("test_a6c_systematic_swaps_fwer_full", {"STUB_A6C_SITES": "0"}),
    ("test_a6c_systematic_swaps_fwer", {"STUB_A6C_SITES": "0"}),
]


@pytest.mark.parametrize(("test", "env"), CONFIRM_SKIPS)
def test_confirm_run_skips_tests_that_bypass_the_harness(
    tmp_path: Path, test: str, env: dict[str, str]
) -> None:
    out = run_gate(tmp_path, test, C.SKIPPED_STATUS, {**CONFIRM, **env})
    assert out.get("skipped") == 1 and not {"failed", "xfailed", "passed"} & set(out), out
    assert out["exit"] == 0, out


@pytest.mark.parametrize(
    "test", ["test_g1_variant_f1_full", "test_a6b_random_swaps_fwer_full",
             "test_a6c_systematic_swaps_fwer_full"],
)  # fmt: skip
def test_confirm_run_still_reaches_the_harness(tmp_path: Path, test: str) -> None:
    """Guard against over-skipping: a confirm run must still call the harness of G1, A6b and A6c
    (that is where block 2 is logged); only then may the test skip (``shard`` status)."""
    calls = tmp_path / "calls.txt"
    out = run_gate(tmp_path, test, C.SHARD_STATUS, {**CONFIRM, "STUB_CALLS": str(calls)})
    assert calls.exists() and calls.read_text().count("harness") == 1, out
    assert out.get("skipped") == 1 and out["exit"] == 0, out


M710 = [
    ("test_variant_class_error.py", "test_class_error_rate_unbiased_with_planted_variants_full"),
    ("test_variant_class_error.py", "test_class_error_rate_sanity_without_variants_full"),
    ("attacks/test_attack_class_error_rates.py", "test_known_rate_mix_full"),
]


@pytest.mark.parametrize(("path", "test"), M710)
def test_confirm_run_skips_m710_full_tests(tmp_path: Path, path: str, test: str) -> None:
    """M7-10's three full tests assert directly (no harness record): under a confirm run they
    must skip (a conftest-level rule; their files are M7-10's), never rerun block 1."""
    out = run_gate(tmp_path, test, C.SKIPPED_STATUS, CONFIRM, target=GATE.parent / path)
    assert out.get("skipped") == 1 and "failed" not in out, out
    assert out["exit"] == 0, out


def test_m710_stub_has_teeth(tmp_path: Path) -> None:
    """Without the confirm variable the stubbed M7-10 body runs and fails at once."""
    path, test = M710[2]
    out = run_gate(tmp_path, test, C.SKIPPED_STATUS, {}, target=GATE.parent / path)
    assert out.get("failed") == 1, out
