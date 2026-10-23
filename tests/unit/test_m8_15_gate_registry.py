"""M8-15 acceptance bookkeeping (docs/specs/M8-15-calibration-gate.md § Acceptance criteria).

- "Three suites registered through ``calibrate``, fast and full variants, listed in the nightly
  array (``calibration_full.sbatch`` task list)": the ``#SBATCH --array`` range of
  ``scripts/slurm/calibration_full.sbatch`` covers the task index of ``ora``, ``gsea`` and
  ``determinants`` (task i = the i-th key of ``replicates_full``).
- "``benchmarks/baseline.json`` gains the fast-suite runtimes; the fast tier is < 10 min":
  proposed schema (test-author; the implementer may rename by lead decision):
  ``baseline["m8_calibration_fast"][<suite>]`` for each of the three suites, with ``test`` (a
  node id in one of the M8-15 gate files), ``secs`` (wall clock of the suite's fast tier,
  < 600), ``cpus`` and ``slurm_job`` (ints).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
SUITES = ("ora", "gsea", "determinants")
GATE_FILES = (
    "tests/calibration/test_function_calibration.py",
    "tests/calibration/test_determinants_calibration.py",
    "tests/calibration/attacks/test_attack_function.py",
)
FAST_BUDGET_S = 600.0


def _array_suites() -> list[str]:
    raw = yaml.safe_load((REPO / "tests" / "calibration" / "tolerances.yaml").read_text())
    return list(raw["replicates_full"])


def test_gate_files_exist() -> None:
    for f in GATE_FILES:
        assert (REPO / f).is_file(), f


@pytest.mark.parametrize("suite", SUITES)
def test_suite_is_in_the_nightly_array(suite: str) -> None:
    text = (REPO / "scripts" / "slurm" / "calibration_full.sbatch").read_text()
    m = re.search(r"^#SBATCH --array=(\d+)-(\d+)", text, re.M)
    assert m, "no #SBATCH --array line"
    lo, hi = int(m.group(1)), int(m.group(2))
    i = _array_suites().index(suite)
    assert lo <= i <= hi, f"{suite} is task {i}, outside --array={lo}-{hi}"


@pytest.mark.parametrize("suite", SUITES)
def test_baseline_records_the_fast_tier_runtime(suite: str) -> None:
    base = json.loads((REPO / "benchmarks" / "baseline.json").read_text())
    assert "m8_calibration_fast" in base, "benchmarks/baseline.json has no m8_calibration_fast"
    entry = base["m8_calibration_fast"][suite]
    node = str(entry["test"])
    assert node.split("::")[0] in GATE_FILES, node
    secs = float(entry["secs"])
    assert 0.0 < secs < FAST_BUDGET_S, f"{suite} fast tier {secs} s (budget {FAST_BUDGET_S} s)"
    assert isinstance(entry["cpus"], int) and entry["cpus"] >= 1
    assert isinstance(entry["slurm_job"], int)
