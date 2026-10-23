"""The nightly calibration array covers every ``replicates_full`` suite (docs/specs/
M7-13-calibration-gate.md § Array extension; ADR-0042; backlog line 614; GIAP_SPEC.md §8.4 "one
suite per task"): ``--array=0-(n-1)`` and the header comment "<n> suites" with n the number of
``replicates_full`` keys (15 today: task 13 ``variants``, 14 ``compare``)."""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests" / "calibration"))
import _criteria as C  # noqa: E402

SBATCH = REPO / "scripts" / "slurm" / "calibration_full.sbatch"


def test_array_spans_every_suite() -> None:
    text = SBATCH.read_text()
    suites = C.array_suites()
    assert len(suites) == 15 and suites[13:] == ("variants", "compare")
    arrays = re.findall(r"^#SBATCH --array=(\d+)-(\d+)\s*$", text, flags=re.M)
    assert arrays == [("0", str(len(suites) - 1))], arrays
    assert f"{len(suites)} suites" in text
