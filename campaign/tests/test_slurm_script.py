"""``campaign/slurm/campaign_array.sbatch`` (campaign/SPEC.md §10 Slurm; R14)."""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "campaign" / "slurm" / "campaign_array.sbatch"


@pytest.fixture(scope="module")
def text() -> str:
    assert SCRIPT.is_file(), f"missing {SCRIPT}"
    return SCRIPT.read_text()


def _directives(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in text.splitlines():
        m = re.match(r"#SBATCH\s+--([\w-]+)(?:=(\S+))?", line)
        if m:
            out[m.group(1)] = m.group(2) or ""
    return out


def _code(text: str) -> str:
    """The script without comment lines, continuation lines joined."""
    lines = [ln for ln in text.splitlines() if not ln.lstrip().startswith("#")]
    return re.sub(r"\\\n\s*", " ", "\n".join(lines))


def test_shebang(text: str) -> None:
    assert text.startswith("#!/bin/bash")


def test_resources(text: str) -> None:
    d = _directives(text)
    assert d["job-name"] == "giap-campaign"
    assert d["cpus-per-task"] == "16"
    assert d["mem"] == "64000M"
    assert d["time"] == "3-00:00:00"
    assert d["output"].endswith("/local/logs/%x.%A_%a.out")
    assert int(d["cpus-per-task"]) * 8 <= 128  # throttle <= 8 (R14)


def test_strict_mode_and_environment(text: str) -> None:
    code = _code(text)
    assert re.search(r"^set -euo pipefail$", code, re.M)
    assert re.search(r"^(source|\.) \S*scripts/farmshare_env\.sh\"?$", code, re.M)
    # set -euo pipefail before anything that can fail
    assert code.index("set -euo pipefail") < code.index("farmshare_env.sh")


def test_runs_the_array_task(text: str) -> None:
    code = _code(text)
    line = next(ln for ln in code.splitlines() if "run-job" in ln)
    assert "uv run --no-sync python -m giap_campaign run-job" in line
    assert '--jobs "$1"' in line
    assert '--index "$SLURM_ARRAY_TASK_ID"' in line
    assert '--db-dir "$2"' in line
    assert "--profile farmshare_local,apptainer" in line
    assert "--threads 16" in line


def test_giap_campaign_is_importable(text: str) -> None:
    # `python -m giap_campaign` needs campaign/src on the path (the package is not installed)
    assert "campaign/src" in _code(text)


def test_bash_syntax() -> None:
    bash = shutil.which("bash")
    assert bash is not None
    proc = subprocess.run(
        [bash, "-n", str(SCRIPT)], capture_output=True, text=True, timeout=30, check=False
    )
    assert proc.returncode == 0, proc.stderr
