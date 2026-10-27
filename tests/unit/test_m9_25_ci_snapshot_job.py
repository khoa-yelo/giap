"""M9-25 work item 4 and index decision 22: one ubuntu-latest / Python 3.12 CI job runs
``pytest tests/snapshot --mpl``; every other job that runs the snapshot tests does so without
``--mpl`` (docs/specs/M9-25-snapshot-baselines.md Q1, default: a separate small ``snapshot``
job; GIPPA_SPEC.md §11.2, §12).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

CI = Path(__file__).resolve().parents[2] / ".github" / "workflows" / "ci.yml"


def _jobs() -> dict[str, Any]:
    doc = yaml.safe_load(CI.read_text(encoding="utf-8"))
    jobs: dict[str, Any] = doc["jobs"]
    return jobs


def _runs(job: dict[str, Any]) -> list[str]:
    return [str(s.get("run", "")) for s in job.get("steps", [])]


def _python(job: dict[str, Any]) -> str | None:
    env = {**(job.get("env") or {})}
    for step in job.get("steps", []):
        env |= step.get("env") or {}
        with_ = step.get("with") or {}
        if "python-version" in with_:
            return str(with_["python-version"])
    return None if "UV_PYTHON" not in env else str(env["UV_PYTHON"])


def _mpl_jobs() -> dict[str, Any]:
    return {n: j for n, j in _jobs().items() if any("--mpl" in r for r in _runs(j))}


def test_one_job_compares_snapshots() -> None:
    assert list(_mpl_jobs()) == ["snapshot"]


def test_snapshot_job_is_ubuntu_312_running_tests_snapshot_with_mpl() -> None:
    job = _jobs()["snapshot"]
    assert job["runs-on"] == "ubuntu-latest"
    assert "strategy" not in job  # one leg
    assert _python(job) == "3.12"
    cmd = next(r for r in _runs(job) if "--mpl" in r)
    assert "pytest" in cmd and "tests/snapshot" in cmd
    assert "if" not in job  # every push and pull request (§12)


def test_no_other_job_passes_mpl() -> None:
    for name, job in _jobs().items():
        if name != "snapshot":
            assert not any("--mpl" in r for r in _runs(job)), name
