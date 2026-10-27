"""The CI ``package`` job (M11-01 item 5, acceptance criterion 3), read from ``ci.yml``."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
CI: dict[str, Any] = yaml.safe_load((REPO / ".github/workflows/ci.yml").read_text())
JOBS: dict[str, Any] = CI["jobs"]


def _job() -> dict[str, Any]:
    assert "package" in JOBS, "ci.yml has no `package` job (M11-01 item 5)"
    job: dict[str, Any] = JOBS["package"]
    return job


def _script(job: dict[str, Any]) -> str:
    return "\n".join(str(s.get("run", "")) for s in job.get("steps", []))


def _as_list(v: Any) -> list[str]:
    return [v] if isinstance(v, str) else list(v or [])


def test_package_job_runs_on_ubuntu_after_lint() -> None:
    job = _job()
    assert "ubuntu" in str(job["runs-on"])
    assert "lint" in _as_list(job.get("needs"))


def test_package_job_uses_python_312() -> None:
    job = _job()
    text = yaml.safe_dump(job)
    assert "3.12" in text


def test_summary_needs_package() -> None:
    assert "package" in _as_list(JOBS["summary"]["needs"])


def test_summary_table_reports_package() -> None:
    assert "needs.package.result" in _script(JOBS["summary"])


def test_build_and_check() -> None:
    job = _job()
    s = _script(job)
    assert "uv build" in s
    assert "scripts/release/check_dist.py dist" in s
    assert s.index("uv build") < s.index("check_dist.py")


def test_clean_venv_outside_checkout() -> None:
    job = _job()
    s = _script(job)
    assert "uv venv /tmp/" in s
    assert "uv pip install" in s and "dist/" in s and ".whl" in s


@pytest.mark.parametrize(
    "cmd",
    [
        "gippa --version",
        "gippa version --json",
        "gippa preflight",
        "workflow_main()",
        "nextflow config",
        "-profile test,stub",
        "-profile farmshare,apptainer",
    ],
)
def test_smoke_commands_present(cmd: str) -> None:
    job = _job()
    assert cmd in _script(job)


def test_preflight_uses_tiny_test_genome() -> None:
    job = _job()
    assert "tests/data/genomes/" in _script(job)


def test_workflow_main_asserted_under_site_packages() -> None:
    job = _job()
    s = _script(job)
    assert "site-packages" in s and "gippa/workflow" in s


def test_smoke_runs_outside_the_checkout() -> None:
    job = _job()
    """The smoke steps run from /tmp (``working-directory`` or ``cd /tmp``), not the checkout."""
    smoke = [
        st
        for st in job.get("steps", [])
        if "gippa --version" in str(st.get("run", "")) or "workflow_main()" in str(st.get("run", ""))
    ]
    assert smoke
    defaults = str(job.get("defaults", {}).get("run", {}).get("working-directory", ""))
    for st in smoke:
        wd = str(st.get("working-directory", defaults))
        assert wd.startswith("/tmp") or "cd /tmp" in str(st["run"]), st
