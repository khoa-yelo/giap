"""The docker-smoke DB-visibility job, parsed statically (M11-14 item 3b, acceptance 11; 14p).

Spec: "``docker-smoke.yml`` gains one job that runs ``BAKTA`` on the 100-kb test genome with a
Bakta light DB fixture, or checks that ``params.db_dir`` is visible and readable inside the BAKTA
container (``test -r <db_dir>/version.json``) ... It is advisory, like the existing jobs."

Readings used here:

- the new job is the one job of ``.github/workflows/docker-smoke.yml`` other than the M11-09
  ``docker-smoke`` job; it names BAKTA;
- it goes through Nextflow's docker profile (``nextflow run``/``nf-test`` with a profile list that
  contains ``docker``) and sets ``db_dir``, so the bind under test is the profile's, not a
  hand-written ``docker run -v``;
- the check reads ``version.json`` (``test -r``) or runs BAKTA on the test genome;
- the BAKTA image comes from ``containers.config`` (by digest), never a literal tag;
- advisory: no job needs it, ``ci.yml``'s ``summary`` does not need it; least privilege and a
  timeout as for the M11-09 job.

Nothing is dispatched; no network, no Docker.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
WORKFLOWS = REPO / ".github" / "workflows"
SMOKE = WORKFLOWS / "docker-smoke.yml"
CI_YML = WORKFLOWS / "ci.yml"
OLD_JOB = "docker-smoke"
READ_ONLY = {"read", "none"}


def _load(path: Path) -> dict[str, Any]:
    data: dict[str, Any] = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert isinstance(data, dict), path
    return data


def _needs(job: dict[str, Any]) -> list[str]:
    n = job.get("needs")
    return [n] if isinstance(n, str) else list(n or [])


def _script(job: dict[str, Any]) -> str:
    return "\n".join(str(s.get("run", "")) for s in job.get("steps") or [])


def _text(job: dict[str, Any]) -> str:
    return yaml.safe_dump(job, sort_keys=False, width=10_000)


@pytest.fixture(scope="module")
def wf() -> dict[str, Any]:
    return _load(SMOKE)


@pytest.fixture(scope="module")
def job_id(wf: dict[str, Any]) -> str:
    new = [j for j in (wf.get("jobs") or {}) if j != OLD_JOB]
    assert len(new) == 1, f"expected one new job beside `{OLD_JOB}` in docker-smoke.yml: {new}"
    return str(new[0])


@pytest.fixture(scope="module")
def job(wf: dict[str, Any], job_id: str) -> dict[str, Any]:
    j = wf["jobs"][job_id]
    assert isinstance(j, dict)
    return j


def test_the_m11_09_job_is_kept(wf: dict[str, Any]) -> None:
    assert OLD_JOB in (wf.get("jobs") or {})


def test_job_targets_bakta(job: dict[str, Any]) -> None:
    assert "BAKTA" in _text(job)


def test_job_goes_through_the_docker_profile(job: dict[str, Any]) -> None:
    script = _script(job)
    profiles = re.findall(r"(?:-profile|--profile)[ =]+(\S+)", script)
    assert any("docker" in p.split(",") for p in profiles), (
        f"no nextflow/nf-test call under a docker profile: {profiles}"
    )
    assert "db_dir" in _text(job), "the job must set params.db_dir"
    assert not re.search(r"docker\s+run\b[^\n]*\s-v\s", script), (
        "a hand-written docker run -v bypasses the profile's DB bind"
    )


def test_job_checks_the_db_is_readable_or_runs_bakta(job: dict[str, Any]) -> None:
    text = _text(job)
    readable = "version.json" in text
    runs_bakta = re.search(r"\bbakta\b.*--db", text) is not None or "light" in text.lower()
    assert readable or runs_bakta, "no DB-visibility check (test -r <db_dir>/version.json)"


def test_bakta_image_comes_from_containers_config(job: dict[str, Any]) -> None:
    text = _text(job)
    assert not re.search(r"oschwengers/bakta:[A-Za-z0-9]", text), "tagged bakta image"
    for ref in re.findall(r"oschwengers/bakta@sha256:[0-9a-f]+", text):
        assert ref in (REPO / "workflow/conf/containers.config").read_text(), ref


def test_advisory_like_the_existing_job(wf: dict[str, Any], job_id: str) -> None:
    for path in WORKFLOWS.glob("*.yml"):
        for name, j in (_load(path).get("jobs") or {}).items():
            assert job_id not in _needs(j or {}), f"{path.name}:{name} needs {job_id}"
    ci_jobs = _load(CI_YML).get("jobs") or {}
    if "summary" in ci_jobs:
        assert job_id not in _needs(ci_jobs["summary"])


def test_runner_timeout_and_least_privilege(wf: dict[str, Any], job: dict[str, Any]) -> None:
    assert str(job.get("runs-on", "")).startswith("ubuntu-")
    assert isinstance(job.get("timeout-minutes"), int) and job["timeout-minutes"] <= 30
    assert "services" not in job and "container" not in job
    for perms in (wf.get("permissions") or {}, job.get("permissions") or {}):
        assert isinstance(perms, dict)
        assert all(v in READ_ONLY for v in perms.values()), perms


def test_third_party_actions_pinned_by_sha(job: dict[str, Any]) -> None:
    for step in job.get("steps") or []:
        uses = step.get("uses")
        if uses and not str(uses).startswith("actions/"):
            assert re.fullmatch(r"[^@\s]+@[0-9a-f]{40}", str(uses)), uses
