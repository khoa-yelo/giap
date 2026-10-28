"""The ``docker-smoke`` Actions job parsed statically (M11-09 definition item 2, examples 1-3,
edge cases, acceptance 3). Nothing is dispatched; no network, no Docker.

Spec: docs/specs/M11-09-platforms-docker-smoke.md. Readings used here:

- The job is found by its id ``docker-smoke`` in any ``.github/workflows/*.yml`` (the spec names
  ci.yml; a dedicated workflow file with native ``paths:`` is accepted too).
- Trigger (Q-K10 default): path-filtered pull requests (``workflow/**``) plus a weekly ``schedule``;
  never on every push. The paths filter is either the workflow's ``on.pull_request.paths`` or a
  ``dorny/paths-filter`` step (in the job or in a job it needs) whose output gates the job.
- Advisory: not in ``summary.needs`` and no other job needs it.
- Pinning (Q-L5): third-party actions by a 40-hex commit SHA; GitHub's own ``actions/*`` may use
  a tag. Images: MINIPROT/BLAST/TBLASTN pinned by ``@sha256`` in containers.config, and every
  literal image the job names carries a digest.
"""

from __future__ import annotations

import re
import shlex
from pathlib import Path
from typing import Any

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
WORKFLOWS = REPO / ".github" / "workflows"
CI_YML = WORKFLOWS / "ci.yml"
CONTAINERS = REPO / "workflow" / "conf" / "containers.config"
SMOKE_TEST = REPO / "workflow" / "tests" / "docker_smoke.nf.test"
JOB = "docker-smoke"
SHA = re.compile(r"^[^@\s]+@[0-9a-f]{40}$")
DIGEST = re.compile(r"@sha256:[0-9a-f]{64}\b")
READ_ONLY = {"read", "none"}
# An image reference with a tag (registry/name:tag or owner/name:tag), not followed by a digest.
TAGGED_IMAGE = re.compile(
    r"(?:docker://)?\b((?:[a-z0-9-]+\.)+[a-z]{2,}(?::\d+)?/[a-z0-9._/-]+|[a-z0-9-]+/[a-z0-9._-]+)"
    r":([A-Za-z0-9][A-Za-z0-9._-]*)(?!\S*@sha256)"
)


# --- helpers ------------------------------------------------------------------------------------


def _load(path: Path) -> dict[str, Any]:
    data: dict[str, Any] = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert isinstance(data, dict), path
    return data


def _on(wf: dict[str, Any]) -> dict[str, Any]:
    on = wf.get("on", wf.get(True))  # PyYAML reads a bare `on:` key as True
    if isinstance(on, str):
        return {on: None}
    if isinstance(on, list):
        return dict.fromkeys(on)
    assert isinstance(on, dict), on
    return on


def _steps(job: dict[str, Any]) -> list[dict[str, Any]]:
    return list(job.get("steps") or [])


def _script(job: dict[str, Any]) -> str:
    return "\n".join(str(s.get("run", "")) for s in _steps(job))


def _uses(job: dict[str, Any]) -> list[str]:
    return [str(s["uses"]) for s in _steps(job) if "uses" in s]


def _needs(job: dict[str, Any]) -> list[str]:
    n = job.get("needs")
    return [n] if isinstance(n, str) else list(n or [])


def _job_text(job: dict[str, Any]) -> str:
    return yaml.safe_dump(job, sort_keys=False, width=10_000)


def _glob_re(pattern: str) -> re.Pattern[str]:
    """GitHub path-filter glob: ``**`` crosses ``/``, ``*`` and ``?`` do not."""
    out, i = [], 0
    while i < len(pattern):
        if pattern.startswith("**/", i):
            out.append("(?:.*/)?")
            i += 3
        elif pattern.startswith("**", i):
            out.append(".*")
            i += 2
        elif pattern[i] == "*":
            out.append("[^/]*")
            i += 1
        elif pattern[i] == "?":
            out.append("[^/]")
            i += 1
        else:
            out.append(re.escape(pattern[i]))
            i += 1
    return re.compile("".join(out) + r"\Z")


def _matches(patterns: list[str], path: str) -> bool:
    hit = False
    for p in patterns:
        if p.startswith("!"):
            if _glob_re(p[1:]).match(path):
                hit = False
        elif _glob_re(p).match(path):
            hit = True
    return hit


def _filter_patterns(filters: Any) -> list[str]:
    """Every pattern in a dorny/paths-filter ``filters`` value (YAML string or mapping)."""
    data = yaml.safe_load(filters) if isinstance(filters, str) else filters
    out: list[str] = []

    def walk(x: Any) -> None:
        if isinstance(x, str):
            out.append(x)
        elif isinstance(x, list):
            for y in x:
                walk(y)
        elif isinstance(x, dict):
            for y in x.values():
                walk(y)

    walk(data)
    return out


@pytest.fixture(scope="module")
def located() -> tuple[Path, dict[str, Any], dict[str, Any]]:
    """(workflow file, workflow, job) for the single job with id ``docker-smoke``."""
    hits = []
    for path in sorted(WORKFLOWS.glob("*.y*ml")):
        wf = _load(path)
        if JOB in (wf.get("jobs") or {}):
            hits.append((path, wf, wf["jobs"][JOB]))
    assert len(hits) == 1, f"expected one `{JOB}` job in .github/workflows, found {len(hits)}"
    return hits[0]


@pytest.fixture(scope="module")
def wf(located: tuple[Path, dict[str, Any], dict[str, Any]]) -> dict[str, Any]:
    return located[1]


@pytest.fixture(scope="module")
def job(located: tuple[Path, dict[str, Any], dict[str, Any]]) -> dict[str, Any]:
    return located[2]


def _filter_jobs(wf: dict[str, Any], job: dict[str, Any]) -> list[dict[str, Any]]:
    """The job itself and the jobs it needs (where a paths-filter step may live)."""
    jobs: dict[str, Any] = wf["jobs"]
    return [job, *(jobs[n] for n in _needs(job) if n in jobs)]


def _pr_patterns(wf: dict[str, Any], job: dict[str, Any]) -> list[str]:
    pr = _on(wf).get("pull_request")
    if isinstance(pr, dict) and pr.get("paths"):
        return [str(p) for p in pr["paths"]]
    for j in _filter_jobs(wf, job):
        for s in _steps(j):
            if "paths-filter" in str(s.get("uses", "")):
                return _filter_patterns((s.get("with") or {}).get("filters"))
    raise AssertionError(
        "no paths filter for pull requests: set on.pull_request.paths or a paths-filter step"
    )


# --- trigger (definition 2, examples 1-2, Q-K10) ------------------------------------------------


def test_runs_on_pull_requests(wf: dict[str, Any]) -> None:
    on = _on(wf)
    assert "pull_request" in on
    assert "pull_request_target" not in on  # forks never get a write token


@pytest.mark.parametrize(
    "path",
    [
        "workflow/conf/containers.config",  # example 1
        "workflow/modules/local/miniprot/main.nf",
        "workflow/nextflow.config",
        "workflow/tests/docker_smoke.nf.test",
    ],
)
def test_paths_filter_selects_workflow_changes(
    wf: dict[str, Any], job: dict[str, Any], path: str
) -> None:
    assert _matches(_pr_patterns(wf, job), path), path


@pytest.mark.parametrize(
    "path",
    ["src/gippa/stats/enrichment.py", "src/gippa/stats/__init__.py", "docs/index.md", "README.md"],
)
def test_paths_filter_skips_other_changes(
    wf: dict[str, Any], job: dict[str, Any], path: str
) -> None:
    """Example 2: a PR touching only src/gippa/stats/ skips the smoke."""
    assert not _matches(_pr_patterns(wf, job), path), path


def test_in_workflow_filter_gates_the_job(wf: dict[str, Any], job: dict[str, Any]) -> None:
    """Without native ``paths:``, the filter's output must gate the job (or its smoke step)."""
    pr = _on(wf).get("pull_request")
    if isinstance(pr, dict) and pr.get("paths"):
        return
    conds = [str(job.get("if", ""))] + [str(s.get("if", "")) for s in _steps(job)]
    assert any("outputs" in c for c in conds), conds


def test_weekly_schedule(wf: dict[str, Any]) -> None:
    sched = _on(wf).get("schedule")
    assert isinstance(sched, list) and sched, "add a weekly `schedule` trigger (Q-K10)"
    crons = [str(s["cron"]).split() for s in sched]
    assert len(crons) == 1, crons
    minute, hour, dom, month, dow = crons[0]
    assert (dom, month) == ("*", "*"), crons[0]
    assert re.fullmatch(r"[0-7]|SUN|MON|TUE|WED|THU|FRI|SAT", dow.upper()), crons[0]
    assert minute.isdigit() and hour.isdigit(), crons[0]  # once a week, not hourly


def test_not_on_every_push(wf: dict[str, Any], job: dict[str, Any]) -> None:
    """Cheapest reading of Q-K10: the smoke does not run on each push to main. A workflow that is
    also triggered by other events must gate the job on ``github.event_name``."""
    on = set(_on(wf))
    extra = on - {"pull_request", "schedule", "workflow_dispatch"}
    if extra:
        cond = str(job.get("if", ""))
        assert "github.event_name" in cond, (sorted(extra), cond)
        assert not re.search(r"==\s*['\"]push['\"]", cond), cond


# --- the run itself (definition 2) --------------------------------------------------------------


def test_ubuntu_runner_and_timeout(job: dict[str, Any]) -> None:
    assert job["runs-on"] == "ubuntu-latest"
    t = job.get("timeout-minutes")
    assert isinstance(t, int) and 0 < t <= 10, t


def test_java_21_and_nextflow_as_in_nf_test_stub(job: dict[str, Any]) -> None:
    java = [s for s in _steps(job) if str(s.get("uses", "")).startswith("actions/setup-java@")]
    assert len(java) == 1
    assert str(java[0]["with"]["java-version"]) == "21"
    s = _script(job)
    assert "get.nextflow.io" in s
    assert "get.nf-test.com" in s


def _nf_test_profile(job: dict[str, Any]) -> tuple[str, set[str]]:
    cmds = [ln for ln in _script(job).replace("\\\n", " ").splitlines() if "nf-test test" in ln]
    assert len(cmds) == 1, cmds
    words = shlex.split(cmds[0].split("nf-test test", 1)[1])
    for i, w in enumerate(words):
        if w == "--profile":
            return cmds[0], set(words[i + 1].split(","))
        if w.startswith("--profile="):
            return cmds[0], set(w.split("=", 1)[1].split(","))
    raise AssertionError(f"no --profile in {cmds[0]!r}")


def test_runs_the_smoke_suite_under_profile_test_docker(job: dict[str, Any]) -> None:
    cmd, profiles = _nf_test_profile(job)
    assert "workflow/tests/docker_smoke.nf.test" in cmd
    assert profiles == {"test", "docker"}, profiles  # not stub, not apptainer, not conda


def test_no_job_container_or_services(job: dict[str, Any]) -> None:
    assert "container" not in job and "services" not in job


# --- images by digest (example 3) ---------------------------------------------------------------


@pytest.mark.parametrize("process", ["MINIPROT", "BLAST", "TBLASTN"])
def test_smoke_images_pinned_by_digest(process: str) -> None:
    text = CONTAINERS.read_text(encoding="utf-8")
    m = re.search(rf"withName:\s*'{process}'\s*\{{\s*container\s*=\s*'([^']+)'", text)
    assert m, process
    assert DIGEST.search(m.group(1)), m.group(1)


def test_job_names_no_tagged_image(job: dict[str, Any]) -> None:
    text = _job_text(job)
    tagged = [m.group(0) for m in TAGGED_IMAGE.finditer(text) if "@" not in m.group(0)]
    assert tagged == [], tagged


def test_pull_retried_once_and_rate_limit_reported(job: dict[str, Any]) -> None:
    """Edge case: ``docker pull`` with one retry, then ``rate-limited`` in the step summary."""
    s = _script(job)
    assert "docker pull" in s
    assert "rate-limited" in s
    assert "GITHUB_STEP_SUMMARY" in s
    for line in s.splitlines():
        m = re.search(r"docker pull\s+(\S+)", line)
        if m and not m.group(1).startswith(("$", '"$', "'$", "-")):
            assert DIGEST.search(m.group(1)), line  # a literal pull is by digest


# --- least privilege and action pins (Q-L5) -----------------------------------------------------


def test_least_privilege(wf: dict[str, Any], job: dict[str, Any]) -> None:
    perms = job.get("permissions", wf.get("permissions"))
    assert isinstance(perms, dict), "set explicit read-only permissions on the job"
    assert perms.get("contents") in ("read", "none"), perms
    bad = {k: v for k, v in perms.items() if str(v) not in READ_ONLY}
    assert bad == {}, bad
    assert "secrets." not in _job_text(job)


def test_third_party_actions_pinned_by_sha(wf: dict[str, Any], job: dict[str, Any]) -> None:
    uses = [u for j in _filter_jobs(wf, job) for u in _uses(j) if not u.startswith("./")]
    assert uses
    third = [u for u in uses if not u.startswith("actions/")]
    assert [u for u in third if not SHA.match(u)] == []
    assert not any(u.startswith("docker://") for u in uses)


# --- advisory, and kept out of the stub run (acceptance 3) --------------------------------------


def test_not_required_by_summary_or_any_job() -> None:
    ci = _load(CI_YML)
    jobs: dict[str, Any] = ci["jobs"]
    assert JOB not in _needs(jobs["summary"])
    assert [n for n, j in jobs.items() if JOB in _needs(j)] == []


def test_stub_run_does_not_pick_up_the_smoke_suite() -> None:
    """``nf-test test workflow/tests --profile test,stub`` would also run docker_smoke.nf.test
    (real tools, no containers): the stub job must exclude it by name or by tag."""
    s = _script(_load(CI_YML)["jobs"]["nf-test-stub"])
    assert "docker_smoke" in s or "--tag" in s, s


# --- the nf-test suite it runs ------------------------------------------------------------------


@pytest.fixture(scope="module")
def smoke_text() -> str:
    assert SMOKE_TEST.is_file(), "workflow/tests/docker_smoke.nf.test does not exist (M11-09)"
    return SMOKE_TEST.read_text(encoding="utf-8")


def test_smoke_suite_runs_real_miniprot_and_blast(smoke_text: str) -> None:
    assert re.search(r'process\s+"MINIPROT"', smoke_text)
    # nf-test runs one suite per file: the second process runs as a setup dependency (lead)
    assert re.search(r'(?:process\s+|run\()"(?:TBLASTN|BLAST)"', smoke_text)
    assert '"-stub"' not in smoke_text and "-stub-run" not in smoke_text


def test_smoke_suite_asserts_outputs(smoke_text: str) -> None:
    assert "process.success" in smoke_text
    assert smoke_text.count("assert") >= 4  # success, non-empty and header for each process
