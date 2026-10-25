"""``.github/workflows/docs.yml`` and the ``ci.yml`` lint steps (M11-07 item 4, acceptance 4,
example 7). Static: nothing is dispatched or deployed."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
DOCS_YML = REPO / ".github" / "workflows" / "docs.yml"
CI_YML = REPO / ".github" / "workflows" / "ci.yml"
SHA = re.compile(r"^[^@\s]+@[0-9a-f]{40}$")
GUARD = "vars.PAGES_ENABLED == 'true'"


@pytest.fixture(scope="module")
def wf() -> dict[str, Any]:
    assert DOCS_YML.is_file(), ".github/workflows/docs.yml does not exist (M11-07 item 4)"
    data: dict[str, Any] = yaml.safe_load(DOCS_YML.read_text(encoding="utf-8"))
    return data


def _on(wf: dict[str, Any]) -> dict[str, Any]:
    on = wf.get("on", wf.get(True))  # PyYAML reads a bare `on:` key as True
    assert isinstance(on, dict), on
    return on


def _script(job: dict[str, Any]) -> str:
    return "\n".join(str(s.get("run", "")) for s in job.get("steps", []))


def _uses(job: dict[str, Any]) -> list[str]:
    return [str(s["uses"]) for s in job.get("steps", []) if "uses" in s]


def _jobs(wf: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], str]:
    """(build job, deploy job, build job id): deploy is the job using actions/deploy-pages."""
    jobs: dict[str, Any] = wf["jobs"]
    deploy = [j for j in jobs.values() if any("actions/deploy-pages@" in u for u in _uses(j))]
    build = [
        (k, j)
        for k, j in jobs.items()
        if any("actions/upload-pages-artifact@" in u for u in _uses(j))
    ]
    assert len(deploy) == 1 and len(build) == 1, list(jobs)
    return build[0][1], deploy[0], build[0][0]


def test_triggers(wf: dict[str, Any]) -> None:
    on = _on(wf)
    assert "workflow_dispatch" in on
    push = on["push"]
    assert push["branches"] == ["main"]
    expected = {
        "docs/**",
        "mkdocs.yml",
        "src/giap/cli.py",
        "src/giap/config/**",
        "src/giap/schemas/**",
        "CHANGELOG.md",
    }
    assert expected <= set(push["paths"])
    assert "pull_request" not in on and "pull_request_target" not in on


def test_build_job_builds_strict_and_uploads(wf: dict[str, Any]) -> None:
    build, _, _ = _jobs(wf)
    s = _script(build)
    assert "mkdocs build --strict" in s
    assert "PAGES_ENABLED" not in str(build.get("if", ""))  # example 7: build always runs


def test_deploy_job_guard(wf: dict[str, Any]) -> None:
    _, deploy, build_id = _jobs(wf)
    cond = str(deploy.get("if", ""))
    assert GUARD in cond.replace('"', "'"), cond
    # deploy only from main, whatever triggered the run (#339 review B2)
    assert "github.ref == 'refs/heads/main'" in cond.replace('"', "'"), cond
    needs = deploy.get("needs")
    assert build_id in ([needs] if isinstance(needs, str) else list(needs or []))


def test_deploy_environment(wf: dict[str, Any]) -> None:
    _, deploy, _ = _jobs(wf)
    env = deploy.get("environment")
    name = env.get("name") if isinstance(env, dict) else env
    assert name == "github-pages"


def test_pages_permissions_only_on_deploy(wf: dict[str, Any]) -> None:
    build, deploy, _ = _jobs(wf)
    assert deploy.get("permissions", {}).get("pages") == "write"
    assert deploy.get("permissions", {}).get("id-token") == "write"
    top = wf.get("permissions")
    assert isinstance(top, dict), "set top-level permissions (least privilege)"
    for scope in (top, build.get("permissions") or {}):
        assert scope.get("pages") in (None, "none", "read")
        assert scope.get("id-token") in (None, "none")
        assert scope.get("contents") in (None, "none", "read")
    for job in wf["jobs"].values():
        if job is not deploy:
            assert (job.get("permissions") or {}).get("pages") != "write"


def test_actions_pinned_by_sha_with_tag_comment(wf: dict[str, Any]) -> None:
    """Q-L5: every third-party action in docs.yml is pinned by a 40-hex commit SHA (tag in a
    comment); GitHub's own ``actions/*`` may use a tag, as in ci.yml (M11-index Q-L5)."""
    uses = [u for j in wf["jobs"].values() for u in _uses(j) if not u.startswith("./")]
    assert uses
    third = [u for u in uses if not u.startswith("actions/")]
    assert [u for u in third if not SHA.match(u)] == []
    for line in DOCS_YML.read_text(encoding="utf-8").splitlines():
        if re.search(r"\buses:\s*[^.\s]", line) and not re.search(r"uses:\s*actions/", line):
            assert "#" in line.split("uses:", 1)[1], line


def test_no_tag_pushes(wf: dict[str, Any]) -> None:
    text = DOCS_YML.read_text(encoding="utf-8")
    assert "git push" not in text and "git tag" not in text


# --- ci.yml > lint ------------------------------------------------------------------------------


def test_ci_lint_checks_references_and_builds_strict() -> None:
    ci: dict[str, Any] = yaml.safe_load(CI_YML.read_text(encoding="utf-8"))
    s = _script(ci["jobs"]["lint"])
    assert "scripts/docs/gen_reference.py --check" in s
    assert "mkdocs build --strict --site-dir /tmp/site" in s
