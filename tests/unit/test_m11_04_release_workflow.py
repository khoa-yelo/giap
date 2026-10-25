"""``.github/workflows/release.yml`` parsed statically (M11-04 definition, acceptance 1, examples
1-9, edge cases). Nothing is dispatched, nothing touches the network.

Spec: docs/specs/M11-04-release-workflow.md. Pinning reading (Q-L5) as in
tests/unit/test_m11_07_docs_workflow.py: third-party actions by a 40-hex commit SHA with the tag in
a comment; GitHub's own ``actions/*`` may use a tag.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
RELEASE_YML = REPO / ".github" / "workflows" / "release.yml"
SHA = re.compile(r"^[^@\s]+@[0-9a-f]{40}$")
JOBS = ("verify", "build", "publish", "github-release")
WRITE_SCOPES = {"write", "write-all"}


@pytest.fixture(scope="module")
def text() -> str:
    assert RELEASE_YML.is_file(), ".github/workflows/release.yml does not exist (M11-04)"
    return RELEASE_YML.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def wf(text: str) -> dict[str, Any]:
    data: dict[str, Any] = yaml.safe_load(text)
    assert isinstance(data, dict)
    return data


@pytest.fixture(scope="module")
def jobs(wf: dict[str, Any]) -> dict[str, Any]:
    j: dict[str, Any] = wf["jobs"]
    missing = [k for k in JOBS if k not in j]
    assert missing == [], f"jobs missing: {missing}"
    return j


def _on(wf: dict[str, Any]) -> dict[str, Any]:
    on = wf.get("on", wf.get(True))  # PyYAML reads a bare `on:` key as True
    assert isinstance(on, dict), on
    return on


def _steps(job: dict[str, Any]) -> list[dict[str, Any]]:
    return list(job.get("steps") or [])


def _script(job: dict[str, Any]) -> str:
    return "\n".join(str(s.get("run", "")) for s in _steps(job))


def _job_text(job: dict[str, Any]) -> str:
    """Everything a job says (runs, withs, envs, ifs) as one string."""
    return yaml.safe_dump(job, sort_keys=False, width=10_000)


def _uses(job: dict[str, Any]) -> list[str]:
    return [str(s["uses"]) for s in _steps(job) if "uses" in s]


def _needs(job: dict[str, Any]) -> list[str]:
    n = job.get("needs")
    return [n] if isinstance(n, str) else list(n or [])


def _cond(job: dict[str, Any]) -> str:
    """Normalise a job ``if:``: no ``${{ }}``, single quotes, no blanks after ``!``, one space."""
    c = str(job.get("if", "")).strip()
    if c.startswith("${{") and c.endswith("}}"):
        c = c[3:-2]
    c = c.replace('"', "'")
    c = re.sub(r"!\s+", "!", c)
    return " ".join(c.split())


def _perms(scope: dict[str, Any]) -> dict[str, Any]:
    p = scope.get("permissions")
    if p is None:
        return {}
    assert isinstance(p, dict), f"permissions must be a mapping, got {p!r}"
    return p


# --- trigger and inputs -------------------------------------------------------------------------


def test_only_workflow_dispatch(wf: dict[str, Any]) -> None:
    assert set(_on(wf)) == {"workflow_dispatch"}


def test_inputs_and_defaults(wf: dict[str, Any]) -> None:
    inputs = _on(wf)["workflow_dispatch"]["inputs"]
    assert set(inputs) == {"version", "dry_run", "target"}
    v = inputs["version"]
    assert v.get("type") == "string" and v.get("required") is True
    assert v.get("default") in (None, "")
    d = inputs["dry_run"]
    assert d.get("type") == "boolean"
    assert d.get("default") is True, "dry_run must default to true"
    t = inputs["target"]
    assert t.get("type") == "choice"
    assert set(t.get("options") or []) == {"testpypi", "pypi"}
    assert t.get("default") == "testpypi"


def test_concurrency_release_not_cancelled(wf: dict[str, Any]) -> None:
    c = wf.get("concurrency")
    assert isinstance(c, dict), c
    assert c.get("group") == "release"
    assert c.get("cancel-in-progress") is False


# --- permissions --------------------------------------------------------------------------------


def test_top_level_contents_read_only(wf: dict[str, Any]) -> None:
    assert _perms(wf) == {"contents": "read"}


def test_id_token_write_only_on_publish(wf: dict[str, Any], jobs: dict[str, Any]) -> None:
    assert _perms(jobs["publish"]).get("id-token") == "write"
    assert _perms(wf).get("id-token") in (None, "none")
    for name, job in jobs.items():
        if name != "publish":
            assert _perms(job).get("id-token") in (None, "none"), name


def test_contents_write_only_on_github_release(wf: dict[str, Any], jobs: dict[str, Any]) -> None:
    assert _perms(jobs["github-release"]).get("contents") == "write"
    for name, job in jobs.items():
        if name != "github-release":
            assert _perms(job).get("contents") in (None, "none", "read"), name


def test_no_other_write_scopes(wf: dict[str, Any], jobs: dict[str, Any]) -> None:
    allowed = {("publish", "id-token"), ("github-release", "contents")}
    assert wf.get("permissions") not in ("write-all", "read-all")
    for name, job in jobs.items():
        assert job.get("permissions") not in ("write-all", "read-all"), name
        for scope, level in _perms(job).items():
            if level in WRITE_SCOPES:
                assert (name, scope) in allowed, (name, scope)


# --- job graph and guards -----------------------------------------------------------------------


def test_needs_chain(jobs: dict[str, Any]) -> None:
    assert _needs(jobs["verify"]) == []
    assert "verify" in _needs(jobs["build"])
    assert "build" in _needs(jobs["publish"])
    assert "publish" in _needs(jobs["github-release"])


@pytest.mark.parametrize("name", ["publish", "github-release"])
def test_publishing_jobs_guarded_by_dry_run_and_actor(jobs: dict[str, Any], name: str) -> None:
    """Examples 4 and 1: another actor or a dry run skips both publishing jobs."""
    c = _cond(jobs[name])
    assert "!inputs.dry_run" in c, c
    assert "github.actor == 'khoa-yelo'" in c, c
    assert "||" not in c, f"an OR can bypass the guard: {c}"
    # exact conjuncts, so an inverted or unconditional guard cannot pass (#342 review B1)
    want = {"!inputs.dry_run", "github.actor == 'khoa-yelo'"}
    if name == "github-release":
        want.add("inputs.target == 'pypi'")
    assert {part.strip() for part in c.split("&&")} == want, c
    for bad in ("always()", "failure()", "cancelled()", "== false", "!= true", "== 'false'"):
        assert bad not in c, (bad, c)


def test_github_release_only_for_pypi(jobs: dict[str, Any]) -> None:
    """Example 5: a TestPyPI upload creates no tag and no GitHub release."""
    c = _cond(jobs["github-release"])
    assert "inputs.target == 'pypi'" in c, c


def test_verify_and_build_run_for_any_actor(jobs: dict[str, Any]) -> None:
    for name in ("verify", "build"):
        c = _cond(jobs[name])
        assert "github.actor" not in c and "dry_run" not in c, (name, c)


def test_publish_environment_is_the_target(jobs: dict[str, Any]) -> None:
    env = jobs["publish"].get("environment")
    name = env.get("name") if isinstance(env, dict) else env
    assert " ".join(str(name).split()) == "${{ inputs.target }}", env


# --- verify -------------------------------------------------------------------------------------


def test_verify_checkout_full_history_at_sha(jobs: dict[str, Any]) -> None:
    co = [
        s for s in _steps(jobs["verify"]) if str(s.get("uses", "")).startswith("actions/checkout@")
    ]
    assert co, "verify checks out the repository"
    w = co[0].get("with") or {}
    assert w.get("fetch-depth") == 0
    assert "github.sha" in str(w.get("ref", ""))


def test_verify_refuses_non_main(jobs: dict[str, Any]) -> None:
    """Example 8."""
    t = _job_text(jobs["verify"])
    assert "refs/heads/main" in t
    assert "release runs only on main" in t


def test_verify_checks_version(jobs: dict[str, Any]) -> None:
    """Example 3: ``version input 1.0.0 != giap 0.0.1``; the input regex is enforced."""
    t = _job_text(jobs["verify"])
    assert "__version__" in t
    assert "version input" in t and "!= giap" in t
    assert "inputs.version" in t
    assert re.search(r"rc(\\d|\[0-9\])\+", t), "verify enforces ^\\d+\\.\\d+\\.\\d+(rc\\d+)?$"


def test_verify_checks_tag_absent(jobs: dict[str, Any]) -> None:
    """Example 9: ``tag v1.0.0 exists``."""
    t = _job_text(jobs["verify"])
    assert "git ls-remote --tags origin" in t
    assert re.search(r"tag v[^\n]*\bexists\b", t), "verify prints `tag v<version> exists`"


def test_verify_runs_readiness(jobs: dict[str, Any]) -> None:
    """Examples 1, 2 and 7: readiness.py verify on validation/release/<version>/readiness.json,
    missing is reported as ``readiness missing``, the overall is read and gated per target."""
    t = _job_text(jobs["verify"])
    assert "scripts/release/readiness.py verify" in t
    assert "validation/release/" in t and "readiness.json" in t
    assert "readiness missing" in t
    assert "overall" in t
    assert "inputs.dry_run" in t and "inputs.target" in t


# --- build --------------------------------------------------------------------------------------


def test_build_builds_checks_and_smokes(jobs: dict[str, Any]) -> None:
    s = _script(jobs["build"])
    assert "uv build" in s
    assert re.search(r"check_dist\.py dist --version\b", s), s
    assert ".whl" in s and "giap --version" in s, "the M11-01 clean-venv smoke"


def test_build_renders_notes_with_draft_on_dry_run(jobs: dict[str, Any]) -> None:
    t = _job_text(jobs["build"])
    assert "scripts/release/notes.py" in t
    assert "--draft" in t and "inputs.dry_run" in t


def test_build_uploads_dist_and_notes(jobs: dict[str, Any]) -> None:
    ups = [
        s
        for s in _steps(jobs["build"])
        if str(s.get("uses", "")).startswith("actions/upload-artifact@")
    ]
    paths = "\n".join(str((s.get("with") or {}).get("path", "")) for s in ups)
    assert "dist" in paths and "notes.md" in paths, paths


# --- publish ------------------------------------------------------------------------------------


def test_publish_uses_pinned_pypi_action(jobs: dict[str, Any]) -> None:
    pub = [u for u in _uses(jobs["publish"]) if u.startswith("pypa/gh-action-pypi-publish@")]
    assert len(pub) == 1, _uses(jobs["publish"])
    assert SHA.match(pub[0]), pub[0]


def test_publish_repository_url_for_testpypi(jobs: dict[str, Any]) -> None:
    step = next(
        s
        for s in _steps(jobs["publish"])
        if str(s.get("uses", "")).startswith("pypa/gh-action-pypi-publish@")
    )
    w = step.get("with") or {}
    assert "repository-url" in w
    assert "test.pypi.org/legacy/" in str(w["repository-url"])
    assert "inputs.target" in str(w["repository-url"])
    assert "password" not in w, "trusted publishing needs no token (R12)"


def test_publish_uploads_build_artefacts_only(jobs: dict[str, Any]) -> None:
    pub = jobs["publish"]
    assert any(u.startswith("actions/download-artifact@") for u in _uses(pub))
    s = _script(pub)
    assert "uv build" not in s and "python -m build" not in s


def test_publish_skip_existing_only_on_matching_sha256(jobs: dict[str, Any]) -> None:
    """Edge case: re-run after a failed publish compares sha256 with the index JSON first."""
    t = _job_text(jobs["publish"])
    assert "skip-existing" in t
    assert "pypi.org/pypi/giap/" in t and "sha256" in t


# --- github-release -----------------------------------------------------------------------------


def test_github_release_command(jobs: dict[str, Any]) -> None:
    s = _script(jobs["github-release"])
    assert re.search(r"gh release create [\"']?v", s), s
    for frag in ("dist/*", "--target", "github.sha", "--title", "GIAP ", "--notes-file notes.md"):
        assert frag in s, frag
    assert "--prerelease" in s
    assert any(u.startswith("actions/download-artifact@") for u in _uses(jobs["github-release"]))


def test_release_created_only_in_github_release(text: str, jobs: dict[str, Any]) -> None:
    """The only place a tag is created (§3.10)."""
    assert text.count("gh release create") == 1
    for name, job in jobs.items():
        if name != "github-release":
            assert "gh release" not in _job_text(job), name


# --- summary ------------------------------------------------------------------------------------


def test_step_summary(text: str) -> None:
    """Item 5 and examples 1 and 4."""
    assert "GITHUB_STEP_SUMMARY" in text
    assert "sha256" in text
    assert "nothing published" in text
    assert "is not khoa-yelo" in text


# --- hygiene (acceptance 1) ---------------------------------------------------------------------


def test_actions_pinned_by_sha_with_tag_comment(wf: dict[str, Any], text: str) -> None:
    uses = [u for j in wf["jobs"].values() for u in _uses(j) if not u.startswith("./")]
    assert uses
    third = [u for u in uses if not u.startswith("actions/")]
    assert third, "pypa/gh-action-pypi-publish at least"
    assert [u for u in third if not SHA.match(u)] == []
    for line in text.splitlines():
        if re.search(r"\buses:\s*[^.\s]", line) and not re.search(r"uses:\s*actions/", line):
            assert "#" in line.split("uses:", 1)[1], line


def test_no_secrets(text: str) -> None:
    assert "secrets." not in text and "secrets[" not in text


def test_no_git_push_or_tag(text: str) -> None:
    assert "git push" not in text and "git tag" not in text
