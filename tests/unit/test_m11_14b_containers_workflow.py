"""``containers.yml`` with key tags, the authenticated check and the publish gate (M11-14b).

Spec docs/specs/M11-14-container-install.md item 2 ("Immutable tags"), item 5 (release flow),
item 6 (prefix from ``pyproject.toml``), examples 5, 5b, 5c, AC 2; ADR-0060 §1-§2, §5.
Nothing is built, pushed or queried: the workflow is parsed as YAML, and the existence check is
run against a fake ``docker`` on ``PATH``.

Interface fixed by these tests (readings chosen where the spec leaves room):

* The existence check of one ref lives in ``scripts/release/ghcr_tag_state.sh REF``. It runs
  ``docker buildx imagetools inspect REF`` and classifies the result (item 2): inspect exit 0 ->
  prints ``exists`` and exits 0; an error whose output matches
  ``manifest unknown|not found|NAME_UNKNOWN`` -> prints ``absent`` and exits 0; any other error ->
  exits non-zero and repeats the error on stderr (never ``absent``). Both the plan job and the
  build job call this script, so the two checks cannot classify differently (AC 2 asks for the
  workflow's own shell to be run against a fake ``docker``; the shared script is that shell).
* The plan job is the job whose output feeds the build matrix (``fromJSON(needs.<plan>...)``).
  It computes tags with ``scripts/release/image_key.py``, logs in to GHCR before the check, has
  ``packages: read``, prints ``publishing off`` when ``vars.GHCR_PUBLISH`` is not ``'true'``, and
  exposes the tools whose tag already exists as an output whose name contains ``exist`` (their
  digests go into the digest PR, example 5); the digest-PR job reads that output.
* The build job's re-check is the step right before ``docker/build-push-action`` and fails with
  ``tag exists: <ref>``.
* The image prefix is not written in the workflow: no ``gippa-`` literal appears in any ``run`` or
  ``with`` value; ``pyproject.toml`` is read in a job the build needs.
"""

from __future__ import annotations

import functools
import os
import re
import shutil
import stat
import subprocess
import tomllib
from pathlib import Path
from typing import Any

import pytest
import yaml
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
WF = REPO / ".github" / "workflows" / "containers.yml"
CHECK = REPO / "scripts" / "release" / "ghcr_tag_state.sh"
REF = "ghcr.io/khoa-yelo/gippa-ostir:1.1.3-3f9c0a1b2d4e"
ABSENT_PATTERNS = ("manifest unknown", "not found", "NAME_UNKNOWN")


@functools.cache
def _text() -> str:
    assert WF.is_file(), ".github/workflows/containers.yml is missing"
    return WF.read_text()


@functools.cache
def _load() -> dict[str, Any]:
    data: dict[str, Any] = yaml.safe_load(_text())
    return data


def _jobs() -> dict[str, Any]:
    jobs: dict[str, Any] = _load()["jobs"]
    return jobs


def _steps(job: dict[str, Any]) -> list[dict[str, Any]]:
    return list(job.get("steps", []))


def _runs(job: dict[str, Any]) -> str:
    return "\n".join(str(s.get("run", "")) for s in _steps(job))


def _needs(job: dict[str, Any]) -> list[str]:
    n = job.get("needs", [])
    return [n] if isinstance(n, str) else list(n)


def _perms(block: Any) -> dict[str, str]:
    if block is None:
        return {}
    if isinstance(block, str):
        return {"*": block}
    return {str(k): str(v) for k, v in block.items()}


def _build() -> tuple[str, dict[str, Any], int]:
    hits = [
        (name, job, i)
        for name, job in _jobs().items()
        for i, s in enumerate(_steps(job))
        if str(s.get("uses", "")).startswith("docker/build-push-action@")
    ]
    assert len(hits) == 1, "exactly one docker/build-push-action step expected"
    return hits[0]


def _plan() -> tuple[str, dict[str, Any]]:
    _, job, _ = _build()
    matrix = yaml.safe_dump(job.get("strategy", {}).get("matrix", {}), width=10_000)
    m = re.search(r"fromJSON\(\s*needs\.([\w-]+)\.outputs\.", matrix)
    assert m, "the build matrix must come from the plan job (fromJSON(needs.<plan>.outputs.x))"
    return m.group(1), _jobs()[m.group(1)]


def _digest_job() -> tuple[str, dict[str, Any]]:
    hits = [(n, j) for n, j in _jobs().items() if "scripts/release/update_digests.py" in _runs(j)]
    assert len(hits) == 1, "exactly one job runs scripts/release/update_digests.py"
    return hits[0]


def _calls_check(step: dict[str, Any]) -> bool:
    return "scripts/release/ghcr_tag_state.sh" in str(step.get("run", ""))


def _login_index(job: dict[str, Any]) -> int:
    for i, s in enumerate(_steps(job)):
        if str(s.get("uses", "")).startswith("docker/login-action@"):
            assert str(s["with"]["registry"]).strip() == "ghcr.io"
            pw = str(s["with"].get("password", ""))
            assert re.search(r"secrets\.GITHUB_TOKEN|github\.token", pw), pw
            return i
    raise AssertionError("no docker/login-action step in the job")


def _tags() -> list[str]:
    _, job, i = _build()
    raw = str(_steps(job)[i]["with"]["tags"])
    return [t.strip() for t in re.split(r"[,\n]", raw) if t.strip()]


# --- immutable tags ---------------------------------------------------------------------------


def test_only_one_tag_is_pushed() -> None:
    assert len(_tags()) == 1, _tags()


def test_the_tag_is_the_key_tag_not_sha7_or_bare_version() -> None:
    (t,) = _tags()
    assert "sha7" not in t and "github.sha" not in t, t
    assert not re.search(r"outputs\.version\s*\}\}\s*$", t), f"bare `:<version>` tag: {t}"
    # either `<image>:${{ …tag… }}` or one expression naming the full ref
    assert re.search(r":\$\{\{[^}]*\btag\b[^}]*\}\}$", t) or re.fullmatch(
        r"\$\{\{[^}]*\b(ref|tag)\b[^}]*\}\}", t
    ), t


def test_plan_computes_tags_with_image_key_py() -> None:
    _, plan = _plan()
    assert "scripts/release/image_key.py" in _runs(plan)


# --- authenticated existence check (examples 5, 5b) -------------------------------------------


def test_plan_job_has_packages_read_and_nothing_written() -> None:
    _, plan = _plan()
    perms = _perms(plan.get("permissions"))
    assert perms.get("packages") == "read", perms
    assert "write" not in perms.values() and perms.get("*") != "write-all", perms


def test_plan_job_logs_in_before_the_check() -> None:
    _, plan = _plan()
    login = _login_index(plan)
    checks = [i for i, s in enumerate(_steps(plan)) if _calls_check(s)]
    assert checks, "the plan job never runs scripts/release/ghcr_tag_state.sh"
    assert login < min(checks)


def test_existing_tools_feed_the_digest_pr() -> None:
    """Example 5: a tool whose tag exists is not built, and its digest goes into the digest PR."""
    plan_name, plan = _plan()
    outs = [o for o in plan.get("outputs", {}) if "exist" in o]
    assert outs, "the plan job needs an output listing the tools whose tag exists"
    _, job = _digest_job()
    assert plan_name in _needs(job)
    blob = yaml.safe_dump(job, width=10_000)
    assert any(f"needs.{plan_name}.outputs.{o}" in blob for o in outs), outs


# --- re-check before push (example 5c) --------------------------------------------------------


def test_build_job_rechecks_right_before_the_push() -> None:
    _, job, i = _build()
    assert i > 0 and _calls_check(_steps(job)[i - 1]), (
        "the step right before docker/build-push-action must run ghcr_tag_state.sh"
    )
    assert _login_index(job) < i - 1
    assert "tag exists: " in str(_steps(job)[i - 1]["run"])


# --- race, gate and environment --------------------------------------------------------------


def test_global_concurrency_group_without_cancel() -> None:
    conc = _load()["concurrency"]
    assert isinstance(conc, dict), conc
    assert str(conc["group"]).strip() == "containers", conc
    assert conc.get("cancel-in-progress") is False
    for name, job in _jobs().items():
        assert "concurrency" not in job, f"job {name} must not override the global group"


def test_build_job_is_gated_on_ghcr_publish() -> None:
    _, job, _ = _build()
    cond = str(job.get("if", ""))
    assert re.search(r"vars\.GHCR_PUBLISH\s*==\s*'true'", cond), cond


def test_build_job_runs_in_the_ghcr_environment() -> None:
    _, job, _ = _build()
    env = job.get("environment")
    name = env.get("name") if isinstance(env, dict) else env
    assert str(name).strip() == "ghcr", env


def test_digest_pr_job_is_also_gated_on_ghcr_publish() -> None:
    """Publishing off: stop after plan; the digest-PR job must not run on a skipped build."""
    _, job = _digest_job()
    assert re.search(r"vars\.GHCR_PUBLISH\s*==\s*'true'", str(job.get("if", ""))), job.get("if")


def test_plan_says_publishing_off() -> None:
    _, plan = _plan()
    blob = yaml.safe_dump(plan, width=10_000)
    assert "vars.GHCR_PUBLISH" in blob
    assert "publishing off" in _runs(plan)


def test_the_workflow_never_sets_ghcr_publish_or_visibility() -> None:
    text = _text()
    assert not re.search(r"gh\s+variable\s+set", text)
    assert not re.search(r"gh\s+api[^\n]*(visibility|environments)", text)
    assert "imagetools create" not in text  # no retagging (ADR-0060 §2)
    assert not re.search(r"gh\s+api[^\n]*-X\s*DELETE", text)


# --- prefix from pyproject (item 6) ------------------------------------------------------------


def test_prefix_is_not_hard_coded() -> None:
    for name, job in _jobs().items():
        for s in _steps(job):
            body = str(s.get("run", "")) + yaml.safe_dump(s.get("with", {}), width=10_000)
            assert not re.search(r"\bgippa-", body), f"{name}: literal `gippa-` image prefix"
    assert not re.search(r"\bgippa-", yaml.safe_dump(_load().get("env", {}), width=10_000))


def test_prefix_is_read_from_pyproject_in_a_job_the_build_needs() -> None:
    _, build, _ = _build()
    readers = [n for n, j in _jobs().items() if "pyproject.toml" in _runs(j)]
    assert set(readers) & set(_needs(build)), readers


def test_digest_pr_commits_images_json() -> None:
    _, job = _digest_job()
    assert "workflow/conf/images.json" in _runs(job)


def test_every_ghcr_ref_uses_the_project_prefix() -> None:
    """Item 6: every GHCR ref pinned in containers.config and images.json is `<owner>/<name>-`."""
    name = tomllib.loads((REPO / "pyproject.toml").read_text())["project"]["name"]
    conf = REPO / "workflow" / "conf"
    text = (conf / "containers.config").read_text()
    if (conf / "images.json").is_file():
        text += (conf / "images.json").read_text()
    for ref in re.findall(r"ghcr\.io/[\w.-]+/[\w.-]+@sha256:[0-9a-f]{64}", text):
        assert ref.startswith(f"ghcr.io/khoa-yelo/{name}-"), ref


# --- the classification script, run against a fake docker -------------------------------------


def _fake_docker(tmp_path: Path, exit_code: int, output: str) -> tuple[dict[str, str], Path]:
    bindir = tmp_path / "bin"
    bindir.mkdir(exist_ok=True)
    log = tmp_path / "docker.log"
    (tmp_path / "out.txt").write_text(output)
    fake = bindir / "docker"
    fake.write_text(
        "#!/bin/sh\n"
        f'printf "%s\\n" "$*" >> "{log}"\n'
        f'cat "{tmp_path / "out.txt"}" >&2\n'
        f"exit {exit_code}\n"
    )
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    env = {"PATH": f"{bindir}{os.pathsep}{os.environ['PATH']}", "HOME": str(tmp_path)}
    return env, log


def _check(
    tmp_path: Path, exit_code: int, output: str
) -> tuple[subprocess.CompletedProcess[str], Path]:
    assert CHECK.is_file(), "scripts/release/ghcr_tag_state.sh is missing"
    env, log = _fake_docker(tmp_path, exit_code, output)
    proc = subprocess.run(
        ["bash", str(CHECK), REF], env=env, capture_output=True, text=True, timeout=30
    )
    return proc, log


def test_example5_inspect_exit_0_is_exists(tmp_path: Path) -> None:
    proc, log = _check(tmp_path, 0, "Name: " + REF)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip().splitlines()[-1] == "exists"
    assert log.read_text().split() == ["buildx", "imagetools", "inspect", REF]


@pytest.mark.parametrize(
    "message",
    [
        f"ERROR: {REF}: not found",
        "ERROR: failed to resolve: manifest unknown",
        'ERROR: unexpected status: 404 {"errors":[{"code":"NAME_UNKNOWN","message":"x"}]}',
    ],
)
def test_registry_404_is_absent(tmp_path: Path, message: str) -> None:
    proc, _ = _check(tmp_path, 1, message)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip().splitlines()[-1] == "absent"


@pytest.mark.parametrize(
    "message",
    [
        "ERROR: unexpected status from HEAD request: 401 Unauthorized",  # example 5b
        "ERROR: failed to authorize: unauthorized: authentication required",
        "ERROR: denied: permission_denied: The token provided does not match expected scopes.",
        "ERROR: dial tcp: lookup ghcr.io on 127.0.0.53:53: server misbehaving",
        "ERROR: toomanyrequests: retry later",
        "",
    ],
)
def test_example5b_any_other_error_fails(tmp_path: Path, message: str) -> None:
    proc, _ = _check(tmp_path, 1, message)
    assert proc.returncode != 0
    assert "absent" not in proc.stdout and "exists" not in proc.stdout
    if message:
        assert message in proc.stderr


@settings(
    max_examples=25,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)
@given(
    noise=st.text(alphabet=st.characters(min_codepoint=32, max_codepoint=126), max_size=40),
    pattern=st.sampled_from(ABSENT_PATTERNS),
    code=st.integers(min_value=1, max_value=125),
)
def test_property_classification(tmp_path: Path, noise: str, pattern: str, code: int) -> None:
    """Absent iff the error names a registry 404; any other failure text never reads as absent."""
    if any(w in noise for w in (*ABSENT_PATTERNS, "absent", "exists")):
        return
    proc, _ = _check(tmp_path, code, f"ERROR: {noise} {pattern}")
    assert (proc.returncode, proc.stdout.strip().splitlines()[-1:]) == (0, ["absent"])
    proc, _ = _check(tmp_path, code, f"ERROR: {noise}")
    assert proc.returncode != 0 and "absent" not in proc.stdout


# --- review round 1 (PR #383) -----------------------------------------------------------------


def test_build_and_digest_pr_run_only_on_main() -> None:
    """Review item 8: a dispatch from a feature branch never reaches a push or the digest PR."""
    _, build, _ = _build()
    _, digest = _digest_job()
    for job in (build, digest):
        assert re.search(r"github\.ref\s*==\s*'refs/heads/main'", str(job.get("if", ""))), job.get(
            "if"
        )


def test_plan_reads_the_revision_label_of_existing_images() -> None:
    """Review item 5: an existing image's provenance comes from its OCI revision label."""
    _, plan = _plan()
    assert "org.opencontainers.image.revision" in _runs(plan)


def _check_without_docker(tmp_path: Path) -> subprocess.CompletedProcess[str]:
    """Run the check with a PATH that has the usual tools but no `docker`."""
    assert CHECK.is_file(), "scripts/release/ghcr_tag_state.sh is missing"
    bindir = tmp_path / "nodocker"
    bindir.mkdir()
    for tool in ("grep", "cat", "sed", "awk", "head", "tail", "tr", "cut", "printf", "env"):
        found = shutil.which(tool)
        if found:
            (bindir / tool).symlink_to(found)
    bash = shutil.which("bash")
    assert bash
    return subprocess.run(
        [bash, str(CHECK), REF],
        env={"PATH": str(bindir), "HOME": str(tmp_path)},
        capture_output=True,
        text=True,
        timeout=30,
    )


def test_missing_docker_is_an_error_not_absent(tmp_path: Path) -> None:
    """Review item 1: bash's `docker: command not found` (rc 127) must not match `not found`."""
    proc = _check_without_docker(tmp_path)
    assert proc.returncode != 0, proc.stdout
    assert "absent" not in proc.stdout and "exists" not in proc.stdout


@pytest.mark.parametrize(
    ("code", "message"),
    [
        (127, "bash: docker: command not found"),
        (126, "bash: /usr/bin/docker: Permission denied"),
        (127, "/usr/libexec/docker/cli-plugins/docker-buildx: not found"),
    ],
)
def test_docker_cannot_run_is_an_error(tmp_path: Path, code: int, message: str) -> None:
    """rc 126/127 means docker (or buildx) could not run: never `absent`, whatever the text."""
    proc, _ = _check(tmp_path, code, message)
    assert proc.returncode != 0, proc.stdout
    assert "absent" not in proc.stdout
