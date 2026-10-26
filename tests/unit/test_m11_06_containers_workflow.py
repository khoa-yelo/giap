"""Static checks of ``.github/workflows/containers.yml`` (M11-06 item 2, examples 6-7, AC 2).

Nothing is built, pushed or dispatched: the workflow file is parsed as YAML and its jobs, triggers,
permissions and steps are inspected. Readings chosen where the spec leaves room:

* the build job is the one job that uses ``docker/build-push-action``; its matrix may be a literal
  list or come from an earlier "plan" job's output (``fromJSON(needs.<job>.outputs.<x>)``);
* the owner in image names is ``khoa-yelo`` literally or ``${{ github.repository_owner }}``;
* the dispatch input ``tools`` defaults to ``all`` or to the five names (comma or space separated).
"""

from __future__ import annotations

import functools
import re
from pathlib import Path
from typing import Any

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
WF = REPO / ".github" / "workflows" / "containers.yml"
TOOLS = {"promotech", "promoter_calculator", "terminatornet", "bactermfinder", "ostir"}
SHA_PIN = re.compile(r"^[\w.-]+/[\w./-]+@[0-9a-f]{40}$")
OWNER = r"(?:khoa-yelo|\$\{\{\s*github\.repository_owner\s*\}\})"


@functools.cache
def _text() -> str:
    assert WF.is_file(), ".github/workflows/containers.yml is missing (M11-06 item 2)"
    return WF.read_text()


@functools.cache
def _load() -> dict[str, Any]:
    data: dict[str, Any] = yaml.safe_load(_text())
    return data


def _on(wf: dict[str, Any]) -> dict[str, Any]:
    raw: dict[Any, Any] = dict(wf)
    on = raw.get("on", raw.get(True))  # PyYAML reads the bare key `on` as True
    assert isinstance(on, dict), "`on:` must be a mapping"
    return on


def _jobs(wf: dict[str, Any]) -> dict[str, Any]:
    jobs: dict[str, Any] = wf["jobs"]
    return jobs


def _steps(job: dict[str, Any]) -> list[dict[str, Any]]:
    return list(job.get("steps", []))


def _runs(job: dict[str, Any]) -> str:
    return "\n".join(str(s.get("run", "")) for s in _steps(job))


def _needs(job: dict[str, Any]) -> list[str]:
    n = job.get("needs", [])
    return [n] if isinstance(n, str) else list(n)


def _uses(wf: dict[str, Any]) -> list[str]:
    out = []
    for job in _jobs(wf).values():
        if "uses" in job:  # reusable workflow call
            out.append(str(job["uses"]))
        for s in _steps(job):
            if "uses" in s:
                out.append(str(s["uses"]))
    return out


def _build(wf: dict[str, Any]) -> tuple[str, dict[str, Any], dict[str, Any]]:
    hits = [
        (name, job, s)
        for name, job in _jobs(wf).items()
        for s in _steps(job)
        if str(s.get("uses", "")).startswith("docker/build-push-action@")
    ]
    assert len(hits) == 1, "exactly one docker/build-push-action step expected"
    return hits[0]


def _perms(block: Any) -> dict[str, str]:
    if block is None:
        return {}
    if isinstance(block, str):  # read-all / write-all
        return {"*": block}
    return {str(k): str(v) for k, v in block.items()}


# --- triggers -------------------------------------------------------------------------------


def test_push_to_main_on_container_and_env_paths() -> None:
    wf = _load()
    push = _on(wf)["push"]
    assert push["branches"] == ["main"]
    paths = set(push["paths"])
    assert {"workflow/containers/**", "workflow/envs/**"} <= paths
    assert "paths-ignore" not in push
    assert "tags" not in push


def test_manual_dispatch_with_tools_input_defaulting_to_all_five() -> None:
    wf = _load()
    inputs = _on(wf)["workflow_dispatch"]["inputs"]
    default = str(inputs["tools"]["default"]).strip()
    assert default == "all" or set(re.split(r"[,\s]+", default)) == TOOLS, default


def test_no_pull_request_target_or_other_triggers() -> None:
    wf = _load()
    on = _on(wf)
    assert "pull_request_target" not in on
    assert set(on) <= {"push", "workflow_dispatch"}, sorted(on)


# --- permissions ----------------------------------------------------------------------------


def test_top_level_permissions_are_read_only() -> None:
    wf = _load()
    assert "permissions" in wf, "set a read-only top-level `permissions:` block"
    top = _perms(wf["permissions"])
    assert top.get("*") != "write-all"
    assert "write" not in top.values(), top


def test_packages_write_only_on_the_build_job() -> None:
    wf = _load()
    name, job, _ = _build(wf)
    assert _perms(job.get("permissions")).get("packages") == "write"
    for other, j in _jobs(wf).items():
        if other != name:
            p = _perms(j.get("permissions"))
            assert p.get("packages") != "write", other
            assert p.get("*") != "write-all", other


# --- the build job --------------------------------------------------------------------------


def test_build_job_is_a_matrix() -> None:
    wf = _load()
    _, job, _ = _build(wf)
    assert "matrix" in job.get("strategy", {}), "the build job needs strategy.matrix"


def test_literal_matrix_lists_exactly_the_five_tools() -> None:
    wf = _load()
    _, job, _ = _build(wf)
    matrix = job["strategy"]["matrix"]
    dynamic = (
        [matrix]
        if isinstance(matrix, str)
        else [v for v in matrix.values() if isinstance(v, str) and "fromJSON" in v]
    )
    if dynamic:  # fromJSON(needs.plan.outputs.x): the plan job is a dependency
        for expr in dynamic:
            m = re.search(r"needs\.([\w-]+)\.outputs\.", expr)
            assert m, expr
            assert m.group(1) in _needs(job)
        if isinstance(matrix, str):
            return
    lists = [v for v in matrix.values() if isinstance(v, list)]
    assert lists or dynamic, matrix
    for v in lists:
        if all(isinstance(x, str) for x in v) and set(v) & (TOOLS | {"uniop"}):
            assert set(v) == TOOLS, v
    assert "uniop" not in yaml.safe_dump(matrix)
    assert "fgsea" not in yaml.safe_dump(matrix)


def test_build_push_action_inputs() -> None:
    wf = _load()
    _, _, step = _build(wf)
    w = step["with"]
    assert str(w["context"]).strip() in {".", "${{ github.workspace }}"}
    assert re.fullmatch(
        r"workflow/containers/\$\{\{\s*matrix\.[\w.]+\s*\}\}/Dockerfile", str(w["file"]).strip()
    ), w["file"]
    assert str(w["platforms"]).strip() == "linux/amd64"
    assert str(w["push"]).strip().lower() == "true"


# test_two_tags_sha7_and_tool_version retired by M11-14 14b: one content-key tag per image
# (tests/unit/test_m11_14b_containers_workflow.py pins the single tag).


def test_oci_labels() -> None:
    wf = _load()
    _, _, step = _build(wf)
    labels = str(step["with"].get("labels", ""))
    for key in ("source", "revision", "licenses"):
        assert f"org.opencontainers.image.{key}=" in labels, key


def test_build_job_logs_in_to_ghcr() -> None:
    wf = _load()
    _, job, _ = _build(wf)
    logins = [s for s in _steps(job) if str(s.get("uses", "")).startswith("docker/login-action@")]
    assert logins, "no docker/login-action step in the build job"
    assert str(logins[0]["with"]["registry"]).strip() == "ghcr.io"


def test_build_job_exposes_the_pushed_digest() -> None:
    wf = _load()
    _, job, step = _build(wf)
    assert "id" in step, "give the build-push step an id to read its digest output"
    assert f"steps.{step['id']}.outputs.digest" in yaml.safe_dump(job)


# --- guards: no UniOP, no fgsea, tools input validated ---------------------------------------


def test_uniop_rejected_before_any_build() -> None:
    """Example 7: dispatch with tools=uniop fails in a job the build job needs."""
    wf = _load()
    name, job, _ = _build(wf)
    jobs = _jobs(wf)
    guards = [n for n, j in jobs.items() if n != name and "uniop is never built here" in _runs(j)]
    assert guards, "no step prints `uniop is never built here`"
    assert set(guards) & set(_needs(job)), "the uniop guard must run in a job the build needs"
    assert any("exit 1" in _runs(jobs[g]) for g in guards)


def test_uniop_never_built_or_pushed() -> None:
    text = _text()
    assert "workflow/containers/uniop" not in text
    assert not re.search(r"giap-uniop", text)


def test_dispatch_input_is_validated_against_the_five() -> None:
    """Every tool name appears in the plan/guard logic, so an unknown name can be rejected."""
    wf = _load()
    name, _, _ = _build(wf)
    plan = "\n".join(_runs(j) for n, j in _jobs(wf).items() if n != name)
    for t in TOOLS:
        assert t in plan, t
    assert "inputs.tools" in plan or "INPUT_TOOLS" in plan or "github.event.inputs.tools" in plan


def test_push_maps_changed_env_or_dockerfile_to_its_tool() -> None:
    """Example 6: a push touching workflow/envs/ostir.yml builds ostir only (paths → tools)."""
    wf = _load()
    name, _, _ = _build(wf)
    plan = "\n".join(_runs(j) for n, j in _jobs(wf).items() if n != name)
    uses_filter = any(u.startswith("dorny/paths-filter@") for u in _uses(wf))
    assert uses_filter or re.search(r"git diff[^\n]*--name-only", plan), (
        "compute the changed files (git diff --name-only or a SHA-pinned paths filter)"
    )
    if not uses_filter:
        assert "envs/" in plan and "containers/" in plan


# --- the digest PR job ---------------------------------------------------------------------


def _digest_job(wf: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    hits = [(n, j) for n, j in _jobs(wf).items() if "scripts/release/update_digests.py" in _runs(j)]
    assert len(hits) == 1, "exactly one job runs scripts/release/update_digests.py"
    return hits[0]


def test_digest_job_runs_after_build() -> None:
    wf = _load()
    build_name, _, _ = _build(wf)
    _, job = _digest_job(wf)
    assert build_name in _needs(job)
    assert "--digests" in _runs(job)


def test_digest_job_opens_a_pr_on_a_ci_branch() -> None:
    wf = _load()
    _, job = _digest_job(wf)
    blob = yaml.safe_dump(job, width=10_000)
    assert "ci/digests-" in blob
    assert "ci: pin GHCR digests (" in blob
    cpr = [
        s
        for s in _steps(job)
        if str(s.get("uses", "")).startswith("peter-evans/create-pull-request@")
    ]
    assert cpr or "gh pr create" in _runs(job)
    for s in cpr:
        assert str(s["with"]["branch"]).startswith("ci/digests-")
        assert str(s.get("with", {}).get("base", "main")) == "main"
    perms = _perms(job.get("permissions"))
    assert perms.get("pull-requests") == "write"
    assert perms.get("contents") == "write"


def test_never_pushes_to_main() -> None:
    wf = _load()
    for name, job in _jobs(wf).items():
        runs = _runs(job)
        assert not re.search(r"git push[^\n]*\b(main|HEAD:main|refs/heads/main)\b", runs), name
        assert "--force" not in runs and "push -f" not in runs, name
        for s in _steps(job):
            br = s.get("with", {}).get("branch") if isinstance(s.get("with"), dict) else None
            assert br != "main", name


# --- action pinning (Q-L5) ------------------------------------------------------------------


def test_every_third_party_action_is_sha_pinned() -> None:
    wf = _load()
    for u in _uses(wf):
        if u.startswith(("actions/", "./")):
            continue
        assert SHA_PIN.match(u), f"{u} is not pinned to a 40-hex commit SHA (Q-L5)"


def test_sha_pins_carry_a_tag_comment() -> None:
    text = _text()
    for line in text.splitlines():
        m = re.search(r"uses:\s*([\w.-]+/[\w./-]+)@([0-9a-f]{40})", line)
        if m and not m.group(1).startswith("actions/"):
            assert "#" in line.split("@", 1)[1], f"tag comment missing: {line.strip()}"


# --- the five tools have what the workflow builds from --------------------------------------


@pytest.mark.parametrize("tool", sorted(TOOLS))
def test_each_tool_has_dockerfile_and_env(tool: str) -> None:
    assert (REPO / "workflow" / "containers" / tool / "Dockerfile").is_file()
    assert (REPO / "workflow" / "envs" / f"{tool}.yml").is_file()
