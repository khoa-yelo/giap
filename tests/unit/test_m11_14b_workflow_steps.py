"""Run ``containers.yml``'s own ``run:`` steps against a fake ``docker`` (M11-14b, AC 2; PR #383).

Spec docs/specs/M11-14-container-install.md item 2 ("Existence check", "Re-check before push",
"Race"), examples 5, 5b, 5c; review of PR #383 items 2, 5 and 9. Offline: the plan job's
``run:`` steps, the build job's re-check and the digest-PR step are taken from the YAML and run
with ``bash -eo pipefail`` (as GitHub does) in a temporary copy of the tree. ``docker``, ``uv``,
``gh`` (and, for the digest-PR step, ``git``) are fakes on ``PATH``; ``${{ … }}`` in a step's
``env:`` is filled from a context the test gives.

Readings fixed here (where the spec and the review leave room):

* The plan job's ``build`` output is a JSON list of matrix entries with ``tool``, ``ref``,
  ``version`` (the tag without ``-<key>``), ``key`` and ``licence``; ``existing`` names the tools
  whose tag already exists and whose digest ``images.json`` does not pin yet. Its exact shape is
  not fixed: the tests feed it to the digest-PR step unchanged.
* On ``push``, the build list is every GHCR tool whose key tag is absent from GHCR, whatever this
  push changed (a pending run cancelled by the global concurrency group must not lose its tools;
  review item 2). ``workflow_dispatch`` with ``tools=<t>`` plans only ``t``.
* For a tag that exists, ``images.json`` ``built_from`` and the ``containers.config`` record line
  carry the image's ``org.opencontainers.image.revision`` label and its ``created`` date, not the
  current commit and today (review item 5). The fake answers
  ``imagetools inspect --format '{{json .Image}}'`` (an image config with ``created`` and
  ``config.Labels``), ``--format '{{json .Manifest}}'`` (with ``digest``) and ``--raw`` (bytes
  whose sha256 is that digest).
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import stat
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml

from gippa.annotate.tools import TOOLS
from gippa.io.image_key import key as image_key
from gippa.io.image_key import tag as image_tag

REPO = Path(__file__).resolve().parents[2]
WF = REPO / ".github" / "workflows" / "containers.yml"
GHCR = sorted(
    ["bactermfinder", "fgsea", "ostir", "promoter_calculator", "promotech", "terminatornet"]
)
NS = "ghcr.io/khoa-yelo"
SHA = "c0ffee" + "0" * 34
OLD_REVISION = "0123456789abcdef0123456789abcdef01234567"
CREATED = "2026-10-20T08:30:00Z"

FAKE_DOCKER = r"""#!{python}
import hashlib, json, os, sys
args = sys.argv[1:]
with open(os.environ["FAKE_DOCKER_LOG"], "a") as log:
    log.write(json.dumps(args) + "\n")
if args[:3] != ["buildx", "imagetools", "inspect"]:
    sys.stderr.write("fake docker: unsupported " + " ".join(args) + "\n")
    sys.exit(2)
rest = args[3:]
ref = [a for a in rest if not a.startswith("--") and "{{{{" not in a][-1]
fmt = next((rest[i + 1] for i, a in enumerate(rest) if a == "--format"), "")
fmt = fmt or next((a.split("=", 1)[1] for a in rest if a.startswith("--format=")), "")
state = json.load(open(os.environ["FAKE_DOCKER_STATE"]))
entry = state.get(ref)
if entry is None:
    sys.stderr.write(f"ERROR: {{ref}}: not found\n")
    sys.exit(1)
if "error" in entry:
    sys.stderr.write(entry["error"] + "\n")
    sys.exit(1)
raw = json.dumps({{"schemaVersion": 2, "ref": ref}}, separators=(",", ":")).encode()
digest = "sha256:" + hashlib.sha256(raw).hexdigest()
if "--raw" in rest:
    sys.stdout.buffer.write(raw)
elif ".Manifest" in fmt:
    print(json.dumps({{"schemaVersion": 2, "mediaType":
        "application/vnd.oci.image.manifest.v1+json", "digest": digest, "size": len(raw)}}))
elif ".Image" in fmt:
    print(json.dumps({{"created": entry["created"], "architecture": "amd64", "os": "linux",
        "config": {{"Labels": {{"org.opencontainers.image.revision": entry["revision"]}}}}}}))
else:
    print(f"Name:      {{ref}}\nMediaType: application/vnd.oci.image.manifest.v1+json\n"
          f"Digest:    {{digest}}")
"""

FAKE_UV = r"""#!/bin/bash
# fake uv: `uv sync …` does nothing; `uv run [--no-sync] [python] ARGS` runs the test's python
if [ "$1" = sync ]; then exit 0; fi
shift
args=()
for a in "$@"; do [ "$a" = --no-sync ] || [ "$a" = --frozen ] || args+=("$a"); done
if [ "${args[0]}" = python ] || [ "${args[0]}" = python3 ]; then args=("${args[@]:1}"); fi
exec "$FAKE_PYTHON" "${args[@]}"
"""


def _digest(ref: str) -> str:
    raw = json.dumps({"schemaVersion": 2, "ref": ref}, separators=(",", ":")).encode()
    return "sha256:" + hashlib.sha256(raw).hexdigest()


# --- the temporary tree and the fakes -------------------------------------------------------


def _git(root: Path, *args: str) -> str:
    env = {
        **os.environ,
        "GIT_AUTHOR_NAME": "t",
        "GIT_AUTHOR_EMAIL": "t@t",
        "GIT_COMMITTER_NAME": "t",
        "GIT_COMMITTER_EMAIL": "t@t",
        "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_SYSTEM": os.devnull,
    }
    out = subprocess.run(
        ["git", "-C", str(root), *args], env=env, capture_output=True, text=True, check=True
    )
    return out.stdout.strip()


def _tree(tmp_path: Path) -> Path:
    root = tmp_path / "repo"
    root.mkdir()
    shutil.copytree(
        REPO / "workflow",
        root / "workflow",
        ignore=shutil.ignore_patterns("tests", "__pycache__", ".nf-test*"),
    )
    shutil.copytree(
        REPO / "scripts" / "release",
        root / "scripts" / "release",
        ignore=shutil.ignore_patterns("__pycache__"),
    )
    shutil.copy(REPO / "pyproject.toml", root / "pyproject.toml")
    (root / "src").symlink_to(REPO / "src")
    (root / "workflow" / "conf" / "images.json").write_text("{}\n")
    _git(root, "init", "-q", "-b", "main")
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "c1")
    return root


def _commit(root: Path, rel: str, extra: str = "# touched\n") -> str:
    p = root / rel
    p.write_text(p.read_text() + extra)
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", f"touch {rel}")
    return _git(root, "rev-parse", "HEAD")


def _exe(path: Path, text: str) -> None:
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _fakes(tmp_path: Path, *, fake_git: bool = False) -> tuple[Path, dict[str, str]]:
    bindir = tmp_path / "bin"
    bindir.mkdir(exist_ok=True)
    _exe(bindir / "docker", FAKE_DOCKER.format(python=sys.executable))
    _exe(bindir / "uv", FAKE_UV)
    _exe(bindir / "python3", f'#!/bin/bash\nexec "{sys.executable}" "$@"\n')
    _exe(bindir / "gh", f'#!/bin/bash\necho "gh $*" >> "{tmp_path / "gh.log"}"\nexit 0\n')
    if fake_git:
        _exe(bindir / "git", f'#!/bin/bash\necho "git $*" >> "{tmp_path / "git.log"}"\nexit 0\n')
    state = tmp_path / "docker_state.json"
    state.write_text("{}")
    env = {
        "PATH": f"{bindir}{os.pathsep}{os.environ['PATH']}",
        "HOME": str(tmp_path),
        "FAKE_DOCKER_STATE": str(state),
        "FAKE_DOCKER_LOG": str(tmp_path / "docker.log"),
        "FAKE_PYTHON": sys.executable,
        "PYTHONPATH": str(REPO / "src"),
        "TMPDIR": str(tmp_path),
    }
    return state, env


def _ref(root: Path, tool: str) -> str:
    k = image_key(tool, root)
    return f"{NS}/gippa-{tool}:{image_tag(TOOLS[tool].version, k)}"


def _set_state(state: Path, entries: dict[str, dict[str, str]]) -> None:
    state.write_text(json.dumps(entries))


def _exists(revision: str = OLD_REVISION) -> dict[str, str]:
    return {"revision": revision, "created": CREATED}


# --- the step runner --------------------------------------------------------------------------


def _wf() -> dict[str, Any]:
    data: dict[str, Any] = yaml.safe_load(WF.read_text())
    return data


def _build_job() -> tuple[str, dict[str, Any], int]:
    for name, job in _wf()["jobs"].items():
        for i, s in enumerate(job.get("steps", [])):
            if str(s.get("uses", "")).startswith("docker/build-push-action@"):
                return name, job, i
    raise AssertionError("no docker/build-push-action step")


def _plan_job() -> tuple[str, dict[str, Any]]:
    _, job, _ = _build_job()
    matrix = yaml.safe_dump(job["strategy"]["matrix"], width=10_000)
    m = re.search(r"needs\.([\w-]+)\.outputs\.", matrix)
    assert m, matrix
    return m.group(1), _wf()["jobs"][m.group(1)]


def _digest_job() -> dict[str, Any]:
    for job in _wf()["jobs"].values():
        if any("update_digests.py" in str(s.get("run", "")) for s in job.get("steps", [])):
            return dict(job)
    raise AssertionError("no digest-PR job")


_EXPR = re.compile(r"\$\{\{\s*(.*?)\s*\}\}")


def _fill(value: Any, ctx: dict[str, str]) -> str:
    def lookup(m: re.Match[str]) -> str:
        expr = m.group(1)
        if expr in ctx:
            return ctx[expr]
        m2 = re.fullmatch(r"github\.event\.inputs\.(\w+)", expr)
        if m2:
            return ctx.get(f"inputs.{m2.group(1)}", "")
        return ""

    return _EXPR.sub(lookup, str(value))


def _parse_outputs(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    lines = path.read_text().splitlines() if path.is_file() else []
    i = 0
    while i < len(lines):
        line = lines[i]
        m = re.fullmatch(r"([\w-]+)<<(\S+)", line)
        if m:
            end = lines.index(m.group(2), i + 1)
            out[m.group(1)] = "\n".join(lines[i + 1 : end])
            i = end + 1
            continue
        if "=" in line:
            k, v = line.split("=", 1)
            out[k] = v
        i += 1
    return out


def _run_steps(
    job: dict[str, Any],
    steps: list[dict[str, Any]],
    root: Path,
    env: dict[str, str],
    ctx: dict[str, str],
) -> tuple[subprocess.CompletedProcess[str], dict[str, str]]:
    """Run the ``run:`` steps in order; returns the last process and the job's outputs."""
    ctx = dict(ctx)
    bash = shutil.which("bash")
    assert bash
    proc: subprocess.CompletedProcess[str] | None = None
    base = {**env, **{k: _fill(v, ctx) for k, v in (_wf().get("env") or {}).items()}}
    base.update({k: _fill(v, ctx) for k, v in (job.get("env") or {}).items()})
    for n, step in enumerate(steps):
        if "run" not in step:
            continue
        out = root.parent / f"out-{n}.txt"
        out.write_text("")
        script = root.parent / f"step-{n}.sh"
        script.write_text(str(step["run"]))
        step_env = {
            **base,
            **{k: _fill(v, ctx) for k, v in (step.get("env") or {}).items()},
            "GITHUB_OUTPUT": str(out),
            "GITHUB_STEP_SUMMARY": str(root.parent / "summary.md"),
            "GITHUB_SHA": ctx["github.sha"],
            "GITHUB_WORKSPACE": str(root),
        }
        proc = subprocess.run(
            [bash, "--noprofile", "--norc", "-eo", "pipefail", str(script)],
            cwd=root,
            env=step_env,
            capture_output=True,
            text=True,
            timeout=120,
        )
        if proc.returncode != 0:
            return proc, {}
        if "id" in step:
            for k, v in _parse_outputs(out).items():
                ctx[f"steps.{step['id']}.outputs.{k}"] = v
    assert proc is not None, "no run: step"
    outputs = {k: _fill(v, ctx) for k, v in (job.get("outputs") or {}).items()}
    return proc, outputs


def _ctx(event: str, before: str, after: str, **extra: str) -> dict[str, str]:
    return {
        "github.event_name": event,
        "github.event.before": before,
        "github.sha": after,
        "github.ref": "refs/heads/main",
        "vars.GHCR_PUBLISH": "true",
        **extra,
    }


def _plan(
    root: Path, env: dict[str, str], ctx: dict[str, str]
) -> tuple[subprocess.CompletedProcess[str], dict[str, str]]:
    _, job = _plan_job()
    return _run_steps(job, list(job["steps"]), root, env, ctx)


def _build_tools(outputs: dict[str, str]) -> list[str]:
    return sorted(e["tool"] for e in json.loads(outputs["build"]))


def _existing_tools(outputs: dict[str, str]) -> list[str]:
    raw = json.loads(outputs.get("existing") or "{}")
    return sorted(raw if isinstance(raw, dict) else (e["tool"] for e in raw))


def _pin(root: Path, tools: dict[str, str]) -> None:
    """images.json pins ``tool -> ref (by digest)`` with the current key."""
    images = {}
    for tool, ref in tools.items():
        k = image_key(tool, root)
        images[tool] = {
            "key": k,
            "tag": image_tag(TOOLS[tool].version, k),
            "ref": ref,
            "processes": [],
            "built_from": OLD_REVISION,
            "licence": TOOLS[tool].licence,
        }
    (root / "workflow" / "conf" / "images.json").write_text(json.dumps(images, indent=2) + "\n")
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "pin")


# --- plan job: what is built (review items 2, 9 M9/M16) ------------------------------------------


def test_plan_matrix_entries_carry_ref_version_key_and_licence(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    _, env = _fakes(tmp_path)
    after = _commit(root, "workflow/envs/ostir.yml")
    proc, out = _plan(root, env, _ctx("workflow_dispatch", "", after, **{"inputs.tools": "all"}))
    assert proc.returncode == 0, proc.stderr
    entries = {e["tool"]: e for e in json.loads(out["build"])}
    assert sorted(entries) == GHCR
    for tool, e in entries.items():
        k = image_key(tool, root)
        t = image_tag(TOOLS[tool].version, k)
        assert e["ref"] == f"{NS}/gippa-{tool}:{t}"
        assert e["key"] == k
        assert e["version"] == t.removesuffix(f"-{k}"), "version is the tag without -<key>"
        assert e["licence"] == TOOLS[tool].licence
    assert entries["ostir"]["version"] == "1.1.3"
    assert entries["promotech"]["version"] == "git-56251ad9"


def test_dispatch_of_one_tool_plans_only_that_tool(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    _, env = _fakes(tmp_path)
    head = _git(root, "rev-parse", "HEAD")
    proc, out = _plan(root, env, _ctx("workflow_dispatch", "", head, **{"inputs.tools": "ostir"}))
    assert proc.returncode == 0, proc.stderr
    assert _build_tools(out) == ["ostir"]


def test_push_plans_tools_of_a_cancelled_pending_run(tmp_path: Path) -> None:
    """Review item 2: push B (promotech) was cancelled while pending; push C (terminatornet) runs.

    promotech's new key has no tag, so C's plan must build it although C did not touch it."""
    root = _tree(tmp_path)
    state, env = _fakes(tmp_path)
    old = {t: _ref(root, t) for t in GHCR}
    _pin(root, {t: f"{old[t].rsplit(':', 1)[0]}@{_digest(old[t])}" for t in GHCR})
    b = _commit(root, "workflow/envs/promotech.yml")
    c = _commit(root, "workflow/envs/terminatornet.yml")
    _set_state(state, {old[t]: _exists() for t in GHCR})  # the c1 tags exist
    proc, out = _plan(root, env, _ctx("push", b, c))
    assert proc.returncode == 0, proc.stderr
    assert _build_tools(out) == ["promotech", "terminatornet"]
    assert _existing_tools(out) == [], "unchanged tools are already pinned to that digest"


def test_push_builds_every_unpublished_tool(tmp_path: Path) -> None:
    """Nothing published yet: a push touching one env file plans all six (no tool is lost)."""
    root = _tree(tmp_path)
    _, env = _fakes(tmp_path)
    before = _git(root, "rev-parse", "HEAD")
    after = _commit(root, "workflow/envs/ostir.yml")
    proc, out = _plan(root, env, _ctx("push", before, after))
    assert proc.returncode == 0, proc.stderr
    assert _build_tools(out) == GHCR


def test_push_changing_only_a_lock_builds_that_tool(tmp_path: Path) -> None:
    """Review item 9 (M16): a lock file is a key input, so its change rebuilds the tool."""
    root = _tree(tmp_path)
    state, env = _fakes(tmp_path)
    old = {t: _ref(root, t) for t in GHCR}
    _pin(root, {t: f"{old[t].rsplit(':', 1)[0]}@{_digest(old[t])}" for t in GHCR})
    before = _git(root, "rev-parse", "HEAD")
    after = _commit(root, "workflow/envs/locks/terminatornet.pip.txt", "\n")
    _set_state(state, {old[t]: _exists() for t in GHCR})
    proc, out = _plan(root, env, _ctx("push", before, after))
    assert proc.returncode == 0, proc.stderr
    assert _build_tools(out) == ["terminatornet"]


# --- plan job: existing tags and re-pinning (review items 5, 9 M11/M12) --------------------------


def test_existing_tag_is_not_built_and_is_reported(tmp_path: Path) -> None:
    """Example 5: inspect exits 0 for ostir's tag: ostir is not built; its digest is reported."""
    root = _tree(tmp_path)
    state, env = _fakes(tmp_path)
    head = _git(root, "rev-parse", "HEAD")
    _set_state(state, {_ref(root, "ostir"): _exists()})
    proc, out = _plan(root, env, _ctx("workflow_dispatch", "", head, **{"inputs.tools": "all"}))
    assert proc.returncode == 0, proc.stderr
    assert "ostir" not in _build_tools(out)
    assert _existing_tools(out) == ["ostir"]
    assert _digest(_ref(root, "ostir")) in out["existing"]


def test_existing_tag_already_pinned_to_that_digest_is_not_repinned(tmp_path: Path) -> None:
    """Review item 9 (M11): re-pin only when images.json names another digest."""
    root = _tree(tmp_path)
    state, env = _fakes(tmp_path)
    o, p = _ref(root, "ostir"), _ref(root, "promotech")
    _pin(
        root,
        {
            "ostir": f"{o.rsplit(':', 1)[0]}@{_digest(o)}",  # the same digest: nothing to do
            "promotech": f"{p.rsplit(':', 1)[0]}@sha256:{'e' * 64}",  # another digest: re-pin
        },
    )
    head = _git(root, "rev-parse", "HEAD")
    _set_state(state, {o: _exists(), p: _exists()})
    proc, out = _plan(root, env, _ctx("workflow_dispatch", "", head, **{"inputs.tools": "all"}))
    assert proc.returncode == 0, proc.stderr
    assert _existing_tools(out) == ["promotech"]
    assert sorted(set(GHCR) - {"ostir", "promotech"}) == _build_tools(out)


def test_example5b_unauthorized_fails_the_plan(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    state, env = _fakes(tmp_path)
    head = _git(root, "rev-parse", "HEAD")
    _set_state(state, {_ref(root, "ostir"): {"error": "ERROR: 401 Unauthorized"}})
    proc, _ = _plan(root, env, _ctx("workflow_dispatch", "", head, **{"inputs.tools": "all"}))
    assert proc.returncode != 0


def test_publishing_off_is_printed(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    _, env = _fakes(tmp_path)
    head = _git(root, "rev-parse", "HEAD")
    ctx = _ctx("workflow_dispatch", "", head, **{"inputs.tools": "ostir", "vars.GHCR_PUBLISH": ""})
    proc, _ = _plan(root, env, ctx)
    assert proc.returncode == 0, proc.stderr
    assert "publishing off" in proc.stdout


# --- build job: the re-check (example 5c; review item 9 M5/M15) ----------------------------------


def _recheck(tmp_path: Path, exists: bool) -> subprocess.CompletedProcess[str]:
    root = _tree(tmp_path)
    state, env = _fakes(tmp_path)
    ref = _ref(root, "ostir")
    _set_state(state, {ref: _exists()} if exists else {})
    _, job, i = _build_job()
    step = job["steps"][i - 1]
    k = image_key("ostir", root)
    ctx = _ctx("push", "", SHA)
    ctx.update(
        {
            "matrix.tool": "ostir",
            "matrix.ref": ref,
            "matrix.key": k,
            "matrix.version": "1.1.3",
            "matrix.licence": TOOLS["ostir"].licence,
        }
    )
    proc, _ = _run_steps(job, [step], root, env, ctx)
    return proc


def test_example5c_recheck_fails_when_the_tag_appeared(tmp_path: Path) -> None:
    proc = _recheck(tmp_path, exists=True)
    assert proc.returncode != 0
    assert "tag exists: " in proc.stdout + proc.stderr


def test_recheck_passes_when_the_tag_is_still_absent(tmp_path: Path) -> None:
    proc = _recheck(tmp_path, exists=False)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    calls = [json.loads(x) for x in (tmp_path / "docker.log").read_text().splitlines()]
    assert any(c[:3] == ["buildx", "imagetools", "inspect"] for c in calls)


# --- digest PR: provenance of existing images (review item 5) ------------------------------------


def test_digest_pr_records_the_existing_images_own_revision_and_date(tmp_path: Path) -> None:
    """promotech's tag exists (built from OLD_REVISION on 2026-10-20); ostir is built now."""
    root = _tree(tmp_path)
    state, env = _fakes(tmp_path)
    head = _git(root, "rev-parse", "HEAD")
    p, o = _ref(root, "promotech"), _ref(root, "ostir")
    _set_state(state, {p: _exists(OLD_REVISION)})
    ctx = _ctx("workflow_dispatch", "", SHA, **{"inputs.tools": "promotech ostir"})
    proc, out = _plan(root, env, ctx)
    assert proc.returncode == 0, proc.stderr
    assert _build_tools(out) == ["ostir"] and _existing_tools(out) == ["promotech"]
    del head

    # the build leg of ostir pushed and uploaded its digest artifact
    digests = root / "digests"
    digests.mkdir()
    (digests / "ostir").write_text(f"{o.rsplit(':', 1)[0]}@sha256:{'a' * 64}\n")

    plan_name, _ = _plan_job()
    job = _digest_job()
    steps = [s for s in job["steps"] if "update_digests.py" in str(s.get("run", ""))]
    for k, v in out.items():
        ctx[f"needs.{plan_name}.outputs.{k}"] = v
    ctx["needs.build.result"] = "success"
    _, genv = _fakes(tmp_path, fake_git=True)
    proc, _ = _run_steps(job, steps, root, genv, ctx)
    assert proc.returncode == 0, proc.stdout + proc.stderr

    images = json.loads((root / "workflow" / "conf" / "images.json").read_text())
    assert images["ostir"]["built_from"] == SHA
    assert images["promotech"]["built_from"] == OLD_REVISION, "the image's revision label"
    assert images["promotech"]["ref"] == f"{p.rsplit(':', 1)[0]}@{_digest(p)}"
    config = (root / "workflow" / "conf" / "containers.config").read_text()
    line = [ln for ln in config.splitlines() if "promotech" in ln and OLD_REVISION in ln]
    assert line, "containers.config records promotech's own revision"
    assert "2026-10-20" in line[0], line[0]


@pytest.mark.parametrize("event", ["push", "workflow_dispatch"])
def test_plan_steps_never_call_docker_push_or_build(tmp_path: Path, event: str) -> None:
    root = _tree(tmp_path)
    _, env = _fakes(tmp_path)
    head = _git(root, "rev-parse", "HEAD")
    after = _commit(root, "workflow/envs/ostir.yml")
    ctx = _ctx(event, head, after, **{"inputs.tools": "all"})
    proc, _ = _plan(root, env, ctx)
    assert proc.returncode == 0, proc.stderr
    calls = [json.loads(x) for x in (tmp_path / "docker.log").read_text().splitlines()]
    assert calls and all(c[:3] == ["buildx", "imagetools", "inspect"] for c in calls), calls
