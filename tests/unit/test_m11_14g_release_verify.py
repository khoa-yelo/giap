"""``release.yml`` verify gains ``image_key.py --check`` and an anonymous ``imagetools inspect`` of
every GHCR pin (M11-14 14g; item 2 "The anonymous check lives only in release.yml verify", item 5
flow, edge case "Private package"; ADR-0060 §5). A CI change: Khoa approves the PR.

Offline: the workflow is parsed statically, and the inspect step's own shell is extracted and run
in a temporary tree with a fake ``docker`` on ``PATH`` (never a registry).

Interface fixed by these tests:

* The verify job runs ``scripts/release/image_key.py --check`` in a step with no ``if:`` (every
  dispatch, dry run or not) and never masks its exit (``||`` on that line is not allowed).
* One verify step (its ``run`` contains ``imagetools inspect``) collects every distinct GHCR ref
  (``ghcr.io/...@sha256:<64 hex>``) from ``workflow/conf/containers.config`` and
  ``workflow/conf/images.json`` (``ref``; fgsea is there only) of the working directory and runs
  ``docker buildx imagetools inspect <ref>`` for each. Any failure exits non-zero and prints
  ``not public: <ref>``. No GHCR pin: exit 0 without calling ``docker``. The step may call
  ``python3``/``python`` and ``uv run [--no-sync] python``.
* Anonymous: the verify job has no ``docker/login-action``, no ``docker login``, reads no
  ``secrets.`` and no ``GITHUB_TOKEN``/``github.token``, and asks for no ``packages`` permission.
"""

from __future__ import annotations

import json
import os
import stat
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
RELEASE_YML = REPO / ".github" / "workflows" / "release.yml"
OWNER = "ghcr.io/khoa-yelo/gippa"
OSTIR = f"{OWNER}-ostir@sha256:{'a' * 64}"
FGSEA = f"{OWNER}-fgsea@sha256:{'b' * 64}"
BAKTA = f"oschwengers/bakta@sha256:{'c' * 64}"


@pytest.fixture(scope="module")
def verify() -> dict[str, Any]:
    wf = yaml.safe_load(RELEASE_YML.read_text(encoding="utf-8"))
    job: dict[str, Any] = wf["jobs"]["verify"]
    return job


def _steps(job: dict[str, Any]) -> list[dict[str, Any]]:
    return list(job.get("steps") or [])


def _inspect_step(job: dict[str, Any]) -> dict[str, Any]:
    hits = [s for s in _steps(job) if "imagetools inspect" in str(s.get("run", ""))]
    assert len(hits) == 1, f"verify: {len(hits)} steps run `imagetools inspect` (want 1)"
    return hits[0]


# ------------------------------------------------------------------ static


def test_verify_runs_image_key_check_unmasked(verify: dict[str, Any]) -> None:
    hits = [
        s for s in _steps(verify) if "scripts/release/image_key.py --check" in str(s.get("run", ""))
    ]
    assert hits, "verify does not run scripts/release/image_key.py --check"
    for s in hits:
        assert "if" not in s, s
        for line in str(s["run"]).splitlines():
            if "image_key.py --check" in line:
                assert "||" not in line, line


def test_image_key_check_runs_after_checkout(verify: dict[str, Any]) -> None:
    steps = _steps(verify)
    co = next(
        i for i, s in enumerate(steps) if str(s.get("uses", "")).startswith("actions/checkout")
    )
    chk = next(i for i, s in enumerate(steps) if "image_key.py --check" in str(s.get("run", "")))
    assert co < chk


def test_inspect_step_runs_on_every_dispatch(verify: dict[str, Any]) -> None:
    step = _inspect_step(verify)
    assert "if" not in step, step
    assert "not public:" in str(step["run"])


def test_verify_is_anonymous(verify: dict[str, Any]) -> None:
    text = yaml.safe_dump(verify, sort_keys=False, width=10_000)
    assert "docker/login-action" not in text
    assert "docker login" not in text
    assert "secrets." not in text
    assert "GITHUB_TOKEN" not in text and "github.token" not in text
    assert "packages" not in (verify.get("permissions") or {}), verify.get("permissions")


# ------------------------------------------------------------------ the step's own shell


FAKE_DOCKER = """#!/usr/bin/env bash
printf '%s\\n' "$*" >> "$DOCKER_LOG"
for a in "$@"; do
  case " $FAIL_REFS " in *" $a "*)
    echo "ERROR: unauthorized: authentication required" >&2; exit 1;;
  esac
done
echo '{"mediaType": "application/vnd.oci.image.index.v1+json"}'
"""
FAKE_UV = """#!/usr/bin/env bash
[ "$1" = run ] && shift
while [ "${1#--}" != "$1" ]; do shift; done
exec "$@"
"""


def _exe(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _tree(tmp: Path, config_refs: dict[str, str], images: dict[str, Any]) -> Path:
    root = tmp / "tree"
    conf = root / "workflow" / "conf"
    conf.mkdir(parents=True)
    lines = [f"  withName: '{p}' {{ container = '{r}' }}" for p, r in config_refs.items()]
    (conf / "containers.config").write_text(
        "// pins (fixture)\nprocess {\n" + "\n".join(lines) + "\n}\n", encoding="utf-8"
    )
    (conf / "images.json").write_text(json.dumps(images, indent=2) + "\n", encoding="utf-8")
    for d in ("src", "scripts"):
        (root / d).symlink_to(REPO / d)
    return root


def _entry(ref: str, processes: list[str]) -> dict[str, Any]:
    return {"key": "0" * 12, "tag": "1.0-" + "0" * 12, "ref": ref, "processes": processes,
            "built_from": "0" * 40, "licence": "MIT"}  # fmt: skip


def _run(
    verify: dict[str, Any], tmp: Path, root: Path, fail: tuple[str, ...] = ()
) -> tuple[subprocess.CompletedProcess[str], list[str]]:
    bin_ = tmp / "bin"
    bin_.mkdir()
    _exe(bin_ / "docker", FAKE_DOCKER)
    _exe(bin_ / "uv", FAKE_UV)
    for name in ("python", "python3"):
        (bin_ / name).symlink_to(sys.executable)
    script = tmp / "step.sh"
    script.write_text(str(_inspect_step(verify)["run"]), encoding="utf-8")
    log, home = tmp / "docker.log", tmp / "home"
    home.mkdir()
    env = {
        **os.environ, "PATH": f"{bin_}{os.pathsep}{os.environ['PATH']}", "HOME": str(home),
        "DOCKER_LOG": str(log), "FAIL_REFS": " ".join(fail), "VERSION": "1.0.0",
        "DRY_RUN": "true", "TARGET": "testpypi", "GITHUB_OUTPUT": str(tmp / "out"),
        "GITHUB_STEP_SUMMARY": str(tmp / "summary"), "GITHUB_REF": "refs/heads/main",
        "PYTHONPATH": str(REPO / "src"),
    }  # fmt: skip
    p = subprocess.run(["bash", "-eo", "pipefail", str(script)], cwd=root, env=env,
                       capture_output=True, text=True, timeout=120, check=False)  # fmt: skip
    calls = log.read_text(encoding="utf-8").splitlines() if log.is_file() else []
    return p, calls


def _inspected(calls: list[str]) -> set[str]:
    return {c.split()[-1] for c in calls if "imagetools inspect" in c}


def test_every_ghcr_pin_is_inspected_and_public_passes(
    verify: dict[str, Any], tmp_path: Path
) -> None:
    images = {"ostir": _entry(OSTIR, ["OSTIR", "OSTIR_NOVEL"]), "fgsea": _entry(FGSEA, [])}
    root = _tree(tmp_path, {"OSTIR": OSTIR, "OSTIR_NOVEL": OSTIR, "BAKTA": BAKTA}, images)
    p, calls = _run(verify, tmp_path, root)
    assert p.returncode == 0, p.stdout + p.stderr
    assert {OSTIR, FGSEA} <= _inspected(calls), calls
    assert all(c.startswith("buildx imagetools inspect") for c in calls if "inspect" in c), calls
    assert BAKTA not in _inspected(calls), "only GHCR pins are checked"


@pytest.mark.parametrize("bad", [OSTIR, FGSEA], ids=["withName-pin", "images-json-only"])
def test_a_private_pin_fails_not_public(verify: dict[str, Any], tmp_path: Path, bad: str) -> None:
    root = _tree(tmp_path, {"OSTIR": OSTIR},
                 {"ostir": _entry(OSTIR, ["OSTIR"]), "fgsea": _entry(FGSEA, [])})  # fmt: skip
    p, _ = _run(verify, tmp_path, root, fail=(bad,))
    assert p.returncode != 0, p.stdout + p.stderr
    assert f"not public: {bad}" in p.stdout + p.stderr, p.stdout + p.stderr


def test_no_ghcr_pin_calls_no_docker(verify: dict[str, Any], tmp_path: Path) -> None:
    root = _tree(tmp_path, {"BAKTA": BAKTA}, {})
    p, calls = _run(verify, tmp_path, root)
    assert p.returncode == 0, p.stdout + p.stderr
    assert _inspected(calls) == set(), calls
