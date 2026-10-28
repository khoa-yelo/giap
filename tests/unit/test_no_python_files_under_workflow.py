"""Guard: no git-tracked ``.py``/``.pyi`` file may ship under ``workflow/``.

pyproject's ``[tool.hatch.build.targets.wheel.sources]`` maps ``"workflow"`` to
``"gippa/workflow"``. That rewrite *changes* a prefix, and hatchling's editable
(dev-mode) build refuses it as soon as a selected ``.py`` file lives under that
prefix::

    Dev mode installations are unsupported when any path rewrite in the
    `sources` option changes a prefix rather than removes it

So ``uv sync`` (an editable install) fails everywhere, CI included. PR #455 hit
this with ``workflow/bin/gippa_terminatornet_shard.py``; the fix was to drop the
suffix. Python helpers under ``workflow/bin/`` are therefore extensionless
executables with a ``#!/usr/bin/env python3`` shebang (see
``workflow/bin/gippa_terminatornet_shard``); tests load them by path.

Scope: only files the wheel selects, i.e. ``only-include`` minus the wheel's
``exclude`` list (``workflow/tests`` and friends are excluded, so hatchling never
sees them). ``.pyi`` stubs are guarded too: they are Python sources in spirit and
would be just as surprising in the workflow tree.
"""

from __future__ import annotations

import fnmatch
import shutil
import subprocess
import tomllib
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
PYPROJECT = REPO / "pyproject.toml"


def _wheel_excludes() -> list[str]:
    cfg = tomllib.loads(PYPROJECT.read_text(encoding="utf-8"))
    wheel = cfg["tool"]["hatch"]["build"]["targets"]["wheel"]
    return [str(p).rstrip("/") for p in wheel.get("exclude", [])]


def _excluded(path: str, patterns: list[str]) -> bool:
    """Gitignore-style approximation: a pattern matches the path or any parent dir."""
    parts = path.split("/")
    prefixes = ["/".join(parts[: i + 1]) for i in range(len(parts))]
    return any(fnmatch.fnmatch(pre, pat) for pat in patterns for pre in prefixes)


def _tracked_workflow_files() -> list[str]:
    if shutil.which("git") is None:
        pytest.skip("git not on PATH")
    proc = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "-z", "--", "workflow"],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        pytest.skip(f"not a git checkout: {proc.stderr.strip()}")
    return [p for p in proc.stdout.split("\0") if p]


def test_sources_rewrite_still_maps_workflow() -> None:
    """If this mapping goes away the guard below can go too; keep them in sync."""
    cfg = tomllib.loads(PYPROJECT.read_text(encoding="utf-8"))
    sources = cfg["tool"]["hatch"]["build"]["targets"]["wheel"]["sources"]
    assert sources.get("workflow") == "gippa/workflow"


def test_exclude_matcher_examples() -> None:
    pats = _wheel_excludes()
    assert _excluded("workflow/tests/bin/bakta", pats)
    assert _excluded("workflow/tests/helper.py", pats)
    assert not _excluded("workflow/bin/gippa_terminatornet_shard", pats)
    assert not _excluded("workflow/bin/foo.py", pats)


def test_no_python_files_under_workflow() -> None:
    files = _tracked_workflow_files()
    assert files, "git ls-files workflow returned nothing; is this the GIPPA repo?"
    pats = _wheel_excludes()
    offenders = sorted(f for f in files if f.endswith((".py", ".pyi")) and not _excluded(f, pats))
    assert not offenders, (
        "Python source files are tracked under workflow/ and selected by the wheel:\n  "
        + "\n  ".join(offenders)
        + "\n\npyproject.toml maps sources 'workflow' -> 'gippa/workflow'; hatchling's "
        "editable build rejects any .py under a prefix-changing rewrite ('Dev mode "
        "installations are unsupported when any path rewrite in the `sources` option "
        "changes a prefix rather than removes it'), so every `uv sync` (and all CI) "
        "fails. Drop the suffix and keep a python3 shebang, as in "
        "workflow/bin/gippa_terminatornet_shard (PR #455); load it in tests by path."
    )
