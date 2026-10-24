"""``giap.commit`` and ``giap.dirty`` of ``<a>.convert.json`` (M10-03 review round 1, B3).

Spec: docs/specs/M10-03-diann-convert.md, "``validation/pilot/diann_convert.py summary``" (the
``giap`` field) and Q9 (running from the worktree is accepted because the commit and the dirty
flag are recorded). Lead decisions on B3:

* ``dirty`` counts tracked edits and untracked files under ``src``, ``validation``, ``tests`` and
  ``scripts`` (the code the job executes), but not ``__pycache__`` or ``*.pyc``;
* a git failure records ``"commit": null, "dirty": null`` (never ``""`` / ``false``);
* the repository is the checkout of the running code (``_REPO``, from the module's own path), not
  ``$GIAP_REPO``.

Each test copies ``validation/pilot/diann_convert.py`` into a temporary checkout (or a plain
directory) and calls its ``_giap()``. No ``.gitignore`` is written, so bytecode is excluded by the
module, not by git.
"""

from __future__ import annotations

import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import pytest
from helpers.diann_convert_m10 import SUMMARY
from helpers.slurm_sandbox import REPO

DIRS = ("src", "validation", "tests", "scripts")


def _git(repo: Path, *args: str) -> str:
    cmd = ["git", "-c", "user.name=t", "-c", "user.email=t@example.invalid", "-C", str(repo)]
    return subprocess.run([*cmd, *args], check=True, capture_output=True, text=True).stdout.strip()


def _load(where: Path, monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    """Import ``where/validation/pilot/diann_convert.py`` (its ``_REPO`` is ``where``)."""
    dst = where / "validation" / "pilot" / "diann_convert.py"
    if not dst.is_file():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(SUMMARY, dst)
    monkeypatch.setattr(sys, "path", list(sys.path))  # the module inserts <where>/src
    spec = importlib.util.spec_from_file_location(f"m10_03_prov_{abs(hash(where))}", dst)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert Path(mod._REPO) == where.resolve(), (mod._REPO, where)
    return mod


@pytest.fixture
def checkout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A committed checkout: one tracked file in each code directory, and ``docs/``."""
    r = (tmp_path / "checkout").resolve()
    r.mkdir()
    _git(r, "init", "-q", "-b", "main")
    for d in DIRS:
        (r / d).mkdir(parents=True, exist_ok=True)
        (r / d / "keep.txt").write_text("x\n")
    (r / "docs").mkdir()
    (r / "docs" / "notes.md").write_text("x\n")
    (r / "validation" / "pilot").mkdir(parents=True)
    shutil.copy2(SUMMARY, r / "validation" / "pilot" / "diann_convert.py")
    _git(r, "add", "-A")
    _git(r, "commit", "-q", "-m", "init")
    monkeypatch.delenv("GIAP_REPO", raising=False)
    return r


def test_clean_checkout(checkout: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    g = _load(checkout, monkeypatch)._giap()
    assert g["commit"] == _git(checkout, "rev-parse", "HEAD")
    assert g["dirty"] is False


def test_tracked_edit_is_dirty(checkout: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (checkout / "src" / "keep.txt").write_text("edited\n")
    assert _load(checkout, monkeypatch)._giap()["dirty"] is True


@pytest.mark.parametrize("d", DIRS)
def test_untracked_code_file_is_dirty(
    checkout: Path, monkeypatch: pytest.MonkeyPatch, d: str
) -> None:
    """Review B3: an uncommitted new module the job imports must not read as clean."""
    (checkout / d / "giap_new").mkdir()
    (checkout / d / "giap_new" / "untracked_new.py").write_text("x = 1\n")
    assert _load(checkout, monkeypatch)._giap()["dirty"] is True


@pytest.mark.parametrize("rel", ["docs/draft.md", "NOTES.txt", "local/pilot/inputs/a1.tsv.gz"])
def test_untracked_file_elsewhere_is_clean(
    checkout: Path, monkeypatch: pytest.MonkeyPatch, rel: str
) -> None:
    p = checkout / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("x\n")
    assert _load(checkout, monkeypatch)._giap()["dirty"] is False


@pytest.mark.parametrize(
    "rel",
    [
        "src/giap/__pycache__/cli.cpython-312.pyc",
        "validation/pilot/__pycache__/diann_convert.cpython-312.pyc",
        "tests/__pycache__/conftest.cpython-312.pyc",
        "scripts/stray.pyc",
    ],
)
def test_bytecode_is_not_dirty(checkout: Path, monkeypatch: pytest.MonkeyPatch, rel: str) -> None:
    """``python -I`` writes bytecode next to the code it runs; that is not an edit."""
    p = checkout / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(b"\0")
    assert _load(checkout, monkeypatch)._giap()["dirty"] is False


def test_repo_is_the_running_checkout_not_giap_repo(
    checkout: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Review B3: ``$GIAP_REPO`` naming another (dirty) checkout does not change the record."""
    other = tmp_path / "other"
    other.mkdir()
    _git(other, "init", "-q", "-b", "main")
    (other / "src").mkdir()
    (other / "src" / "a.txt").write_text("x\n")
    _git(other, "add", "-A")
    _git(other, "commit", "-q", "-m", "other")
    (other / "src" / "a.txt").write_text("dirty\n")
    monkeypatch.setenv("GIAP_REPO", str(other))
    g = _load(checkout, monkeypatch)._giap()
    assert g["commit"] == _git(checkout, "rev-parse", "HEAD") != _git(other, "rev-parse", "HEAD")
    assert g["dirty"] is False


@pytest.mark.parametrize("giap_repo", ["unset", "valid"])
def test_git_failure_records_null(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, giap_repo: str
) -> None:
    """Review B3: not a git checkout gives ``commit: null, dirty: null``, not ``""``/``false``,
    also when ``$GIAP_REPO`` names a valid checkout."""
    plain = (tmp_path / "plain").resolve()
    plain.mkdir()
    monkeypatch.setenv("GIT_CEILING_DIRECTORIES", str(tmp_path.resolve()))
    if giap_repo == "valid":
        monkeypatch.setenv("GIAP_REPO", str(REPO))
    else:
        monkeypatch.delenv("GIAP_REPO", raising=False)
    g = _load(plain, monkeypatch)._giap()
    assert g["commit"] is None, g
    assert g["dirty"] is None, g
