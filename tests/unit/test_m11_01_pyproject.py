"""Static packaging config (issue M11-01 items 1-3, acceptance criteria 1 and 4).

Reads ``pyproject.toml``, ``.gitignore`` and ``src/gippa/__init__.py``; never builds (building needs
hatchling from PyPI, so the real build runs in the CI ``package`` job only).
"""

from __future__ import annotations

import fnmatch
import re
import subprocess
import tomllib
from pathlib import Path
from typing import Any

import pytest

import gippa

REPO = Path(__file__).resolve().parents[2]
PYPROJECT: dict[str, Any] = tomllib.loads((REPO / "pyproject.toml").read_text())
HATCH: dict[str, Any] = PYPROJECT.get("tool", {}).get("hatch", {})
BUILD: dict[str, Any] = HATCH.get("build", {})
WHEEL: dict[str, Any] = BUILD.get("targets", {}).get("wheel", {})
SDIST: dict[str, Any] = BUILD.get("targets", {}).get("sdist", {})
SDIST_ALLOWED = {
    "src",
    "workflow",
    "pyproject.toml",
    "README.md",
    "LICENSE",
    "NOTICE.md",
    "CHANGELOG.md",
    "CITATION.cff",
}


def _norm(p: str) -> str:
    return p.strip().strip("/").removeprefix("./")


def _excludes(target: dict[str, Any]) -> list[str]:
    return [*BUILD.get("exclude", []), *target.get("exclude", [])]


def _gitignore() -> list[str]:
    return [
        ln.strip()
        for ln in (REPO / ".gitignore").read_text().splitlines()
        if ln.strip() and not ln.startswith("#")
    ]


def _matches(pattern: str, rel: str) -> bool:
    """Loose gitignore-style match of ``pattern`` against a repo-relative path ``rel``."""
    raw = pattern.strip().removeprefix("./")
    anchored = raw.startswith("/") or "/" in raw.strip("/")
    pat = _norm(raw)
    if pat.startswith("**/"):
        pat, anchored = pat[3:], "/" in pat[3:]
    if pat.endswith("/**"):
        pat = pat[:-3]
    parts = rel.split("/")
    if anchored:
        return any(fnmatch.fnmatch("/".join(parts[: i + 1]), pat) for i in range(len(parts)))
    return any(fnmatch.fnmatch(part, pat) for part in parts)


def _excluded(target: dict[str, Any], rel: str, *, gitignore: bool = True) -> bool:
    pats = _excludes(target) + (_gitignore() if gitignore else [])
    return any(_matches(p, rel) for p in pats)


def _workflow_mapped(target: dict[str, Any]) -> bool:
    """``workflow`` lands at ``gippa/workflow`` through force-include or sources."""
    fi = {_norm(k): _norm(v) for k, v in target.get("force-include", {}).items()}
    src = target.get("sources", {})
    src = {_norm(k): _norm(v) for k, v in src.items()} if isinstance(src, dict) else {}
    return fi.get("workflow") == "gippa/workflow" or src.get("workflow") == "gippa/workflow"


# --- item 1: single-source version ------------------------------------------------------------


def test_project_has_no_literal_version() -> None:
    assert "version" not in PYPROJECT["project"]


def test_version_is_dynamic() -> None:
    assert "version" in PYPROJECT["project"].get("dynamic", [])


def test_hatch_version_path_is_init() -> None:
    assert HATCH.get("version", {}).get("path") == "src/gippa/__init__.py"


def test_init_holds_the_only_literal() -> None:
    text = (REPO / "src/gippa/__init__.py").read_text()
    found = re.findall(r'^__version__\s*=\s*"([^"]+)"\s*$', text, flags=re.M)
    assert found == [gippa.__version__]
    others = [
        p
        for p in (REPO / "src/gippa").rglob("*.py")
        if p.name != "__init__.py" or p.parent != REPO / "src/gippa"
        if re.search(r'^__version__\s*=\s*["\']', p.read_text(), flags=re.M)
    ]
    assert others == []


# --- item 2: the wheel carries workflow/ as gippa/workflow -------------------------------------


def test_wheel_maps_workflow_to_gippa_workflow() -> None:
    assert _workflow_mapped(WHEEL), WHEEL


def test_wheel_still_ships_the_python_package() -> None:
    pkgs = [_norm(p) for p in WHEEL.get("packages", [])]
    srcs = WHEEL.get("sources", {})
    stripped = isinstance(srcs, dict) and _norm(srcs.get("src", "x")) == ""
    stripped = stripped or (isinstance(srcs, list) and "src" in [_norm(x) for x in srcs])
    assert "src/gippa" in pkgs or stripped, WHEEL


@pytest.mark.parametrize(
    "rel",
    [
        "workflow/tests/main.nf.test",
        "workflow/tests/real/assets",
        "workflow/.nf-test/tests/x",
        "workflow/.nf-test.log",
        "workflow/.nextflow/history",
        "workflow/.nextflow.log",
        "workflow/work/ab/cdef/.command.sh",
    ],
)
def test_wheel_excludes_test_and_runtime_dirs(rel: str) -> None:
    assert _excluded(WHEEL, rel), (rel, _excludes(WHEEL))


def test_wheel_excludes_workflow_tests_in_hatch_config_itself() -> None:
    """``workflow/tests`` is tracked by git, so .gitignore cannot drop it: hatch must."""
    assert _excluded(WHEEL, "workflow/tests/main.nf.test", gitignore=False), _excludes(WHEEL)


@pytest.mark.parametrize(
    "rel",
    [
        "workflow/bin/gippa_fgsea.R",
        "workflow/assets/ps/ps",
        "workflow/main.nf",
        "workflow/nextflow.config",
        "workflow/envs/bakta.yml",
        "workflow/conf/test.config",
        "workflow/modules/local/bakta/main.nf",
        "workflow/modules/local/uniop/main.nf",
        "workflow/subworkflows",
    ],
)
def test_wheel_keeps_workflow_content(rel: str) -> None:
    assert not _excluded(WHEEL, rel), (rel, _excludes(WHEEL))


def test_executables_are_tracked_with_mode_755() -> None:
    """Hatch copies the file mode; git must keep +x on workflow/bin/* and assets/ps/ps."""
    out = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "-s", "workflow/bin", "workflow/assets/ps/ps"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.splitlines()
    assert out, "no tracked executables"
    assert all(ln.split()[0] == "100755" for ln in out), out


def test_no_symlinks_in_shipped_workflow() -> None:
    links = [
        p
        for p in (REPO / "workflow").rglob("*")
        if p.is_symlink() and not _excluded(WHEEL, p.relative_to(REPO).as_posix())
    ]
    assert links == []


# --- item 3: sdist content --------------------------------------------------------------------


def test_sdist_selection_is_restricted() -> None:
    chosen = {_norm(p) for p in SDIST.get("only-include", SDIST.get("include", []))}
    assert {"src", "workflow"} <= chosen, SDIST
    assert chosen <= SDIST_ALLOWED, chosen - SDIST_ALLOWED


def test_sdist_names_every_top_level_file_or_hatch_adds_it() -> None:
    """README (``readme``), LICENSE (``license-files``) and pyproject are always added by hatch;
    NOTICE.md, CHANGELOG.md and CITATION.cff must be selected explicitly."""
    chosen = {_norm(p) for p in SDIST.get("only-include", SDIST.get("include", []))}
    assert {"NOTICE.md", "CHANGELOG.md", "CITATION.cff"} <= chosen, chosen


def test_sdist_excludes_workflow_tests() -> None:
    assert _excluded(SDIST, "workflow/tests/main.nf.test", gitignore=False), _excludes(SDIST)
    assert _excluded(SDIST, "workflow/work/ab/x")


def test_no_new_runtime_dependency_for_packaging() -> None:
    deps = {
        re.split(r"[<>=!~ \[;]", d, maxsplit=1)[0].lower()
        for d in PYPROJECT["project"]["dependencies"]
    }
    assert not deps & {"hatchling", "build", "setuptools", "importlib-resources"}
