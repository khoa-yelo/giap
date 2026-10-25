"""``docs/user-guide/platforms.md`` (M11-09 definition item 1, example 4, acceptance 1).

Spec: docs/specs/M11-09-platforms-docker-smoke.md. Readings used here:

- The matrix is the first Markdown table whose header names "Platform"; rows are found by their
  first cell (Apptainer, conda, Docker on Linux, macOS, Windows); each status cell starts with
  ``tested``, ``smoke-tested`` or ``untested`` (backticks/bold allowed).
- Example 4: a ``tested`` cell names at least one CI job id from ``.github/workflows/*.yml`` or a
  Slurm script under ``scripts/slurm/`` (by file name or stem, e.g. ``nightly``); every
  ``*.sbatch`` or ``scripts/...`` path the page mentions exists.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "user-guide" / "platforms.md"
MKDOCS = REPO / "mkdocs.yml"
WORKFLOWS = REPO / ".github" / "workflows"
SLURM = REPO / "scripts" / "slurm"
STATUSES = ("smoke-tested", "untested", "tested")
ROWS = ("apptainer", "conda", "docker", "macos", "windows")
# definition 1: (Python package, Full pipeline) per row
EXPECTED = {
    "apptainer": ("tested", "tested"),
    "conda": ("tested", "tested"),
    "docker": ("tested", "smoke-tested"),
    "macos": ("tested", "untested"),
    "windows": ("untested", "untested"),
}


@pytest.fixture(scope="module")
def text() -> str:
    assert PAGE.is_file(), "docs/user-guide/platforms.md does not exist"
    return PAGE.read_text(encoding="utf-8")


def _cells(line: str) -> list[str]:
    return [c.strip() for c in line.strip().strip("|").split("|")]


def _table(text: str) -> tuple[list[str], list[list[str]]]:
    lines = text.splitlines()
    for i, ln in enumerate(lines):
        if ln.lstrip().startswith("|") and "platform" in ln.lower():
            header = _cells(ln)
            rows = []
            for row in lines[i + 2 :]:
                if not row.lstrip().startswith("|"):
                    break
                rows.append(_cells(row))
            return header, rows
    raise AssertionError("no Markdown table with a 'Platform' column")


def _key(first: str) -> str | None:
    f = first.lower()
    if "macos" in f:
        return "macos"
    if "windows" in f:
        return "windows"
    if "apptainer" in f:
        return "apptainer"
    if "conda" in f:
        return "conda"
    if "docker" in f and "linux" in f:
        return "docker"
    return None


@pytest.fixture(scope="module")
def matrix(text: str) -> dict[str, tuple[str, str]]:
    header, rows = _table(text)
    assert len(header) == 3, header
    assert "python package" in header[1].lower(), header
    assert "full pipeline" in header[2].lower(), header
    out: dict[str, tuple[str, str]] = {}
    for r in rows:
        assert len(r) == 3, r
        k = _key(r[0])
        assert k is not None, f"unrecognised platform row {r[0]!r}"
        assert k not in out, f"duplicate row for {k}"
        out[k] = (r[1], r[2])
    assert sorted(out) == sorted(ROWS), sorted(out)
    return out


def _status(cell: str) -> str:
    c = cell.lstrip("`*_ ").lower()
    for s in STATUSES:
        if c.startswith(s) and not c[len(s) : len(s) + 1].isalnum():
            return s
    raise AssertionError(f"cell does not start with tested/smoke-tested/untested: {cell!r}")


def _job_ids() -> set[str]:
    ids: set[str] = set()
    for p in WORKFLOWS.glob("*.y*ml"):
        data: dict[str, Any] = yaml.safe_load(p.read_text(encoding="utf-8"))
        ids |= set((data or {}).get("jobs") or {})
    return ids


def _slurm_names() -> set[str]:
    names: set[str] = set()
    for p in SLURM.iterdir():
        if p.is_file():
            names |= {p.name, p.stem, f"scripts/slurm/{p.name}"}
    return names


def test_page_is_filled(text: str) -> None:
    assert "Filled by M11-09" not in text
    assert '!!! note "Stub"' not in text


def test_in_nav_under_user_guide() -> None:
    """Read as text: mkdocs.yml may carry !!python tags that yaml.safe_load rejects."""
    lines = MKDOCS.read_text(encoding="utf-8").splitlines()
    start = next(i for i, ln in enumerate(lines) if re.match(r"^\s*- User guide:\s*$", ln))
    indent = len(lines[start]) - len(lines[start].lstrip())
    section = []
    for ln in lines[start + 1 :]:
        if ln.strip() and len(ln) - len(ln.lstrip()) <= indent:
            break
        section.append(ln.strip())
    assert "- user-guide/platforms.md" in section, section


@pytest.mark.parametrize("row", ROWS)
def test_matrix_statuses(matrix: dict[str, tuple[str, str]], row: str) -> None:
    py, full = matrix[row]
    assert (_status(py), _status(full)) == EXPECTED[row], matrix[row]


def test_tested_cells_name_an_existing_job_or_slurm_script(
    matrix: dict[str, tuple[str, str]],
) -> None:
    """Example 4."""
    known = _job_ids() | _slurm_names()
    for row, cells in matrix.items():
        for cell in cells:
            if _status(cell) != "tested":
                continue
            rest = cell.lstrip("`*_ ")[len("tested") :]
            tokens = {t.strip(".,;:`*()") for t in re.split(r"[\s,;()`]+", rest)}
            assert tokens & known, f"{row}: {cell!r} names no CI job or Slurm script"


def test_docker_smoke_cell_names_the_job(matrix: dict[str, tuple[str, str]]) -> None:
    cell = matrix["docker"][1]
    assert "docker-smoke" in cell, cell
    assert "docker-smoke" in _job_ids(), "the docker-smoke job does not exist"


def test_referenced_scripts_exist(text: str) -> None:
    refs = set(re.findall(r"(?<![\w/])((?:scripts/)[\w./-]+\w|[\w-]+\.sbatch)\b", text))
    for ref in refs:
        path = REPO / ref if ref.startswith("scripts/") else SLURM / ref
        assert path.exists(), ref


def test_matrix_cells_name_existing_profiles(matrix: dict[str, tuple[str, str]]) -> None:
    """``-profile a,b`` in a cell names profiles defined in workflow/nextflow.config."""
    cfg = (REPO / "workflow" / "nextflow.config").read_text(encoding="utf-8")
    for cells in matrix.values():
        for cell in cells:
            for spec in re.findall(r"-profile[ =]([\w,]+)", cell):
                for prof in spec.split(","):
                    assert re.search(rf"^\s*{prof}\s*\{{", cfg, re.M), prof


def test_notes(text: str) -> None:
    low = text.lower()
    for word in ("uniop", "fgsea", "dia-nn", "wsl2", "antismash", "bakta"):
        assert word in low, word
    assert "conda" in low and "docker" in low
    assert "amd64" in low and "emulat" in low  # macOS arm64 runs amd64 images under emulation
