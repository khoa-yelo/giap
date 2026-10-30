"""Runtime dependencies of the Bioconda recipe agree with ``pyproject.toml`` (M11-10 definition
item 2, examples 1, 2 and 5, edge cases "extras", "environment markers" and "version pins
duplicated", acceptance 1 and 2).

The recipe is rendered from the real ``pyproject.toml`` (cwd = repository root), and for the
modified-pyproject cases from a throw-away copy of the repository layout in ``tmp_path``
(``pyproject.toml``, ``LICENSE``, ``scripts/release/bioconda_recipe.py``,
``packaging/bioconda/``), so the script finds the modified file whether it reads it from the
working directory or relative to its own path. No network.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path
from typing import Any

import pytest
import yaml
from packaging.requirements import Requirement

REPO = Path(__file__).resolve().parents[2]
SHA = "fedcba9876543210" * 4

# Definition item 2: the explicit PyPI -> conda names that differ or need noting; others identical.
EXPECTED_CONDA_NAME = {
    "pyahocorasick": "pyahocorasick",
    "codon-bias": "codon-bias",
    "scikit-learn": "scikit-learn",
    "pydantic": "pydantic",
    "biopython": "biopython",
}
EXTRA_RUN = {"python", "nextflow", "openjdk"}


def pyproject_deps(path: Path = REPO / "pyproject.toml") -> list[Requirement]:
    data = tomllib.loads(path.read_text())
    return [Requirement(d) for d in data["project"]["dependencies"]]


def copy_layout(dst: Path) -> Path:
    (dst / "scripts" / "release").mkdir(parents=True)
    shutil.copy2(REPO / "scripts" / "release" / "bioconda_recipe.py", dst / "scripts" / "release")
    shutil.copytree(REPO / "packaging" / "bioconda", dst / "packaging" / "bioconda")
    for name in ("pyproject.toml", "LICENSE", "README.md"):
        shutil.copy2(REPO / name, dst / name)
    return dst


def add_dependency(root: Path, dep: str) -> None:
    py = root / "pyproject.toml"
    text = py.read_text()
    marker = "\ndependencies = [\n"
    assert text.count(marker) == 1
    new = text.replace(marker, f'{marker}  "{dep}",\n')
    py.write_text(new)


def run_render(root: Path, out: Path, *extra: str) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([str(REPO / "src"), env.get("PYTHONPATH", "")])
    return subprocess.run(
        [
            sys.executable,
            str(root / "scripts" / "release" / "bioconda_recipe.py"),
            "--version",
            "1.0.0",
            "--sha256",
            SHA,
            *extra,
            "--out",
            str(out),
        ],
        cwd=root,
        env=env,
        capture_output=True,
        text=True,
        timeout=120,
    )


def run_entries(out: Path) -> list[str]:
    data: dict[str, Any] = yaml.safe_load((out / "gippa" / "meta.yaml").read_text())
    return [" ".join(str(r).split()) for r in data["requirements"]["run"]]


def split_entry(entry: str) -> tuple[str, str]:
    name, _, spec = entry.partition(" ")
    return name, spec.replace(" ", "")


@pytest.fixture(scope="module")
def main_run(tmp_path_factory: pytest.TempPathFactory) -> list[str]:
    out = tmp_path_factory.mktemp("main")
    proc = run_render(REPO, out)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return run_entries(out)


def test_example_1_run_count(main_run: list[str]) -> None:
    """Example 1: every runtime dependency (21 on main) + python + nextflow + openjdk."""
    deps = pyproject_deps()
    assert len(main_run) == len(deps) + 3
    names = [split_entry(e)[0] for e in main_run]
    assert len(names) == len(set(names)), "duplicate run entries"
    assert set(names) >= EXTRA_RUN


def test_acceptance_2_every_dependency_with_same_lower_bound(main_run: list[str]) -> None:
    by_name = dict(split_entry(e) for e in main_run)
    for req in pyproject_deps():
        conda = EXPECTED_CONDA_NAME.get(req.name, req.name)
        assert conda in by_name, f"{req.name} missing from run:"
        assert by_name[conda] == str(req.specifier).replace(" ", ""), (
            f"{req.name}: pyproject {req.specifier!s} vs recipe {by_name[conda]!r}"
        )


def test_run_has_nothing_beyond_dependencies(main_run: list[str]) -> None:
    allowed = {EXPECTED_CONDA_NAME.get(r.name, r.name) for r in pyproject_deps()} | EXTRA_RUN
    assert {split_entry(e)[0] for e in main_run} == allowed


def test_codon_bias_present(main_run: list[str]) -> None:
    """codon-bias is bioconda-only, but a bioconda recipe may depend on it."""
    assert "codon-bias" in [split_entry(e)[0] for e in main_run]


def test_example_5_conda_spacing(main_run: list[str]) -> None:
    raw = set(main_run)
    assert "pydantic >=2" in raw
    assert "click >=8.1" in raw
    assert "pydantic>=2" not in raw and "click>=8.1" not in raw


def test_extras_not_in_recipe(main_run: list[str]) -> None:
    data = tomllib.loads((REPO / "pyproject.toml").read_text())
    extras = {
        Requirement(d).name
        for group in data["project"]["optional-dependencies"].values()
        for d in group
    }
    runtime = {r.name for r in pyproject_deps()}
    names = {split_entry(e)[0] for e in main_run}
    assert not (extras - runtime) & names


def test_example_2_unmapped_dependency_fails(tmp_path: Path) -> None:
    root = copy_layout(tmp_path / "repo")
    add_dependency(root, "newlib>=1")
    out = tmp_path / "out"
    proc = run_render(root, out)
    assert proc.returncode != 0
    assert "no conda name for newlib (add it to the map)" in proc.stdout + proc.stderr
    assert not (out / "gippa" / "meta.yaml").exists()


def test_pin_change_propagates(tmp_path: Path) -> None:
    """Edge case "version pins duplicated": a bound changed in pyproject shows up in the recipe."""
    root = copy_layout(tmp_path / "repo")
    py = root / "pyproject.toml"
    text = py.read_text()
    assert '"numpy<2.6",' in text
    py.write_text(text.replace('"numpy<2.6",', '"numpy>=1.26",', 1))
    out = tmp_path / "out"
    proc = run_render(root, out)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "numpy >=1.26" in run_entries(out)


def test_environment_marker_fails(tmp_path: Path) -> None:
    """Edge case: a marker cannot become a selector in a noarch: python recipe."""
    root = copy_layout(tmp_path / "repo")
    py = root / "pyproject.toml"
    text = py.read_text()
    py.write_text(text.replace('"numpy<2.6",', "\"numpy; python_version < '3.12'\",", 1))
    assert "python_version" in py.read_text()
    out = tmp_path / "out"
    proc = run_render(root, out)
    assert proc.returncode != 0
    assert "numpy" in proc.stdout + proc.stderr
    assert not (out / "gippa" / "meta.yaml").exists()
