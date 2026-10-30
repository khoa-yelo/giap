"""A fresh install resolves the newest release of every dependency, so the numeric stack is capped.

pandas 3.1 removed a name pandarallel 1.6.5 imports, and codon-bias imports pandarallel: an uncapped
wheel crashed at `gippa run` while every test on the locked environment passed.
"""

from __future__ import annotations

import tomllib
from pathlib import Path

import pytest
from packaging.requirements import Requirement

REPO = Path(__file__).resolve().parents[2]


def _requirements() -> dict[str, Requirement]:
    project = tomllib.loads((REPO / "pyproject.toml").read_text())["project"]
    reqs = (Requirement(d) for d in project["dependencies"])
    return {r.name: r for r in reqs}


def test_pandas_excludes_the_release_that_breaks_pandarallel() -> None:
    spec = _requirements()["pandas"].specifier
    assert "3.0.6" in spec
    assert "3.1.0" not in spec


@pytest.mark.parametrize(
    "name",
    ["numpy", "pandas", "pyarrow", "numba", "scipy", "scikit-learn", "statsmodels", "lightgbm"],
)
def test_numeric_dependency_has_an_upper_bound(name: str) -> None:
    assert any(s.operator in {"<", "<=", "~=", "=="} for s in _requirements()[name].specifier)
