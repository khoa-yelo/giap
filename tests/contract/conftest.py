"""Contract-test configuration.

``requires_r`` (M8-17a, docs/specs/M8-17-r-fgsea-engine.md §6) marks tests that run the real
``gippa-fgsea`` R environment; they run only on Slurm (``scripts/slurm/test_branch.sbatch -m
requires_r``). Registered here until it is added to ``[tool.pytest.ini_options] markers``.
"""

from __future__ import annotations

import pytest


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers", "requires_r: runs R fgsea in the gippa-fgsea env (Slurm only, never in CI)"
    )
