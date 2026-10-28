"""Confirm-run harness detection (tests/calibration/conftest.py ``HARNESS_NAMES`` /
``reaches_harness``; review of #271, B1). Unmarked: runs on every commit.

A confirm run (``GIPPA_CALIBRATION_CONFIRM``) skips every calibration test that does not reach a
harness entry point, so an entry point missing from ``HARNESS_NAMES`` silently drops block 2 of
its checks and ``--merge`` ends in MERGE_ERROR. Every public ``calibrate*`` function of
``_criteria`` must be listed, and the M7-13b FDR-count tests must be found as harness tests.
"""

from __future__ import annotations

import importlib.util
import inspect
import types
from pathlib import Path
from typing import Any

import _criteria as C
import pytest

HERE = Path(__file__).parent


def _conftest(config: pytest.Config) -> types.ModuleType:
    want = (HERE / "conftest.py").resolve()
    for p in config.pluginmanager.get_plugins():
        f = getattr(p, "__file__", None)
        if isinstance(p, types.ModuleType) and f and Path(f).resolve() == want:
            return p
    raise RuntimeError("tests/calibration/conftest.py is not loaded")


def _load(path: Path, name: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def public_calibrators() -> set[str]:
    return {n for n, f in vars(C).items() if n.startswith("calibrate") and inspect.isfunction(f)
            and f.__module__ == C.__name__}  # fmt: skip


def test_every_public_calibrate_entry_point_is_a_harness_name(
    pytestconfig: pytest.Config,
) -> None:
    names = public_calibrators()
    assert {"calibrate", "calibrate_pooled", "calibrate_fdr_counts"} <= names
    assert names <= set(_conftest(pytestconfig).HARNESS_NAMES), names


def test_fdr_count_tests_reach_the_harness(pytestconfig: pytest.Config) -> None:
    cf = _conftest(pytestconfig)
    gate = _load(HERE / "test_novel_orfs_gate.py", "_novel_gate_for_harness")
    attack = _load(HERE / "attacks" / "test_attack_m7_gate.py", "_attack_m7_for_harness")
    for fn in (gate.test_g3_novel_orf_f1, gate.test_g4_class_fdr, gate.test_g4_class_fdr_full,
               gate.test_g4k3_clustered_false_orfs, gate.test_g4k3_clustered_false_orfs_full,
               attack._gc, attack.test_attack_gc_poor_intergenic_fdp_at_level,
               attack.test_attack_gc_poor_intergenic_fdp_at_level_full):  # fmt: skip
        assert cf.reaches_harness([fn]), fn.__name__
