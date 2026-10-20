"""The 15 % benchmark gate of ``benchmarks/test_annotation_steps.py`` (M2-21 lead decisions 10, 13).

Decision 10/13: the slowdown gate (``benchmark_slowdown_max`` of ``ci/gates.yaml``) applies only to
steps whose baseline is at least 1 s; shorter steps (``context_layers_sim100k`` at 0.15 s on wheat)
are recorded, not gated.

The benchmark module is imported from its file and its ``_bench(benchmark, name, run)`` helper is
driven with a fake ``benchmark`` fixture (``pedantic`` calls the function ``rounds +
warmup_rounds`` times) and a fake clock, so every round takes exactly the chosen median; the
baseline lookup is replaced. No step runs. (Test author: the gate stays inside ``_bench``; a
refactor that moves it must keep ``_bench``'s behaviour.)
"""

from __future__ import annotations

import importlib.util
import json
import sys
from collections.abc import Callable
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Any

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
BENCH = REPO / "benchmarks" / "test_annotation_steps.py"
SLOWDOWN = float(yaml.safe_load((REPO / "ci" / "gates.yaml").read_text())["benchmark_slowdown_max"])


@pytest.fixture(scope="module")
def mod() -> ModuleType:
    spec = importlib.util.spec_from_file_location("m2_21_bench_annotation_steps", BENCH)
    assert spec and spec.loader
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


class FakeBenchmark:
    def __init__(self) -> None:
        self.extra_info: dict[str, Any] = {}

    def pedantic(self, fn: Callable[[], object], rounds: int = 1, iterations: int = 1,
                 warmup_rounds: int = 0, **_: Any) -> None:  # fmt: skip
        for _ in range((rounds + warmup_rounds) * iterations):
            fn()


def _clock(step_s: float) -> SimpleNamespace:
    """A clock advancing ``step_s`` per call: a round timed by two calls lasts ``step_s``."""
    t = [0.0]

    def pc() -> float:
        t[0] += step_s
        return t[0]

    return SimpleNamespace(perf_counter=pc, perf_counter_ns=lambda: int(pc() * 1e9),
                           monotonic=pc, time=pc, sleep=lambda s: None)  # fmt: skip


def _bench(mod: ModuleType, monkeypatch: pytest.MonkeyPatch, base: float | None,
           median: float) -> FakeBenchmark:  # fmt: skip
    """Run ``_bench`` with every round lasting ``median`` seconds."""
    monkeypatch.setattr(mod, "baseline", lambda name, klass: base)
    monkeypatch.setattr(mod, "time", _clock(median))
    b = FakeBenchmark()
    mod._bench(b, "step_under_test", lambda i: None)
    return b


def test_fake_clock_gives_the_chosen_median(
    mod: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    b = _bench(mod, monkeypatch, None, 2.5)
    assert b.extra_info.get("median_s") == pytest.approx(2.5, abs=1e-3), b.extra_info


@pytest.mark.parametrize(
    ("base", "median"),
    [
        (0.1519, 1.0),  # context_layers_sim100k on wheat: recorded, not gated
        (0.5, 50.0),
        (0.999, 2.0),
    ],
)
def test_short_baselines_are_recorded_not_gated(
    mod: ModuleType, monkeypatch: pytest.MonkeyPatch, base: float, median: float
) -> None:
    b = _bench(mod, monkeypatch, base, median)
    assert b.extra_info.get("median_s") == pytest.approx(median, abs=1e-3), b.extra_info


@pytest.mark.parametrize(
    ("base", "factor"),
    [
        (1.0, 1 + SLOWDOWN + 0.05),  # exactly 1 s is gated (">= 1 s")
        (1.2724, 1 + SLOWDOWN + 0.05),  # producer_chain_sim100k on wheat
        (30.0, 2.0),
    ],
)
def test_long_baselines_fail_beyond_the_slowdown(
    mod: ModuleType, monkeypatch: pytest.MonkeyPatch, base: float, factor: float
) -> None:
    with pytest.raises(AssertionError):
        _bench(mod, monkeypatch, base, base * factor)


@pytest.mark.parametrize("base", [1.0, 1.2724, 30.0])
def test_long_baselines_pass_within_the_slowdown(
    mod: ModuleType, monkeypatch: pytest.MonkeyPatch, base: float
) -> None:
    _bench(mod, monkeypatch, base, base * (1 + SLOWDOWN - 0.05))


def test_recorded_wheat_baselines_split_at_one_second() -> None:
    """The two committed wheat entries sit on either side of the 1 s line (decision 10)."""
    wheat = json.loads((REPO / "benchmarks" / "baseline.json").read_text())["annotation_steps"][
        "wheat"
    ]
    assert wheat["context_layers_sim100k"]["median_s"] < 1.0
    assert wheat["producer_chain_sim100k"]["median_s"] >= 1.0
