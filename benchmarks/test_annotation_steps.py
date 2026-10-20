"""Runtime of the in-process annotation steps on SIM100K (M2-21; GIAP_SPEC.md §12 benchmarks).

Two benchmarks, both through ``giap.steps.registry.run_step`` and no external tool (the tool
outputs are the recorded / recorded-format files of ``tests/helpers/layers_chain.py``):

* ``producer_chain_sim100k``: ``layers_chain.build_chain`` (giapsim SIM100K, then ``parse_bakta``,
  ``uniop_input``, ``operon_pairs``, ``codon_usage``, ``terminator_layer``, ``region_layer``,
  ``promoter_windows``, ``promoter_layer``, ``tx_layer``, ``ostir_input``, ``rbs_layer``);
* ``context_layers_sim100k``: the ``context_layers`` step on that chain.

Median of 5 rounds after one warm-up round. ``benchmarks/baseline.json`` > ``annotation_steps``
> ``<node class>`` holds the accepted medians per node class (the class is the node-name
prefix: wheat, barley, rye, oat; the nightly pins wheat, lead decision 1); a median more than
``benchmark_slowdown_max`` (``ci/gates.yaml``, 15 %) above the baseline of the same class
fails. Without a baseline for the class, or with a baseline under ``GATE_MIN_S`` (1 s; lead
decisions 10 and 13: timer noise dominates such steps), the medians are only reported
(``extra_info``, ``--benchmark-json``).
Run with ``pytest benchmarks/test_annotation_steps.py --benchmark-only``.
"""

from __future__ import annotations

import json
import socket
import statistics
import sys
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))

from helpers.layers_chain import Chain, build_chain, fresh_steps  # noqa: E402

from giap.config import GiapConfig  # noqa: E402

ROUNDS = 5
GATE_MIN_S = 1.0  # baselines below this are recorded, not gated (lead decisions 10, 13)
SLOWDOWN_MAX = float(
    yaml.safe_load((ROOT / "ci" / "gates.yaml").read_text())["benchmark_slowdown_max"]
)


def node_class(host: str | None = None) -> str:
    """``wheat-05`` → ``wheat``."""
    name = (host or socket.gethostname()).split(".")[0]
    return name.split("-")[0].lower()


def baseline(name: str, klass: str) -> float | None:
    """Accepted median (s) of ``name`` on node class ``klass``, if recorded."""
    data = json.loads((ROOT / "benchmarks" / "baseline.json").read_text())
    entry = data.get("annotation_steps", {}).get(klass, {}).get(name)
    return float(entry["median_s"]) if entry else None


def _bench(benchmark: Any, name: str, run: Callable[[int], object]) -> None:
    times: list[float] = []
    count = iter(range(1_000_000))

    def once() -> None:
        i = next(count)
        t0 = time.perf_counter()
        run(i)
        times.append(time.perf_counter() - t0)

    benchmark.pedantic(once, rounds=ROUNDS, iterations=1, warmup_rounds=1)
    median = statistics.median(times[-ROUNDS:])  # the warm-up round is not a measurement
    klass = node_class()
    benchmark.extra_info.update({"median_s": round(median, 4), "node_class": klass,
                                 "node": socket.gethostname().split(".")[0]})  # fmt: skip
    base = baseline(name, klass)
    if base is not None and base >= GATE_MIN_S:
        limit = base * (1 + SLOWDOWN_MAX)
        assert median <= limit, (
            f"{name}: median {median:.3f} s > {limit:.3f} s (baseline {base:.3f} s on {klass})"
        )


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return build_chain(tmp_path_factory.mktemp("bench_chain"))


def test_producer_chain_sim100k(benchmark: Any, tmp_path: Path) -> None:
    _bench(benchmark, "producer_chain_sim100k", lambda i: build_chain(tmp_path / f"r{i}"))


def test_context_layers_sim100k(benchmark: Any, chain: Chain, tmp_path: Path) -> None:
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)  # registry reload outside the timed rounds

        def run(i: int) -> None:
            out = tmp_path / f"r{i}"
            outs = {"annotation": out / "annotation", "layers": out / "layers.parquet"}
            reg.run_step("context_layers", GiapConfig(seed=1), dict(chain.inputs), outs)
            assert outs["layers"].is_file()

        _bench(benchmark, "context_layers_sim100k", run)
