"""M9-26 work item 1 and acceptance "Benchmarks recorded": ``benchmarks/test_plots_bench.py`` times
the ``plots`` and ``report`` steps on a 5-Mb gippasim run (≈ 4,500 genes, ≈ 50,000 placements) and
``benchmarks/baseline.json`` holds its wall time and peak RSS (docs/specs/M9-26-scale-and-
determinism.md; GIPPA_SPEC.md §11.2, §12 "benchmarks against baseline.json on one node class").

Interface choices (M9-26 test author; lead decision 4, review nits N1-N4):

- two flat entries, keyed ``plots_5mb_sim4500`` and ``report_5mb_sim4500``, each naming
  ``benchmarks/test_plots_bench.py::<test>`` in its ``benchmark`` field, with top-level ``secs``
  (wall clock) and ``maxrss_mb`` (MB), ``source`` (the M9-26 spec), ``cpus``, ``commit`` (a bare
  short or full SHA, no prose), ``budget_s`` and a ``note`` saying interpreter start-up is
  excluded; the run size as ``n_genes`` and ``n_placements``; where it ran (``node`` or
  ``node_class``);
- budgets: 120 s for ``plots``, 30 s for ``report`` (the 600-s line stays the K-12 one); the
  benchmark module states them as a module-level dict literal ``BUDGETS_S = {"plots": 120.0,
  "report": 30.0}``;
- the placement lower bound is 35,000 (sim4500 gives ≈ 42k), stated in the benchmark module as
  ``MIN_PLACEMENTS = 35_000``.
"""

from __future__ import annotations

import ast
import json
import re
from pathlib import Path
from typing import Any

import pytest

REPO = Path(__file__).resolve().parents[3]
BENCH = REPO / "benchmarks" / "test_plots_bench.py"
BASELINE = REPO / "benchmarks" / "baseline.json"
KEYS = {"plots": "plots_5mb_sim4500", "report": "report_5mb_sim4500"}
BUDGETS = {"plots": 120.0, "report": 30.0}
MIN_PLACEMENTS = 35_000
STEPS = tuple(KEYS)


def _data() -> dict[str, Any]:
    data: dict[str, Any] = json.loads(BASELINE.read_text(encoding="utf-8"))
    return data


def _entry(step: str) -> dict[str, Any]:
    data = _data()
    assert KEYS[step] in data, f"no baseline.json entry {KEYS[step]!r}"
    entry: dict[str, Any] = data[KEYS[step]]
    return entry


def _bench_constant(name: str) -> Any:
    """A module-level literal of the benchmark module, read without importing it."""
    tree = ast.parse(BENCH.read_text(encoding="utf-8"))
    for node in tree.body:
        targets = [node.target] if isinstance(node, ast.AnnAssign) else getattr(node, "targets", [])
        if any(isinstance(t, ast.Name) and t.id == name for t in targets) and node.value:
            return ast.literal_eval(node.value)
    raise AssertionError(f"{BENCH.relative_to(REPO)} defines no module-level {name}")


def test_plots_benchmark_module_exists() -> None:
    assert BENCH.is_file(), f"{BENCH.relative_to(REPO)} missing"
    text = BENCH.read_text(encoding="utf-8")
    assert "plots" in text and "report" in text


def test_exactly_the_two_flat_entries_name_the_benchmark_module() -> None:
    found = {
        k
        for k, v in _data().items()
        if isinstance(v, dict)
        and str(v.get("benchmark", "")).startswith("benchmarks/test_plots_bench.py::")
    }
    assert found == set(KEYS.values()), found


@pytest.mark.parametrize("step", STEPS)
def test_baseline_names_an_existing_benchmark_test(step: str) -> None:
    name = str(_entry(step)["benchmark"]).split("::", 1)[1]
    assert f"def {name.split('[')[0]}(" in BENCH.read_text(encoding="utf-8"), name


@pytest.mark.parametrize("step", STEPS)
def test_baseline_records_flat_wall_time_and_peak_rss(step: str) -> None:
    entry = _entry(step)
    assert not isinstance(entry.get("plots"), dict) and not isinstance(entry.get("report"), dict)
    assert float(entry["secs"]) > 0
    assert float(entry["maxrss_mb"]) > 0


@pytest.mark.parametrize("step", STEPS)
def test_baseline_records_source_cpus_and_a_real_commit(step: str) -> None:
    entry = _entry(step)
    assert "docs/specs/M9-26-scale-and-determinism.md" in str(entry.get("source", "")), entry
    assert isinstance(entry.get("cpus"), int) and entry["cpus"] >= 1, entry.get("cpus")
    assert re.fullmatch(r"[0-9a-f]{7,40}", str(entry.get("commit", ""))), entry.get("commit")


@pytest.mark.parametrize("step", STEPS)
def test_baseline_note_says_interpreter_start_up_is_excluded(step: str) -> None:
    note = str(_entry(step).get("note", ""))
    assert re.search(r"interpreter start-?up", note, re.IGNORECASE), note


@pytest.mark.parametrize("step", STEPS)
def test_budget_is_120_s_for_plots_and_30_s_for_report(step: str) -> None:
    entry = _entry(step)
    assert float(entry["budget_s"]) == BUDGETS[step]
    assert float(entry["secs"]) <= BUDGETS[step], entry["secs"]
    budgets = _bench_constant("BUDGETS_S")
    assert {k: float(v) for k, v in budgets.items()} == BUDGETS, budgets


def test_benchmark_placement_lower_bound_is_35000() -> None:
    assert int(_bench_constant("MIN_PLACEMENTS")) == MIN_PLACEMENTS


@pytest.mark.parametrize("step", STEPS)
def test_baseline_records_the_5_mb_run_size_and_node(step: str) -> None:
    entry = _entry(step)
    assert 4_000 <= int(entry["n_genes"]) <= 5_500, entry["n_genes"]
    assert MIN_PLACEMENTS <= int(entry["n_placements"]) <= 60_000, entry["n_placements"]
    assert entry.get("node") or entry.get("node_class"), "where the benchmark ran is not recorded"
