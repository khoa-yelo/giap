"""M9-26 work item 1 and acceptance "Benchmarks recorded": ``benchmarks/test_plots_bench.py`` times
the ``plots`` and ``report`` steps on a 5-Mb giapsim run (≈ 4,500 genes, ≈ 50,000 placements) and
``benchmarks/baseline.json`` holds its wall time and peak RSS (docs/specs/M9-26-scale-and-
determinism.md; GIAP_SPEC.md §11.2, §12 "benchmarks against baseline.json on one node class").

Interface choices (M9-26 test author; the spec names no keys): the entry is found by its
``benchmark`` field (``benchmarks/test_plots_bench.py::<test>``, as every other entry), not by its
key. It has one sub-record per step, ``plots`` and ``report``, each with ``secs`` (wall clock) and
a peak-RSS figure (``maxrss_mb``, ``peak_rss_mb`` or ``peak_mb``, MB); the run size as
``n_genes`` and ``n_placements``; and where it ran (``node`` or ``node_class``).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

REPO = Path(__file__).resolve().parents[3]
BENCH = REPO / "benchmarks" / "test_plots_bench.py"
BASELINE = REPO / "benchmarks" / "baseline.json"
RSS_KEYS = ("maxrss_mb", "peak_rss_mb", "peak_mb")


def _entry() -> dict[str, Any]:
    data = json.loads(BASELINE.read_text(encoding="utf-8"))
    found = [
        v
        for v in data.values()
        if isinstance(v, dict)
        and str(v.get("benchmark", "")).startswith("benchmarks/test_plots_bench.py::")
    ]
    assert len(found) == 1, f"{len(found)} baseline.json entries for benchmarks/test_plots_bench.py"
    return found[0]


def test_plots_benchmark_module_exists() -> None:
    assert BENCH.is_file(), f"{BENCH.relative_to(REPO)} missing"
    text = BENCH.read_text(encoding="utf-8")
    assert "plots" in text and "report" in text


def test_baseline_names_an_existing_benchmark_test() -> None:
    entry = _entry()
    name = str(entry["benchmark"]).split("::", 1)[1]
    assert f"def {name.split('[')[0]}(" in BENCH.read_text(encoding="utf-8"), name


@pytest.mark.parametrize("step", ["plots", "report"])
def test_baseline_records_wall_time_and_peak_rss_per_step(step: str) -> None:
    rec = _entry().get(step)
    assert isinstance(rec, dict), f"no {step!r} sub-record"
    assert float(rec["secs"]) > 0
    rss = [k for k in RSS_KEYS if k in rec]
    assert rss, f"{step}: no peak RSS among {RSS_KEYS}"
    assert float(rec[rss[0]]) > 0


def test_baseline_records_the_5_mb_run_size_and_node() -> None:
    entry = _entry()
    assert 4_000 <= int(entry["n_genes"]) <= 5_500, entry["n_genes"]
    assert 40_000 <= int(entry["n_placements"]) <= 60_000, entry["n_placements"]
    assert entry.get("node") or entry.get("node_class"), "where the benchmark ran is not recorded"
