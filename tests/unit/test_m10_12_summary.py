"""``scripts/dev/benchmark_summary.py`` (docs/specs/M10-12-benchmark-standard.md, criteria 6-7).

Pure functions on ``bench_run.json`` records (examples 3, 4, 5, 10, the tie case, ``no_data``,
``inputs_differ``) and the thin CLI (``summary.json``/``summary.md``, ``--write`` into a copy of
the real ``benchmarks/baseline.json``). The numbers are the spec's hand-made ones.

Conventions settled by the test author (the spec names the functions, not every argument):

* ``cell_summary(runs, estimate)`` takes **every** record of one cell (all repeats and
  attempts) and does the attempt selection itself; ``estimate`` is ``ESTIMATES[variant]``, a
  module-level mapping keyed by variant (``standard``: 2,100 s, range 960-3,300; ``core``:
  1,500 s, no range). The result is a dict with at least ``status``, ``median_s``,
  ``mean_of_valid_s``, ``spread``, ``noisy``, ``ratio_to_estimate``, ``within_estimate_range``
  and ``superseded`` (run-directory names: ``r2`` for attempt ``a``, ``r2.b`` for ``b``).
* ``baseline_entry(cell)`` of a ``complete`` cell has exactly the documented keys.
* ``markdown_rows(cells)`` returns the rows as one string or a list of strings.
* ``summary.json`` holds the cells under ``cells``, keyed ``<genome>_<variant>`` or as a list of
  dicts with ``genome`` and ``variant``.
"""

from __future__ import annotations

import functools
import importlib.util
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
import yaml
from helpers.bench_m10_12 import (
    BASELINE,
    BASELINE_KEYS,
    GATES,
    GIB,
    SUMMARY,
    bench_run,
    write_runs,
)
from hypothesis import given, settings
from hypothesis import strategies as st

CELLS = ("k12_standard", "k12_core", "sim5m_standard", "sim5m_core")


@functools.cache
def _mod() -> ModuleType:
    assert SUMMARY.is_file(), "scripts/dev/benchmark_summary.py does not exist (M10-12)"
    spec = importlib.util.spec_from_file_location("benchmark_summary", SUMMARY)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules["benchmark_summary"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def bs() -> ModuleType:
    return _mod()


def _cell(bs: ModuleType, runs: list[dict[str, Any]], variant: str = "standard") -> dict[str, Any]:
    return bs.cell_summary(runs, bs.ESTIMATES[variant])


def _three(walls: list[float], genome: str = "k12", variant: str = "standard", **kw: Any):  # type: ignore[no-untyped-def]
    return [bench_run(genome, variant, k + 1, wall_s=w, **kw) for k, w in enumerate(walls)]


def _md(bs: ModuleType, cells: list[dict[str, Any]]) -> str:
    out = bs.markdown_rows(cells)
    return out if isinstance(out, str) else "\n".join(out)


# ================================================================== the gate and the estimates


def test_noisy_threshold_is_read_from_the_gate_file(bs: ModuleType) -> None:
    gate = yaml.safe_load(GATES.read_text())["benchmark_slowdown_max"]
    assert gate == 0.15
    src = SUMMARY.read_text()
    assert "benchmark_slowdown_max" in src and "gates.yaml" in src
    assert not re.search(r"(?<![\d.])0?\.15(?!\d)", src), "read, not hard-coded"


def test_estimates_are_the_spec_13_2_numbers(bs: ModuleType) -> None:
    assert set(bs.ESTIMATES) >= {"standard", "core"}
    std = _cell(bs, _three([2000, 2000, 2000]), "standard")
    assert std["ratio_to_estimate"] == pytest.approx(2000 / 2100)
    assert std["within_estimate_range"] is True
    core = _cell(bs, _three([1500, 1500, 1500], variant="core"), "core")
    assert core["ratio_to_estimate"] == pytest.approx(1.0)
    assert core["within_estimate_range"] is None


# ================================================================== examples 3, 4, 5, 10


def test_example3_k12_standard(bs: ModuleType) -> None:
    cell = _cell(bs, _three([3912, 3700, 4106]))
    assert cell["status"] == "complete"
    assert cell["median_s"] == 3912
    assert cell["spread"] == pytest.approx((4106 - 3700) / 3912, abs=1e-4)
    assert cell["spread"] == pytest.approx(0.1038, abs=1e-4)
    assert cell["noisy"] is False
    assert cell["ratio_to_estimate"] == pytest.approx(1.863, abs=1e-3)
    assert cell["within_estimate_range"] is False
    entry = bs.baseline_entry(cell)
    assert set(entry) == BASELINE_KEYS, set(entry) ^ BASELINE_KEYS
    assert entry["median_s"] == 3912 and entry["runs_s"] == [3912, 3700, 4106]
    assert entry["noisy"] is False and entry["estimate_s"] == 2100
    assert entry["ratio_to_estimate"] == pytest.approx(1.863, abs=1e-3)
    assert entry["within_estimate_range"] is False
    assert entry["cpu"] == "Intel Xeon Gold 6148 (CPU_SKU:6148)"
    assert (entry["threads"], entry["mem_mb"], entry["task_memory"]) == (16, 64000, "56.GB")
    assert entry["profile"] == "farmshare_local,apptainer"
    assert entry["rounds"] == 3 and len(entry["nodes"]) == 3 and len(entry["slurm_jobs"]) == 3
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", str(entry["date"])), entry["date"]
    assert entry["commit"] and "inputs_sha256" in entry and "adr" in entry


def test_example4_k12_core_noisy(bs: ModuleType) -> None:
    cell = _cell(bs, _three([1500, 1520, 1800], variant="core"), "core")
    assert cell["status"] == "complete" and cell["median_s"] == 1520
    assert cell["spread"] == pytest.approx(300 / 1520, abs=1e-4)
    assert cell["noisy"] is True
    assert cell["ratio_to_estimate"] == pytest.approx(1.013, abs=1e-3)
    assert cell["within_estimate_range"] is None
    entry = bs.baseline_entry(cell)
    assert entry["noisy"] is True and entry["median_s"] == 1520
    assert "noisy (spread 19.7 % > 15 %)" in _md(bs, [cell])


def _example5(with_b: bool) -> list[dict[str, Any]]:
    runs = [
        bench_run("sim5m", "core", 1, wall_s=1490),
        bench_run("sim5m", "core", 2, wall_s=1200, status="failed", exit_code=4, n_failed=1),
        bench_run("sim5m", "core", 3, wall_s=1510),
    ]
    if with_b:
        runs.append(bench_run("sim5m", "core", 2, "b", wall_s=1505))
    return runs


def test_example5_incomplete_without_the_replacement(bs: ModuleType) -> None:
    cell = _cell(bs, _example5(False), "core")
    assert cell["status"] == "incomplete"
    assert cell["median_s"] is None
    assert cell["mean_of_valid_s"] == 1500
    assert "2 of 3 repeats" in _md(bs, [cell])


def test_example5_complete_after_r2b(bs: ModuleType) -> None:
    cell = _cell(bs, _example5(True), "core")
    assert cell["status"] == "complete"
    assert cell["median_s"] == 1505
    assert list(cell["superseded"]) == ["r2"]
    assert bs.baseline_entry(cell)["runs_s"] == [1490, 1505, 1510]


def test_example10_per_step_medians(bs: ModuleType) -> None:
    rt, peak = [600, 590, 610], [20.1, 20.6, 20.4]
    runs = [
        bench_run("k12", "standard", k + 1, wall_s=w,
                  steps={"bakta": {"realtime_s": rt[k], "peak_rss_bytes": int(peak[k] * GIB)},
                         "operons": {"realtime_s": 100 + k, "peak_rss_bytes": GIB}})
        for k, w in enumerate([3912, 3700, 4106])
    ]  # fmt: skip
    steps = bs.baseline_entry(_cell(bs, runs))["steps"]
    assert list(steps) == ["bakta", "operons"], "first-start order of the median repeat"
    assert steps["bakta"]["median_s"] == 600
    assert steps["bakta"]["peak_rss_max_bytes"] == int(20.6 * GIB)
    assert steps["operons"]["median_s"] == 101


def test_companion_numbers(bs: ModuleType) -> None:
    runs = [
        bench_run(wall_s=1000, repeat=1, max_rss_bytes=20 * GIB, task_peak=10 * GIB,
                  total_cpu_s=0.70 * 1000 * 16, node="wheat-05"),
        bench_run(wall_s=1100, repeat=2, max_rss_bytes=None, task_peak=12 * GIB,
                  total_cpu_s=0.80 * 1100 * 16, node="wheat-07"),
        bench_run(wall_s=1050, repeat=3, max_rss_bytes=22 * GIB, task_peak=11 * GIB,
                  total_cpu_s=0.75 * 1050 * 16, node="wheat-05"),
    ]  # fmt: skip
    entry = bs.baseline_entry(_cell(bs, runs))
    assert entry["peak_rss_sacct_max_bytes"] == 22 * GIB, "max over repeats; null skipped"
    assert entry["task_peak_rss_max_bytes"] == 12 * GIB
    assert entry["cpu_efficiency_median"] == pytest.approx(0.75)
    assert entry["nodes"] == ["wheat-05", "wheat-07", "wheat-05"]


# ================================================================== edge cases


def test_ties(bs: ModuleType) -> None:
    cell = _cell(bs, _three([1500, 1500, 1800], variant="core"), "core")
    assert cell["median_s"] == 1500 and cell["spread"] == pytest.approx(0.2)


def test_spread_at_the_gate_is_not_noisy(bs: ModuleType) -> None:
    cell = _cell(bs, _three([1000, 1000, 1150]))
    assert cell["spread"] == pytest.approx(0.15) and cell["noisy"] is False


def test_no_data(bs: ModuleType) -> None:
    cell = _cell(bs, [])
    assert cell["status"] == "no_data" and cell["median_s"] is None


@pytest.mark.parametrize("status", ["invalid", "timeout", "failed"])
def test_non_ok_runs_are_excluded(bs: ModuleType, status: str) -> None:
    runs = _three([1000, 1000, 1000])
    reasons = ["cached_tasks:3"] if status == "invalid" else None
    runs[1] = bench_run(repeat=2, wall_s=99999, status=status, invalid_reasons=reasons)
    cell = _cell(bs, runs)
    assert cell["status"] == "incomplete" and cell["mean_of_valid_s"] == 1000


def test_inputs_differ(bs: ModuleType) -> None:
    runs = _three([1000, 1010, 1020], genome="sim5m")
    runs[2]["inputs"]["genome"]["sha256"] = "b" * 64
    cell = _cell(bs, runs)
    assert cell["status"] == "incomplete"
    assert "inputs_differ" in json.dumps(cell)


@pytest.mark.parametrize(
    ("attempts", "wall"),
    [([("a", "ok", 1500), ("b", "ok", 1600)], 1600),
     ([("a", "ok", 1500), ("b", "failed", 1600)], 1500),
     ([("a", "failed", 1), ("b", "invalid", 2), ("c", "ok", 1550)], 1550)],
    ids=["latest_valid_wins", "invalid_b_keeps_a", "c_after_two"],
)  # fmt: skip
def test_latest_valid_attempt(bs: ModuleType, attempts: list[tuple[str, str, int]],
                              wall: int) -> None:  # fmt: skip
    runs = [bench_run(repeat=1, wall_s=1000), bench_run(repeat=3, wall_s=2000)]
    runs += [bench_run(repeat=2, attempt=a, wall_s=w, status=s) for a, s, w in attempts]
    cell = _cell(bs, runs)
    assert cell["status"] == "complete"
    assert bs.baseline_entry(cell)["runs_s"] == [1000, wall, 2000]


# ============================================================ review of PR #300 (lead decisions)


def test_step_order_comes_from_the_median_repeat(bs: ModuleType) -> None:
    """When the repeats' first-start orders differ, the median repeat's order wins (finding 12,
    M10). The median (3,912 s) is r2, not the first record."""

    def steps(first: str, second: str) -> dict[str, dict[str, float]]:
        return {first: {"realtime_s": 600.0, "peak_rss_bytes": GIB, "order": 0},
                second: {"realtime_s": 180.0, "peak_rss_bytes": GIB, "order": 1}}  # fmt: skip

    runs = [bench_run(repeat=1, wall_s=3700, steps=steps("operons", "bakta")),
            bench_run(repeat=2, wall_s=3912, steps=steps("bakta", "operons")),
            bench_run(repeat=3, wall_s=4106, steps=steps("operons", "bakta"))]  # fmt: skip
    cell = _cell(bs, runs)
    assert cell["median_s"] == 3912
    assert list(bs.baseline_entry(cell)["steps"]) == ["bakta", "operons"]


def test_mixed_node_class(bs: ModuleType) -> None:
    """Repeats on two node classes are never mixed into one median (M11)."""
    runs = _three([1000, 1010, 1020])
    runs[1]["cpu_sku"] = "CPU_SKU:9384X"
    cell = _cell(bs, runs)
    assert cell["status"] == "incomplete" and cell["median_s"] is None, cell
    assert "mixed_node_class" in cell["reasons"], cell["reasons"]
    with pytest.raises(ValueError):
        bs.baseline_entry(cell)


def test_an_off_class_cell_is_never_summarised(bs: ModuleType) -> None:
    """Off-class runs are invalid and never summarised (lead, 2026-10-24, nit 20): even three
    records labelled ok on another SKU give no complete cell and no baseline entry."""
    runs = _three([1000, 1010, 1020])
    for r in runs:
        r["cpu_sku"] = "CPU_SKU:9384X"
    cell = _cell(bs, runs)
    assert cell["status"] != "complete" and cell["median_s"] is None, cell


def _set_flag(run: dict[str, Any], name: str, value: str) -> None:
    argv = run["argv"]
    if name in argv:
        argv[argv.index(name) + 1] = value
    else:
        argv[2:2] = [name, value]


@pytest.mark.parametrize(
    "what", ["commit", "profile", "memory", "modules"],
)  # fmt: skip
def test_repeats_that_differ_in_inputs_are_inputs_differ(bs: ModuleType, what: str) -> None:
    """Valid repeats run on another commit, profile, task memory or module set are not one
    measurement (finding 8; lead: reason inputs_differ)."""
    runs = _three([1500, 1520, 1540], variant="core")
    odd = runs[2]
    if what == "commit":
        odd["commit"] = "2" * 40
    elif what == "profile":
        _set_flag(odd, "--profile", "farmshare_local,conda")
    elif what == "memory":
        _set_flag(odd, "--memory", "48.GB")
    else:
        _set_flag(odd, "--modules", "core,bgc")
    cell = _cell(bs, runs, "core")
    assert cell["status"] == "incomplete" and cell["median_s"] is None, cell
    assert "inputs_differ" in cell["reasons"], cell["reasons"]


def test_repeats_differing_only_in_run_paths_and_jobs_are_complete(bs: ModuleType) -> None:
    runs = _three([1500, 1520, 1540], variant="core")
    assert len({json.dumps(r["argv"]) for r in runs}) == 3, "outdirs differ per repeat"
    assert _cell(bs, runs, "core")["status"] == "complete"


def test_not_valid_names_the_latest_attempt(bs: ModuleType) -> None:
    """A repeat without a valid attempt is reported by its latest attempt (finding 12, M14)."""
    runs = [bench_run(repeat=1, wall_s=1000), bench_run(repeat=3, wall_s=1020),
            bench_run(repeat=2, wall_s=900, status="failed", exit_code=4, n_failed=1),
            bench_run(repeat=2, attempt="b", wall_s=950, status="invalid",
                      invalid_reasons=["cached_tasks:3"])]  # fmt: skip
    cell = _cell(bs, runs)
    assert [nv["run"] for nv in cell["not_valid"]] == ["r2.b"], cell["not_valid"]
    assert cell["not_valid"][0]["status"] == "invalid"
    md = _md(bs, [cell])
    assert "r2.b invalid" in md and "cached_tasks:3" in md, md
    assert "r2 failed" not in md, md


def test_cli_write_sets_the_adr(tmp_path: Path, baseline: Path) -> None:
    """``--write --adr ADR-NNNN`` fills ``adr``; a same-commit rewrite without ``--adr`` keeps it
    (finding 9)."""
    d = _complete_bench(tmp_path)
    res = _cli(d, "--baseline", baseline, "--write", "--adr", "ADR-0061")
    assert res.returncode == 0, res.stderr
    wheat = json.loads(baseline.read_text())["standard_run"]["wheat"]
    assert {e["adr"] for e in wheat.values()} == {"ADR-0061"}, wheat
    res = _cli(d, "--baseline", baseline, "--write")
    assert res.returncode == 0, res.stderr
    wheat = json.loads(baseline.read_text())["standard_run"]["wheat"]
    assert {e["adr"] for e in wheat.values()} == {"ADR-0061"}, "kept on a same-commit rewrite"


def test_cli_write_without_adr_leaves_it_null(tmp_path: Path, baseline: Path) -> None:
    assert _cli(_complete_bench(tmp_path), "--baseline", baseline, "--write").returncode == 0
    wheat = json.loads(baseline.read_text())["standard_run"]["wheat"]
    assert {e["adr"] for e in wheat.values()} == {None}


# ================================================================== properties


@settings(max_examples=60, deadline=None)
@given(st.lists(st.integers(1, 20_000), min_size=3, max_size=3), st.permutations(range(3)))
def test_median_spread_noisy_property(walls: list[int], order: list[int]) -> None:
    bs = _mod()
    runs = _three([float(w) for w in walls])
    cell = _cell(bs, [runs[i] for i in order])
    med = sorted(walls)[1]
    spread = (max(walls) - min(walls)) / med
    assert cell["median_s"] == med
    assert cell["spread"] == pytest.approx(spread, rel=1e-6, abs=1e-6)
    gate = yaml.safe_load(GATES.read_text())["benchmark_slowdown_max"]
    if abs(spread - gate) > 1e-9:
        assert cell["noisy"] is (spread > gate)


_STATUS = st.sampled_from(["ok", "failed", "timeout", "invalid"])


@settings(max_examples=60, deadline=None)
@given(st.lists(_STATUS, min_size=1, max_size=3), st.randoms(use_true_random=False))
def test_attempt_selection_property(statuses: list[str], rnd: Any) -> None:
    """Per repeat, the latest valid attempt wins, whatever the order of the records."""
    bs = _mod()
    letters = "abc"
    r2 = [bench_run(repeat=2, attempt=letters[i], wall_s=1100 + 10 * i, status=s)
          for i, s in enumerate(statuses)]  # fmt: skip
    runs = [bench_run(repeat=1, wall_s=1000), bench_run(repeat=3, wall_s=1200), *r2]
    rnd.shuffle(runs)
    cell = _cell(bs, runs)
    ok = [i for i, s in enumerate(statuses) if s == "ok"]
    if not ok:
        assert cell["status"] == "incomplete" and cell["median_s"] is None
    else:
        assert cell["status"] == "complete"
        assert bs.baseline_entry(cell)["runs_s"] == [1000, 1100 + 10 * ok[-1], 1200]


# ================================================================== CLI


def _complete_bench(tmp: Path, commit: str = "1f567901" + "0" * 32) -> Path:
    d = tmp / "bench" / "2026-10-30_1f56790"
    walls = {"k12_standard": [3912, 3700, 4106], "k12_core": [1500, 1520, 1800],
             "sim5m_standard": [5000, 5100, 5200], "sim5m_core": [1490, 1505, 1510]}  # fmt: skip
    runs = []
    for cell, ws in walls.items():
        g, v = cell.split("_")
        runs += [bench_run(g, v, k + 1, wall_s=w, commit=commit) for k, w in enumerate(ws)]
    write_runs(d, runs)
    return d


def _cli(*args: str | Path) -> subprocess.CompletedProcess[str]:
    _mod()
    return subprocess.run([sys.executable, str(SUMMARY), *map(str, args)], capture_output=True,
                          text=True, timeout=120)  # fmt: skip


def _cells(summary: dict[str, Any]) -> dict[str, dict[str, Any]]:
    cells = summary["cells"]
    if isinstance(cells, dict):
        return cells
    return {f"{c['genome']}_{c['variant']}": c for c in cells}


@pytest.fixture
def baseline(tmp_path: Path) -> Path:
    p = tmp_path / "baseline.json"
    shutil.copyfile(BASELINE, p)
    return p


def test_cli_writes_summary_files(tmp_path: Path) -> None:
    d = _complete_bench(tmp_path)
    res = _cli(d)
    assert res.returncode == 0, res.stderr
    summary = json.loads((d / "summary.json").read_text())
    cells = _cells(summary)
    assert set(cells) == set(CELLS)
    assert all(c["status"] == "complete" for c in cells.values())
    md = (d / "summary.md").read_text()
    assert "gippa run standard" in md and "gippa run core" in md, "PROGRESS rows"
    assert "noisy (spread 19.7 % > 15 %)" in md


def test_cli_cells_in_fixed_order(tmp_path: Path) -> None:
    d = _complete_bench(tmp_path)
    assert _cli(d).returncode == 0
    summary = json.loads((d / "summary.json").read_text())
    assert list(_cells(summary)) == list(CELLS)


def test_cli_without_write_leaves_the_baseline(tmp_path: Path, baseline: Path) -> None:
    before = baseline.read_bytes()
    assert _cli(_complete_bench(tmp_path), "--baseline", baseline).returncode == 0
    assert baseline.read_bytes() == before


def test_cli_write_adds_complete_cells_only(tmp_path: Path, baseline: Path) -> None:
    d = _complete_bench(tmp_path)
    (d / "r3" / "sim5m_core" / "bench_run.json").unlink()  # sim5m/core: 2 of 3
    before_text = baseline.read_text()
    before = json.loads(before_text)
    res = _cli(d, "--baseline", baseline, "--write")
    assert res.returncode == 0, res.stderr
    text = baseline.read_text()
    after = json.loads(text)
    assert set(after["standard_run"]) == {"wheat"}
    assert set(after["standard_run"]["wheat"]) == {"k12_standard", "k12_core", "sim5m_standard"}
    for entry in after["standard_run"]["wheat"].values():
        assert set(entry) == BASELINE_KEYS
    assert after["standard_run"]["wheat"]["k12_standard"]["median_s"] == 3912
    rest = {k: v for k, v in after.items() if k != "standard_run"}
    assert json.dumps(rest, indent=2, ensure_ascii=False) + "\n" == before_text, (
        "every other key byte-for-byte, in its order"
    )
    for k in before:
        assert json.dumps(after[k], ensure_ascii=False) == json.dumps(before[k], ensure_ascii=False)
    assert text == json.dumps(after, indent=2, ensure_ascii=False) + "\n", "the file's own format"


def test_cli_refuses_another_commit_without_replace(tmp_path: Path, baseline: Path) -> None:
    assert (
        _cli(_complete_bench(tmp_path / "one"), "--baseline", baseline, "--write").returncode == 0
    )
    written = baseline.read_bytes()
    other = _complete_bench(tmp_path / "two", commit="2" * 40)
    res = _cli(other, "--baseline", baseline, "--write")
    assert res.returncode == 1, res.stdout + res.stderr
    assert baseline.read_bytes() == written, "file unchanged"
    res = _cli(other, "--baseline", baseline, "--write", "--replace")
    assert res.returncode == 0, res.stderr
    entry = json.loads(baseline.read_text())["standard_run"]["wheat"]["k12_standard"]
    assert str(entry["commit"]).startswith("2222222"), entry["commit"]


def test_cli_same_commit_rewrite_is_allowed(tmp_path: Path, baseline: Path) -> None:
    d = _complete_bench(tmp_path)
    assert _cli(d, "--baseline", baseline, "--write").returncode == 0
    first = baseline.read_bytes()
    assert _cli(d, "--baseline", baseline, "--write").returncode == 0
    assert baseline.read_bytes() == first, "idempotent"


def test_cli_no_data(tmp_path: Path, baseline: Path) -> None:
    d = tmp_path / "bench" / "2026-10-30_1f56790"
    d.mkdir(parents=True)
    before = baseline.read_bytes()
    res = _cli(d, "--baseline", baseline, "--write")
    assert res.returncode == 0, res.stderr
    cells = _cells(json.loads((d / "summary.json").read_text()))
    assert set(cells) == set(CELLS) and {c["status"] for c in cells.values()} == {"no_data"}
    assert baseline.read_bytes() == before, "no baseline write"


@pytest.mark.parametrize("damage", ["not_json", "wrong_schema", "missing_key"])
def test_cli_schema_invalid_record_exits_2(tmp_path: Path, damage: str) -> None:
    d = _complete_bench(tmp_path)
    p = d / "r1" / "k12_core" / "bench_run.json"
    rec = json.loads(p.read_text())
    if damage == "not_json":
        p.write_text("{not json")
    elif damage == "wrong_schema":
        rec["schema"] = "gippa-bench-run/0"
        p.write_text(json.dumps(rec))
    else:
        del rec["wall_s"]
        p.write_text(json.dumps(rec))
    assert _cli(d).returncode == 2
