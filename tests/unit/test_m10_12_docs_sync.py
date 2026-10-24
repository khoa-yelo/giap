"""User-guide runtime table vs ``benchmarks/baseline.json`` (M10-12 criterion 10).

docs/specs/M10-12-benchmark-standard.md, lead decision L9: ``docs/user-guide/runtimes.md`` holds
the ``k12`` rows (variant, median h:mm, range of the 3 runs, peak memory GiB, node class, date,
commit). Skips while ``standard_run`` is absent from ``baseline.json`` (the evidence PR adds
both); from then on every ``k12`` entry must have a matching row.
"""

from __future__ import annotations

import json
from typing import Any

import pytest
from helpers.bench_m10_12 import BASELINE
from helpers.slurm_sandbox import REPO

GUIDE = REPO / "docs" / "user-guide" / "runtimes.md"


def _hmm(seconds: float) -> set[str]:
    """``h:mm`` of a median, floor or rounded to the minute."""
    out = set()
    for m in (int(seconds // 60), round(seconds / 60)):
        out.add(f"{m // 60}:{m % 60:02d}")
    return out


def _entries() -> dict[str, Any]:
    data = json.loads(BASELINE.read_text())
    if "standard_run" not in data:
        pytest.skip("benchmarks/baseline.json has no standard_run yet (M10-12 evidence PR)")
    return data["standard_run"]


def test_guide_table_matches_the_baseline_entries() -> None:
    classes = _entries()
    assert GUIDE.is_file(), "docs/user-guide/runtimes.md (lead decision L9)"
    rows = [ln for ln in GUIDE.read_text().splitlines() if ln.lstrip().startswith("|")]
    checked = 0
    for node_class, cells in classes.items():
        for cell, entry in cells.items():
            if not cell.startswith("k12_"):
                continue
            variant = cell.split("_", 1)[1]
            commit = str(entry["commit"])[:7]
            hits = [r for r in rows if variant in r and node_class in r and commit in r]
            assert hits, f"no guide row for {node_class}/{cell} at {commit}"
            median = float(entry["median_s"])
            assert any(any(h in r for h in _hmm(median)) for r in hits), (cell, hits)
            checked += 1
    assert checked, "standard_run has no k12 entry"


def test_guide_shows_no_sim5m_rows() -> None:
    _entries()
    text = GUIDE.read_text() if GUIDE.is_file() else ""
    rows = [ln for ln in text.splitlines() if ln.lstrip().startswith("|")]
    assert not [r for r in rows if "sim5m" in r.lower()], "the guide shows only k12 (L8)"


def test_guide_peak_memory_is_the_step_max_rss() -> None:
    """How much memory a user needs is the run's sacct MaxRSS (JVM, head and concurrent tasks),
    not the largest single task (PR #300 review finding 10). ``guide_rows`` moved to the docs
    evidence PR with this test; it skips until then."""
    classes = _entries()
    assert GUIDE.is_file(), "docs/user-guide/runtimes.md (lead decision L9)"
    rows = [ln for ln in GUIDE.read_text().splitlines() if ln.lstrip().startswith("|")]
    for node_class, cells in classes.items():
        for cell, entry in cells.items():
            if not cell.startswith("k12_") or entry.get("peak_rss_sacct_max_bytes") is None:
                continue
            variant, commit = cell.split("_", 1)[1], str(entry["commit"])[:7]
            hits = [r for r in rows if variant in r and node_class in r and commit in r]
            gib = f"{entry['peak_rss_sacct_max_bytes'] / 1024**3:.1f}"
            assert any(gib in r for r in hits), (cell, gib, hits)
