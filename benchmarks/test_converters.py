"""Performance budget for the DIA-NN converter: 500k-row, 60+-column parquet report <= 15 s.

Moved from tests/unit/test_convert_review.py (which keeps a loose 10x guard). Not part of the
per-commit suite; run with ``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ._unit_helpers import load_test_module, timed

_T = load_test_module("unit/test_convert_review.py")


def test_diann_500k_parquet_budget(benchmark: Any, tmp_path: Path) -> None:
    path = tmp_path / "report.parquet"
    _T._big_report(path)

    def run() -> float:
        secs, res = timed(lambda: _T._diann([path]))
        assert not _T._errors(res.issues)
        assert len(res.peptides) == 500_000
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs < 15.0, f"conversion took {secs:.1f} s"
