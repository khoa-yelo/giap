"""Performance budget for peptide-table ingest: 200k rows read+validate+aggregate <= 10 s.

Moved from tests/unit/test_peptide_table_review.py (which keeps a loose 10x guard). Not part of
the per-commit suite; run with ``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from giap.config.model import InputConfig
from giap.io.peptides import aggregate_peptides, read_peptide_table, validate_peptides

from ._unit_helpers import load_test_module, timed

_T = load_test_module("unit/test_peptide_table_review.py")


def test_peptide_table_200k_rows_budget(benchmark: Any, tmp_path: Path) -> None:
    path = tmp_path / "big.tsv"
    _T._synthetic_table(path)

    def body() -> Any:
        df_raw, _ = read_peptide_table(path)
        df, _issues = validate_peptides(df_raw, InputConfig())
        return df, aggregate_peptides(df)

    def run() -> float:
        secs, (df, _agg) = timed(body)
        assert len(df) == 200_000
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs < 10.0, f"read+validate+aggregate took {secs:.1f} s"
