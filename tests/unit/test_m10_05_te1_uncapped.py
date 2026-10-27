"""``validation/proxies/te_proxies.te1_result(..., max_exceptions=...)`` (M10-05 review T1).

The pilot's TE1 needs every miss row (``te1_exceptions.tsv`` is uncapped, spec "Outputs"); the
proxy keeps its 50-row cap by default. Lead decision (M10-05 review round): a keyword
``max_exceptions: int | None = MAX_EXCEPTIONS`` on ``te1_result``; ``None`` means no cap. The
counts never depend on the cap. The proxy's own tests (tests/unit/test_te_proxies.py, on main) are
untouched; this file covers only the new keyword.
"""

from __future__ import annotations

from typing import Any

import pandas as pd
import pytest
from helpers.pilot_metrics_m10 import load_te_proxies

from gippa.model.coords import Interval

HIT_COLS = ("peptide", "route", "locus_class", "gene", "contig", "start", "end", "strand")


@pytest.fixture(scope="module")
def te() -> Any:
    return load_te_proxies()


def _world(n_good: int, n_miss: int) -> tuple[dict[str, set[str]], pd.DataFrame, dict[str, Any]]:
    """``n_good`` peptides exact in frame at G1; ``n_miss`` peptides with no hit at all."""
    sources: dict[str, set[str]] = {}
    rows = []
    for i in range(n_good):
        p = f"GOOD{i:04d}K"
        sources[p] = {"G1"}
        rows.append({"peptide": p, "route": "exact", "locus_class": "cds_inframe", "gene": "G1",
                     "contig": "chr", "start": 1, "end": 30, "strand": "+"})  # fmt: skip
    for i in range(n_miss):
        sources[f"MISS{i:04d}K"] = {"G1"}
    hits = pd.DataFrame(rows, columns=list(HIT_COLS))
    genes = {"G1": Interval("chr", 0, 60000, "+")}
    return sources, hits, genes


def _misses(res: Any) -> list[dict[str, Any]]:
    return [r for r in res.exceptions if r["kind"] == "miss"]


def test_default_keeps_the_50_row_cap(te: Any) -> None:
    sources, hits, genes = _world(10, 60)
    res = te.te1_result(sources, {}, hits, genes, {"chr": None})
    assert te.MAX_EXCEPTIONS == 50
    assert len(res.exceptions) == 50 and res.counts["misses"] == 60


def test_max_exceptions_none_returns_every_miss(te: Any) -> None:
    sources, hits, genes = _world(10, 60)
    res = te.te1_result(sources, {}, hits, genes, {"chr": None}, max_exceptions=None)
    rows = _misses(res)
    assert len(rows) == 60 == res.counts["misses"]
    assert {r["peptide"] for r in rows} == {p for p in sources if p.startswith("MISS")}


def test_max_exceptions_int_caps_and_counts_do_not_change(te: Any) -> None:
    sources, hits, genes = _world(10, 60)
    capped = te.te1_result(sources, {}, hits, genes, {"chr": None}, max_exceptions=3)
    full = te.te1_result(sources, {}, hits, genes, {"chr": None}, max_exceptions=None)
    assert len(capped.exceptions) == 3
    assert capped.counts == full.counts
    assert capped.shares == full.shares and capped.passed == full.passed
    assert capped.exceptions == full.exceptions[:3]  # the cap truncates, it does not reorder
