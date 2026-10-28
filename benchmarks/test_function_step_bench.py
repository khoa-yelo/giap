"""Runtime of the ``function`` step on ``sim4500`` (docs/specs/M8-06-function-step-and-run.md,
M8-06a acceptance "Benchmark ``function_step_4500_genes``").

One ``build_m4_fixture(size="sim4500")`` work directory, the synthetic function DB of
``helpers.function_db`` built from its own C3 ``xrefs`` (lead decision 1: no real GO or KEGG
content), ``gene_table`` and ``detectability`` (with ``tu_folds`` wired) are prepared once and not
timed. Only ``run_step("function", ...)`` is timed, with every resource installed and the shipped
``stats.bootstrap`` = 1,000. The measured time belongs under ``function_step_4500_genes`` in
``benchmarks/baseline.json``. The budget follows the M8 index feasibility line ("seconds to
minutes"); it is about 3 x the recorded time.

Run with ``pytest benchmarks/test_function_step_bench.py --benchmark-only`` (needs ``tests/sim``
on PYTHONPATH for gippasim).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from ._unit_helpers import timed

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tests" / "helpers"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tests"))
import m4_chain as mc
from helpers.function_db import write_function_db

SEED = 20261022
BUDGET_S = 45.0  # about 3 x the measured 14.1 s (numba compilation included)
TABLES = ("category_summary", "enrichment_detection", "enrichment_gsea")


@pytest.fixture(scope="module")
def prepared(tmp_path_factory: pytest.TempPathFactory) -> tuple[dict[str, Path], Path, Any]:
    from gippa.steps import load_step_modules
    from gippa.steps.gene_table import run_gene_table
    from gippa.steps.registry import run_step

    load_step_modules()
    root = tmp_path_factory.mktemp("function_bench")
    w = mc.build(root / "work", seed=SEED, size="sim4500")
    db = write_function_db(root / "db", pd.read_parquet(w / "xrefs.parquet"))
    cfg = mc.config(SEED)
    gt = w / "genetable"
    with mc.db_env(db.root):
        run_gene_table(w, cfg)
    run_step("detectability", cfg, {"genes_pre": gt / "genes_pre.parquet"},
             {"gene_table": gt / "genes.parquet",
              "detectability_summary": gt / "detection_summary_parts" / "detectability.parquet",
              "tu_folds": gt / "tu_folds.parquet"})  # fmt: skip
    ins = {"gene_table": gt / "genes.parquet", "replicons": w / "replicons.parquet",
           "tu_folds": gt / "tu_folds.parquet"} | db.inputs()  # fmt: skip
    return ins, w / "tables", cfg


def test_function_step_4500_genes_budget(
    benchmark: Any, prepared: tuple[dict[str, Path], Path, Any]
) -> None:
    from gippa.steps.registry import run_step

    ins, out, cfg = prepared
    assert cfg.stats.bootstrap == 1000
    outs = {k: out / f"{k}.tsv" for k in TABLES}

    def run() -> float:
        secs, _ = timed(lambda: run_step("function", cfg, ins, outs))
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    rows = {k: len(outs[k].read_text(encoding="utf-8").splitlines()) - 1 for k in TABLES}
    benchmark.extra_info["rows"] = rows
    assert rows["category_summary"] > 0 and rows["enrichment_detection"] > 0, rows
    assert secs <= BUDGET_S, f"function on sim4500 took {secs:.1f} s (budget {BUDGET_S} s)"
