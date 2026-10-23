"""M8-06a item 1 / example 7: ``detectability`` persists M4-13's TU folds as the optional output
``tu_folds`` (``genetable/tu_folds.parquet``; docs/specs/M8-06-function-step-and-run.md; ADR-0035).

Pinned here:

* ``tu_folds`` is in the step's ``outputs`` and ``optional_outputs``; left unwired, nothing is
  written and the other outputs are unchanged.
* Columns ``tu_id`` (string) and ``fold`` (int64), rows exactly
  ``DetectabilityResult.fold_of_tu``, sorted by natural TU ID, written with
  ``write_parquet_deterministic`` (byte-reproducible).
* A constant model (fewer than 50 detected or undetected eligible genes) gives a header-only file.
* Wiring the output does not move the step's random stream: ``gene_table`` is byte-identical with
  and without it.

The non-constant fixture is the ~600-gene M4-13 synthetic layout of
``helpers.function_step_data.world_pre``.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq  # type: ignore[import-untyped]
import pytest
from helpers.function_step_data import natural_key, toy4, world_pre, write_gene_table

from giap.config import GiapConfig


@pytest.fixture(autouse=True, scope="module")
def _simple_gsea_engine() -> Iterator[None]:
    """The step runs on ``fgsea_simple``, never R: unit tests use fakes (ADR-0038 §4, M8-17b 10).
    Module scope, so module-scoped fixtures that run the step see it too."""
    import giap.steps.function as fn
    from giap.function.gsea import SIMPLE_ENGINE

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(fn, "_LazyRFgsea", lambda log: SIMPLE_ENGINE)
        yield


SEED = 7


def _run(pre: Path, out: Path, *, folds: bool, seed: int = SEED) -> dict[str, Path]:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    outputs = {"gene_table": out / "genetable" / "genes.parquet",
               "detectability_summary": out / "genetable" / "det.parquet"}  # fmt: skip
    if folds:
        outputs["tu_folds"] = out / "genetable" / "tu_folds.parquet"
    reg.run_step("detectability", GiapConfig(seed=seed), {"genes_pre": pre}, outputs)
    return outputs


def _fit(pre: Path, seed: int = SEED):  # type: ignore[no-untyped-def]
    from giap.genetable.covariates import COVARIATE_COLUMNS
    from giap.genetable.detectability import fit_detectability
    from giap.io.genes_tsv import read_genes_parquet
    from giap.rng import for_step

    p = read_genes_parquet(pre)
    idx = pd.Index(p["locus_tag"].astype(str))
    det = pd.Series(p["detected"].astype(bool).to_numpy(), index=idx, dtype="boolean")
    det[p["partial"].astype(bool).to_numpy()] = pd.NA
    cov = p[["locus_tag", *COVARIATE_COLUMNS]].astype(dict.fromkeys(COVARIATE_COLUMNS, "float64"))
    cfg = GiapConfig(seed=seed)
    return fit_detectability(det, cov, pd.Series(p["tu_id"].to_numpy(), index=idx),
                             cfg.detectability, for_step("detectability", seed))  # fmt: skip


@pytest.fixture(scope="module")
def pre600(tmp_path_factory: pytest.TempPathFactory) -> Path:
    root = tmp_path_factory.mktemp("det_folds")
    pre, _ = world_pre()
    return write_gene_table(root / "genes_pre.parquet", pre)


@pytest.fixture(scope="module")
def wired(pre600: Path, tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    return _run(pre600, tmp_path_factory.mktemp("wired"), folds=True)


def test_tu_folds_is_an_optional_output_of_detectability() -> None:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    spec = reg.get("detectability")
    assert set(spec.outputs) == {"gene_table", "detectability_summary", "tu_folds"}
    assert spec.optional_outputs == ("tu_folds",)
    assert spec.inputs == ("genes_pre",)


def test_example_7_rows_equal_fold_of_tu(pre600: Path, wired: dict[str, Path]) -> None:
    fit = _fit(pre600)
    assert not fit.constant, "fixture meant to give a non-constant model"
    table = pq.read_table(wired["tu_folds"])
    assert table.schema.names == ["tu_id", "fold"]
    assert str(table.schema.field("tu_id").type) == "string"
    assert str(table.schema.field("fold").type) == "int64"
    df = table.to_pandas()
    got = dict(zip(df["tu_id"].astype(str), df["fold"].astype(int), strict=True))
    assert got == dict(fit.fold_of_tu)
    assert len(df) == len(fit.fold_of_tu) > 0


def test_example_7_rows_sorted_by_natural_tu_id(wired: dict[str, Path]) -> None:
    ids = pd.read_parquet(wired["tu_folds"])["tu_id"].astype(str).tolist()
    assert ids == sorted(ids, key=natural_key)
    assert ids != sorted(ids), "fixture meant to have TU_2 < TU_10 (natural, not lexical)"


def test_tu_folds_byte_reproducible(pre600: Path, wired: dict[str, Path], tmp_path: Path) -> None:
    again = _run(pre600, tmp_path, folds=True)
    assert again["tu_folds"].read_bytes() == wired["tu_folds"].read_bytes()


def test_wiring_tu_folds_does_not_move_the_stream(
    pre600: Path, wired: dict[str, Path], tmp_path: Path
) -> None:
    plain = _run(pre600, tmp_path, folds=False)
    for key in ("gene_table", "detectability_summary"):
        assert plain[key].read_bytes() == wired[key].read_bytes(), key


def test_unwired_tu_folds_is_not_written(pre600: Path, tmp_path: Path) -> None:
    out = _run(pre600, tmp_path, folds=False)
    assert out["gene_table"].is_file()
    assert not list(tmp_path.rglob("*tu_folds*"))


def test_constant_model_gives_a_header_only_file(tmp_path: Path) -> None:
    """Item 1: ``detectability_constant`` (here 4 genes) -> header-only, so ``function`` sees
    ``{}`` and records ``tu_recal_skipped`` (sim100k's documented outcome)."""
    genes = toy4(length_aa=[300, 200, 250, 150], mw_da=[33000.0, 22000.0, 27500.0, 16500.0],
                 pi=[6.0, 7.0, 5.5, 8.0], gravy=[-0.2, 0.1, -0.4, 0.0],
                 tm_segments=[0, 1, 0, 0], n_observable_peptides=[12, 8, 10, 6])  # fmt: skip
    pre = write_gene_table(tmp_path / "genes_pre.parquet", genes)
    assert _fit(pre).constant
    out = _run(pre, tmp_path / "out", folds=True)
    table = pq.read_table(out["tu_folds"])
    assert table.num_rows == 0
    assert table.schema.names == ["tu_id", "fold"]
    assert str(table.schema.field("tu_id").type) == "string"
    assert str(table.schema.field("fold").type) == "int64"


def test_example_7_folds_feed_function_without_skip(
    pre600: Path, wired: dict[str, Path], tmp_path: Path
) -> None:
    """Example 7 end to end: detectability's own ``tu_folds.parquet`` into ``function`` gives no
    ``tu_recal_skipped``."""
    from helpers.function_db import write_function_db
    from helpers.function_step_data import write_replicons

    from giap.steps import registry as reg

    _, xrefs = world_pre()
    db = write_function_db(tmp_path / "db", xrefs)
    ins = {"gene_table": wired["gene_table"], "tu_folds": wired["tu_folds"],
           "replicons": write_replicons(tmp_path / "replicons.parquet")}  # fmt: skip
    ins |= db.inputs(["cog2020_fun"])
    outs = {k: tmp_path / "tables" / f"{k}.tsv"
            for k in ("category_summary", "enrichment_detection", "enrichment_gsea")}  # fmt: skip
    cfg = GiapConfig.model_validate({"seed": SEED, "stats": {"bootstrap": 50},
                                     "function": {"gsea": {"n_perm_simple": 100}}})  # fmt: skip
    import logging

    seen: list[logging.LogRecord] = []

    class H(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            seen.append(record)

    log = logging.getLogger("giap.steps.function")
    h, old = H(level=logging.DEBUG), log.level
    log.addHandler(h)
    log.setLevel(logging.DEBUG)
    try:
        reg.run_step("function", cfg, ins, outs)
    finally:
        log.removeHandler(h)
        log.setLevel(old)
    events = [r for r in seen if getattr(r, "event", None) == "function.run_flags"]
    assert len(events) == 1
    assert "tu_recal_skipped" not in events[0].data["flags"]  # type: ignore[attr-defined]
