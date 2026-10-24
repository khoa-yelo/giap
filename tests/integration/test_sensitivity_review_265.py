"""M8-16 review round at step level (lead decision 5; code review 6372795283 B2a, B2b, N3):

* B2a: the step reads ``tu_folds`` and hands them to the ``enriched_categories`` ORA as
  ``folds`` (``None`` without the file), as the ``function`` step does;
* B2b: ``operon_flags`` and ``determinants`` both name a ``features`` keyword (C3 vs M8-07's TSV);
  in one step call each family gets its own table, and the operon rows equal the run without
  ``determinants_features``;
* N3: a determinants variant below the 100-row / 25-TU minimum writes no row and logs the event
  ``sensitivity.determinants_variant_skipped`` naming the variant; the step does not crash and the
  other variants keep their rows.

Families are replaced (monkeypatched ``FAMILIES``) by spies around the real ones or by small
fakes, so each test pins the step's binding and bookkeeping, not the statistics.
"""

from __future__ import annotations

import inspect
import logging
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import determinants_models_data as dm
from helpers import sensitivity_data as sd
from helpers import sensitivity_function_data as sf
from helpers import sensitivity_operons_data as so
from helpers.layers_chain import Chain, build_chain

RESOURCES = ("go_basic", "goslim_generic", "kegg_ko_pathway", "kegg_pathway_names", "cog2020_fun")
REQUIRED = ("replicons", "regions", "region_genes", "tus")
VARIANT_EVENT = "sensitivity.determinants_variant_skipped"


def _S() -> Any:
    from giap.context import sensitivity as S

    return S


class _Records(logging.Handler):
    def __init__(self) -> None:
        super().__init__(logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


def _run_step(
    cfg: Any, ins: dict[str, Path], out: Path, families: dict[str, Any]
) -> tuple[pd.DataFrame, list[logging.LogRecord]]:
    from giap.steps import load_step_modules
    from giap.steps.registry import run_step

    load_step_modules()
    h = _Records()
    root = logging.getLogger("giap")
    old = root.level
    root.addHandler(h)
    root.setLevel(logging.DEBUG)
    try:
        with pytest.MonkeyPatch.context() as mp:
            mp.setattr(_S(), "FAMILIES", families)
            run_step("sensitivity", cfg, ins, {"sensitivity": out})
    finally:
        root.removeHandler(h)
        root.setLevel(old)
    return sd.read_sensitivity(out), h.records


def _spy(family: Callable[..., Any], store: list[dict[str, Any]]) -> Callable[..., Any]:
    """``family`` recording the keyword tables of each call; its signature kept so the step's
    ``_bind`` binds the same tables."""

    def wrapped(genes: pd.DataFrame, cfg: Any, rng: np.random.Generator, **kw: Any) -> Any:
        store.append(dict(kw))
        return family(genes, cfg, rng, **kw)

    wrapped.__signature__ = inspect.signature(family)  # type: ignore[attr-defined]
    return wrapped


def _det_features_tsv(genes: pd.DataFrame, path: Path) -> Path:
    """An M8-07 features TSV for ``genes`` (noise feature columns)."""
    from giap.io.determinants_tables import write_determinants_features

    g = genes[genes["tu_id"].notna()].reset_index(drop=True)
    n = len(g)
    r = np.random.default_rng(3)
    cols = {c: r.normal(0.0, 1.0, n).tolist() for c in dm.PLANTED_BLOCKS}
    ab = pd.to_numeric(g["abundance_log2"], errors="coerce").astype(float).tolist()
    feats = dm.features(cols, tu=g["tu_id"].astype(str).tolist(),
                        tags=g["locus_tag"].astype(str).tolist(),
                        detected=g["detected"].astype(bool).tolist(), abundance=ab,
                        imputed=g["abundance_imputed"].astype(bool).tolist())  # fmt: skip
    path.parent.mkdir(parents=True, exist_ok=True)
    write_determinants_features(path, feats)
    return path


# --------------------------------------------------------------------------- B2a: ORA folds


@dataclass
class World:
    root: Path
    inputs: dict[str, Path]
    genes: pd.DataFrame


def _function_columns(g: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    g = g.copy()
    n = len(g)
    g["cog_category"] = [("J", "K", "L")[i % 3] for i in range(n)]
    g["cog_id"] = [f"COG{i % 3:04d}" for i in range(n)]
    g["go_terms"] = [[f"GO:{i % 12 + 1:07d}"] for i in range(n)]
    g["kegg_ko"] = [[f"K{i % 5 + 1:05d}"] for i in range(n)]
    rows = [(t, db, x) for t, gos, kos, c in zip(g["locus_tag"], g["go_terms"], g["kegg_ko"],
                                                 g["cog_id"], strict=True)
            for db, x in [("GO", gos[0]), ("KEGG", kos[0]), ("COG", c)]]  # fmt: skip
    return g, pd.DataFrame(rows, columns=["locus_tag", "db", "id"])


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> World:
    from helpers.function_db import write_function_db
    from helpers.function_step_data import write_tu_folds

    from giap.io.genes_tsv import read_genes_parquet, write_genes_parquet

    root = tmp_path_factory.mktemp("sens_review")
    work = sd.sim_work_dir(root / "work")
    g, xrefs = _function_columns(read_genes_parquet(work / "genetable" / "genes.parquet"))
    gt = root / "in" / "genes.parquet"
    gt.parent.mkdir(parents=True, exist_ok=True)
    write_genes_parquet(gt, g)
    db = write_function_db(root / "db", xrefs)
    tus = sorted(set(g["tu_id"].dropna().astype(str)))
    folds = write_tu_folds(root / "in" / "tu_folds.parquet", {t: i % 5 for i, t in enumerate(tus)})
    ins = {"gene_table": gt, **{k: work / f"{k}.parquet" for k in REQUIRED},
           **db.inputs(RESOURCES), "tu_folds": folds}  # fmt: skip
    return World(root, ins, read_genes_parquet(gt))


def _ora_folds(w: World, out: Path, *, drop: tuple[str, ...] = ()) -> list[Any]:
    import giap.function.sensitivity as FS
    from giap.function.gsea import SIMPLE_ENGINE

    seen: list[Any] = []
    real = FS.detection_enrichment

    def spy(*a: Any, **k: Any) -> Any:
        seen.append(k.get("folds"))
        return real(*a, **k)

    fam = sf.with_engine(FS.enriched_categories_family, SIMPLE_ENGINE)
    ins = {k: v for k, v in w.inputs.items() if k not in drop}
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(FS, "detection_enrichment", spy)
        rows, _ = _run_step(sf.config(), ins, out, {sf.EC: fam})
    assert set(rows["family"]) == {sf.EC}
    return seen


def test_step_passes_tu_folds_to_the_ora(world: World, tmp_path: Path) -> None:
    """B2a: every ORA call (primary and each variant) gets the ``tu_folds`` file's folds."""
    from giap.steps.function import read_folds

    want = read_folds(world.inputs["tu_folds"])
    seen = _ora_folds(world, tmp_path / "s.tsv")
    assert len(seen) == 3, f"primary + missing=min + missing=exclude expected: {len(seen)}"
    assert all(s is not None and dict(s) == want for s in seen), "ORA ran without the TU folds"


def test_step_without_tu_folds_runs_the_ora_uncalibrated(world: World, tmp_path: Path) -> None:
    seen = _ora_folds(world, tmp_path / "s.tsv", drop=("tu_folds",))
    assert seen and all(s is None for s in seen)


# --------------------------------------------------------------------------- B2b: features


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return build_chain(Path(tmp_path_factory.mktemp("sens_review_operons")))


def _fake_det(store: list[dict[str, Any]]) -> Callable[..., Any]:
    from giap.determinants.sensitivity import DeterminantsSensitivity

    def determinants(genes: pd.DataFrame, cfg: Any, rng: np.random.Generator, *,
                     features: pd.DataFrame) -> Any:  # fmt: skip
        store.append({"features": features})
        return DeterminantsSensitivity(frozenset({"g1"}), {"g1": "codon"}, 0.5)

    return determinants


def test_operon_flags_and_determinants_each_get_their_features(chain: Chain,
                                                              tmp_path: Path) -> None:  # fmt: skip
    """B2b: one step call with pairs_raw + C2 / C3 and ``determinants_features``: operon_flags
    gets C3, determinants M8-07's TSV, both families write rows, and the operon rows equal the
    run without ``determinants_features``."""
    from giap.io.determinants_tables import read_determinants_features

    S = _S()
    cfg = so.config()
    gt = so.gene_table(so.tables(chain)["genes"], imputed=(1, 2, 3))
    ins = so.write_inputs(chain, tmp_path / "in", gt)
    det_tsv = _det_features_tsv(gt, tmp_path / "in" / "determinants_features.tsv")
    op_calls: list[dict[str, Any]] = []
    det_calls: list[dict[str, Any]] = []
    fams = {so.FAMILY: _spy(S.operon_flags_family, op_calls), sf.DET: _fake_det(det_calls)}
    both, _ = _run_step(cfg, {**ins, "determinants_features": det_tsv}, tmp_path / "a.tsv", fams)

    c3 = pd.read_parquet(ins["features"])
    m807 = read_determinants_features(det_tsv)
    assert op_calls and det_calls
    for kw in op_calls:
        pd.testing.assert_frame_equal(kw["features"].reset_index(drop=True), c3)
    for kw in det_calls:
        pd.testing.assert_frame_equal(kw["features"].reset_index(drop=True), m807)
    assert set(both["family"]) == {so.FAMILY, sf.DET}

    alone, _ = _run_step(cfg, ins, tmp_path / "b.tsv",
                         {so.FAMILY: S.operon_flags_family, sf.DET: _fake_det([])})  # fmt: skip
    assert set(alone["family"]) == {so.FAMILY}
    op = both[both["family"] == so.FAMILY].reset_index(drop=True)
    pd.testing.assert_frame_equal(op, alone.reset_index(drop=True))


# --------------------------------------------------------------------------- N3: variant skip


def test_variant_below_the_minimum_gets_no_row_and_an_event(chain: Chain, tmp_path: Path) -> None:
    """N3: a determinants variant whose rerun is below the minimum (the fake family returns
    ``None`` under ``missing=exclude``) writes no row and logs
    ``sensitivity.determinants_variant_skipped`` naming the variant; ``missing=min`` keeps its
    row and the step does not crash."""
    from giap.determinants.sensitivity import DeterminantsSensitivity

    def determinants(genes: pd.DataFrame, cfg: Any, rng: np.random.Generator, *,
                     features: pd.DataFrame) -> Any:  # fmt: skip
        if cfg.abundance.missing == "exclude":
            return None
        return DeterminantsSensitivity(frozenset({"g1"}), {"g1": "codon"}, 0.5)

    cfg = so.config(abundance={"missing": "median", "sensitivity": ["min", "exclude"]})
    gt = so.gene_table(so.tables(chain)["genes"])
    ins = so.write_inputs(chain, tmp_path / "in", gt, drop=("genes", "features"))
    ins.pop("pairs_raw")
    det_tsv = _det_features_tsv(gt, tmp_path / "in" / "determinants_features.tsv")
    rows, records = _run_step(cfg, {**ins, "determinants_features": det_tsv},
                              tmp_path / "s.tsv", {sf.DET: determinants})  # fmt: skip
    assert list(zip(rows["family"], rows["variant"], strict=True)) == [(sf.DET, "missing=min")]
    hits = [r for r in records if getattr(r, "event", None) == VARIANT_EVENT]
    assert len(hits) == 1, [getattr(r, "event", None) for r in records]
    text = hits[0].getMessage() + " " + repr(getattr(hits[0], "data", None))
    assert "missing=exclude" in text, text
    assert "missing=min" not in text, text
