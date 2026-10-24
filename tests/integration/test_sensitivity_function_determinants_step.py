"""M8-16 item 3 at step level: the registered ``sensitivity`` step with the ``enriched_categories``
and ``determinants`` families on the M4-19 giapsim work directory (docs/specs/M8-16-sensitivity-
function-determinants.md items 1-3, degenerate cases, acceptance "a run with function and
determinants on writes the two families"; GIAP_SPEC.md §8.6, §6.2 ``sensitivity`` row).

A new file rather than an extension of ``test_sensitivity_step.py`` (the test author never edits
an existing test). Interface assumed (test author):

* the step declares the five ``FUNCTION_RESOURCES`` keys and ``determinants_features``
  (``tables/determinants_features.tsv``, M8-10's writer) as optional inputs; its required inputs
  are unchanged; it builds the collections from the resources (``read_resources`` +
  ``build_collections`` on the primary gene table) and reads the features TSV;
* the families are called with the step's family streams ``for_step("sensitivity/<family>",
  seed)`` and the base config (``stats.permutations = 199``), as M4-19's;
* a family whose inputs are absent (no resource, ``function`` not in ``modules``, no features
  table, fewer than 100 model rows) writes no row and logs a record with an ``event`` attribute
  whose message names the family; the determinants CV R² of every run goes to an event whose
  ``event`` name contains ``determinants`` and whose ``data`` names every variant;
* the GSEA engine is injected by wrapping the registered family (``with_engine``:
  ``SIMPLE_ENGINE``, recorded), so no R is needed.

Fixture: ``sensitivity_data.sim_work_dir`` (1,000 genes, 60 imputed) with synthetic function
columns (the 40 most abundant observed genes in COG J, the rest in K / L; GO / KEGG IDs by
index, the M8-06 run's scheme), the ``helpers.function_db`` DB for those IDs and a features
table whose ``cai`` tracks abundance. Only the two M8-16 families are left in ``FAMILIES`` (the
others have their own tests).
"""

from __future__ import annotations

import json
import logging
import re
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import determinants_models_data as dm
from helpers import sensitivity_data as sd
from helpers import sensitivity_function_data as sf

REPO = Path(__file__).resolve().parents[2]
WORKFLOW = REPO / "workflow"
RESOURCES = ("go_basic", "goslim_generic", "kegg_ko_pathway", "kegg_pathway_names", "cog2020_fun")
REQUIRED = {"gene_table", "replicons", "regions", "region_genes", "tus"}


def _spec() -> Any:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    return reg.get("sensitivity")


# --------------------------------------------------------------------------- fixture


@dataclass
class World:
    root: Path
    inputs: dict[str, Path]
    genes: pd.DataFrame
    features: pd.DataFrame


def _function_columns(g: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    g = g.copy()
    n = len(g)
    ab = g["abundance_log2"].astype(float).to_numpy()
    obs = g["detected"].astype(bool).to_numpy() & ~g["abundance_imputed"].astype(bool).to_numpy()
    order = np.argsort(-np.where(obs & np.isfinite(ab), ab, -np.inf), kind="stable")
    top = set(order[:40].tolist())
    letters = ("K", "L")  # 3 COG terms: BH keeps J's floor p (~1/(B+1)) under 0.05
    g["cog_category"] = ["J" if i in top else letters[i % 2] for i in range(n)]
    g["cog_id"] = [f"COG{i % 3:04d}" for i in range(n)]
    g["go_terms"] = [[f"GO:{i % 12 + 1:07d}"] for i in range(n)]
    g["kegg_ko"] = [[f"K{i % 5 + 1:05d}"] for i in range(n)]
    rows = [(t, db, x) for t, gos, kos, c in zip(g["locus_tag"], g["go_terms"], g["kegg_ko"],
                                                 g["cog_id"], strict=True)
            for db, x in [("GO", gos[0]), ("KEGG", kos[0]), ("COG", c)]]  # fmt: skip
    return g, pd.DataFrame(rows, columns=["locus_tag", "db", "id"])


def _features(g: pd.DataFrame, keep: int | None = None) -> pd.DataFrame:
    el = g[~g["partial"].astype(bool) & g["tu_id"].notna()].reset_index(drop=True)
    if keep is not None:
        el = el.iloc[:keep]
    n = len(el)
    r = np.random.default_rng(17)
    ab = el["abundance_log2"].astype(float).to_numpy()
    cai = np.where(np.isfinite(ab), ab, 0.0) + r.normal(0.0, 1.0, n)
    cols = {c: r.normal(0.0, 1.0, n).tolist() for c in dm.PLANTED_BLOCKS}
    cols["cai"] = cai.tolist()
    return dm.features(
        cols, tu=el["tu_id"].astype(str).tolist(), tags=el["locus_tag"].astype(str).tolist(),
        detected=el["detected"].astype(bool).tolist(), abundance=ab.tolist(),
        imputed=el["abundance_imputed"].astype(bool).tolist(),
    )  # fmt: skip


def _build(root: Path, **sim: Any) -> World:
    from helpers.function_db import write_function_db

    from giap.io.determinants_tables import write_determinants_features
    from giap.io.genes_tsv import read_genes_parquet, write_genes_parquet

    work = sd.sim_work_dir(root / "work", **sim)
    g, xrefs = _function_columns(read_genes_parquet(work / "genetable" / "genes.parquet"))
    gt = root / "in" / "genes.parquet"
    gt.parent.mkdir(parents=True, exist_ok=True)
    write_genes_parquet(gt, g)
    db = write_function_db(root / "db", xrefs)
    feats = _features(g)
    ft = root / "in" / "determinants_features.tsv"
    write_determinants_features(ft, feats)
    ins = {"gene_table": gt, **{k: work / f"{k}.parquet" for k in REQUIRED - {"gene_table"}},
           **db.inputs(RESOURCES), "determinants_features": ft}  # fmt: skip
    return World(root, ins, read_genes_parquet(gt), feats)


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> World:
    return _build(tmp_path_factory.mktemp("sens_m8"))


@pytest.fixture(scope="module")
def world_2pep(tmp_path_factory: pytest.TempPathFactory) -> World:
    return _build(tmp_path_factory.mktemp("sens_m8_2pep"),
                  evidence={"emit": frozenset({"intensity"})})  # fmt: skip


class _Engine:
    """``SIMPLE_ENGINE`` recording each call's parameters."""

    method = "fgsea_simple"

    def __init__(self) -> None:
        self.params: list[Any] = []

    def __call__(self, call: Any) -> Any:
        from giap.function.gsea import SIMPLE_ENGINE

        self.params.append(call.params)
        return SIMPLE_ENGINE(call)


def _families(mp: pytest.MonkeyPatch, engine: _Engine, *, ec: bool = True,
              det: bool = True) -> None:  # fmt: skip
    from giap.context import sensitivity as S
    from giap.determinants.sensitivity import determinants_family
    from giap.function.sensitivity import enriched_categories_family

    fams: dict[str, Any] = {}
    if ec:
        fams[sf.EC] = sf.with_engine(enriched_categories_family, engine)
    if det:
        fams[sf.DET] = determinants_family
    mp.setattr(S, "FAMILIES", fams)


class _Records(logging.Handler):
    def __init__(self) -> None:
        super().__init__(logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


@dataclass
class Run:
    out: Path
    rows: pd.DataFrame
    records: list[logging.LogRecord]
    engine: _Engine


def _run(
    w: World,
    out: Path,
    cfg: Any,
    *,
    drop: tuple[str, ...] = (),
    replace: dict[str, Path] | None = None,
    ec: bool = True,
    det: bool = True,
) -> Run:
    from giap.steps.registry import run_step

    _spec()
    ins = {k: v for k, v in w.inputs.items() if k not in drop} | dict(replace or {})
    eng = _Engine()
    h = _Records()
    root = logging.getLogger("giap")
    old = root.level
    root.addHandler(h)
    root.setLevel(logging.DEBUG)
    try:
        with pytest.MonkeyPatch.context() as mp:
            _families(mp, eng, ec=ec, det=det)
            run_step("sensitivity", cfg, ins, {"sensitivity": out})
    finally:
        root.removeHandler(h)
        root.setLevel(old)
    return Run(out, sd.read_sensitivity(out), h.records, eng)


@pytest.fixture(scope="module")
def full(world: World) -> Iterator[Run]:
    yield _run(world, world.root / "out" / "sensitivity.tsv", sf.config())


def _events(records: list[logging.LogRecord]) -> list[logging.LogRecord]:
    return [r for r in records if getattr(r, "event", None)]


def _validate(tsv: Path) -> None:
    frictionless = pytest.importorskip("frictionless")
    import giap

    path = Path(giap.__file__).resolve().parent / "schemas" / "tables" / "sensitivity.schema.json"
    d = {k: v for k, v in json.loads(path.read_text(encoding="utf-8")).items()
         if k != "coordinates"}  # fmt: skip
    res = frictionless.Resource(path=tsv.name, basepath=str(tsv.parent),
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    report = res.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


# --------------------------------------------------------------------------- wiring


def test_step_declares_the_optional_inputs() -> None:
    """Item 3: the function resources and the determinants features are optional inputs; the
    required inputs stay M4-19's."""
    spec = _spec()
    assert set(RESOURCES) <= set(spec.optional_inputs)
    assert "determinants_features" in spec.optional_inputs
    assert set(spec.inputs) == REQUIRED
    assert spec.outputs == ("sensitivity",)


def _sensitivity_call() -> str:
    for nf in sorted(WORKFLOW.rglob("*.nf")):
        text = re.sub(r"(?m)//.*$", "", nf.read_text(encoding="utf-8"))
        m = re.search(r"stepIn\(\s*'sensitivity'\s*,(.*?)\n\s*SENSITIVITY\w*\(", text, re.S)
        if m:
            return m.group(1)
    raise AssertionError("no stepIn('sensitivity', ...) followed by SENSITIVITY(...)")


def test_workflow_passes_the_function_and_determinants_channels() -> None:
    """Item 3: the ``sensitivity`` step's inputs include ``determinants_features`` and the
    function resources (``functionResources()``, a channel named ``function_resources`` or the
    item keys). DETERMINANTS runs after CONTEXT, so the call has to move after it (static check;
    the run-time check is the e2e)."""
    call = _sensitivity_call()
    assert "determinants_features" in call
    assert re.search(r"functionResources|function_resources|go_basic", call), call


# --------------------------------------------------------------------------- the two families


def test_writes_both_families(full: Run) -> None:
    rows = full.rows
    got = sorted(zip(rows["family"], rows["variant"], strict=True))
    assert got == sorted(
        (f, v) for f in (sf.DET, sf.EC) for v in ("missing=min", "missing=exclude")
    )
    _validate(full.out)
    for r in rows.itertuples(index=False):
        n_p, n_v, n_s = int(r.n_primary), int(r.n_variant), int(r.n_shared)
        assert n_s <= min(n_p, n_v)
        if n_p + n_v == 0:
            assert r.jaccard == "NA"
        else:
            assert float(r.jaccard) == pytest.approx(n_s / (n_p + n_v - n_s), abs=1e-6)
    assert (rows.loc[rows["family"] == sf.EC, "agreement"] == "NA").all()


def test_enriched_categories_n_primary_is_the_primary_rerun(world: World, full: Run) -> None:
    """n_primary = the family on the primary gene table with the step's collections, stream and
    base config (cog:J is flagged: its 40 genes are the most abundant)."""
    from giap.function.collections import build_collections
    from giap.function.gsea import SIMPLE_ENGINE
    from giap.function.sensitivity import enriched_categories_family
    from giap.genetable.eligibility import eligible
    from giap.io.function_resources import read_resources
    from giap.rng import for_step

    res = read_resources({k: world.inputs[k] for k in RESOURCES})
    colls, _ = build_collections(world.genes, eligible(world.genes), res)
    want = enriched_categories_family(
        world.genes, sf.base_config(), for_step(f"sensitivity/{sf.EC}", sf.SEED),
        collections=colls, tus=pd.read_parquet(world.inputs["tus"]), engine=SIMPLE_ENGINE,
    )  # fmt: skip
    assert "cog:J:gsea" in want, "fixture: COG J not enriched"
    ec = full.rows[full.rows["family"] == sf.EC]
    assert set(ec["n_primary"].astype(int)) == {len(want)}


def test_gsea_reruns_at_b_199(full: Run) -> None:
    assert full.engine.params, "GSEA was never run"
    assert {p.n_perm_simple for p in full.engine.params} == {sf.SENS_B}
    assert {p.sample_size for p in full.engine.params} == {sf.SAMPLE_SIZE}


def test_determinants_rows_min_exclude_agree(world: World, full: Run) -> None:
    """Imputed genes never enter the determinants models (M8-08 rule 1), so with common random
    numbers ``min`` / ``exclude`` reproduce the primary: agreement 1, n_variant = n_shared =
    n_primary, which is the primary rerun's flagged count."""
    from giap.determinants.sensitivity import determinants_family
    from giap.rng import for_step

    det = full.rows[full.rows["family"] == sf.DET]
    want = determinants_family(world.genes, sf.base_config(),
                               for_step(f"sensitivity/{sf.DET}", sf.SEED),
                               features=world.features)  # fmt: skip
    assert want is not None and np.isfinite(float(want.cv_r2))
    assert set(det["n_primary"].astype(int)) == {len(want.flagged)}
    assert (det["n_variant"] == det["n_primary"]).all()
    assert (det["n_shared"] == det["n_primary"]).all()
    assert (det["agreement"].astype(float) == 1.0).all()


def test_determinants_cv_r2_goes_to_an_event(full: Run) -> None:
    """Item 2: the CV R² of each run is an event (the §8.6 table has no R² column)."""
    hits = [r for r in _events(full.records) if "determinants" in str(r.event)]  # type: ignore[attr-defined]
    assert hits, "no determinants event"
    text = " ".join(json.dumps(getattr(r, "data", None), default=str) for r in hits)
    for v in ("missing=min", "missing=exclude"):
        assert v in text, f"no CV R² of {v} in the events"


def test_same_seed_byte_identical(world: World, full: Run) -> None:
    again = _run(world, world.root / "again" / "sensitivity.tsv", sf.config())
    assert again.out.read_bytes() == full.out.read_bytes()


def test_two_peptide_adds_the_single_peptide_variant(world_2pep: World) -> None:
    run = _run(world_2pep, world_2pep.root / "out" / "sensitivity.tsv", sf.config(), det=False)
    assert set(run.rows["variant"]) == {"missing=min", "missing=exclude",
                                        "single_peptide_detected"}  # fmt: skip
    assert set(run.rows["family"]) == {sf.EC}


# --------------------------------------------------------------------------- absent inputs


def _skipped(run: Run, family: str) -> None:
    assert family not in set(run.rows["family"]), "rows written for a family without inputs"
    assert any(family in r.getMessage() for r in _events(run.records)), \
        f"no event naming {family}"  # fmt: skip


def test_no_function_resources_no_rows_and_event(world: World, tmp_path: Path) -> None:
    run = _run(world, tmp_path / "s.tsv", sf.config(), drop=RESOURCES, det=False)
    assert len(run.rows) == 0
    _skipped(run, sf.EC)


def test_function_module_off_no_rows(world: World, tmp_path: Path) -> None:
    cfg = sf.config(modules=["core"])
    run = _run(world, tmp_path / "s.tsv", cfg, det=False)
    assert len(run.rows) == 0
    _skipped(run, sf.EC)


def test_no_determinants_features_no_rows_and_event(world: World, tmp_path: Path) -> None:
    run = _run(world, tmp_path / "s.tsv", sf.config(), drop=("determinants_features",), ec=False)
    assert len(run.rows) == 0
    _skipped(run, sf.DET)


def test_determinants_below_the_minimum_no_rows_and_event(world: World, tmp_path: Path) -> None:
    from giap.io.determinants_tables import write_determinants_features

    small = tmp_path / "determinants_features.tsv"
    write_determinants_features(small, _features(world.genes, keep=60))
    run = _run(world, tmp_path / "s.tsv", sf.config(), replace={"determinants_features": small},
               ec=False)  # fmt: skip
    assert len(run.rows) == 0
    _skipped(run, sf.DET)
