"""M8-06a: the ``function`` step and its three tables (docs/specs/M8-06-function-step-and-run.md,
items 1-5, examples 4-9, the 06a degenerate cases and lead decisions 1-3; M8-04 decision 13;
GIPPA_SPEC.md §6.2, §6.3, §6.6, §6.7, §7.7, §9.1-§9.4).

Interface pinned here (test-author reading of the spec):

* ``@step("function", module="function", inputs={gene_table, replicons, tu_folds},
  optional_inputs=FUNCTION_RESOURCES, outputs={category_summary, enrichment_detection,
  enrichment_gsea})``; ``run_function(genes, replicons, folds, res, cfg)`` is the pure core and
  returns ``(summary, detection, gsea, run_flags)``.
* Spec item 2's order. The summary draws from ``for_step("function.bootstrap", seed)`` and GSEA
  from ``for_step("function.gsea", seed)``; ORA draws nothing. The oracle in
  ``helpers.function_step_data.oracle_tables`` calls exactly that chain, so the step's tables
  must equal it byte for byte.
* Lead decision 3 (M8-04 decision 13): a header-only ``tu_folds.parquet`` means the step passes
  ``folds=None`` explicitly. Tests marked **[#179]** need the #179 fix round: every ORA row
  carries ``ora_uncalibrated``, and a skipped recalibration puts ``tu_recal_skipped`` on every
  ORA row too.
* Item 3: each run flag is logged once as the WARNING line ``function: run flag <flag>``, then
  one event ``function.run_flags`` ``{"flags": [...]}`` and one event ``function.collections``
  ``{<collection>: {"n_terms", "coverage", "n_unknown_ids"}}`` for each built collection
  (``extra={"event": ..., "data": ...}`` on ``ctx.log``, as ``operons.missed_operons``).
* No eligible gene: three header-only files, run flags ``build_collections`` flags +
  ``["no_eligible_genes"]``, and the summary, ORA and GSEA are never called.

Example 4's toy gene table and the ~600-gene world (non-constant detectability, synthetic function
DB from ``helpers.function_db``, lead decision 1) are in ``helpers.function_step_data``.
"""

from __future__ import annotations

import json
import logging
from collections.abc import Iterator, Mapping
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.function_db import FUNCTION_RESOURCES, FunctionDb, write_function_db
from helpers.function_step_data import (
    oracle_tables,
    read_tu_folds,
    toy4,
    world_pre,
    write_gene_table,
    write_replicons,
    write_tu_folds,
)

from gippa.config import GippaConfig


@pytest.fixture(autouse=True, scope="module")
def _simple_gsea_engine() -> Iterator[None]:
    """The step runs on ``fgsea_simple``, never R: unit tests use fakes (ADR-0038 §4, M8-17b 10).
    Module scope, so module-scoped fixtures that run the step see it too."""
    import gippa.steps.function as fn
    from gippa.function.gsea import SIMPLE_ENGINE

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(fn, "_LazyRFgsea", lambda log: SIMPLE_ENGINE)
        yield


SCHEMA_DIR = Path(__file__).resolve().parents[3] / "src" / "gippa" / "schemas" / "tables"
REGISTRY_JSON = (Path(__file__).resolve().parents[3] / "workflow" / "tests" / "data" / "annotate"
                 / "step_registry.json")  # fmt: skip
OUTPUTS = {
    "category_summary": "category_summary.tsv",
    "enrichment_detection": "enrichment_detection.tsv",
    "enrichment_gsea": "enrichment_gsea.tsv",
}
REQUIRED = {"gene_table", "replicons", "tu_folds"}
FLAGS_EVENT = "function.run_flags"
COLL_EVENT = "function.collections"
J_NAME = "Translation, ribosomal structure and biogenesis"
RUN_ONLY = ("layer_missing:cog", "layer_missing:kegg", "layer_missing:go",
            "layer_missing:abundance", "bakta_db_light", "tu_recal_skipped",
            "too_few_ranked_tus", "no_eligible_genes")  # fmt: skip


def _cfg(seed: int = 7, **over: Any) -> GippaConfig:
    raw: dict[str, Any] = {"seed": seed, "stats": {"bootstrap": 100},
                           "function": {"gsea": {"n_perm_simple": 200}}}  # fmt: skip
    for k, v in over.items():
        raw[k] = {**raw.get(k, {}), **v} if isinstance(v, dict) else v
    return GippaConfig.model_validate(raw)


def _columns() -> dict[str, tuple[str, ...]]:
    from gippa.function import gsea, ora, summary

    return {"category_summary": tuple(summary.COLUMNS),
            "enrichment_detection": tuple(ora.COLUMNS),
            "enrichment_gsea": tuple(gsea.COLUMNS)}  # fmt: skip


# ----------------------------------------------------------------------------- running the step


class _Recorder(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


class Run:
    def __init__(self, out: Path, records: list[logging.LogRecord], inputs: dict[str, Path],
                 cfg: GippaConfig) -> None:  # fmt: skip
        self.out = out
        self.records = records
        self.inputs = inputs
        self.cfg = cfg
        self.paths = {k: out / "tables" / v for k, v in OUTPUTS.items()}

    def text(self, key: str) -> str:
        return self.paths[key].read_text(encoding="utf-8")

    def frame(self, key: str) -> pd.DataFrame:
        return pd.read_csv(self.paths[key], sep="\t", dtype=str, keep_default_na=False)

    def events(self, name: str) -> list[dict[str, Any]]:
        out = []
        for r in self.records:
            if getattr(r, "event", None) == name:
                data = getattr(r, "data", None)
                assert isinstance(data, dict), data
                out.append(data)
        return out

    def run_flags(self) -> list[str]:
        found = self.events(FLAGS_EVENT)
        assert len(found) == 1, f"expected one {FLAGS_EVENT} event, got {len(found)}"
        flags = found[0]["flags"]
        assert isinstance(flags, list), flags
        return [str(f) for f in flags]

    def warnings(self) -> list[str]:
        return [r.getMessage() for r in self.records if r.levelno >= logging.WARNING]


def run_function_step(inputs: Mapping[str, Path], out: Path, cfg: GippaConfig) -> Run:
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    outputs = {k: out / "tables" / v for k, v in OUTPUTS.items()}
    rec = _Recorder()
    log = logging.getLogger("gippa.steps.function")
    old = log.level
    log.addHandler(rec)
    log.setLevel(logging.DEBUG)
    try:
        reg.run_step("function", cfg, dict(inputs), outputs)
    finally:
        log.removeHandler(rec)
        log.setLevel(old)
    return Run(out, rec.records, dict(inputs), cfg)


def toy_inputs(root: Path, genes: pd.DataFrame, resources: Mapping[str, Path],
               folds: Mapping[str, int] | None = None) -> dict[str, Path]:  # fmt: skip
    ins = {
        "gene_table": write_gene_table(root / "genetable" / "genes.parquet", genes),
        "replicons": write_replicons(root / "replicons.parquet"),
        "tu_folds": write_tu_folds(root / "genetable" / "tu_folds.parquet", folds or {}),
    }
    return ins | dict(resources)


def _empty_xrefs() -> pd.DataFrame:
    return pd.DataFrame({"locus_tag": [], "db": [], "id": []}, dtype=object)


@pytest.fixture(scope="module")
def toy_db(tmp_path_factory: pytest.TempPathFactory) -> FunctionDb:
    return write_function_db(tmp_path_factory.mktemp("toy_db"), _empty_xrefs())


def _validate(path: Path, stem: str) -> Any:
    frictionless = pytest.importorskip("frictionless")
    d = json.loads((SCHEMA_DIR / f"{stem}.schema.json").read_text(encoding="utf-8"))
    d = {k: v for k, v in d.items() if k != "coordinates"}
    res = frictionless.Resource(path=path.name, basepath=str(path.parent),
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    return res.validate()


def _assert_valid(run: Run) -> None:
    for key, path in run.paths.items():
        report = _validate(path, key)
        assert report.valid, (key, report.flatten(["rowNumber", "fieldName", "type", "note"]))


def _header_only(run: Run, key: str) -> None:
    assert run.text(key) == "\t".join(_columns()[key]) + "\n", key


def _tokens(cell: str) -> set[str]:
    return set() if cell == "NA" else set(cell.split(";"))


# ----------------------------------------------------------------------------- registry


def test_function_step_is_registered_with_the_spec_keys() -> None:
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    spec = reg.get("function")
    assert spec.kind == "py"
    assert spec.module == "function"
    assert spec.precomputed is None
    assert set(spec.inputs) == REQUIRED
    assert set(spec.optional_inputs) == set(FUNCTION_RESOURCES)
    assert set(spec.outputs) == {*OUTPUTS, "function_support"}  # M9-17 plot support
    assert spec.optional_outputs == ("function_support",)


def test_function_resources_constant_is_the_five_setup_items() -> None:
    from gippa.setupdb.installed import FUNCTION_RESOURCES as LIVE

    assert LIVE == ("go_basic", "goslim_generic", "kegg_ko_pathway", "kegg_pathway_names",
                    "cog2020_fun")  # fmt: skip


def test_function_runs_after_detectability_which_produces_tu_folds() -> None:
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    det = reg.get("detectability")
    assert "tu_folds" in det.outputs
    order = reg.topo_order([s.name for s in reg.all_steps()])
    assert order.index("detectability") < order.index("function")


def test_function_is_switched_by_modules() -> None:
    """§6.3 / item 6: planned only when ``"function" in cfg.modules`` (default on)."""
    from gippa.runner.local import _enabled
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    spec = reg.get("function")
    assert "function" in GippaConfig(seed=1).modules
    assert _enabled(GippaConfig(seed=1), spec)
    assert not _enabled(GippaConfig.model_validate({"seed": 1, "modules": ["core"]}), spec)


def test_step_registry_json_pins_function_and_tu_folds() -> None:
    """The registry pin (``test_nf_annotate_registry.py`` checks the whole file)."""
    got = json.loads(REGISTRY_JSON.read_text(encoding="utf-8"))
    assert got["function"] == {
        "inputs": sorted(REQUIRED),
        "optional_inputs": sorted(FUNCTION_RESOURCES),
        "outputs": sorted([*OUTPUTS, "function_support"]),
    }
    assert got["detectability"]["outputs"] == sorted(
        ["gene_table", "detectability_summary", "tu_folds"]
    )


# ----------------------------------------------------------------------------- example 4


@pytest.fixture(scope="module")
def ex4(toy_db: FunctionDb, tmp_path_factory: pytest.TempPathFactory) -> Run:
    root = tmp_path_factory.mktemp("ex4")
    ins = toy_inputs(root, toy4(), toy_db.inputs(["cog2020_fun"]))
    return run_function_step(ins, root / "out", _cfg())


def test_example_4_three_tables_written_with_module_columns(ex4: Run) -> None:
    for key in OUTPUTS:
        assert ex4.paths[key].is_file(), key
        assert tuple(ex4.text(key).split("\n", 1)[0].split("\t")) == _columns()[key], key


def test_example_4_category_summary_row_is_exact(ex4: Run) -> None:
    lines = ex4.text("category_summary").splitlines()
    assert len(lines) == 2, lines
    want = ["cog", "J", J_NAME, "3", "2", "1.8", "0.2", "NA", "NA", "NA", "0.75", "ci_few_units"]
    assert lines[1].split("\t") == want


def test_example_4_detection_and_gsea_header_only(ex4: Run) -> None:
    """3 TUs < ``ora_size`` 5; no abundance -> GSEA header-only."""
    _header_only(ex4, "enrichment_detection")
    _header_only(ex4, "enrichment_gsea")


def test_example_4_run_flag_union_order_is_pinned(ex4: Run) -> None:
    """Item 2 step 6: collections' flags, then ORA's, then GSEA's, first-seen order."""
    assert ex4.run_flags() == ["layer_missing:kegg", "layer_missing:go", "tu_recal_skipped",
                               "layer_missing:abundance"]  # fmt: skip


def test_example_4_each_run_flag_logged_once_as_a_warning(ex4: Run) -> None:
    lines = ex4.warnings()
    for flag in ex4.run_flags():
        assert lines.count(f"function: run flag {flag}") == 1, (flag, lines)
    logged = [m for m in lines if m.startswith("function: run flag ")]
    assert len(logged) == len(ex4.run_flags()), logged


def test_example_4_collections_event(ex4: Run) -> None:
    found = ex4.events(COLL_EVENT)
    assert len(found) == 1, found
    data = found[0]
    assert set(data) == {"cog"}
    assert data["cog"]["n_terms"] == 1
    assert data["cog"]["coverage"] == pytest.approx(0.75)
    assert data["cog"]["n_unknown_ids"] == 0
    json.dumps(data)  # §6.7 events are plain JSON
    json.dumps(ex4.events(FLAGS_EVENT))


def test_example_4_run_function_is_the_pure_core(toy_db: FunctionDb) -> None:
    from gippa.io.function_resources import read_resources
    from gippa.steps.function import run_function

    reps = pd.DataFrame({"contig": ["chr1"], "length": [1_000_000], "topology": ["circular"],
                         "kind": ["chromosome"]})  # fmt: skip
    res = read_resources(toy_db.inputs(["cog2020_fun"]))
    summary, detection, gsea, flags = run_function(toy4(), reps, {}, res, _cfg())
    assert flags == ["layer_missing:kegg", "layer_missing:go", "tu_recal_skipped",
                     "layer_missing:abundance"]  # fmt: skip
    assert len(summary) == 1 and summary.iloc[0]["category"] == "J"
    assert len(detection) == 0 and len(gsea) == 0


def test_example_4_bakta_light_adds_the_run_flag(toy_db: FunctionDb, tmp_path: Path) -> None:
    """Item 3: ``bakta_db_light`` is a run flag (from ``build_collections``); tables as usual."""
    genes = toy4(annotation_source=["bakta_light"] * 4)
    run = run_function_step(toy_inputs(tmp_path, genes, toy_db.inputs(["cog2020_fun"])),
                            tmp_path / "out", _cfg())  # fmt: skip
    assert run.run_flags() == ["layer_missing:kegg", "layer_missing:go", "bakta_db_light",
                               "tu_recal_skipped", "layer_missing:abundance"]  # fmt: skip
    assert len(run.text("category_summary").splitlines()) == 2


def test_example_4_validates(ex4: Run) -> None:
    _assert_valid(ex4)


# ----------------------------------------------------------------------------- example 5


@pytest.fixture(scope="module")
def ex5(tmp_path_factory: pytest.TempPathFactory) -> Run:
    root = tmp_path_factory.mktemp("ex5")
    return run_function_step(toy_inputs(root, toy4(), {}), root / "out", _cfg())


def test_example_5_no_resource_three_header_only_files(ex5: Run) -> None:
    for key in OUTPUTS:
        _header_only(ex5, key)


def test_example_5_run_flags(ex5: Run) -> None:
    flags = ex5.run_flags()
    assert flags[:3] == ["layer_missing:cog", "layer_missing:kegg", "layer_missing:go"]
    assert flags == ["layer_missing:cog", "layer_missing:kegg", "layer_missing:go",
                     "tu_recal_skipped", "layer_missing:abundance"]  # fmt: skip


def test_example_5_validates(ex5: Run) -> None:
    _assert_valid(ex5)


# ----------------------------------------------------------------------------- example 6


def _forbid_steps_3_to_5(monkeypatch: pytest.MonkeyPatch) -> None:
    """The guard must not call the summary, ORA or GSEA (spec item 2, "No eligible gene")."""
    import importlib

    def boom(*_: Any, **__: Any) -> Any:
        raise AssertionError("called despite no eligible gene")

    names = {"gippa.function.summary": "category_summary",
             "gippa.function.ora": "detection_enrichment",
             "gippa.function.gsea": "gsea_table"}  # fmt: skip
    step_mod = importlib.import_module("gippa.steps.function")
    for mod, fn in names.items():
        monkeypatch.setattr(importlib.import_module(mod), fn, boom)
        monkeypatch.setattr(step_mod, fn, boom, raising=False)


@pytest.mark.parametrize("case", ["all_partial", "empty"])
def test_example_6_no_eligible_gene(
    case: str, toy_db: FunctionDb, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippa.steps import load_step_modules

    load_step_modules()
    _forbid_steps_3_to_5(monkeypatch)
    genes = toy4(partial=[True] * 4) if case == "all_partial" else toy4().iloc[0:0]
    run = run_function_step(toy_inputs(tmp_path, genes, toy_db.inputs(["cog2020_fun"])),
                            tmp_path / "out", _cfg())  # fmt: skip
    for key in OUTPUTS:
        _header_only(run, key)
    assert run.run_flags() == ["layer_missing:kegg", "layer_missing:go", "no_eligible_genes"]
    assert "function: run flag no_eligible_genes" in run.warnings()
    _assert_valid(run)


# ----------------------------------------------------------------------------- failing closed


def test_unparsable_resource_fails_the_step(toy_db: FunctionDb, tmp_path: Path) -> None:
    """Item 2 step 1: a bad OBO raises (exit 4 = ``StepError``); no fallback, no table left."""
    from gippa.steps.registry import StepError

    bad = tmp_path / "db" / "go" / "go-basic.obo"
    bad.parent.mkdir(parents=True)
    bad.write_text("format-version: 1.2\n\n[Term]\nid: GO:0000001\nname: x\n"
                   "namespace: biological_process\nis_a: GO:0009999 ! dangling\n",
                   encoding="utf-8")  # fmt: skip
    ins = toy_inputs(tmp_path, toy4(), toy_db.inputs(["cog2020_fun"])) | {"go_basic": bad}
    with pytest.raises(StepError):
        run_function_step(ins, tmp_path / "out", _cfg())
    for name in OUTPUTS.values():
        assert not (tmp_path / "out" / "tables" / name).exists()


def _toy5(**g5: Any) -> pd.DataFrame:
    base = toy4()
    extra = toy4().iloc[[0]].copy()
    extra["locus_tag"] = "g5"
    extra["tu_id"] = "TU_5"
    for k, v in g5.items():
        extra[k] = pd.Series([v], index=extra.index, dtype=object)
    return pd.concat([base, extra], ignore_index=True)


@pytest.mark.parametrize("bad", [{"p_detect": float("nan")}, {"tu_id": None}])
def test_eligible_gene_without_p_detect_or_tu_fails_closed(
    bad: dict[str, Any], toy_db: FunctionDb, tmp_path: Path
) -> None:
    """Degenerate case: ``ValueError`` from ``eligible_genes`` -> exit 4, naming the gene."""
    from gippa.steps.registry import StepError

    ins = toy_inputs(tmp_path, _toy5(**bad), toy_db.inputs(["cog2020_fun"]))
    with pytest.raises(StepError) as err:
        run_function_step(ins, tmp_path / "out", _cfg())
    cause = err.value.__cause__
    assert isinstance(cause, ValueError), cause
    assert "g5" in str(cause)


# ----------------------------------------------------------------------------- the 600-gene world


class World:
    def __init__(self, root: Path) -> None:
        from gippa.genetable.covariates import COVARIATE_COLUMNS
        from gippa.genetable.detectability import fit_detectability
        from gippa.io.genes_tsv import read_genes_parquet
        from gippa.rng import for_step
        from gippa.steps import load_step_modules
        from gippa.steps import registry as reg

        load_step_modules()
        self.root = root
        _, xrefs = world_pre()
        self.db = write_function_db(root / "db", xrefs)
        pre, _ = world_pre(root / "db")
        self.pre = pre
        pre_path = write_gene_table(root / "genetable" / "genes_pre.parquet", pre)
        cfg = _cfg()
        self.gene_table = root / "genetable" / "genes.parquet"
        reg.run_step("detectability", cfg, {"genes_pre": pre_path},
                     {"gene_table": self.gene_table,
                      "detectability_summary": root / "genetable" / "det.parquet"})  # fmt: skip
        # the folds the detectability step draws (its own stream), as item 1 persists them
        p = read_genes_parquet(pre_path)
        idx = pd.Index(p["locus_tag"].astype(str))
        det = pd.Series(p["detected"].astype(bool).to_numpy(), index=idx, dtype="boolean")
        cov = p[["locus_tag", *COVARIATE_COLUMNS]].astype(
            dict.fromkeys(COVARIATE_COLUMNS, "float64")
        )
        fit = fit_detectability(det, cov, pd.Series(p["tu_id"].to_numpy(), index=idx),
                                cfg.detectability, for_step("detectability", cfg.seed))  # fmt: skip
        assert not fit.constant, "the world must give a non-constant model"
        self.folds = dict(fit.fold_of_tu)
        self.tu_folds = write_tu_folds(root / "genetable" / "tu_folds.parquet", self.folds)
        self.empty_folds = write_tu_folds(root / "genetable" / "tu_folds_empty.parquet", {})
        self.replicons = write_replicons(root / "replicons.parquet")
        self.linear = write_replicons(root / "replicons_linear.parquet", topology="linear")

    def inputs(self, items: Any = FUNCTION_RESOURCES, **over: Path) -> dict[str, Path]:
        ins = {"gene_table": self.gene_table, "replicons": self.replicons,
               "tu_folds": self.tu_folds}  # fmt: skip
        return ins | self.db.inputs(items) | over

    def run(self, name: str, cfg: GippaConfig | None = None, **kw: Any) -> Run:
        return run_function_step(self.inputs(**kw), self.root / name, cfg or _cfg())


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> World:
    return World(tmp_path_factory.mktemp("function_world"))


@pytest.fixture(scope="module")
def full(world: World) -> Run:
    return world.run("full")


@pytest.fixture(scope="module")
def skipped(world: World) -> Run:
    """Example 7's second half: the same run with a header-only ``tu_folds.parquet``."""
    return world.run("skipped", tu_folds=world.empty_folds)


def test_world_tables_equal_the_item_2_chain(world: World, full: Run, tmp_path: Path) -> None:
    paths, union = oracle_tables(world.gene_table, world.replicons, world.folds,
                                 world.db.inputs(), full.cfg, tmp_path / "oracle")  # fmt: skip
    for key in OUTPUTS:
        assert full.paths[key].read_bytes() == paths[key].read_bytes(), key
    assert full.run_flags() == union


def test_world_tables_have_rows_and_validate(full: Run) -> None:
    for key in OUTPUTS:
        assert len(full.frame(key)) > 0, key
    _assert_valid(full)


def test_world_tables_carry_no_coordinates(full: Run) -> None:
    for key in OUTPUTS:
        header = set(full.text(key).split("\n", 1)[0].split("\t"))
        assert not header & {"contig", "start", "end", "strand"}, key


def test_world_no_run_only_token_in_a_summary_or_gsea_row(full: Run, skipped: Run) -> None:
    for run in (full, skipped):
        for key in ("category_summary", "enrichment_gsea"):
            toks = set().union(*map(_tokens, run.frame(key)["flags"]))
            assert not toks & set(RUN_ONLY), (key, toks)


def test_world_collections_event_matches_the_collections(world: World, full: Run) -> None:
    from gippa.function.collections import COLLECTIONS, build_collections
    from gippa.genetable.eligibility import eligible
    from gippa.io.function_resources import read_resources
    from gippa.io.genes_tsv import read_genes_parquet

    genes = read_genes_parquet(world.gene_table)
    colls, _ = build_collections(genes, eligible(genes), read_resources(world.db.inputs()))
    found = full.events(COLL_EVENT)
    assert len(found) == 1
    data = found[0]
    assert list(data) == [c for c in COLLECTIONS if c in {g.collection for g in colls}]
    for c in colls:
        got = data[c.collection]
        assert got["n_terms"] == len(c.sets)
        assert got["coverage"] == pytest.approx(c.coverage, rel=1e-12)
        assert got["n_unknown_ids"] == c.n_unknown_ids
    json.dumps(data)


# ----------------------------------------------------------------------------- example 7, function


def test_example_7_fitted_folds_no_tu_recal_skipped(full: Run) -> None:
    assert "tu_recal_skipped" not in full.run_flags()
    rows = full.frame("enrichment_detection")
    assert len(rows) > 0
    assert not rows["flags"].map(lambda s: "tu_recal_skipped" in _tokens(s)).any()


def test_example_7_header_only_folds_skip_the_recalibration(full: Run, skipped: Run) -> None:
    assert "tu_recal_skipped" in skipped.run_flags()
    a, b = full.frame("enrichment_detection"), skipped.frame("enrichment_detection")
    assert a[["collection", "term", "direction"]].equals(b[["collection", "term", "direction"]])
    assert (a["p_value"] != b["p_value"]).any()


def test_decision_13_skip_puts_tu_recal_skipped_on_every_ora_row(skipped: Run) -> None:
    """[#179] M8-04 decision 13: a skipped recalibration flags every ORA row."""
    rows = skipped.frame("enrichment_detection")
    assert len(rows) > 0
    assert rows["flags"].map(lambda s: "tu_recal_skipped" in _tokens(s)).all()


def test_decision_13_every_ora_row_is_ora_uncalibrated(full: Run, skipped: Run) -> None:
    """[#179] ``CALIBRATION_FLAGS = ("ora_uncalibrated",)`` in pass 1 (ADR-0035 decision 3)."""
    for run in (full, skipped):
        rows = run.frame("enrichment_detection")
        assert rows["flags"].map(lambda s: "ora_uncalibrated" in _tokens(s)).all()


def test_decision_3_folds_passed_explicitly(
    world: World, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Lead decision 3: ``folds`` is a keyword; a header-only file gives ``folds=None``, a full
    one the file's ``{tu_id: fold}``."""
    import importlib

    from gippa.steps import load_step_modules

    load_step_modules()
    ora = importlib.import_module("gippa.function.ora")
    step_mod = importlib.import_module("gippa.steps.function")
    real = ora.detection_enrichment
    calls: list[dict[str, Any]] = []

    def spy(*args: Any, **kwargs: Any) -> Any:
        calls.append(dict(kwargs))
        return real(*args, **kwargs)

    monkeypatch.setattr(ora, "detection_enrichment", spy)
    monkeypatch.setattr(step_mod, "detection_enrichment", spy, raising=False)
    items = ["cog2020_fun"]
    world.run("spy_empty", items=items, tu_folds=world.empty_folds)
    world.run("spy_full", items=items)
    assert len(calls) == 2, calls
    assert "folds" in calls[0] and calls[0]["folds"] is None
    assert dict(calls[1]["folds"]) == read_tu_folds(world.tu_folds) == world.folds


# ----------------------------------------------------------------------------- example 8


def test_example_8_same_inputs_and_seed_byte_identical(world: World, full: Run) -> None:
    again = world.run("again")
    for key in OUTPUTS:
        assert again.paths[key].read_bytes() == full.paths[key].read_bytes(), key
    assert again.run_flags() == full.run_flags()


def test_example_8_seed_8_moves_only_the_random_parts(world: World, full: Run) -> None:
    seed8 = world.run("seed8", cfg=_cfg(seed=8))
    same = seed8.paths["enrichment_detection"].read_bytes()
    assert same == full.paths["enrichment_detection"].read_bytes()
    assert (
        seed8.paths["category_summary"].read_bytes() != full.paths["category_summary"].read_bytes()
    )
    a, b = full.frame("enrichment_gsea"), seed8.frame("enrichment_gsea")
    assert a[["collection", "term"]].equals(b[["collection", "term"]])
    assert (a["p_value"] != b["p_value"]).any()


def test_streams_are_independent(world: World, full: Run) -> None:
    """Item 2: separate sub-streams, so more bootstrap draws never move GSEA and more GSEA
    permutations never move the summary; ORA stays fixed in both."""
    more_boot = world.run("more_boot", cfg=_cfg(stats={"bootstrap": 150}))
    assert (
        more_boot.paths["enrichment_gsea"].read_bytes()
        == full.paths["enrichment_gsea"].read_bytes()
    )
    more_perm = world.run("more_perm", cfg=_cfg(function={"gsea": {"n_perm_simple": 300}}))
    assert (more_perm.paths["category_summary"].read_bytes()
            == full.paths["category_summary"].read_bytes())  # fmt: skip
    for run in (more_boot, more_perm):
        assert (run.paths["enrichment_detection"].read_bytes()
                == full.paths["enrichment_detection"].read_bytes())  # fmt: skip


# ----------------------------------------------------------------------------- degenerate cases


def test_only_cog_installed(world: World) -> None:
    run = world.run("cog_only", items=["cog2020_fun"])
    assert run.run_flags()[:2] == ["layer_missing:kegg", "layer_missing:go"]
    for key in OUTPUTS:
        df = run.frame(key)
        assert len(df) > 0, key
        assert set(df["collection"]) == {"cog"}, key


def test_go_basic_without_slim_has_no_go_slim(world: World) -> None:
    items = [i for i in FUNCTION_RESOURCES if i != "goslim_generic"]
    run = world.run("no_slim", items=items)
    assert set(run.frame("category_summary")["collection"]) == {"cog", "kegg"}
    for key in ("enrichment_detection", "enrichment_gsea"):
        colls = set(run.frame(key)["collection"])
        assert "go_slim" not in colls, key
        assert "go_bp" in colls, key


def test_linear_replicon_changes_only_genomic_cluster(world: World, tmp_path: Path) -> None:
    """Degenerate case "circular vs linear": only ``genomic_cluster`` reads the topology. This world
    has no origin-spanning term, so the hand-off of C1 to GSEA is pinned in
    test_function_step_review_181.py (#181 review B1)."""
    items = ["cog2020_fun", "kegg_ko_pathway", "kegg_pathway_names"]
    circ = world.run("circ", items=items)
    lin = world.run("lin", items=items, replicons=world.linear)
    paths, _ = oracle_tables(world.gene_table, world.linear, world.folds, world.db.inputs(items),
                             lin.cfg, tmp_path / "oracle")  # fmt: skip
    for key in OUTPUTS:
        assert lin.paths[key].read_bytes() == paths[key].read_bytes(), key
    for key in ("category_summary", "enrichment_detection"):
        assert lin.paths[key].read_bytes() == circ.paths[key].read_bytes(), key
    a, b = circ.frame("enrichment_gsea"), lin.frame("enrichment_gsea")
    rest = [c for c in a.columns if c != "flags"]
    assert a[rest].equals(b[rest])
    drop = lambda s: _tokens(s) - {"genomic_cluster"}  # noqa: E731
    assert a["flags"].map(drop).tolist() == b["flags"].map(drop).tolist()


def test_minus_strand_has_no_effect(world: World, full: Run) -> None:
    """Degenerate case: the tables have no coordinates; TU membership and ``tu_position`` come
    from C5, so flipping every strand changes nothing."""
    from gippa.io.genes_tsv import read_genes_parquet

    genes = read_genes_parquet(world.gene_table)
    genes["strand"] = "-"
    flipped = write_gene_table(world.root / "minus" / "genes.parquet", genes)
    run = world.run("minus_out", gene_table=flipped)
    for key in OUTPUTS:
        assert run.paths[key].read_bytes() == full.paths[key].read_bytes(), key


def test_no_abundance_world(world: World) -> None:
    """Degenerate case "no abundance": summary medians NA, ORA as usual, GSEA header-only."""
    from gippa.io.genes_tsv import read_genes_parquet

    genes = read_genes_parquet(world.gene_table)
    genes["abundance_log2"] = float("nan")
    genes["abundance_scale"] = None
    genes["abundance_type"] = None
    path = write_gene_table(world.root / "noab" / "genes.parquet", genes)
    items = ["cog2020_fun"]
    run = world.run("noab_out", items=items, gene_table=path)
    s = run.frame("category_summary")
    assert len(s) > 0 and set(s["median_abundance_log2"]) == {"NA"}
    assert len(run.frame("enrichment_detection")) > 0
    _header_only(run, "enrichment_gsea")
    assert "layer_missing:abundance" in run.run_flags()


def test_fewer_than_two_ranked_tus(toy_db: FunctionDb, tmp_path: Path) -> None:
    """Degenerate case "< 2 ranked TUs": GSEA header-only and ``too_few_ranked_tus``."""
    genes = toy4(abundance_log2=[21.0, float("nan"), float("nan"), float("nan")],
                 abundance_scale=["absolute"] * 4)  # fmt: skip
    run = run_function_step(toy_inputs(tmp_path, genes, toy_db.inputs(["cog2020_fun"])),
                            tmp_path / "out", _cfg())  # fmt: skip
    _header_only(run, "enrichment_gsea")
    flags = run.run_flags()
    assert "too_few_ranked_tus" in flags and "layer_missing:abundance" not in flags
