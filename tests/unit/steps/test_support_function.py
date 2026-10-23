"""M9-17: ``plots/support/gsea_running.tsv[.gz]`` from the ``function`` step (V14(a); docs/specs/
M9-17-plot-support-analysis.md "Support files", "Definitions", example 1 (on the step), example 5,
the "function off / no abundance" degenerate case and the "No change to any tables/ file bytes"
acceptance; lead decision 1).

The step runs on ``fgsea_simple`` (never R), as in tests/unit/steps/test_function_step.py. The
~600-gene world is the M8-06a one (``helpers.function_step_data``). Contract: tests/helpers/
support_analysis_data.py (output key ``function_support``, a directory).

Reading of the spec pinned here: the top 6 sets are the ``enrichment_gsea.tsv`` rows with a
``p_value``, by ``q_value``, then |``nes``| descending, then ``collection``, ``term``; each set has
one row per ranked TU in the engine's rank order (``rank_position`` 1..n, ``score`` the centred
TU score of ``giap.function.gsea.tu_ranking``, non-increasing), ``in_set`` true on the set's
TUs, and ``running_es`` fgsea's running sum, whose max deviation is the row's ``es``.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from helpers import support_analysis_data as sa
from helpers.function_db import FUNCTION_RESOURCES, write_function_db
from helpers.function_step_data import world_pre, write_gene_table, write_replicons, write_tu_folds

from giap.config import GiapConfig

KEY = sa.SUPPORT_KEYS["function"]
NAME = "gsea_running"
TABLES = {
    "category_summary": "category_summary.tsv",
    "enrichment_detection": "enrichment_detection.tsv",
    "enrichment_gsea": "enrichment_gsea.tsv",
}


@pytest.fixture(autouse=True, scope="module")
def _simple_gsea_engine() -> Iterator[None]:
    import giap.steps.function as fn
    from giap.function.gsea import SIMPLE_ENGINE

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(fn, "_LazyRFgsea", lambda log: SIMPLE_ENGINE)
        yield


def _cfg(seed: int = 7) -> GiapConfig:
    return GiapConfig.model_validate({"seed": seed, "stats": {"bootstrap": 100},
                                      "function": {"gsea": {"n_perm_simple": 200}}})  # fmt: skip


class World:
    """The M8-06a world: detectability-fitted gene table, folds, replicons, function DB."""

    def __init__(self, root: Path) -> None:
        from giap.steps import load_step_modules
        from giap.steps import registry as reg

        load_step_modules()
        self.root = root
        _, xrefs = world_pre()
        self.db = write_function_db(root / "db", xrefs)
        pre, _ = world_pre(root / "db")
        pre_path = write_gene_table(root / "genetable" / "genes_pre.parquet", pre)
        self.gene_table = root / "genetable" / "genes.parquet"
        reg.run_step("detectability", _cfg(), {"genes_pre": pre_path},
                     {"gene_table": self.gene_table,
                      "detectability_summary": root / "genetable" / "det.parquet"})  # fmt: skip
        self.tu_folds = write_tu_folds(root / "genetable" / "tu_folds.parquet", {})
        self.replicons = write_replicons(root / "replicons.parquet")

    def inputs(self, **over: Path) -> dict[str, Path]:
        ins = {"gene_table": self.gene_table, "replicons": self.replicons,
               "tu_folds": self.tu_folds}  # fmt: skip
        return ins | self.db.inputs(FUNCTION_RESOURCES) | over

    def run(self, name: str, *, support: bool = True, cfg: GiapConfig | None = None,
            **over: Path) -> dict[str, Path]:  # fmt: skip
        from giap.steps import registry as reg

        out = self.root / name
        outputs = {k: out / "tables" / v for k, v in TABLES.items()}
        if support:
            outputs[KEY] = out / KEY
        reg.run_step("function", cfg or _cfg(), self.inputs(**over), outputs)
        return outputs


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> World:
    return World(tmp_path_factory.mktemp("m917_function_world"))


@pytest.fixture(scope="module")
def full(world: World) -> dict[str, Path]:
    return world.run("full")


@pytest.fixture(scope="module")
def running(full: dict[str, Path]) -> pd.DataFrame:
    return sa.read(sa.need(full[KEY], NAME))


@pytest.fixture(scope="module")
def gsea(full: dict[str, Path]) -> pd.DataFrame:
    return pd.read_csv(full["enrichment_gsea"], sep="\t", dtype=str, keep_default_na=False)


def _top6(gsea: pd.DataFrame) -> list[tuple[str, str]]:
    t = gsea[gsea["p_value"] != "NA"].copy()
    t["_q"] = sa.nums(t["q_value"])
    t["_a"] = -np.abs(sa.nums(t["nes"]))
    t = t.sort_values(["_q", "_a", "collection", "term"], kind="stable")
    return list(zip(t["collection"], t["term"], strict=True))[: sa.GSEA_TOP]


def _sets(running: pd.DataFrame) -> list[tuple[str, str]]:
    return list(dict.fromkeys(zip(running["collection"], running["term"], strict=True)))


# ----------------------------------------------------------------------------- wiring


def test_registry_declares_the_optional_support_output() -> None:
    spec = sa.step_spec("function")
    assert KEY in spec.outputs and KEY in spec.optional_outputs
    for key in TABLES:
        assert key in spec.outputs and key not in spec.optional_outputs


def test_without_the_support_key_the_step_still_runs(world: World) -> None:
    paths = world.run("no_support", support=False)
    assert all(paths[k].is_file() for k in TABLES)


def test_tables_bytes_unchanged_by_the_support_writer(world: World, full: dict[str, Path]) -> None:
    """Acceptance: the tables with and without the support writer are byte-identical."""
    bare = world.run("bare", support=False)
    assert sa.tables_bytes({k: full[k] for k in TABLES}) == sa.tables_bytes(
        {k: bare[k] for k in TABLES}
    )


def test_published_under_plots_support_and_read_by_run_view(world: World, tmp_path: Path) -> None:
    work = tmp_path / "work"
    sa.run_with_spec_outputs("function", _cfg(), world.inputs(), work)
    view = sa.publish_and_view(work, tmp_path / "OUT")
    assert sa.find(tmp_path / "OUT" / "plots" / "support", NAME) is not None
    got = view.support(NAME)
    assert got is not None and list(got.columns) == list(sa.COLUMNS[NAME])


# ----------------------------------------------------------------------------- schema


def test_header_and_written_form(full: dict[str, Path]) -> None:
    path = sa.need(full[KEY], NAME)
    assert sa.header(path) == list(sa.COLUMNS[NAME])
    sa.check_written_form(path)
    assert sa.names(full[KEY]) == [NAME], "the function step writes gsea_running only"


def test_validates_against_its_schema(full: dict[str, Path]) -> None:
    import json

    frictionless = pytest.importorskip("frictionless")
    import giap.schemas

    path = sa.need(full[KEY], NAME)
    d = json.loads((Path(giap.schemas.__file__).parent / "plot_support" / f"{NAME}.json")
                   .read_text(encoding="utf-8"))  # fmt: skip
    d = {k: v for k, v in d.items() if k != "coordinates"}
    res = frictionless.Resource(path=path.name, basepath=str(path.parent), format="tsv",
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    report = res.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def test_cells(running: pd.DataFrame) -> None:
    assert set(running["in_set"]) <= {"true", "false"}
    assert (running["score"] != "NA").all() and (running["running_es"] != "NA").all()
    assert running["rank_position"].str.fullmatch(r"[1-9][0-9]*").all()


# ----------------------------------------------------------------------------- the top 6 sets


def test_the_top_6_sets_in_rule_order(running: pd.DataFrame, gsea: pd.DataFrame) -> None:
    want = _top6(gsea)
    assert len(want) == sa.GSEA_TOP, "the world must test at least 6 sets (vacuous otherwise)"
    assert _sets(running) == want


def test_one_row_per_ranked_tu_per_set(running: pd.DataFrame, world: World) -> None:
    from giap.function.gsea import tu_ranking
    from giap.io.genes_tsv import read_genes_parquet

    rank = tu_ranking(read_genes_parquet(world.gene_table))
    n = len(rank.tu_ids)
    score_of = dict(zip(rank.tu_ids, rank.scores.tolist(), strict=True))
    first: list[str] | None = None
    for _, grp in running.groupby(["collection", "term"], sort=False):
        assert len(grp) == n
        assert grp["rank_position"].astype(int).tolist() == list(range(1, n + 1))
        tus = grp["tu_id"].tolist()
        assert sorted(tus) == sorted(rank.tu_ids)
        first = first or tus
        assert tus == first, "every set uses the same ranking"
        s = sa.nums(grp["score"])
        assert np.all(np.diff(s) <= 1e-12), "score is the ranked (non-increasing) centred score"
        for tu, v in zip(tus, s, strict=True):
            assert sa.close(v, score_of[tu]), (tu, v, score_of[tu])


def test_in_set_counts_and_running_sum(running: pd.DataFrame, gsea: pd.DataFrame) -> None:
    """``running_es`` is the fgsea walk over the rows' |score| and ``in_set``; its max deviation
    is the table's ``es`` (spec: an assertion in the step, 1e-9 in memory; %.6g here)."""
    es_of = {(c, t): sa.num(e) for c, t, e in zip(gsea["collection"], gsea["term"], gsea["es"],
                                                   strict=True)}  # fmt: skip
    for (coll, term), grp in running.groupby(["collection", "term"], sort=False):
        hit = (grp["in_set"] == "true").to_numpy()
        assert 0 < hit.sum() < len(grp)
        w = np.abs(sa.nums(grp["score"]))
        nr, nm = w[hit].sum(), (~hit).sum()
        step = np.where(hit, w / nr if nr > 0 else 1.0 / hit.sum(), -1.0 / nm)
        walk = np.cumsum(step)
        got = sa.nums(grp["running_es"])
        np.testing.assert_allclose(got, walk, rtol=0, atol=1e-5, err_msg=f"{coll} {term}")
        hi, lo = max(0.0, got.max()), min(0.0, got.min())
        dev = hi if hi > -lo else lo
        assert sa.close(dev, es_of[(coll, term)], rel=1e-4, abs_=1e-5), (coll, term)


def test_in_set_tus_are_the_size_of_the_set(running: pd.DataFrame, gsea: pd.DataFrame) -> None:
    size = {(c, t): int(s) for c, t, s in zip(gsea["collection"], gsea["term"],
                                               gsea["size_tus"], strict=True)}  # fmt: skip
    for (coll, term), grp in running.groupby(["collection", "term"], sort=False):
        assert int((grp["in_set"] == "true").sum()) == size[(coll, term)], (coll, term)


# ----------------------------------------------------------------------------- determinism


def test_example_5_same_run_twice_byte_identical(world: World, full: dict[str, Path]) -> None:
    again = world.run("again")
    assert sa.tree_bytes(full[KEY]) == sa.tree_bytes(again[KEY])


# ----------------------------------------------------------------------------- degenerate


def test_no_abundance_no_gsea_running_file(world: World) -> None:
    from giap.io.genes_tsv import read_genes_parquet

    genes = read_genes_parquet(world.gene_table)
    genes["abundance_log2"] = float("nan")
    genes["abundance_scale"] = None
    genes["abundance_type"] = None
    path = write_gene_table(world.root / "noab" / "genes.parquet", genes)
    out = world.run("noab_out", gene_table=path)
    assert out[KEY].is_dir(), "the declared support directory is always created"
    assert sa.find(out[KEY], NAME) is None


def test_function_module_off_plans_no_support() -> None:
    """Function off: the step is not planned, so nothing writes ``gsea_running``."""
    from giap.runner.local import _enabled

    spec = sa.step_spec("function")
    off = GiapConfig.model_validate({"seed": 1, "modules": ["core"]})
    assert not _enabled(off, spec)
    assert KEY in spec.outputs
