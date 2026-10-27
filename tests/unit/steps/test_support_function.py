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
TU score of ``gippa.function.gsea.tu_ranking``, non-increasing), ``in_set`` true on the set's
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
from hypothesis import example, given, settings
from hypothesis import strategies as st

from gippa.config import GippaConfig

KEY = sa.SUPPORT_KEYS["function"]
NAME = "gsea_running"
TABLES = {
    "category_summary": "category_summary.tsv",
    "enrichment_detection": "enrichment_detection.tsv",
    "enrichment_gsea": "enrichment_gsea.tsv",
}


@pytest.fixture(autouse=True, scope="module")
def _simple_gsea_engine() -> Iterator[None]:
    import gippa.steps.function as fn
    from gippa.function.gsea import SIMPLE_ENGINE

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(fn, "_LazyRFgsea", lambda log: SIMPLE_ENGINE)
        yield


def _cfg(seed: int = 7) -> GippaConfig:
    return GippaConfig.model_validate({"seed": seed, "stats": {"bootstrap": 100},
                                      "function": {"gsea": {"n_perm_simple": 200}}})  # fmt: skip


class World:
    """The M8-06a world: detectability-fitted gene table, folds, replicons, function DB."""

    def __init__(self, root: Path) -> None:
        from gippa.steps import load_step_modules
        from gippa.steps import registry as reg

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

    def run(self, name: str, *, support: bool = True, cfg: GippaConfig | None = None,
            **over: Path) -> dict[str, Path]:  # fmt: skip
        from gippa.steps import registry as reg

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
    import gippa.schemas

    path = sa.need(full[KEY], NAME)
    d = json.loads((Path(gippa.schemas.__file__).parent / "plot_support" / f"{NAME}.json")
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
    from gippa.function.gsea import tu_ranking
    from gippa.io.genes_tsv import read_genes_parquet

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
    from gippa.io.genes_tsv import read_genes_parquet

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
    from gippa.runner.local import _enabled

    spec = sa.step_spec("function")
    off = GippaConfig.model_validate({"seed": 1, "modules": ["core"]})
    assert not _enabled(off, spec)
    assert KEY in spec.outputs


# ----------------------------------------------------------------------------- review round
# Lead decision 8 (PR #272): at a maxP = -minP near-tie the walk and the engine may land on
# opposite sides (|maxP + minP| ~ 1e-16); the check also accepts |dev| = |es| there (statistician
# B1, review N3; kills mutant M17). Publishing unlinks the other-suffix sibling (review B1b).

#: the brief's near-tie: centred scores already in rank order, set at ranks 2, 3, 4, 6, 7, 8, 9
NEAR_TIE_SCORES = (1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0)
NEAR_TIE_MEMBERS = (2, 3, 4, 6, 7, 8, 9)


def _ranking_inputs(
    scores: np.ndarray, members: np.ndarray, es: float | None = None
) -> tuple[object, list[object], pd.DataFrame]:
    """(TuRanking, collections, enrichment_gsea) for ``gsea_running`` on one set: TU_i holds gene
    g_i and score ``scores[i - 1]``; the table's ``es`` is the engine's (``gsea_es``) unless
    given."""
    from types import SimpleNamespace

    from gippa.function.gsea import TuRanking
    from gippa.stats.fgsea import gsea_es

    n = len(scores)
    ids = tuple(f"TU_{i}" for i in range(1, n + 1))
    rank = TuRanking(ids, np.asarray(scores, dtype=np.float64),
                     {t: (f"g_{t[3:]}",) for t in ids})  # fmt: skip
    sets = {"T1": frozenset(f"g_{int(m)}" for m in members)}
    colls = [SimpleNamespace(collection="go", sets=sets)]
    if es is None:
        es, _ = gsea_es(np.asarray(scores, float), np.asarray(members, np.int64))
    gsea = pd.DataFrame({"collection": ["go"], "term": ["T1"], "p_value": [0.01],
                         "q_value": [0.02], "nes": [1.5], "es": [float(es)]})  # fmt: skip
    gsea.attrs["rank_labels"] = list(ids)
    return rank, colls, gsea


def _call(
    monkeypatch: pytest.MonkeyPatch, scores: np.ndarray, members: np.ndarray,
    es: float | None = None,
) -> pd.DataFrame | None:  # fmt: skip
    import gippa.function.running as R

    rank, colls, gsea = _ranking_inputs(scores, members, es)
    monkeypatch.setattr(R, "tu_ranking", lambda genes: rank)
    return R.gsea_running(pd.DataFrame(), colls, gsea)  # type: ignore[arg-type]


def test_near_tie_fixture_is_a_near_tie() -> None:
    """The brief's example: the walk's max deviation sits on the maxP = -minP tie and its sign
    is the opposite of the engine's ES (fails if the fixture loses its teeth)."""
    from gippa.function.running import _max_dev
    from gippa.stats.fgsea import gsea_es, running_es

    s, S = np.array(NEAR_TIE_SCORES), np.array(NEAR_TIE_MEMBERS, dtype=np.int64)
    es, _ = gsea_es(s, S)
    dev, tie = _max_dev(running_es(s, S))
    assert tie and es != 0.0 and np.sign(dev) != np.sign(es)
    assert abs(abs(dev) - abs(es)) <= 1e-9


def test_gsea_running_accepts_the_engine_side_of_a_near_tie(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Decision 8 tie rule: no raise; the rows are the walk and the set's ranks."""
    from gippa.stats.fgsea import running_es

    s, S = np.array(NEAR_TIE_SCORES), np.array(NEAR_TIE_MEMBERS, dtype=np.int64)
    out = _call(monkeypatch, s, S)
    assert out is not None and list(out.columns) == list(sa.COLUMNS[NAME])
    assert out["rank_position"].tolist() == list(range(1, 11))
    assert out.loc[out["in_set"], "rank_position"].tolist() == list(NEAR_TIE_MEMBERS)
    np.testing.assert_allclose(out["running_es"].to_numpy(float), running_es(s, S), atol=1e-12)


def test_gsea_running_accepts_an_exact_tie_with_es_zero(monkeypatch: pytest.MonkeyPatch) -> None:
    """fgsea's exact tie (ES 0, side 0): accepted, as before (kills M17 with the near tie)."""
    from gippa.stats.fgsea import gsea_es

    s, S = np.ones(4), np.array([1, 4], dtype=np.int64)  # walk 0.5, 0, -0.5, 0
    assert gsea_es(s, S)[0] == 0.0
    assert _call(monkeypatch, s, S) is not None


@pytest.mark.parametrize("es", [0.2, -0.2, 1 / 3 + 1e-6])
def test_gsea_running_still_raises_on_a_wrong_es_at_a_tie(
    monkeypatch: pytest.MonkeyPatch, es: float
) -> None:
    """The tie clause accepts |es| = |dev| only: another value at the tie still raises."""
    s, S = np.array(NEAR_TIE_SCORES), np.array(NEAR_TIE_MEMBERS, dtype=np.int64)
    with pytest.raises(ValueError, match="gsea_running"):
        _call(monkeypatch, s, S, es=es)


def test_gsea_running_still_raises_off_a_tie(monkeypatch: pytest.MonkeyPatch) -> None:
    """Example 1 (ES 0.8, no tie) with the opposite sign in the table raises."""
    s, S = np.array([2.0, 1.0, 0.5, -0.5, -1.0, -2.0]), np.array([1, 3], dtype=np.int64)
    assert _call(monkeypatch, s, S) is not None
    with pytest.raises(ValueError, match="gsea_running"):
        _call(monkeypatch, s, S, es=-0.8)


@st.composite
def _discrete(draw: st.DrawFn) -> tuple[np.ndarray, np.ndarray]:
    n = draw(st.integers(4, 30))
    raw = np.array(draw(st.lists(st.sampled_from([-1.0, 0.0, 1.0, 0.1, -0.2, 0.3]),
                                 min_size=n, max_size=n)))  # fmt: skip
    s = np.sort(raw - np.median(raw))[::-1].copy()
    m = draw(st.integers(1, n - 1))
    S = np.array(sorted(draw(st.sets(st.integers(1, n), min_size=m, max_size=m))), np.int64)
    return s, S


@settings(max_examples=300, deadline=None)
@given(_discrete())
@example((np.array(NEAR_TIE_SCORES), np.array(NEAR_TIE_MEMBERS, dtype=np.int64)))
def test_property_gsea_running_never_rejects_the_engine_es(
    case: tuple[np.ndarray, np.ndarray],
) -> None:
    """Invariant: on any ranking with discrete (tied) centred scores, the step's check accepts
    the ES the engine computed on that same ranking."""
    s, S = case
    with pytest.MonkeyPatch.context() as mp:
        out = _call(mp, s, S)
    assert out is not None and len(out) == len(s)


def test_function_step_writes_support_at_a_near_tie(
    world: World, monkeypatch: pytest.MonkeyPatch
) -> None:
    """On the step: every top set's walk is forced onto the maxP = -minP tie (an extra dip to
    -maxP, the side the engine did not take for a positive ES); the step neither raises nor
    drops the file."""
    import gippa.function.running as R

    real = R.running_es
    tied: list[bool] = []

    def tie_walk(stats: np.ndarray, members: np.ndarray) -> np.ndarray:
        w = np.array(real(stats, members), dtype=np.float64)
        hi, lo = max(0.0, float(w.max())), min(0.0, float(w.min()))
        if hi > -lo:
            w[int(np.argmin(w))] = -hi
            tied.append(True)
        return w

    monkeypatch.setattr(R, "running_es", tie_walk)
    out = world.run("near_tie")
    assert tied, "vacuous: no top set with a positive ES"
    rows = sa.read(sa.need(out[KEY], NAME))
    assert len(_sets(rows)) == sa.GSEA_TOP


def _publish_twice(root: Path, first: int, second: int) -> tuple[Path, object]:
    from gippa.io.plot_support import write_support

    out = root / "OUT"
    for tag, n in (("run1", first), ("run2", second)):
        d = root / f"work_{tag}" / KEY
        d.mkdir(parents=True)
        write_support(NAME, pd.DataFrame({"tag": [tag] * n, "x": range(n)}), d)
        view = sa.publish_and_view(root / f"work_{tag}", out)
    return out / "plots" / "support", view


@pytest.mark.parametrize(("first", "second"), [(10, 200_000), (200_000, 10)],
                         ids=["tsv_then_gz", "gz_then_tsv"])  # fmt: skip
def test_publish_twice_leaves_one_file_and_run_view_reads_the_new_one(
    tmp_path: Path, first: int, second: int
) -> None:
    """Review B1b: two runs into one OUT whose support file flips suffix (.tsv <-> .tsv.gz)
    leave exactly one file, and ``RunView.support`` reads the second run's."""
    support, view = _publish_twice(tmp_path, first, second)
    files = sorted(p.name for p in support.iterdir() if p.name.startswith(NAME))
    assert len(files) == 1, files
    assert files[0] == (f"{NAME}.tsv.gz" if second > first else f"{NAME}.tsv")
    got = view.support(NAME)  # type: ignore[attr-defined]
    assert got is not None and set(got["tag"]) == {"run2"} and len(got) == second
