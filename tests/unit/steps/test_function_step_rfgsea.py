"""M8-17b T9, A8 and A11: the ``function`` step on the R fgsea engine, without R
(docs/specs/M8-17-r-fgsea-engine.md §4 "Provenance", §9, examples 7-9; lead decisions 2 and 5).

Interface pinned here (test-author reading of §9):

- ``function_step`` builds ``r_fgsea_engine(resolve_launcher(os.environ), tmp_root=Path.cwd())``
  lazily: only when ``gsea_table`` would call it. A run without abundance, or without a testable
  set, never resolves the launcher, so it succeeds even when ``GIPPA_FGSEA_PREFIX`` points
  nowhere.
- An engine error propagates: ``gippa _step function`` exits 4 (``ExitCode.STEP_FAILED``) and no
  ``tables/enrichment_*.tsv`` is written.
- Each engine call logs one §6.7 event ``function.gsea_engine`` with data
  ``{"call": "main" | "retest", "jobs": n, "r": ..., "fgsea": "1.36.2", "launcher": "prefix",
  "warnings": [...]}``.
- ``run_function(..., *, engine=SIMPLE_ENGINE)``.

The R side is the fake front end of ``helpers.fgsea_r_data`` (every sent set gets ES 0.5, p 0.01,
so GO parents with tested descendants are significant and a re-test call follows). Tests never
read the real ``$GIPPA_ROOT``: it is pointed at an empty directory.
"""

from __future__ import annotations

import inspect
import logging
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from click.testing import CliRunner
from helpers.fgsea_r_data import fake_prefix
from helpers.fgsea_switch_data import RecordingEngine, SimWorld, sim_cfg, sim_world
from helpers.function_db import FUNCTION_RESOURCES
from helpers.function_step_data import oracle_tables, write_gene_table, write_tu_folds

from gippa.exitcodes import ExitCode

OUTPUTS = {
    "category_summary": "tables/category_summary.tsv",
    "enrichment_detection": "tables/enrichment_detection.tsv",
    "enrichment_gsea": "tables/enrichment_gsea.tsv",
}
EVENT = "function.gsea_engine"
EVENT_KEYS = {"call", "jobs", "r", "fgsea", "launcher", "warnings"}


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> SimWorld:
    return sim_world(tmp_path_factory.mktemp("rfgsea_step"))


@pytest.fixture(scope="module")
def folds(world: SimWorld) -> Path:
    return write_tu_folds(world.gene_table.parent / "tu_folds_empty.parquet", {})


@pytest.fixture
def isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """cwd = a fresh directory; ``$GIPPA_ROOT`` empty; no inherited fgsea variables."""
    root = tmp_path / "empty_gippa_root"
    root.mkdir()
    monkeypatch.setenv("GIPPA_ROOT", str(root))
    monkeypatch.delenv("GIPPA_FGSEA_PREFIX", raising=False)
    monkeypatch.delenv("GIPPA_FGSEA_IMAGE", raising=False)
    cwd = tmp_path / "task"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    return cwd


def _inputs(world: SimWorld, folds: Path, gene_table: Path | None = None) -> dict[str, Path]:
    ins = {
        "gene_table": gene_table or world.gene_table,
        "replicons": world.replicons_path,
        "tu_folds": folds,
    }
    return ins | dict(world.db.inputs(FUNCTION_RESOURCES))


def _no_abundance(world: SimWorld, where: Path) -> Path:
    from gippa.io.genes_tsv import read_genes_parquet

    genes = read_genes_parquet(world.gene_table)
    genes["abundance_log2"] = float("nan")
    genes["abundance_scale"] = None
    genes["abundance_type"] = None
    return write_gene_table(where / "noab" / "genes.parquet", genes)


def _cli(ins: dict[str, Path], sets: tuple[str, ...] = ()) -> Any:
    from gippa.cli import cli

    argv = ["_step", "function"]
    for k, v in ins.items():
        argv += ["--in", f"{k}={v}"]
    for k, v in OUTPUTS.items():
        argv += ["--out", f"{k}={v}"]
    for s in (
        "seed=7",
        "stats.bootstrap=100",
        "function.gsea.n_perm_simple=200",
        *sets,
    ):
        argv += ["--set", s]
    return CliRunner().invoke(cli, argv)


class _Recorder(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


def _run_step(ins: dict[str, Path], out: Path, cfg: Any = None) -> list[dict[str, Any]]:
    """``run_step("function")``; returns the data of every ``function.gsea_engine`` event."""
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    rec = _Recorder()
    log = logging.getLogger("gippa.steps.function")
    old = log.level
    log.addHandler(rec)
    log.setLevel(logging.DEBUG)
    try:
        reg.run_step(
            "function", cfg or sim_cfg(), dict(ins), {k: out / v for k, v in OUTPUTS.items()}
        )
    finally:
        log.removeHandler(rec)
        log.setLevel(old)
    return [r.data for r in rec.records if getattr(r, "event", None) == EVENT]  # type: ignore[attr-defined]


def _gsea(out: Path) -> pd.DataFrame:
    return pd.read_csv(out / OUTPUTS["enrichment_gsea"], sep="\t", dtype=str, keep_default_na=False)


# ---------------------------------------------------------------- A8: failures


def test_missing_prefix_with_abundance_exits_4_without_tables(
    world: SimWorld, folds: Path, isolated: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GIPPA_FGSEA_PREFIX", str(isolated.parent / "nowhere"))
    res = _cli(_inputs(world, folds))
    assert res.exit_code == ExitCode.STEP_FAILED == 4, res.output
    for k in ("enrichment_detection", "enrichment_gsea", "category_summary"):
        assert not (isolated / OUTPUTS[k]).exists(), k
    log = (isolated / "logs" / "function.log").read_text(encoding="utf-8")
    assert "FgseaEngineUnavailable" in log


def test_unbuilt_env_under_gippa_root_exits_4(world: SimWorld, folds: Path, isolated: Path) -> None:
    """No ``GIPPA_FGSEA_PREFIX`` and an empty ``$GIPPA_ROOT``: launcher rule 4, exit 4."""
    res = _cli(_inputs(world, folds))
    assert res.exit_code == 4, res.output
    log = (isolated / "logs" / "function.log").read_text(encoding="utf-8")
    assert "sbatch scripts/slurm/build_envs.sbatch fgsea" in log
    assert not (isolated / OUTPUTS["enrichment_gsea"]).exists()


def test_r_failure_exits_4_without_tables(
    world: SimWorld, folds: Path, isolated: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Example 9: the script exits 1 with an fgseaMultilevel error."""
    prefix = fake_prefix(
        isolated.parent / "env",
        {"mode": "fail", "code": 1, "stderr": "Error in fgseaMultilevel(...) : x\n"},
    )
    monkeypatch.setenv("GIPPA_FGSEA_PREFIX", str(prefix))
    res = _cli(_inputs(world, folds))
    assert res.exit_code == 4, res.output
    log = (isolated / "logs" / "function.log").read_text(encoding="utf-8")
    assert "R fgsea failed (exit 1" in log and "Error in fgseaMultilevel" in log
    assert not list(isolated.glob("tables/enrichment_*.tsv"))
    assert len((prefix / "calls.txt").read_text().splitlines()) == 1  # never retried
    assert not list(isolated.glob("gippa-fgsea-*")), "temporary call directory left behind"


# ---------------------------------------------------------------- A8: no R needed


def test_no_abundance_runs_without_r(
    world: SimWorld, folds: Path, isolated: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Example 7: zero engine calls, the engine is never constructed."""
    monkeypatch.setenv("GIPPA_FGSEA_PREFIX", str(isolated.parent / "nowhere"))
    res = _cli(_inputs(world, folds, _no_abundance(world, isolated.parent)))
    assert res.exit_code == 0, res.output
    gsea = (isolated / OUTPUTS["enrichment_gsea"]).read_text(encoding="utf-8")
    assert len(gsea.splitlines()) == 1


def test_no_abundance_never_resolves_the_launcher(
    world: SimWorld, folds: Path, isolated: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import gippa.io.fgsea_r as fr
    import gippa.steps.function as step_mod

    def boom(*_: Any, **__: Any) -> Any:
        raise AssertionError("the R engine was built without abundance")

    for name in ("resolve_launcher", "r_fgsea_engine"):
        monkeypatch.setattr(fr, name, boom)
        monkeypatch.setattr(step_mod, name, boom, raising=False)
    events = _run_step(_inputs(world, folds, _no_abundance(world, isolated.parent)), isolated)
    assert events == []


def test_no_testable_set_runs_without_r(
    world: SimWorld, folds: Path, isolated: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Example 8 through the step: every set is below ``function.gsea_size[0]``."""
    monkeypatch.setenv("GIPPA_FGSEA_PREFIX", str(isolated.parent / "nowhere"))
    res = _cli(_inputs(world, folds), sets=("function.gsea_size=[400,500]",))
    assert res.exit_code == 0, res.output
    assert len((isolated / OUTPUTS["enrichment_gsea"]).read_text().splitlines()) == 1


# ---------------------------------------------------------------- events and the R engine


def test_one_event_per_engine_call(
    world: SimWorld, folds: Path, isolated: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """With ``gsea_size`` [20, 500] go_bp's GO:0990001 (19 TUs) is untested, so its parent
    GO:0990100 minus its significant child GO:0990000 keeps TUs: a re-test call follows."""
    from gippa.config import GippaConfig

    prefix = fake_prefix(isolated.parent / "env")
    monkeypatch.setenv("GIPPA_FGSEA_PREFIX", str(prefix))
    cfg = GippaConfig.model_validate(
        sim_cfg().model_dump()
        | {"function": {"gsea_size": [20, 500], "gsea": {"n_perm_simple": 200}}}
    )
    events = _run_step(_inputs(world, folds), isolated, cfg)
    calls = (prefix / "calls.txt").read_text().splitlines()
    assert len(calls) == 2, "expected a main and a re-test call"
    assert len(events) == len(calls)
    assert [e["call"] for e in events] == ["main", "retest"]
    for e in events:
        assert set(e) >= EVENT_KEYS, e
        assert e["fgsea"] == "1.36.2"
        assert e["launcher"] == "prefix"
        assert str(e["r"]).startswith("R version 4.5.3")
        assert isinstance(e["jobs"], int) and e["jobs"] >= 1
        assert isinstance(e["warnings"], list)
    # the main call has one job per collection with an in-range set (six fixture collections)
    assert 1 <= events[0]["jobs"] <= 6


def test_event_carries_the_r_warnings(
    world: SimWorld, folds: Path, isolated: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    msg = "There are ties in the preranked stats (3.1% of the list)."
    prefix = fake_prefix(isolated.parent / "env", {"warnings": [["main:cog", msg]]})
    monkeypatch.setenv("GIPPA_FGSEA_PREFIX", str(prefix))
    events = _run_step(_inputs(world, folds), isolated)
    assert msg in str(events[0]["warnings"])


def test_r_engine_table_and_tmp_root(
    world: SimWorld, folds: Path, isolated: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prefix = fake_prefix(isolated.parent / "env")
    monkeypatch.setenv("GIPPA_FGSEA_PREFIX", str(prefix))
    _run_step(_inputs(world, folds), isolated)
    df = _gsea(isolated)
    assert len(df) > 0
    assert set(df["method"]) == {"fgsea_multilevel"}
    assert set(df["n_perm"]) == {"200"}
    for d in (prefix / "calls.txt").read_text().splitlines():
        assert Path(d).parent == isolated, "tmp_root must be the task directory (Path.cwd())"
    assert not list(isolated.glob("gippa-fgsea-*"))


def test_r_engine_step_is_byte_identical(
    world: SimWorld, folds: Path, isolated: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prefix = fake_prefix(isolated.parent / "env")
    monkeypatch.setenv("GIPPA_FGSEA_PREFIX", str(prefix))
    a, b = isolated / "a", isolated / "b"
    _run_step(_inputs(world, folds), a)
    _run_step(_inputs(world, folds), b)
    for k in OUTPUTS.values():
        assert (a / k).read_bytes() == (b / k).read_bytes(), k


def test_a11_detection_table_unchanged_by_the_engine(
    world: SimWorld, folds: Path, isolated: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A11: the ORA path is untouched; ``enrichment_detection.tsv`` equals the M8-06 chain."""
    prefix = fake_prefix(isolated.parent / "env")
    monkeypatch.setenv("GIPPA_FGSEA_PREFIX", str(prefix))
    _run_step(_inputs(world, folds), isolated)
    paths, _ = oracle_tables(
        world.gene_table,
        world.replicons_path,
        None,
        world.db.inputs(FUNCTION_RESOURCES),
        sim_cfg(),
        isolated.parent / "oracle",
    )
    got = (isolated / OUTPUTS["enrichment_detection"]).read_bytes()
    assert got == paths["enrichment_detection"].read_bytes()
    got = (isolated / OUTPUTS["category_summary"]).read_bytes()
    assert got == paths["category_summary"].read_bytes()


# ---------------------------------------------------------------- run_function


def test_run_function_engine_keyword_defaults_to_simple() -> None:
    from gippa.function.gsea import SIMPLE_ENGINE
    from gippa.steps.function import run_function

    p = inspect.signature(run_function).parameters["engine"]
    assert p.kind is inspect.Parameter.KEYWORD_ONLY
    assert p.default is SIMPLE_ENGINE


def test_run_function_passes_the_engine(world: SimWorld) -> None:
    from gippa.io.function_resources import read_resources
    from gippa.steps.function import run_function

    eng = RecordingEngine()
    res = read_resources(dict(world.db.inputs(FUNCTION_RESOURCES)))
    _, _, gsea, _ = run_function(world.genes, world.replicons, {}, res, sim_cfg(), engine=eng)
    assert len(eng.calls) >= 1
    assert len(gsea) > 0 and set(gsea["method"]) == {"fgsea_multilevel"}
    # the stream is for_step("function.gsea", seed) (§5)
    from gippa.rng import for_step

    rng = for_step("function.gsea", 7)
    rng.permutation(len(eng.calls[0].tu_ids))
    assert eng.calls[0].seed == int(rng.integers(1, 2**31 - 1))
