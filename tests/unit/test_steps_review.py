"""Lead decisions after code review (docs/specs/M2-02-step-registry-precomputed.md, items 1-9).

Assumptions pinned here (the spec leaves the exact mechanics open):
- The step-module loader discovers modules with ``pkgutil``-style iteration over
  ``gippa.steps.__path__`` (so a test can point that path at a temporary directory).
- ``run_step`` hands the step temporary output paths and renames them into place on success.
- Opt-in step ``bactermfinder`` is enabled by ``annotation.terminators.bactermfinder``.
- Precomputed parse errors count rows from the first data row (header excluded) as row 1.
"""

from __future__ import annotations

import itertools
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import pytest

from gippa.config import GippaConfig, resolve_config
from gippa.exitcodes import ExitCode
from gippa.io.tsv import write_tsv
from gippa.model.coords import Replicon


def _exit_code(argv: list[str]) -> int:
    from gippa.cli import main

    with pytest.raises(SystemExit) as info:
        main(argv)
    code = info.value.code
    return 0 if code is None else int(code)


@contextmanager
def _fresh() -> Iterator[ModuleType]:
    from gippa.steps import registry as reg

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(reg, "_REGISTRY", {})
        yield reg


@pytest.fixture
def genome(tmp_path: Path) -> Path:
    path = tmp_path / "g.fa"
    path.write_text(">chr1 [topology=circular]\n" + "ACGT" * 50 + "\n")
    return path


@pytest.fixture
def no_step_modules(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    """Point the loader at an empty package directory (no real step modules get imported)."""
    import gippa.steps

    d = tmp_path / "steps_pkg"
    d.mkdir()
    monkeypatch.setattr(gippa.steps, "__path__", [str(d)])
    return d


@pytest.fixture(autouse=True)
def _test_db(tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch) -> None:
    """``gippa run`` runs preflight: a test ``$GIPPA_DB_DIR``, never the real one (M3-17 rev.
    12)."""
    from helpers.mapping_steps import write_test_db

    monkeypatch.setenv("GIPPA_DB_DIR", str(write_test_db(tmp_path_factory.mktemp("test_db"))))


def _run_argv(genome: Path, outdir: Path, *extra: str) -> list[str]:
    """``gippa run --executor local`` with one peptide as evidence (a run without evidence exits 1,
    M3-17 review decision 7)."""
    pep = genome.parent / "evidence_peptides.tsv"
    pep.write_text("peptide\tprotein_ids\nPEPTIDEK\tP1\n")
    return ["run", "--genome", str(genome), "--peptides", str(pep), "--outdir", str(outdir),
            "--executor", "local", *extra]  # fmt: skip


# ---------------------------------------------------------------- 1. loader; empty registry


def test_run_local_with_no_registered_steps_exits_5(
    registry: ModuleType,
    no_step_modules: Path,
    genome: Path,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    code = _exit_code(_run_argv(genome, tmp_path / "o"))
    assert code == ExitCode.INTERNAL == 5
    assert "not implemented" in capsys.readouterr().err


def test_run_local_loader_imports_step_modules(
    registry: ModuleType,
    no_step_modules: Path,
    genome: Path,
    tmp_path: Path,
) -> None:
    modname = "gippa.steps.zz_review_toy"
    (no_step_modules / "zz_review_toy.py").write_text(
        "from gippa.steps.registry import step\n\n"
        "@step('zz_review_toy', inputs=[], outputs=['zz_review_toy.txt'])\n"
        "def _toy(ctx):\n"
        "    ctx.outputs['zz_review_toy.txt'].write_text('loaded')\n"
    )
    try:
        code = _exit_code(_run_argv(genome, tmp_path / "o"))
        assert modname in sys.modules
        assert "zz_review_toy" in [s.name for s in registry.all_steps()]
        assert code == ExitCode.SUCCESS
        found = list((tmp_path / "o").rglob("zz_review_toy.txt"))
        assert found and found[0].read_text() == "loaded"
    finally:
        sys.modules.pop(modname, None)


# ---------------------------------------------------------------- 2. AMR rows with NA coords

REPLICONS = {"chr1": Replicon("chr1", 1000, "circular", "chromosome")}
REGION_COLS = [
    "region_id",
    "region_type",
    "contig",
    "start",
    "end",
    "source_tool",
    "score",
    "genes",
    "notes",
]


def _load(key: str, path: Path) -> pd.DataFrame:
    from gippa.annotate.precomputed import load_precomputed

    return load_precomputed(key, path, REPLICONS)


def test_regions_amr_row_with_na_start_end_loads(tmp_path: Path) -> None:
    path = tmp_path / "regions.tsv"
    rows = [
        ["R0001", "prophage", "chr1", 11, 400, "genomad", 0.97, ["G_0001"], None],
        ["R0002", "amr", "chr1", None, None, "amrfinderplus", None, ["G_0005"], "blaTEM-1"],
    ]
    write_tsv(path, REGION_COLS, rows)
    df = _load("regions", path)
    assert df["region_type"].tolist() == ["prophage", "amr"]
    assert int(df["start"].iloc[0]) == 10 and int(df["end"].iloc[0]) == 400
    assert pd.isna(df["start"].iloc[1]) and pd.isna(df["end"].iloc[1])
    assert df["genes"].iloc[1] == ["G_0005"]


def test_regions_non_amr_row_with_na_start_still_rejected(tmp_path: Path) -> None:
    path = tmp_path / "regions.tsv"
    write_tsv(path, REGION_COLS, [["R1", "prophage", "chr1", None, 40, "genomad", 1.0, [], None]])
    with pytest.raises(ValueError, match="start"):
        _load("regions", path)


# ---------------------------------------------------------------- 3. stale outputs; atomic writes


def test_noop_step_with_stale_output_fails_output_check(
    registry: ModuleType, tmp_path: Path
) -> None:
    registry.step("toy", inputs=[], outputs=["b"])(lambda ctx: None)
    out = tmp_path / "b.txt"
    out.write_text("stale from a previous run")
    with pytest.raises(registry.WiringError, match="toy"):
        registry.run_step("toy", GippaConfig(), {}, {"b": out})
    assert not out.exists()


def test_output_removed_before_step_runs_and_written_via_temp_path(
    registry: ModuleType, tmp_path: Path
) -> None:
    out = tmp_path / "b.txt"
    out.write_text("stale")
    seen: dict[str, Any] = {}

    def toy(ctx: Any) -> None:
        seen["final_exists"] = out.exists()
        seen["path"] = ctx.outputs["b"]
        ctx.outputs["b"].write_text("fresh")

    registry.step("toy", inputs=[], outputs=["b"])(toy)
    registry.run_step("toy", GippaConfig(), {}, {"b": out})
    assert seen["final_exists"] is False
    assert seen["path"] != out  # temporary path, renamed into place on success
    assert out.read_text() == "fresh"
    assert not Path(seen["path"]).exists()


def test_failing_step_leaves_no_partial_output(registry: ModuleType, tmp_path: Path) -> None:
    def toy(ctx: Any) -> None:
        ctx.outputs["b"].write_text("partial")
        raise RuntimeError("boom")

    registry.step("toy", inputs=[], outputs=["b"])(toy)
    out = tmp_path / "b.txt"
    out.write_text("stale")
    with pytest.raises(registry.StepError):
        registry.run_step("toy", GippaConfig(), {}, {"b": out})
    assert not out.exists()
    assert [p.name for p in tmp_path.iterdir() if "partial" in p.read_text()] == []


# ---------------------------------------------------------------- 4. opt-in steps follow config


def _opt_in_graph(reg: ModuleType) -> None:
    def write(key: str) -> Any:
        def fn(ctx: Any) -> None:
            ctx.outputs[key].write_text(key)

        return fn

    reg.step("bactermfinder", inputs=[], outputs=["bactermfinder.txt"], module="opt_in")(
        write("bactermfinder.txt")
    )
    reg.step(
        "terminator_layer",
        inputs=[],
        optional_inputs=["bactermfinder.txt"],
        outputs=["terminators.txt"],
    )(write("terminators.txt"))


def test_opt_in_step_not_planned_when_config_off_even_if_targeted(registry: ModuleType) -> None:
    from gippa.runner.local import plan_local

    _opt_in_graph(registry)
    cfg = resolve_config(None, ["annotation.terminators.bactermfinder=false"])
    assert plan_local(cfg, ["bactermfinder", "terminator_layer"]) == ["terminator_layer"]


def test_opt_in_step_planned_when_config_on(registry: ModuleType) -> None:
    from gippa.runner.local import plan_local

    _opt_in_graph(registry)
    cfg = resolve_config(None, ["annotation.terminators.bactermfinder=true"])
    assert plan_local(cfg, ["terminator_layer"]) == ["bactermfinder", "terminator_layer"]


@pytest.mark.parametrize("flag", [False, True])
def test_cli_add_bactermfinder_controls_opt_in_step(
    flag: bool,
    registry: ModuleType,
    no_step_modules: Path,
    genome: Path,
    tmp_path: Path,
) -> None:
    _opt_in_graph(registry)
    outdir = tmp_path / "o"
    extra = ["--add-bactermfinder"] if flag else []
    assert _exit_code(_run_argv(genome, outdir, *extra)) == ExitCode.SUCCESS
    assert list(outdir.rglob("terminators.txt"))
    assert bool(list(outdir.rglob("bactermfinder.txt"))) is flag


# ---------------------------------------------------------------- 5. precomputed vs stale key


def test_precomputed_layer_wins_over_stale_workdir_key(
    registry: ModuleType, tmp_path: Path
) -> None:
    from gippa.runner.local import run_local

    registry.external("term_tool", inputs=[], outputs=["term_tool.tsv"], precomputed="terminators")
    registry.step(
        "term_layer",
        inputs=["term_tool.tsv"],
        outputs=["terminators.parquet"],
        precomputed="terminators",
    )(lambda ctx: ctx.outputs["terminators.parquet"].write_text("computed"))

    def use(ctx: Any) -> None:
        ctx.outputs["use.txt"].write_text(ctx.inputs["terminators.parquet"].read_text())

    registry.step("use", inputs=["terminators.parquet"], outputs=["use.txt"])(use)
    layer = tmp_path / "terminators.tsv"
    layer.write_text("contig\tstart\tend\tstrand\tconfidence\tsource\n")
    work = tmp_path / "work"
    work.mkdir()
    (work / "terminators.parquet").write_text("STALE from an earlier run")
    cfg = resolve_config(None, [f"annotation.precomputed.terminators={layer}"])
    run_local(cfg, work, ["use"])
    assert (work / "use.txt").read_text() == layer.read_text()


# ---------------------------------------------------------------- 6. duplicate output keys


def test_two_py_steps_with_same_output_key_rejected(registry: ModuleType) -> None:
    registry.step("a", inputs=[], outputs=["x.parquet"])(lambda ctx: None)
    with pytest.raises(ValueError, match=r"x\.parquet"):
        registry.step("b", inputs=[], outputs=["y", "x.parquet"])(lambda ctx: None)
    assert [s.name for s in registry.all_steps()] == ["a"]


def test_py_and_ext_step_with_same_output_key_rejected(registry: ModuleType) -> None:
    registry.external("tool", inputs=[], outputs=["out/"])
    with pytest.raises(ValueError, match="out/"):
        registry.step("parse", inputs=[], outputs=["out/"])(lambda ctx: None)


# ---------------------------------------------------------------- 7. failing step log


def test_step_failure_writes_log_with_traceback(
    registry: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    import yaml

    def toy(ctx: Any) -> None:
        raise RuntimeError("boom-in-toy")

    registry.step("toy", inputs=[], outputs=["b"])(toy)
    resolved = tmp_path / "resolved.yaml"
    resolved.write_text(yaml.safe_dump(GippaConfig().model_dump(mode="json"), sort_keys=True))
    monkeypatch.chdir(tmp_path)
    argv = ["_step", "toy", "--config", str(resolved), "--out", f"b={tmp_path / 'b.txt'}"]
    assert _exit_code(argv) == ExitCode.STEP_FAILED
    assert "logs/toy.log" in capsys.readouterr().err
    log = (tmp_path / "logs" / "toy.log").read_text()
    assert "Traceback" in log
    assert "RuntimeError: boom-in-toy" in log


def test_run_local_step_failure_writes_log_and_points_to_it(
    registry: ModuleType,
    no_step_modules: Path,
    genome: Path,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    def toy(ctx: Any) -> None:
        raise RuntimeError("boom-local")

    registry.step("toy", inputs=[], outputs=["b"])(toy)
    outdir = tmp_path / "o"
    assert _exit_code(_run_argv(genome, outdir)) == ExitCode.STEP_FAILED
    assert "toy.log" in capsys.readouterr().err
    logs = [p for p in outdir.rglob("toy.log") if p.parent.name == "logs"]
    assert logs, "no logs/toy.log under outdir"
    text = logs[0].read_text()
    assert "Traceback" in text and "boom-local" in text


# ---------------------------------------------------------------- 8. repeated --in -> list


def test_repeated_in_key_collects_list(registry: ModuleType, tmp_path: Path) -> None:
    import yaml

    seen: dict[str, Any] = {}

    def toy(ctx: Any) -> None:
        seen["chunks"] = ctx.inputs["chunks"]
        seen["one"] = ctx.inputs["one"]
        ctx.outputs["b"].write_text("".join(p.read_text() for p in ctx.inputs["chunks"]))

    registry.step("toy", inputs=["chunks", "one"], outputs=["b"])(toy)
    paths = []
    for i in range(3):
        p = tmp_path / f"chunk{i}.txt"
        p.write_text(str(i))
        paths.append(p)
    one = tmp_path / "one.txt"
    one.write_text("1")
    resolved = tmp_path / "resolved.yaml"
    resolved.write_text(yaml.safe_dump(GippaConfig().model_dump(mode="json"), sort_keys=True))
    argv = ["_step", "toy", "--config", str(resolved), "--in", f"one={one}"]
    for p in paths:
        argv += ["--in", f"chunks={p}"]
    argv += ["--out", f"b={tmp_path / 'b.txt'}"]
    assert _exit_code(argv) == ExitCode.SUCCESS
    assert list(seen["chunks"]) == paths
    assert isinstance(seen["one"], Path)
    assert (tmp_path / "b.txt").read_text() == "012"


# ---------------------------------------------------------------- 9. parse errors name key/col/row

TERM_COLS = ["contig", "start", "end", "strand", "confidence", "source"]


@pytest.mark.parametrize(
    ("bad_row", "column"),
    [
        (["chr1", 5, 40, "+", "high", "X"], "confidence"),
        (["chr1", 5, 40, "?", 0.5, "X"], "strand"),
        (["chrZ", 5, 40, "+", 0.5, "X"], "contig"),
        (["chr1", None, 40, "+", 0.5, "X"], "start"),
        (["chr1", "five", 40, "+", 0.5, "X"], "start"),
    ],
)
def test_precomputed_parse_error_names_key_column_and_row(
    bad_row: list[Any], column: str, tmp_path: Path
) -> None:
    path = tmp_path / "terminators.tsv"
    good = ["chr1", 100, 160, "-", 0.82, "TerminatorNet"]
    write_tsv(path, TERM_COLS, [good, good, bad_row])
    with pytest.raises(ValueError) as info:
        _load("terminators", path)
    msg = str(info.value)
    assert "precomputed terminators" in msg
    assert column in msg
    assert "row 3" in msg


# ---------------------------------------------------------------- determinism of run_local


def _chain(reg: ModuleType, order: list[str]) -> None:
    def x(ctx: Any) -> None:
        ctx.outputs["b.txt"].write_text(repr(ctx.rng.integers(0, 10**9, size=5).tolist()))

    def y(ctx: Any) -> None:
        ctx.outputs["c.txt"].write_text(ctx.inputs["b.txt"].read_text() + f"|{ctx.rng.random()!r}")

    def z(ctx: Any) -> None:
        ctx.outputs["d.txt"].write_text(f"z{ctx.rng.random()!r}")

    def w(ctx: Any) -> None:
        text = ctx.inputs["c.txt"].read_text() + ctx.inputs["d.txt"].read_text()
        ctx.outputs["e.txt"].write_text(text)

    defs = {
        "x": ([], ["b.txt"], x),
        "y": (["b.txt"], ["c.txt"], y),
        "z": ([], ["d.txt"], z),
        "w": (["c.txt", "d.txt"], ["e.txt"], w),
    }
    for name in order:
        ins, outs, fn = defs[name]

        def rec(ctx: Any, _fn: Any = fn, _name: str = name) -> None:
            ctx.log.info("run %s", _name)
            _fn(ctx)

        reg.step(name, inputs=ins, outputs=outs)(rec)


KEYS = ("b.txt", "c.txt", "d.txt", "e.txt")


def test_run_local_rerun_same_workdir_is_byte_identical(
    registry: ModuleType, tmp_path: Path
) -> None:
    from gippa.runner.local import run_local

    _chain(registry, ["x", "y", "z", "w"])
    work = tmp_path / "work"
    run_local(GippaConfig(seed=9), work, ["w"])
    first = {k: (work / k).read_bytes() for k in KEYS}
    run_local(GippaConfig(seed=9), work, ["w"])
    assert {k: (work / k).read_bytes() for k in KEYS} == first


def test_run_local_independent_of_registration_order(tmp_path: Path) -> None:
    from gippa.runner import local

    results: list[tuple[list[str], dict[str, bytes]]] = []
    for i, order in enumerate(itertools.permutations(["x", "y", "z", "w"])):
        with _fresh() as reg:
            _chain(reg, list(order))
            calls: list[str] = []
            real = reg.run_step

            def spy(name: str, *a: Any, _real: Any = real, _calls: list[str] = calls) -> None:
                _calls.append(name)
                _real(name, *a)

            with pytest.MonkeyPatch.context() as mp:
                mp.setattr(reg, "run_step", spy)
                work = tmp_path / f"w{i}"
                local.run_local(GippaConfig(seed=4), work, ["w"])
            results.append((calls, {k: (work / k).read_bytes() for k in KEYS}))
    first_calls, first_bytes = results[0]
    assert first_calls == ["x", "y", "z", "w"]
    for calls, data in results[1:]:
        assert calls == first_calls
        assert data == first_bytes
