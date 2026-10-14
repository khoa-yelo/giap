"""Local executor (M2-02 example 8; GIAP_SPEC.md §6.4, §6.6).

Convention assumed here: ``run_local`` stores every step key at ``workdir / key`` (keys are the
work-directory file names of the M2-index contract table).
"""

from __future__ import annotations

from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from giap.config import GiapConfig, resolve_config
from giap.exitcodes import ExitCode

MSG = "annotation.precomputed.terminators is required with --executor local"


def _local() -> ModuleType:
    from giap.runner import local

    return local


def _chain(reg: ModuleType, calls: list[str]) -> None:
    """x: () -> b.txt ; y: b.txt -> c.txt ; registered in reverse order on purpose."""

    def y(ctx: Any) -> None:
        calls.append("y")
        ctx.outputs["c.txt"].write_text(ctx.inputs["b.txt"].read_text() + "|y")

    def x(ctx: Any) -> None:
        calls.append("x")
        ctx.outputs["b.txt"].write_text(f"x:{ctx.rng.integers(0, 10**9)}")

    reg.step("y", inputs=["b.txt"], outputs=["c.txt"])(y)
    reg.step("x", inputs=[], outputs=["b.txt"])(x)


def _terminator_graph(reg: ModuleType, calls: list[str]) -> None:
    def rec(name: str, out: str) -> Any:
        def fn(ctx: Any) -> None:
            calls.append(name)
            ctx.outputs[out].write_text(name)

        return fn

    reg.external("term_tool", inputs=[], outputs=["term_tool.tsv"], precomputed="terminators")
    reg.step(
        "term_layer",
        inputs=["term_tool.tsv"],
        outputs=["terminators.parquet"],
        precomputed="terminators",
    )(rec("term_layer", "terminators.parquet"))
    reg.step("first", inputs=[], outputs=["first.txt"])(rec("first", "first.txt"))
    reg.step("use", inputs=["first.txt", "terminators.parquet"], outputs=["use.txt"])(
        rec("use", "use.txt")
    )


def test_run_local_runs_python_steps_in_topological_order(
    registry: ModuleType, tmp_path: Path
) -> None:
    calls: list[str] = []
    _chain(registry, calls)
    _local().run_local(GiapConfig(seed=3), tmp_path / "work", ["y"])
    assert calls == ["x", "y"]
    assert (tmp_path / "work" / "c.txt").read_text().endswith("|y")


def test_plan_local_includes_upstream_py_steps(registry: ModuleType) -> None:
    _chain(registry, [])
    assert _local().plan_local(GiapConfig(), ["y"]) == ["x", "y"]


def test_plan_local_skips_steps_of_disabled_modules(registry: ModuleType) -> None:
    def noop(ctx: Any) -> None:
        return None

    registry.step("core_step", inputs=[], outputs=["core.txt"])(noop)
    registry.step("bgc_step", inputs=[], outputs=["bgc.txt"], module="bgc")(noop)
    cfg = GiapConfig(modules=["core"])
    assert _local().plan_local(cfg, ["core_step", "bgc_step"]) == ["core_step"]
    cfg_all = GiapConfig(modules=["core", "bgc"])
    assert _local().plan_local(cfg_all, ["core_step", "bgc_step"]) == ["bgc_step", "core_step"]


def test_plan_local_replaces_steps_with_precomputed_layer(
    registry: ModuleType, tmp_path: Path
) -> None:
    _terminator_graph(registry, [])
    layer = tmp_path / "terminators.tsv"
    layer.write_text("contig\tstart\tend\tstrand\tconfidence\tsource\n")
    cfg = resolve_config(None, [f"annotation.precomputed.terminators={layer}"])
    plan = _local().plan_local(cfg, ["use"])
    assert plan == ["first", "use"]
    assert "term_tool" not in plan and "term_layer" not in plan


def test_missing_precomputed_layer_raises_layer_missing(
    registry: ModuleType, tmp_path: Path
) -> None:
    calls: list[str] = []
    _terminator_graph(registry, calls)
    cfg = resolve_config(None, ["annotation.precomputed.terminators=null"])
    local = _local()
    with pytest.raises(local.LayerMissing) as info:
        local.run_local(cfg, tmp_path / "work", ["use"])
    assert info.value.key == "annotation.precomputed.terminators"
    assert MSG in str(info.value)
    assert calls == []  # validated before any step runs


def test_run_local_outputs_byte_identical(registry: ModuleType, tmp_path: Path) -> None:
    _chain(registry, [])
    local = _local()
    local.run_local(GiapConfig(seed=5), tmp_path / "w1", ["y"])
    local.run_local(GiapConfig(seed=5), tmp_path / "w2", ["y"])
    for key in ("b.txt", "c.txt"):
        assert (tmp_path / "w1" / key).read_bytes() == (tmp_path / "w2" / key).read_bytes()


def test_cli_run_local_missing_layer_exits_1(
    registry: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    import giap.cli
    from giap.cli import main

    # a registered step: with none, --executor local exits 5 (review decision 1)
    _terminator_graph(registry, [])
    local = _local()

    def boom(*args: Any, **kwargs: Any) -> Any:
        raise local.LayerMissing("annotation.precomputed.terminators")

    for name in ("plan_local", "run_local"):
        monkeypatch.setattr(local, name, boom)
        monkeypatch.setattr(giap.cli, name, boom, raising=False)
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1 [topology=circular]\n" + "ACGT" * 50 + "\n")
    argv = ["run", "--genome", str(genome), "--outdir", str(tmp_path / "o"), "--executor", "local"]
    with pytest.raises(SystemExit) as info:
        main(argv)
    assert info.value.code == ExitCode.INVALID_INPUT
    assert MSG in capsys.readouterr().err
