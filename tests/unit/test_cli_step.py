"""Hidden command ``gippa _step NAME --config resolved.yaml --in k=path --out k=path``.

M2-02 examples 5-7; GIPPA_SPEC.md §5.1 (exit 4 = a pipeline step failed, 5 = internal error).
The toy steps are registered in this process, so ``gippa.cli.main`` is called in-process.
"""

from __future__ import annotations

from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
import yaml

from gippa.config import GippaConfig
from gippa.exitcodes import ExitCode
from gippa.rng import for_step


def _exit_code(argv: list[str]) -> int:
    from gippa.cli import main

    with pytest.raises(SystemExit) as info:
        main(argv)
    code = info.value.code
    if code is None:
        return 0
    assert isinstance(code, int), f"non-integer exit code {code!r}"
    return code


@pytest.fixture
def resolved(tmp_path: Path) -> Path:
    path = tmp_path / "resolved.yaml"
    path.write_text(yaml.safe_dump(GippaConfig(seed=11).model_dump(mode="json"), sort_keys=True))
    return path


@pytest.fixture
def infile(tmp_path: Path) -> Path:
    path = tmp_path / "in.txt"
    path.write_text("payload\n")
    return path


def _argv(name: str, resolved: Path, ins: dict[str, Path], outs: dict[str, Path]) -> list[str]:
    argv = ["_step", name, "--config", str(resolved)]
    for k, p in ins.items():
        argv += ["--in", f"{k}={p}"]
    for k, p in outs.items():
        argv += ["--out", f"{k}={p}"]
    return argv


def _copy(ctx: Any) -> None:
    ctx.outputs["b"].write_bytes(ctx.inputs["a"].read_bytes())


def test_step_runs_and_writes_output(
    registry: ModuleType, resolved: Path, infile: Path, tmp_path: Path
) -> None:
    registry.step("toy", inputs=["a"], outputs=["b"])(_copy)
    out = tmp_path / "out" / "out.txt"
    assert _exit_code(_argv("toy", resolved, {"a": infile}, {"b": out})) == ExitCode.SUCCESS
    assert out.read_text() == "payload\n"


def test_step_loads_resolved_config_and_seeds_rng(
    registry: ModuleType, resolved: Path, tmp_path: Path
) -> None:
    def toy(ctx: Any) -> None:
        ctx.outputs["b"].write_text(f"{ctx.cfg.seed} {ctx.rng.random()!r}")

    registry.step("toy", inputs=[], outputs=["b"])(toy)
    out = tmp_path / "out.txt"
    assert _exit_code(_argv("toy", resolved, {}, {"b": out})) == ExitCode.SUCCESS
    assert out.read_text() == f"11 {for_step('toy', 11).random()!r}"


def test_step_repeatable_in_and_out(registry: ModuleType, resolved: Path, tmp_path: Path) -> None:
    def toy(ctx: Any) -> None:
        for k in ("c", "d"):
            text = ctx.inputs["a"].read_text() + ctx.inputs["b"].read_text()
            ctx.outputs[k].write_text(text)

    registry.step("toy", inputs=["a", "b"], outputs=["c", "d"])(toy)
    a, b = tmp_path / "a.txt", tmp_path / "b.txt"
    a.write_text("A")
    b.write_text("B")
    outs = {"c": tmp_path / "c.txt", "d": tmp_path / "d.txt"}
    assert _exit_code(_argv("toy", resolved, {"a": a, "b": b}, outs)) == ExitCode.SUCCESS
    assert outs["c"].read_text() == outs["d"].read_text() == "AB"


def test_step_optional_input_may_be_absent(
    registry: ModuleType, resolved: Path, infile: Path, tmp_path: Path
) -> None:
    def toy(ctx: Any) -> None:
        ctx.outputs["b"].write_text(",".join(sorted(ctx.inputs)))

    registry.step("toy", inputs=["a"], outputs=["b"], optional_inputs=["opt"])(toy)
    out = tmp_path / "o.txt"
    assert _exit_code(_argv("toy", resolved, {"a": infile}, {"b": out})) == ExitCode.SUCCESS
    assert out.read_text() == "a"


def test_step_directory_input(registry: ModuleType, resolved: Path, tmp_path: Path) -> None:
    def toy(ctx: Any) -> None:
        ctx.outputs["b"].write_text(str(ctx.inputs["genomad/"].is_dir()))

    registry.step("toy", inputs=["genomad/"], outputs=["b"])(toy)
    d = tmp_path / "genomad"
    d.mkdir()
    out = tmp_path / "o.txt"
    assert _exit_code(_argv("toy", resolved, {"genomad/": d}, {"b": out})) == ExitCode.SUCCESS
    assert out.read_text() == "True"


def test_failing_step_exits_4(
    registry: ModuleType,
    resolved: Path,
    infile: Path,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    def toy(ctx: Any) -> None:
        raise RuntimeError("boom")

    registry.step("toy", inputs=["a"], outputs=["b"])(toy)
    code = _exit_code(_argv("toy", resolved, {"a": infile}, {"b": tmp_path / "o.txt"}))
    assert code == ExitCode.STEP_FAILED == 4
    err = capsys.readouterr().err
    assert "step 'toy' failed" in err
    assert "boom" in err
    assert "logs/toy.log" in err


def test_unknown_step_exits_5(
    registry: ModuleType,
    resolved: Path,
    infile: Path,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    code = _exit_code(_argv("nope", resolved, {"a": infile}, {"b": tmp_path / "o.txt"}))
    assert code == ExitCode.INTERNAL == 5
    err = capsys.readouterr().err
    assert "nope" in err
    assert "not implemented" not in err.lower()


def test_undeclared_input_key_exits_5(
    registry: ModuleType,
    resolved: Path,
    infile: Path,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    registry.step("toy", inputs=["a"], outputs=["b"])(_copy)
    argv = _argv("toy", resolved, {"a": infile, "zz": infile}, {"b": tmp_path / "o.txt"})
    assert _exit_code(argv) == ExitCode.INTERNAL
    err = capsys.readouterr().err
    assert "toy" in err and "zz" in err


def test_missing_declared_output_exits_5(
    registry: ModuleType,
    resolved: Path,
    infile: Path,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    def toy(ctx: Any) -> None:
        return None  # forgets to write "b"

    registry.step("toy", inputs=["a"], outputs=["b"])(toy)
    out = tmp_path / "forgotten.txt"
    assert _exit_code(_argv("toy", resolved, {"a": infile}, {"b": out})) == ExitCode.INTERNAL
    err = capsys.readouterr().err
    assert "toy" in err and "b" in err
    assert not out.exists()


def test_missing_required_input_exits_5(
    registry: ModuleType, resolved: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    registry.step("toy", inputs=["a"], outputs=["b"])(_copy)
    missing = tmp_path / "missing.txt"
    argv = _argv("toy", resolved, {"a": missing}, {"b": tmp_path / "o.txt"})
    assert _exit_code(argv) == ExitCode.INTERNAL
    assert "toy" in capsys.readouterr().err


def test_step_outputs_are_byte_identical_across_runs(
    registry: ModuleType, resolved: Path, infile: Path, tmp_path: Path
) -> None:
    def toy(ctx: Any) -> None:
        draws = ctx.rng.normal(size=10)
        ctx.outputs["b"].write_text(ctx.inputs["a"].read_text() + repr(draws.tolist()))

    registry.step("toy", inputs=["a"], outputs=["b"])(toy)
    o1, o2 = tmp_path / "r1" / "b.txt", tmp_path / "r2" / "b.txt"
    assert _exit_code(_argv("toy", resolved, {"a": infile}, {"b": o1})) == 0
    assert _exit_code(_argv("toy", resolved, {"a": infile}, {"b": o2})) == 0
    assert o1.read_bytes() == o2.read_bytes()
