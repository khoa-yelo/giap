"""Genome-wide Promotech: valid config, rejected at run time (M2-13 Q2, revised lead decision 4).

``annotation.promoters.mode: genome`` (GIPPA_SPEC.md §5.2; ``gippa run --promoters``) loads and
passes ``gippa preflight`` (resource checks: tests/unit/preflight/test_run_review.py). ``gippa run``
exits 1 before launching anything with "annotation.promoters.mode: genome is not implemented", and
the ``promoter_windows`` step raises the same error if it is ever reached.
"""

from __future__ import annotations

import contextlib
import importlib
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from click.testing import CliRunner, Result

GENOME_SET = ["--set", "annotation.promoters.mode=genome"]


def _cfg(tmp_path: Path) -> Path:
    p = tmp_path / "gippa.yaml"
    p.write_text("annotation:\n  promoters:\n    mode: genome\n")
    return p


def _fasta(tmp_path: Path) -> Path:
    p = tmp_path / "g.fa"
    p.write_text(">c1\n" + "ACGT" * 50 + "\n")
    return p


def _text(r: Result) -> str:
    out = r.output
    with contextlib.suppress(ValueError, AttributeError):  # stderr may be mixed into output
        out += r.stderr
    return out


def _says_genome_not_implemented(text: str) -> bool:
    """One line naming the mode (not the generic ``gippa run: not implemented`` stub)."""
    return any("genome" in ln and "not implemented" in ln for ln in text.splitlines())


# ------------------------------------------------------------------- config load and preflight


def test_model_accepts_genome_mode() -> None:
    from gippa.config import GippaConfig

    cfg = GippaConfig.model_validate({"annotation": {"promoters": {"mode": "genome"}}})
    assert cfg.annotation.promoters.mode == "genome"


def test_resolve_config_accepts_genome_mode_from_file_and_set(tmp_path: Path) -> None:
    from gippa.config import resolve_config

    assert resolve_config(_cfg(tmp_path)).annotation.promoters.mode == "genome"
    assert resolve_config(None, GENOME_SET[1:]).annotation.promoters.mode == "genome"


def test_preflight_treats_genome_mode_like_targeted(tmp_path: Path) -> None:
    from gippa.cli import cli

    g = _fasta(tmp_path)
    base = CliRunner().invoke(cli, ["preflight", "--genome", str(g)])
    gm = CliRunner().invoke(cli, ["preflight", "--genome", str(g), *GENOME_SET])
    assert gm.exit_code == base.exit_code, (_text(base), _text(gm))
    assert not _says_genome_not_implemented(_text(gm))
    assert "invalid configuration" not in _text(gm)


# ------------------------------------------------------------------------------ gippa run


@pytest.fixture
def launches(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Record (and block) every way ``gippa run`` could start work."""
    calls: list[str] = []

    def _stub(name: str) -> Any:
        def f(*a: Any, **k: Any) -> Any:
            calls.append(name)
            raise AssertionError(f"{name} launched in genome mode")

        return f

    from gippa.runner import local, nextflow
    from gippa.steps import registry

    monkeypatch.setattr(local, "plan_local", _stub("plan_local"))
    monkeypatch.setattr(local, "run_local", _stub("run_local"))
    monkeypatch.setattr(nextflow, "run_nextflow", _stub("run_nextflow"))
    monkeypatch.setattr(registry, "run_step", _stub("run_step"))
    import subprocess

    monkeypatch.setattr(subprocess, "run", _stub("subprocess.run"))
    monkeypatch.setattr(subprocess, "Popen", _stub("subprocess.Popen"))
    return calls


@pytest.mark.parametrize("executor", ["local", "nextflow"])
@pytest.mark.parametrize("how", ["flag", "set", "config"])
def test_run_genome_mode_exits_1_before_launching(
    tmp_path: Path, launches: list[str], executor: str, how: str
) -> None:
    from gippa.cli import cli

    mode = {
        "flag": ["--promoters", "genome"],
        "set": GENOME_SET,
        "config": ["--config", str(_cfg(tmp_path))],
    }[how]
    args = ["run", "--genome", str(_fasta(tmp_path)), "--outdir", str(tmp_path / "out")]
    r = CliRunner().invoke(cli, [*args, "--executor", executor, *mode])
    assert r.exit_code == 1, (r.exit_code, _text(r), r.exception)
    assert _says_genome_not_implemented(_text(r)), _text(r)
    assert "annotation.promoters.mode" in _text(r)
    assert launches == []


# ---------------------------------------------------------------------- promoter_windows step


def test_promoter_windows_step_rejects_genome_mode(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippa.config import GippaConfig
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    names = load_step_modules()
    monkeypatch.setattr(reg, "_REGISTRY", {})
    for n in names:
        if n != "registry":
            importlib.reload(importlib.import_module(f"gippa.steps.{n}"))
    spec = reg.get("promoter_windows")
    assert spec.func is not None
    inputs: dict[str, Any] = {}
    for k in (*spec.inputs, *spec.optional_inputs):
        p = tmp_path / "in" / k
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"")
        inputs[k] = p
    outputs = {k: tmp_path / "out" / k for k in spec.outputs}
    cfg = GippaConfig.model_validate({"annotation": {"promoters": {"mode": "genome"}}})
    ctx = reg.StepContext(cfg, inputs, outputs, np.random.default_rng(0), logging.getLogger("test"))
    with pytest.raises(ValueError, match=r"annotation\.promoters\.mode: genome is not implemented"):
        spec.func(ctx)
