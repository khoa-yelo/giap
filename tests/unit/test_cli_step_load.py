"""``giap _step NAME`` finds every registered step (M2-07 follow-up, code review item 5 of PR #65).

Nextflow calls ``giap _step <name>`` in a fresh process, where nothing has imported the step
modules yet: the command must call ``giap.steps.load_step_modules()`` before looking the name up.
An unknown name stays a wiring error, exit code 5 (M2-02 example 7; GIAP_SPEC.md §5.1).

Each call runs in a fresh interpreter (as under Nextflow), so the in-process registry state of
other tests cannot hide the bug. The steps get no ``--in``/``--out``: a found step fails its
key check ("required input … not given", also exit 5), an unfound one says "unknown step".
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from giap.config import GiapConfig
from giap.exitcodes import ExitCode

REPO = Path(__file__).resolve().parents[2]
_MAIN = "import sys; from giap.cli import main; main(sys.argv[1:])"


def _env() -> dict[str, str]:
    env = dict(os.environ)
    paths = [str(REPO / "src"), str(REPO / "tests" / "sim")]
    env["PYTHONPATH"] = os.pathsep.join([*paths, env.get("PYTHONPATH", "")])
    return env


def _giap(argv: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-c", _MAIN, *argv],
        cwd=cwd, env=_env(), capture_output=True, text=True, timeout=120, check=False,
    )  # fmt: skip


def _registered() -> list[tuple[str, str]]:
    code = (
        "import json; from giap.steps import load_step_modules, registry as r; "
        "load_step_modules(); print(json.dumps([[s.name, s.kind] for s in r.all_steps()]))"
    )
    res = subprocess.run(
        [sys.executable, "-c", code], env=_env(), capture_output=True, text=True, check=True
    )
    return [(n, k) for n, k in json.loads(res.stdout)]


@pytest.fixture
def resolved(tmp_path: Path) -> Path:
    path = tmp_path / "resolved.yaml"
    path.write_text(yaml.safe_dump(GiapConfig().model_dump(mode="json"), sort_keys=True))
    return path


@pytest.mark.parametrize("name", ["parse_bakta", "uniop_input", "operon_pairs"])
def test_step_finds_a_registered_step_in_a_fresh_process(
    name: str, resolved: Path, tmp_path: Path
) -> None:
    res = _giap(["_step", name, "--config", str(resolved)], tmp_path)
    err = res.stdout + res.stderr
    assert "unknown step" not in err, err
    assert res.returncode == ExitCode.INTERNAL, err
    assert "required input" in err, err


def test_step_finds_every_registered_step(resolved: Path, tmp_path: Path) -> None:
    steps = _registered()
    assert {"parse_bakta", "uniop_input"} <= {n for n, _ in steps}, steps
    for name, _kind in steps:
        res = _giap(["_step", name, "--config", str(resolved)], tmp_path)
        err = res.stdout + res.stderr
        assert "unknown step" not in err, f"{name}: {err}"


def test_step_unknown_name_is_still_exit_5(resolved: Path, tmp_path: Path) -> None:
    res = _giap(["_step", "no_such_step", "--config", str(resolved)], tmp_path)
    err = res.stdout + res.stderr
    assert res.returncode == ExitCode.INTERNAL, err
    assert "unknown step 'no_such_step'" in err, err
