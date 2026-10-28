"""``gippa _step NAME`` finds every registered step (M2-07 follow-up, code review item 5 of PR #65).

Nextflow calls ``gippa _step <name>`` in a fresh process, where nothing has imported the step
modules yet: the command must call ``gippa.steps.load_step_modules()`` before looking the name up.
An unknown name stays a wiring error, exit code 5 (M2-02 example 7; GIPPA_SPEC.md §5.1).

Each test runs in a fresh interpreter (as under Nextflow), so the in-process registry state of
other tests cannot hide the bug. The steps get no ``--in``/``--out``: a found step fails its
key check ("required input … not given", also exit 5), an unfound one says "unknown step".

The every-step check runs all ``gippa _step`` calls in ONE fresh interpreter (backlog: one process
per step took ~55 s and grew with every step). Its first call happens before anything has
imported a step module, as under Nextflow; the child then lists the registry itself and calls
``main(["_step", name, ...])`` for each name, the same entry point the console script uses.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from gippa.config import GippaConfig
from gippa.exitcodes import ExitCode

REPO = Path(__file__).resolve().parents[2]
_MAIN = "import sys; from gippa.cli import main; main(sys.argv[1:])"
# One subprocess per call can be slow on a loaded batch node (NFS imports; Slurm 1849060).
_TIMEOUT = 300

# Child for the every-step check. argv: FIRST CONFIG. Prints JSON
# {"fresh": [...], "registered": [[name, kind], ...], "calls": {name: [exit_code, stderr]}}.
_ALL_STEPS = """
import contextlib, io, json, sys
from gippa.cli import main

first, config = sys.argv[1], sys.argv[2]
fresh = sorted(m for m in sys.modules if m.startswith("gippa.steps."))


def call(name):
    buf = io.StringIO()
    code = None
    with contextlib.redirect_stderr(buf), contextlib.redirect_stdout(buf):
        try:
            main(["_step", name, "--config", config])
        except SystemExit as exc:
            code = exc.code
    return [code, buf.getvalue()]


calls = {first: call(first)}  # before anything else imports a step module
from gippa.steps import load_step_modules, registry as r

load_step_modules()
registered = [[s.name, s.kind] for s in r.all_steps()]
for name, _kind in registered:
    if name not in calls:
        calls[name] = call(name)
print(json.dumps({"fresh": fresh, "registered": registered, "calls": calls}))
"""


def _env() -> dict[str, str]:
    env = dict(os.environ)
    paths = [str(REPO / "src"), str(REPO / "tests" / "sim")]
    env["PYTHONPATH"] = os.pathsep.join([*paths, env.get("PYTHONPATH", "")])
    return env


def _gippa(argv: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-c", _MAIN, *argv],
        cwd=cwd, env=_env(), capture_output=True, text=True, timeout=_TIMEOUT, check=False,
    )  # fmt: skip


@pytest.fixture
def resolved(tmp_path: Path) -> Path:
    path = tmp_path / "resolved.yaml"
    path.write_text(yaml.safe_dump(GippaConfig().model_dump(mode="json"), sort_keys=True))
    return path


@pytest.mark.parametrize("name", ["parse_bakta", "uniop_input", "operon_pairs"])
def test_step_finds_a_registered_step_in_a_fresh_process(
    name: str, resolved: Path, tmp_path: Path
) -> None:
    res = _gippa(["_step", name, "--config", str(resolved)], tmp_path)
    err = res.stdout + res.stderr
    assert "unknown step" not in err, err
    assert res.returncode == ExitCode.INTERNAL, err
    assert "required input" in err, err


def test_step_finds_every_registered_step(resolved: Path, tmp_path: Path) -> None:
    res = subprocess.run(
        [sys.executable, "-c", _ALL_STEPS, "parse_bakta", str(resolved)],
        cwd=tmp_path, env=_env(), capture_output=True, text=True, timeout=_TIMEOUT, check=False,
    )  # fmt: skip
    assert res.returncode == 0, res.stdout + res.stderr
    out = json.loads(res.stdout.strip().splitlines()[-1])
    assert out["fresh"] == [], f"step modules imported before the first call: {out['fresh']}"
    names = [n for n, _ in out["registered"]]
    assert {"parse_bakta", "uniop_input"} <= set(names), names
    assert set(out["calls"]) == set(names), sorted(set(out["calls"]) ^ set(names))
    first_code, first_err = out["calls"]["parse_bakta"]
    assert first_code == ExitCode.INTERNAL and "required input" in first_err, first_err
    for name in names:
        code, err = out["calls"][name]
        assert "unknown step" not in err, f"{name}: {err}"
        assert code is not None, f"{name}: main() returned without an exit code: {err}"


def test_step_unknown_name_is_still_exit_5(resolved: Path, tmp_path: Path) -> None:
    res = _gippa(["_step", "no_such_step", "--config", str(resolved)], tmp_path)
    err = res.stdout + res.stderr
    assert res.returncode == ExitCode.INTERNAL, err
    assert "unknown step 'no_such_step'" in err, err
