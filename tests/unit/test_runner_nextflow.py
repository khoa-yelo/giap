"""Nextflow launcher (M2-06 examples 1-4; GIPPA_SPEC.md §6.3 "Launch", "Error strategy").

Pure Python: no Nextflow is started. ``run_nextflow`` is exercised with a stand-in command (the
current Python interpreter) so only its exit-code mapping and trace parsing are tested.

Conventions assumed here (M2-06 public interface):
- ``write_params`` writes ``OUT/work/params.json`` and a file named ``resolved_config.yaml`` under
  ``OUT``; ``resolve_config(resolved_config.yaml)`` gives back the same configuration (the file is
  what ``GIPPA_STEP`` passes to ``gippa _step --config``).
- ``params.json`` carries ``threads``, ``memory``, ``outdir``, ``modules`` and ``annotation``
  (the ``nextflow.config`` params), and every input path as given (absolute).
- ``failed_step`` returns the last FAILED row whose task has no later COMPLETED row (review
  decision 2) and names the step by the process name's last ``:`` component, lower-cased
  (``ANNOTATE:TERMINATORNET (chr1)`` -> ``terminatornet``); for the generic ``GIPPA_STEP`` process
  the tag (the Python step name) is the step. The work directory is the absolute task directory
  ``OUT/work/<hash dir>``.
"""

from __future__ import annotations

import importlib
import json
import sys
import time
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.config import GippaConfig, resolve_config
from gippa.exitcodes import ExitCode


def _nf() -> ModuleType:
    return importlib.import_module("gippa.runner.nextflow")


W = Path("/opt/gippa/workflow/main.nf")


# ---------------------------------------------------------------- launch_command (examples 1, 2)


def _expected(out: Path) -> list[str]:
    return [
        "nextflow",
        "run",
        str(W),
        "-params-file",
        str(out / "work/params.json"),
        "-profile",
        "farmshare,apptainer",
        "-work-dir",
        str(out / "work"),
        "-with-trace",
        str(out / "logs/trace.tsv"),
        "-with-timeline",
        str(out / "logs/timeline.html"),
        "-with-report",
        str(out / "logs/nf_report.html"),
        "-with-dag",
        str(out / "logs/dag.html"),
    ]


def test_launch_command_example_1(tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    cmd = _nf().launch_command(
        W, out / "work/params.json", "farmshare,apptainer", out, resume=False
    )
    assert cmd == _expected(out)
    assert all(isinstance(c, str) for c in cmd)


def test_launch_command_example_2_resume(tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    cmd = _nf().launch_command(W, out / "work/params.json", "farmshare,apptainer", out, resume=True)
    assert cmd == [*_expected(out), "-resume"]


def test_launch_command_never_uses_plugins_or_user_params(tmp_path: Path) -> None:
    """§6.3: Nextflow never parses user parameters itself; no plugins in version 1."""
    out = tmp_path / "OUT"
    cmd = _nf().launch_command(W, out / "work/params.json", "test", out, resume=True)
    assert not any(c.startswith("--") for c in cmd)
    assert "-plugins" not in cmd


_profile = st.lists(
    st.sampled_from(
        ["farmshare", "farmshare_local", "apptainer", "conda", "docker", "test", "standard"]
    ),
    min_size=1,
    max_size=3,
    unique=True,
).map(",".join)


@settings(max_examples=50, deadline=None)
@given(profile=_profile, resume=st.booleans())
def test_launch_command_shape_property(profile: str, resume: bool) -> None:
    out = Path("/scratch/OUT")
    params = out / "work/params.json"
    cmd = _nf().launch_command(W, params, profile, out, resume)
    assert cmd[:3] == ["nextflow", "run", str(W)]
    assert cmd[cmd.index("-profile") + 1] == profile
    assert cmd[cmd.index("-params-file") + 1] == str(params)
    assert cmd[cmd.index("-work-dir") + 1] == str(out / "work")
    assert ("-resume" in cmd) is resume
    if resume:
        assert cmd[-1] == "-resume"
    base = _nf().launch_command(W, params, profile, out, False)
    assert cmd[: len(base)] == base


# ---------------------------------------------------------------- write_params (example 3)


def _inputs(tmp_path: Path) -> dict[str, Path]:
    d = tmp_path / "in"
    d.mkdir(exist_ok=True)
    genome = d / "genome.fa"
    genome.write_text(">c1\nACGT\n")
    peptides = d / "peptides.tsv"
    peptides.write_text("peptide\nPEPTIDE\n")
    return {"genome": genome, "peptides": peptides}


def _strings(obj: Any) -> list[str]:
    if isinstance(obj, str):
        return [obj]
    if isinstance(obj, dict):
        return [s for v in obj.values() for s in _strings(v)]
    if isinstance(obj, list):
        return [s for v in obj for s in _strings(v)]
    return []


def _assert_sorted_keys(text: str) -> None:
    def hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        keys = [k for k, _ in pairs]
        assert keys == sorted(keys), f"params.json keys not sorted: {keys}"
        return dict(pairs)

    json.loads(text, object_pairs_hook=hook)


def test_write_params_locations_and_content(tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    cfg = resolve_config(overrides=["threads=3", "modules=[core, mobile]"])
    inputs = _inputs(tmp_path)
    params, resolved = _nf().write_params(out, cfg, inputs)
    assert params == out / "work" / "params.json"
    assert resolved.name == "resolved_config.yaml"
    assert out in resolved.parents
    data = json.loads(params.read_text())
    for key in ("threads", "memory", "outdir", "modules", "annotation"):
        assert key in data, f"params.json lacks params.{key}"
    assert data["threads"] == 3
    assert data["modules"] == ["core", "mobile"]
    assert Path(data["outdir"]) == out
    values = _strings(data)
    for path in inputs.values():
        assert str(path) in values, f"input {path} missing from params.json"


def test_write_params_resolved_config_round_trips(tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    cfg = resolve_config(overrides=["seed=7", "annotation.terminators.bactermfinder=true"])
    _, resolved = _nf().write_params(out, cfg, _inputs(tmp_path))
    assert resolve_config(resolved) == cfg


def test_write_params_example_3_byte_identical(tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    cfg = resolve_config()
    inputs = _inputs(tmp_path)
    p1, r1 = _nf().write_params(out, cfg, inputs)
    first = (p1.read_bytes(), r1.read_bytes())
    time.sleep(1.05)  # a timestamp with one-second resolution would change
    p2, r2 = _nf().write_params(out, cfg, inputs)
    assert (p2, r2) == (p1, r1)
    assert (p2.read_bytes(), r2.read_bytes()) == first
    _assert_sorted_keys(p1.read_text())


@settings(max_examples=15, deadline=None)
@given(
    threads=st.integers(min_value=1, max_value=128),
    modules=st.lists(st.sampled_from(["mobile", "bgc", "function"]), unique=True),
    seed=st.integers(min_value=0, max_value=2**31 - 1),
)
def test_write_params_deterministic_property(
    tmp_path_factory: pytest.TempPathFactory, threads: int, modules: list[str], seed: int
) -> None:
    base = tmp_path_factory.mktemp("wp")
    cfg = GippaConfig.model_validate(
        {**resolve_config().model_dump(), "threads": threads, "modules": modules, "seed": seed}
    )
    inputs = _inputs(base)
    a = _nf().write_params(base / "A", cfg, inputs)[0].read_text()
    b = _nf().write_params(base / "A", cfg, inputs)[0].read_text()
    assert a == b
    _assert_sorted_keys(a)
    data = json.loads(a)
    assert data["threads"] == threads
    assert data["modules"] == cfg.modules


# ---------------------------------------------------------------- failed_step, run_nextflow (ex. 4)

_FIELDS = [
    "task_id",
    "hash",
    "native_id",
    "name",
    "status",
    "exit",
    "submit",
    "duration",
    "realtime",
    "%cpu",
    "peak_rss",
    "peak_vmem",
    "rchar",
    "wchar",
    "workdir",
]


def _row(out: Path, task: int, hash_: str, name: str, status: str, exit_: str) -> str:
    a, b = hash_.split("/")
    workdir = out / "work" / a / (b + "0123456789abcdef0123456789")
    workdir.mkdir(parents=True, exist_ok=True)
    vals = [str(task), hash_, str(1000 + task), name, status, exit_, "2026-10-14 10:00:00.000"]
    vals += ["1s", "1s", "100.0%", "1 MB", "1 MB", "1 KB", "1 KB", str(workdir)]
    return "\t".join(vals)


def _trace(out: Path, rows: list[str]) -> Path:
    path = out / "logs" / "trace.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\t".join(_FIELDS) + "\n" + "\n".join(rows) + "\n")
    return path


def _workdir(out: Path, hash_: str) -> str:
    a, b = hash_.split("/")
    return str(out / "work" / a / (b + "0123456789abcdef0123456789"))


def test_failed_step_example_4(tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    trace = _trace(
        out,
        [
            _row(out, 1, "3a/1b2c3d", "ANNOTATE:PARSE_BAKTA", "COMPLETED", "0"),
            _row(out, 2, "7f/9e8d7c", "ANNOTATE:TERMINATORNET (chr1)", "FAILED", "1"),
        ],
    )
    assert _nf().failed_step(trace) == ("terminatornet", _workdir(out, "7f/9e8d7c"))


def test_failed_step_none_when_all_completed(tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    trace = _trace(
        out,
        [
            _row(out, 1, "3a/1b2c3d", "ANNOTATE:PARSE_BAKTA", "COMPLETED", "0"),
            _row(out, 2, "4b/2c3d4e", "ANNOTATE:OSTIR (chr1)", "CACHED", "0"),
        ],
    )
    assert _nf().failed_step(trace) is None


def test_failed_step_last_unrecovered_failed_row(tmp_path: Path) -> None:
    """Review decision 2 (was: the first FAILED/ABORTED row, i.e. genomad)."""
    out = tmp_path / "OUT"
    trace = _trace(
        out,
        [
            _row(out, 1, "3a/1b2c3d", "ANNOTATE:PARSE_BAKTA", "COMPLETED", "0"),
            _row(out, 2, "5c/aaaaaa", "ANNOTATE:GENOMAD", "ABORTED", "-"),
            _row(out, 3, "6d/bbbbbb", "ANNOTATE:TERMINATORNET (chr1)", "FAILED", "1"),
        ],
    )
    assert _nf().failed_step(trace) == ("terminatornet", _workdir(out, "6d/bbbbbb"))


def test_failed_step_gippa_step_uses_tag(tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    trace = _trace(
        out, [_row(out, 1, "8e/cccccc", "ANNOTATE:GIPPA_STEP (parse_bakta)", "FAILED", "4")]
    )
    assert _nf().failed_step(trace) == ("parse_bakta", _workdir(out, "8e/cccccc"))


def _py(code: str) -> list[str]:
    return [sys.executable, "-c", code]


def test_run_nextflow_success_returns_0(tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    out.mkdir()
    assert _nf().run_nextflow(_py("import sys; sys.exit(0)"), out) == ExitCode.SUCCESS


def test_run_nextflow_example_4_returns_4_and_names_step(
    tmp_path: Path, capfd: pytest.CaptureFixture[str]
) -> None:
    out = tmp_path / "OUT"
    _trace(
        out,
        [
            _row(out, 1, "3a/1b2c3d", "ANNOTATE:PARSE_BAKTA", "COMPLETED", "0"),
            _row(out, 2, "7f/9e8d7c", "ANNOTATE:TERMINATORNET (chr1)", "FAILED", "1"),
        ],
    )
    code = _nf().run_nextflow(_py("import sys; sys.exit(1)"), out)
    assert code == ExitCode.STEP_FAILED
    assert int(code) == 4
    cap = capfd.readouterr()
    msg = f"step terminatornet failed; log: {_workdir(out, '7f/9e8d7c')}/.command.log"
    assert msg in cap.err  # review decision 6: messages on stderr
