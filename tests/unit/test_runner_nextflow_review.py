"""Runner review fixes (M2-06 "Lead decisions after code review", items 1, 2, 5, 6).

Pure Python: no Nextflow is started (the Nextflow-side checks live in
``tests/integration/test_nf_review.py``).

- (1) A trace written with Nextflow's DEFAULT field set (no ``workdir`` column; header taken from a
  local Nextflow 26.10 run) still names the failed step; the work directory is rebuilt from the
  abbreviated ``hash`` under the run's work dir ``OUT/work`` (the trace lives at
  ``OUT/logs/trace.tsv``) and must be the existing task directory.
- (2) ``failed_step`` = the last FAILED row whose task name has no later COMPLETED row; a
  retried-then-succeeded attempt never counts. For ``GIPPA_STEP`` the tag is the step, so two
  ``GIPPA_STEP`` tasks with different tags are different steps.
- (6) Input paths in ``params.json`` are absolute even when given relative (Nextflow's cwd is OUT);
  a failed run without a trace still exits 4; runner messages go to stderr only.
"""

from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path
from types import ModuleType

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.config import resolve_config
from gippa.exitcodes import ExitCode


def _nf() -> ModuleType:
    return importlib.import_module("gippa.runner.nextflow")


# Nextflow's default trace fields (no workdir), verified with `nextflow run ... -with-trace`.
DEFAULT_FIELDS = [
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
]
_FULL_SUFFIX = "0123456789abcdef0123456789"


def _taskdir(out: Path, hash_: str) -> Path:
    a, b = hash_.split("/")
    return out / "work" / a / (b + _FULL_SUFFIX)


def _default_row(out: Path, task: int, hash_: str, name: str, status: str, exit_: str) -> str:
    _taskdir(out, hash_).mkdir(parents=True, exist_ok=True)
    vals = [str(task), hash_, str(1000 + task), name, status, exit_, "2026-10-14 10:00:00.000"]
    vals += ["1s", "1s", "100.0%", "1 MB", "1 MB", "1 KB", "1 KB"]
    return "\t".join(vals)


def _default_trace(out: Path, rows: list[str]) -> Path:
    path = out / "logs" / "trace.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\t".join(DEFAULT_FIELDS) + "\n" + "\n".join(rows) + "\n")
    return path


def _full_trace(out: Path, rows: list[tuple[int, str, str, str, str]]) -> Path:
    """Trace with the configured field set of decision 1 (incl. tag and workdir)."""
    fields = ["task_id", "hash", "name", "tag", "status", "exit", "attempt", "workdir"]
    lines = ["\t".join(fields)]
    for task, hash_, name, status, exit_ in rows:
        d = _taskdir(out, hash_)
        d.mkdir(parents=True, exist_ok=True)
        tag = name.split(" (", 1)[1].rstrip(")") if " (" in name else "-"
        lines.append("\t".join([str(task), hash_, name, tag, status, exit_, "1", str(d)]))
    path = out / "logs" / "trace.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


# ---------------------------------------------------------------- decision 1: default trace fields


def test_failed_step_with_default_trace_fields(tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    trace = _default_trace(
        out,
        [
            _default_row(out, 1, "3a/1b2c3d", "ANNOTATE:PARSE_BAKTA", "COMPLETED", "0"),
            _default_row(out, 2, "7f/9e8d7c", "ANNOTATE:TERMINATORNET (chr1)", "FAILED", "1"),
        ],
    )
    hit = _nf().failed_step(trace)
    assert hit is not None, "no failed step found in a default-field trace"
    step, workdir = hit
    assert step == "terminatornet"
    assert Path(workdir) == _taskdir(out, "7f/9e8d7c"), workdir


def test_run_nextflow_default_trace_names_step_and_log(
    tmp_path: Path, capfd: pytest.CaptureFixture[str]
) -> None:
    out = tmp_path / "OUT"
    _default_trace(
        out, [_default_row(out, 1, "8e/cccccc", "ANNOTATE:GIPPA_STEP (parse_bakta)", "FAILED", "4")]
    )
    code = _nf().run_nextflow([sys.executable, "-c", "import sys; sys.exit(1)"], out)
    assert code == ExitCode.STEP_FAILED
    err = capfd.readouterr().err
    assert f"step parse_bakta failed; log: {_taskdir(out, '8e/cccccc')}/.command.log" in err


# ---------------------------------------------------------------- decision 2: retries


def test_retry_then_success_of_a_then_failure_of_b(tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    trace = _full_trace(
        out,
        [
            (1, "11/aaaaaa", "ANNOTATE:BAKTA (chr1)", "FAILED", "137"),
            (2, "22/bbbbbb", "ANNOTATE:BAKTA (chr1)", "COMPLETED", "0"),
            (3, "33/cccccc", "ANNOTATE:GENOMAD (chr1)", "FAILED", "1"),
        ],
    )
    assert _nf().failed_step(trace) == ("genomad", str(_taskdir(out, "33/cccccc")))


def test_retry_then_success_of_a_then_failure_of_b_default_fields(tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    trace = _default_trace(
        out,
        [
            _default_row(out, 1, "11/aaaaaa", "ANNOTATE:BAKTA (chr1)", "FAILED", "137"),
            _default_row(out, 3, "33/cccccc", "ANNOTATE:GENOMAD (chr1)", "FAILED", "1"),
            _default_row(out, 2, "22/bbbbbb", "ANNOTATE:BAKTA (chr1)", "COMPLETED", "0"),
        ],
    )
    hit = _nf().failed_step(trace)
    assert hit is not None
    assert hit[0] == "genomad"
    assert Path(hit[1]) == _taskdir(out, "33/cccccc")


def test_retried_then_succeeded_only_is_none(tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    trace = _full_trace(
        out,
        [
            (1, "11/aaaaaa", "ANNOTATE:BAKTA (chr1)", "FAILED", "143"),
            (2, "22/bbbbbb", "ANNOTATE:BAKTA (chr1)", "COMPLETED", "0"),
        ],
    )
    assert _nf().failed_step(trace) is None


def test_gippa_step_tags_are_distinct_steps(tmp_path: Path) -> None:
    """A later COMPLETED GIPPA_STEP with another tag does not hide the failed one."""
    out = tmp_path / "OUT"
    trace = _full_trace(
        out,
        [
            (1, "11/aaaaaa", "ANNOTATE:GIPPA_STEP (parse_bakta)", "FAILED", "1"),
            (2, "22/bbbbbb", "ANNOTATE:GIPPA_STEP (operons)", "COMPLETED", "0"),
        ],
    )
    assert _nf().failed_step(trace) == ("parse_bakta", str(_taskdir(out, "11/aaaaaa")))


_NAMES = [
    "ANNOTATE:BAKTA (chr1)",
    "ANNOTATE:GENOMAD (chr1)",
    "ANNOTATE:TERMINATORNET (chr1)",
    "ANNOTATE:GIPPA_STEP (parse_bakta)",
    "ANNOTATE:GIPPA_STEP (operons)",
]
_STEP = {
    "ANNOTATE:BAKTA (chr1)": "bakta",
    "ANNOTATE:GENOMAD (chr1)": "genomad",
    "ANNOTATE:TERMINATORNET (chr1)": "terminatornet",
    "ANNOTATE:GIPPA_STEP (parse_bakta)": "parse_bakta",
    "ANNOTATE:GIPPA_STEP (operons)": "operons",
}


def _reference(events: list[tuple[str, str]]) -> str | None:
    for i in range(len(events) - 1, -1, -1):
        name, status = events[i]
        if status == "FAILED" and not any(
            n == name and s == "COMPLETED" for n, s in events[i + 1 :]
        ):
            return _STEP[name]
    return None


@settings(max_examples=60, deadline=None)
@given(
    events=st.lists(
        st.tuples(st.sampled_from(_NAMES), st.sampled_from(["COMPLETED", "FAILED", "CACHED"])),
        max_size=8,
    )
)
def test_failed_step_last_unrecovered_failure_property(
    tmp_path_factory: pytest.TempPathFactory, events: list[tuple[str, str]]
) -> None:
    out = tmp_path_factory.mktemp("tr") / "OUT"
    rows = [
        (i + 1, f"{i:02x}/{i:06x}", name, status, "1" if status == "FAILED" else "0")
        for i, (name, status) in enumerate(events)
    ]
    hit = _nf().failed_step(_full_trace(out, rows))
    expected = _reference(events)
    assert (hit[0] if hit else None) == expected


# ---------------------------------------------------------------- decision 6


def test_write_params_relative_inputs_become_absolute(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "genome.fa").write_text(">c1\nACGT\n")
    (tmp_path / "in" / "peptides.tsv").write_text("peptide\nPEPTIDE\n")
    inputs = {"genome": Path("in/genome.fa"), "peptides": Path("in/peptides.tsv")}
    params, _ = _nf().write_params(Path("OUT"), resolve_config(), inputs)
    data = json.loads(params.read_text())
    assert set(data["inputs"]) == {"genome", "peptides"}
    for key, rel in inputs.items():
        got = Path(data["inputs"][key])
        assert got.is_absolute(), f"input {key} not absolute: {got}"
        assert got == tmp_path / rel
    assert Path(data["outdir"]).is_absolute()
    assert Path(data["outdir"]) == tmp_path / "OUT"
    assert Path(data["config"]).is_absolute()


def test_run_nextflow_without_trace_exits_4_on_stderr(
    tmp_path: Path, capfd: pytest.CaptureFixture[str]
) -> None:
    out = tmp_path / "OUT"
    out.mkdir()
    code = _nf().run_nextflow([sys.executable, "-c", "import sys; sys.exit(3)"], out)
    assert code == ExitCode.STEP_FAILED
    cap = capfd.readouterr()
    assert "nextflow exited 3" in cap.err
    assert str(out / ".nextflow.log") in cap.err
    assert "gippa run" not in cap.out


def test_run_nextflow_step_message_on_stderr_only(
    tmp_path: Path, capfd: pytest.CaptureFixture[str]
) -> None:
    out = tmp_path / "OUT"
    _full_trace(out, [(1, "7f/9e8d7c", "ANNOTATE:TERMINATORNET (chr1)", "FAILED", "1")])
    code = _nf().run_nextflow([sys.executable, "-c", "import sys; sys.exit(1)"], out)
    assert code == ExitCode.STEP_FAILED
    cap = capfd.readouterr()
    msg = f"step terminatornet failed; log: {_taskdir(out, '7f/9e8d7c')}/.command.log"
    assert msg in cap.err
    assert msg not in cap.out


def test_runner_docstring_documents_minimum_of_user_and_profile_limits() -> None:
    """Decision 5: the --threads/--memory vs profile-limit rule is documented in the runner."""
    doc = (_nf().__doc__ or "").lower()
    assert "minimum" in doc or "min(" in doc, "runner docstring does not document the minimum rule"


def test_write_params_relative_precomputed_layers_become_absolute(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A relative ``annotation.precomputed.<key>`` must reach both params.json and the resolved
    config absolute: region_layer reads it inside a Nextflow task directory (review of #146)."""
    import yaml

    monkeypatch.chdir(tmp_path)
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "regions.tsv").write_text("region_id\n")
    cfg = resolve_config(overrides=["annotation.precomputed.regions=in/regions.tsv"])
    params, resolved = _nf().write_params(Path("OUT"), cfg, {})
    want = str(tmp_path / "in" / "regions.tsv")
    assert json.loads(params.read_text())["annotation"]["precomputed"]["regions"] == want
    assert yaml.safe_load(resolved.read_text())["annotation"]["precomputed"]["regions"] == want
    assert json.loads(params.read_text())["annotation"]["precomputed"]["operon_pairs"] is None
