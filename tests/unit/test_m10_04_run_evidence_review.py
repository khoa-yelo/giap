"""``validation/pilot/run_evidence.py arm``: the code review's fixes (review of 45efc4c9, lead
decisions B1, B2, S1-S4, N1-N3; docs/specs/M10-04-pilot-runs.md rule set A).

* B1: A6 fails ``cpu_cap:unknown:<n>`` when ``n`` trace tasks have no ``sacct`` allocation row
  (no native IDs at all: ``n`` is the task count). The CLI and the pure function share one ID
  filter, so they agree on any trace.
* B2: ``--expect-commit`` matches ``^[0-9a-f]{7,40}$`` (else exit 2); the only match is
  ``commit.startswith(want)``; a null or missing manifest commit is ``commit_mismatch``. The
  manifest holds the full 40-character sha (``FULL_COMMIT``), ``--expect-commit`` the short one.
* S1: every ``--expect-input role=sha`` is checked: ``input_missing:<role>`` when the manifest has
  no such role, ``input_mismatch:<role>`` when the sha differs.
* S2: an unparseable ``run_manifest.json`` is ``manifest_invalid:/`` (exit 1, evidence written); a
  trace ``attempt`` of ``-`` is not a traceback.
* S3: the evidence copies ``skipped`` and ``layers_missing`` from the manifest (no failure rule).
* S4: tests that kill the reviewer's surviving mutants (``mutants.log``: 3-7, 12, 15,
  17-19, 24, 38).
* N1: no head ``sacct`` row is ``run_failed:no_sacct``; N2: ``--expect-input`` without ``=`` is an
  argparse error (exit 2); N3: a ``sacct`` line with the wrong field count is a parse error.
"""

from __future__ import annotations

import json
import re
import tempfile
from pathlib import Path
from typing import Any

import pytest
from helpers.pilot_runs_m10_04 import (
    COMMIT,
    FULL_COMMIT,
    HEAD_JOB,
    SHA_ANNOTATION,
    SHA_GENOME,
    SHA_PEPTIDES,
    TASK_JOBS,
    arm_cli_args,
    arm_doc,
    default_trace,
    edit_manifest,
    sacct_rows,
    trace_task,
    write_arm,
)
from helpers.pilot_runs_m10_04 import run_evidence_cli as cli
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

EXPECTED_INPUTS = {"genome": SHA_GENOME, "annotation": SHA_ANNOTATION, "peptides": SHA_PEPTIDES}


def _arm_cli(tmp: Path, out: Path, rows: list[dict[str, str]] | None = None, **kw: str):  # type: ignore[no-untyped-def]
    dest = tmp / "evidence" / "a2.json"
    return cli(tmp, *arm_cli_args(out, dest, **kw), rows=rows), dest


def _head_rows(**kw: Any) -> list[dict[str, str]]:
    """Only the head job's rows (no task-job rows)."""
    return [r for r in sacct_rows(**kw) if r["JobID"].split(".")[0] == HEAD_JOB]


def _no_traceback(res: Any) -> None:
    assert "Traceback" not in res.stderr, res.stderr


# ------------------------------------------------------------------ B1: CPU count unknown


def test_a6_fails_when_no_task_has_a_sacct_row(tmp_path: Path) -> None:
    """The reviewer's probe: only the head's rows; ``max_alloc_cpus`` 0 used to pass."""
    out = write_arm(tmp_path / "out_a2")
    assert arm_doc(out, _head_rows())["failures"] == ["cpu_cap:unknown:2"]


def test_a6_counts_the_tasks_without_a_row(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2")
    rows = sacct_rows(tasks={TASK_JOBS[0]: 14})
    assert arm_doc(out, rows)["failures"] == ["cpu_cap:unknown:1"]


def test_a6_fails_when_the_trace_has_no_native_ids(tmp_path: Path) -> None:
    """Local executor or wrong profile: ``native_id`` is ``-`` everywhere; n = the task count."""
    trace = [trace_task("GIAP:CORE:GENE_TABLE (gene_table)", "-", task_id=1),
             trace_task("GIAP:CORE:DETECTION (detection)", "-", task_id=2),
             trace_task("GIAP:CORE:OPERONS (operons)", "", task_id=3)]  # fmt: skip
    out = write_arm(tmp_path / "out_a2", trace=trace)
    assert arm_doc(out)["failures"] == ["cpu_cap:unknown:3"]


def test_a6_unknown_through_the_cli(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2")
    res, dest = _arm_cli(tmp_path, out, rows=_head_rows())
    assert res.returncode == 1, res.stdout + res.stderr
    assert json.loads(dest.read_text())["failures"] == ["cpu_cap:unknown:2"]


def test_cli_and_pure_function_share_the_id_filter(tmp_path: Path) -> None:
    """An array-style ``native_id`` (``1900002_4``) with its own ``sacct`` row: whatever the filter
    decides, the CLI (which queries ``sacct``) and ``arm_evidence`` (given every row) agree."""
    trace = [trace_task("GIAP:CORE:GENE_TABLE (gene_table)", "1900002_4", task_id=1),
             trace_task("GIAP:CORE:DETECTION (detection)", TASK_JOBS[1], task_id=2)]  # fmt: skip
    out = write_arm(tmp_path / "out_a2", trace=trace)
    rows = sacct_rows(tasks={"1900002_4": 16, TASK_JOBS[1]: 14})
    pure = arm_doc(out, rows)["failures"]
    res, dest = _arm_cli(tmp_path, out, rows=rows)
    _no_traceback(res)
    assert json.loads(dest.read_text())["failures"] == pure
    assert pure in (["cpu_cap:16"], ["cpu_cap:unknown:1"]), pure


@settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(st.lists(st.booleans(), min_size=1, max_size=5))
def test_a6_unknown_count_property(has_row: list[bool]) -> None:
    """Property: ``cpu_cap:unknown:<k>`` exactly when k > 0 tasks lack an allocation row."""
    ids = [f"19100{i:02d}" for i in range(len(has_row))]
    trace = [trace_task(f"GIAP:CORE:GENE_TABLE (g{i})", j, task_id=i + 1)
             for i, j in enumerate(ids)]  # fmt: skip
    rows = sacct_rows(tasks={j: 4 for j, ok in zip(ids, has_row, strict=True) if ok})
    k = has_row.count(False)
    with tempfile.TemporaryDirectory() as d:
        out = write_arm(Path(d) / "out_a2", trace=trace)
        failures = arm_doc(out, rows)["failures"]
    assert failures == ([f"cpu_cap:unknown:{k}"] if k else [])


# ------------------------------------------------------------------ B2: commit prefix


def test_full_manifest_sha_matches_its_short_prefix(tmp_path: Path) -> None:
    """Real runs: the manifest holds ``rev-parse HEAD`` (40), ``--expect-commit`` the short sha."""
    out = write_arm(tmp_path / "out_a2")
    assert len(FULL_COMMIT) == 40 and FULL_COMMIT.startswith(COMMIT) and len(COMMIT) == 7
    doc = arm_doc(out, commit=COMMIT)
    assert doc["failures"] == [] and doc["manifest"]["git_commit"] == FULL_COMMIT
    assert arm_doc(out, commit=FULL_COMMIT)["failures"] == []


def test_wrong_short_prefix_is_a_mismatch(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2")
    assert arm_doc(out, commit="0123abd")["failures"] == [f"commit_mismatch:{FULL_COMMIT}"]


def test_truncated_manifest_commit_is_a_mismatch(tmp_path: Path) -> None:
    """Only ``commit.startswith(want)``: a manifest commit shorter than ``want`` never matches."""
    out = write_arm(tmp_path / "out_a2", commit="0123ab")
    assert arm_doc(out, commit=COMMIT)["failures"] == ["commit_mismatch:0123ab"]


def test_null_manifest_commit_is_a_mismatch(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2", commit=None)
    assert json.loads((out / "run_manifest.json").read_text())["git_commit"] is None
    failures = arm_doc(out)["failures"]
    assert len(failures) == 1 and failures[0].startswith("commit_mismatch"), failures


def test_missing_manifest_commit_is_a_mismatch(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2")
    edit_manifest(out, git_commit=None)  # deletes the key (also a schema violation)
    failures = arm_doc(out)["failures"]
    assert any(f.startswith("commit_mismatch") for f in failures), failures


def test_empty_expected_commit_never_matches(tmp_path: Path) -> None:
    """``arm_evidence`` called without a commit (the CLI rejects it): a mismatch, not a pass."""
    out = write_arm(tmp_path / "out_a2")
    assert arm_doc(out, commit="")["failures"] == [f"commit_mismatch:{FULL_COMMIT}"]
    out = write_arm(tmp_path / "out_null", commit=None)
    failures = arm_doc(out, commit="")["failures"]
    assert len(failures) == 1 and failures[0].startswith("commit_mismatch"), failures


@pytest.mark.parametrize(
    "bad", ["0", "0123ab", "0123ABC", "g123abc", " 0123abc", FULL_COMMIT + "0", "", "HEAD"]
)
def test_cli_rejects_a_malformed_expect_commit(tmp_path: Path, bad: str) -> None:
    """``--expect-commit`` must match ``^[0-9a-f]{7,40}$``: argparse error, exit 2, no evidence."""
    out = write_arm(tmp_path / "out_a2")
    res, dest = _arm_cli(tmp_path, out, commit=bad)
    assert res.returncode == 2, (bad, res.stdout + res.stderr)
    _no_traceback(res)
    assert not dest.exists()


def test_cli_full_sha_against_short_expect_passes(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2")
    res, dest = _arm_cli(tmp_path, out, commit=COMMIT)
    assert res.returncode == 0, res.stdout + res.stderr
    assert json.loads(dest.read_text())["manifest"]["git_commit"] == FULL_COMMIT
    res, dest = _arm_cli(tmp_path, out, commit="0123abd")
    assert res.returncode == 1, res.stdout + res.stderr
    assert json.loads(dest.read_text())["failures"] == [f"commit_mismatch:{FULL_COMMIT}"]


_HEX40 = st.text(alphabet="0123456789abcdef", min_size=40, max_size=40)


@settings(max_examples=25, deadline=None)
@given(_HEX40, st.integers(min_value=7, max_value=40), st.data())
def test_commit_prefix_property(sha: str, k: int, data: st.DataObject) -> None:
    """Property: every 7-40 character prefix of the manifest sha passes; changing any one of its
    characters fails ``commit_mismatch:<sha>``."""
    with tempfile.TemporaryDirectory() as d:
        out = write_arm(Path(d) / "out_a2", commit=sha)
        assert arm_doc(out, commit=sha[:k])["failures"] == []
        i = data.draw(st.integers(min_value=0, max_value=k - 1))
        c = data.draw(st.sampled_from([h for h in "0123456789abcdef" if h != sha[i]]))
        wrong = sha[:i] + c + sha[i + 1 : k]
        assert arm_doc(out, commit=wrong)["failures"] == [f"commit_mismatch:{sha}"]


# ------------------------------------------------------------------ S1: every expected input


@pytest.mark.parametrize("role", ["genome", "annotation", "peptides"])
def test_every_expected_input_is_checked(tmp_path: Path, role: str) -> None:
    """S1 (annotation used to be ignored) and mutant 24 (genome dropped from the checked roles)."""
    out = write_arm(tmp_path / "out_a2")
    inputs = {**EXPECTED_INPUTS, role: "f" * 64}
    assert arm_doc(out, inputs=inputs)["failures"] == [f"input_mismatch:{role}"]


@pytest.mark.parametrize("role", ["genome", "annotation"])
def test_an_expected_role_absent_from_the_manifest(tmp_path: Path, role: str) -> None:
    out = write_arm(tmp_path / "out_a2")
    doc = json.loads((out / "run_manifest.json").read_text())
    edit_manifest(out, inputs=[r for r in doc["inputs"] if r["role"] != role])
    failures = arm_doc(out)["failures"]
    assert f"input_missing:{role}" in failures, failures
    assert f"input_mismatch:{role}" not in failures, failures


def test_an_extra_expected_role_is_missing(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2")
    inputs = {**EXPECTED_INPUTS, "reference_peptides": "4" * 64}
    assert arm_doc(out, inputs=inputs)["failures"] == ["input_missing:reference_peptides"]


def test_cli_checks_the_annotation_sha(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2")
    res, dest = _arm_cli(tmp_path, out, annotation="9" * 64)
    assert res.returncode == 1, res.stdout + res.stderr
    assert json.loads(dest.read_text())["failures"] == ["input_mismatch:annotation"]
    res, dest = _arm_cli(tmp_path, out, genome="9" * 64)
    assert res.returncode == 1, res.stdout + res.stderr
    assert json.loads(dest.read_text())["failures"] == ["input_mismatch:genome"]


# ------------------------------------------------------------------ S2: unreadable manifest


@pytest.mark.parametrize("data", [b'{"status": "comp', b"", b"\x00\xff garbage"])
def test_unparseable_manifest_is_manifest_invalid(tmp_path: Path, data: bytes) -> None:
    """Edge case "Head job killed before the tail": a truncated manifest fails, never crashes."""
    out = write_arm(tmp_path / "out_a2")
    (out / "run_manifest.json").write_bytes(data)
    failures = arm_doc(out)["failures"]
    assert "manifest_invalid:/" in failures, failures


def test_cli_unparseable_manifest_exits_1_with_evidence(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2")
    (out / "run_manifest.json").write_text('{"status": "comp', encoding="utf-8")
    res, dest = _arm_cli(tmp_path, out)
    assert res.returncode == 1, res.stdout + res.stderr
    _no_traceback(res)
    assert "manifest_invalid:/" in json.loads(dest.read_text())["failures"]


def test_attempt_dash_is_not_a_traceback(tmp_path: Path) -> None:
    """A trace ``attempt`` of ``-`` (no retry information) counts as not retried."""
    trace = default_trace()
    trace[0] = {**trace[0], "attempt": "-"}
    out = write_arm(tmp_path / "out_a2", trace=trace)
    doc = arm_doc(out)
    assert doc["failures"] == [] and doc["tasks"]["retried"] == 0
    res, _ = _arm_cli(tmp_path, out)
    assert res.returncode == 0, res.stdout + res.stderr
    _no_traceback(res)


# ------------------------------------------------------------------ S3: skipped, layers_missing


def test_skipped_and_layers_missing_are_copied(tmp_path: Path) -> None:
    out = write_arm(tmp_path / "out_a2")
    skipped = [{"name": "novel_orfs", "reason": "not_run"},
               {"name": "alt_sites", "reason": "opt_in_off"}]  # fmt: skip
    edit_manifest(out, skipped=skipped, layers_missing=["terminators"])
    doc = arm_doc(out)
    assert doc["manifest"]["skipped"] == skipped
    assert doc["manifest"]["layers_missing"] == ["terminators"]
    assert doc["failures"] == [], "no step_not_run rule for now (backlog)"


def test_empty_skipped_is_copied_as_empty(tmp_path: Path) -> None:
    doc = arm_doc(write_arm(tmp_path / "out_a2"))
    assert doc["manifest"]["skipped"] == [] and doc["manifest"]["layers_missing"] == []


# ------------------------------------------------------------------ S4: surviving mutants (A)


@pytest.mark.parametrize(("state", "code"), [("COMPLETED", "1:0"), ("CANCELLED", "0:0"),
                                             ("COMPLETED", "0:15")])  # fmt: skip
def test_head_job_state_and_exit_code_both_count(tmp_path: Path, state: str, code: str) -> None:
    """Mutants 3 and 4: A2 needs ``COMPLETED`` *and* ``0:0``."""
    out = write_arm(tmp_path / "out_a2")
    failures = arm_doc(out, sacct_rows(state=state, exit_code=code))["failures"]
    assert len(failures) == 1 and failures[0].startswith("run_failed:"), failures


def test_a_failed_step_with_exit_0_fails_a2(tmp_path: Path) -> None:
    """Mutant 5: manifest ``exit_code`` 0, ``status`` completed, but a step ``failed``."""
    out = write_arm(tmp_path / "out_a2")
    doc = json.loads((out / "run_manifest.json").read_text())
    steps = [{**s, "status": "failed"} if s["name"] == "detection" else s for s in doc["steps"]]
    edit_manifest(out, steps=steps)
    assert arm_doc(out)["failures"] == ["run_failed:detection"]


def test_a_nonzero_manifest_exit_code_fails_a2(tmp_path: Path) -> None:
    """Mutant 6: ``status`` completed but ``exit_code`` 4 (and no failed step)."""
    out = write_arm(tmp_path / "out_a2")
    edit_manifest(out, exit_code=4)
    assert "run_failed:4" in arm_doc(out)["failures"]


def test_git_dirty_null_is_dirty(tmp_path: Path) -> None:
    """Mutant 7: ``git_dirty: null`` (not a git checkout) is not ``false``."""
    out = write_arm(tmp_path / "out_a2", dirty=None)
    assert json.loads((out / "run_manifest.json").read_text())["git_dirty"] is None
    assert arm_doc(out)["failures"] == ["dirty_checkout"]


@pytest.mark.parametrize(("cpus", "want"), [(14, []), (15, ["cpu_cap:15"]), (16, ["cpu_cap:16"])])
def test_cpu_cap_boundary(tmp_path: Path, cpus: int, want: list[str]) -> None:
    """Mutant 12: 15 CPUs breaches the cap of 14."""
    out = write_arm(tmp_path / "out_a2")
    assert arm_doc(out, sacct_rows(tasks={TASK_JOBS[0]: cpus, TASK_JOBS[1]: 2}))["failures"] == want


def test_head_rows_do_not_take_a_task_with_a_longer_id(tmp_path: Path) -> None:
    """Mutant 15: task ``19000011`` starts with head ``1900001`` but is not one of its steps."""
    trace = [*default_trace(), trace_task("GIAP:CORE:OPERONS (operons)", "19000011", task_id=3)]
    out = write_arm(tmp_path / "out_a2", trace=trace)
    rows = sacct_rows()
    rows.append({"JobID": "19000011", "State": "COMPLETED", "ExitCode": "0:0",
                 "Elapsed": "00:01:00", "MaxRSS": "9999999K", "NodeList": "oat-04",
                 "AllocCPUS": "2"})  # fmt: skip
    doc = arm_doc(out, rows)
    assert doc["failures"] == []
    assert doc["job"]["max_rss"] == "2104832K"


def test_head_max_rss_and_task_peak_rss(tmp_path: Path) -> None:
    """Mutants 17 and 18: the head MaxRSS is the largest step's; ``K`` is KiB (1024)."""
    out = write_arm(tmp_path / "out_a2")
    rows = sacct_rows()
    rows.insert(2, {"JobID": f"{HEAD_JOB}.extern", "State": "COMPLETED", "ExitCode": "0:0",
                    "Elapsed": "05:59:12", "MaxRSS": "1024K", "NodeList": "wheat-07",
                    "AllocCPUS": "2"})  # fmt: skip
    doc = arm_doc(out, rows)
    assert doc["job"]["max_rss"] == "2104832K"
    assert doc["tasks"]["peak_max_rss_bytes"] == 1258291 * 1024


@pytest.mark.parametrize(("text", "want"), [
    ("2104832K", 2104832 * 1024), ("1.5G", int(1.5 * 1024**3)), ("12M", 12 * 1024**2),
    ("2T", 2 * 1024**4), ("12345", 12345), ("0", 0), ("", None), (None, None),
])  # fmt: skip
def test_rss_bytes_units(text: str | None, want: int | None) -> None:
    from helpers.pilot_runs_m10_04 import load_run_evidence

    assert load_run_evidence().rss_bytes(text) == want


def test_failed_then_retried_task_is_counted(tmp_path: Path) -> None:
    """Mutant 19: a FAILED attempt 1 then a COMPLETED attempt 2 of the same task."""
    op = "GIAP:CORE:OPERONS (operons)"
    trace = [*default_trace(), trace_task(op, "1900004", status="FAILED", task_id=3),
             trace_task(op, "1900005", attempt=2, task_id=4)]  # fmt: skip
    out = write_arm(tmp_path / "out_a2", trace=trace)
    rows = sacct_rows(tasks={**dict.fromkeys(TASK_JOBS, 14), "1900004": 4, "1900005": 4})
    tasks = arm_doc(out, rows)["tasks"]
    assert tasks["count"] == 4 and tasks["failed"] == 1 and tasks["retried"] == 1


def test_cli_queries_sacct_in_chunks(tmp_path: Path) -> None:
    """Mutant 38: 250 task jobs; the one breaching the cap is past the first chunk of 200."""
    ids = [str(1910000 + i) for i in range(250)]
    trace = [trace_task(f"GIAP:CORE:GENE_TABLE (g{i})", j, task_id=i + 1)
             for i, j in enumerate(ids)]  # fmt: skip
    out = write_arm(tmp_path / "out_a2", trace=trace)
    rows = sacct_rows(tasks={j: (16 if i == 240 else 4) for i, j in enumerate(ids)})
    res, dest = _arm_cli(tmp_path, out, rows=rows)
    assert res.returncode == 1, res.stdout + res.stderr
    doc = json.loads(dest.read_text())
    assert doc["failures"] == ["cpu_cap:16"] and doc["tasks"]["count"] == 250


# ------------------------------------------------------------------ nits N1-N3


def test_no_head_sacct_rows_is_run_failed_no_sacct(tmp_path: Path) -> None:
    """N1: never ``run_failed:None None``."""
    out = write_arm(tmp_path / "out_a2")
    rows = [r for r in sacct_rows() if r["JobID"].split(".")[0] != HEAD_JOB]
    failures = arm_doc(out, rows)["failures"]
    assert failures == ["run_failed:no_sacct"], failures


def test_expect_input_without_equals_is_a_usage_error(tmp_path: Path) -> None:
    """N2: argparse error (exit 2), no traceback, no evidence."""
    out = write_arm(tmp_path / "out_a2")
    dest = tmp_path / "evidence" / "a2.json"
    args = arm_cli_args(out, dest)
    args[args.index(f"genome={SHA_GENOME}")] = "genome"
    res = cli(tmp_path, *args)
    assert res.returncode == 2, res.stdout + res.stderr
    _no_traceback(res)
    assert re.search(r"(?i)usage|error", res.stderr), res.stderr
    assert not dest.exists()


_SHORT_SACCT = r"""#!{python}
import sys
# one field short for every query (``-o`` lists 6 or 4 fields)
fields = sys.argv[sys.argv.index("-o") + 1].split(",")
print("|".join(["1900001"] + ["X"] * (len(fields) - 2)))
"""


def test_sacct_line_with_the_wrong_field_count_is_a_parse_error(tmp_path: Path) -> None:
    """N3: ``zip(strict=False)`` used to accept a short line silently."""
    out = write_arm(tmp_path / "out_a2")
    dest = tmp_path / "evidence" / "a2.json"
    res = cli(tmp_path, *arm_cli_args(out, dest), sacct_script=_SHORT_SACCT)
    assert res.returncode == 2, res.stdout + res.stderr
    _no_traceback(res)
    assert not dest.exists()
