"""``scripts/dev/bench_record.py``: the pure functions behind ``bench_run.json`` (M10-12).

docs/specs/M10-12-benchmark-standard.md ("Run validity", "Interface", edge cases "sacct gaps"
and "Nextflow task retries"); review of PR #300, findings 1-4, 7, 12, 13 and nit 18, with the
lead's decisions. The sbatch heredoc moved into this module (finding 13); the sandbox tests in
``test_m10_12_sbatch.py`` run it through the script, these call it directly.

Interface (settled by the test author from the reviewer's names; ``main`` is the implementer's):

* ``parse_sacct(row: str | None) -> dict | None``: one line of ``sacct -P -n -o
  JobIDRaw,ElapsedRaw,MaxRSS,TotalCPU,AllocCPUS,State,ExitCode`` → ``{elapsed_s: int,
  max_rss_bytes: int | None, total_cpu_s: float | None, alloc_cpus: int, state: str, ...}``
  (``state`` is the first word: ``CANCELLED by 12`` → ``CANCELLED``). ``None`` for no row and
  for an **incomplete** row (fewer than 7 fields, or ElapsedRaw, AllocCPUS or State empty or
  not a number where one is due): the caller retries. Never raises. MaxRSS empty or 0 → ``None``
  (never 0); units K/M/G/T are binary.
* ``trace_stats(rows) -> {"trace": {...}, "steps": {...}}``: ``rows`` are the
  ``csv.DictReader`` rows of ``trace.tsv``. ``trace`` is the record's ``trace`` object
  (``n_tasks n_cached n_failed n_retried hosts task_peak_rss_max_bytes failed_processes``);
  ``steps`` maps the ``measure_runtimes.py`` process name to ``{realtime_s, peak_rss_bytes,
  order}``, ``order`` ranking the steps by their first trace ``submit``. Never raises on a
  malformed row.
* ``invalid_reasons(*, variant, commit, node, sku, trace_rows, manifest, sacct, gippa_other,
  caches_before, caches_after) -> (reasons, warnings)``: ``manifest`` is the parsed
  ``run_manifest.json`` or ``None``; ``sacct`` the ``parse_sacct`` result or ``None``;
  ``gippa_other`` the record's ``cotenancy.gippa_other`` (``{job, name, cpus: int, node}``);
  the caches are the before/after file lists. Reasons as the spec names them, plus (lead,
  2026-10-24) ``manifest_missing``, ``commit_missing`` (empty or absent ``git_commit``) and
  ``sacct_missing``; ``contended_other_nodes`` is a warning.
* ``run_status(rc, sacct_state, reasons, forced=None) -> str``: ``forced`` > ``failed``
  (rc ∉ {0, 124, 137}, or an ``OUT_OF_MEMORY`` state whatever rc) > ``timeout`` (124, 137) >
  ``invalid`` (any reason) > ``ok``.

No statistical claim lives here (counts and maxima only), so there is no calibration test.
"""

from __future__ import annotations

import functools
import importlib.util
import sys
from types import ModuleType
from typing import Any

import pytest
from helpers.slurm_sandbox import REPO
from hypothesis import given, settings
from hypothesis import strategies as st

RECORD = REPO / "scripts" / "dev" / "bench_record.py"
COMMIT = "1f567901" + "ab" * 16
NODE = "wheat-05"
GIB = 1024**3


@functools.cache
def _mod() -> ModuleType:
    assert RECORD.is_file(), "scripts/dev/bench_record.py does not exist (PR #300, finding 13)"
    spec = importlib.util.spec_from_file_location("bench_record", RECORD)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules["bench_record"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def br() -> ModuleType:
    return _mod()


def trow(name: str, status: str = "COMPLETED", attempt: int | str = 1, *,
         submit: str = "2026-10-30 10:00:01.000", realtime: str = "3m",
         peak_rss: str = "2 GB", hostname: str = NODE) -> dict[str, str]:  # fmt: skip
    """One ``trace.tsv`` row as ``csv.DictReader`` gives it."""
    return {
        "task_id": "1", "hash": "ab/cdef01", "native_id": "1001", "name": name, "status": status,
        "exit": "0" if status == "COMPLETED" else "1", "attempt": str(attempt), "submit": submit,
        "duration": realtime, "realtime": realtime, "%cpu": "95.0%", "peak_rss": peak_rss,
        "peak_vmem": peak_rss, "rchar": "1 GB", "wchar": "1 GB", "hostname": hostname,
    }  # fmt: skip


SACCT_OK = {"elapsed_s": 3912, "max_rss_bytes": 21_600_000 * 1024, "total_cpu_s": 46944.0,
            "alloc_cpus": 16, "state": "COMPLETED"}  # fmt: skip


def good(**over: Any) -> dict[str, Any]:
    """The arguments of a clean core run; ``over`` replaces some."""
    kw: dict[str, Any] = {
        "variant": "core",
        "commit": COMMIT,
        "node": NODE,
        "sku": "CPU_SKU:6148",
        "trace_rows": [trow("MAPPING:MAP_PEPTIDES (k12)"), trow("OPERONS:OPERONS (k12)")],
        "manifest": {"git_commit": COMMIT, "git_dirty": False, "resume": False,
                     "steps": [{"name": "operons", "status": "completed"}]},
        "sacct": dict(SACCT_OK),
        "gippa_other": [],
        "caches_before": ["/c/containers/bakta.img", "/c/envs/env-old"],
        "caches_after": ["/c/containers/bakta.img", "/c/envs/env-old"],
    }  # fmt: skip
    kw.update(over)
    return kw


def reasons_of(br: ModuleType, **over: Any) -> list[str]:
    return list(br.invalid_reasons(**good(**over))[0])


def warnings_of(br: ModuleType, **over: Any) -> list[str]:
    return list(br.invalid_reasons(**good(**over))[1])


# ================================================================== run_status (finding 1, M1)


@pytest.mark.parametrize(
    ("rc", "state", "reasons", "forced", "want"),
    [
        (0, "COMPLETED", [], None, "ok"),
        (0, "COMPLETED", ["cached_tasks:3"], None, "invalid"),
        (124, "COMPLETED", [], None, "timeout"),
        (137, "CANCELLED", [], None, "timeout"),
        (124, "COMPLETED", ["cached_tasks:3"], None, "timeout"),
        (4, "FAILED", ["cached_tasks:3"], None, "failed"),
        (1, "FAILED", [], None, "failed"),
        # an OOM kill is failed whatever the exit code, 137 and 124 included (spec, lead)
        (137, "OUT_OF_MEMORY", [], None, "failed"),
        (124, "OUT_OF_MEMORY", [], None, "failed"),
        (0, "OUT_OF_MEMORY", [], None, "failed"),
        (137, "OUT_OF_MEMORY", ["cached_tasks:3"], None, "failed"),
        # no sacct row: no state, so the exit code decides; the reasons make rc 0 invalid
        (0, None, ["sacct_missing"], None, "invalid"),
        (137, None, ["sacct_missing"], None, "timeout"),
        # the SIGTERM trap forces timeout
        (143, None, [], "timeout", "timeout"),
        (0, "COMPLETED", [], "timeout", "timeout"),
    ],
    ids=["ok", "invalid", "timeout_124", "timeout_137", "timeout_over_invalid",
         "failed_over_invalid", "failed", "oom_137_is_failed", "oom_124_is_failed",
         "oom_0_is_failed", "oom_over_invalid", "no_state_rc0", "no_state_137", "forced_term",
         "forced_over_ok"],
)  # fmt: skip
def test_run_status(br: ModuleType, rc: int, state: str | None, reasons: list[str],
                    forced: str | None, want: str) -> None:  # fmt: skip
    assert br.run_status(rc, state, reasons, forced) == want


def test_failed_beats_timeout(br: ModuleType) -> None:
    """Criterion 5 precedence on a run with two causes: the timeout's 137 and an OOM state."""
    assert br.run_status(137, "OUT_OF_MEMORY", []) == "failed"
    assert br.run_status(124, "OUT_OF_MEMORY", ["image_or_env_created"]) == "failed"


_STATES = st.sampled_from([None, "COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY",
                           "NODE_FAIL"])  # fmt: skip
_RCS = st.sampled_from([0, 1, 2, 4, 124, 125, 137, 143, 255])


@settings(max_examples=200, deadline=None)
@given(_RCS, _STATES, st.lists(st.sampled_from(["cached_tasks:1", "alloc_cpus:18"]), max_size=2),
       st.sampled_from([None, "timeout"]))  # fmt: skip
def test_run_status_precedence_property(rc: int, state: str | None, reasons: list[str],
                                        forced: str | None) -> None:  # fmt: skip
    """forced > failed (rc not 0/124/137, or OOM) > timeout (124/137) > invalid > ok."""
    got = _mod().run_status(rc, state, reasons, forced)
    if forced:
        want = forced
    elif state == "OUT_OF_MEMORY" or rc not in (0, 124, 137):
        want = "failed"
    elif rc:
        want = "timeout"
    else:
        want = "invalid" if reasons else "ok"
    assert got == want


# ================================================================== parse_sacct (3, 4, nit 18)


def test_parse_sacct_full_row(br: ModuleType) -> None:
    s = br.parse_sacct("5100.0|3912|21600000K|13:02:24|16|COMPLETED|0:0\n")
    assert s is not None
    assert s["elapsed_s"] == 3912 and isinstance(s["elapsed_s"], int)
    assert s["max_rss_bytes"] == 21_600_000 * 1024
    assert s["total_cpu_s"] == pytest.approx(13 * 3600 + 2 * 60 + 24)
    assert s["alloc_cpus"] == 16 and s["state"] == "COMPLETED"


@pytest.mark.parametrize(
    ("rss", "want"),
    [("", None), ("0", None), ("0K", None), ("1024K", 1024**2), ("20.5G", int(20.5 * GIB)),
     ("3M", 3 * 1024**2)],
)  # fmt: skip
def test_parse_sacct_max_rss(br: ModuleType, rss: str, want: int | None) -> None:
    """MaxRSS empty (no 60-s sample yet) or 0 is null, never 0 (spec "sacct gaps", nit 18)."""
    s = br.parse_sacct(f"5100.0|3912|{rss}|13:02:24|16|COMPLETED|0:0")
    assert s is not None and s["max_rss_bytes"] == want, s


@pytest.mark.parametrize(
    ("cpu", "want"),
    [("13:02:24", 46944.0), ("02:24.500", 144.5), ("1-02:03:04", 93784.0),
     ("00:00:00", 0.0)],
)  # fmt: skip
def test_parse_sacct_total_cpu(br: ModuleType, cpu: str, want: float) -> None:
    s = br.parse_sacct(f"5100.0|3912|1000K|{cpu}|16|COMPLETED|0:0")
    assert s is not None and s["total_cpu_s"] == pytest.approx(want), s


def test_parse_sacct_empty_total_cpu_does_not_raise(br: ModuleType) -> None:
    s = br.parse_sacct("5100.0|3912|1000K||16|COMPLETED|0:0")
    assert s is None or s["total_cpu_s"] is None, s


def test_parse_sacct_state_is_the_first_word(br: ModuleType) -> None:
    s = br.parse_sacct("5100.0|60|1000K|00:01:00|16|CANCELLED by 4242|0:15")
    assert s is not None and s["state"] == "CANCELLED"
    s = br.parse_sacct("5100.0|60|1000K|00:01:00|16|OUT_OF_MEMORY|0:125")
    assert s is not None and s["state"] == "OUT_OF_MEMORY"


@pytest.mark.parametrize(
    "row",
    [None, "", "\n", "5100.0", "5100.0|3912", "5100.0||||||",
     "5100.0||21600000K|13:02:24|16|COMPLETING|0:0",
     "5100.0|3912|21600000K|13:02:24||COMPLETED|0:0",
     "5100.0|3912|21600000K|13:02:24|16||0:0",
     "5100.0|abc|21600000K|13:02:24|16|COMPLETED|0:0"],
    ids=["none", "empty", "newline", "one_field", "short", "all_empty", "elapsed_empty",
         "alloc_empty", "state_empty", "elapsed_not_a_number"],
)  # fmt: skip
def test_parse_sacct_incomplete_row_is_none(br: ModuleType, row: str | None) -> None:
    """No row or a row still being filled in: None (the caller retries), never an exception."""
    assert br.parse_sacct(row) is None


@settings(max_examples=300, deadline=None)
@given(st.lists(st.text(alphabet="0123456789:.-KMGT ABCDEFGHIJ_|", max_size=12), max_size=9))
def test_parse_sacct_never_raises(fields: list[str]) -> None:
    s = _mod().parse_sacct("|".join(fields))
    if s is not None:
        assert isinstance(s["elapsed_s"], int) and isinstance(s["alloc_cpus"], int)
        assert s["state"] and " " not in s["state"]
        assert s["max_rss_bytes"] is None or s["max_rss_bytes"] > 0


# ================================================================== trace_stats (findings 3, 7, 12)


def test_trace_stats_numbers(br: ModuleType) -> None:
    rows = [trow("ANNOTATE:BAKTA (k12)", realtime="10m", peak_rss="20.6 GB"),
            trow("OPERONS:OPERONS (k12)", realtime="3m", peak_rss="2 GB",
                 submit="2026-10-30 10:10:00.000")]  # fmt: skip
    out = br.trace_stats(rows)
    t, s = out["trace"], out["steps"]
    assert (t["n_tasks"], t["n_cached"], t["n_failed"], t["n_retried"]) == (2, 0, 0, 0)
    assert t["hosts"] == [NODE] and t["failed_processes"] == []
    assert t["task_peak_rss_max_bytes"] == pytest.approx(20.6 * GIB, rel=1e-6)
    assert set(s) == {"bakta", "operons"}
    assert s["bakta"]["realtime_s"] == 600 and s["operons"]["realtime_s"] == 180
    assert s["bakta"]["peak_rss_bytes"] == pytest.approx(20.6 * GIB, rel=1e-6)


def test_trace_stats_sums_the_tasks_of_a_step(br: ModuleType) -> None:
    rows = [trow("ANNOTATE:PROMOTER (TU_1)", realtime="1m", peak_rss="1 GB"),
            trow("ANNOTATE:PROMOTER (TU_2)", realtime="2m", peak_rss="3 GB")]  # fmt: skip
    s = br.trace_stats(rows)["steps"]["promoter"]
    assert s["realtime_s"] == 180
    assert s["peak_rss_bytes"] == pytest.approx(3 * GIB, rel=1e-6)


def _ordered(steps: dict[str, dict[str, Any]]) -> list[str]:
    return sorted(steps, key=lambda n: steps[n]["order"])


def test_step_order_is_the_first_submit_not_the_finish_order(br: ModuleType) -> None:
    """Nextflow writes a trace row when a task finishes; the order is by ``submit`` (finding 7,
    M13). PROMOTER's first task is submitted at 10:02, before OPERONS, though it finished last."""
    rows = [
        trow("OPERONS:OPERONS (k12)", submit="2026-10-30 10:05:00.000"),
        trow("ANNOTATE:BAKTA (k12)", submit="2026-10-30 10:00:00.000"),
        trow("ANNOTATE:PROMOTER (TU_9)", submit="2026-10-30 10:30:00.000"),
        trow("REPORT:REPORT (k12)", submit="2026-10-30 10:20:00.000"),
        trow("MAPPING:MAP_PEPTIDES (k12)", submit="2026-10-30 10:01:00.000"),
        trow("ANNOTATE:PROMOTER (TU_1)", submit="2026-10-30 10:02:00.000"),
    ]
    steps = br.trace_stats(rows)["steps"]
    assert _ordered(steps) == ["bakta", "map_peptides", "promoter", "operons", "report"]
    orders = [s["order"] for s in steps.values()]
    assert len(set(orders)) == len(orders) and all(isinstance(o, int) for o in orders), orders


def test_a_replaced_failed_attempt_is_not_a_failed_process(br: ModuleType) -> None:
    """The memory x 2 retry: attempt 1 FAILED, attempt 2 COMPLETED (finding 12, M3)."""
    rows = [trow("MAPPING:MAP_PEPTIDES (k12)"),
            trow("OPERONS:OPERONS (k12)", "FAILED", 1),
            trow("OPERONS:OPERONS (k12)", "COMPLETED", 2)]  # fmt: skip
    t = br.trace_stats(rows)["trace"]
    assert t["failed_processes"] == [] and t["n_failed"] == 0, t
    assert t["n_retried"] == 1 and t["n_tasks"] == 3, t


def test_a_task_failed_at_its_last_attempt_is_listed_once(br: ModuleType) -> None:
    rows = [trow("MAPPING:MAP_PEPTIDES (k12)"),
            trow("OPERONS:OPERONS (k12)", "FAILED", 1),
            trow("OPERONS:OPERONS (k12)", "FAILED", 2)]  # fmt: skip
    t = br.trace_stats(rows)["trace"]
    assert t["failed_processes"] == ["operons"] and t["n_failed"] == 1, t


def test_cached_rows_are_counted_not_measured(br: ModuleType) -> None:
    rows = [trow("OPERONS:OPERONS (k12)", realtime="3m", peak_rss="2 GB")]
    rows += [trow("MAPPING:DIGEST (k12)", "CACHED", realtime="-", peak_rss="-")] * 3
    out = br.trace_stats(rows)
    assert out["trace"]["n_cached"] == 3 and out["trace"]["n_tasks"] == 4
    assert set(out["steps"]) == {"operons"}
    assert out["trace"]["task_peak_rss_max_bytes"] == pytest.approx(2 * GIB, rel=1e-6)


def test_hosts_are_distinct_and_dashes_dropped(br: ModuleType) -> None:
    rows = [trow("A:X (1)", hostname=NODE), trow("A:Y (1)", hostname=NODE),
            trow("A:Z (1)", hostname="-"), trow("A:W (1)", hostname="wheat-07")]  # fmt: skip
    assert br.trace_stats(rows)["trace"]["hosts"] == [NODE, "wheat-07"]


def test_empty_trace(br: ModuleType) -> None:
    """gippa failed before Nextflow started: n_tasks 0, steps {} (spec "Failed run")."""
    out = br.trace_stats([])
    assert out["steps"] == {}
    t = out["trace"]
    assert (t["n_tasks"], t["n_cached"], t["n_failed"], t["n_retried"]) == (0, 0, 0, 0)
    assert t["hosts"] == [] and t["failed_processes"] == []
    assert t["task_peak_rss_max_bytes"] is None


@pytest.mark.parametrize(
    "rows",
    [[{}], [{"name": ""}], [{"status": "COMPLETED"}],
     [trow("OPERONS:OPERONS (k12)", attempt="-")], [trow("OPERONS:OPERONS (k12)", attempt="")],
     [trow("OPERONS:OPERONS (k12)", submit="-")], [trow("OPERONS:OPERONS (k12)", peak_rss="-")],
     [{k: v for k, v in trow("OPERONS:OPERONS (k12)").items() if k != "name"}]],
    ids=["empty_row", "empty_name", "no_name", "attempt_dash", "attempt_empty", "submit_dash",
         "rss_dash", "name_key_missing"],
)  # fmt: skip
def test_trace_stats_malformed_rows_do_not_raise(
    br: ModuleType, rows: list[dict[str, str]]
) -> None:
    """A crash here used to leave no bench_run.json at all (finding 3)."""
    out = br.trace_stats(rows)
    assert out["trace"]["n_tasks"] == len(rows)


_NAMES = st.sampled_from(["OPERONS:OPERONS (k12)", "ANNOTATE:BAKTA (k12)", "A:P (TU_1)",
                          "A:P (TU_2)", "", "X"])  # fmt: skip
_ROW = st.builds(
    trow,
    _NAMES,
    st.sampled_from(["COMPLETED", "FAILED", "CACHED", "ABORTED", ""]),
    st.sampled_from([1, 2, 3, "-", ""]),
    submit=st.sampled_from([f"2026-10-30 10:00:0{i}.000" for i in range(10)] + ["-"]),
    realtime=st.sampled_from(["1m", "30s", "-", ""]),
    peak_rss=st.sampled_from(["1 GB", "500 MB", "-", ""]),
    hostname=st.sampled_from([NODE, "-", "wheat-07"]),
)


@settings(max_examples=150, deadline=None)
@given(st.lists(_ROW, max_size=8))
def test_trace_stats_invariants(rows: list[dict[str, str]]) -> None:
    out = _mod().trace_stats(rows)
    t, s = out["trace"], out["steps"]
    assert t["n_tasks"] == len(rows)
    assert t["n_cached"] == sum(r["status"] == "CACHED" for r in rows)
    assert t["n_failed"] == len(t["failed_processes"]) <= sum(r["status"] == "FAILED" for r in rows)
    assert t["hosts"] == sorted(set(t["hosts"])) and "-" not in t["hosts"]
    orders = [v["order"] for v in s.values()]
    assert len(set(orders)) == len(orders)


# ================================================================== invalid_reasons


def test_clean_run_has_no_reason_and_no_warning(br: ModuleType) -> None:
    rs, ws = br.invalid_reasons(**good())
    assert list(rs) == [] and list(ws) == [], (rs, ws)


CASES: list[tuple[str, dict[str, Any], str]] = [
    ("cached_tasks", {"trace_rows": [trow("A:X (1)")] + [trow("A:D (1)", "CACHED")] * 3},
     "cached_tasks:3"),
    ("cached_steps", {"manifest": {"git_commit": COMMIT, "git_dirty": False, "resume": False,
                                   "steps": [{"name": "operons", "status": "cached"}]}},
     "cached_steps:1"),
    ("manifest_resume", {"manifest": {"git_commit": COMMIT, "git_dirty": False, "resume": True,
                                      "steps": []}}, "manifest_resume"),
    ("manifest_missing", {"manifest": None}, "manifest_missing"),
    ("commit_empty", {"manifest": {"git_commit": "", "git_dirty": False, "resume": False,
                                   "steps": []}}, "commit_missing"),
    ("commit_null", {"manifest": {"git_commit": None, "git_dirty": False, "resume": False,
                                  "steps": []}}, "commit_missing"),
    ("commit_absent", {"manifest": {"git_dirty": False, "resume": False, "steps": []}},
     "commit_missing"),
    ("commit_other", {"manifest": {"git_commit": "f" * 40, "git_dirty": False, "resume": False,
                                   "steps": []}}, "commit_mismatch"),
    ("commit_dirty", {"manifest": {"git_commit": COMMIT, "git_dirty": True, "resume": False,
                                   "steps": []}}, "commit_mismatch"),
    ("host_mismatch", {"trace_rows": [trow("A:X (1)", hostname="wheat-07")]}, "host_mismatch"),
    ("bakta_ran_in_core", {"trace_rows": [trow("ANNOTATE:BAKTA (k12)")]}, "bakta_ran_in_core"),
    ("same_node_gippa", {"gippa_other": [{"job": "777", "name": "gippa-run-pilot-a1", "cpus": 2,
                                        "node": NODE}]}, "same_node_gippa:777"),
    ("image_created", {"caches_after": ["/c/containers/bakta.img", "/c/envs/env-old",
                                        "/c/containers/new-tool.img"]}, "image_or_env_created"),
    ("alloc_cpus", {"sacct": dict(SACCT_OK, alloc_cpus=18)}, "alloc_cpus:18"),
    ("sacct_missing", {"sacct": None}, "sacct_missing"),
]  # fmt: skip


@pytest.mark.parametrize(("over", "reason"), [c[1:] for c in CASES], ids=[c[0] for c in CASES])
def test_each_invalid_reason(br: ModuleType, over: dict[str, Any], reason: str) -> None:
    """Criterion 5 and the lead's additions (findings 2, 4; M5: an empty git_commit)."""
    assert reason in reasons_of(br, **over)


def test_node_class_outside_the_pinned_sku(br: ModuleType) -> None:
    """Off-class runs are invalid, never summarised (lead, 2026-10-24, nit 20)."""
    rs = reasons_of(br, sku="CPU_SKU:9384X")
    assert any(r.startswith("node_class:") and "9384X" in r for r in rs), rs
    assert any(r.startswith("node_class:") for r in reasons_of(br, sku="")), "unknown SKU"


def test_bakta_in_a_standard_run_is_fine(br: ModuleType) -> None:
    assert reasons_of(br, variant="standard", trace_rows=[trow("ANNOTATE:BAKTA (k12)")]) == []


def test_a_short_or_full_hostname_is_the_same_node(br: ModuleType) -> None:
    assert reasons_of(br, trace_rows=[trow("A:X (1)", hostname=NODE + ".stanford.edu")]) == []
    assert reasons_of(br, node=NODE + ".stanford.edu") == []


def test_alloc_16_is_fine_and_no_sacct_is_not_an_alloc_reason(br: ModuleType) -> None:
    assert reasons_of(br, sacct=dict(SACCT_OK, alloc_cpus=16)) == []
    rs = reasons_of(br, sacct=None)
    assert not [r for r in rs if r.startswith("alloc_cpus")], rs


def test_all_reasons_are_recorded(br: ModuleType) -> None:
    rs = reasons_of(br, trace_rows=[trow("A:D (1)", "CACHED")] * 2,
                    manifest={"git_commit": COMMIT, "git_dirty": False, "resume": True,
                              "steps": []})  # fmt: skip
    assert {"cached_tasks:2", "manifest_resume"} <= set(rs), rs


@pytest.mark.parametrize(
    ("cpus", "flagged"), [(16, True), (8, True), (7, False), (2, False)]
)  # fmt: skip
def test_contended_other_nodes_from_8_cpus(br: ModuleType, cpus: int, flagged: bool) -> None:
    """A gippa-* job with >= 8 CPUs on another node: a warning, valid (example 7b; M4 at 8)."""
    other = [{"job": "778", "name": "gippa-diann-a2", "cpus": cpus, "node": "wheat-07"}]
    rs, ws = br.invalid_reasons(**good(gippa_other=other))
    assert list(rs) == [], rs
    assert ("contended_other_nodes" in ws) is flagged, ws


def test_same_node_job_of_any_size_invalidates(br: ModuleType) -> None:
    full = NODE + ".stanford.edu"
    other = [{"job": "779", "name": "gippa-x", "cpus": 16, "node": NODE},
             {"job": "780", "name": "gippa-y", "cpus": 1, "node": full}]  # fmt: skip
    rs = reasons_of(br, gippa_other=other)
    assert {"same_node_gippa:779", "same_node_gippa:780"} <= set(rs), rs
