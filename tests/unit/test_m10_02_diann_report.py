"""``validation/pilot/diann_report.py`` (``runs``, ``check``, ``sample``, ``provenance``) and
``scripts/slurm/diann_report.sbatch`` (docs/specs/M10-02-diann-searches.md).

Spec: "validation/pilot/diann_report.py", "check rules", "Row selection", "summary.json",
examples 10-12, edge cases "Partial report", "Empty input", "Ties", "Stratum quota",
acceptance criteria 6-8; lead decision 2 (stratified parquet sample). Khoa's decision of
2026-10-23: the sample is never committed; ``diann_report.sbatch`` writes it under
``$GIAP_ROOT/pilot/recorded/diann/2.3.1/<arm>/``.

Fixtures are parquet files of a few dozen rows written in ``tmp_path`` with pyarrow, using the
converter's column names (``helpers.pilot_kit_m10.report_rows``). No DIA-NN output is read.

Conventions settled by the test author:

* The CLI runs as ``python -I validation/pilot/diann_report.py <command> …``.
* The selection rule is also a pure function,
  ``select_rows(run, decoy, precursor_id, n=1000, d=100) -> list[int]``: three equal-length
  sequences (one entry per report row in file order; ``decoy`` truthy for a decoy) and the
  0-based file-order indices of the sample, ascending. The Hypothesis properties call it.
* ``check.json`` ``failures`` are strings; each names its rule: ``footer`` (rule 1, as in
  example 10's "parquet footer unreadable"), the missing column (rule 2), ``decoy`` (rule 4),
  the run label or ``Run`` value concerned (rule 3), ``sha256`` (rule 5).
* ``schema.tsv`` may start with the header ``column\\tarrow_type\\tnullable``.

The sample is a deterministic selection, not an estimator: there is no statistical claim and so
no calibration test (spec §8.4).
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
from collections import Counter
from pathlib import Path
from typing import Any

import pytest
from helpers.pilot_kit_m10 import (
    DIANN_REPORT,
    M10_SCRIPTS,
    RUNS,
    M10Sandbox,
    acquisition_doc,
    load_diann_report_module,
    m10_sandbox,
    report_rows,
    run_report,
    sha256,
    write_json,
    write_report,
    write_sums,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

pa = pytest.importorskip("pyarrow")
pq = pytest.importorskip("pyarrow.parquet")

THREE = {r: 6 for r in RUNS}
TWO = {r: 2 for r in RUNS}


def _key(run: str, pid: str) -> str:
    return hashlib.sha256(f"{run}\t{pid}".encode()).hexdigest()


def _reference(run: list[str], decoy: list[int], pid: list[str], n: int, d: int) -> list[int]:
    """The spec's rule, written out independently (an oracle for the property tests)."""
    runs = sorted(set(run))
    r = len(runs)
    quota: dict[tuple[str, int], int] = {}
    for kind, total in ((0, n - d), (1, d)):
        base, rem = divmod(total, r) if r else (0, 0)
        for j, name in enumerate(runs):
            quota[(name, kind)] = base + (1 if j < rem else 0)
    keyed = sorted((_key(run[i], pid[i]), i) for i in range(len(run)))
    taken: set[int] = set()
    for (name, kind), q in quota.items():
        stratum = [i for _, i in keyed if run[i] == name and int(bool(decoy[i])) == kind]
        taken.update(stratum[:q])
    for _, i in keyed:
        if len(taken) >= n:
            break
        taken.add(i)
    return sorted(taken)


def _acq(tmp: Path, runs: tuple[str, ...] | list[str] = RUNS, **kw: Any) -> Path:
    return write_json(tmp / "acquisition.json", acquisition_doc(list(runs), **kw))


def _good(tmp: Path, rows: list[dict[str, Any]] | None = None, **drop: Any) -> tuple[Path, Path]:
    report = write_report(tmp / "A1_report.parquet", rows or report_rows(THREE, TWO), **drop)
    return report, write_sums(tmp / "outputs.sha256", report)


def _check(tmp: Path, report: Path, sums: Path, acq: Path | None = None):  # type: ignore[no-untyped-def]
    out = tmp / "check.json"
    res = run_report("check", str(report), "--acquisition", str(acq or _acq(tmp)),
                     "--sums", str(sums), "--out", str(out))  # fmt: skip
    doc = json.loads(out.read_text()) if out.is_file() else None
    return res, doc


def _failures(doc: dict[str, Any] | None) -> list[str]:
    assert doc is not None, "check.json was not written"
    return [str(f) for f in doc["failures"]]


# ------------------------------------------------------------ runs (criterion 8)


def test_runs_prints_the_labels_in_file_order(tmp_path: Path) -> None:
    order = ["wt_rep3", "wt_rep1", "wt_rep2"]
    res = run_report("runs", str(_acq(tmp_path, order)))
    assert res.returncode == 0, res.stderr
    assert res.stdout.split() == order


@pytest.mark.parametrize("mode", ["ddaPASEF", "mixed", "ms1_only", "unknown", "inconsistent"])
def test_runs_exits_3_for_other_modes(tmp_path: Path, mode: str) -> None:
    res = run_report("runs", str(_acq(tmp_path, mode=mode)))
    assert res.returncode == 3, (res.returncode, res.stdout, res.stderr)
    msg = res.stdout + res.stderr
    assert "§17" in msg and mode in msg and "ADR-0051" in msg, msg
    assert not res.stdout.strip() or mode in res.stdout, "no run labels on a refusal"


@pytest.mark.parametrize("variant", ["synchro-PASEF", "slice-PASEF", "midia-PASEF"])
def test_runs_exits_3_for_a_variant_in_notes(tmp_path: Path, variant: str) -> None:
    res = run_report("runs", str(_acq(tmp_path, notes={"wt_rep1": f"{variant} scheme"})))
    assert res.returncode == 3, (res.returncode, res.stdout, res.stderr)
    assert variant.split("-")[0].lower() in (res.stdout + res.stderr).lower()


@pytest.mark.parametrize("kind", ["missing", "unparsable", "not-an-object", "no-runs"])
def test_runs_exits_2_for_a_bad_file(tmp_path: Path, kind: str) -> None:
    path = tmp_path / "acquisition.json"
    if kind == "unparsable":
        path.write_text("{")
    elif kind == "not-an-object":
        path.write_text("[1, 2]")
    elif kind == "no-runs":
        path.write_text(json.dumps({"dataset_mode": "diaPASEF"}))
    res = run_report("runs", str(path))
    assert res.returncode == 2, (res.returncode, res.stdout, res.stderr)


def test_no_arguments_is_a_usage_error() -> None:
    assert run_report().returncode == 2
    assert run_report("frobnicate").returncode == 2


# ------------------------------------------------------------ check (criterion 6, example 10)


def test_check_accepts_a_well_formed_report(tmp_path: Path) -> None:
    report, sums = _good(tmp_path)
    res, doc = _check(tmp_path, report, sums)
    assert res.returncode == 0, (res.stdout, res.stderr, doc)
    assert doc and doc["ok"] is True and doc["failures"] == []
    assert doc["rows"] == 3 * 8
    assert set(doc["columns"]) >= {"Run", "Decoy", "Precursor.Id"}
    assert doc["runs"] == {r: 8 for r in RUNS}
    assert doc["decoy_rows"] == 6
    rows = report_rows(THREE, TWO)
    assert math.isclose(doc["max_q"], max(r["Q.Value"] for r in rows))
    assert doc["pyarrow"] == pa.__version__


def test_example10_a_truncated_report_fails_on_its_footer(tmp_path: Path) -> None:
    report, _ = _good(tmp_path)
    data = report.read_bytes()
    report.write_bytes(data[:-100])
    sums = write_sums(tmp_path / "outputs.sha256", report)
    res, doc = _check(tmp_path, report, sums)
    assert res.returncode == 1, (res.stdout, res.stderr)
    assert doc and doc["ok"] is False
    assert any("footer" in f for f in _failures(doc)), doc


def test_an_empty_or_zero_row_report_fails_rule1(tmp_path: Path) -> None:
    empty = tmp_path / "A1_report.parquet"
    empty.write_bytes(b"")
    res, doc = _check(tmp_path, empty, write_sums(tmp_path / "outputs.sha256", empty))
    assert res.returncode == 1 and doc and doc["ok"] is False
    zero = write_report(tmp_path / "A1_report.parquet", [])
    res, doc = _check(tmp_path, zero, write_sums(tmp_path / "outputs.sha256", zero))
    assert res.returncode == 1 and doc and doc["ok"] is False


@pytest.mark.parametrize("column", ["Decoy", "PG.MaxLFQ", "Protein.Group", "Modified.Sequence"])
def test_a_missing_required_column_is_named(tmp_path: Path, column: str) -> None:
    report, sums = _good(tmp_path, drop=[column])
    res, doc = _check(tmp_path, report, sums)
    assert res.returncode == 1 and doc and doc["ok"] is False
    assert any(column in f for f in _failures(doc)), doc


def test_either_q_value_column_satisfies_rule2(tmp_path: Path) -> None:
    report, sums = _good(tmp_path, drop=["Global.Q.Value"])
    res, doc = _check(tmp_path, report, sums)
    assert res.returncode == 0, doc
    report, sums = _good(tmp_path, drop=["Global.Q.Value", "Q.Value"])
    res, doc = _check(tmp_path, report, sums)
    assert res.returncode == 1 and any("Q.Value" in f for f in _failures(doc)), doc


def test_a_report_with_no_decoy_rows_fails_rule4(tmp_path: Path) -> None:
    report, sums = _good(tmp_path, report_rows(THREE, {}))
    res, doc = _check(tmp_path, report, sums)
    assert res.returncode == 1 and doc and doc["ok"] is False
    (only,) = _failures(doc)
    assert "decoy" in only.lower(), only


def test_a_run_missing_from_the_report_is_named(tmp_path: Path) -> None:
    rows = report_rows({"wt_rep1": 6, "wt_rep2": 6}, {"wt_rep1": 2, "wt_rep2": 2})
    report, sums = _good(tmp_path, rows)
    res, doc = _check(tmp_path, report, sums)
    assert res.returncode == 1
    (only,) = _failures(doc)
    assert "wt_rep3" in only, only


def test_an_unexpected_run_value_is_named(tmp_path: Path) -> None:
    rows = report_rows({**THREE, "blank_07": 3}, TWO)
    report, sums = _good(tmp_path, rows)
    res, doc = _check(tmp_path, report, sums)
    assert res.returncode == 1
    (only,) = _failures(doc)
    assert "blank_07" in only, only


def test_a_run_value_with_a_path_and_d_suffix_matches_its_label(tmp_path: Path) -> None:
    rows = report_rows(THREE, TWO)
    for r in rows:
        r["Run"] = f"/data/pxd078459/{r['Run']}.d"
    report, sums = _good(tmp_path, rows)
    res, doc = _check(tmp_path, report, sums)
    assert res.returncode == 0, doc


def test_a_report_that_disagrees_with_outputs_sha256_fails_rule5(tmp_path: Path) -> None:
    report, sums = _good(tmp_path)
    sums.write_text(f"{'0' * 64}  {report.name}\n")
    res, doc = _check(tmp_path, report, sums)
    assert res.returncode == 1
    (only,) = _failures(doc)
    assert "sha256" in only, only


def test_a_report_absent_from_outputs_sha256_fails_rule5(tmp_path: Path) -> None:
    report, sums = _good(tmp_path)
    sums.write_text("")
    res, doc = _check(tmp_path, report, sums)
    assert res.returncode == 1 and any("sha256" in f for f in _failures(doc))


def test_check_usage_errors_are_exit_2(tmp_path: Path) -> None:
    report, _ = _good(tmp_path)
    assert run_report("check", str(report)).returncode == 2


# ------------------------------------------------------------ sample: the rule (criterion 7)


def _hand_rows(skip_wt3_targets: bool) -> list[dict[str, Any]]:
    """File order: for j in 0..6, runs WT_2, WT_3, WT_1; rows j >= 5 are decoys; id PEP<j>."""
    rows = []
    for j in range(7):
        for run in ("WT_2", "WT_3", "WT_1"):
            if skip_wt3_targets and run == "WT_3" and 1 <= j <= 4:
                continue
            rows.append({"Run": run, "Decoy": int(j >= 5), "Precursor.Id": f"PEP{j}"})
    return rows


def _full(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for i, r in enumerate(rows):
        seq = f"PEPTIDE{r['Precursor.Id'][3:]}K"
        out.append({
            **r, "Stripped.Sequence": seq, "Modified.Sequence": seq, "Precursor.Charge": 2,
            "Protein.Ids": "P1", "Protein.Group": "P1", "Q.Value": 0.001 * i,
            "Global.Q.Value": 0.001 * i, "PG.MaxLFQ": 1.0 + i, "Fixture.Row": i,
        })  # fmt: skip
    return out


# Written out by hand from the rule (sha256 prefixes in the comments). n = 10, d = 2, 3 runs:
# targets 8 = 3/3/2, decoys 2 = 1/1/0 (remainders to WT_1, WT_2, the first runs sorted).
#   WT_1 targets 2:18cf 5:3f0c 8:7b32 11:223e 14:cf9f -> 2, 11, 5
#   WT_2 targets 0:e8c2 3:24ae 6:7e82 9:5ba4 12:401b  -> 3, 12, 9
#   WT_3 targets 1:5eef 4:3588 7:94c0 10:0228 13:9db1 -> 10, 4
#   WT_1 decoys 17:5f70 20:ff12 -> 17;  WT_2 decoys 15:a2ea 18:3d28 -> 18;  WT_3 decoys: quota 0
HAND_EXPECTED = [2, 3, 4, 5, 9, 10, 11, 12, 17, 18]
# Example 12 in small: WT_3 has one target (1:5eef) for a quota of 2. The shortfall of 1 is the
# smallest key not yet taken over all strata: 12:20ce (a WT_3 decoy).
HAND_SHORTFALL_EXPECTED = [1, 2, 3, 4, 7, 8, 9, 12, 13, 14]


@pytest.mark.parametrize(
    ("skip", "want"), [(False, HAND_EXPECTED), (True, HAND_SHORTFALL_EXPECTED)],
    ids=["quotas", "shortfall"],
)  # fmt: skip
def test_select_rows_matches_the_hand_example(skip: bool, want: list[int]) -> None:
    mod = load_diann_report_module()
    rows = _hand_rows(skip)
    got = mod.select_rows([r["Run"] for r in rows], [r["Decoy"] for r in rows],
                          [r["Precursor.Id"] for r in rows], n=10, d=2)  # fmt: skip
    assert list(got) == want


@pytest.mark.parametrize(
    ("skip", "want"), [(False, HAND_EXPECTED), (True, HAND_SHORTFALL_EXPECTED)],
    ids=["quotas", "shortfall"],
)  # fmt: skip
def test_sample_cli_writes_the_hand_example_in_file_order(
    tmp_path: Path, skip: bool, want: list[int]
) -> None:
    report = write_report(tmp_path / "A1_report.parquet", _full(_hand_rows(skip)))
    out = tmp_path / "sample"
    res = run_report("sample", str(report), "--out", str(out), "--n", "10", "--decoys", "2")
    assert res.returncode == 0, (res.stdout, res.stderr)
    table = pq.read_table(out / "report.sample.parquet")
    assert table.column("Fixture.Row").to_pylist() == want
    assert table.schema.equals(pq.read_schema(report)), "the sample keeps the source schema"
    assert (
        pq.ParquetFile(out / "report.sample.parquet")
        .metadata.row_group(0)
        .column(0)
        .compression.lower()
        == "zstd"
    )


def test_example11_the_full_size_quotas(tmp_path: Path) -> None:
    """Example 11: 3 runs, 600,000 target and 12,000 decoy rows -> 300 + 34/33/33 per run."""
    n_t, n_d = 200_000, 4_000
    run, decoy, pid = [], [], []
    for j in range(n_t + n_d):
        for r in ("WT_3", "WT_1", "WT_2"):
            run.append(r)
            decoy.append(int(j >= n_t))
            pid.append(f"PEP{j}")
    table = pa.table({"Run": run, "Decoy": decoy, "Precursor.Id": pid,
                      "Stripped.Sequence": pid, "Modified.Sequence": pid,
                      "Precursor.Charge": [2] * len(run), "Protein.Ids": ["P1"] * len(run),
                      "Protein.Group": ["P1"] * len(run), "Global.Q.Value": [0.0] * len(run),
                      "PG.MaxLFQ": [1.0] * len(run)})  # fmt: skip
    report = tmp_path / "A1_report.parquet"
    pq.write_table(table, report)
    out = tmp_path / "sample"
    res = run_report("sample", str(report), "--out", str(out))
    assert res.returncode == 0, (res.stdout, res.stderr)
    s = pq.read_table(out / "report.sample.parquet", columns=["Run", "Decoy"]).to_pylist()
    assert len(s) == 1000
    counts = Counter((r["Run"], r["Decoy"]) for r in s)
    assert counts == {("WT_1", 0): 300, ("WT_2", 0): 300, ("WT_3", 0): 300,
                      ("WT_1", 1): 34, ("WT_2", 1): 33, ("WT_3", 1): 33}  # fmt: skip
    summary = json.loads((out / "summary.json").read_text())
    assert summary["sample_rows"] == 1000 and summary["duplicate_keys"] == 0
    assert summary["rows"] == 3 * (n_t + n_d) and summary["decoy_rows"] == 3 * n_d


def test_example12_a_short_stratum_is_filled_from_the_others(tmp_path: Path) -> None:
    rows = report_rows(
        {"WT_1": 400, "WT_2": 400, "WT_3": 120}, {"WT_1": 50, "WT_2": 50, "WT_3": 50}
    )
    report = write_report(tmp_path / "A1_report.parquet", rows)
    out = tmp_path / "sample"
    assert run_report("sample", str(report), "--out", str(out)).returncode == 0
    got = pq.read_table(out / "report.sample.parquet").column("Fixture.Row").to_pylist()
    assert len(got) == 1000
    assert sum(1 for i in got if rows[i]["Run"] == "WT_3" and not rows[i]["Decoy"]) == 120
    want = _reference([r["Run"] for r in rows], [r["Decoy"] for r in rows],
                      [r["Precursor.Id"] for r in rows], 1000, 100)  # fmt: skip
    assert got == want


def test_a_report_under_n_rows_is_returned_whole(tmp_path: Path) -> None:
    rows = report_rows(THREE, TWO)
    report = write_report(tmp_path / "A1_report.parquet", rows)
    out = tmp_path / "sample"
    assert run_report("sample", str(report), "--out", str(out)).returncode == 0
    got = pq.read_table(out / "report.sample.parquet")
    assert got.column("Fixture.Row").to_pylist() == list(range(len(rows)))
    assert json.loads((out / "summary.json").read_text())["sample_rows"] == len(rows)


def test_sampling_twice_gives_identical_rows(tmp_path: Path) -> None:
    report = write_report(tmp_path / "A1_report.parquet", report_rows({r: 60 for r in RUNS}, TWO))
    a, b = tmp_path / "a", tmp_path / "b"
    for d in (a, b):
        assert run_report("sample", str(report), "--out", str(d), "--n", "30",
                          "--decoys", "3").returncode == 0  # fmt: skip
    assert pq.read_table(a / "report.sample.parquet").equals(
        pq.read_table(b / "report.sample.parquet")
    )


def test_duplicate_keys_are_counted_and_broken_by_file_order(tmp_path: Path) -> None:
    rows = report_rows(THREE, TWO)
    rows[5]["Precursor.Id"] = rows[2]["Precursor.Id"]
    rows[5]["Run"] = rows[2]["Run"]
    report = write_report(tmp_path / "A1_report.parquet", rows)
    out = tmp_path / "sample"
    assert run_report("sample", str(report), "--out", str(out)).returncode == 0
    assert json.loads((out / "summary.json").read_text())["duplicate_keys"] == 1


def test_the_precursor_id_fallback_is_deterministic(tmp_path: Path) -> None:
    rows = [{k: v for k, v in r.items() if k != "Precursor.Id"}
            for r in report_rows({r: 40 for r in RUNS}, TWO)]  # fmt: skip
    report = write_report(tmp_path / "A1_report.parquet", rows)
    out = tmp_path / "sample"
    res = run_report("sample", str(report), "--out", str(out), "--n", "20", "--decoys", "3")
    assert res.returncode == 0, res.stderr
    s = pq.read_table(out / "report.sample.parquet").to_pylist()
    assert len(s) == 20 and sum(r["Decoy"] for r in s) == 3


def test_sample_writes_the_schema_and_summary(tmp_path: Path) -> None:
    rows = report_rows(THREE, TWO)
    report = write_report(tmp_path / "A1_report.parquet", rows)
    out = tmp_path / "sample"
    assert run_report("sample", str(report), "--out", str(out)).returncode == 0
    schema = pq.read_schema(report)
    lines = (out / "schema.tsv").read_text().splitlines()
    if lines and lines[0].split("\t") == ["column", "arrow_type", "nullable"]:
        lines = lines[1:]
    cells = [ln.split("\t") for ln in lines]
    assert [c[0] for c in cells] == schema.names
    assert [c[1] for c in cells] == [str(f.type) for f in schema]
    assert [c[2].lower() for c in cells] == [str(f.nullable).lower() for f in schema]
    s = json.loads((out / "summary.json").read_text())
    assert s["rows"] == len(rows) and s["columns"] == len(schema.names)
    assert s["rows_per_run"] == {r: 8 for r in RUNS}
    assert s["decoy_rows"] == 6
    assert s["distinct_precursors"] == len({r["Precursor.Id"] for r in rows})
    q01 = {r["Stripped.Sequence"] for r in rows if not r["Decoy"] and r["Global.Q.Value"] <= 0.01}
    assert s["distinct_stripped_q01"] == len(q01)
    assert math.isclose(s["max_q"], max(r["Q.Value"] for r in rows))
    assert s["report_sha256"] == sha256(report) and s["report_bytes"] == report.stat().st_size
    assert s["sample_rule"] == "strata(run,decoy) sha256(Run\tPrecursor.Id) n=1000 d=100"
    assert s["pyarrow"] == pa.__version__


def test_sample_refuses_a_structurally_invalid_report(tmp_path: Path) -> None:
    report, _ = _good(tmp_path)
    report.write_bytes(report.read_bytes()[:-100])
    out = tmp_path / "sample"
    res = run_report("sample", str(report), "--out", str(out))
    assert res.returncode == 1, (res.stdout, res.stderr)
    assert not (out / "report.sample.parquet").exists()


def test_diann_report_writes_nothing_into_the_repository(tmp_path: Path) -> None:
    report, sums = _good(tmp_path)
    before = os.popen(f"git -C {DIANN_REPORT.parent} status --porcelain").read()
    _check(tmp_path, report, sums)
    run_report("sample", str(report), "--out", str(tmp_path / "s"))
    assert os.popen(f"git -C {DIANN_REPORT.parent} status --porcelain").read() == before


# ------------------------------------------------------------ Hypothesis: invariants of the rule

_strata = st.lists(
    st.tuples(
        st.sampled_from(["WT_1", "WT_2", "WT_3", "B_9"]),
        st.booleans(),
        st.integers(min_value=0, max_value=40),
    ),
    min_size=1,
    max_size=60,
)


@settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(rows=_strata, n=st.integers(min_value=1, max_value=40), dfrac=st.floats(0, 1))
def test_select_rows_invariants(rows: list[tuple[str, bool, int]], n: int, dfrac: float) -> None:
    mod = load_diann_report_module()
    run = [r for r, _, _ in rows]
    decoy = [int(d) for _, d, _ in rows]
    pid = [f"P{k}" for _, _, k in rows]
    d = int(dfrac * n)
    got = list(mod.select_rows(run, decoy, pid, n=n, d=d))
    assert got == sorted(set(got)), "ascending file order, no repeats"
    assert all(0 <= i < len(rows) for i in got)
    assert len(got) == min(n, len(rows)), "n rows, or the whole report"
    assert got == _reference(run, decoy, pid, n, d)
    # Within each stratum the taken rows are the smallest (key, index): a prefix.
    for s in {(r, x) for r, x in zip(run, decoy, strict=True)}:
        members = sorted((_key(run[i], pid[i]), i) for i in range(len(rows))
                         if (run[i], decoy[i]) == s)  # fmt: skip
        taken = [i for _, i in members if i in set(got)]
        assert taken == [i for _, i in members][: len(taken)], s
    assert list(mod.select_rows(run, decoy, pid, n=n, d=d)) == got, "deterministic"


@settings(max_examples=40, deadline=None)
@given(r=st.integers(1, 7), n=st.integers(1, 600), dfrac=st.floats(0, 1))
def test_quotas_split_evenly_with_the_remainder_to_the_first_runs(
    r: int, n: int, dfrac: float
) -> None:
    """With every stratum large enough, each gets exactly its quota (rule step 2)."""
    mod = load_diann_report_module()
    d = int(dfrac * n)
    names = [f"R{j:02d}" for j in range(r)]
    big = n  # every stratum can fill its quota
    run, decoy, pid = [], [], []
    for name in reversed(names):  # file order unrelated to the sorted order
        for kind in (0, 1):
            for k in range(big):
                run.append(name)
                decoy.append(kind)
                pid.append(f"K{k}")
    got = mod.select_rows(run, decoy, pid, n=n, d=d)
    counts = Counter((run[i], decoy[i]) for i in got)
    for kind, total in ((0, n - d), (1, d)):
        base, rem = divmod(total, r)
        for j, name in enumerate(names):
            assert counts.get((name, kind), 0) == base + (1 if j < rem else 0), (name, kind)


# ------------------------------------------------------------ diann_report.sbatch in the sandbox


@pytest.fixture
def m10(tmp_path: Path) -> M10Sandbox:
    return m10_sandbox(tmp_path)


def _searched(
    m10: M10Sandbox, arm: str, rows: list[dict[str, Any]], truncate: bool = False
) -> Path:
    out = m10.out(arm)
    report = write_report(out / f"{arm}_report.parquet", rows)
    if truncate:
        report.write_bytes(report.read_bytes()[:-100])
    write_sums(out / "outputs.sha256", report)
    return out


def _sample_dir(m10: M10Sandbox, arm: str) -> Path:
    return m10.root / "pilot" / "recorded" / "diann" / "2.3.1" / arm


@pytest.mark.parametrize("arm", ["A1", "A2"])
def test_report_job_checks_then_samples_under_giap_root(m10: M10Sandbox, arm: str) -> None:
    out = _searched(m10, arm, report_rows({r: 30 for r in RUNS}, TWO))
    res = m10.run("diann_report", arm)
    assert res.returncode == 0, f"{res.stdout}\n{res.stderr}"
    check = json.loads((out / "check.json").read_text())
    assert check["ok"] is True, check
    sample = _sample_dir(m10, arm)
    for name in ("report.sample.parquet", "schema.tsv", "summary.json"):
        assert (sample / name).is_file(), f"{name} missing under {sample} (Khoa 2026-10-23)"
    assert m10.repo_is_clean() == "", "the sample is never written into the checkout"


def test_example10_report_job_on_a_truncated_report(m10: M10Sandbox) -> None:
    out = _searched(m10, "A1", report_rows(THREE, TWO), truncate=True)
    res = m10.run("diann_report", "A1")
    assert res.returncode == 1, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    check = json.loads((out / "check.json").read_text())
    assert check["ok"] is False and any("footer" in f for f in check["failures"]), check
    assert not _sample_dir(m10, "A1").exists() or not any(_sample_dir(m10, "A1").iterdir())
    assert not (out / "sample").exists()


def test_report_job_without_a_report_fails(m10: M10Sandbox) -> None:
    res = m10.run("diann_report", "A1")
    assert res.returncode != 0
    assert not _sample_dir(m10, "A1").exists()


def test_report_job_house_rules() -> None:
    path = M10_SCRIPTS["diann_report"]
    assert path.is_file(), f"{path} does not exist (M10-02 deliverable 5)"
    text = path.read_text()

    def directive(key: str) -> str | None:
        m = re.search(rf"^#SBATCH\s+--{key}(?:[= ](\S+))?", text, re.M)
        return (m.group(1) or "") if m else None

    assert (directive("job-name") or "").startswith("giap-diann-report")
    assert directive("partition") == "normal" and directive("qos") == "normal"
    assert directive("cpus-per-task") == "4"
    assert (directive("mem") or "").upper().rstrip("B") == "16000M"  # 4 x 4000 MB (S1)
    assert directive("time") == "01:00:00"
    out = directive("output") or ""
    assert "/local/logs/" in out or "$GIAP_ROOT" in out
    code = "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("#"))
    assert "scripts/farmshare_env.sh" in code
    assert re.search(r"^set -euo pipefail$", code, re.M)
    assert "tests/data/recorded" not in code, "Khoa 2026-10-23: samples stay on FarmShare"


def test_the_kit_stays_within_the_r14_cpu_cap() -> None:
    def cpus(name: str) -> int:
        m = re.search(r"^#SBATCH\s+--cpus-per-task[= ](\d+)", M10_SCRIPTS[name].read_text(), re.M)
        assert m, name
        return int(m.group(1))

    total = 2 * cpus("diann_pilot") + 2 * cpus("diann_report") + cpus("diann_image")
    assert total <= 128, total
