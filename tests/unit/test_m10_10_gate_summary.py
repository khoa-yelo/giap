"""``validation/pilot/cascadia_convert.py`` ``gate`` and ``summary`` (M10-10a criteria 1 and 2).

Spec: docs/specs/M10-10-a3-run.md, "Entry gate", "``cascadia_convert.py``", example 1 and the
summary rules 1-4. The 12-row fixture is converted with the real M10-09 converter
(``convert_cascadia``) and every ``counts`` field is worked out by hand below. Interface
conventions: ``tests/helpers/a3_run_m10_10.py``.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.a3_run_m10_10 import (
    ACQUISITION,
    FAIL_REASONS,
    RUNS,
    SEQ_LIKE,
    TWELVE,
    WEIGHTS_SHA,
    feasibility_doc,
    load_cascadia_convert,
    sha256,
    strings,
    write_assets,
    write_feasibility,
    write_ssl,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.config.model import InputConfig
from giap.io.converters.cascadia import convert_cascadia
from giap.io.issues import Issue

# ------------------------------------------------------------------ criterion 1: the gate


def _gate(tmp: Path, doc: dict[str, Any] | str | None, *, weights: str = WEIGHTS_SHA,
          capsys: pytest.CaptureFixture[str] | None = None) -> tuple[int, str]:  # fmt: skip
    mod = load_cascadia_convert()
    feas = tmp / "a3_feasibility.json"
    if doc is not None:
        write_feasibility(feas, doc)
    assets = write_assets(tmp / "a3_assets.tsv", weights)
    rc = mod.main(["gate", "--feasibility", str(feas), "--assets", str(assets)])
    text = ""
    if capsys is not None:
        out = capsys.readouterr()
        text = out.out + out.err
    return rc, text


def test_gate_passes_on_a_pass_verdict(tmp_path: Path) -> None:
    assert _gate(tmp_path, feasibility_doc())[0] == 0


@pytest.mark.parametrize("reason", FAIL_REASONS)
def test_example1_gate_exits_3_on_every_fail_reason(
    tmp_path: Path, reason: str, capsys: pytest.CaptureFixture[str]
) -> None:
    rc, text = _gate(tmp_path, feasibility_doc(verdict="fail", reason=reason), capsys=capsys)
    assert rc == 3, text
    assert f"M10-08 verdict fail ({reason}): M10-10 is skipped; open M10-11a" in text, text


def test_gate_missing_file_exits_2(tmp_path: Path) -> None:
    assert _gate(tmp_path, None)[0] == 2


@pytest.mark.parametrize(
    "case",
    ["schema_version_2", "schema_version_str", "run_wt2", "not_json", "list", "no_verdict",
     "verdict_unknown", "pass_with_reason", "no_weights"],
)  # fmt: skip
def test_gate_malformed_exits_2(tmp_path: Path, case: str) -> None:
    doc: dict[str, Any] | str = feasibility_doc()
    if case == "schema_version_2":
        doc = feasibility_doc(schema_version=2)
    elif case == "schema_version_str":
        doc = feasibility_doc(schema_version="1")
    elif case == "run_wt2":
        doc = feasibility_doc(run="WT_2")
    elif case == "not_json":
        doc = '{"schema_version": 1, "verdict": "pass"'
    elif case == "list":
        doc = "[]\n"
    elif case == "no_verdict":
        assert isinstance(doc, dict)
        doc.pop("verdict")
    elif case == "verdict_unknown":
        doc = feasibility_doc(verdict="maybe")
    elif case == "pass_with_reason":
        doc = feasibility_doc(reason="too_few_peptides")
    elif case == "no_weights":
        assert isinstance(doc, dict)
        doc.pop("weights_sha256")
    assert _gate(tmp_path, doc)[0] == 2


@pytest.mark.parametrize("weights", ["0" * 64, ""], ids=["other", "empty_assets_row"])
def test_gate_weights_hash_differs_from_assets_exits_2(tmp_path: Path, weights: str) -> None:
    """The committed a3_assets.tsv carries no weights sha256 until M10-08b: that is not a pass."""
    assert _gate(tmp_path, feasibility_doc(), weights=weights)[0] == 2


def test_gate_assets_without_a_weights_row_exits_2(tmp_path: Path) -> None:
    mod = load_cascadia_convert()
    feas = write_feasibility(tmp_path / "f.json", feasibility_doc())
    assets = tmp_path / "a3_assets.tsv"
    write_assets(assets)
    lines = [ln for ln in assets.read_text().splitlines() if not ln.startswith("cascadia_weights")]
    assets.write_text("\n".join(lines) + "\n")
    assert mod.main(["gate", "--feasibility", str(feas), "--assets", str(assets)]) == 2


_VERDICTS = st.sampled_from(["pass", "fail", "PASS", "", "unknown"])
_REASONS = st.sampled_from([None, *FAIL_REASONS, "other"])


@settings(
    max_examples=60, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(verdict=_VERDICTS, reason=_REASONS, schema=st.sampled_from([1, 2, 0]),
       same_weights=st.booleans())  # fmt: skip
def test_property_gate_exit_code(
    tmp_path_factory: pytest.TempPathFactory, verdict: str, reason: str | None, schema: int,
    same_weights: bool,
) -> None:  # fmt: skip
    """0 only for a well-formed pass; 3 only for a well-formed fail with a known reason."""
    tmp = tmp_path_factory.mktemp("gate")
    doc = feasibility_doc(verdict=verdict, reason=reason, schema_version=schema)
    rc, _ = _gate(tmp, doc, weights=WEIGHTS_SHA if same_weights else "0" * 64)
    well_formed = schema == 1
    if well_formed and same_weights and verdict == "pass" and reason is None:
        assert rc == 0
    elif well_formed and verdict == "fail" and reason in FAIL_REASONS:
        assert rc in (3, 2)  # a fail verdict never passes; a weights mismatch may come first
        if same_weights:
            assert rc == 3
    else:
        assert rc == 2


# ------------------------------------------------------------------ criterion 2: summarize


def _twelve(tmp: Path) -> tuple[pd.DataFrame, list[Issue], list[Path]]:
    paths = [write_ssl(tmp / f"{r}.ssl", r, TWELVE[r]) for r in RUNS]
    res = convert_cascadia(paths, search_db=None, cfg=InputConfig())
    assert len(res.peptides) == 12  # one output row per input row (M10-09 D4)
    return res.peptides, list(res.issues), paths


SSL_ROWS = {"WT_1": 4, "WT_2": 4, "WT_3": 4}


def _summarize(peptides: pd.DataFrame, issues: list[Issue], *, ssl_rows: dict[str, int] = SSL_ROWS,
               runs: tuple[str, ...] = RUNS) -> dict[str, Any]:  # fmt: skip
    mod = load_cascadia_convert()
    doc = mod.summarize(dict(ssl_rows), peptides, list(issues), list(runs))
    assert isinstance(doc, dict)
    return json.loads(json.dumps(doc))  # JSON-serialisable; keys as written to a3.convert.json


def test_summarize_counts_by_hand(tmp_path: Path) -> None:
    peptides, issues, _ = _twelve(tmp_path)
    doc = _summarize(peptides, issues)
    c = doc["counts"]
    assert c["rows"] == 12
    assert c["rows_per_run"] == {"WT_1": 4, "WT_2": 4, "WT_3": 4}
    # PEPTLDEK MPEPTLDEKR SAMPLLER ACDEFGHK LLLKR CAMPLLERK V*31 GGGGGGG
    assert c["distinct_peptides"] == 8
    assert c["distinct_peptides_7_30"] == 6  # not LLLKR (5) and not V*31 (31)
    assert c["modified_rows"] == 2  # M[+15.995]… and C[+57.021]…
    assert {str(k): v for k, v in c["charge_counts"].items()} == {"2": 7, "3": 4, "4": 1}
    assert c["score_min"] == pytest.approx(0.85)
    assert c["score_median"] == pytest.approx(0.915)  # (0.91 + 0.92) / 2 over the 12 rows
    assert c["score_max"] == pytest.approx(0.99)
    assert c["score_at_or_below_default_cutoff"] == 0
    assert doc["failures"] == []


def test_summarize_counts_rows_at_the_cutoff(tmp_path: Path) -> None:
    rows = {r: list(v) for r, v in TWELVE.items()}
    rows["WT_2"][0] = ("ACDEFGHK", 2, 0.8)  # kept by the converter, counted here
    paths = [write_ssl(tmp_path / f"{r}.ssl", r, rows[r]) for r in RUNS]
    res = convert_cascadia(paths, search_db=None, cfg=InputConfig())
    doc = _summarize(res.peptides, list(res.issues))
    assert doc["counts"]["score_at_or_below_default_cutoff"] == 1
    assert doc["counts"]["score_min"] == pytest.approx(0.8)


def test_summarize_issue_table(tmp_path: Path) -> None:
    peptides, issues, _ = _twelve(tmp_path)
    doc = _summarize(peptides, issues)
    assert doc["issues"]["stripped"] == {"severity": "info", "rows": 2}
    assert doc["issues"]["short_peptide"]["severity"] == "warn"
    assert doc["failures"] == [], "warn codes are listed, not failed (rule 3)"


def _has_rule(doc: dict[str, Any], n: int) -> bool:
    return any(str(f).startswith(f"rule {n}") for f in doc["failures"])


def test_rule1_table_rows_differ_from_ssl_rows(tmp_path: Path) -> None:
    peptides, issues, _ = _twelve(tmp_path)
    doc = _summarize(peptides, issues, ssl_rows={"WT_1": 5, "WT_2": 4, "WT_3": 4})
    assert _has_rule(doc, 1), doc["failures"]
    doc = _summarize(peptides.iloc[1:], issues)
    assert _has_rule(doc, 1), doc["failures"]


@pytest.mark.parametrize("case", ["missing_run", "extra_run", "renamed_run"])
def test_rule2_runs_differ_from_acquisition(tmp_path: Path, case: str) -> None:
    peptides, issues, _ = _twelve(tmp_path)
    if case == "missing_run":
        doc = _summarize(peptides, issues, runs=(*RUNS, "WT_4"))
    elif case == "extra_run":
        doc = _summarize(peptides, issues, runs=RUNS[:2])
    else:
        p = peptides.copy()
        p["run"] = p["run"].replace({"WT_3": "WT_3.mzML"})
        doc = _summarize(p, issues)
    assert _has_rule(doc, 2), doc["failures"]


def test_rule3_an_error_issue_fails(tmp_path: Path) -> None:
    peptides, issues, _ = _twelve(tmp_path)
    doc = _summarize(peptides, [*issues, Issue("run_in_several_inputs", "error", "x")])
    assert _has_rule(doc, 3), doc["failures"]
    assert doc["issues"]["run_in_several_inputs"]["severity"] == "error"


@pytest.mark.parametrize("case", ["database_row", "q_value", "pep", "is_decoy", "intensity"])
def test_rule4_database_fields_fail(tmp_path: Path, case: str) -> None:
    peptides, issues, _ = _twelve(tmp_path)
    p = peptides.copy()
    if case == "database_row":
        p["source"] = p["source"].astype(object)
        p.loc[p.index[3], "source"] = "database"
    else:
        p[case] = 0.0 if case != "is_decoy" else False
    assert _has_rule(_summarize(p, issues), 4)


def test_summary_json_has_no_sequence_and_is_stable(tmp_path: Path) -> None:
    peptides, issues, _ = _twelve(tmp_path)
    mod = load_cascadia_convert()
    a = mod.summarize(dict(SSL_ROWS), peptides, list(issues), list(RUNS))
    b = mod.summarize(dict(SSL_ROWS), peptides.copy(), list(issues), list(RUNS))
    ta = json.dumps(a, indent=2, sort_keys=True)
    assert ta == json.dumps(b, indent=2, sort_keys=True)
    bad = [s for s in strings(a) if re.fullmatch(SEQ_LIKE, s)]
    assert bad == [], f"sequence-like strings in the summary: {bad}"


# ------------------------------------------------------------------ the summary CLI


def _cli(tmp: Path, *, sabotage: str | None = None) -> tuple[int, dict[str, Any], dict[str, Path]]:
    from click.testing import CliRunner

    from giap.cli import cli

    a3 = tmp / "A3"
    paths = [write_ssl(a3 / r / f"{r}.ssl", r, TWELVE[r]) for r in RUNS]
    sums = a3 / "ssl.sha256"
    sums.write_text("".join(f"{sha256(p)}  {p}\n" for p in paths))
    pep, log = tmp / "inputs" / "a3.tsv.gz", tmp / "inputs" / "a3.convert.log"
    pep.parent.mkdir()
    argv = ["convert", "--from", "cascadia", *map(str, paths), "--out", str(pep),
            "--source", "denovo"]  # fmt: skip
    res = CliRunner().invoke(cli, argv)
    assert res.exit_code == 0, res.output
    log.write_text("".join(ln + "\n" for ln in res.output.splitlines()
                           if ln.startswith("giap convert:")))  # fmt: skip
    if sabotage == "extra_ssl_row":
        with paths[0].open("a") as fh:
            fh.write("\t".join(["/x/WT_1.mzML", "9", "2", "PEPTLDEK", "UNKNOWN", "0.9", "1", "1",
                                "1"]) + "\n")  # fmt: skip
        sums.write_text("".join(f"{sha256(p)}  {p}\n" for p in paths))
    out = tmp / "a3.convert.json"
    mod = load_cascadia_convert()
    rc = mod.main(["summary", "--ssl-sums", str(sums), "--peptides", str(pep), "--log", str(log),
                   "--acquisition", str(ACQUISITION), "--out", str(out)])  # fmt: skip
    doc = json.loads(out.read_text()) if out.is_file() else {}
    return rc, doc, {"ssl": a3, "peptides": pep, "log": log, "out": out, "sums": sums}


def test_summary_cli_writes_a3_convert_json(tmp_path: Path) -> None:
    rc, doc, p = _cli(tmp_path)
    assert rc == 0, doc.get("failures")
    text = p["out"].read_text()
    assert text == json.dumps(json.loads(text), indent=2, sort_keys=True) + "\n"
    assert doc["failures"] == []
    for r in RUNS:
        entry = doc["ssl"][r]
        path = p["ssl"] / r / f"{r}.ssl"
        assert entry["sha256"] == sha256(path)
        assert entry["bytes"] == path.stat().st_size
        assert entry["rows"] == 4
    assert doc["outputs"]["peptides"]["sha256"] == sha256(p["peptides"])
    assert doc["outputs"]["peptides"]["rows"] == 12, "data rows of the table, short ones included"
    assert doc["outputs"]["log"]["sha256"] == sha256(p["log"])
    assert doc["counts"]["rows_per_run"] == {"WT_1": 4, "WT_2": 4, "WT_3": 4}
    assert doc["issues"]["stripped"] == {"severity": "info", "rows": 2}
    assert set(doc["giap"]) >= {"version", "commit", "dirty"}
    assert isinstance(doc["command"], list) and "cascadia" in doc["command"]
    assert set(doc["job"]) >= {"id", "node", "start", "end"}
    bad = [s for s in strings(doc) if re.fullmatch(SEQ_LIKE, s)]
    assert bad == [], f"sequence-like strings in a3.convert.json: {bad}"


def test_summary_cli_is_byte_stable_outside_job(tmp_path: Path) -> None:
    _, a, _ = _cli(tmp_path / "a")
    _, b, _ = _cli(tmp_path / "b")
    for d in (a, b):
        d.pop("job")
        for e in [*d["ssl"].values(), *d["outputs"].values()]:
            e.pop("path")  # tmp_path differs between the two runs, the rest must not
        d.pop("command")
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)


def test_summary_cli_exits_1_on_a_rule_failure(tmp_path: Path) -> None:
    rc, doc, _ = _cli(tmp_path, sabotage="extra_ssl_row")
    assert rc == 1
    assert any(str(f).startswith("rule 1") for f in doc["failures"]), doc["failures"]


# ------------------------------------------------------------------ properties

_SEQ = st.text(alphabet="ACDEFGHKLMNPQRSTVWY", min_size=3, max_size=34)
_ROW = st.tuples(_SEQ, st.integers(1, 5), st.floats(0.801, 1.0, allow_nan=False))


@settings(
    max_examples=40, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(rows=st.fixed_dictionaries({r: st.lists(_ROW, min_size=1, max_size=8) for r in RUNS}))
def test_property_summary_invariants(
    tmp_path_factory: pytest.TempPathFactory, rows: dict[str, list[tuple[str, int, float]]]
) -> None:
    tmp = tmp_path_factory.mktemp("prop")
    paths = [write_ssl(tmp / f"{r}.ssl", r, rows[r]) for r in RUNS]
    res = convert_cascadia(paths, search_db=None, cfg=InputConfig())
    n = {r: len(v) for r, v in rows.items()}
    doc = _summarize(res.peptides, list(res.issues), ssl_rows=n)
    c = doc["counts"]
    assert c["rows"] == sum(n.values()) == sum(c["rows_per_run"].values())
    assert c["rows_per_run"] == n
    assert c["distinct_peptides_7_30"] <= c["distinct_peptides"] <= c["rows"]
    assert c["score_min"] <= c["score_median"] <= c["score_max"]
    assert sum(c["charge_counts"].values()) == c["rows"]
    assert c["score_at_or_below_default_cutoff"] == 0
    assert not _has_rule(doc, 1) and not _has_rule(doc, 2) and not _has_rule(doc, 4)
    assert not [s for s in strings(doc) if re.fullmatch(SEQ_LIKE, s)]
