"""``validation/pilot/diann_convert.py`` summary (M10-03 criterion 6).

Spec: docs/specs/M10-03-diann-convert.md, "validation/pilot/diann_convert.py summary" and the
edge case "q-values at the cut-off".

Interface fixed here (the spec names ``summarize(report_cols, peptides, abundance, issues)``):

* ``report_cols``: a pandas DataFrame of report columns, at least ``Run``, ``Decoy`` and
  ``Global.Q.Value`` (as ``pq.read_table(report, columns=…).to_pandas()`` gives them);
* ``peptides``: the converted table **read back** from the written TSV (``read_peptide_table``);
* ``abundance``: the MaxLFQ table read back as a DataFrame (``protein_id``, ``maxlfq``);
* ``issues``: a list of ``gippa.io.issues.Issue``;
* returns a dict with ``counts`` (the fields below), ``issues`` (code → ``{"severity", "rows"}``)
  and ``failures`` (a list of strings; empty when every summary rule holds).

The CLI: ``python -I diann_convert.py summary --arm A --report R --peptides P --abundance B
--log L --out J`` writes ``J`` (exit 0) or exits 1 when ``failures`` is not empty.

The 12-row fixture is hand-made; each expected count below is computed by hand in the comments.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as pq
import pytest
from helpers.diann_convert_m10 import (
    SUMMARY,
    gippa_convert,
    load_summary_module,
    row,
    sha256,
    table,
    write_parquet,
)

from gippa.io.issues import Issue
from gippa.io.peptides import read_peptide_table


def _r(
    run: str,
    pep: str,
    ids: str,
    q: float | None,
    decoy: int = 0,
    inten: float | None = 1000.0,
    lfq: float = 500.0,
) -> dict[str, Any]:
    return row(
        Run=run,
        Stripped_Sequence=pep,
        Modified_Sequence=pep,
        Protein_Ids=ids,
        Protein_Group=ids,
        Global_Q_Value=q,
        PEP=0.01,
        Decoy=decoy,
        Precursor_Charge=2,
        Precursor_Normalised=inten,
        PG_MaxLFQ=lfq,
    )


ROWS = [
    _r("WT_1", "PEPTMIDEK", "P0A7V0", 0.0021, inten=15234.5),  # 1 target, q ≤ 0.01
    _r("WT_1", "PLPTMIDDK", "P0A7V0", 0.4, decoy=1, lfq=900.0),  # 2 decoy
    _r("WT_1", "KEDIMTPEP", "P0A7V0", None, decoy=1),  # 3 decoy, q NA
    _r("WT_1", "LLGNVLVR", "Cont_P00761", 0.001),  # 4 contaminant only
    _r("WT_2", "AVLGDLVR", "Cont_P00722;P00722", 0.005),  # 5 mixed contaminant
    _r("WT_2", "QCROSSINGK", "P0A7V1", 0.01000004),  # 6 > 0.01, written "0.01": a crossing
    _r("WT_2", "NOCROSSK", "P0A7V2", 0.0100001),  # 7 written "0.0100001": no crossing
    _r("WT_2", "ATCUTOFFK", "P0A7V3", 0.01),  # 8 at the cut-off already
    _r("WT_3", "PEPTMIDEK", "P0A7V0", 0.002, inten=0.0),  # 9 intensity 0 → NA; repeat peptide
    _r("WT_3", "ABOVEQK", "P0A7V1", 0.02),  # 10 q > 0.01
    _r("WT_3", "NOINTENSK", "JMJJOG_01234", 0.003, inten=None, lfq=0.0),  # 11 intensity NA
    _r("WT_3", "DECOYTHREEK", "P0A7V9", 0.9, decoy=1, lfq=900.0),  # 12 decoy
]
COUNTS = {
    "rows": 12,
    "target_rows": 9,
    "decoy_rows": 3,  # rows 2, 3, 12
    "decoy_rows_q_na": 1,  # row 3
    "contaminant_only_rows": 1,  # row 4
    "mixed_contaminant_rows": 1,  # row 5
    "rows_per_run": {"WT_1": 4, "WT_2": 4, "WT_3": 4},
    # targets with written q ≤ 0.01: rows 1, 4, 5, 6, 8, 9, 11; row 9 repeats row 1's peptide
    "distinct_peptides_target_q01": 6,
    "q_crossings": 1,  # row 6
    "intensity_na_rows": 2,  # rows 9, 11
    # non-decoy groups with a positive value: P0A7V0, Cont_P00761, Cont_P00722;P00722, P0A7V1,
    # P0A7V2, P0A7V3 (JMJJOG_01234 has 0; P0A7V9 is decoy-only)
    "abundance_groups": 6,
}


@pytest.fixture(scope="module")
def mod():  # type: ignore[no-untyped-def]
    return load_summary_module()


@pytest.fixture
def converted(tmp_path: Path) -> dict[str, Any]:
    """The fixture report converted by this checkout's CLI and read back."""
    rep = write_parquet(tmp_path / "A1_report.parquet", table(ROWS))
    pep, ab = tmp_path / "a1.tsv.gz", tmp_path / "a1_abundance.tsv.gz"
    res = gippa_convert(rep, pep, ab, "other")
    assert res.returncode == 0, res.stderr
    log = tmp_path / "a1.convert.log"
    log.write_text(res.stderr)
    peptides, issues = read_peptide_table(pep)
    assert not [i for i in issues if i.severity == "error"], issues
    abundance = pd.read_csv(ab, sep="\t")
    report = pq.read_table(rep).to_pandas()
    return {
        "report": report,
        "peptides": peptides,
        "abundance": abundance,
        "paths": {"report": rep, "peptides": pep, "abundance": ab, "log": log},
    }


def _summ(mod, c: dict[str, Any], issues: list[Issue] | None = None, **over: Any):  # type: ignore[no-untyped-def]
    args = {
        "report_cols": c["report"],
        "peptides": c["peptides"],
        "abundance": c["abundance"],
        "issues": issues or [],
        **over,
    }
    return mod.summarize(args["report_cols"], args["peptides"], args["abundance"], args["issues"])


def test_the_summary_script_exists() -> None:
    assert SUMMARY.is_file(), "validation/pilot/diann_convert.py does not exist (M10-03)"


@pytest.mark.parametrize("field", sorted(COUNTS))
def test_counts_by_hand(mod, converted: dict[str, Any], field: str) -> None:  # type: ignore[no-untyped-def]
    out = _summ(mod, converted)
    assert out["counts"][field] == COUNTS[field], (field, out["counts"].get(field))


def test_the_good_fixture_has_no_failure(mod, converted: dict[str, Any]) -> None:  # type: ignore[no-untyped-def]
    assert _summ(mod, converted)["failures"] == []


def test_q_crossing_is_the_writer_rounding(converted: dict[str, Any]) -> None:
    """The premise of ``q_crossings``: row 6 reads back as exactly 0.01."""
    df = converted["peptides"].set_index("peptide")
    assert df.loc["QCROSSINGK", "q_value"] == 0.01
    assert df.loc["NOCROSSK", "q_value"] > 0.01


# ------------------------------------------------------------------ each rule fails when broken


def test_rule1_rows_differ_from_the_report(mod, converted: dict[str, Any]) -> None:  # type: ignore[no-untyped-def]
    short = converted["peptides"].iloc[:-1].reset_index(drop=True)
    assert _summ(mod, converted, peptides=short)["failures"]


def test_rule2_decoy_count_differs_from_the_report(mod, converted: dict[str, Any]) -> None:  # type: ignore[no-untyped-def]
    rep = converted["report"].copy()
    rep.loc[0, "Decoy"] = 1  # the report now has 4 decoys, the table 3
    assert _summ(mod, converted, report_cols=rep)["failures"]


def test_rule2_no_decoy_fails(mod, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    rows = [{**r, "Decoy": 0} for r in ROWS]
    rep = write_parquet(tmp_path / "r.parquet", table(rows))
    pep, ab = tmp_path / "p.tsv", tmp_path / "ab.tsv"
    assert gippa_convert(rep, pep, ab, "other").returncode == 0
    peptides, _ = read_peptide_table(pep)
    out = mod.summarize(pq.read_table(rep).to_pandas(), peptides, pd.read_csv(ab, sep="\t"), [])
    assert out["counts"]["decoy_rows"] == 0
    assert out["failures"]


def test_rule3_an_error_issue_fails_and_a_warn_is_listed(mod, converted: dict[str, Any]) -> None:  # type: ignore[no-untyped-def]
    warn = Issue("decoy_without_protein_ids", "warn", "x (2 rows)", 2, ("1", "2"))
    out = _summ(mod, converted, issues=[warn])
    assert out["failures"] == []
    assert out["issues"] == {"decoy_without_protein_ids": {"severity": "warn", "rows": 2}}
    err = Issue("bad_type:PEP", "error", "PEP: invalid value (1 rows)", 1, ("3",))
    assert _summ(mod, converted, issues=[warn, err])["failures"]


def test_rule4_run_keys_differ_from_the_report(mod, converted: dict[str, Any]) -> None:  # type: ignore[no-untyped-def]
    rep = converted["report"].copy()
    rep["Run"] = rep["Run"].astype(str)
    rep.loc[0, "Run"] = "WT_4"
    assert _summ(mod, converted, report_cols=rep)["failures"]


# ------------------------------------------------------------------ the CLI


def _cli(c: dict[str, Any], out: Path, arm: str = "A1") -> subprocess.CompletedProcess[str]:
    p = c["paths"]
    cmd = [
        sys.executable,
        "-I",
        str(SUMMARY),
        "summary",
        "--arm",
        arm,
        "--report",
        str(p["report"]),
        "--peptides",
        str(p["peptides"]),
        "--abundance",
        str(p["abundance"]),
        "--log",
        str(p["log"]),
        "--out",
        str(out),
    ]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=300, check=False)


def test_cli_writes_convert_json(converted: dict[str, Any], tmp_path: Path) -> None:
    assert SUMMARY.is_file(), "validation/pilot/diann_convert.py does not exist (M10-03)"
    p = converted["paths"]
    p["log"].write_text(
        p["log"].read_text() + "gippa convert: warn: decoy_without_protein_ids: "
        "decoy rows without protein IDs (2 rows) (e.g. rows 1, 2)\n"
    )
    out = tmp_path / "a1.convert.json"
    res = _cli(converted, out)
    assert res.returncode == 0, res.stderr
    doc = json.loads(out.read_text())
    for key in (
        "arm",
        "search_db",
        "report",
        "outputs",
        "gippa",
        "command",
        "job",
        "issues",
        "counts",
        "failures",
    ):
        assert key in doc, key
    assert doc["arm"] == "A1" and doc["search_db"] == "other"
    assert doc["failures"] == []
    assert doc["counts"] == {**doc["counts"], **COUNTS}
    assert doc["report"]["sha256"] == sha256(p["report"])
    assert doc["report"]["rows"] == 12
    assert doc["report"]["bytes"] == p["report"].stat().st_size
    assert doc["issues"]["decoy_without_protein_ids"] == {"severity": "warn", "rows": 2}
    shas = json.dumps(doc["outputs"])
    for name in ("peptides", "abundance", "log"):
        assert sha256(p[name]) in shas, name
    assert '"PEPTMIDEK"' not in out.read_text(), "convert.json holds no row (Q-K1)"


def test_cli_exits_1_on_a_failed_rule(converted: dict[str, Any], tmp_path: Path) -> None:
    assert SUMMARY.is_file(), "validation/pilot/diann_convert.py does not exist (M10-03)"
    p = converted["paths"]
    p["log"].write_text("gippa convert: error: bad_type:PEP: PEP: invalid value (1 rows)\n")
    out = tmp_path / "a1.convert.json"
    res = _cli(converted, out)
    assert res.returncode == 1, (res.returncode, res.stderr)
