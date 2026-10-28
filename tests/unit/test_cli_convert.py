"""``gippa convert --from generic|diann`` end to end (issue M1-07; GIPPA_SPEC.md §5.1, §5.4)."""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from gippa.cli import cli
from gippa.exitcodes import ExitCode

DIANN_ROWS = [
    {
        "Run": run,
        "Protein.Group": "P1;P2",
        "Protein.Ids": "P1;P2",
        "Stripped.Sequence": "PEPTIDEK",
        "Modified.Sequence": "PEPTIDEK",
        "Precursor.Charge": 2,
        "Q.Value": 0.002,
        "Global.Q.Value": 0.004,
        "PEP": 0.01,
        "Precursor.Quantity": 1400.0,
        "Precursor.Normalised": norm,
        "PG.MaxLFQ": lfq,
    }
    for run, norm, lfq in (("r1", 1500.0, 10.0), ("r2", 0.0, 30.0), ("r3", 900.0, 0.0))
]


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


def _diann(tmp_path: Path, kind: str = "tsv", rows: list[dict[str, object]] | None = None) -> Path:
    df = pd.DataFrame(rows if rows is not None else DIANN_ROWS)
    if kind == "parquet":
        p = tmp_path / "report.parquet"
        df.to_parquet(p, engine="pyarrow", index=False)
    else:
        p = tmp_path / "report.tsv"
        df.to_csv(p, sep="\t", index=False, lineterminator="\n")
    return p


def _generic(tmp_path: Path, text: str) -> Path:
    p = tmp_path / "in.tsv"
    p.write_text(text, encoding="utf-8")
    return p


def test_example9_generic_end_to_end(runner: CliRunner, tmp_path: Path) -> None:
    src = _generic(tmp_path, "peptide\tis_decoy\nPEPT[+79.966]IDEK\t1\n")
    out = tmp_path / "peptides.tsv"
    result = runner.invoke(cli, ["convert", "--from", "generic", str(src), "--out", str(out)])
    assert result.exit_code == ExitCode.SUCCESS, result.output
    assert out.read_text() == "peptide\tis_decoy\nPEPTIDEK\ttrue\n"
    assert "stripped" in result.stderr


def test_generic_error_exits_1_and_writes_nothing(runner: CliRunner, tmp_path: Path) -> None:
    src = _generic(tmp_path, "peptide\tq_value\nPEPTIDEK\tabc\n")
    out = tmp_path / "peptides.tsv"
    result = runner.invoke(cli, ["convert", "--from", "generic", str(src), "--out", str(out)])
    assert result.exit_code == ExitCode.INVALID_INPUT
    assert not out.exists()
    assert "bad_type:q_value" in result.stderr


def test_diann_tsv_end_to_end(runner: CliRunner, tmp_path: Path) -> None:
    out = tmp_path / "peptides.tsv"
    ab = tmp_path / "abundance.tsv"
    argv = ["convert", "--from", "diann", str(_diann(tmp_path)), "--out", str(out)]
    result = runner.invoke(cli, [*argv, "--abundance-out", str(ab)])
    assert result.exit_code == ExitCode.SUCCESS, result.output
    assert out.read_text() == (
        "peptide\tmodified_peptide\tprotein_ids\tq_value\tpep\tintensity\trun\tcharge\tsource\n"
        "PEPTIDEK\tPEPTIDEK\tP1;P2\t0.004\t0.01\t1500\tr1\t2\tdatabase\n"
        "PEPTIDEK\tPEPTIDEK\tP1;P2\t0.004\t0.01\tNA\tr2\t2\tdatabase\n"
        "PEPTIDEK\tPEPTIDEK\tP1;P2\t0.004\t0.01\t900\tr3\t2\tdatabase\n"
    )
    assert ab.read_text() == "protein_id\tmaxlfq\nP1;P2\t20\n"
    assert "--abundance-type maxlfq" in result.output


def test_diann_parquet_and_tsv_give_identical_bytes(runner: CliRunner, tmp_path: Path) -> None:
    outs = []
    for kind in ("tsv", "parquet"):
        d = tmp_path / kind
        d.mkdir()
        out = d / "peptides.tsv"
        argv = ["convert", "--from", "diann", str(_diann(d, kind)), "--out", str(out)]
        result = runner.invoke(cli, argv)
        assert result.exit_code == ExitCode.SUCCESS, result.output
        outs.append(out.read_bytes())
    assert outs[0] == outs[1]


def test_output_is_deterministic(runner: CliRunner, tmp_path: Path) -> None:
    src = _diann(tmp_path)
    blobs = []
    for i in range(2):
        out = tmp_path / f"o{i}.tsv.gz"
        ab = tmp_path / f"a{i}.tsv"
        argv = ["convert", "--from", "diann", str(src), "--out", str(out)]
        result = runner.invoke(cli, [*argv, "--abundance-out", str(ab)])
        assert result.exit_code == ExitCode.SUCCESS, result.output
        blobs.append((out.read_bytes(), ab.read_bytes()))
    assert blobs[0] == blobs[1]


def test_gz_output_is_gzipped(runner: CliRunner, tmp_path: Path) -> None:
    out = tmp_path / "peptides.tsv.gz"
    argv = ["convert", "--from", "diann", str(_diann(tmp_path)), "--out", str(out)]
    result = runner.invoke(cli, argv)
    assert result.exit_code == ExitCode.SUCCESS, result.output
    text = gzip.decompress(out.read_bytes()).decode()
    assert text.startswith("peptide\tmodified_peptide\t")


def test_diann_missing_column_exits_1_no_output(runner: CliRunner, tmp_path: Path) -> None:
    rows = [{k: v for k, v in r.items() if k != "Stripped.Sequence"} for r in DIANN_ROWS]
    out = tmp_path / "peptides.tsv"
    argv = ["convert", "--from", "diann", str(_diann(tmp_path, rows=rows)), "--out", str(out)]
    result = runner.invoke(cli, argv)
    assert result.exit_code == ExitCode.INVALID_INPUT
    assert not out.exists()
    assert "missing_column:Stripped.Sequence" in result.stderr


def test_diann_empty_report_exits_1(runner: CliRunner, tmp_path: Path) -> None:
    src = tmp_path / "report.tsv"
    pd.DataFrame(columns=list(DIANN_ROWS[0])).to_csv(src, sep="\t", index=False)
    out = tmp_path / "peptides.tsv"
    result = runner.invoke(cli, ["convert", "--from", "diann", str(src), "--out", str(out)])
    assert result.exit_code == ExitCode.INVALID_INPUT
    assert not out.exists()


def test_diann_with_source_denovo_is_error(runner: CliRunner, tmp_path: Path) -> None:
    out = tmp_path / "peptides.tsv"
    argv = ["convert", "--from", "diann", str(_diann(tmp_path)), "--out", str(out)]
    result = runner.invoke(cli, [*argv, "--source", "denovo"])
    assert result.exit_code == ExitCode.INVALID_INPUT
    assert not out.exists()


def test_diann_search_db_flag_adds_column(runner: CliRunner, tmp_path: Path) -> None:
    out = tmp_path / "peptides.tsv"
    argv = ["convert", "--from", "diann", str(_diann(tmp_path)), "--out", str(out)]
    result = runner.invoke(cli, [*argv, "--search-db", "other"])
    assert result.exit_code == ExitCode.SUCCESS, result.output
    header, first = out.read_text().splitlines()[:2]
    assert header.endswith("\tsource\tsearch_db")
    assert first.endswith("\tdatabase\tother")


@pytest.mark.parametrize("engine", ["casanovo", "instanovo"])  # cascadia: M10-09
def test_other_engines_not_implemented(runner: CliRunner, tmp_path: Path, engine: str) -> None:
    src = _generic(tmp_path, "x\n")
    out = tmp_path / "peptides.tsv"
    result = runner.invoke(cli, ["convert", "--from", engine, str(src), "--out", str(out)])
    assert result.exit_code == ExitCode.INTERNAL
    assert "not implemented" in result.output.lower()
    assert not out.exists()
