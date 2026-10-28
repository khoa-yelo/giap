"""``gippa convert --from cascadia`` through the CLI (M10-09; GIPPA_SPEC.md §5.1, §5.4, §6.6).

Spec: docs/specs/M10-09-cascadia-convert.md, decision D7, example 12, criteria 1 and 3.
Fixtures are hand-made 0.0.7-shaped ``.ssl`` text written into ``tmp_path``.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest
from click.testing import CliRunner, Result

from gippa.cli import cli
from gippa.exitcodes import ExitCode

W1 = "file\tscan\tcharge\tsequence\tscore-type\tscore\tretention-time\tstart-time\tend-time"
EX1_LINES = [
    W1,
    "/x/mzml/WT_1.mzML\t0\t2\tPEPTLDEK\tUNKNOWN\t0.95123\t4.5\t4.4\t4.6",
    "/x/mzml/WT_1.mzML\t3\t3\tM[+15.99]PEPTLDEKR\tUNKNOWN\t0.88120001\t4.5\t4.4\t4.6",
]
EX1_OUT = (
    "peptide\tmodified_peptide\tscore\trun\tcharge\tspectrum_id\tsource\n"
    "MPEPTLDEKR\tM[+15.99]PEPTLDEKR\t0.8812\tWT_1\t3\tcascadia:3\tdenovo\n"
    "PEPTLDEK\tPEPTLDEK\t0.95123\tWT_1\t2\tcascadia:0\tdenovo\n"
)


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


def _ssl(tmp_path: Path, name: str = "WT_1.ssl", lines: list[str] = EX1_LINES) -> Path:
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def _convert(runner: CliRunner, src: Path, out: Path, *extra: str) -> Result:
    return runner.invoke(
        cli, ["convert", "--from", "cascadia", str(src), "--out", str(out), *extra]
    )


# ------------------------------------------------------------------ criterion 1


def test_criterion1_cascadia_no_longer_exits_5(runner: CliRunner, tmp_path: Path) -> None:
    out = tmp_path / "p.tsv"
    result = _convert(runner, _ssl(tmp_path), out)
    assert result.exit_code == ExitCode.SUCCESS, result.output
    assert out.read_text() == EX1_OUT


@pytest.mark.parametrize("engine", ["casanovo", "instanovo"])
def test_criterion1_other_de_novo_engines_still_exit_5(
    runner: CliRunner, tmp_path: Path, engine: str
) -> None:
    out = tmp_path / "p.tsv"
    argv = ["convert", "--from", engine, str(_ssl(tmp_path)), "--out", str(out)]
    result = runner.invoke(cli, argv)
    assert result.exit_code == ExitCode.INTERNAL
    assert "not implemented" in result.output.lower()
    assert not out.exists()


# ------------------------------------------------------------------ example 12 / D7 flags


def test_source_denovo_is_accepted(runner: CliRunner, tmp_path: Path) -> None:
    out = tmp_path / "p.tsv"
    result = _convert(runner, _ssl(tmp_path), out, "--source", "denovo")
    assert result.exit_code == ExitCode.SUCCESS, result.output
    assert out.read_text() == EX1_OUT


def test_source_database_is_rejected(runner: CliRunner, tmp_path: Path) -> None:
    out = tmp_path / "p.tsv"
    result = _convert(runner, _ssl(tmp_path), out, "--source", "database")
    assert result.exit_code == ExitCode.INVALID_INPUT, result.output
    assert "--source database" in result.stderr
    assert "cascadia" in result.stderr
    assert not out.exists()


@pytest.mark.parametrize("db", ["same", "other"])
def test_search_db_is_rejected(runner: CliRunner, tmp_path: Path, db: str) -> None:
    out = tmp_path / "p.tsv"
    result = _convert(runner, _ssl(tmp_path), out, "--search-db", db)
    assert result.exit_code == ExitCode.INVALID_INPUT, result.output
    assert "search_db_invalid" in result.stderr
    assert not out.exists()


def test_abundance_out_is_rejected(runner: CliRunner, tmp_path: Path) -> None:
    out, ab = tmp_path / "p.tsv", tmp_path / "a.tsv"
    result = _convert(runner, _ssl(tmp_path), out, "--abundance-out", str(ab))
    assert result.exit_code == ExitCode.INVALID_INPUT, result.output
    assert "no_abundance" in result.stderr
    assert not out.exists()
    assert not ab.exists()


def test_gzip_crlf_bom_input_gives_the_same_bytes(runner: CliRunner, tmp_path: Path) -> None:
    src = tmp_path / "WT_1.ssl.gz"
    text = "﻿" + "\r\n".join(EX1_LINES) + "\r\n"
    src.write_bytes(gzip.compress(text.encode("utf-8"), mtime=0))
    out = tmp_path / "p.tsv"
    result = _convert(runner, src, out)
    assert result.exit_code == ExitCode.SUCCESS, result.output
    assert out.read_bytes() == EX1_OUT.encode()


def test_stripped_info_is_reported(runner: CliRunner, tmp_path: Path) -> None:
    result = _convert(runner, _ssl(tmp_path), tmp_path / "p.tsv")
    assert result.exit_code == ExitCode.SUCCESS, result.output
    assert "info: stripped" in result.stderr


# ------------------------------------------------------------------ error cases exit 1


@pytest.mark.parametrize(
    ("bad_line", "code"),
    [
        ("/x/mzml/WT_1.mzML\t5\t2\tSAMPLLER\tUNKNOWN\tnan\t4.5\t4.4\t4.6", "bad_type:score"),
        ("/x/mzml/WT_1.mzML\t5\t2\tSAMPLLER\tUNKNOWN\t1.2\t4.5\t4.4\t4.6", "out_of_range:score"),
        (
            "/x/mzml/WT_1.mzML\t5\t2\tSAMPLLER\tPERCOLATOR QVALUE\t0.9\t4.5\t4.4\t4.6",
            "unsupported_score_type",
        ),
        ("/x/mzml/WT_1.mzML\t-5\t2\tSAMPLLER\tUNKNOWN\t0.9\t4.5\t4.4\t4.6", "bad_type:scan"),
        ("/x/mzml/WT_1.mzML\t0\t2\tSAMPLLER\tUNKNOWN\t0.9\t4.5\t4.4\t4.6", "duplicate_spectrum_id"),
        ("/x/mzml/WT_1.mzML\t5\t0\tSAMPLLER\tUNKNOWN\t0.9\t4.5\t4.4\t4.6", "out_of_range:charge"),
    ],
)
def test_error_exits_1_and_writes_nothing(
    runner: CliRunner, tmp_path: Path, bad_line: str, code: str
) -> None:
    out = tmp_path / "p.tsv"
    result = _convert(runner, _ssl(tmp_path, lines=[*EX1_LINES, bad_line]), out)
    assert result.exit_code == ExitCode.INVALID_INPUT, result.output
    assert code in result.stderr
    assert not out.exists()


def test_header_only_alone_exits_1(runner: CliRunner, tmp_path: Path) -> None:
    out = tmp_path / "p.tsv"
    result = _convert(runner, _ssl(tmp_path, lines=[W1]), out)
    assert result.exit_code == ExitCode.INVALID_INPUT, result.output
    assert "empty_table" in result.stderr
    assert not out.exists()


def test_same_run_twice_exits_1(runner: CliRunner, tmp_path: Path) -> None:
    src = _ssl(tmp_path)
    out = tmp_path / "p.tsv"
    argv = ["convert", "--from", "cascadia", str(src), str(src), "--out", str(out)]
    result = runner.invoke(cli, argv)
    assert result.exit_code == ExitCode.INVALID_INPUT, result.output
    assert "run_in_several_inputs" in result.stderr
    assert not out.exists()


# ------------------------------------------------------------------ criterion 3 determinism


@pytest.mark.parametrize("suffix", [".tsv", ".tsv.gz"])
def test_two_runs_are_byte_identical(runner: CliRunner, tmp_path: Path, suffix: str) -> None:
    src = _ssl(tmp_path)
    blobs = []
    for i in range(2):
        out = tmp_path / f"o{i}{suffix}"
        result = _convert(runner, src, out)
        assert result.exit_code == ExitCode.SUCCESS, result.output
        blobs.append(out.read_bytes())
    assert blobs[0] == blobs[1]
    text = gzip.decompress(blobs[0]) if suffix.endswith(".gz") else blobs[0]
    assert text == EX1_OUT.encode()


def test_argument_order_gives_the_same_bytes(runner: CliRunner, tmp_path: Path) -> None:
    wt1 = _ssl(tmp_path)
    wt2 = _ssl(
        tmp_path,
        "WT_2.ssl",
        [W1, "/x/mzml/WT_2.mzML\t0\t2\tPEPTLDEK\tUNKNOWN\t0.9\t4.5\t4.4\t4.6"],
    )
    blobs = []
    for i, order in enumerate(((wt1, wt2), (wt2, wt1))):
        out = tmp_path / f"o{i}.tsv"
        argv = ["convert", "--from", "cascadia", *map(str, order), "--out", str(out)]
        result = runner.invoke(cli, argv)
        assert result.exit_code == ExitCode.SUCCESS, result.output
        blobs.append(out.read_bytes())
    assert blobs[0] == blobs[1]
    assert blobs[0].decode().splitlines()[-1].split("\t")[3] == "WT_2"
