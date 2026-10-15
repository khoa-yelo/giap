"""Exit codes of the console entry point ``giap.cli.main`` (GIAP_SPEC.md §5.1).

§5.1 reserves exit code 2 for environment errors, so click usage errors (missing file, unknown
option, bad choice) must surface as 1 (INVALID_INPUT) through ``main``. Every subcommand also
accepts ``--version``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import giap
from giap.exitcodes import ExitCode

SUBCOMMANDS = ["setup", "preflight", "run", "makedb", "convert", "compare", "report"]


def _exit_code(argv: list[str]) -> int:
    from giap.cli import main

    with pytest.raises(SystemExit) as info:
        main(argv)
    code = info.value.code
    if code is None:
        return 0
    assert isinstance(code, int), f"non-integer exit code {code!r}"
    return code


@pytest.fixture
def genome(tmp_path: Path) -> Path:
    path = tmp_path / "g.fa"
    path.write_text(">chr1 [topology=circular]\n" + "ACGT" * 50 + "\n")
    return path


@pytest.fixture
def peptides(tmp_path: Path) -> Path:
    path = tmp_path / "peptides.tsv"
    path.write_text("peptide\tpsm_count\nPEPTIDEK\t1\n")
    return path


def test_missing_input_file_is_invalid_input(tmp_path: Path) -> None:
    missing = tmp_path / "missing.fa"
    assert _exit_code(["preflight", "--genome", str(missing)]) == ExitCode.INVALID_INPUT


def test_unknown_option_is_invalid_input() -> None:
    assert _exit_code(["run", "--bogus"]) == ExitCode.INVALID_INPUT


def test_bad_abundance_type_choice_is_invalid_input(genome: Path, tmp_path: Path) -> None:
    argv = ["run", "--genome", str(genome), "--outdir", str(tmp_path / "out")]
    assert _exit_code([*argv, "--abundance-type", "nope"]) == ExitCode.INVALID_INPUT


def test_bad_genome_relation_choice_is_invalid_input(genome: Path, tmp_path: Path) -> None:
    argv = ["run", "--genome", str(genome), "--outdir", str(tmp_path / "out")]
    assert _exit_code([*argv, "--genome-relation", "bogus"]) == ExitCode.INVALID_INPUT


@pytest.mark.parametrize("argv", [["--help"], ["preflight", "--help"], ["-h"]])
def test_help_exits_zero(argv: list[str], capsys: pytest.CaptureFixture[str]) -> None:
    assert _exit_code(argv) == 0
    assert "Usage" in capsys.readouterr().out


def test_no_arguments_is_neither_environment_nor_internal_error() -> None:
    code = _exit_code([])
    assert code in (0, ExitCode.INVALID_INPUT)
    assert code not in (2, ExitCode.INTERNAL)


@pytest.mark.parametrize("command", SUBCOMMANDS)
def test_version_option_on_every_subcommand(
    command: str, capsys: pytest.CaptureFixture[str]
) -> None:
    assert _exit_code([command, "--version"]) == 0
    assert giap.__version__ in capsys.readouterr().out


def test_modules_with_spaces_reach_the_runner(
    genome: Path, peptides: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # M2-07: `--modules "core, mobile"` parses and reaches the launcher (no nextflow on PATH → 2).
    monkeypatch.setenv("PATH", str(tmp_path / "no_nextflow_here"))
    argv = [
        "run",
        "--genome",
        str(genome),
        "--peptides",
        str(peptides),
        "--outdir",
        str(tmp_path / "out"),
        "--modules",
        "core, mobile",
    ]
    assert _exit_code(argv) == ExitCode.ENVIRONMENT


def test_preflight_is_not_a_stub(
    genome: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """M1-12: preflight runs its rules; warnings alone exit 0, ``--strict`` turns them into 3."""
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIAP_DB_DIR", str(tmp_path / "db_empty"))
    prot = tmp_path / "p.faa"
    prot.write_text(">P1\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA\n")
    argv = ["preflight", "--genome", str(genome), "--proteins", str(prot)]
    assert _exit_code(argv) == ExitCode.SUCCESS
    assert _exit_code([*argv, "--strict"]) == ExitCode.WARNINGS_STRICT
    assert _exit_code(["preflight", "--genome", str(genome)]) == ExitCode.INVALID_INPUT
