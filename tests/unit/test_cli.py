"""Command-line interface skeleton (GIPPA_SPEC.md §5.1)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

SUBCOMMANDS = ["setup", "preflight", "run", "makedb", "convert", "compare", "report"]


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def genome(tmp_path: Path) -> Path:
    path = tmp_path / "g.fa"
    path.write_text(">chr1 [topology=circular]\n" + "ACGT" * 50 + "\n")
    return path


def test_main_entry_point_exists() -> None:
    from gippa.cli import main

    assert callable(main)


def test_version_option_prints_package_version(runner: CliRunner) -> None:
    import gippa
    from gippa.cli import cli

    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert gippa.__version__ in result.output


def test_help_lists_public_subcommands(runner: CliRunner) -> None:
    from gippa.cli import cli

    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    for name in SUBCOMMANDS:
        assert name in result.output


def test_help_hides_hidden_commands(runner: CliRunner) -> None:
    from gippa.cli import cli

    result = runner.invoke(cli, ["--help"])
    assert "_step" not in result.output
    commands_section = result.output.split("Commands:", 1)[1]
    assert not any(line.split()[:1] == ["version"] for line in commands_section.splitlines())


def test_hidden_step_command_is_registered() -> None:
    from gippa.cli import cli

    assert "_step" in cli.commands
    assert cli.commands["_step"].hidden


def test_version_json_has_gippa_and_python(runner: CliRunner) -> None:
    import gippa
    from gippa.cli import cli

    result = runner.invoke(cli, ["version", "--json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert data["gippa"] == gippa.__version__
    assert "python" in data


def test_exit_codes() -> None:
    from gippa.exitcodes import ExitCode

    assert ExitCode.SUCCESS == 0
    assert ExitCode.INVALID_INPUT == 1
    assert ExitCode.ENVIRONMENT == 2
    assert ExitCode.WARNINGS_STRICT == 3
    assert ExitCode.STEP_FAILED == 4
    assert ExitCode.INTERNAL == 5


@pytest.mark.parametrize("name", SUBCOMMANDS)
def test_subcommands_accept_shared_options(runner: CliRunner, name: str) -> None:
    from gippa.cli import cli

    result = runner.invoke(cli, [name, "--help"])
    assert result.exit_code == 0
    for opt in ("--config", "--set", "--log-level"):
        assert opt in result.output


@pytest.fixture
def proteins(tmp_path: Path) -> Path:
    path = tmp_path / "p.faa"
    path.write_text(">P1\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA\n")
    return path


@pytest.fixture(autouse=True)
def _no_db_manifest(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # preflight (M1-12): a missing db_manifest.json is info, so exit codes come from the inputs
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "db_empty"))
    monkeypatch.setenv("GIPPA_ROOT", str(tmp_path / "db_empty"))  # never the real $GIPPA_ROOT


def test_preflight_is_implemented(runner: CliRunner, genome: Path, proteins: Path) -> None:
    from gippa.cli import cli
    from gippa.exitcodes import ExitCode

    argv = ["preflight", "--profile", "none", "--genome", str(genome), "--proteins", str(proteins)]
    result = runner.invoke(cli, argv)
    # a 200 bp genome only warns (total length, contiguity); warnings exit 0 without --strict
    assert result.exit_code == ExitCode.SUCCESS, result.output
    assert "not implemented" not in result.output.lower()


def test_unknown_set_key_exits_invalid_input(runner: CliRunner, genome: Path) -> None:
    from gippa.cli import cli
    from gippa.exitcodes import ExitCode

    result = runner.invoke(
        cli, ["preflight", "--genome", str(genome), "--set", "stats.bogus_key=1"]
    )
    assert result.exit_code == ExitCode.INVALID_INPUT
    assert "bogus_key" in result.output


def test_valid_set_is_accepted_repeatedly(runner: CliRunner, genome: Path, proteins: Path) -> None:
    from gippa.cli import cli
    from gippa.exitcodes import ExitCode

    args = ["preflight", "--profile", "none", "--genome", str(genome), "--proteins", str(proteins)]
    result = runner.invoke(cli, [*args, "--set", "seed=1", "--set", "stats.permutations=5"])
    # config resolves; preflight (M1-12) then succeeds with warnings only
    assert result.exit_code == ExitCode.SUCCESS, result.output


def test_config_file_with_unknown_key_exits_invalid_input(
    runner: CliRunner, genome: Path, tmp_path: Path
) -> None:
    from gippa.cli import cli
    from gippa.exitcodes import ExitCode

    cfg = tmp_path / "c.yaml"
    cfg.write_text("not_a_key: 1\n")
    result = runner.invoke(cli, ["preflight", "--genome", str(genome), "--config", str(cfg)])
    assert result.exit_code == ExitCode.INVALID_INPUT


def test_run_rejects_bad_genome_relation(runner: CliRunner, genome: Path, tmp_path: Path) -> None:
    from gippa.cli import cli

    result = runner.invoke(
        cli,
        [
            "run",
            "--genome",
            str(genome),
            "--outdir",
            str(tmp_path / "o"),
            "--genome-relation",
            "bogus",
        ],
    )
    assert result.exit_code == 2  # click usage error


def test_run_rejects_bad_abundance_missing(runner: CliRunner, genome: Path, tmp_path: Path) -> None:
    from gippa.cli import cli

    result = runner.invoke(
        cli,
        [
            "run",
            "--genome",
            str(genome),
            "--outdir",
            str(tmp_path / "o"),
            "--abundance-missing",
            "knn",
        ],
    )
    assert result.exit_code == 2


def test_run_requires_outdir(runner: CliRunner, genome: Path) -> None:
    from gippa.cli import cli

    result = runner.invoke(cli, ["run", "--genome", str(genome)])
    assert result.exit_code == 2


def test_run_lists_spec_options(runner: CliRunner) -> None:
    from gippa.cli import cli

    result = runner.invoke(cli, ["run", "--help"])
    assert result.exit_code == 0
    for opt in (
        "--genome",
        "--peptides",
        "--proteins",
        "--abundance",
        "--annotation",
        "--reference-peptides",
        "--replicon-depth",
        "--outdir",
        "--db-dir",
        "--threads",
        "--memory",
        "--modules",
        "--genome-relation",
        "--abundance-missing",
        "--abundance-type",
        "--abundance-scale",
        "--promoters",
        "--terminators",
        "--add-bactermfinder",
        "--publication",
        "--seed",
        "--parquet",
        "--sensitivity",
        "--strict",
        "--profile",
        "--resume",
        "--executor",
    ):
        assert opt in result.output, opt


def test_run_valid_options_reach_the_runner(
    runner: CliRunner, genome: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # M2-07 wires `gippa run`: valid options get past option parsing and reach preflight, which
    # (run-wiring decisions 1 and 3d) runs before the `nextflow` check and stops a genome-only run
    # (no peptides or proteins) with exit 1 — not a usage error, not "not implemented".
    from gippa.cli import cli
    from gippa.exitcodes import ExitCode

    monkeypatch.setenv("PATH", str(tmp_path / "no_nextflow_here"))

    result = runner.invoke(
        cli,
        [
            "run",
            "--genome",
            str(genome),
            "--outdir",
            str(tmp_path / "o"),
            "--genome-relation",
            "same",
            "--seed",
            "3",
            "--publication",
        ],
    )
    assert result.exit_code == ExitCode.INVALID_INPUT, result.output
    assert "not implemented" not in result.output.lower()
