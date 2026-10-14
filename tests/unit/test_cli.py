"""Command-line interface skeleton (GIAP_SPEC.md §5.1)."""

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
    from giap.cli import main

    assert callable(main)


def test_version_option_prints_package_version(runner: CliRunner) -> None:
    from giap.cli import cli

    import giap

    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert giap.__version__ in result.output


def test_help_lists_public_subcommands(runner: CliRunner) -> None:
    from giap.cli import cli

    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    for name in SUBCOMMANDS:
        assert name in result.output


def test_help_hides_hidden_commands(runner: CliRunner) -> None:
    from giap.cli import cli

    result = runner.invoke(cli, ["--help"])
    assert "_step" not in result.output
    commands_section = result.output.split("Commands:", 1)[1]
    assert not any(line.split()[:1] == ["version"] for line in commands_section.splitlines())


def test_hidden_step_command_is_registered() -> None:
    from giap.cli import cli

    assert "_step" in cli.commands
    assert cli.commands["_step"].hidden


def test_version_json_has_giap_and_python(runner: CliRunner) -> None:
    from giap.cli import cli

    import giap

    result = runner.invoke(cli, ["version", "--json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert data["giap"] == giap.__version__
    assert "python" in data


def test_exit_codes() -> None:
    from giap.exitcodes import ExitCode

    assert ExitCode.SUCCESS == 0
    assert ExitCode.INVALID_INPUT == 1
    assert ExitCode.ENVIRONMENT == 2
    assert ExitCode.WARNINGS_STRICT == 3
    assert ExitCode.STEP_FAILED == 4
    assert ExitCode.INTERNAL == 5


@pytest.mark.parametrize("name", SUBCOMMANDS)
def test_subcommands_accept_shared_options(runner: CliRunner, name: str) -> None:
    from giap.cli import cli

    result = runner.invoke(cli, [name, "--help"])
    assert result.exit_code == 0
    for opt in ("--config", "--set", "--log-level"):
        assert opt in result.output


def test_preflight_not_implemented_exits_internal(runner: CliRunner, genome: Path) -> None:
    from giap.cli import cli
    from giap.exitcodes import ExitCode

    result = runner.invoke(cli, ["preflight", "--genome", str(genome)])
    assert result.exit_code == ExitCode.INTERNAL
    assert "not implemented" in result.output.lower()


def test_unknown_set_key_exits_invalid_input(runner: CliRunner, genome: Path) -> None:
    from giap.cli import cli
    from giap.exitcodes import ExitCode

    result = runner.invoke(
        cli, ["preflight", "--genome", str(genome), "--set", "stats.bogus_key=1"]
    )
    assert result.exit_code == ExitCode.INVALID_INPUT
    assert "bogus_key" in result.output


def test_valid_set_is_accepted_repeatedly(runner: CliRunner, genome: Path) -> None:
    from giap.cli import cli
    from giap.exitcodes import ExitCode

    result = runner.invoke(
        cli,
        ["preflight", "--genome", str(genome), "--set", "seed=1", "--set", "stats.permutations=5"],
    )
    # config resolves; command then reports it is not implemented yet (M0)
    assert result.exit_code == ExitCode.INTERNAL


def test_config_file_with_unknown_key_exits_invalid_input(
    runner: CliRunner, genome: Path, tmp_path: Path
) -> None:
    from giap.cli import cli
    from giap.exitcodes import ExitCode

    cfg = tmp_path / "c.yaml"
    cfg.write_text("not_a_key: 1\n")
    result = runner.invoke(cli, ["preflight", "--genome", str(genome), "--config", str(cfg)])
    assert result.exit_code == ExitCode.INVALID_INPUT


def test_run_rejects_bad_genome_relation(runner: CliRunner, genome: Path, tmp_path: Path) -> None:
    from giap.cli import cli

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
    from giap.cli import cli

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
    from giap.cli import cli

    result = runner.invoke(cli, ["run", "--genome", str(genome)])
    assert result.exit_code == 2


def test_run_lists_spec_options(runner: CliRunner) -> None:
    from giap.cli import cli

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


def test_run_valid_options_reach_not_implemented(
    runner: CliRunner, genome: Path, tmp_path: Path
) -> None:
    from giap.cli import cli
    from giap.exitcodes import ExitCode

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
    assert result.exit_code == ExitCode.INTERNAL
    assert "not implemented" in result.output.lower()
