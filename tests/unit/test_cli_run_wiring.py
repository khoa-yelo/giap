"""``gippa run`` wires preflight output and the database directory into Nextflow (M2-07 follow-up).

docs/specs/M2-07-nextflow-tool-processes.md, "Follow-up PR (stacked, branch m2/run-wiring)" and code
review item 1 of PR #65: ``workflow/main.nf`` takes the ANNOTATE branch only when
``params.inputs.genome`` and ``params.inputs.replicons`` are both set; ``gippa run --genome G`` must
therefore run preflight (GIPPA_SPEC.md §5.1: errors stop the run, exit 1) and pass the normalised
``replicons.tsv`` (``id id type c|l -``, preflight format) plus ``params.db_dir`` (``--db-dir``,
else ``$GIPPA_DB_DIR``) to Nextflow.

A fake ``nextflow`` executable on ``PATH`` records its argv; nothing real is launched. Every test
provides a database directory whose ``db_manifest.json`` lists every requirement, so the run-mode
environment check of preflight passes whichever way the implementation probes it, and a protein
list (preflight errors without any evidence: ``inputs.evidence_present``).
"""

from __future__ import annotations

import json
import os
import random
import stat
from pathlib import Path

import pytest
from click.testing import CliRunner
from helpers.fake_env import fake_gippa_root, write_setup_manifest

from gippa.exitcodes import ExitCode
from gippa.runner.nextflow import launch_command

FAKE_NF = """#!/bin/sh
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
exit 0
"""
ALL_MODULES = ["core", "mobile", "bgc", "function"]


@pytest.fixture(autouse=True)
def _fake_gippa_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Every pinned image and prebuilt env present (run-wiring decision 4a); never the real root."""
    root = fake_gippa_root(tmp_path / "gippa_root")
    monkeypatch.setenv("GIPPA_ROOT", str(root))
    return root


def _write_db(path: Path) -> Path:
    path.mkdir(parents=True)
    write_setup_manifest(path)  # real `gippa setup` format (decision 4c)
    return path


def _fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            fh.writelines(seq[i : i + 60] + "\n" for i in range(0, len(seq), 60))
    return path


def _seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture
def genome(tmp_path: Path) -> Path:
    rng = random.Random(3)
    return _fasta(tmp_path / "in" / "genome.fa", [("chr1 [topology=circular]", _seq(rng, 6000))])


@pytest.fixture
def proteins(tmp_path: Path) -> Path:
    rng = random.Random(7)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    recs = [(f"prot{i}", "M" + "".join(rng.choice(aa) for _ in range(150))) for i in range(3)]
    return _fasta(tmp_path / "in" / "proteins.faa", recs)


@pytest.fixture
def db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    path = _write_db(tmp_path / "db_env")
    monkeypatch.setenv("GIPPA_DB_DIR", str(path))
    return path


@pytest.fixture
def fake_nf(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(FAKE_NF)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    args = tmp_path / "nf_args.txt"
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_ARGS", str(args))
    return args


def _run(argv: list[str], proteins: Path) -> tuple[int, str]:
    from gippa.cli import cli

    result = CliRunner().invoke(cli, ["run", *argv, "--proteins", str(proteins)])
    out = result.output + (result.stderr if result.stderr_bytes is not None else "")
    if result.exception is not None and not isinstance(result.exception, SystemExit):
        raise result.exception
    return result.exit_code, out


def _params(out: Path) -> dict[str, object]:
    path = out / "work" / "params.json"
    assert path.is_file(), "gippa run wrote no params.json"
    data = json.loads(path.read_text())
    assert isinstance(data, dict)
    return data


def _inputs(out: Path) -> dict[str, str]:
    inputs = _params(out)["inputs"]
    assert isinstance(inputs, dict)
    return inputs


def _fasta_ids(path: Path) -> list[str]:
    return [ln[1:].split()[0] for ln in path.read_text().splitlines() if ln.startswith(">")]


def _replicon_rows(path: Path) -> list[list[str]]:
    return [ln.split("\t") for ln in path.read_text().splitlines() if ln.strip()]


# ---------------------------------------------------------------- inputs.replicons (review item 1)


def test_run_passes_replicons_so_main_nf_takes_annotate(
    proteins: Path, genome: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    code, text = _run(["--genome", str(genome), "--outdir", str(out)], proteins)
    assert code == ExitCode.SUCCESS, text
    inputs = _inputs(out)
    # main.nf: `if (inputs.genome && inputs.replicons) ANNOTATE(...) else smoke`
    assert inputs.get("genome") and inputs.get("replicons"), f"smoke branch: {inputs}"
    replicons = Path(inputs["replicons"])
    assert replicons.is_absolute() and replicons.is_file(), replicons
    rows = _replicon_rows(replicons)
    assert rows == [["chr1", "chr1", "chromosome", "c", "-"]], rows


def test_run_genome_and_replicon_ids_agree(
    proteins: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    """Bakta reads ``--replicons`` against the genome it gets: the IDs must be the same ones
    (a header preflight sanitises, ``chr|1`` → ``chr_1``, must not split them)."""
    rng = random.Random(5)
    genome = _fasta(
        tmp_path / "in" / "g.fa",
        [
            ("chr|1 [topology=circular]", _seq(rng, 6000)),
            ("p|2 [topology=linear]", _seq(rng, 3000)),
        ],
    )
    out = tmp_path / "OUT"
    code, text = _run(["--genome", str(genome), "--outdir", str(out)], proteins)
    assert code == ExitCode.SUCCESS, text
    inputs = _inputs(out)
    g = Path(inputs["genome"])
    assert g.is_absolute() and g.is_file(), g
    rows = _replicon_rows(Path(inputs["replicons"]))
    assert all(len(r) == 5 and r[3] in {"c", "l"} for r in rows), rows
    assert [r[0] for r in rows] == _fasta_ids(g)


def test_run_writes_the_preflight_report(
    proteins: Path, genome: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    code, text = _run(["--genome", str(genome), "--outdir", str(out)], proteins)
    assert code == ExitCode.SUCCESS, text
    assert (out / "preflight" / "preflight.json").is_file()


def test_run_preflight_error_exits_1_and_never_launches_nextflow(
    proteins: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    bad = tmp_path / "in" / "bad.fa"
    bad.parent.mkdir(exist_ok=True)
    bad.write_text("this is not a FASTA file\n")
    out = tmp_path / "OUT"
    code, text = _run(["--genome", str(bad), "--outdir", str(out)], proteins)
    assert code == ExitCode.INVALID_INPUT, text
    assert not fake_nf.exists(), "nextflow launched after a preflight error"


def test_run_launch_argv_points_at_the_written_params(
    proteins: Path, genome: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    code, text = _run(["--genome", str(genome), "--outdir", str(out)], proteins)
    assert code == ExitCode.SUCCESS, text
    args = fake_nf.read_text().splitlines()
    params = out.absolute() / "work" / "params.json"
    want = launch_command(Path(args[1]), params, "farmshare,apptainer", out.absolute(), False)
    assert ["nextflow", *args] == want
    assert _inputs(out).get("replicons"), "the launched params file has no inputs.replicons"


# ---------------------------------------------------------------- params.db_dir


def test_run_db_dir_flag_reaches_params(
    proteins: Path, genome: Path, fake_nf: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv("GIPPA_DB_DIR", raising=False)
    flag = _write_db(tmp_path / "db_flag")
    out = tmp_path / "OUT"
    code, text = _run(
        ["--genome", str(genome), "--outdir", str(out), "--db-dir", str(flag)], proteins
    )
    assert code == ExitCode.SUCCESS, text
    assert _params(out).get("db_dir") == str(flag.absolute())


def test_run_db_dir_defaults_to_gippa_db_dir(
    proteins: Path, genome: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    code, text = _run(["--genome", str(genome), "--outdir", str(out)], proteins)
    assert code == ExitCode.SUCCESS, text
    assert _params(out).get("db_dir") == str(db.absolute())


def test_run_db_dir_flag_wins_over_the_environment(
    proteins: Path, genome: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    flag = _write_db(tmp_path / "db_flag")
    out = tmp_path / "OUT"
    code, text = _run(
        ["--genome", str(genome), "--outdir", str(out), "--db-dir", str(flag)], proteins
    )
    assert code == ExitCode.SUCCESS, text
    assert _params(out).get("db_dir") == str(flag.absolute())


def test_run_db_dir_relative_flag_is_made_absolute(
    proteins: Path,
    genome: Path,
    db: Path,
    fake_nf: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _write_db(tmp_path / "db_rel")
    monkeypatch.chdir(tmp_path)
    out = tmp_path / "OUT"
    code, text = _run(
        ["--genome", str(genome), "--outdir", str(out), "--db-dir", "db_rel"], proteins
    )
    assert code == ExitCode.SUCCESS, text
    assert _params(out).get("db_dir") == str((tmp_path / "db_rel").absolute())
