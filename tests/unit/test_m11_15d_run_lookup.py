"""M11-15d: how ``gippa run`` finds Nextflow and checks Java (spec item 5).

docs/specs/M11-15-portable-clusters.md item 5 "Lookup in ``gippa run``/``demo --full``":
``nextflow`` on ``PATH``, else ``<db_dir>/bin/nextflow`` (installed by ``gippa setup``); Java
missing or older than 17 is an error (exit 2) before Nextflow starts. Examples 5.1 (a later
``gippa run`` uses the launcher setup installed), 5.3, 5.4, 5.7; AC "run lookup-order test".

Pinned interface: ``gippa.runner.prereqs.find_nextflow(environ, db_dir) -> Path | None``
(``environ["PATH"]`` first; then ``<db_dir>/bin/nextflow`` when it is an executable file).
``gippa run`` launches the path found (``argv[0]``), so a launcher only in the DB directory runs.
``gippa demo --full`` is 15b's command (not on this branch); its Java/lookup tests come with it.

Nothing real runs: ``nextflow`` and ``java`` are fake scripts, the DB directory and
``$GIPPA_ROOT`` are fakes in ``tmp_path`` (as in ``test_cli_run_nextflow.py``).
"""

from __future__ import annotations

import os
import random
import tomllib
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.fake_env import fake_gippa_root, write_setup_manifest

from gippa.exitcodes import ExitCode

pytestmark = pytest.mark.real_prereqs  # these tests build their own PATH

REPO = Path(__file__).resolve().parents[2]
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
JAVA_HEAD = "Java 17 or newer is needed for Nextflow"
ADVICE = f"micromamba create -n {PKG}-java -c conda-forge openjdk=21"

FAKE_NF = """#!/bin/sh
printf '%s\\n' "$0" > "$FAKE_NF_ARGS.argv0"
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
exit 0
"""


def java_script(stderr: str, code: int = 0) -> str:
    lines = "".join(f"printf '%s\\n' '{ln}' >&2\n" for ln in stderr.splitlines())
    return f"#!/bin/sh\n{lines}exit {code}\n"


JAVA_21 = java_script('openjdk version "21.0.1" 2023-10-17\nOpenJDK Runtime Environment')


def make_exe(path: Path, text: str, mode: int = 0o755) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    path.chmod(mode)
    return path


def host_path_without(*names: str) -> list[str]:
    return [
        d
        for d in os.environ.get("PATH", "").split(os.pathsep)
        if d and not any((Path(d) / n).exists() for n in names)
    ]


def _prereqs() -> Any:
    try:
        from gippa.runner import prereqs
    except ImportError as exc:
        pytest.fail(f"gippa.runner.prereqs is missing (M11-15d lookup): {exc}")
    return prereqs


# ---------------------------------------------------------------- find_nextflow


def test_path_wins_over_the_db_dir(tmp_path: Path) -> None:
    on_path = make_exe(tmp_path / "bin" / "nextflow", FAKE_NF)
    make_exe(tmp_path / "db" / "bin" / "nextflow", FAKE_NF)
    got = _prereqs().find_nextflow({"PATH": str(tmp_path / "bin")}, tmp_path / "db")
    assert got is not None and Path(got) == on_path


def test_db_dir_launcher_when_none_on_path(tmp_path: Path) -> None:
    (tmp_path / "bin").mkdir()
    in_db = make_exe(tmp_path / "db" / "bin" / "nextflow", FAKE_NF)
    got = _prereqs().find_nextflow({"PATH": str(tmp_path / "bin")}, tmp_path / "db")
    assert got is not None and Path(got).absolute() == in_db.absolute()


@pytest.mark.parametrize("case", ["no_db_dir", "empty_db", "not_executable", "a_directory"])
def test_nothing_found(tmp_path: Path, case: str) -> None:
    (tmp_path / "bin").mkdir()
    db: Path | None = tmp_path / "db"
    if case == "no_db_dir":
        db = None
    elif case == "empty_db":
        (tmp_path / "db").mkdir()
    elif case == "not_executable":
        make_exe(tmp_path / "db" / "bin" / "nextflow", FAKE_NF, mode=0o644)
    else:
        (tmp_path / "db" / "bin" / "nextflow").mkdir(parents=True)
    assert _prereqs().find_nextflow({"PATH": str(tmp_path / "bin")}, db) is None


# ---------------------------------------------------------------- gippa run (CLI)


@pytest.fixture(autouse=True)
def _env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("GIPPA_ROOT", str(fake_gippa_root(tmp_path / "gippa_root")))
    for var in ("JAVA_HOME", "JAVA_CMD", "GIPPA_CONTAINER_RUNTIME", "NXF_APPTAINER_CACHEDIR"):
        monkeypatch.delenv(var, raising=False)


@pytest.fixture
def inputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    rng = random.Random(3)
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    genome = tmp_path / "in" / "genome.fa"
    genome.parent.mkdir()
    genome.write_text(
        ">chr1 [topology=circular]\n"
        + "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60))
        + "\n"
    )
    aa = "ACDEFGHIKLMNPQRSTVWY"
    proteins = tmp_path / "in" / "proteins.faa"
    proteins.write_text(
        "".join(f">prot{i}\nM{''.join(rng.choice(aa) for _ in range(150))}\n" for i in range(3))
    )
    db = tmp_path / "db"
    db.mkdir()
    write_setup_manifest(db)
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    return {"genome": genome, "proteins": proteins, "db": db}


def run_cli(inp: dict[str, Path], out: Path) -> tuple[int, str]:
    from gippa.cli import cli

    argv = ["run", "--genome", str(inp["genome"]), "--proteins", str(inp["proteins"])]
    res = CliRunner().invoke(cli, [*argv, "--outdir", str(out)])
    text = res.output + (res.stderr if res.stderr_bytes is not None else "")
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res.exit_code, text


def set_path(monkeypatch: pytest.MonkeyPatch, fakebin: Path) -> None:
    rest = host_path_without("nextflow")
    monkeypatch.setenv("PATH", os.pathsep.join([str(fakebin), *rest]))


def test_lookup_order_in_gippa_run(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Example 5.1's later run: only ``<db_dir>/bin/nextflow`` → it is launched; with a
    ``nextflow`` on ``PATH`` too, the ``PATH`` one is launched."""
    fakebin = tmp_path / "fakebin"
    make_exe(fakebin / "java", JAVA_21)
    in_db = make_exe(inputs["db"] / "bin" / "nextflow", FAKE_NF)
    set_path(monkeypatch, fakebin)
    args = tmp_path / "db_nf_args.txt"
    monkeypatch.setenv("FAKE_NF_ARGS", str(args))
    code, text = run_cli(inputs, tmp_path / "OUT1")
    assert code == ExitCode.SUCCESS, text
    assert args.is_file(), "the launcher in <db_dir>/bin was not run"
    assert Path(args.with_name(args.name + ".argv0").read_text().strip()).resolve() == (
        in_db.resolve()
    )
    assert args.read_text().splitlines()[0] == "run"

    on_path = make_exe(fakebin / "nextflow", FAKE_NF)
    args2 = tmp_path / "path_nf_args.txt"
    monkeypatch.setenv("FAKE_NF_ARGS", str(args2))
    code, text = run_cli(inputs, tmp_path / "OUT2")
    assert code == ExitCode.SUCCESS, text
    assert Path(args2.with_name(args2.name + ".argv0").read_text().strip()).resolve() == (
        on_path.resolve()
    )


def test_no_nextflow_anywhere_names_both_places(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fakebin = tmp_path / "fakebin"
    make_exe(fakebin / "java", JAVA_21)
    set_path(monkeypatch, fakebin)
    code, text = run_cli(inputs, tmp_path / "OUT")
    assert code == ExitCode.ENVIRONMENT, text
    assert str(inputs["db"] / "bin" / "nextflow") in text, text
    assert f"{PKG} setup" in text


@pytest.mark.parametrize(
    ("java", "found"),
    [
        (java_script('openjdk version "11.0.2" 2019-01-15\nOpenJDK Runtime Environment'), "11"),
        (java_script('java version "1.8.0_292"\nJava(TM) SE Runtime Environment'), "8"),
        (java_script("Unable to locate a Java Runtime.", 1), None),  # macOS stub, example 5.7
        (java_script("garbage"), None),
    ],
    ids=["5.3-java11", "5.4-java8", "5.7-macos-stub", "garbage"],
)
def test_examples_5_3_5_4_5_7_run_refuses_old_or_missing_java(
    inputs: dict[str, Path],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    java: str,
    found: str | None,
) -> None:
    fakebin = tmp_path / "fakebin"
    exe = make_exe(fakebin / "java", java)
    make_exe(fakebin / "nextflow", FAKE_NF)
    set_path(monkeypatch, fakebin)
    args = tmp_path / "nf_args.txt"
    monkeypatch.setenv("FAKE_NF_ARGS", str(args))
    code, text = run_cli(inputs, tmp_path / "OUT")
    assert code == ExitCode.ENVIRONMENT, text
    assert JAVA_HEAD in text and ADVICE in text, text
    if found is None:
        assert "not found" in text
    else:
        assert f"(found {found} at {exe})" in text
    assert not args.exists(), "Nextflow must not start without Java 17+"


def test_java_home_is_the_java_checked_by_run(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``JAVA_HOME`` names the Java Nextflow runs: an old one there fails even with Java 21 on
    ``PATH``."""
    fakebin = tmp_path / "fakebin"
    make_exe(fakebin / "java", JAVA_21)
    make_exe(fakebin / "nextflow", FAKE_NF)
    old = make_exe(
        tmp_path / "jdk11" / "bin" / "java",
        java_script('openjdk version "11.0.2" 2019-01-15'),
    )
    monkeypatch.setenv("JAVA_HOME", str(tmp_path / "jdk11"))
    set_path(monkeypatch, fakebin)
    monkeypatch.setenv("FAKE_NF_ARGS", str(tmp_path / "nf_args.txt"))
    code, text = run_cli(inputs, tmp_path / "OUT")
    assert code == ExitCode.ENVIRONMENT, text
    assert f"(found 11 at {old})" in text
