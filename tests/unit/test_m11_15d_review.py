"""M11-15d review B-items and nits on PR #437 (the Nextflow launcher and the Java check).

- **B1.** The Java check finds Java the way the Nextflow launcher does (``nextflow`` lines
  303-320 of the bundled 25.10/26.10 launcher, quoted in :data:`LAUNCHER_JAVA_SNIPPET`):

  1. ``NXF_JAVA_HOME`` (non-empty) replaces ``JAVA_HOME`` and unsets ``JAVA_CMD``;
  2. ``JAVA_CMD`` when it is executable;
  3. ``$JAVA_HOME/bin/java`` when ``JAVA_HOME`` is a directory (whether or not that file exists);
  4. macOS ``/usr/libexec/java_home -v 17+`` (not tested here: Linux);
  5. ``java`` on ``PATH``.

  ``gippa.runner.prereqs.check_java(environ)`` and the exit code of ``gippa run`` both follow it.
  The oracle test runs the launcher's own shell lines over every combination of the four
  variables and compares the Java it picks with ``check_java``.
- **B2 (lead decision).** ``run_manifest.json`` records the launcher a run used under
  ``launch.nextflow = {"path": <absolute>, "source": "path" | "db_dir"}`` (``version`` optional).
  Why ``launch``: M11-15a adds the top-level ``launch`` object (how ``gippa run`` chose the
  profile: ``profile``, ``profile_source``, ``profile_reason``), volatile (``/launch``) and
  outside the provenance core. The launcher is the same kind of fact: a machine-specific launch
  choice that must not change ``provenance_digest``. ``OUT/work/provenance.json`` is the core
  only (no paths), so it is not the place. The tests read only ``launch.nextflow``, so they hold
  with and without 15a's profile keys.
- **B3.** ``gippa setup --check`` with no ``nextflow`` on ``PATH`` never plans the launcher.
- **Nits.** 1: a relative ``PATH`` entry still launches an absolute path; 2: an adopted launcher
  that lost its x bit gets it back; 3: an adopted launcher keeps its original ``retrieved``.

Nothing real runs: ``java`` and ``nextflow`` are ``/bin/sh`` scripts in ``tmp_path``; the network
is ``helpers.uniop_14u1.FakeNet`` behind sockets/DNS blocked (``no_network``).
"""

from __future__ import annotations

import dataclasses
import hashlib
import itertools
import json
import os
import random
import shutil
import stat
import subprocess
import sys
import tomllib
from collections.abc import Mapping
from pathlib import Path
from types import MappingProxyType
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.fake_env import fake_gippa_root, write_setup_manifest
from helpers.setup_opt_outs import NO_RESOURCE_FLAGS, opt_outs
from helpers.uniop_14u1 import FakeNet, Route, install_fake_net, isolate, manifest_items, setup_cli

from gippa.exitcodes import ExitCode

pytestmark = [pytest.mark.usefixtures("no_network"), pytest.mark.real_prereqs]

REPO = Path(__file__).resolve().parents[2]
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
JAVA_HEAD = "Java 17 or newer is needed for Nextflow"
JAVA_VARS = ("NXF_JAVA_HOME", "JAVA_HOME", "JAVA_CMD")

# The Java lookup of the Nextflow launcher, verbatim (nextflow 25.10/26.10 launcher script,
# Apache-2.0, lines 303-320 of local/micromamba/envs/giap-wf/bin/nextflow).
LAUNCHER_JAVA_SNIPPET = r"""
# use nextflow custom java home path
if [[ "$NXF_JAVA_HOME" ]]; then
  JAVA_HOME="$NXF_JAVA_HOME"
  unset JAVA_CMD
fi
# Determine the Java command to use to start the JVM.
if [ ! -x "$JAVA_CMD" ] ; then
    if [ -d "$JAVA_HOME" ] ; then
        if [ -x "$JAVA_HOME/jre/sh/java" ] ; then
            # IBM's JDK on AIX uses strange locations for the executables
            JAVA_CMD="$JAVA_HOME/jre/sh/java"
        else
            JAVA_CMD="$JAVA_HOME/bin/java"
        fi
    elif [ -x /usr/libexec/java_home ]; then
        JAVA_CMD="$(/usr/libexec/java_home -v 17+ 2>/dev/null)/bin/java" || JAVA_CMD=java
    else
        JAVA_CMD="$(command -v java)" || JAVA_CMD=java
    fi
fi
printf '%s' "$JAVA_CMD"
"""

FAKE_NF = """#!/bin/sh
case "$1" in
-v) echo "nextflow version 26.10.0.0"; exit 0;;
esac
printf '%s\\n' "$0" > "$FAKE_NF_ARGS.argv0"
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
exit 0
"""


def java_script(major: int) -> str:
    """A fake ``java`` whose ``-version`` banner (stderr) reports ``major``."""
    v = f"{major}.0.1" if major >= 9 else f"1.{major}.0_292"
    return (
        "#!/bin/sh\n"
        f"printf '%s\\n' 'openjdk version \"{v}\" 2023-10-17' >&2\n"
        f"printf '%s\\n' 'OpenJDK Runtime Environment (build {v}+12)' >&2\n"
        "exit 0\n"
    )


def make_exe(path: Path, text: str, mode: int = 0o755) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    path.chmod(mode)
    return path


def make_jdk(home: Path, major: int) -> Path:
    """``home/bin/java`` reporting ``major``; returns ``home``."""
    make_exe(home / "bin" / "java", java_script(major))
    return home


def host_path_without(*names: str) -> list[str]:
    return [
        d
        for d in os.environ.get("PATH", "").split(os.pathsep)
        if d and not any((Path(d) / n).exists() for n in names)
    ]


def _prereqs() -> Any:
    try:
        from gippa.runner import prereqs
    except ImportError as exc:  # pragma: no cover - the module exists on this branch
        pytest.fail(f"gippa.runner.prereqs is missing: {exc}")
    return prereqs


@pytest.fixture(autouse=True)
def _env(monkeypatch: pytest.MonkeyPatch) -> None:
    for var in (*JAVA_VARS, "GIPPA_CONTAINER_RUNTIME", "NXF_APPTAINER_CACHEDIR", "NXF_HOME"):
        monkeypatch.delenv(var, raising=False)


# ======================================== B1: check_java follows the launcher


@pytest.fixture
def jdks(tmp_path: Path) -> dict[str, Path]:
    """Fake JDKs with distinct majors, so the version found names the source that won."""
    return {
        "jdk25": make_jdk(tmp_path / "jdk25", 25),
        "jdk21": make_jdk(tmp_path / "jdk21", 21),
        "jdk17": make_jdk(tmp_path / "jdk17", 17),
        "jdk11": make_jdk(tmp_path / "jdk11", 11),
        "pathbin11": make_exe(tmp_path / "pathbin11" / "java", java_script(11)).parent,
        "pathbin21": make_exe(tmp_path / "pathbin21" / "java", java_script(21)).parent,
        "emptybin": (tmp_path / "emptybin").resolve(),
    }


def check(environ: Mapping[str, str]) -> tuple[int | None, str | None]:
    got = _prereqs().check_java(dict(environ))
    return got.major, got.path


def same_file(a: str | None, b: Path) -> bool:
    return a is not None and Path(a).resolve() == b.resolve()


def test_nxf_java_home_overrides_java_home(jdks: dict[str, Path]) -> None:
    """Review B1 row 1: ``NXF_JAVA_HOME``=21, ``JAVA_HOME``=11 -> Nextflow runs 21."""
    env = {"NXF_JAVA_HOME": str(jdks["jdk21"]), "JAVA_HOME": str(jdks["jdk11"]),
           "PATH": str(jdks["pathbin11"])}  # fmt: skip
    major, path = check(env)
    assert major == 21, (major, path)
    assert same_file(path, jdks["jdk21"] / "bin" / "java"), path


def test_nxf_java_home_ignores_java_cmd(jdks: dict[str, Path], tmp_path: Path) -> None:
    """The launcher unsets ``JAVA_CMD`` when ``NXF_JAVA_HOME`` is set."""
    cmd = make_exe(tmp_path / "cmd11" / "java", java_script(11))
    env = {"NXF_JAVA_HOME": str(jdks["jdk21"]), "JAVA_CMD": str(cmd),
           "PATH": str(jdks["pathbin11"])}  # fmt: skip
    major, path = check(env)
    assert major == 21, (major, path)
    assert same_file(path, jdks["jdk21"] / "bin" / "java"), path


def test_nxf_java_home_nonexistent_falls_through_to_path(jdks: dict[str, Path]) -> None:
    """``NXF_JAVA_HOME`` that is no directory: ``JAVA_HOME`` is replaced by it (so the real
    ``JAVA_HOME`` is ignored), not a directory, so ``java`` on ``PATH``."""
    env = {"NXF_JAVA_HOME": "/nonexistent/jdk", "JAVA_HOME": str(jdks["jdk11"]),
           "PATH": str(jdks["pathbin21"])}  # fmt: skip
    major, path = check(env)
    assert major == 21, (major, path)
    assert same_file(path, jdks["pathbin21"] / "java"), path


def test_empty_nxf_java_home_is_ignored(jdks: dict[str, Path], tmp_path: Path) -> None:
    cmd = make_exe(tmp_path / "cmd21" / "java", java_script(21))
    env = {"NXF_JAVA_HOME": "", "JAVA_CMD": str(cmd), "PATH": str(jdks["pathbin11"])}
    major, path = check(env)
    assert major == 21, (major, path)
    assert same_file(path, cmd), path


def test_java_cmd_wins_over_path(jdks: dict[str, Path], tmp_path: Path) -> None:
    """Review B1 row 2: ``JAVA_CMD``=21, ``PATH`` java 11 -> 21."""
    cmd = make_exe(tmp_path / "cmd21" / "java", java_script(21))
    major, path = check({"JAVA_CMD": str(cmd), "PATH": str(jdks["pathbin11"])})
    assert major == 21, (major, path)
    assert same_file(path, cmd), path


def test_java_cmd_wins_over_java_home(jdks: dict[str, Path], tmp_path: Path) -> None:
    cmd = make_exe(tmp_path / "cmd21" / "java", java_script(21))
    env = {"JAVA_CMD": str(cmd), "JAVA_HOME": str(jdks["jdk11"]), "PATH": str(jdks["pathbin11"])}
    major, path = check(env)
    assert major == 21, (major, path)
    assert same_file(path, cmd), path


@pytest.mark.parametrize("home", ["unset", "jdk17"])
def test_non_executable_java_cmd_falls_through(
    jdks: dict[str, Path], tmp_path: Path, home: str
) -> None:
    """A ``JAVA_CMD`` without the x bit is skipped: ``JAVA_HOME`` (a directory), else ``PATH``."""
    cmd = make_exe(tmp_path / "noexec" / "java", java_script(25), mode=0o644)
    env = {"JAVA_CMD": str(cmd), "PATH": str(jdks["pathbin21"])}
    if home == "jdk17":
        env["JAVA_HOME"] = str(jdks["jdk17"])
    major, path = check(env)
    want, exe = (17, jdks["jdk17"] / "bin" / "java") if home == "jdk17" else (
        21, jdks["pathbin21"] / "java")  # fmt: skip
    assert major == want, (major, path)
    assert same_file(path, exe), path


def test_missing_java_cmd_falls_through(jdks: dict[str, Path]) -> None:
    major, path = check({"JAVA_CMD": "/nonexistent/java", "PATH": str(jdks["pathbin21"])})
    assert major == 21, (major, path)


@pytest.mark.parametrize("kind", ["nonexistent", "a_file"])
def test_java_home_not_a_directory_falls_back_to_path(
    jdks: dict[str, Path], tmp_path: Path, kind: str
) -> None:
    """Review B1 row 3: ``JAVA_HOME=/nonexistent``, ``PATH`` java 21 -> 21, not "not found"."""
    home = "/nonexistent"
    if kind == "a_file":
        (tmp_path / "afile").write_text("not a JDK\n", encoding="utf-8")
        home = str(tmp_path / "afile")
    major, path = check({"JAVA_HOME": home, "PATH": str(jdks["pathbin21"])})
    assert major == 21, (major, path)
    assert same_file(path, jdks["pathbin21"] / "java"), path


def test_java_home_directory_without_java_does_not_fall_back(
    jdks: dict[str, Path], tmp_path: Path
) -> None:
    """A ``JAVA_HOME`` directory is used even without ``bin/java`` (the launcher then fails), so
    the check reports not found rather than the ``PATH`` java Nextflow would never run."""
    (tmp_path / "jdk_broken").mkdir()
    major, _path = check({"JAVA_HOME": str(tmp_path / "jdk_broken"),
                          "PATH": str(jdks["pathbin21"])})  # fmt: skip
    assert major is None


# ---------------------------------------------------------------- oracle: the launcher's own lines


_NXF = ("unset", "empty", "jdk25", "nonexistent", "dir_no_java")
_CMD = ("unset", "exec21", "noexec", "nonexistent")
_HOME = ("unset", "jdk17", "nonexistent", "a_file", "dir_no_java")
_PATHJ = ("none", "java11")


def _value(kind: str, tree: Mapping[str, Path], var: str) -> str | None:
    return {
        "unset": None,
        "empty": "",
        "jdk25": str(tree["jdk25"]),
        "jdk17": str(tree["jdk17"]),
        "nonexistent": f"/nonexistent/{var.lower()}",
        "dir_no_java": str(tree["dir_no_java"]),
        "a_file": str(tree["a_file"]),
        "exec21": str(tree["exec21"]),
        "noexec": str(tree["noexec"]),
    }[kind]


@pytest.fixture(scope="module")
def tree(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("jdk_tree")
    (root / "dir_no_java").mkdir()
    (root / "a_file").write_text("x\n", encoding="utf-8")
    (root / "nojava").mkdir()
    return {
        "jdk25": make_jdk(root / "jdk25", 25),
        "jdk17": make_jdk(root / "jdk17", 17),
        "exec21": make_exe(root / "cmd" / "java", java_script(21)),
        "noexec": make_exe(root / "noexec" / "java", java_script(8), mode=0o644),
        "path11": make_exe(root / "path11" / "java", java_script(11)).parent,
        "nojava": root / "nojava",
        "dir_no_java": root / "dir_no_java",
        "a_file": root / "a_file",
    }


_MAJOR_OF = {"jdk25": 25, "jdk17": 17, "exec21": 21, "path11": 11}


def _oracle(environ: Mapping[str, str]) -> str:
    bash = shutil.which("bash") or "/bin/bash"
    proc = subprocess.run(
        [bash, "--noprofile", "--norc", "-c", LAUNCHER_JAVA_SNIPPET],
        capture_output=True, text=True, env=dict(environ), check=True, timeout=10,
    )  # fmt: skip
    return proc.stdout


@pytest.mark.skipif(
    sys.platform == "darwin" or os.access("/usr/libexec/java_home", os.X_OK),
    reason="the macOS java_home branch (step 4) is not tested on Linux",
)
@pytest.mark.parametrize(
    ("nxf", "cmd", "home", "pathj"), list(itertools.product(_NXF, _CMD, _HOME, _PATHJ))
)
def test_check_java_agrees_with_the_launcher(
    tree: dict[str, Path], nxf: str, cmd: str, home: str, pathj: str
) -> None:
    """Every combination of ``NXF_JAVA_HOME``/``JAVA_CMD``/``JAVA_HOME``/``PATH``: the Java
    ``check_java`` runs is the ``JAVA_CMD`` the launcher's own lines pick (same major; same file
    when that file is a fake JDK)."""
    env = {"PATH": str(tree["path11"] if pathj == "java11" else tree["nojava"])}
    for var, kind in (("NXF_JAVA_HOME", nxf), ("JAVA_CMD", cmd), ("JAVA_HOME", home)):
        value = _value(kind, tree, var)
        if value is not None:
            env[var] = value
    picked = _oracle(env)
    exes = {"jdk25": tree["jdk25"] / "bin" / "java", "jdk17": tree["jdk17"] / "bin" / "java",
            "exec21": tree["exec21"], "path11": tree["path11"] / "java"}  # fmt: skip
    want = None
    for name, exe in exes.items():
        if picked and os.path.isabs(picked) and Path(picked).resolve() == exe.resolve():
            want = (_MAJOR_OF[name], exe)
    major, path = check(env)
    if want is None:
        assert major is None, f"launcher picks {picked!r} (no usable Java); check: {major} {path}"
    else:
        assert major == want[0], f"launcher picks {picked!r}; check_java ran {path} ({major})"
        assert same_file(path, want[1]), (picked, path)


# ---------------------------------------------------------------- B1 through gippa run (exit codes)


@pytest.fixture
def inputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    monkeypatch.setenv("GIPPA_ROOT", str(fake_gippa_root(tmp_path / "gippa_root")))
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


def set_path(monkeypatch: pytest.MonkeyPatch, *front: Path | str) -> None:
    """``front`` first, then the host ``PATH`` without any ``java`` or ``nextflow``."""
    rest = host_path_without("nextflow", "java")
    monkeypatch.setenv("PATH", os.pathsep.join([*(str(f) for f in front), *rest]))


@pytest.fixture
def nf_args(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    args = tmp_path / "nf_args.txt"
    monkeypatch.setenv("FAKE_NF_ARGS", str(args))
    return args


def _run_case(
    case: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path | None, int | None]:
    """Set up one B1 case; returns (the java Nextflow would run, its major) for the refusals."""
    fakebin = tmp_path / "fakebin"
    make_exe(fakebin / "nextflow", FAKE_NF)
    if case == "nxf_java_home_21_over_java_home_11":
        make_exe(fakebin / "java", java_script(11))
        monkeypatch.setenv("NXF_JAVA_HOME", str(make_jdk(tmp_path / "jdk21", 21)))
        monkeypatch.setenv("JAVA_HOME", str(make_jdk(tmp_path / "jdk11", 11)))
    elif case == "java_cmd_21_over_path_11":
        make_exe(fakebin / "java", java_script(11))
        monkeypatch.setenv("JAVA_CMD", str(make_exe(tmp_path / "cmd" / "java", java_script(21))))
    elif case == "java_home_nonexistent_path_21":
        make_exe(fakebin / "java", java_script(21))
        monkeypatch.setenv("JAVA_HOME", "/nonexistent")
    elif case == "non_executable_java_cmd_falls_through_to_path_21":
        make_exe(fakebin / "java", java_script(21))
        bad = make_exe(tmp_path / "cmd" / "java", java_script(11), mode=0o644)
        monkeypatch.setenv("JAVA_CMD", str(bad))
    elif case == "non_executable_java_cmd_falls_through_to_path_11":
        old = make_exe(fakebin / "java", java_script(11))
        monkeypatch.setenv("JAVA_CMD", str(make_exe(tmp_path / "cmd" / "java", java_script(21),
                                                    mode=0o644)))  # fmt: skip
        set_path(monkeypatch, fakebin)
        return old, 11
    elif case == "java_cmd_11_over_path_21":
        make_exe(fakebin / "java", java_script(21))
        old = make_exe(tmp_path / "cmd" / "java", java_script(11))
        monkeypatch.setenv("JAVA_CMD", str(old))
        set_path(monkeypatch, fakebin)
        return old, 11
    elif case == "nxf_java_home_11_unsets_java_cmd_21":
        make_exe(fakebin / "java", java_script(21))
        jdk = make_jdk(tmp_path / "jdk11", 11)
        monkeypatch.setenv("NXF_JAVA_HOME", str(jdk))
        monkeypatch.setenv("JAVA_CMD", str(make_exe(tmp_path / "cmd" / "java", java_script(21))))
        set_path(monkeypatch, fakebin)
        return jdk / "bin" / "java", 11
    else:  # pragma: no cover
        raise AssertionError(case)
    set_path(monkeypatch, fakebin)
    return None, None


@pytest.mark.parametrize(
    "case",
    [
        "nxf_java_home_21_over_java_home_11",
        "java_cmd_21_over_path_11",
        "java_home_nonexistent_path_21",
        "non_executable_java_cmd_falls_through_to_path_21",
    ],
)
def test_run_starts_nextflow_when_the_launcher_would_find_java_17(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, nf_args: Path,
    case: str,
) -> None:  # fmt: skip
    """Review B1: the three false refusals (and a non-executable ``JAVA_CMD``) exit 0."""
    _run_case(case, tmp_path, monkeypatch)
    code, text = run_cli(inputs, tmp_path / "OUT")
    assert code == ExitCode.SUCCESS, text
    assert JAVA_HEAD not in text, text
    assert nf_args.is_file(), "Nextflow was not started"


@pytest.mark.parametrize(
    "case",
    [
        "non_executable_java_cmd_falls_through_to_path_11",
        "java_cmd_11_over_path_21",
        "nxf_java_home_11_unsets_java_cmd_21",
    ],
)
def test_run_refuses_the_old_java_the_launcher_would_run(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, nf_args: Path,
    case: str,
) -> None:  # fmt: skip
    """The other direction: the Java the launcher picks is 11 although a 21 is around: exit 2,
    naming that Java, and Nextflow never starts."""
    exe, major = _run_case(case, tmp_path, monkeypatch)
    code, text = run_cli(inputs, tmp_path / "OUT")
    assert code == ExitCode.ENVIRONMENT, text
    assert JAVA_HEAD in text, text
    assert f"(found {major} at {exe})" in text, text
    assert not nf_args.exists(), "Nextflow must not start without Java 17+"


def test_setup_warning_follows_java_cmd(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """``gippa setup`` checks the same Java: ``JAVA_CMD``=21 with ``PATH`` java 11 -> no warning."""
    isolate(tmp_path, monkeypatch)
    fakebin = tmp_path / "fakebin"
    make_exe(fakebin / "java", java_script(11))
    make_exe(fakebin / "nextflow", FAKE_NF)
    monkeypatch.setenv("JAVA_CMD", str(make_exe(tmp_path / "cmd" / "java", java_script(21))))
    set_path(monkeypatch, fakebin)
    install_fake_net(monkeypatch)
    code, out = setup_cli(["--db-dir", str(tmp_path / "db"), *SKIP_DBS, "--no-nextflow"])
    assert code == 0, out
    assert JAVA_HEAD not in out, out


def test_conftest_fake_prereqs_unsets_nxf_java_home() -> None:
    """The autouse ``_fake_prereqs`` fixture also drops ``NXF_JAVA_HOME`` (it would override the
    fake ``java`` of every other unit test)."""
    text = (REPO / "tests" / "unit" / "conftest.py").read_text(encoding="utf-8")
    fixture = text.split("def _fake_prereqs(", 1)[1]
    assert '"NXF_JAVA_HOME"' in fixture, "tests/unit/conftest.py _fake_prereqs keeps NXF_JAVA_HOME"


# ======================================== B2: the manifest records the launcher


def run_manifest(out: Path) -> dict[str, Any]:
    path = out / "run_manifest.json"
    assert path.is_file(), f"no {path}"
    data = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def launcher_record(m: Mapping[str, Any]) -> dict[str, Any]:
    launch = m.get("launch")
    assert isinstance(launch, dict), f"run_manifest.json has no `launch` object: {launch!r}"
    rec = launch.get("nextflow")
    assert isinstance(rec, dict), f"run_manifest.json has no `launch.nextflow`: {launch!r}"
    assert set(rec) <= {"path", "source", "version"}, rec
    return rec


def test_manifest_records_a_path_launcher(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, nf_args: Path
) -> None:
    from helpers.manifest_data import validate

    fakebin = tmp_path / "fakebin"
    make_exe(fakebin / "java", java_script(21))
    on_path = make_exe(fakebin / "nextflow", FAKE_NF)
    make_exe(inputs["db"] / "bin" / "nextflow", FAKE_NF)  # PATH wins over it
    set_path(monkeypatch, fakebin)
    code, text = run_cli(inputs, tmp_path / "OUT")
    assert code == ExitCode.SUCCESS, text
    m = run_manifest(tmp_path / "OUT")
    rec = launcher_record(m)
    assert rec["source"] == "path", rec
    assert Path(rec["path"]).is_absolute() and Path(rec["path"]) == on_path.absolute(), rec
    validate(m)


def test_manifest_records_a_db_dir_launcher(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, nf_args: Path
) -> None:
    from helpers.manifest_data import validate

    fakebin = tmp_path / "fakebin"
    make_exe(fakebin / "java", java_script(21))
    in_db = make_exe(inputs["db"] / "bin" / "nextflow", FAKE_NF)
    set_path(monkeypatch, fakebin)
    code, text = run_cli(inputs, tmp_path / "OUT")
    assert code == ExitCode.SUCCESS, text
    m = run_manifest(tmp_path / "OUT")
    rec = launcher_record(m)
    assert rec["source"] == "db_dir", rec
    assert Path(rec["path"]).is_absolute() and Path(rec["path"]) == in_db.absolute(), rec
    validate(m)


def test_manifest_launcher_path_is_absolute_for_a_relative_path_entry(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, nf_args: Path
) -> None:
    fakebin = tmp_path / "relbin"
    make_exe(fakebin / "java", java_script(21))
    exe = make_exe(fakebin / "nextflow", FAKE_NF)
    monkeypatch.chdir(tmp_path)
    set_path(monkeypatch, "relbin")
    code, text = run_cli(inputs, tmp_path / "OUT")
    assert code == ExitCode.SUCCESS, text
    rec = launcher_record(run_manifest(tmp_path / "OUT"))
    assert rec["source"] == "path"
    assert Path(rec["path"]).is_absolute(), rec
    assert Path(rec["path"]).resolve() == exe.resolve(), rec


def test_launcher_record_is_volatile_and_outside_the_provenance_digest(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, nf_args: Path
) -> None:
    """Which launcher ran is machine-specific: ``stable_view`` drops it and two runs that differ
    only in the launcher have the same ``provenance_digest``."""
    from gippa.manifest import stable_view

    fakebin = tmp_path / "fakebin"
    make_exe(fakebin / "java", java_script(21))
    make_exe(inputs["db"] / "bin" / "nextflow", FAKE_NF)
    set_path(monkeypatch, fakebin)
    code, text = run_cli(inputs, tmp_path / "OUT1")
    assert code == ExitCode.SUCCESS, text
    make_exe(fakebin / "nextflow", FAKE_NF)
    code, text = run_cli(inputs, tmp_path / "OUT2")
    assert code == ExitCode.SUCCESS, text
    m1, m2 = run_manifest(tmp_path / "OUT1"), run_manifest(tmp_path / "OUT2")
    assert launcher_record(m1)["source"] == "db_dir" and launcher_record(m2)["source"] == "path"
    assert m1["provenance_digest"] == m2["provenance_digest"]
    stable = stable_view(m1, m1["volatile"])
    assert "nextflow" not in (stable.get("launch") or {}), stable.get("launch")


def test_schema_pins_the_launcher_source_vocabulary(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, nf_args: Path
) -> None:
    from helpers.manifest_data import schema_errors

    fakebin = tmp_path / "fakebin"
    make_exe(fakebin / "java", java_script(21))
    make_exe(fakebin / "nextflow", FAKE_NF)
    set_path(monkeypatch, fakebin)
    code, text = run_cli(inputs, tmp_path / "OUT")
    assert code == ExitCode.SUCCESS, text
    m = run_manifest(tmp_path / "OUT")
    launcher_record(m)
    assert schema_errors(m) == []
    for bad in ({"source": "PATH"}, {"source": "elsewhere"}, {"path": 3}, {"extra": "x"}):
        doc = json.loads(json.dumps(m))
        doc["launch"]["nextflow"].update(bad)
        assert schema_errors(doc), f"the schema accepts launch.nextflow with {bad}"
    doc = json.loads(json.dumps(m))
    del doc["launch"]["nextflow"]["source"]
    assert schema_errors(doc), "launch.nextflow.source is required"


# ================================================================ B3 and nits 2-3: gippa setup

NF_VERSION = "26.10.0"
NF_PATH = f"/nextflow-io/nextflow/releases/download/v{NF_VERSION}/nextflow-{NF_VERSION}-dist"
ASSET_PATH = "/github-production-release-asset/43446421/0f5e7c1a-26-10-0-dist"
# M11-15c2 turns antiSMASH, the resource files and the images on by default; these tests are
# about Java and the Nextflow launcher, so everything else is opted out (opt_outs: only flags the
# current setup accepts). Without --no-images the real apptainer would pull every image.
SKIP_DBS = ["--bakta", "skip", "--no-genomad", "--no-promotech-models",
            *opt_outs(*NO_RESOURCE_FLAGS, "images")]  # fmt: skip
FAKE_DIST = b"#!/bin/sh\n# fake nextflow-26.10.0-dist (test bytes, not Nextflow)\nexit 0\n"


def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def pin_nextflow(monkeypatch: pytest.MonkeyPatch, hex_digest: str) -> None:
    from gippa.setupdb import items

    original = items.SOURCES
    patched = MappingProxyType({**original, "nextflow": dataclasses.replace(
        original["nextflow"], checksum=f"sha256:{hex_digest}")})  # fmt: skip
    for name, mod in list(sys.modules.items()):
        if (name == "gippa" or name.startswith("gippa.")) and getattr(
            mod, "SOURCES", None
        ) is original:
            monkeypatch.setattr(mod, "SOURCES", patched)


def serve_nextflow(net: FakeNet) -> None:
    from gippa.setupdb.items import SOURCES

    host = SOURCES["nextflow"].redirect_hosts[0]
    loc = f"https://{host}{ASSET_PATH}?sp=r&sig=x"
    net.routes[("github.com", NF_PATH)] = Route(302, {"Location": loc})
    net.routes[(host, ASSET_PATH)] = Route(200, {"Content-Type": "application/octet-stream"},
                                           FAKE_DIST)  # fmt: skip


@pytest.fixture
def setup_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> FakeNet:
    """No ``nextflow`` anywhere on ``PATH``, a fake Java 21, the launcher asset served."""
    isolate(tmp_path, monkeypatch)
    fakebin = tmp_path / "fakebin"
    make_exe(fakebin / "java", java_script(21))
    set_path(monkeypatch, fakebin)
    pin_nextflow(monkeypatch, sha256_hex(FAKE_DIST))
    net = install_fake_net(monkeypatch)
    serve_nextflow(net)
    old = os.umask(0o022)
    yield net
    os.umask(old)


def test_b3_check_never_plans_the_launcher(setup_env: FakeNet, tmp_path: Path) -> None:
    db = tmp_path / "db"
    _code, out = setup_cli(["--db-dir", str(db), *SKIP_DBS, "--check"])
    assert setup_env.contacted == [], f"--check contacted {setup_env.contacted}: {out}"
    assert not (db / "bin" / "nextflow").exists(), out
    assert "nextflow" not in manifest_items(db), out


def test_nit2_adopted_launcher_gets_its_x_bit_back(setup_env: FakeNet, tmp_path: Path) -> None:
    db = tmp_path / "db"
    exe = db / "bin" / "nextflow"
    exe.parent.mkdir(parents=True)
    exe.write_bytes(FAKE_DIST)
    exe.chmod(0o644)
    code, out = setup_cli(["--db-dir", str(db), *SKIP_DBS])
    assert code == 0, out
    assert setup_env.contacted == [], "a verified launcher is adopted, not downloaded"
    assert stat.S_IMODE(exe.stat().st_mode) == 0o755
    assert manifest_items(db)["nextflow"]["status"] == "ok"


def test_nit3_adopted_launcher_keeps_its_retrieved_date(setup_env: FakeNet, tmp_path: Path) -> None:
    db = tmp_path / "db"
    code, out = setup_cli(["--db-dir", str(db), *SKIP_DBS])
    assert code == 0, out
    assert (db / "bin" / "nextflow").read_bytes() == FAKE_DIST
    path = db / "db_manifest.json"
    doc = json.loads(path.read_text(encoding="utf-8"))
    for item in doc["items"]:
        if item["name"] == "nextflow":
            item["retrieved"] = "2001-02-03"
    path.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    setup_env.contacted.clear()
    code, out = setup_cli(["--db-dir", str(db), *SKIP_DBS])
    assert code == 0, out
    assert setup_env.contacted == []
    rec = manifest_items(db)["nextflow"]
    assert (rec["status"], rec["retrieved"]) == ("ok", "2001-02-03"), rec


# ================================================================ nit 1: an absolute argv[0]


def test_nit1_relative_path_entry_launches_an_absolute_path(
    inputs: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, nf_args: Path
) -> None:
    """``PATH=relbin:...``: Nextflow runs in ``OUT``, so ``argv[0]`` must be absolute."""
    fakebin = tmp_path / "relbin"
    make_exe(fakebin / "java", java_script(21))
    exe = make_exe(fakebin / "nextflow", FAKE_NF)
    monkeypatch.chdir(tmp_path)
    set_path(monkeypatch, "relbin")
    code, text = run_cli(inputs, tmp_path / "OUT")
    assert code == ExitCode.SUCCESS, text
    argv0 = nf_args.with_name(nf_args.name + ".argv0").read_text(encoding="utf-8").strip()
    assert os.path.isabs(argv0), argv0
    assert Path(argv0).resolve() == exe.resolve()
