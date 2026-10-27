"""M11-15d: the Java check for Nextflow (docs/specs/M11-15-portable-clusters.md item 5).

Java is never downloaded. ``java -version`` (``$JAVA_HOME/bin/java`` when ``JAVA_HOME`` is set,
the launcher Nextflow itself runs; else ``java`` on ``PATH``) is run with a 10 s timeout and the
major version parsed. Any nonzero exit (the macOS ``/usr/bin/java`` stub without a JDK) or
unparseable output counts as "not found". Examples 5.3, 5.4, 5.7; AC "Java version parser tests
(``"21.0.1"``, ``"17"``, ``"1.8.0_292"``, ``"11.0.2"``, the macOS stub's exit 1, garbage → not
found)".

Interface pinned here (new module ``gippa.runner.prereqs``, the prerequisites of a Nextflow
launch, used by ``gippa setup`` (warning) and ``gippa run`` / ``gippa demo --full`` (exit 2)):

- ``JAVA_MIN == 17``, ``JAVA_TIMEOUT_S == 10``;
- ``parse_java_major(text) -> int | None`` on the combined stdout + stderr of ``java -version``;
- ``JavaCheck(major: int | None, path: str | None)`` (frozen dataclass);
- ``check_java(environ) -> JavaCheck`` runs ``[java, "-version"]`` through ``subprocess.run``
  with ``timeout=JAVA_TIMEOUT_S``; ``environ`` gives ``PATH`` and ``JAVA_HOME``;
- ``java_problem(check) -> str | None``: ``None`` for Java >= 17, else the spec's message.

Never the real ``java``: every check runs a fake ``java`` script in ``tmp_path``.
"""

from __future__ import annotations

import os
import stat
import subprocess
import tomllib
from pathlib import Path
from typing import Any

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

pytestmark = pytest.mark.real_prereqs  # these tests build their own PATH

REPO = Path(__file__).resolve().parents[2]
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
ADVICE = f"install it, e.g.: micromamba create -n {PKG}-java -c conda-forge openjdk=21"
HEAD = "Java 17 or newer is needed for Nextflow ("

# Real `java -version` outputs (stderr), one per release line the spec names.
OPENJDK_21 = (
    'openjdk version "21.0.1" 2023-10-17\n'
    "OpenJDK Runtime Environment (build 21.0.1+12-29)\n"
    "OpenJDK 64-Bit Server VM (build 21.0.1+12-29, mixed mode, sharing)\n"
)
OPENJDK_17 = (
    'openjdk version "17" 2021-09-14\n'
    "OpenJDK Runtime Environment (build 17+35-2724)\n"
    "OpenJDK 64-Bit Server VM (build 17+35-2724, mixed mode, sharing)\n"
)
ORACLE_8 = (
    'java version "1.8.0_292"\n'
    "Java(TM) SE Runtime Environment (build 1.8.0_292-b10)\n"
    "Java HotSpot(TM) 64-Bit Server VM (build 25.292-b10, mixed mode)\n"
)
OPENJDK_11 = (
    'openjdk version "11.0.2" 2019-01-15\n'
    "OpenJDK Runtime Environment 18.9 (build 11.0.2+9)\n"
    "OpenJDK 64-Bit Server VM 18.9 (build 11.0.2+9, mixed mode)\n"
)
MACOS_STUB = (
    "The operation couldn't be completed. Unable to locate a Java Runtime.\n"
    "Please visit http://www.java.com for information on installing Java.\n"
)


def _prereqs() -> Any:
    try:
        from gippa.runner import prereqs
    except ImportError as exc:  # 15d adds the module
        pytest.fail(f"gippa.runner.prereqs is missing (M11-15d Java check): {exc}")
    return prereqs


def fake_java(bindir: Path, stderr: str, code: int = 0, stdout: str = "") -> Path:
    """A ``java`` that prints ``stderr`` (and ``stdout``) for ``-version`` and exits ``code``;
    every call is logged to ``<bindir>/java.calls``."""
    bindir.mkdir(parents=True, exist_ok=True)
    (bindir / "java.err").write_text(stderr, encoding="utf-8")
    (bindir / "java.out").write_text(stdout, encoding="utf-8")
    exe = bindir / "java"
    exe.write_text(
        "#!/bin/sh\n"
        f"printf '%s\\n' \"$*\" >> '{bindir}/java.calls'\n"
        f"while IFS= read -r l; do printf '%s\\n' \"$l\"; done < '{bindir}/java.out'\n"
        f"while IFS= read -r l; do printf '%s\\n' \"$l\" >&2; done < '{bindir}/java.err'\n"
        f"exit {code}\n",
        encoding="utf-8",
    )
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return exe


def env_for(*dirs: Path, java_home: Path | None = None) -> dict[str, str]:
    env = {"PATH": os.pathsep.join(str(d) for d in dirs)}
    if java_home is not None:
        env["JAVA_HOME"] = str(java_home)
    return env


# ---------------------------------------------------------------- the parser (AC)


def test_constants() -> None:
    p = _prereqs()
    assert p.JAVA_MIN == 17
    assert p.JAVA_TIMEOUT_S == 10


@pytest.mark.parametrize(
    ("text", "major"),
    [
        (OPENJDK_21, 21),
        (OPENJDK_17, 17),
        (ORACLE_8, 8),
        (OPENJDK_11, 11),
        ('openjdk version "21.0.1"', 21),
        ('openjdk version "17"', 17),
        ('java version "1.8.0_292"', 8),
        ('openjdk version "11.0.2"', 11),
        ('openjdk version "1.7.0_80"', 7),
        ('openjdk version "25" 2025-09-16 LTS\n', 25),
        # a JAVA_TOOL_OPTIONS line comes first on many clusters
        ("Picked up JAVA_TOOL_OPTIONS: -Xmx2g\n" + OPENJDK_21, 21),
    ],
)
def test_parse_java_major(text: str, major: int) -> None:
    assert _prereqs().parse_java_major(text) == major


@pytest.mark.parametrize(
    "text",
    [
        "",
        MACOS_STUB,
        "garbage",
        "java: command not found",
        'openjdk version ""',
        'openjdk version "x.y.z"',
        "version 21.0.1",  # no quoted version string
    ],
)
def test_garbage_is_not_found(text: str) -> None:
    assert _prereqs().parse_java_major(text) is None


@settings(max_examples=200, deadline=None)
@given(
    vendor=st.sampled_from(["openjdk", "java"]),
    major=st.integers(min_value=9, max_value=60),
    minor=st.integers(min_value=0, max_value=30),
    patch=st.integers(min_value=0, max_value=30),
    form=st.sampled_from(["major", "minor", "patch"]),
    tail=st.sampled_from(["", " 2024-01-16", " 2024-01-16 LTS"]),
)
def test_modern_version_property(
    vendor: str, major: int, minor: int, patch: int, form: str, tail: str
) -> None:
    """Java 9+ (JEP 223): the major is the first component."""
    v = {"major": f"{major}", "minor": f"{major}.{minor}", "patch": f"{major}.{minor}.{patch}"}
    text = f'{vendor} version "{v[form]}"{tail}\nOpenJDK Runtime Environment\n'
    assert _prereqs().parse_java_major(text) == major


@settings(max_examples=100, deadline=None)
@given(minor=st.integers(min_value=1, max_value=8), update=st.integers(min_value=0, max_value=999))
def test_legacy_version_property(minor: int, update: int) -> None:
    """Java <= 8 reports ``1.<major>.0_<update>``: the major is the second component."""
    text = f'java version "1.{minor}.0_{update}"\n'
    assert _prereqs().parse_java_major(text) == minor


@settings(max_examples=300, deadline=None)
@given(st.text())
def test_parser_never_raises(text: str) -> None:
    got = _prereqs().parse_java_major(text)
    assert got is None or (isinstance(got, int) and got >= 1)


# ---------------------------------------------------------------- running java (fake)


@pytest.fixture(autouse=True)
def _no_host_java(monkeypatch: pytest.MonkeyPatch) -> None:
    for var in ("JAVA_HOME", "JAVA_CMD", "JAVA_TOOL_OPTIONS", "_JAVA_OPTIONS", "JDK_JAVA_OPTIONS"):
        monkeypatch.delenv(var, raising=False)


@pytest.mark.parametrize(
    ("stderr", "major"), [(OPENJDK_21, 21), (OPENJDK_17, 17), (OPENJDK_11, 11), (ORACLE_8, 8)]
)
def test_check_java_on_path(tmp_path: Path, stderr: str, major: int) -> None:
    exe = fake_java(tmp_path / "bin", stderr)
    got = _prereqs().check_java(env_for(tmp_path / "bin"))
    assert (got.major, got.path) == (major, str(exe))
    assert (tmp_path / "bin" / "java.calls").read_text().splitlines() == ["-version"]


def test_version_on_stdout_is_read_too(tmp_path: Path) -> None:
    fake_java(tmp_path / "bin", "", stdout=OPENJDK_21)
    assert _prereqs().check_java(env_for(tmp_path / "bin")).major == 21


def test_no_java_on_path_is_not_found(tmp_path: Path) -> None:
    (tmp_path / "empty").mkdir()
    got = _prereqs().check_java(env_for(tmp_path / "empty"))
    assert got.major is None
    assert got.path is None


def test_example_5_7_macos_stub_exit_1_is_not_found(tmp_path: Path) -> None:
    fake_java(tmp_path / "usr" / "bin", MACOS_STUB, code=1)
    assert _prereqs().check_java(env_for(tmp_path / "usr" / "bin")).major is None


def test_nonzero_exit_with_a_version_line_is_not_found(tmp_path: Path) -> None:
    """Any nonzero exit counts as not found, whatever was printed."""
    fake_java(tmp_path / "bin", OPENJDK_21, code=1)
    assert _prereqs().check_java(env_for(tmp_path / "bin")).major is None


def test_unparseable_output_is_not_found(tmp_path: Path) -> None:
    fake_java(tmp_path / "bin", "garbage\n")
    assert _prereqs().check_java(env_for(tmp_path / "bin")).major is None


def test_java_home_wins_over_path(tmp_path: Path) -> None:
    """``$JAVA_HOME/bin/java`` is the launcher Nextflow runs when ``JAVA_HOME`` is set."""
    fake_java(tmp_path / "path_bin", OPENJDK_11)
    home_java = fake_java(tmp_path / "jdk21" / "bin", OPENJDK_21)
    got = _prereqs().check_java(env_for(tmp_path / "path_bin", java_home=tmp_path / "jdk21"))
    assert (got.major, got.path) == (21, str(home_java))
    assert not (tmp_path / "path_bin" / "java.calls").exists()


def test_check_java_uses_the_10_s_timeout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fake_java(tmp_path / "bin", OPENJDK_21)
    real = subprocess.run
    seen: list[Any] = []

    def spy(*args: Any, **kw: Any) -> Any:
        seen.append(kw.get("timeout"))
        return real(*args, **kw)

    monkeypatch.setattr(subprocess, "run", spy)
    assert _prereqs().check_java(env_for(tmp_path / "bin")).major == 21
    assert seen == [10]


def test_a_hanging_java_is_not_found(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fake_java(tmp_path / "bin", OPENJDK_21)

    def hang(argv: Any, *_a: Any, **kw: Any) -> Any:
        raise subprocess.TimeoutExpired(argv, kw.get("timeout") or 10)

    monkeypatch.setattr(subprocess, "run", hang)
    assert _prereqs().check_java(env_for(tmp_path / "bin")).major is None


def test_a_java_that_cannot_start_is_not_found(tmp_path: Path) -> None:
    bindir = tmp_path / "bin"
    bindir.mkdir()
    exe = bindir / "java"
    exe.write_text("#!/nonexistent/interpreter\n")
    exe.chmod(0o755)
    assert _prereqs().check_java(env_for(bindir)).major is None


# ---------------------------------------------------------------- the message


@pytest.mark.parametrize("major", [17, 21, 25])
def test_new_enough_java_has_no_problem(major: int) -> None:
    p = _prereqs()
    assert p.java_problem(p.JavaCheck(major=major, path="/opt/jdk/bin/java")) is None


@pytest.mark.parametrize(("major", "path"), [(11, "/usr/bin/java"), (8, "/opt/jdk8/bin/java")])
def test_examples_5_3_5_4_old_java_message(major: int, path: str) -> None:
    p = _prereqs()
    assert p.java_problem(p.JavaCheck(major=major, path=path)) == (
        f"Java 17 or newer is needed for Nextflow (found {major} at {path}); {ADVICE}"
    )


@pytest.mark.parametrize("path", [None, "/usr/bin/java"])
def test_not_found_message(path: str | None) -> None:
    p = _prereqs()
    msg = p.java_problem(p.JavaCheck(major=None, path=path))
    assert msg is not None
    assert msg.startswith(HEAD) and msg.endswith(ADVICE), msg
    assert "not found" in msg
    assert "found None" not in msg


def test_end_to_end_old_java_message(tmp_path: Path) -> None:
    exe = fake_java(tmp_path / "bin", OPENJDK_11)
    p = _prereqs()
    msg = p.java_problem(p.check_java(env_for(tmp_path / "bin")))
    assert msg == f"Java 17 or newer is needed for Nextflow (found 11 at {exe}); {ADVICE}"
