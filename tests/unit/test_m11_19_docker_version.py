"""M11-19: the Docker / podman CLI version that ``docker run --pull=never`` needs.

docs/specs/M11-19-nextflow-docker-pull-never.md item 2 (``parse_docker_version``,
``docker_client_version``, ``pull_never_problem``, ``DOCKER_MIN``/``PODMAN_MIN``), item 4
(podman-docker), examples 1-7, acceptance criterion 1, Q2 (an unknown version never blocks) and
Q3 (podman floor 3.0).

Interface pinned here (implementer, ``src/gippa/preflight/environment.py``):

- ``parse_docker_version(text) -> ("docker" | "podman", (major, minor)) | None``: pure; the first
  line of ``docker --version`` stdout must match ``^(Docker|podman) version (\\d+)\\.(\\d+)``
  (case-sensitive), anything else is ``None``.
- ``docker_client_version()``: runs ``docker --version`` (found with ``shutil.which("docker")``)
  once; ``None`` without docker on ``PATH``, on a failed call or unparseable output.
- ``pull_never_problem(version) -> str``: ``""`` when supported or unknown, else the rule 1
  reason (``<family> <X.Y> is older than <MIN>, ...``).

No network, no real docker: a logging fake ``docker`` sits first on ``PATH`` (after the autouse
container guard of tests/conftest.py), and ``$GIPPA_ROOT``/``$GIPPA_DB_DIR`` are unset.
"""

from __future__ import annotations

import json
import os
import stat
import sys
from pathlib import Path
from typing import Any

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

TOO_OLD = (
    "{family} {x}.{y} is older than {min}, the first version whose docker run has --pull=never, "
    "and runs never pull images; upgrade Docker (Docker >= 20.10, or podman >= 3.0 with "
    "podman-docker)"
)

FAKE_DOCKER = r'''#!{python} -I
"""Fake docker (tests only): ``--version`` prints FAKE_DOCKER_VERSION, exits FAKE_DOCKER_RC."""
import json, os, sys
argv = sys.argv[1:]
with open(os.environ["FAKE_DOCKER_LOG"], "a") as fh:
    fh.write(json.dumps(argv) + "\n")
if argv == ["--version"]:
    line = os.environ.get("FAKE_DOCKER_VERSION", "")
    if line.startswith("podman"):
        sys.stderr.write("Emulate Docker CLI using podman. Create /etc/containers/nodocker to "
                         "quiet msg.\n")
    if line:
        sys.stdout.write(line + "\n")
    sys.exit(int(os.environ.get("FAKE_DOCKER_RC", "0")))
sys.stderr.write("fake docker: unexpected call\n")
sys.exit(3)
'''


def _env() -> Any:
    from gippa.preflight import environment

    return environment


def _fn(name: str) -> Any:
    fn = getattr(_env(), name, None)
    assert fn is not None, f"gippa.preflight.environment.{name} is not implemented (M11-19)"
    return fn


@pytest.fixture(autouse=True)
def fake_docker(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A logging fake ``docker`` first on ``PATH``; returns its log."""
    for var in ("GIPPA_DB_DIR", "GIPPA_ROOT", "FAKE_DOCKER_VERSION", "FAKE_DOCKER_RC"):
        monkeypatch.delenv(var, raising=False)
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "docker"
    exe.write_text(FAKE_DOCKER.format(python=sys.executable))
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    log = tmp_path / "docker.log"
    monkeypatch.setenv("FAKE_DOCKER_LOG", str(log))
    return log


def _calls(log: Path) -> list[list[str]]:
    return [json.loads(x) for x in log.read_text().splitlines()] if log.exists() else []


# ---------------------------------------------------------------- constants


def test_floors() -> None:
    assert _env().DOCKER_MIN == (20, 10)
    assert _env().PODMAN_MIN == (3, 0)


def test_setup_docker_hint() -> None:
    assert _env().SETUP_DOCKER == "run: gippa setup --images --profile docker"


# ---------------------------------------------------------------- examples 1-7


@pytest.mark.parametrize(
    ("text", "want"),
    [
        ("Docker version 27.3.1, build ce12230", ("docker", (27, 3))),  # 1
        ("Docker version 20.10.24+dfsg1, build 297e128", ("docker", (20, 10))),  # 2
        ("Docker version 19.03.15, build 99e3ed8919", ("docker", (19, 3))),  # 3
        ("Docker version 1.13.1, build 7d71120/1.13.1", ("docker", (1, 13))),  # 4
        ("podman version 4.9.3", ("podman", (4, 9))),  # 5
        ("podman version 2.2.1", ("podman", (2, 2))),  # 5
        ("Docker version 24.0.7, build 24.0.7-0ubuntu4.1", ("docker", (24, 0))),  # 6
        ("nerdctl version 1.7.6", None),  # 7
        ("", None),  # 7
        ("docker: command not found", None),  # 7
    ],
)
def test_parse_docker_version_examples(text: str, want: Any) -> None:
    assert _fn("parse_docker_version")(text) == want


@pytest.mark.parametrize(
    "text",
    [
        "docker version 27.3.1",
        "Podman version 4.9.3",
        "DOCKER VERSION 27.3.1",
        "Docker version x.y",
    ],
)
def test_parse_is_case_sensitive_and_needs_numbers(text: str) -> None:
    assert _fn("parse_docker_version")(text) is None


def test_parse_reads_the_first_line_only() -> None:
    parse = _fn("parse_docker_version")
    assert parse("Docker version 27.3.1, build ce12230\nsomething else\n") == ("docker", (27, 3))
    assert parse("Emulate Docker CLI using podman\npodman version 4.9.3\n") is None


@pytest.mark.parametrize(
    ("version", "want"),
    [
        (("docker", (27, 3)), ""),  # 1
        (("docker", (20, 10)), ""),  # 2: the floor itself passes
        (("docker", (24, 0)), ""),  # 6
        (("podman", (4, 9)), ""),  # 5
        (None, ""),  # 7: unknown never blocks (Q2)
        (("docker", (19, 3)), TOO_OLD.format(family="docker", x=19, y=3, min="20.10")),  # 3
        (("docker", (1, 13)), TOO_OLD.format(family="docker", x=1, y=13, min="20.10")),  # 4
        (("podman", (2, 2)), TOO_OLD.format(family="podman", x=2, y=2, min="3.0")),  # 5
    ],
)
def test_pull_never_problem_examples(version: Any, want: str) -> None:
    assert _fn("pull_never_problem")(version) == want


def test_example_3_renders_the_parsed_ints() -> None:
    """``19.03`` reads ``19.3``: X.Y comes from the parsed ints."""
    v = _fn("parse_docker_version")("Docker version 19.03.15, build 99e3ed8919")
    got = _fn("pull_never_problem")(v)
    assert got.startswith("docker 19.3 is older than 20.10"), got


def test_example_5_podman_is_not_compared_with_the_docker_floor() -> None:
    v = _fn("parse_docker_version")("podman version 4.9.3")
    assert _fn("pull_never_problem")(v) == ""
    old = _fn("pull_never_problem")(_fn("parse_docker_version")("podman version 2.2.1"))
    assert old.startswith("podman 2.2 is older than 3.0"), old


# ---------------------------------------------------------------- AC 1: property


@settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(
    family=st.sampled_from(["docker", "podman"]),
    major=st.integers(min_value=0, max_value=400),
    minor=st.integers(min_value=0, max_value=400),
)
def test_property_floor(family: str, major: int, minor: int) -> None:
    floor = (20, 10) if family == "docker" else (3, 0)
    got = _fn("pull_never_problem")((family, (major, minor)))
    if (major, minor) >= floor:
        assert got == ""
    else:
        lo = "20.10" if family == "docker" else "3.0"
        assert got == TOO_OLD.format(family=family, x=major, y=minor, min=lo)


@settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(
    family=st.sampled_from(["Docker", "podman"]),
    major=st.integers(min_value=0, max_value=400),
    minor=st.integers(min_value=0, max_value=400),
    patch=st.integers(min_value=0, max_value=400),
    tail=st.sampled_from(["", ", build abc1234", "+dfsg1, build 297e128", "-rc.1"]),
)
def test_property_parse_round_trip(
    family: str, major: int, minor: int, patch: int, tail: str
) -> None:
    text = f"{family} version {major}.{minor}.{patch}{tail}"
    assert _fn("parse_docker_version")(text) == (family.lower(), (major, minor))


# ---------------------------------------------------------------- docker_client_version


@pytest.mark.parametrize(
    ("line", "want"),
    [
        ("Docker version 27.3.1, build ce12230", ("docker", (27, 3))),
        ("Docker version 19.03.15, build 99e3ed8919", ("docker", (19, 3))),
        ("podman version 4.9.3", ("podman", (4, 9))),
        ("nerdctl version 1.7.6", None),
        ("", None),
    ],
)
def test_client_version_runs_docker_version_once(
    fake_docker: Path, monkeypatch: pytest.MonkeyPatch, line: str, want: Any
) -> None:
    monkeypatch.setenv("FAKE_DOCKER_VERSION", line)
    assert _fn("docker_client_version")() == want
    assert _calls(fake_docker) == [["--version"]], "only the client version is asked, once"


def test_client_version_podman_notice_on_stderr_is_ignored(
    fake_docker: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_DOCKER_VERSION", "podman version 3.0.1")
    assert _fn("docker_client_version")() == ("podman", (3, 0))


def test_client_version_nonzero_exit_is_unknown(
    fake_docker: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_DOCKER_VERSION", "Docker version 27.3.1, build ce12230")
    monkeypatch.setenv("FAKE_DOCKER_RC", "1")
    assert _fn("docker_client_version")() is None
    assert _calls(fake_docker) == [["--version"]]


def test_client_version_without_docker_on_path(
    tmp_path: Path, fake_docker: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    empty = tmp_path / "nobin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    assert _fn("docker_client_version")() is None
    assert _calls(fake_docker) == []
