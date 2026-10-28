"""build_envs.sh builds and locks an env that has no pip (job 1902432, miniprot).

miniprot's env (workflow/envs/miniprot.yml) holds no Python. In lock mode the script ran
``micromamba run -n gippa-miniprot pip freeze``, which printed ``exec: pip: not found``; the tool
was reported FAILED with the binary installed but no locks and no marker. An env without
``$prefix/bin/pip`` gets only ``<tool>.explicit.txt`` and the marker; an old ``<tool>.pip.txt``
goes with the other old locks. A ``pip freeze`` that fails in an env that has pip still fails the
build.

The harness is the one of test_build_envs_swap.py, with a micromamba stub that models pip: the fake
``create`` installs ``bin/pip`` unless told ``no_pip_<env>``, and ``run ... pip freeze`` in an env
without ``bin/pip`` fails the way micromamba does.
"""

from __future__ import annotations

import contextlib
import os
import signal
import subprocess
import time
from collections.abc import Iterator
from pathlib import Path

import pytest
from unit.annotate.test_build_envs_swap import (
    FAKE_MICROMAMBA,
    OLD_LOCK,
    OLD_PIP,
    Box,
    _survivors,
)

_CREATE = '    echo new > "$envs/$n/bin/payload"\n'
_PIP = '      *"pip freeze"*)\n'
assert _CREATE in FAKE_MICROMAMBA and _PIP in FAKE_MICROMAMBA
FAKE_MICROMAMBA_PIP = FAKE_MICROMAMBA.replace(
    _CREATE,
    _CREATE
    + '    if [ ! -e "$FAKE_DIR/no_pip_$n" ]; then\n'
    + '      printf \'#!/bin/sh\\nexit 0\\n\' > "$envs/$n/bin/pip"; chmod +x "$envs/$n/bin/pip"\n'
    + "    fi\n",
).replace(
    _PIP,
    _PIP + '        [ -x "$envs/$n/bin/pip" ] || { echo "exec: pip: not found" >&2; exit 1; }\n',
)


@pytest.fixture
def box(tmp_path: Path) -> Iterator[Box]:
    b = Box(tmp_path)
    stub = b.root / "bin" / "micromamba"
    stub.write_text(FAKE_MICROMAMBA_PIP)
    stub.chmod(0o755)
    b.version("bakta", "bakta 1.12.2")
    b.version("miniprot", "0.18")
    b.flag("no_pip", "miniprot")
    yield b
    deadline = time.monotonic() + 10
    while (left := _survivors(b.fake)) and time.monotonic() < deadline:
        time.sleep(0.1)
    for pid in left:
        with contextlib.suppress(ProcessLookupError):
            os.kill(pid, signal.SIGKILL)
    assert not left, f"processes outlived the test: {left}"


def _ok(proc: subprocess.CompletedProcess[str]) -> None:
    assert proc.returncode == 0, proc.stdout + proc.stderr


def _assert_built_without_pip(box: Box, tool: str = "miniprot") -> None:
    prefix = box.prefix(tool)
    assert (prefix / "bin" / "payload").read_text() == "new\n"
    assert not (prefix / "bin" / "pip").exists()
    assert (prefix / ".gippa-built").read_text().strip() == box.yml_sha(tool)
    assert "new-2.0-0.conda" in (box.locks / f"{tool}.explicit.txt").read_text()
    assert not (box.locks / f"{tool}.pip.txt").exists(), "an env without pip has no pip lock"
    left = sorted(p.name for p in box.locks.glob(f"{tool}.*") if p.name.endswith((".tmp", ".prev")))
    assert left == []
    assert not Path(f"{prefix}.gippa-old").exists()


# ---------------------------------------------------------------- the job 1902432 case


def test_first_build_of_env_without_pip_writes_explicit_lock_and_marker(box: Box) -> None:
    """No committed lock yet (lock mode): the build succeeds with one lock and the marker."""
    proc = box.run("miniprot")
    _ok(proc)
    assert "FAILED" not in proc.stderr
    assert "pip: not found" not in proc.stdout + proc.stderr
    _assert_built_without_pip(box)


def test_rebuild_of_env_without_pip_removes_the_old_pip_lock(box: Box) -> None:
    """An old <tool>.pip.txt from an earlier (pip-bearing) build goes with the other old locks."""
    box.old_env("miniprot")
    assert (box.locks / "miniprot.pip.txt").read_text() == OLD_PIP
    _ok(box.run("miniprot"))
    _assert_built_without_pip(box)
    assert (box.locks / "miniprot.explicit.txt").read_text() != OLD_LOCK
    assert sorted(p.name for p in box.envs.glob("gippa-miniprot.gippa-*")) == []


def test_existing_env_without_pip_heals_missing_locks(box: Box) -> None:
    """A verified env (marker matches) whose locks never moved: the next run re-exports them."""
    prefix = box.prefix("miniprot")
    (prefix / "bin").mkdir(parents=True)
    (prefix / "bin" / "payload").write_text("new\n")
    (prefix / ".gippa-built").write_text(box.yml_sha("miniprot") + "\n")
    proc = box.run("miniprot")
    _ok(proc)
    assert box.creates() == [], "a verified env is not rebuilt"
    _assert_built_without_pip(box)


def test_lock_left_aside_by_a_killed_swap_does_not_return_for_env_without_pip(box: Box) -> None:
    """A pip.txt.prev of an earlier killed swap is not put back as the lock of a pip-less env."""
    box.old_env("miniprot")
    (box.locks / "miniprot.pip.txt").rename(box.locks / "miniprot.pip.txt.prev")
    _ok(box.run("miniprot"))
    _assert_built_without_pip(box)


def test_pip_tmp_left_by_a_killed_run_does_not_become_the_lock(box: Box) -> None:
    """A SIGKILL can leave pip.txt.tmp behind; a later pip-less build must not install it."""
    box.old_env("miniprot")
    (box.locks / "miniprot.pip.txt.tmp").write_text("stale==1\n")
    _ok(box.run("miniprot"))
    _assert_built_without_pip(box)
    assert not (box.locks / "miniprot.pip.txt.tmp").exists()


def test_env_without_pip_next_to_env_with_pip(box: Box) -> None:
    """One run, two tools: the pip-less env gets no pip lock, the other keeps its pip lock."""
    box.old_env("bakta")
    proc = box.run("miniprot", "bakta")
    _ok(proc)
    _assert_built_without_pip(box)
    assert (box.locks / "bakta.pip.txt").read_text() == "newpkg==2.0\n"
    assert "new-2.0-0.conda" in (box.locks / "bakta.explicit.txt").read_text()
    assert (box.prefix("bakta") / ".gippa-built").read_text().strip() == box.yml_sha("bakta")


# ---------------------------------------------------------------- unchanged: envs with pip


def test_env_with_pip_still_writes_pip_lock(box: Box) -> None:
    _ok(box.run("bakta"))
    assert (box.locks / "bakta.pip.txt").read_text() == "newpkg==2.0\n"


@pytest.mark.parametrize("old", [False, True])
def test_pip_freeze_failure_in_env_with_pip_still_fails(box: Box, old: bool) -> None:
    """A real pip failure is not mistaken for "no pip": no marker, the old locks are kept."""
    if old:
        box.old_env("bakta")
    box.flag("fail_pip", "bakta")
    proc = box.run("bakta")
    assert proc.returncode == 1, proc.stdout + proc.stderr
    assert "FAILED: bakta" in proc.stderr
    assert not list(box.locks.glob("bakta.*.tmp"))
    if old:
        prefix = box.prefix("bakta")
        assert (prefix / "bin" / "payload").read_text() == "old\n"
        assert (prefix / ".gippa-built").read_text() == "0" * 64 + "\n"
        assert (box.locks / "bakta.explicit.txt").read_text() == OLD_LOCK
        assert (box.locks / "bakta.pip.txt").read_text() == OLD_PIP
    else:
        assert not (box.prefix("bakta") / ".gippa-built").exists()
        assert not (box.locks / "bakta.pip.txt").exists()


def test_no_pip_decision_follows_bin_pip_not_the_tool(box: Box) -> None:
    """miniprot built with a bin/pip whose freeze breaks: that is a real failure, not "no pip"."""
    box.flag("fail_pip", "miniprot")
    (box.fake / "no_pip_gippa-miniprot").unlink()
    proc = box.run("miniprot")
    assert proc.returncode == 1, proc.stdout + proc.stderr
    assert not (box.prefix("miniprot") / ".gippa-built").exists()


# ---------------------------------------------------------------- frozen mode never uses pip


def test_frozen_build_never_calls_pip(box: Box) -> None:
    lock = (
        "@EXPLICIT\nhttps://conda.anaconda.org/conda-forge/linux-64/new-2.0-0.conda#gippa-fgsea\n"
    )
    (box.locks / "fgsea.explicit.txt").write_text(lock)
    box.version("fgsea", "[1] \u20181.36.2\u2019")
    box.flag("no_pip", "fgsea")
    _ok(box.run("fgsea"))
    calls = (box.fake / "calls.log").read_text().splitlines()
    assert not [c for c in calls if "pip" in c.split()], calls
    assert (box.locks / "fgsea.explicit.txt").read_text() == lock
    assert not (box.locks / "fgsea.pip.txt").exists()
