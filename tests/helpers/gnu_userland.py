"""Host probes for tests that execute the FarmShare Slurm scripts (GIAP_SPEC.md G19).

``scripts/slurm/*.sbatch`` and friends run only on FarmShare and use the GNU userland:
``sha256sum``, ``flock``, ``realpath -e`` and ``stat -c``. On a host without it (macOS: BSD
``realpath``/``stat``, no ``sha256sum``/``flock``) the scripts cannot run, so such tests skip:

* the ``gnu_userland`` marker (``tests/conftest.py`` skips marked tests at collection);
* ``require_gnu_userland()``, called by the shared script runners (``Sandbox.run`` in
  ``helpers.slurm_sandbox`` and ``M10Sandbox.popen``), so every test that runs a script through
  them, including future ones, skips without needing the marker.

The probe runs once per session (per xdist worker) against ``/usr/bin:/bin``, the system part
of the PATH the sandboxed scripts get.

``x87_longdouble`` is the matching probe for ``scripts/golden/make_wallenius_exact.py``, which
refuses to run unless ``np.longdouble`` is the 80-bit x87 format (absent on e.g. macOS arm64).
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from functools import cache

import pytest

GNU_USERLAND_REASON = "FarmShare/GNU userland only (sha256sum, flock, GNU realpath/stat)"
X87_LONGDOUBLE_REASON = "needs 80-bit x87 long double (absent on e.g. macOS arm64)"


# the system part of the sandbox PATH (``helpers.slurm_sandbox``: ``bin:/usr/bin:/bin``); the
# scripts see these tools, not the host PATH (a Homebrew gnubin there would not reach them)
SCRIPT_PATH = "/usr/bin:/bin"


def _probes(target: str) -> dict[str, list[str]]:
    return {
        "sha256sum": ["sha256sum", target],
        "flock": ["flock", "-n", target, "true"],
        "realpath -e": ["realpath", "-e", "/"],
        "stat -c": ["stat", "-c", "%s", target],
    }


def _works(argv: list[str], path: str) -> bool:
    exe = shutil.which(argv[0], path=path)
    if exe is None:
        return False
    try:
        proc = subprocess.run(
            [exe, *argv[1:]], capture_output=True, timeout=30, check=False, env={"PATH": path}
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return proc.returncode == 0


def probe_gnu_userland(path: str) -> tuple[str, ...]:
    """The GNU userland pieces missing (absent or not GNU) on ``path``; empty when present.

    The target is a private temporary file, so ``flock -n`` never meets another lock holder
    (parallel probes on a shared file fail each other: PR #318 review B1)."""
    fd, target = tempfile.mkstemp(prefix="giap-gnu-probe-")
    try:
        os.close(fd)
        return tuple(name for name, argv in _probes(target).items() if not _works(argv, path))
    finally:
        os.unlink(target)


@cache
def missing_gnu_userland() -> tuple[str, ...]:
    """``probe_gnu_userland`` of ``SCRIPT_PATH``, probed once per session."""
    return probe_gnu_userland(SCRIPT_PATH)


def have_gnu_userland() -> bool:
    return not missing_gnu_userland()


def require_gnu_userland() -> None:
    """Skip the calling test when the host lacks the GNU userland the scripts need."""
    if not have_gnu_userland():
        pytest.skip(GNU_USERLAND_REASON)


def have_x87_longdouble() -> bool:
    import numpy as np

    return bool(np.finfo(np.longdouble).nmant >= 63)
