"""Prerequisites of a Nextflow launch: Java 17+ and the ``nextflow`` launcher (M11-15 item 5).

Java is never downloaded: ``$JAVA_HOME/bin/java`` (the Java the Nextflow launcher runs when
``JAVA_HOME`` is set), else ``java`` on ``PATH``, is run as ``java -version`` with a 10 s timeout
and its major version parsed. A nonzero exit (the macOS ``/usr/bin/java`` stub without a JDK),
a timeout, a launcher that cannot start or unparseable output count as "not found". ``gippa
setup`` warns; ``gippa run`` exits 2. ``nextflow`` is looked up on ``PATH``, else
``<db_dir>/bin/nextflow`` (the ``-dist`` launcher ``gippa setup`` installs).
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "JAVA_MIN",
    "JAVA_TIMEOUT_S",
    "LAUNCHER",
    "JavaCheck",
    "check_java",
    "find_nextflow",
    "java_problem",
    "nextflow_version",
    "parse_java_major",
]

JAVA_MIN = 17
JAVA_TIMEOUT_S = 10
LAUNCHER = "bin/nextflow"  # under the DB directory (the db_sources.yaml ``nextflow`` row)
_PKG = __name__.split(".", 1)[0]
_ADVICE = f"install it, e.g.: micromamba create -n {_PKG}-java -c conda-forge openjdk=21"
# ``openjdk version "21.0.1" 2023-10-17`` / ``java version "1.8.0_292"``; bounded digit runs
_JAVA_VERSION = re.compile(r'\bversion "(\d{1,9})(?:\.(\d{1,9}))?[^"\n]*"')
_NF_VERSION = re.compile(r"\bversion (\d+\.\d+\.\d+(?:\.\d+)?)")


@dataclass(frozen=True, slots=True)
class JavaCheck:
    """The Java found: its major version (``None`` = not found or unusable) and its path."""

    major: int | None
    path: str | None


def parse_java_major(text: str) -> int | None:
    """Major version from ``java -version`` output: ``1.<n>...`` is ``n`` (Java <= 8), else the
    first component (JEP 223); ``None`` when no quoted version is found."""
    m = _JAVA_VERSION.search(text)
    if m is None:
        return None
    first = int(m.group(1))
    major = int(m.group(2)) if first == 1 and m.group(2) is not None else first
    return major if major >= 1 else None


def _java_exe(environ: Mapping[str, str]) -> str | None:
    home = environ.get("JAVA_HOME")
    if home:
        return str(Path(home) / "bin" / "java")
    return shutil.which("java", path=environ.get("PATH", ""))


def check_java(environ: Mapping[str, str]) -> JavaCheck:
    """Run ``java -version`` (``environ`` gives ``JAVA_HOME`` and ``PATH``) and parse it."""
    exe = _java_exe(environ)
    if exe is None:
        return JavaCheck(None, None)
    try:
        proc = subprocess.run(
            [exe, "-version"],
            capture_output=True,
            text=True,
            errors="replace",
            timeout=JAVA_TIMEOUT_S,
            check=False,
            env=dict(environ),
        )
    except (OSError, subprocess.SubprocessError):
        return JavaCheck(None, exe)
    if proc.returncode != 0:
        return JavaCheck(None, exe)
    return JavaCheck(parse_java_major(f"{proc.stdout}\n{proc.stderr}"), exe)


def java_problem(check: JavaCheck) -> str | None:
    """``None`` for Java >= 17, else the message naming what was found and how to install it."""
    if check.major is not None and check.major >= JAVA_MIN:
        return None
    if check.major is not None:
        found = f"found {check.major} at {check.path}"
    elif check.path is not None:
        found = f"not found: {check.path} -version gave no usable version"
    else:
        found = "not found"
    return f"Java {JAVA_MIN} or newer is needed for Nextflow ({found}); {_ADVICE}"


def find_nextflow(environ: Mapping[str, str], db_dir: Path | None) -> Path | None:
    """``nextflow`` on ``environ["PATH"]``, else ``<db_dir>/bin/nextflow`` when it is an
    executable file, else ``None``."""
    on_path = shutil.which("nextflow", path=environ.get("PATH", ""))
    if on_path is not None:
        return Path(on_path)
    if db_dir is None:
        return None
    exe = db_dir / LAUNCHER
    return exe if exe.is_file() and os.access(exe, os.X_OK) else None


def nextflow_version(exe: Path, environ: Mapping[str, str]) -> str:
    """The version ``<exe> -v`` reports (``nextflow version 25.10.0.10289``), else ``""``."""
    try:
        proc = subprocess.run(
            [str(exe), "-v"],
            capture_output=True,
            text=True,
            errors="replace",
            timeout=JAVA_TIMEOUT_S * 6,
            check=False,
            env=dict(environ),
        )
    except (OSError, subprocess.SubprocessError):
        return ""
    m = _NF_VERSION.search(proc.stdout)
    return m.group(1) if proc.returncode == 0 and m is not None else ""
