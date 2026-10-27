"""Fake Nextflow prerequisites for the unit tests (M11-15d, spec item 5).

From 15d on, ``gippa setup`` runs ``java -version`` (a warning) and downloads the pinned Nextflow
``-dist`` launcher when ``nextflow`` is not on ``PATH``; ``gippa run`` runs ``java -version``
(exit 2 below Java 17) and looks up ``nextflow``. Unit tests must never run the host's ``java``
or ``nextflow`` nor download anything, so ``tests/unit/conftest.py`` puts the fakes written here
first on ``PATH`` for every test (opt out: ``@pytest.mark.real_prereqs``; tests marked
``needs_nextflow`` keep the real tools).

- ``java -version`` prints an OpenJDK 21 banner on stderr and exits 0.
- ``nextflow -v`` / ``-version`` print a 26.10.0 version; any other call prints a hint on stderr
  and exits 1 (a test that launches Nextflow brings its own fake, first on ``PATH``).

``FAKE_PREREQS_ENV`` names the directory, so ``setup_guard`` lets only these two fakes start.
Tests that replace ``PATH`` wholesale call :func:`link_prereqs` on their own ``bin`` directory.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

FAKE_PREREQS_ENV = "GIPPA_TEST_FAKE_PREREQS"
NAMES = ("java", "nextflow")
JAVA_VERSION = "21.0.1"
NEXTFLOW_VERSION = "26.10.0"

FAKE_JAVA = f"""#!/bin/sh
echo 'openjdk version "{JAVA_VERSION}" 2023-10-17' >&2
echo 'OpenJDK Runtime Environment (build {JAVA_VERSION}+12-29)' >&2
echo 'OpenJDK 64-Bit Server VM (build {JAVA_VERSION}+12-29, mixed mode, sharing)' >&2
exit 0
"""

FAKE_NEXTFLOW = f"""#!/bin/sh
case "$1" in
-v) echo "nextflow version {NEXTFLOW_VERSION}.0"; exit 0;;
-version|--version)
    echo ""
    echo "      N E X T F L O W"
    echo "      version {NEXTFLOW_VERSION} build 0"
    exit 0;;
esac
echo "fake nextflow (tests/helpers/prereqs.py): called with '$*'; bring your own fake" >&2
exit 1
"""


def write_prereqs(bindir: Path) -> Path:
    """Write the fake ``java`` and ``nextflow`` into ``bindir`` (mode 0755)."""
    bindir.mkdir(parents=True, exist_ok=True)
    for name, text in (("java", FAKE_JAVA), ("nextflow", FAKE_NEXTFLOW)):
        exe = bindir / name
        exe.write_text(text, encoding="utf-8")
        exe.chmod(0o755)
    return bindir


def link_prereqs(bindir: Path, names: tuple[str, ...] = NAMES) -> None:
    """Copy the fakes of ``$GIPPA_TEST_FAKE_PREREQS`` (else fresh ones) into ``bindir``."""
    src = os.environ.get(FAKE_PREREQS_ENV)
    bindir.mkdir(parents=True, exist_ok=True)
    if not src:
        write_prereqs(bindir)
        return
    for name in names:
        shutil.copy2(Path(src) / name, bindir / name)


def is_fake_prereq(exe: object) -> bool:
    """``exe`` (a path or a bare name resolved on ``PATH``) is one of the conftest fakes."""
    src = os.environ.get(FAKE_PREREQS_ENV)
    if not src or not isinstance(exe, (str, os.PathLike)):
        return False
    path = str(exe)
    if os.sep not in path:
        path = shutil.which(path) or ""
    return bool(path) and Path(path).resolve().parent == Path(src).resolve()
