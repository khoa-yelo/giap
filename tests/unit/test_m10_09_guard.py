"""Guard: no ``.ssl`` file is tracked and recorded Cascadia data hold only the provenance note.

Spec: docs/specs/M10-09-cascadia-convert.md, "What is committed", criterion 9; Khoa's rule of
2026-10-23 (row-level search output is never committed, not even a synthetic ``.ssl``: CI
fixtures are written into ``tmp_path``).

These tests read ``git ls-files`` and pass while nothing forbidden is tracked; they exist to fail
the commit that would add a ``.ssl`` file.
"""

from __future__ import annotations

import fnmatch
import subprocess
from pathlib import Path, PurePosixPath

REPO = Path(__file__).resolve().parents[2]
FORBIDDEN = ("*.ssl", "*.ssl.gz")
PROVENANCE = "tests/data/recorded/cascadia/0.0.7/PROVENANCE.md"


def _tracked(*paths: str) -> list[str]:
    res = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "-z", "--", *paths],
        capture_output=True,
        text=True,
        check=True,
    )
    return [p for p in res.stdout.split("\0") if p]


def _forbidden(name: str) -> bool:
    return any(fnmatch.fnmatchcase(name.lower(), pat) for pat in FORBIDDEN)


def test_no_ssl_file_is_tracked() -> None:
    bad = [p for p in _tracked() if _forbidden(PurePosixPath(p).name)]
    assert not bad, f"no .ssl file is ever committed (M10-09 criterion 9): {bad}"


def test_recorded_cascadia_holds_only_the_provenance_note() -> None:
    extra = [p for p in _tracked("tests/data/recorded/cascadia/") if p != PROVENANCE]
    assert not extra, f"only {PROVENANCE} may be tracked there: {extra}"


def test_the_patterns_catch_cascadia_outputs() -> None:
    """Self-check against the names Cascadia and the sample step write (W1, M10-08)."""
    for name in ("WT_1.ssl", "WT_1.ssl.gz", "a3.SSL", "x.Ssl.GZ"):
        assert _forbidden(name), name
    for name in ("WT_1.ssl.md", "PROVENANCE.md", "ssl.py"):
        assert not _forbidden(name), name
