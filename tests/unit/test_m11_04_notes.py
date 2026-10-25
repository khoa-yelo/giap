"""``scripts/release/notes.py`` (M11-04 build item 2, acceptance 2, example 1, edge case "empty
notes section").

Interface read from the spec: ``python scripts/release/notes.py <version> [--draft]`` run from the
repository root reads ``CHANGELOG.md`` there and writes ``notes.md`` there; exit 0 on success, 1
when the ``## GIAP <version> (...)`` section is missing (without ``--draft``) or empty. Each test
runs the script in a throw-away copy (``tmp_path``) holding ``CHANGELOG.md``, ``pyproject.toml``
and ``changelog.d/``; the real files are never written and nothing touches the network.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
NOTES = REPO / "scripts" / "release" / "notes.py"

HEADER = """# Changelog

All notable changes are recorded here by towncrier from `changelog.d/` fragments.

<!-- towncrier release notes start -->

"""

CHANGELOG = (
    HEADER
    + """## GIAP 1.1.0 (2027-01-10)

### Features

- Alpha feature in the first section.

## GIAP 1.0.10 (2026-12-01)

### Fixes

- Beta fix that only 1.0.10 carries.

## GIAP 1.0.1 (2026-11-15)

### Fixes

- Gamma fix in the middle section.

## GIAP 1.0.0 (2026-11-01)

### Features

- Delta feature in the last section.
"""
)
BODIES = {
    "1.1.0": "Alpha feature in the first section.",
    "1.0.10": "Beta fix that only 1.0.10 carries.",
    "1.0.1": "Gamma fix in the middle section.",
    "1.0.0": "Delta feature in the last section.",
}
DRAFT_TEXT = "Draft-only fragment text for the rehearsal."


def _setup(root: Path, changelog: str = CHANGELOG, fragment: bool = True) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    (root / "CHANGELOG.md").write_text(changelog, encoding="utf-8")
    shutil.copy2(REPO / "pyproject.toml", root / "pyproject.toml")
    (root / "src" / "giap").mkdir(parents=True, exist_ok=True)
    shutil.copy2(REPO / "src" / "giap" / "__init__.py", root / "src" / "giap" / "__init__.py")
    frag = root / "changelog.d"
    frag.mkdir(exist_ok=True)
    (frag / ".gitkeep").write_text("", encoding="utf-8")
    if fragment:
        (frag / "+rehearsal.feature.md").write_text(DRAFT_TEXT + "\n", encoding="utf-8")
    return root


def _run(root: Path, *args: str) -> subprocess.CompletedProcess[str]:
    assert NOTES.is_file(), "scripts/release/notes.py does not exist (M11-04 build item 2)"
    env = dict(os.environ)
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env.get("PATH", "")
    return subprocess.run(
        [sys.executable, str(NOTES), *args],
        cwd=root,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
        env=env,
    )


def _notes(root: Path) -> str:
    p = root / "notes.md"
    assert p.is_file(), "notes.md not written"
    return p.read_text(encoding="utf-8")


def _first_line(text: str) -> str:
    return next(line for line in text.splitlines() if line.strip())


# --- section extraction -------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("version", "where"), [("1.1.0", "first"), ("1.0.1", "middle"), ("1.0.0", "last")]
)
def test_extracts_one_section(tmp_path: Path, version: str, where: str) -> None:
    root = _setup(tmp_path / "r")
    proc = _run(root, version)
    assert proc.returncode == 0, proc.stderr
    notes = _notes(root)
    assert BODIES[version] in notes, where
    for other, body in BODIES.items():
        if other != version:
            assert body not in notes, (where, other)
    assert "towncrier release notes start" not in notes
    assert "All notable changes" not in notes
    assert "DRAFT" not in notes


def test_version_match_is_exact(tmp_path: Path) -> None:
    """``1.0.1`` must not pick up the ``1.0.10`` section (nor the reverse)."""
    root = _setup(tmp_path / "r")
    assert _run(root, "1.0.1").returncode == 0
    assert BODIES["1.0.10"] not in _notes(root)
    assert _run(root, "1.0.10").returncode == 0
    n = _notes(root)
    assert BODIES["1.0.10"] in n and BODIES["1.0.1"] not in n


def test_subheadings_kept(tmp_path: Path) -> None:
    """The ``###`` type headings belong to the section; extraction stops at the next ``##``."""
    root = _setup(tmp_path / "r")
    assert _run(root, "1.1.0").returncode == 0
    n = _notes(root)
    assert "### Features" in n and "### Fixes" not in n


# --- missing and empty sections -----------------------------------------------------------------


def test_missing_section_exits_1(tmp_path: Path) -> None:
    root = _setup(tmp_path / "r")
    proc = _run(root, "2.0.0")
    assert proc.returncode == 1, (proc.stdout, proc.stderr)
    assert "2.0.0" in proc.stdout + proc.stderr


def test_missing_section_unbuilt_changelog_exits_1(tmp_path: Path) -> None:
    """Example 1's state (header and marker only) without ``--draft``: no fallback."""
    root = _setup(tmp_path / "r", changelog=HEADER)
    proc = _run(root, "0.0.1")
    assert proc.returncode == 1, (proc.stdout, proc.stderr)
    p = root / "notes.md"
    assert not p.exists() or DRAFT_TEXT not in p.read_text(encoding="utf-8")


@pytest.mark.parametrize("draft", [False, True])
def test_empty_section_exits_1(tmp_path: Path, draft: bool) -> None:
    """Edge case: a release with no notes is a mistake, also under ``--draft``."""
    log = HEADER + "## GIAP 2.0.0 (2027-02-01)\n\n\n## GIAP 1.0.0 (2026-11-01)\n\n- Old.\n"
    root = _setup(tmp_path / "r", changelog=log)
    proc = _run(root, "2.0.0", *(["--draft"] if draft else []))
    assert proc.returncode == 1, (proc.stdout, proc.stderr)


# --- draft fallback -----------------------------------------------------------------------------


def test_draft_fallback_uses_towncrier_draft(tmp_path: Path) -> None:
    """Example 1: no ``## GIAP 0.0.1`` section, ``--draft`` renders the towncrier draft, headed
    ``DRAFT``."""
    root = _setup(tmp_path / "r", changelog=HEADER)
    proc = _run(root, "0.0.1", "--draft")
    assert proc.returncode == 0, proc.stderr
    notes = _notes(root)
    assert "DRAFT" in _first_line(notes), notes
    assert DRAFT_TEXT in notes
    assert (root / "CHANGELOG.md").read_text(encoding="utf-8") == HEADER, "draft writes nothing"
    assert (root / "changelog.d" / "+rehearsal.feature.md").is_file(), "fragments kept"


def test_draft_flag_prefers_the_built_section(tmp_path: Path) -> None:
    """``--draft`` is a fallback only: an existing section wins and is not marked DRAFT."""
    root = _setup(tmp_path / "r")
    proc = _run(root, "1.0.0", "--draft")
    assert proc.returncode == 0, proc.stderr
    notes = _notes(root)
    assert BODIES["1.0.0"] in notes
    assert "DRAFT" not in notes and DRAFT_TEXT not in notes


# --- release candidates -------------------------------------------------------------------------


def test_rc_reads_base_section_with_prefix(tmp_path: Path) -> None:
    root = _setup(tmp_path / "r")
    proc = _run(root, "1.0.0rc1")
    assert proc.returncode == 0, proc.stderr
    notes = _notes(root)
    assert "Release candidate" in _first_line(notes), notes
    assert BODIES["1.0.0"] in notes
    assert all(b not in notes for v, b in BODIES.items() if v != "1.0.0")


def test_rc_without_base_section_exits_1(tmp_path: Path) -> None:
    root = _setup(tmp_path / "r")
    assert _run(root, "3.0.0rc2").returncode == 1


# --- property: every section of a generated changelog round-trips -------------------------------

_word = st.text(alphabet="abcdefghijklmnopqrstuvwxyz", min_size=3, max_size=10)


@settings(
    max_examples=10,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture, HealthCheck.too_slow],
)
@given(
    versions=st.lists(
        st.tuples(st.integers(0, 12), st.integers(0, 12), st.integers(0, 12)).map(
            lambda t: f"{t[0]}.{t[1]}.{t[2]}"
        ),
        min_size=1,
        max_size=5,
        unique=True,
    ),
    pick=st.integers(0, 4),
    words=st.lists(_word, min_size=5, max_size=5),
)
def test_property_extracts_exactly_the_picked_section(
    tmp_path_factory: pytest.TempPathFactory, versions: list[str], pick: int, words: list[str]
) -> None:
    pick %= len(versions)
    bodies = {v: f"Entry {i} {words[i]} for version {v} only." for i, v in enumerate(versions)}
    log = HEADER + "".join(
        f"## GIAP {v} (2026-11-{i + 1:02d})\n\n### Features\n\n- {bodies[v]}\n\n"
        for i, v in enumerate(versions)
    )
    root = _setup(tmp_path_factory.mktemp("prop"), changelog=log, fragment=False)
    target = versions[pick]
    proc = _run(root, target)
    assert proc.returncode == 0, proc.stderr
    notes = _notes(root)
    assert bodies[target] in notes
    assert [v for v, b in bodies.items() if v != target and b in notes] == []
