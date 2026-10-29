"""M11-02: towncrier rendering rules and the retroactive milestone fragments.

Spec: docs/specs/M11-02-changelog-citation.md (Definition items 1-2, examples 1-2, edge cases,
acceptance criteria 1, 2, 4). towncrier runs from the venv on a copy of ``changelog.d/`` and
``pyproject.toml`` in ``tmp_path``; no network, the real ``CHANGELOG.md`` is never written.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
FRAGMENTS = ROOT / "changelog.d"
PYPROJECT = ROOT / "pyproject.toml"

REPO_URL = "https://github.com/khoa-yelo/gippa"
TYPES = [
    ("feature", "Features"),
    ("bugfix", "Fixes"),
    ("doc", "Documentation"),
    ("removal", "Removals and breaking changes"),
    ("misc", "Other"),
]
NAME_RE = re.compile(r"^(\+[a-z0-9-]+|\d+)\.(feature|bugfix|doc|removal|misc)\.md$")
MILESTONE_RE = re.compile(r"^\+m(\d{2})-[a-z0-9-]+\.feature\.md$")
TITLE = "## GIPPA 1.0.0 (2026-11-01)"
# Internal references that must not reach users: issue ids (M6-11, M10-15b), PR numbers (#317).
INTERNAL_REF_RE = re.compile(r"\bM\d{1,2}-\d{2}[a-z]?\b|#\d+")
# The 1.0.0 release-prep build (M11-13) consumes the fragments into CHANGELOG.md: the fragment
# checks then run on the built section instead (test_released_*).
RELEASED_RE = re.compile(r"^## GIPPA 1\.0\.0 \(\d{4}-\d{2}-\d{2}\)$", re.MULTILINE)
RELEASED = bool(RELEASED_RE.search((ROOT / "CHANGELOG.md").read_text(encoding="utf-8")))
pre_release = pytest.mark.skipif(RELEASED, reason="fragments consumed by the 1.0.0 build")
MILESTONE_MARKERS = [f"<!-- m{k:02d} -->" for k in range(6)] + [
    "<!-- m06a -->",
    "<!-- m06b -->",
    *(f"<!-- m{k:02d} -->" for k in range(7, 11)),
]


def _towncrier_cfg() -> dict:
    return tomllib.loads(PYPROJECT.read_text(encoding="utf-8"))["tool"]["towncrier"]


def _fragment_files() -> list[Path]:
    return sorted(p for p in FRAGMENTS.iterdir() if p.is_file() and p.name != ".gitkeep")


def _milestone_fragments() -> list[Path]:
    return [p for p in _fragment_files() if MILESTONE_RE.match(p.name)]


def _draft(tmp_path: Path, extra: dict[str, str] | None = None) -> str:
    """Render ``towncrier build --draft`` on a copy of the repo's changelog inputs."""
    shutil.copytree(FRAGMENTS, tmp_path / "changelog.d")
    shutil.copy2(PYPROJECT, tmp_path / "pyproject.toml")
    shutil.copy2(ROOT / "CHANGELOG.md", tmp_path / "CHANGELOG.md")
    (tmp_path / "src" / "gippa").mkdir(parents=True)
    shutil.copy2(ROOT / "src" / "gippa" / "__init__.py", tmp_path / "src" / "gippa" / "__init__.py")
    for name, text in (extra or {}).items():
        (tmp_path / "changelog.d" / name).write_text(text, encoding="utf-8")
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "towncrier",
            "build",
            "--draft",
            "--version",
            "1.0.0",
            "--date",
            "2026-11-01",
            "--dir",
            str(tmp_path),
            "--config",
            str(tmp_path / "pyproject.toml"),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr + proc.stdout
    out = proc.stdout
    start = out.find(TITLE)
    assert start >= 0, f"draft has no {TITLE!r} heading:\n{out}"
    return out[start:]


def _sections(draft: str) -> dict[str, list[str]]:
    """Map each ``### <type>`` heading to its top-level bullets (in order)."""
    sections: dict[str, list[str]] = {}
    current: str | None = None
    for line in draft.splitlines():
        if line.startswith("### "):
            current = line[4:].strip()
            sections[current] = []
        elif line.startswith("- ") and current is not None:
            sections[current].append(line)
    return sections


# ------------------------------------------------------------------ [tool.towncrier] (item 1)


def test_towncrier_name_and_title_format() -> None:
    cfg = _towncrier_cfg()
    assert cfg.get("name") == "GIPPA"
    assert cfg.get("title_format") == "## GIPPA {version} ({project_date})"
    assert cfg.get("directory") == "changelog.d"
    assert cfg.get("filename") == "CHANGELOG.md"


def test_towncrier_issue_format_links_pull_requests() -> None:
    assert _towncrier_cfg().get("issue_format") == f"[#{{issue}}]({REPO_URL}/pull/{{issue}})"


def test_towncrier_types_in_order_with_content() -> None:
    types = _towncrier_cfg().get("type")
    assert isinstance(types, list), "[[tool.towncrier.type]] tables are missing"
    assert [(t["directory"], t["name"]) for t in types] == TYPES
    assert all(t.get("showcontent", True) is True for t in types)


def test_towncrier_orphan_prefix_is_plus() -> None:
    # towncrier's default orphan prefix is "+"; the milestone fragment names rely on it.
    assert _towncrier_cfg().get("orphan_prefix", "+") == "+"


# ------------------------------------------------------------------ fragment names (edge cases)


def test_fragment_directory_exists() -> None:
    assert FRAGMENTS.is_dir()


@pre_release
def test_fragment_directory_is_nonempty() -> None:
    assert _fragment_files()


def test_every_fragment_name_matches_the_rule() -> None:
    bad = [p.name for p in _fragment_files() if not NAME_RE.match(p.name)]
    assert bad == [], f"fragments towncrier would ignore or misprint: {bad}"


@pre_release
def test_old_fragment_names_are_renamed() -> None:
    names = {p.name for p in _fragment_files()}
    assert "m0-scaffold.feature.md" not in names
    assert "m6-11-regulondb.feature.md" not in names
    assert "+m00-scaffold.feature.md" in names
    assert "+m06-regulondb.feature.md" in names


def test_no_fragment_is_empty() -> None:
    empty = [p.name for p in _fragment_files() if not p.read_text(encoding="utf-8").strip()]
    assert empty == [], f"empty fragments render a blank bullet: {empty}"


def test_name_rule_rejects_known_bad_names() -> None:
    # self-test of the rule the directory check applies
    for good in ("+m00-scaffold.feature.md", "331.bugfix.md", "+known-limitations.misc.md"):
        assert NAME_RE.match(good), good
    for bad in (
        "x.feat.md",
        "m0-scaffold.feature.md",
        "+M01-x.feature.md",
        "331.bugfix.txt",
        "README.md",
    ):
        assert not NAME_RE.match(bad), bad


# ------------------------------------------------------------------ milestone fragments (item 2)


@pre_release
@pytest.mark.parametrize("k", range(11))
def test_every_milestone_has_a_feature_fragment(k: int) -> None:
    hits = [p.name for p in _milestone_fragments() if int(MILESTONE_RE.match(p.name).group(1)) == k]  # type: ignore[union-attr]
    assert hits, f"no +m{k:02d}-*.feature.md fragment for milestone M{k}"


@pre_release
def test_milestone_six_has_operons_and_regulondb() -> None:
    m06 = [p.name for p in _milestone_fragments() if p.name.startswith("+m06-")]
    assert len(m06) == 2, m06
    assert "+m06-regulondb.feature.md" in m06


@pre_release
def test_twelve_milestone_fragments() -> None:
    assert len(_milestone_fragments()) == 12, [p.name for p in _milestone_fragments()]


@pre_release
@pytest.mark.parametrize("k", range(1, 11))
def test_new_milestone_fragments_are_three_to_eight_lines(k: int) -> None:
    for p in _milestone_fragments():
        if p.name.startswith(f"+m{k:02d}-"):
            lines = [ln for ln in p.read_text(encoding="utf-8").splitlines() if ln.strip()]
            assert 3 <= len(lines) <= 8, f"{p.name}: {len(lines)} lines"


def test_orphan_fragments_carry_no_internal_issue_numbers() -> None:
    offenders = {}
    for p in _fragment_files():
        if p.name.startswith("+"):
            hits = INTERNAL_REF_RE.findall(p.read_text(encoding="utf-8"))
            if hits:
                offenders[p.name] = hits
    assert offenders == {}


@pre_release
def test_known_limitations_fragment_lists_section_16_deferrals() -> None:
    p = FRAGMENTS / "+known-limitations.misc.md"
    assert p.is_file()
    text = p.read_text(encoding="utf-8")
    for needle in ("sample", "V4", "V9", "RO-Crate", "validation"):
        assert needle in text, f"known limitations does not mention {needle!r}"


# ------------------------------------------------------------------ example 1: the 1.0 draft


@pre_release
def test_example1_draft_title_and_feature_level(tmp_path: Path) -> None:
    draft = _draft(tmp_path)
    lines = [ln for ln in draft.splitlines() if ln.strip()]
    assert lines[0] == TITLE
    # towncrier renders the type headings one level below a "## " title_format
    assert lines[1] == "### Features"


@pre_release
def test_example1_twelve_feature_bullets_in_milestone_order(tmp_path: Path) -> None:
    draft = _draft(tmp_path)
    bullets = _sections(draft)["Features"]
    feature_files = [p for p in _fragment_files() if p.name.endswith(".feature.md")]
    assert len(bullets) == len(feature_files)
    milestone = _milestone_fragments()
    assert len(milestone) == 12
    # each milestone fragment's first line appears exactly once, in file-name (= M0 -> M10) order
    positions = []
    for p in milestone:
        first = p.read_text(encoding="utf-8").strip().splitlines()[0].strip()
        assert draft.count(first) == 1, f"{p.name}: first line not unique in the draft"
        positions.append(draft.index(first))
    order = [p.name for _, p in sorted(zip(positions, milestone, strict=True))]
    assert order == [p.name for p in milestone], "Features bullets are not in M0 -> M10 order"


def test_example1_no_fragment_name_suffix(tmp_path: Path) -> None:
    draft = _draft(tmp_path)
    assert "(m0-scaffold)" not in draft
    assert "(m6-11-regulondb)" not in draft
    assert not re.search(r"\(\+?m\d+[a-z0-9-]*\)\s*$", draft, flags=re.MULTILINE)
    assert "known-limitations" not in draft


@pre_release
def test_example1_other_section_has_known_limitations(tmp_path: Path) -> None:
    draft = _draft(tmp_path)
    sections = _sections(draft)
    assert next(iter(sections)) == "Features"
    assert "Other" in sections
    assert list(sections).index("Other") > list(sections).index("Features")
    first = (
        (FRAGMENTS / "+known-limitations.misc.md")
        .read_text(encoding="utf-8")
        .strip()
        .splitlines()[0]
    )
    other_block = draft[draft.index("### Other") :]
    assert first.strip() in other_block


# ------------------------------------------------------------------ example 2: a PR fragment


def test_example2_pr_fragment_renders_a_pull_link(tmp_path: Path) -> None:
    # one fragment of each neighbouring type, so the order check holds once changelog.d is empty
    extra = {"331.bugfix.md": "Fix X.\n", "+x.feature.md": "Add X.\n", "+x.misc.md": "Note X.\n"}
    draft = _draft(tmp_path, extra)
    sections = _sections(draft)
    assert "Fixes" in sections
    # the real changelog.d may hold other bugfix fragments; this one renders exactly once
    mine = [ln for ln in sections["Fixes"] if "#331" in ln]
    assert mine == [f"- Fix X. ([#331]({REPO_URL}/pull/331))"]
    names = list(sections)
    assert names.index("Features") < names.index("Fixes") < names.index("Other")


def test_draft_does_not_write_changelog(tmp_path: Path) -> None:
    before = (ROOT / "CHANGELOG.md").read_bytes()
    _draft(tmp_path)
    assert (ROOT / "CHANGELOG.md").read_bytes() == before
    assert (tmp_path / "CHANGELOG.md").read_bytes() == before


# ------------------------------------------------------------------ after the 1.0.0 build


def _released_section() -> str:
    text = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    m = RELEASED_RE.search(text)
    assert m is not None
    end = text.find("\n## ", m.end())
    return text[m.start() : end if end >= 0 else len(text)]


def _block(section: str, heading: str) -> str:
    start = section.index(f"### {heading}")
    end = section.find("\n### ", start + 1)
    return section[start : end if end >= 0 else len(section)]


@pytest.mark.skipif(not RELEASED, reason="before the 1.0.0 build")
def test_released_section_title_and_type_order() -> None:
    section = _released_section()
    lines = [ln for ln in section.splitlines() if ln.strip()]
    assert lines[1] == "### Features"
    names = list(_sections(section))
    assert names[0] == "Features"
    assert names.index("Features") < names.index("Fixes") < names.index("Other")


@pytest.mark.skipif(not RELEASED, reason="before the 1.0.0 build")
def test_released_milestone_bullets_in_order() -> None:
    features = _block(_released_section(), "Features")
    for marker in MILESTONE_MARKERS:
        assert features.count(f"- {marker}") == 1, marker
    positions = [features.index(m) for m in MILESTONE_MARKERS]
    assert positions == sorted(positions), "milestone bullets are not in M0 -> M10 order"


@pytest.mark.skipif(not RELEASED, reason="before the 1.0.0 build")
def test_released_other_section_has_known_limitations() -> None:
    section = _released_section()
    other = _block(section, "Other")
    assert "- Known limitations of version 1:" in other
    for needle in ("sample", "V4", "V9", "RO-Crate", "validation"):
        assert needle in other, f"known limitations does not mention {needle!r}"
    assert "known-limitations" not in section
    assert not re.search(r"\(\+?m\d+[a-z0-9-]*\)\s*$", section, flags=re.MULTILINE)
