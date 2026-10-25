"""``mkdocs build --strict`` on this checkout (M11-07 acceptance 1 and 5, examples 5 and 6).

The site is built once per module into ``tmp_path`` (offline: no plugin here fetches anything).
"""

# ruff: noqa: RUF001  (the U+2212 minus of example 5 is deliberate)
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("mkdocs", reason="mkdocs is not installed (extra docs)")
pytest.importorskip("material", reason="mkdocs-material is not installed (extra docs)")

REPO = Path(__file__).resolve().parents[2]
CONFIG = REPO / "mkdocs.yml"


def _build(config: Path, site: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "mkdocs",
            "build",
            "--strict",
            "-f",
            str(config),
            "--site-dir",
            str(site),
        ],
        cwd=config.parent,
        capture_output=True,
        text=True,
        timeout=600,
    )


@pytest.fixture(scope="module")
def site(tmp_path_factory: pytest.TempPathFactory) -> Path:
    assert CONFIG.is_file(), "mkdocs.yml does not exist (M11-07 item 1)"
    out = tmp_path_factory.mktemp("site")
    proc = _build(CONFIG, out)
    assert proc.returncode == 0, proc.stdout[-4000:] + proc.stderr[-4000:]
    return out


def _html(site: Path) -> list[Path]:
    return sorted(site.rglob("*.html"))


def test_strict_build_passes(site: Path) -> None:
    assert (site / "index.html").is_file()
    assert (site / "reference" / "cli" / "index.html").is_file()


def test_example6_specs_and_backlog_not_built(site: Path) -> None:
    assert not (site / "specs").exists()
    assert not (site / "backlog").exists() and not (site / "backlog.html").exists()
    assert not list(site.rglob("M3-10-placement-fdr*"))


_REMOTE = re.compile(
    r"<(script|link)\b[^>]*?\b(src|href)\s*=\s*[\"']?((?:https?:)?//[^\"'\s>]+)", re.I
)


def test_no_third_party_host_but_jsdelivr(site: Path) -> None:
    """Acceptance 5: ``<script src>`` / ``<link href>`` reach only jsDelivr (MathJax)."""
    bad: set[tuple[str, str]] = set()
    for page in _html(site):
        text = page.read_text(encoding="utf-8", errors="replace")
        for m in _REMOTE.finditer(text):
            url = m.group(3)
            host = re.sub(r"^(?:https?:)?//", "", url).split("/", 1)[0].lower()
            if host != "cdn.jsdelivr.net":
                bad.add((page.relative_to(site).as_posix(), url))
    assert sorted(bad)[:20] == []


def test_no_google_fonts(site: Path) -> None:
    for page in _html(site):
        text = page.read_text(encoding="utf-8", errors="replace")
        assert "fonts.googleapis.com" not in text and "fonts.gstatic.com" not in text, page


def test_mathjax_loaded_and_math_rendered(site: Path) -> None:
    assert "cdn.jsdelivr.net/npm/mathjax@3" in (site / "index.html").read_text(encoding="utf-8")
    pages = [p.read_text(encoding="utf-8") for p in (site / "methods").rglob("index.html")]
    assert any('class="arithmatex"' in t for t in pages)


def test_changelog_page_embeds_changelog(site: Path) -> None:
    line = "All notable changes are recorded here by towncrier"
    assert line in (REPO / "CHANGELOG.md").read_text(encoding="utf-8")
    assert any(line in p.read_text(encoding="utf-8") for p in _html(site))


def test_citation_page_renders_citation_cff(site: Path) -> None:
    title = "GIAP: genome-integrated analysis of proteomics data"
    assert title in (REPO / "CITATION.cff").read_text(encoding="utf-8")
    pages = [p for p in _html(site) if "citation" in p.relative_to(site).as_posix().lower()]
    assert pages, "no citation page in the site"
    assert any(title in p.read_text(encoding="utf-8") for p in pages)


def test_methods_and_adr_pages_built(site: Path) -> None:
    for p in (REPO / "docs" / "methods").glob("*.md"):
        stem = "" if p.stem == "index" else p.stem
        assert (site / "methods" / stem / "index.html").is_file(), p.name
    for p in (REPO / "docs" / "adr").glob("ADR-*.md"):
        assert (site / "adr" / p.stem / "index.html").is_file(), p.name


# --- example 5: a link-shaped interval fails the strict build ------------------------------------


def _mini(tmp: Path, body: str) -> subprocess.CompletedProcess[str]:
    """A one-page site inheriting this repository's mkdocs.yml (its validation and extensions)."""
    (tmp / "d").mkdir()
    (tmp / "d" / "index.md").write_text(f"# Note\n\n{body}\n", encoding="utf-8")
    cfg = tmp / "mkdocs.yml"
    cfg.write_text(f"INHERIT: {CONFIG}\ndocs_dir: d\nnav:\n  - Note: index.md\n", encoding="utf-8")
    return _build(cfg, tmp / "site")


def test_example5_unescaped_interval_fails(tmp_path: Path) -> None:
    assert CONFIG.is_file(), "mkdocs.yml does not exist (M11-07 item 1)"
    proc = _mini(tmp_path, "The interval [0, 1](1 − q) holds.")
    assert proc.returncode != 0
    assert "1 − q" in proc.stdout + proc.stderr


@pytest.mark.parametrize(
    "body", ["The interval \\[0, 1\\](1 − q) holds.", "The interval `[0, 1](1 − q)` holds."]
)
def test_example5_escaped_interval_builds(tmp_path: Path, body: str) -> None:
    assert CONFIG.is_file(), "mkdocs.yml does not exist (M11-07 item 1)"
    proc = _mini(tmp_path, body)
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]
