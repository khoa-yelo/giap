"""``mkdocs.yml`` and its nav (M11-07 items 1 and 2), read through MkDocs' own config loader."""

# ruff: noqa: RUF002  (the U+2212 minus of example 5 is deliberate)
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pytest

mkdocs_config = pytest.importorskip("mkdocs.config", reason="mkdocs is not installed (extra docs)")

REPO = Path(__file__).resolve().parents[2]
CONFIG = REPO / "mkdocs.yml"
DOCS = REPO / "docs"


@pytest.fixture(scope="module")
def cfg() -> Any:
    assert CONFIG.is_file(), "mkdocs.yml does not exist (M11-07 item 1)"
    return mkdocs_config.load_config(config_file=str(CONFIG))


def _files(nav: Any) -> list[str]:
    """Every page path in a nav (sub)tree, in order."""
    if isinstance(nav, str):
        return [nav]
    if isinstance(nav, list):
        return [f for item in nav for f in _files(item)]
    if isinstance(nav, dict):
        return [f for v in nav.values() for f in _files(v)]
    return []


def _section(cfg: Any, title: str) -> Any:
    for item in cfg["nav"]:
        if isinstance(item, dict) and title in item:
            return item[title]
    raise AssertionError(f"nav has no section {title!r}")


def test_site_basics(cfg: Any) -> None:
    assert cfg["site_name"] == "GIPPA"
    assert Path(cfg["docs_dir"]).resolve() == DOCS.resolve()
    assert cfg["theme"].name == "material"
    assert cfg["theme"]["font"] is False  # no Google Fonts request


@pytest.mark.parametrize(
    ("path", "excluded"),
    [
        ("specs/M3-10-placement-fdr.md", True),
        ("specs/M11-index.md", True),
        ("backlog.md", True),
        ("adr/ADR-0001-import-spec-decisions.md", False),
        ("methods/placement-fdr.md", False),
        ("index.md", False),
    ],
)
def test_exclude_docs(cfg: Any, path: str, excluded: bool) -> None:
    spec = cfg["exclude_docs"]
    assert spec is not None, "mkdocs.yml has no exclude_docs (Q-L6)"
    assert bool(spec.match_file(path)) is excluded


@pytest.mark.parametrize("key", ["not_found", "anchors", "unrecognized_links"])
def test_link_validation_at_warn(cfg: Any, key: str) -> None:
    """Strict turns these warnings into errors; ``unrecognized_links`` is what catches example 5's
    ``[0, 1](1 − q)`` (MkDocs logs it at INFO by default)."""
    assert cfg["validation"]["links"][key] >= logging.WARNING


def test_markdown_extensions(cfg: Any) -> None:
    exts = set(cfg["markdown_extensions"])
    assert {"tables", "toc", "admonition", "pymdownx.superfences", "pymdownx.arithmatex"} <= exts
    mdx = cfg["mdx_configs"]
    assert mdx["toc"].get("permalink")
    assert mdx["pymdownx.arithmatex"].get("generic") is True


def test_mathjax3_from_jsdelivr_only(cfg: Any) -> None:
    scripts = [str(getattr(s, "path", s)) for s in cfg["extra_javascript"]]
    assert any("cdn.jsdelivr.net/npm/mathjax@3" in s for s in scripts), scripts
    remote = [s for s in scripts if "//" in s]
    assert all(s.startswith("https://cdn.jsdelivr.net/") for s in remote), remote
    css = [str(getattr(s, "path", s)) for s in cfg["extra_css"]]
    assert not [s for s in css if "//" in s]


# --- nav ----------------------------------------------------------------------------------------

TOP = [
    "Home",
    "Get started",
    "User guide",
    "Tutorials",
    "Reference",
    "Methods",
    "Decisions",
    "Developer guide",
    "Changelog",
    "Citation",
]


def test_nav_top_level_order(cfg: Any) -> None:
    titles = [next(iter(item)) for item in cfg["nav"] if isinstance(item, dict)]
    assert titles == TOP
    assert _section(cfg, "Home") == "index.md"


def test_every_nav_page_exists(cfg: Any) -> None:
    missing = [f for f in _files(cfg["nav"]) if "://" not in f and not (DOCS / f).is_file()]
    assert missing == []


@pytest.mark.parametrize(
    ("title", "pages"),
    [
        (
            "Get started",
            [
                "get-started/install.md",
                "get-started/upgrade.md",
                "get-started/setup-databases.md",
                "get-started/manual-databases.md",
                "get-started/quickstart.md",
            ],
        ),
        (
            "User guide",
            [
                f"user-guide/{n}.md"
                for n in [
                    "inputs",
                    "run",
                    "outputs",
                    "interpreting-results",
                    "plots-and-report",
                    "compare",
                    "makedb",
                    "convert",
                    "licences",
                    "platforms",
                    "runtimes",
                ]
            ],
        ),
        ("Tutorials", ["tutorials/ecoli-k12.md", "tutorials/nissle-pilot.md"]),
        (
            "Reference",
            [
                "reference/cli.md",
                "reference/configuration.md",
                "reference/tables.md",
                "reference/flags.md",
            ],
        ),
        ("Developer guide", ["developer.md", "release.md"]),
    ],
)
def test_section_pages(cfg: Any, title: str, pages: list[str]) -> None:
    assert _files(_section(cfg, title)) == pages


def test_methods_section_has_the_index_and_every_note(cfg: Any) -> None:
    files = _files(_section(cfg, "Methods"))
    assert files[0] == "methods/index.md"
    notes = sorted(
        f"methods/{p.name}" for p in (DOCS / "methods").glob("*.md") if p.name != "index.md"
    )
    assert len(notes) >= 29
    assert sorted(files[1:]) == notes


def test_decisions_section_has_the_index_and_every_adr(cfg: Any) -> None:
    files = _files(_section(cfg, "Decisions"))
    assert files[0] == "adr/index.md"
    adrs = sorted(f"adr/{p.name}" for p in (DOCS / "adr").glob("ADR-*.md"))
    assert sorted(files[1:]) == adrs
    assert "adr/TEMPLATE.md" not in files


@pytest.mark.parametrize("title", ["Changelog", "Citation"])
def test_single_page_sections(cfg: Any, title: str) -> None:
    files = _files(_section(cfg, title))
    assert len(files) == 1 and files[0].endswith(".md")


def test_no_specs_or_backlog_in_nav(cfg: Any) -> None:
    files = _files(cfg["nav"])
    assert not [f for f in files if f.startswith("specs/") or f == "backlog.md"]
