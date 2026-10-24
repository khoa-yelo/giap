"""M9-16 behaviour 2, example 4: the report's methods snippets are generated from marked blocks
of ``docs/methods/<topic>.md`` and a test fails when a snippet and its source differ
(docs/specs/M9-16-report-content.md; GIAP_SPEC.md §10.3 section 12).

Interface pinned here (M9-16 test author, beyond the spec text):
- ``scripts/dev/sync_report_snippets.py [--check] [--docs DIR] [--out DIR]`` (defaults
  ``docs/methods`` and ``src/giap/report/snippets``): writes ``<topic>.html`` for every
  ``<topic>.md`` holding a ``<!-- report:begin -->`` … ``<!-- report:end -->`` block; ``--check``
  writes nothing and exits 1, naming the topic, when a snippet differs from its source, is
  missing, or has no source block (an orphan); exit 0 when all agree;
- the generated HTML: one ``<p>`` per blank-line-separated paragraph (lines joined by a space),
  one ``<ul>`` of ``<li>`` per run of ``- `` lines; text HTML-escaped; nothing else;
- the snippets are package data of ``giap.report`` (``importlib.resources``); the methods
  section renders every snippet; ``islands.html`` exists and holds the validity-domain sentence;
- no Markdown library is imported at run time.
"""

from __future__ import annotations

import ast
import shutil
import subprocess
import sys
from importlib import resources
from pathlib import Path

import pytest
from helpers.report_content import ISLANDS_VALIDITY, page, section_text, words
from helpers.report_fixtures import parse, run_out

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "dev" / "sync_report_snippets.py"
DOCS = ROOT / "docs" / "methods"
SNIPPETS = ROOT / "src" / "giap" / "report" / "snippets"
BEGIN, END = "<!-- report:begin -->", "<!-- report:end -->"


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True,
                          cwd=ROOT, check=False)  # fmt: skip


def _sourced_topics() -> list[str]:
    return sorted(p.stem for p in DOCS.glob("*.md") if BEGIN in p.read_text(encoding="utf-8"))


def test_repo_snippets_are_in_sync() -> None:
    res = _run("--check")
    assert res.returncode == 0, res.stdout + res.stderr


def test_every_marked_doc_has_a_snippet_and_back() -> None:
    topics = _sourced_topics()
    assert "islands" in topics, "docs/methods/islands.md has no report: block"
    assert sorted(p.stem for p in SNIPPETS.glob("*.html")) == topics


def _copy(tmp_path: Path) -> tuple[Path, Path]:
    docs, out = tmp_path / "docs", tmp_path / "snippets"
    shutil.copytree(DOCS, docs)
    shutil.copytree(SNIPPETS, out)
    return docs, out


def test_example_4_edited_snippet_fails_the_check(tmp_path: Path) -> None:
    docs, out = _copy(tmp_path)
    target = out / "islands.html"
    target.write_text(target.read_text(encoding="utf-8") + "<p>edited by hand</p>\n")
    res = _run("--check", "--docs", str(docs), "--out", str(out))
    assert res.returncode == 1, res.stdout + res.stderr
    assert "islands" in res.stdout + res.stderr


def test_edited_source_fails_the_check_until_regenerated(tmp_path: Path) -> None:
    docs, out = _copy(tmp_path)
    md = docs / "islands.md"
    md.write_text(md.read_text(encoding="utf-8").replace(END, "One more line.\n\n" + END, 1))
    assert _run("--check", "--docs", str(docs), "--out", str(out)).returncode == 1
    assert _run("--docs", str(docs), "--out", str(out)).returncode == 0
    assert _run("--check", "--docs", str(docs), "--out", str(out)).returncode == 0
    assert "One more line." in (out / "islands.html").read_text(encoding="utf-8")


def test_missing_and_orphan_snippets_fail_the_check(tmp_path: Path) -> None:
    docs, out = _copy(tmp_path)
    (out / "islands.html").unlink()
    assert _run("--check", "--docs", str(docs), "--out", str(out)).returncode == 1
    docs, out = _copy(tmp_path / "b")
    (out / "no_such_topic.html").write_text("<p>x</p>\n")
    assert _run("--check", "--docs", str(docs), "--out", str(out)).returncode == 1


def test_generated_html_paragraphs_lists_escaping(tmp_path: Path) -> None:
    docs, out = tmp_path / "docs", tmp_path / "snippets"
    docs.mkdir()
    out.mkdir()
    (docs / "toy.md").write_text(
        "# Toy\n\nOutside the block.\n\n"
        f"{BEGIN}\nFirst paragraph,\nsecond line.\n\n- item a\n- item b < c & d\n\n"
        f"Last paragraph.\n{END}\n\nAlso outside.\n",
        encoding="utf-8",
    )
    (docs / "plain.md").write_text("# No block here\n", encoding="utf-8")
    res = _run("--docs", str(docs), "--out", str(out))
    assert res.returncode == 0, res.stdout + res.stderr
    assert sorted(p.name for p in out.iterdir()) == ["toy.html"]
    text = (out / "toy.html").read_text(encoding="utf-8")
    assert "&lt;" in text and "&amp;" in text
    root = parse(text).root
    assert [words(n.text()) for n in root.by_tag("p")] == [
        "First paragraph, second line.",
        "Last paragraph.",
    ]
    assert [words(n.text()) for n in root.by_tag("li")] == ["item a", "item b < c & d"]
    assert len(root.by_tag("ul")) == 1
    assert "Outside" not in text and "outside" not in text
    assert "<!--" not in text


def test_snippets_are_package_data() -> None:
    files = resources.files("giap.report").joinpath("snippets")
    names = sorted(f.name for f in files.iterdir() if f.name.endswith(".html"))
    assert "islands.html" in names


def test_islands_snippet_holds_the_validity_domain_sentence() -> None:
    text = words(parse((SNIPPETS / "islands.html").read_text(encoding="utf-8")).root.text())
    assert ISLANDS_VALIDITY in text


def test_methods_section_renders_every_snippet(tmp_path: Path) -> None:
    methods = section_text(page(run_out(tmp_path)), "methods")
    for path in sorted(SNIPPETS.glob("*.html")):
        body = words(parse(path.read_text(encoding="utf-8")).root.text())
        assert body and body in methods, path.name


@pytest.mark.parametrize("path", sorted((ROOT / "src" / "giap" / "report").rglob("*.py")))
def test_no_markdown_library_at_run_time(path: Path) -> None:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods = [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            mods = [node.module or ""]
        else:
            continue
        assert not any(m.split(".")[0] in {"markdown", "mistune", "markdown_it"} for m in mods)
