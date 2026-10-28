"""M9-16 behaviour 7, example 5 and the wording acceptance item: ``gippa.report.checks.wording``
finds "express" (any case, at a word start) outside quoted citations; every template, snippet,
report and plot string passes it (docs/specs/M9-16-report-content.md; GIPPA_SPEC.md §1 rule 6,
§10.3 "Wording"; M9-index decision 25).

Interface pinned here (M9-16 test author, beyond the spec text):
- ``wording(html: str) -> list[str]``: one entry per hit, in document order, each the text
  fragment (at least the hit word and the word after it, when there is one) around the hit;
  ``[]`` when clean;
- only text content is scanned; text inside ``<q>``, ``<blockquote>`` and ``<cite>`` (quoted
  citations) is skipped;
- the snapshot fixture (M9-25) does not exist yet, so the rendered-report scan uses a rich
  synthetic ``OUT/`` tree, and the plot-title scan reads every string literal (docstrings
  excepted) of ``src/gippa/viz`` and ``src/gippa/report`` (``checks.py``, which holds the pattern,
  excepted).
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest
from helpers.report_content import (
    det_rows,
    html,
    write_category_summary,
    write_det_model,
    write_diagnostics,
    write_genes,
    write_gsea,
    write_islands,
    write_manifest,
    write_missed_event,
    write_pairs,
    write_sensitivity,
    write_tier_fdr,
    write_variants,
)
from helpers.report_fixtures import run_out

ROOT = Path(__file__).resolve().parents[3]
REPORT = ROOT / "src" / "gippa" / "report"


def _wording(text: str) -> list[str]:
    from gippa.report.checks import wording

    return wording(text)


def test_example_5_expressed_genes() -> None:
    hits = _wording("<p>Top 20 expressed genes</p>")
    assert len(hits) == 1 and "expressed genes" in hits[0], hits


@pytest.mark.parametrize(
    "text",
    ["<h2>Expression levels</h2>", "<p>EXPRESSED</p>", "<td>co-expressed pairs</td>",
     "<p>genes (expressed)</p>", "<p>We express it.</p>"],
)  # fmt: skip
def test_any_case_word_start(text: str) -> None:
    assert _wording(text), text


@pytest.mark.parametrize(
    "text",
    ["<p>genes detected</p>", "<p>an inexpressible idea</p>", "<p>nonexpressed</p>", ""],
)
def test_clean_text(text: str) -> None:
    assert _wording(text) == [], text


@pytest.mark.parametrize(
    "text",
    ["<p>As <q>expressed genes</q> in the cited paper.</p>",
     "<blockquote>Highly expressed proteins</blockquote>",
     "<p><cite>Gene expression atlas</cite></p>"],
)  # fmt: skip
def test_quoted_citations_are_skipped(text: str) -> None:
    assert _wording(text) == [], text


def test_hits_outside_a_quote_still_count() -> None:
    hits = _wording("<p><q>expressed</q> and expression</p>")
    assert len(hits) == 1 and "expression" in hits[0], hits


def test_two_hits_in_order() -> None:
    hits = _wording("<p>expressed first</p><p>then expression second</p>")
    assert len(hits) == 2 and "first" in hits[0] and "second" in hits[1], hits


# ---------------------------------------------------------------- scans


def _templates() -> list[Path]:
    return sorted(p for p in (REPORT / "templates").rglob("*") if p.is_file()) + sorted(
        (REPORT / "snippets").glob("*.html")
    )


def test_snippet_directory_exists() -> None:
    assert sorted((REPORT / "snippets").glob("*.html")), "no snippets to scan"


@pytest.mark.parametrize("path", _templates(), ids=lambda p: p.name)
def test_templates_and_snippets_pass(path: Path) -> None:
    assert _wording(path.read_text(encoding="utf-8")) == [], path


def _rich(tmp_path: Path) -> Path:
    out = run_out(tmp_path)
    write_genes(out, 12, 3, n_2pep=8)
    write_tier_fdr(out, {"0.01": 0.0098})
    write_islands(out, [0.001, 0.3])
    write_diagnostics(out, [{"contig": "chr1", "scanned": True, "block_rule": "auto",
                             "l_hat": 4, "block_len": 9, "morans_i_lag1": 0.42}])  # fmt: skip
    write_pairs(out)
    write_missed_event(out)
    write_variants(out, [("V1", "high", "0.01")])
    write_category_summary(out, {"cog": 0.62, "kegg": 0.41})
    write_gsea(out, ["genomic_cluster"])
    write_det_model(out, det_rows("elastic_net", 0.21, 0.03, 0.001, "null_uncalibrated"))
    write_sensitivity(out, [("islands", "v1", 0, 0, 0, None)])
    write_manifest(out)
    return out


def test_rendered_report_passes(tmp_path: Path) -> None:
    assert _wording(html(_rich(tmp_path))) == []


def _literals(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    docs: set[int] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Module | ast.ClassDef | ast.FunctionDef | ast.AsyncFunctionDef):
            body = node.body
            if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant):
                docs.add(id(body[0].value))
    strings = [n for n in ast.walk(tree) if isinstance(n, ast.Constant)]
    return [n.value for n in strings if isinstance(n.value, str) and id(n) not in docs]


#: ``checks.py`` holds the pattern itself and is left out
# checklist.py (M9-24) states the "express" rule itself, as report/checks.py does
_PY = sorted(
    p for p in (ROOT / "src" / "gippa" / "viz").rglob("*.py") if p.name != "checklist.py"
) + sorted(p for p in REPORT.rglob("*.py") if p.name != "checks.py")


@pytest.mark.parametrize("path", _PY, ids=lambda p: str(p.relative_to(ROOT)))
def test_plot_titles_and_report_strings_pass(path: Path) -> None:
    from html import escape

    bad = [s for s in _literals(path) if _wording(f"<p>{escape(s)}</p>")]
    assert not bad, (path, bad)
