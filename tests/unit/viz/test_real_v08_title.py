"""K-12 plot survey, batch 1 item 7: V08's title and n/N labels.

K-12: "A (RNA processing and modification) has the highest median abundance (3.3 log2, 1 genes)".
Required: the "highest median" claim needs at least 5 measured genes in the category (else the
next category with n >= 5 is named, or no such claim is made); counts read "1 gene" / "2 genes";
the "n_detected/n_genes" labels do not overlap (rotated or thinned) at the K-12 scale (25
categories with 3-digit counts).
"""

from __future__ import annotations

import itertools
import re
from collections.abc import Sequence
from pathlib import Path

import pytest
from helpers.viz_out_tree import gene_row, make_out, write_schema_table
from matplotlib.backends.backend_agg import FigureCanvasAgg

from gippa.config import resolve_config
from gippa.io.run_view import RunView
from gippa.viz import base, style
from gippa.viz.base import Instance, Skip

MIN_N = 5
#: M9-06 as amended on 2026-10-24 (#322 review S5): the claim names its eligibility rule
AMONG = "among categories with ≥ 5 genes"
LETTERS = "JAKLBDVTMNUOCGEFHIPQRSWYZ"  # 25 COG letters


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


def g(tag: str, cog: str | None, log2: float | None, detected: bool = True) -> dict[str, object]:
    return gene_row(locus_tag=tag, feature_type="cds", cog_category=cog, abundance_log2=log2,
                    abundance_imputed=False, abundance_type="ibaq" if log2 is not None else "none",
                    detected=detected, detection_status="detected" if detected else "not_detected",
                    partial=False)  # fmt: skip


def _render(tmp: Path, genes: Sequence[dict[str, object]]) -> Instance:
    import gippa.viz.detection  # noqa: F401  (registers V08)

    out = make_out(tmp)
    write_schema_table(out / "tables" / "genes.tsv", "genes", genes)
    result = base.registered()["V08"].fn(RunView.from_outdir(out), resolve_config())
    assert not isinstance(result, Skip), result
    return result[0]


def _one_genes(text: str) -> list[str]:
    return re.findall(r"\b1 genes\b", text)


def test_a_single_gene_category_is_not_named_highest(tmp_path: Path) -> None:
    genes = [g("a1", "A", 30.0)]
    genes += [g(f"j{i}", "J", 20.0 + i) for i in range(6)]  # median 22.5, n = 6
    genes += [g(f"c{i}", "C", 18.0) for i in range(5)]  # n = 5
    title = _render(tmp_path, genes).title
    assert not title.startswith("A "), title
    assert title.startswith("J ("), title
    assert "highest median" in title
    assert "6 genes" in title, title
    assert AMONG in title, title


def test_exactly_five_measured_genes_may_carry_the_claim(tmp_path: Path) -> None:
    genes = [g(f"c{i}", "C", 25.0) for i in range(MIN_N)]
    genes += [g(f"j{i}", "J", 20.0) for i in range(8)]
    title = _render(tmp_path, genes).title
    assert title.startswith("C (") and "highest median" in title, title
    assert AMONG in title, title


def test_no_category_with_five_measured_genes_makes_no_highest_claim(tmp_path: Path) -> None:
    genes = [g("a1", "A", 30.0), *[g(f"j{i}", "J", 20.0 + i) for i in range(3)]]
    title = _render(tmp_path, genes).title
    assert "highest" not in title.lower(), title
    assert re.search(r"\d", title), title  # still states a number (checklist title item)
    assert _one_genes(title) == [], title


def test_singular_gene_without_a_cog_letter(tmp_path: Path) -> None:
    title = _render(tmp_path, [g("n1", None, 20.0)]).title
    assert _one_genes(title) == [], title
    assert "1 gene" in title, title


def test_base_plural_is_the_shared_helper() -> None:
    assert base.plural(1, "gene") == "1 gene"
    assert base.plural(6, "gene") == "6 genes"


def test_n_over_n_labels_do_not_overlap_at_k12_scale(tmp_path: Path) -> None:
    genes = []
    for k, letter in enumerate(LETTERS):
        n, det = 130 + k * 4, 92 + k * 3  # 3-digit labels such as "92/130"
        genes += [g(f"{letter}{i:03d}", letter, 20.0 + (i % 7) * 0.1) for i in range(det)]
        genes += [g(f"{letter}u{i:03d}", letter, None, detected=False) for i in range(n - det)]
    fig = _render(tmp_path, genes).figure
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    renderer = canvas.get_renderer()
    labels = [t for ax in fig.axes for t in ax.texts
              if t.get_visible() and re.fullmatch(r"\d+/\d+", t.get_text().strip())]  # fmt: skip
    assert len(labels) >= 2
    boxes = [(t.get_text(), t.get_window_extent(renderer)) for t in labels]
    bad = [(a, b) for (a, ea), (b, eb) in itertools.combinations(boxes, 2)
           if min(ea.x1, eb.x1) - max(ea.x0, eb.x0) > 0.5
           and min(ea.y1, eb.y1) - max(ea.y0, eb.y0) > 0.5]  # fmt: skip
    assert bad == [], bad[:5]


def test_the_claim_is_scoped_even_when_no_smaller_category_beats_it(tmp_path: Path) -> None:
    """The wording does not depend on whether a small category has a higher median."""
    genes = [g(f"j{i}", "J", 22.0 + i) for i in range(6)] + [g("c1", "C", 18.0)]
    title = _render(tmp_path, genes).title
    assert title.startswith("J (") and AMONG in title, title
