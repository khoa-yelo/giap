"""The M9-15 twins in the report: every rendered plot with a registered twin gets one, built from
its data TSV; a skipped plot or a header-only data TSV gets none; ``report.interactive: false``
(``--no-interactive``) leaves no twin and no plotly.js; the page stays inside the size budget
(M9-15 degenerate rows, M9-14 behaviour 2, 3 and 5; GIAP_SPEC.md §10.3 "< 25 MB").

The real plot registry and the real twins are used (no swapped registries): a small ``OUT/``
with ``genes.tsv``, ``islands.tsv``, ``autocorrelation.tsv`` and a preflight JSON renders V2, V8,
V10, V17 and V18; V5, V6, V7 and V11 skip (their tables are absent). The ``plots`` and ``report``
steps run in-process (``helpers.interactive_fixtures.run_steps``). Never touches the real
``$GIAP_DB_DIR``; no network.

Interface assumed (M9-15 test author): the twin div of an instance without a suffix is
``plotly-<pid>``, its JSON sits in ``<script type="application/json" data-plotly=...>`` (M9-14).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from helpers import twins as tw
from helpers.interactive_fixtures import BANNER, page_of, run_steps, twin_divs
from helpers.report_fixtures import Page

from giap.io.plot_files import svg_title

pytestmark = pytest.mark.filterwarnings("ignore::DeprecationWarning")

RENDERED = {"V02", "V08", "V10", "V17", "V18"}
SKIPPED = {"V05", "V06", "V07", "V11"}
BUDGET = 25_000_000


@pytest.fixture(autouse=True)
def _no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIAP_DB_DIR", str(tmp_path / "db_empty"))


def _out(tmp: Path, *, genes: bool = True) -> Path:
    from helpers.viz_landscape import correlogram_rows, gene, island, landscape_out

    rows = []
    for i, (mid, log2) in enumerate(zip([10_000, 30_000, 50_000, 70_000, 90_000],
                                        [1.0, 2.0, 3.0, 2.0, 1.0], strict=True)):  # fmt: skip
        g = gene(f"g{i + 1}", mid, log2)
        g["product"] = f"protein {i + 1}"
        rows.append(g)
    rows.append(gene("g9", 60_000, None, status="not_detected"))
    return landscape_out(
        tmp,
        genes=rows if genes else None,
        islands=[island("ISL1", "high", 40_001, 50_000, 0.01)],
        autocorr=correlogram_rows(),
    )


def _specs(page: Page) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for s in page.root.by_tag("script"):
        if s.attrs.get("type") == "application/json" and "data-plotly" in s.attrs:
            out[s.attrs["data-plotly"]] = json.loads(s.text())
    return out


def _data_tsvs(plots: Path) -> dict[str, Path]:
    return {p.stem: p for p in sorted((plots / "data").glob("V*.tsv"))}


@pytest.fixture
def built(tmp_path: Path) -> tuple[Path, Page]:
    html = run_steps(tmp_path, out=_out(tmp_path))
    return html, page_of(html)


def test_every_rendered_plotly_plot_has_its_twin(built: tuple[Path, Page], tmp_path: Path) -> None:
    _, page = built
    tsvs = _data_tsvs(tmp_path / "work" / "plots")
    rendered = {stem[:3] for stem, p in tsvs.items()
                if len(tw.read(p.read_text(encoding="utf-8")))}  # fmt: skip
    assert rendered >= RENDERED, rendered
    want = {f"plotly-{pid}" for pid in rendered if pid in tw.M915}
    ids = {d.attrs["id"] for d in twin_divs(page.root)}
    assert ids == want
    assert set(_specs(page)) == want


def test_skipped_plots_have_no_twin(built: tuple[Path, Page]) -> None:
    _, page = built
    ids = {d.attrs["id"] for d in twin_divs(page.root)}
    assert {f"plotly-{p}" for p in RENDERED} <= ids, ids
    assert not {f"plotly-{p}" for p in SKIPPED} & ids


def test_page_twins_equal_figure_json_of_the_data_tsvs(
    built: tuple[Path, Page], tmp_path: Path
) -> None:
    _, page = built
    specs = _specs(page)
    assert {f"plotly-{p}" for p in RENDERED} <= set(specs), set(specs)
    for stem, path in _data_tsvs(tmp_path / "work" / "plots").items():
        pid = stem[:3]
        if f"plotly-{pid}" not in specs:
            continue
        # decision 7: the page title is the finding title stored in the static SVG
        svg = path.parents[1] / f"{stem}.svg"
        title = (svg_title(svg.read_text(encoding="utf-8")) if svg.exists() else None) or stem
        text = tw.figure_text(pid, tw.read(path.read_text(encoding="utf-8")), title)
        assert text is not None
        assert specs[f"plotly-{pid}"] == json.loads(text), pid


def test_v2_twin_in_the_page_shows_the_gene_product(built: tuple[Path, Page]) -> None:
    _, page = built
    hover = tw.all_hovers(_specs(page)["plotly-V02"])
    assert any("g2" in h and "protein 2" in h for h in hover)


def test_report_with_twins_is_inside_the_size_budget(built: tuple[Path, Page]) -> None:
    html, _ = built
    text = html.read_text(encoding="utf-8")
    assert text.count(BANNER) == 1
    assert html.stat().st_size < BUDGET


def test_no_interactive_leaves_no_twins(tmp_path: Path) -> None:
    html = run_steps(tmp_path, "report.interactive=false", out=_out(tmp_path))
    text = html.read_text(encoding="utf-8")
    assert "Plotly.newPlot" not in text and "data-plotly" not in text
    assert "plotly.js" not in text.lower()
    assert twin_divs(page_of(html).root) == []


def test_header_only_data_tsv_gets_no_twin(tmp_path: Path) -> None:
    """A plot in its empty state (header-only data TSV) shows only the static figure."""
    from helpers.viz_landscape import landscape_out

    assert tw.load().has_twin("V02") and tw.load().has_twin("V10")
    out = landscape_out(tmp_path, genes=[], autocorr=[])
    html = run_steps(tmp_path, out=out)
    tsvs = _data_tsvs(tmp_path / "work" / "plots")
    empty = {stem[:3] for stem, p in tsvs.items()
             if not len(tw.read(p.read_text(encoding="utf-8")))}  # fmt: skip
    assert empty & {"V02", "V08", "V10"}, empty
    ids = {d.attrs["id"] for d in twin_divs(page_of(html).root)}
    assert not {f"plotly-{p}" for p in empty} & ids
