"""``CMP/report/compare.html`` review fixes for PR #254 (M9-23 lead decision 5, B5 page items).

- ``report.interactive=false``: no V16 twin and no Plotly.js anywhere in the page (§10.3; the
  inlined bundle is ~4.8 MB);
- the runs section gives each run directory relative to ``CMP/`` ("no absolute path", the
  ``build_compare`` contract), never the absolute path;
- the detection section shows at most ``TOP_INTERSECTIONS`` = 30 intersection rows, the 30
  largest (the table links the full TSV).

Same runs as ``test_compare_html.py``; never touches the real ``$GIAP_DB_DIR``; no network.
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner
from helpers import viz_compare as vc
from helpers.compare_data import STABILITY_COLUMNS, gene, write_run
from helpers.report_fixtures import Page, parse, table_rows

from giap.exitcodes import ExitCode

LABELS = ("A1", "A2", "A3")
N = 12


def _genes(shift: float, detected_upto: int) -> list[dict[str, object]]:
    vals = np.sin(np.arange(N) * 1.7) * 3 + np.arange(N) * 0.4
    return [gene(f"g{i:02d}", start=1000 * i, detected=i < detected_upto,
                 log2=float(v + shift * (i % 3)), tu=f"TU_{i // 2}")
            for i, v in enumerate(vals)]  # fmt: skip


def _run_compare(tmp: Path, name: str, *extra: str) -> Path:
    from giap.cli import cli

    (tmp / "db_empty").mkdir(exist_ok=True)
    runs = [tmp / f"run{lab}" for lab in LABELS]
    if not runs[0].exists():
        for lab, s, d in zip(LABELS, (0.0, 0.5, 1.0), (12, 10, 6), strict=True):
            write_run(tmp / f"run{lab}", _genes(s, d))
    out = tmp / name
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GIAP_DB_DIR", str(tmp / "db_empty"))
        mp.setenv("GIAP_ROOT", str(tmp / "db_empty"))
        res = CliRunner().invoke(cli, ["compare", *map(str, runs), "--outdir", str(out),
                                       "--labels", ",".join(LABELS),
                                       "--set", "stats.bootstrap=50", *extra])  # fmt: skip
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == ExitCode.SUCCESS, res.output
    return out


@pytest.fixture(scope="module")
def runs_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return tmp_path_factory.mktemp("compare_html_254")


@pytest.fixture(scope="module")
def off(runs_dir: Path) -> Iterator[Path]:
    yield _run_compare(runs_dir, "CMP_off", "--set", "report.interactive=false")


@pytest.fixture(scope="module")
def on(runs_dir: Path) -> Iterator[Path]:
    yield _run_compare(runs_dir, "CMP_on")


def _html(cmp: Path) -> str:
    return (cmp / "report" / "compare.html").read_text(encoding="utf-8")


# --- report.interactive off: no twin, no Plotly.js --------------------------------------------


def test_interactive_off_has_no_twin(off: Path) -> None:
    page = parse(_html(off))
    assert page.root.by_attr("data-plotly") == []


def test_interactive_off_has_no_plotly_js(off: Path) -> None:
    html = _html(off)
    assert re.search("(?i)plotly", html) is None, "Plotly.js or a twin in a static page"
    assert len(html.encode()) < 1_000_000, len(html.encode())


def test_interactive_off_still_inlines_v16(off: Path) -> None:
    page = parse(_html(off))
    assert page.section("detection").by_tag("svg")


def test_interactive_on_is_the_contrast(on: Path) -> None:
    """The same runs with the default ``report.interactive`` do carry the twin (so the off test
    is not vacuous)."""
    assert parse(_html(on)).section("abundance").by_attr("data-plotly")


# --- run directories are relative to CMP/ -----------------------------------------------------


def _runs_rows(page: Page) -> list[list[str]]:
    tables = page.section("runs").by_tag("table")
    assert len(tables) == 1
    header, rows = table_rows(tables[0])
    assert "directory" in header, header
    i = header.index("directory")
    return [[r[0], r[i]] for r in rows]


def test_run_directories_are_relative(off: Path, runs_dir: Path) -> None:
    rows = _runs_rows(parse(_html(off)))
    assert rows == [[lab, f"../run{lab}"] for lab in LABELS], rows


def test_page_has_no_absolute_path(off: Path, on: Path, runs_dir: Path) -> None:
    for cmp in (off, on):
        html = _html(cmp)
        for root in {str(runs_dir), str(runs_dir.resolve())}:
            assert root not in html, (cmp.name, root)


# --- at most 30 intersection rows -------------------------------------------------------------


def _build(genes: dict[str, pd.DataFrame]) -> tuple[Page, pd.DataFrame]:
    from giap.report.compare import RunInfo, build_compare

    det, ab = vc.inputs(genes, B=5)
    runs = [RunInfo(lab, f"../run{lab}", "a" * 64, "NA") for lab in genes]
    tables = {"compare_detection": det, "compare_abundance": ab,
              "compare_stability": pd.DataFrame(columns=list(STABILITY_COLUMNS))}  # fmt: skip
    return parse(build_compare(runs, None, tables)), vc.intersection_rows(det)


def test_at_most_30_intersection_rows_are_shown() -> None:
    page, inter = _build(vc.many_intersections(40, 6))
    assert len(inter) > 30
    table = page.table("compare_detection")
    assert table is not None
    header, rows = table_rows(table)
    assert len(rows) == 30, len(rows)
    i, n = header.index("intersection"), header.index("n_genes")
    assert [r[i] for r in rows] == [str(x) for x in inter["intersection"].head(30)]
    assert [int(r[n]) for r in rows] == list(range(40, 10, -1))


def test_fewer_than_30_intersections_are_all_shown() -> None:
    page, inter = _build(vc.many_intersections(12, 4))
    table = page.table("compare_detection")
    assert table is not None
    _, rows = table_rows(table)
    assert len(rows) == len(inter) <= 30


def test_the_full_intersection_table_is_linked() -> None:
    page, _ = _build(vc.many_intersections(40, 6))
    hrefs = [n.attrs.get("href") for n in page.section("detection").by_tag("a")]
    assert "../tables/compare_detection.tsv" in hrefs, hrefs
