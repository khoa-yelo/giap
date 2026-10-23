# ruff: noqa: RUF001, RUF003  (ρ is the spec's own character)
"""``CMP/report/compare.html`` and the V16 files written by ``giap compare`` (M9-23
"report/compare.html", example 4, acceptance criterion 2; GIAP_SPEC.md §7.9.3, §10.2 V16
"written only by ``giap compare``", §10.3 self-contained rules).

Structure checks with html.parser only (``helpers.report_fixtures.parse``), the M9-14 external-load
checker (``helpers.interactive_fixtures.external_loads``) and the run-report link rule
(tests/unit/test_cli_report.py): every non-fragment link is relative and names a file under
``CMP/``. Never touches the real ``$GIAP_DB_DIR``; no network.

Interface assumed where the spec leaves it open (M9-23 test author):

- ``giap compare`` writes, besides ``CMP/tables/``, ``CMP/plots/V16_compare.svg`` (``svg`` is in
  the default ``report.formats``), ``CMP/plots/data/V16_compare.tsv`` and
  ``CMP/report/compare.html``;
- the page's sections are exactly ``<section id="runs|detection|abundance|stability|
  provenance">`` in that order (the run report's ``<section id>`` convention, M9-02);
- tables are ``<table data-table="compare_detection|compare_abundance|compare_stability">`` in
  the ``detection``, ``abundance`` and ``stability`` sections; the abundance table's header holds
  the Bland-Altman columns;
- V16 has a Plotly twin (``giap.viz.interactive.has_twin("V16")`` after importing
  ``giap.viz.compare``); its figure JSON sits in a ``<script data-plotly=...>`` inside the
  ``abundance`` section (the M9-14 embedding);
- the ``runs`` section names every label and the genome SHA-256.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner
from helpers.compare_data import SHA_A, gene, write_run
from helpers.interactive_fixtures import EXTERNAL, external_loads
from helpers.report_fixtures import Page, parse, table_rows

from giap.exitcodes import ExitCode

SECTIONS = ["runs", "detection", "abundance", "stability", "provenance"]
LABELS = ("A1", "A2", "A3")
N = 12


def _genes(shift: float, detected_upto: int) -> list[dict[str, object]]:
    vals = np.sin(np.arange(N) * 1.7) * 3 + np.arange(N) * 0.4
    return [gene(f"g{i:02d}", start=1000 * i, detected=i < detected_upto,
                 log2=float(v + shift * (i % 3)), tu=f"TU_{i // 2}")
            for i, v in enumerate(vals)]  # fmt: skip


@pytest.fixture(scope="module")
def cmp(tmp_path_factory: pytest.TempPathFactory) -> Iterator[Path]:
    from giap.cli import cli

    tmp = tmp_path_factory.mktemp("compare_html")
    (tmp / "db_empty").mkdir()
    runs = [write_run(tmp / f"run{lab}", _genes(s, d))
            for lab, s, d in zip(LABELS, (0.0, 0.5, 1.0), (12, 10, 6), strict=True)]  # fmt: skip
    out = tmp / "CMP"
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GIAP_DB_DIR", str(tmp / "db_empty"))
        mp.setenv("GIAP_ROOT", str(tmp / "db_empty"))
        res = CliRunner().invoke(
            cli,
            [
                "compare",
                *map(str, runs),
                "--outdir",
                str(out),
                "--labels",
                ",".join(LABELS),
                "--set",
                "stats.bootstrap=50",
            ],
        )
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == ExitCode.SUCCESS, res.output
    yield out


@pytest.fixture(scope="module")
def page(cmp: Path) -> Page:
    path = cmp / "report" / "compare.html"
    assert path.is_file(), sorted(p.relative_to(cmp).as_posix() for p in cmp.rglob("*"))
    return parse(path.read_text(encoding="utf-8"))


def test_v16_files_are_written(cmp: Path) -> None:
    assert (cmp / "plots" / "V16_compare.svg").is_file()
    assert (cmp / "plots" / "data" / "V16_compare.tsv").is_file()


def test_v16_svg_carries_its_title(cmp: Path) -> None:
    svg = (cmp / "plots" / "V16_compare.svg").read_text(encoding="utf-8")
    assert "detected in all 3 runs; median pairwise ρ" in svg


def test_section_ids(page: Page) -> None:
    assert [s.attrs.get("id") for s in page.sections()] == SECTIONS


def test_example_4_no_external_loads(page: Page) -> None:
    assert external_loads(page.html) == []


def test_example_4_relative_links_resolve(cmp: Path, page: Page) -> None:
    links = page.links()
    assert "../plots/data/V16_compare.tsv" in links
    for name in ("compare_detection", "compare_abundance", "compare_stability"):
        assert f"../tables/{name}.tsv" in links, links
    for url in links:
        assert not EXTERNAL.match(url), url
        if url.startswith("#"):
            continue
        assert not url.startswith("/"), url
        assert (cmp / "report" / url).resolve().is_file(), url


def test_tables_sit_in_their_sections(page: Page) -> None:
    for sid, key in (("detection", "compare_detection"), ("abundance", "compare_abundance"),
                     ("stability", "compare_stability")):  # fmt: skip
        found = page.section(sid).find_all(
            lambda n, k=key: n.tag == "table" and n.attrs.get("data-table") == k
        )
        assert len(found) == 1, (sid, key)


def test_abundance_table_has_bland_altman_columns(page: Page) -> None:
    table = page.table("compare_abundance")
    assert table is not None
    header, rows = table_rows(table)
    for col in ("spearman_rho", "ci_low", "ci_high", "ba_bias_log2", "ba_loa_low", "ba_loa_high"):
        assert col in header, header
    assert len(rows) == 3


def test_detection_section_embeds_v16(page: Page) -> None:
    assert page.section("detection").by_tag("svg"), "no inline V16 SVG in the detection section"


def test_abundance_section_has_the_concordance_twin(page: Page) -> None:
    import giap.viz.compare  # noqa: F401  (registers V16 and its twin)
    from giap.viz.interactive import has_twin

    assert has_twin("V16")
    scripts = page.section("abundance").by_attr("data-plotly")
    assert scripts, "no Plotly figure JSON in the abundance section"


def test_runs_section_names_labels_and_genome(page: Page) -> None:
    text = page.section("runs").text()
    for lab in LABELS:
        assert lab in text
    assert SHA_A in text


def test_rerun_is_byte_identical(cmp: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Same runs, same seed, a sibling outdir: the same page and V16 data (§7.9.4 determinism)."""
    from giap.cli import cli

    monkeypatch.setenv("GIAP_DB_DIR", str(cmp.parent / "db_empty"))
    monkeypatch.setenv("GIAP_ROOT", str(cmp.parent / "db_empty"))
    runs = [cmp.parent / f"run{lab}" for lab in LABELS]
    again = cmp.parent / "CMP2"
    res = CliRunner().invoke(
        cli,
        [
            "compare",
            *map(str, runs),
            "--outdir",
            str(again),
            "--labels",
            ",".join(LABELS),
            "--set",
            "stats.bootstrap=50",
        ],
    )
    assert res.exit_code == ExitCode.SUCCESS, res.output
    for rel in ("report/compare.html", "plots/data/V16_compare.tsv"):
        assert (again / rel).read_bytes() == (cmp / rel).read_bytes(), rel
