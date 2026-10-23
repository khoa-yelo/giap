"""``giap.report.build`` and the ``report`` step on tiny hand-built ``OUT/`` trees (M9-02
behaviour 2-6, examples 6 and 7, acceptance "cards omitted when their table is absent; every card
carries an error estimate or ``unestimated``"; test idea "two SVGs whose internal ids collide").
Structure checks only (html.parser); no ``--mpl`` baselines and no HTML snapshot (M9-25).

Interface assumed beyond the spec text (M9-02 test author; see tests/helpers/report_fixtures.py):
- ``SECTIONS`` entries have the attributes ``id``, ``heading``, ``plots`` (two-digit plot IDs, in
  display order) and ``tables`` (``giap.io.run_view.PATHS`` keys);
- ``build_report(view, {}, {}, cfg)`` is valid (no plot rendered, no skip recorded);
  ``skipped`` maps a plot ID to its M9-01 skip reason (``input_absent:pairs``);
- HTML hooks: ``<section id>``, ``data-plot="Vnn"`` on each listed plot's slot,
  ``class="card" data-card=<name>`` with ``.value`` and ``.estimate`` children,
  ``<table data-table=<key>>``;
- the ``report`` step's output ``report`` is a directory holding ``report.html``; its links point
  at ``../plots/...`` and ``../tables/...``.
"""

from __future__ import annotations

import re
from pathlib import Path
from types import ModuleType

import pytest
from helpers.report_fixtures import (
    INPUTS,
    SECTION_IDS,
    SECTION_PLOTS,
    Page,
    cfg,
    island_rows,
    local_refs,
    parse,
    register_fake,
    run_out,
    sha256,
    swap_registry,
    table_rows,
)
from helpers.viz_out_tree import gene_row, schema_columns, write_schema_table

import giap
from giap.io.run_view import RunView
from giap.steps import registry as reg

EXTERNAL = re.compile(r"^\s*(https?:|//|file:)", re.I)


@pytest.fixture
def fake_registry(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    return swap_registry(monkeypatch)


def _report(out: Path, *, skipped: dict[str, str] | None = None, **kw: object) -> Page:
    from giap.report.build import build_report

    view = RunView.from_outdir(out)
    return parse(build_report(view, {}, skipped or {}, cfg(*[f"{k}={v}" for k, v in kw.items()])))


def _write_genes(out: Path, detected: list[bool]) -> None:
    rows = [
        gene_row(
            locus_tag=f"g{i}",
            detected=d,
            detection_status="detected" if d else "not_detected",
            detection_tier="decoy_picked",
        )
        for i, d in enumerate(detected)
    ]
    write_schema_table(out / "tables" / "genes.tsv", "genes", rows)


def _write_tier_fdr(out: Path, fdr: float | None) -> None:
    rows = [
        {"section": "tier", "group_type": "q_threshold", "group": g, "metric": m, "value": v}
        for g in ("0.01", "0.05")
        for m, v in (("fdr", fdr), ("n_detected", 7))
    ]
    write_schema_table(out / "tables" / "detection_summary.tsv", "detection_summary", rows)


# ---------------------------------------------------------------- SECTIONS


def test_sections_ids_order_and_plots() -> None:
    from giap.report.build import SECTIONS

    assert tuple(s.id for s in SECTIONS) == SECTION_IDS
    for s in SECTIONS:
        assert tuple(s.plots) == SECTION_PLOTS[s.id], s.id
        assert s.heading.strip(), s.id
    by_id = {s.id: s for s in SECTIONS}
    assert {"islands", "regions_stats"} <= set(by_id["genome-context"].tables)
    assert "sensitivity" in set(by_id["sensitivity"].tables)


def test_section_tables_are_published_tables() -> None:
    from giap.io.run_view import PATHS
    from giap.report.build import SECTIONS

    for s in SECTIONS:
        for key in s.tables:
            assert PATHS[key].startswith("tables/") and PATHS[key].endswith(".tsv"), (s.id, key)


# ---------------------------------------------------------------- example 6, structure


def test_example_6_thirteen_sections_in_order(tmp_path: Path) -> None:
    page = _report(run_out(tmp_path))
    assert [s.attrs.get("id") for s in page.sections()] == list(SECTION_IDS)
    assert page.html.lstrip().lower().startswith("<!doctype html")


def test_every_section_exists_on_an_empty_run(tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    (out / "preflight").mkdir(parents=True)
    (out / "preflight" / "preflight.json").write_text('{"inputs": []}')
    page = _report(out)
    assert [s.attrs.get("id") for s in page.sections()] == list(SECTION_IDS)


def test_every_listed_plot_has_a_slot_in_its_section(tmp_path: Path) -> None:
    page = _report(run_out(tmp_path))
    for sid, pids in SECTION_PLOTS.items():
        section = page.section(sid)
        got = [n.attrs["data-plot"] for n in section.by_attr("data-plot")]
        for pid in pids:
            assert pid in got, (sid, pid, got)


def test_css_is_inlined_and_nothing_is_fetched(tmp_path: Path) -> None:
    page = _report(run_out(tmp_path, islands=2))
    styles = page.root.by_tag("style")
    assert styles and styles[0].text().strip(), "no inline <style>"
    assert not [n for n in page.root.by_tag("link") if "stylesheet" in n.attrs.get("rel", "")]
    assert not [n for n in page.root.by_tag("script") if n.attrs.get("src")]
    bad = [u for u in page.links() if EXTERNAL.match(u)]
    assert not bad, bad
    assert "@import" not in page.html


def test_no_absolute_path_no_timestamp_no_express(tmp_path: Path) -> None:
    page = _report(run_out(tmp_path, islands=2))
    assert str(tmp_path) not in page.html
    assert not re.search(r"\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}", page.html)
    assert "express" not in page.html.lower()  # §1 rule 6 ("detected", never "expressed")


def test_build_report_is_deterministic(tmp_path: Path) -> None:
    out = run_out(tmp_path, islands=3)
    assert _report(out).html == _report(out).html


# ---------------------------------------------------------------- tables (behaviour 4)


def test_table_shows_first_20_rows_as_written_with_a_link(tmp_path: Path) -> None:
    out = run_out(tmp_path, islands=25)
    page = _report(out)
    table = page.table("islands")
    assert table is not None, "no <table data-table='islands'>"
    header, body = table_rows(table)
    assert header == schema_columns("islands")
    assert len(body) == 20
    assert [r[0] for r in body] == [f"ISL{i + 1:04d}" for i in range(20)]
    # values as written in the TSV: 1-based starts, not RunView's 0-based frame
    assert body[0][header.index("start")] == "1001"
    assert body[0][header.index("p_value")] == "0.001"
    assert body[0][header.index("ci_low")] == "NA"
    section = page.section("genome-context")
    assert "../tables/islands.tsv" in [n.attrs.get("href") for n in section.by_tag("a")]


def test_absent_table_says_not_produced(tmp_path: Path) -> None:
    page = _report(run_out(tmp_path))
    assert page.table("sensitivity") is None
    assert "not produced in this run" in page.section("sensitivity").text()


def test_header_only_table_is_shown_with_no_rows(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_schema_table(out / "tables" / "sensitivity.tsv", "sensitivity")
    table = _report(out).table("sensitivity")
    assert table is not None
    header, body = table_rows(table)
    assert header == schema_columns("sensitivity") and body == []


# ---------------------------------------------------------------- skip notes (behaviour 2)


def test_skipped_plot_note_names_the_reason(tmp_path: Path) -> None:
    page = _report(
        run_out(tmp_path),
        skipped={"V11": "input_absent:pairs", "V14": "module_off:function"},
    )
    v11 = " ".join(page.slot("V11").text().split())
    assert "not run" in v11 and "pairs" in v11, v11
    v14 = " ".join(page.slot("V14").text().split())
    assert "not run" in v14 and "function" in v14 and "off" in v14, v14
    assert page.slot("V11") in page.section("operons").iter()


# ---------------------------------------------------------------- cards (behaviour 5)


def _card_text(page: Page, name: str, cls: str) -> str:
    card = page.cards()[name]
    found = card.by_class(cls)
    assert found, f"card {name}: no .{cls}"
    return " ".join(found[0].text().split())


def test_cards_omitted_when_their_table_is_absent(tmp_path: Path) -> None:
    page = _report(run_out(tmp_path))
    assert page.cards() == {}


def test_genes_detected_card_with_tier_fdr(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    _write_genes(out, [True] * 7 + [False] * 2)
    _write_tier_fdr(out, 0.0117)
    page = _report(out)
    cards = page.cards()
    assert "genes_detected" in cards
    assert cards["genes_detected"] in page.section("summary").iter()
    assert re.search(r"\b7\b", _card_text(page, "genes_detected", "value"))
    est = _card_text(page, "genes_detected", "estimate")
    assert est and "unestimated" not in est


def test_genes_detected_card_unestimated_without_tier_fdr(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    _write_genes(out, [True, True, False])
    _write_tier_fdr(out, None)
    assert "unestimated" in _card_text(_report(out), "genes_detected", "estimate")
    (out / "tables" / "detection_summary.tsv").unlink()
    assert "unestimated" in _card_text(_report(out), "genes_detected", "estimate")


def test_islands_card_counts_p_at_alpha_context_fwer(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    rows = island_rows(3, p_values=[0.001, 0.04, 0.2])
    write_schema_table(out / "tables" / "islands.tsv", "islands", rows)
    page = _report(out)
    assert re.search(r"\b2\b", _card_text(page, "islands", "value"))
    assert "FWER" in _card_text(page, "islands", "estimate")
    tight = _report(out, **{"stats.alpha.context": 0.01})
    assert re.search(r"\b1\b", _card_text(tight, "islands", "value"))


def test_operon_flags_card_counts_q_at_alpha(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    base = {"gene_a": "g1", "gene_b": "g2", "flag": "split"}
    rows = [
        {**base, "pair_id": "P1", "q_value": 0.05},
        {**base, "pair_id": "P2", "q_value": 0.5},
        {**base, "pair_id": "P3", "q_value": None},
    ]
    write_schema_table(out / "tables" / "operon_flags.tsv", "operon_flags", rows)
    page = _report(out)
    assert re.search(r"\b1\b", _card_text(page, "operon_flags", "value"))
    assert _card_text(page, "operon_flags", "estimate")


def test_novel_orfs_card_unestimated_class_fdr(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    rows = [
        {"orf_id": f"N{i}", "class": "intergenic", "called": c, "class_fdr": "unestimated"}
        for i, c in enumerate([True, True, False])
    ]
    write_schema_table(out / "tables" / "novel_orfs.tsv", "novel_orfs", rows)
    page = _report(out)
    assert re.search(r"\b2\b", _card_text(page, "novel_orfs", "value"))
    assert "unestimated" in _card_text(page, "novel_orfs", "estimate")


def test_every_card_has_a_value_and_an_estimate(tmp_path: Path) -> None:
    out = run_out(tmp_path, islands=2)
    _write_genes(out, [True, False])
    _write_tier_fdr(out, None)
    write_schema_table(out / "mapping" / "peptide_hits.tsv", "peptide_hits",
                       [{"peptide": "PEPTIDEK", "placement_q": 0.001}])  # fmt: skip
    write_schema_table(out / "mapping" / "mapping_fdr.tsv", "mapping_fdr")
    page = _report(out)
    cards = page.cards()
    assert {"genes_detected", "islands", "peptides_placed"} <= set(cards)
    for name in cards:
        assert _card_text(page, name, "value"), name
        assert _card_text(page, name, "estimate"), name


# ---------------------------------------------------------------- provenance (behaviour 6)


def test_provenance_version_inputs_and_config_digest(tmp_path: Path) -> None:
    from helpers.viz_out_tree import write_resolved_config

    out = run_out(tmp_path)
    resolved = write_resolved_config(out / "resolved_config.yaml", ["core"])
    text = _report(out).section("provenance").text()
    assert giap.__version__ in text
    for item in INPUTS:
        assert str(item["file"]) in text and str(item["sha256"]) in text
    assert sha256(resolved) in text


def test_provenance_without_resolved_config(tmp_path: Path) -> None:
    text = _report(run_out(tmp_path)).section("provenance").text()
    assert giap.__version__ in text
    assert str(INPUTS[0]["sha256"]) in text


# ---------------------------------------------------------------- the report step (SVG embedding)


def _steps(tmp_path: Path, c: object | None = None) -> Page:
    out = run_out(tmp_path)
    work = tmp_path / "work"
    pre = out / "preflight" / "preflight.json"
    c = c or cfg()
    reg.run_step("plots", c, {"preflight": pre}, {"plots": work / "plots"})  # type: ignore[arg-type]
    reg.run_step("report", c, {"plots": work / "plots", "preflight": pre},  # type: ignore[arg-type]
                 {"report": work / "report"})  # fmt: skip
    html = work / "report" / "report.html"
    assert html.is_file()
    assert [p.name for p in (work / "report").iterdir()] == ["report.html"]
    return parse(html.read_text(encoding="utf-8"))


def test_report_step_embeds_rendered_svgs_without_prolog(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    register_fake("V07")
    register_fake("V17")
    page = _steps(tmp_path)
    assert "<?xml" not in page.html and not page.pis
    for pid, sid in (("V07", "detection"), ("V17", "inputs-qc")):
        slot = page.slot(pid)
        assert slot in page.section(sid).iter()
        assert slot.by_tag("svg"), f"{pid}: no inline <svg>"
        assert "not run" not in slot.text()


def test_two_svgs_with_colliding_ids_get_unique_prefixed_ids(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    register_fake("V07")
    register_fake("V17")
    register_fake("V03", suffixes=("a", "b"))
    page = _steps(tmp_path)
    ids = page.ids()
    dup = sorted({i for i in ids if ids.count(i) > 1})
    assert not dup, dup
    for pid in ("V07", "V17", "V03"):
        svg_ids = [
            n.attrs["id"] for s in page.slot(pid).by_tag("svg") for n in s.iter() if "id" in n.attrs
        ]
        assert svg_ids, f"{pid}: the fake SVG has no ids (matplotlib clip paths expected)"
        assert all(i.startswith(pid) for i in svg_ids), (pid, svg_ids[:5])
    missing = sorted(set(local_refs(page)) - set(ids))
    assert not missing, f"references to ids not in the page: {missing[:5]}"


def test_plot_links_follow_the_formats(fake_registry: ModuleType, tmp_path: Path) -> None:
    register_fake("V07")
    page = _steps(tmp_path, cfg("report.formats=[svg, png, pdf]"))
    hrefs = [n.attrs.get("href") for n in page.slot("V07").by_tag("a")]
    for want in ("../plots/V07_obs_vs_exp.png", "../plots/V07_obs_vs_exp.pdf",
                 "../plots/data/V07_obs_vs_exp.tsv"):  # fmt: skip
        assert want in hrefs, hrefs


def test_svg_only_formats_link_no_png_or_pdf(fake_registry: ModuleType, tmp_path: Path) -> None:
    register_fake("V07")
    page = _steps(tmp_path)
    hrefs = [str(h) for h in (n.attrs.get("href") for n in page.root.by_tag("a")) if h]
    assert not [h for h in hrefs if h.endswith((".png", ".pdf"))], hrefs
    assert "../plots/data/V07_obs_vs_exp.tsv" in hrefs


def test_report_step_notes_a_skipped_plot(fake_registry: ModuleType, tmp_path: Path) -> None:
    register_fake("V11", needs=["pairs"])
    page = _steps(tmp_path)
    text = " ".join(page.slot("V11").text().split())
    assert "not run" in text and "pairs" in text, text
    assert not page.slot("V11").by_tag("svg")


@pytest.mark.parametrize("pid", ["V07"])
def test_example_7_report_step_links_resolve_in_work_layout(
    fake_registry: ModuleType, tmp_path: Path, pid: str
) -> None:
    register_fake(pid)
    page = _steps(tmp_path, cfg("report.formats=[svg, png]"))
    base = tmp_path / "work" / "report"
    for url in page.links():
        if url.startswith("#") or "/tables/" in url:
            continue  # tables are linked against OUT/, checked in test_cli_report.py
        assert not EXTERNAL.match(url), url
        assert (base / url).resolve().is_file(), url
