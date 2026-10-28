"""PR #205 review round, the summary cards and unreadable tables (M9-02 behaviour 5; GIPPA_SPEC.md
§10.3 "Wording", R4; review B1, B2, B5 and the lead decisions of that round).

* B1: the novel-ORF card's FDR is over CALLED rows only, the max ``class_fdr`` among them; any
  called row ``unestimated`` makes the card ``unestimated``; uncalled ORFs never affect it.
* B2: the operon card counts rows with ``q ≤ alpha`` and says so: "flagged at q ≤ alpha
  (discordant or co-abundant)", never "flagged operon pairs ... FDR-controlled"; a concordant
  ``conflict_terminator_in_operon`` row (q > alpha by definition, §7.4.5) is not counted.
* B5: a schema-invalid or unreadable table never makes ``build_report`` (or the ``report``
  step) raise; the affected card / table preview carries a note naming the file, the rest of the
  page is still built.

HTML hooks as tests/helpers/report_fixtures.py; card texts are compared whitespace-normalised.
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path
from types import ModuleType

import pytest
from helpers.report_fixtures import Page, cfg, parse, run_out, swap_registry
from helpers.viz_out_tree import gene_row, write_schema_table
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.io.run_view import RunView
from gippa.io.tsv import format_value
from gippa.steps import registry as reg

NOVEL_CLASSES = ("nterm_extension", "cterm_extension", "intergenic", "sorf", "antisense")


# no surrogates: a lone "\ud800" cannot be UTF-8 encoded into the corrupt table
_CELL_CHARS = st.characters(blacklist_characters="\t\r\n", blacklist_categories=("Cs",))


@pytest.fixture
def fake_registry(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    return swap_registry(monkeypatch)


def _page(out: Path, *overrides: str) -> Page:
    from gippa.report.build import build_report

    return parse(build_report(RunView.from_outdir(out), {}, {}, cfg(*overrides)))


def _words(text: str) -> str:
    return " ".join(text.split())


def _card(page: Page, name: str) -> str:
    cards = page.cards()
    assert name in cards, f"no {name} card; cards: {sorted(cards)}"
    return _words(cards[name].text())


def _part(page: Page, name: str, cls: str) -> str:
    found = page.cards()[name].by_class(cls)
    assert found, f"card {name}: no .{cls}"
    return _words(found[0].text())


# ---------------------------------------------------------------- B1 novel-ORF card


def _novel(out: Path, rows: list[tuple[str, str, bool, str]]) -> None:
    body = [{"orf_id": i, "class": c, "called": k, "class_fdr": f} for i, c, k, f in rows]
    write_schema_table(out / "tables" / "novel_orfs.tsv", "novel_orfs", body)


def test_b1_uncalled_orf_in_another_class_does_not_lend_its_fdr(tmp_path: Path) -> None:
    """The reviewer's probe: two called ORFs with class FDR ``unestimated`` and one uncalled ORF
    of another class with class FDR 0.02 -> ``unestimated``, never "class FDR ≤ 0.02"."""
    out = run_out(tmp_path)
    _novel(out, [("N1", "intergenic", True, "unestimated"),
                 ("N2", "intergenic", True, "unestimated"),
                 ("N3", "nterm_extension", False, "0.02")])  # fmt: skip
    page = _page(out)
    assert re.search(r"\b2\b", _part(page, "novel_orfs", "value"))
    est = _part(page, "novel_orfs", "estimate")
    assert "unestimated" in est, est
    assert "0.02" not in _card(page, "novel_orfs"), _card(page, "novel_orfs")


def test_b1_max_class_fdr_over_called_rows_only(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    _novel(out, [("N1", "intergenic", True, "0.03"),
                 ("N2", "sorf", True, "0.08"),
                 ("N3", "antisense", False, "0.5"),
                 ("N4", "antisense", False, "0.5")])  # fmt: skip
    page = _page(out)
    assert re.search(r"\b2\b", _part(page, "novel_orfs", "value"))
    est = _part(page, "novel_orfs", "estimate")
    assert "0.08" in est, est
    assert "0.5" not in est and "unestimated" not in est, est


def test_b1_uncalled_unestimated_class_does_not_make_the_card_unestimated(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    _novel(out, [("N1", "intergenic", True, "0.03"),
                 ("N2", "sorf", False, "unestimated")])  # fmt: skip
    est = _part(_page(out), "novel_orfs", "estimate")
    assert "0.03" in est and "unestimated" not in est, est


def test_b1_one_called_unestimated_row_makes_the_card_unestimated(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    _novel(out, [("N1", "intergenic", True, "0.03"),
                 ("N2", "sorf", True, "unestimated")])  # fmt: skip
    est = _part(_page(out), "novel_orfs", "estimate")
    assert "unestimated" in est, est


_FDRS = st.sampled_from(["unestimated", "0.01", "0.02", "0.05", "0.125", "0.3"])


@st.composite
def _novel_rows(draw: st.DrawFn) -> list[tuple[str, str, bool, str]]:
    """Rows of a novel_orfs table; one class FDR per class (as the step writes it); ≥ 1 called."""
    fdr = {c: draw(_FDRS) for c in NOVEL_CLASSES}
    n = draw(st.integers(1, 8))
    rows = []
    for i in range(n):
        c = draw(st.sampled_from(NOVEL_CLASSES))
        rows.append((f"N{i + 1}", c, draw(st.booleans()), fdr[c]))
    if not any(k for _, _, k, _ in rows):
        oid, c0, _, f0 = rows[0]
        rows[0] = (oid, c0, True, f0)
    return rows


@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(rows=_novel_rows())
def test_b1_property_estimate_is_the_called_rows_bound(
    rows: list[tuple[str, str, bool, str]],
) -> None:
    """For any table with ≥ 1 called ORF: ``unestimated`` iff some called row is; else the max
    class FDR over called rows; dropping every uncalled row changes nothing on the card."""
    called = [f for _, _, k, f in rows if k]
    with tempfile.TemporaryDirectory() as d:
        full, only = Path(d) / "a", Path(d) / "b"
        out_full, out_only = run_out(full), run_out(only)
        _novel(out_full, rows)
        _novel(out_only, [r for r in rows if r[2]])
        est = _part(_page(out_full), "novel_orfs", "estimate")
        assert est == _part(_page(out_only), "novel_orfs", "estimate"), rows
        assert _part(_page(out_full), "novel_orfs", "value") == str(len(called))
    if "unestimated" in called:
        assert "unestimated" in est, (rows, est)
    else:
        bound = format_value(max(float(f) for f in called))
        assert bound in est and "unestimated" not in est, (rows, est)


# ---------------------------------------------------------------- B2 operon card

_PAIR = {"gene_a": "g1", "gene_b": "g2", "uniop_class": "OP", "promoter_between_yes": False,
         "terminator_between_yes": False, "lfdr": 0.5, "evidence_note": "x"}  # fmt: skip


def _flags(out: Path) -> None:
    rows = [
        # discordant, exactly at the default alpha = 0.1 (a boundary row: q ≤ alpha counts)
        {**_PAIR, "pair_id": "PAIR_1", "abundance_pattern": "discordant",
         "flag": "internal_promoter_supported", "promoter_between_yes": True, "q_value": 0.1},
        # co-abundant despite a terminator
        {**_PAIR, "pair_id": "PAIR_2", "abundance_pattern": "coabundant",
         "flag": "coabundant_despite_terminator", "terminator_between_yes": True,
         "q_value": 0.02},
        # concordant pairs flagged for a terminator inside an operon: q > alpha, never counted
        {**_PAIR, "pair_id": "PAIR_3", "abundance_pattern": "concordant",
         "flag": "conflict_terminator_in_operon", "terminator_between_yes": True,
         "q_value": 0.5},
        {**_PAIR, "pair_id": "PAIR_4", "abundance_pattern": "concordant",
         "flag": "conflict_terminator_in_operon", "terminator_between_yes": True,
         "q_value": None},
    ]  # fmt: skip
    write_schema_table(out / "tables" / "operon_flags.tsv", "operon_flags", rows)


def test_b2_operon_card_counts_q_at_alpha_and_says_so(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    _flags(out)
    page = _page(out)
    assert re.fullmatch(r"\D*2\D*", _part(page, "operon_flags", "value")), _card(
        page, "operon_flags"
    )
    text = _card(page, "operon_flags")
    assert "flagged at q ≤ 0.1 (discordant or co-abundant)" in text, text
    assert "flagged operon pairs" not in text, text


def test_b2_operon_card_label_follows_alpha(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    _flags(out)
    page = _page(out, "stats.alpha.operon_flags=0.05")
    assert re.fullmatch(r"\D*1\D*", _part(page, "operon_flags", "value"))
    text = _card(page, "operon_flags")
    assert "flagged at q ≤ 0.05 (discordant or co-abundant)" in text, text
    assert "0.1" not in text, text


def test_b2_only_concordant_conflicts_count_zero(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    rows = [{**_PAIR, "pair_id": f"PAIR_{i}", "abundance_pattern": "concordant",
             "flag": "conflict_terminator_in_operon", "terminator_between_yes": True,
             "q_value": q} for i, q in ((1, 0.4), (2, 0.9))]  # fmt: skip
    write_schema_table(out / "tables" / "operon_flags.tsv", "operon_flags", rows)
    page = _page(out)
    assert re.fullmatch(r"\D*0\D*", _part(page, "operon_flags", "value"))
    assert "flagged at q ≤ 0.1 (discordant or co-abundant)" in _card(page, "operon_flags")


# ---------------------------------------------------------------- B5 unreadable tables

GARBAGE = b"\xff\xfe\x00\x81binary\x00\x9c\xd0\n\x00\x00\xff"


def _good_genes(out: Path) -> None:
    rows = [gene_row(locus_tag=f"g{i}", detected=d, detection_status="detected" if d else
                     "not_detected", detection_tier="decoy_picked")
            for i, d in enumerate([True, False])]  # fmt: skip
    write_schema_table(out / "tables" / "genes.tsv", "genes", rows)


def _break(out: Path, case: str) -> str:
    """Write one broken table; returns its file name."""
    t = out / "tables"
    if case == "genes_missing_column":
        (t / "genes.tsv").write_text("locus_tag\tdetected\ng1\ttrue\n", encoding="utf-8")
        return "genes.tsv"
    if case == "genes_non_boolean_detected":
        _good_genes(out)
        text = (t / "genes.tsv").read_text(encoding="utf-8").splitlines()
        header = text[0].split("\t")
        col = header.index("detected")
        rows = [r.split("\t") for r in text[1:]]
        for r in rows:
            r[col] = "maybe"
        body = "\n".join(["\t".join(header), *("\t".join(r) for r in rows)]) + "\n"
        (t / "genes.tsv").write_text(body, encoding="utf-8")
        return "genes.tsv"
    if case == "novel_missing_class_fdr":
        (t / "novel_orfs.tsv").write_text("orf_id\tclass\tcalled\nN1\tsorf\ttrue\n")
        return "novel_orfs.tsv"
    if case == "operon_flags_missing_q":
        (t / "operon_flags.tsv").write_text("pair_id\tgene_a\tgene_b\tflag\nPAIR_1\ta\tb\tother\n")
        return "operon_flags.tsv"
    if case == "peptide_hits_missing_q":
        (out / "mapping").mkdir(exist_ok=True)
        (out / "mapping" / "peptide_hits.tsv").write_text("peptide\nPEPTIDEK\n")
        return "peptide_hits.tsv"
    if case == "detection_summary_garbage":
        _good_genes(out)
        (t / "detection_summary.tsv").write_bytes(GARBAGE)
        return "detection_summary.tsv"
    if case == "islands_binary":
        (t / "islands.tsv").write_bytes(GARBAGE)
        return "islands.tsv"
    raise AssertionError(case)


CASES = [
    "genes_missing_column",
    "genes_non_boolean_detected",
    "novel_missing_class_fdr",
    "operon_flags_missing_q",
    "peptide_hits_missing_q",
    "detection_summary_garbage",
    "islands_binary",
]


@pytest.mark.parametrize("case", CASES)
def test_b5_broken_table_is_a_note_not_a_crash(tmp_path: Path, case: str) -> None:
    out = run_out(tmp_path)
    if case != "islands_binary":
        from helpers.report_fixtures import island_rows

        write_schema_table(out / "tables" / "islands.tsv", "islands", island_rows(2))
    name = _break(out, case)
    page = _page(out)  # must not raise
    assert [s.attrs.get("id") for s in page.sections()][:1] == ["summary"]
    summary = _words(page.section("summary").text())
    assert name in summary, f"no note naming {name} in the summary: {summary}"
    if case != "islands_binary":
        assert "islands" in page.cards(), "a broken table removed an unrelated card"


def test_b5_unreadable_table_preview_is_a_note_not_not_produced(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    _break(out, "islands_binary")
    page = _page(out)
    ctx = _words(page.section("genome-context").text())
    assert "islands.tsv" in ctx, ctx
    assert "islands.tsv: not produced" not in ctx, "the file exists; it could not be read"
    assert page.table("islands") is None


@pytest.mark.parametrize("case", ["genes_missing_column", "islands_binary"])
def test_b5_report_step_writes_the_report_on_a_broken_table(
    fake_registry: ModuleType, tmp_path: Path, case: str
) -> None:
    out = run_out(tmp_path)
    name = _break(out, case)
    work = tmp_path / "work"
    pre = out / "preflight" / "preflight.json"
    tables = {"genes_tsv": out / "tables" / "genes.tsv", "islands": out / "tables" / "islands.tsv"}
    inputs = {"preflight": pre, **{k: p for k, p in tables.items() if p.is_file()}}
    c = cfg()
    reg.run_step("plots", c, dict(inputs), {"plots": work / "plots"})
    reg.run_step("report", c, {**inputs, "plots": work / "plots"}, {"report": work / "report"})
    html = (work / "report" / "report.html").read_text(encoding="utf-8")
    assert name in _words(parse(html).section("summary").text())


_CARD_TABLES = {
    "genes.tsv": "tables/genes.tsv",
    "detection_summary.tsv": "tables/detection_summary.tsv",
    "islands.tsv": "tables/islands.tsv",
    "operon_flags.tsv": "tables/operon_flags.tsv",
    "novel_orfs.tsv": "tables/novel_orfs.tsv",
    "peptide_hits.tsv": "mapping/peptide_hits.tsv",
    "unmapped.tsv": "mapping/unmapped.tsv",
}


@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    which=st.sampled_from(sorted(_CARD_TABLES)),
    body=st.one_of(
        st.binary(max_size=200),
        st.lists(
            st.lists(st.text(_CELL_CHARS, max_size=6),
                     min_size=1, max_size=4),
            max_size=4,
        ).map(lambda rows: "\n".join("\t".join(r) for r in rows).encode("utf-8")),
    ),
)  # fmt: skip
def test_b5_property_build_report_never_raises_on_a_corrupt_table(which: str, body: bytes) -> None:
    with tempfile.TemporaryDirectory() as d:
        out = run_out(Path(d), islands=0 if which == "islands.tsv" else 1)
        if which != "genes.tsv":
            _good_genes(out)
        path = out / _CARD_TABLES[which]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(body)
        page = _page(out)
        assert len(page.sections()) == 13
