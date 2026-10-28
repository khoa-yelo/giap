# ruff: noqa: RUF001, RUF002, RUF003  (the report's wording is Unicode: –, ℓ̂, −∞, σ70)
"""M9-16 summary cards: examples 1-3, every row of the card table, the degenerate rows and lead
decision 3 ("every number on a card is traceable to a named table and column")
(docs/specs/M9-16-report-content.md; GIPPA_SPEC.md §10.3 "Summary cards", "Wording"; R4).

Interface pinned here (M9-16 test author, beyond the spec text):
- card names (``data-card``): ``genes_detected``, ``peptides_placed``, ``proteins_placed``,
  ``islands``, ``operon_flags``, ``alt_calls``, ``novel_orfs``, ``model_r2`` (the M9-02 names
  kept; the Mapping row of the card table is the two cards ``peptides_placed`` and
  ``proteins_placed``);
- every card element carries ``data-source``: ``;``-joined ``<OUT-relative path>:<column>``
  entries naming each table column its number and its error estimate were read from (decision 3);
  the card dicts of ``gippa.report.build.cards`` carry the same list under ``"sources"``;
- the whole card text (whitespace-normalised) holds the example wording; the ``.estimate`` child
  holds the error estimate or the word ``unestimated``;
- ``alt_calls`` counts ``variants.tsv`` rows per ``confidence`` (``high``, ``candidate``), its
  estimate the ``class_error_rate`` range ``<min>–<max>`` (en dash), or ``unestimated`` when no row
  has a numeric rate;
- ``model_r2`` reads the ``elastic_net`` / ``all`` rows of ``tables/determinants_model.tsv``
  (ADR-0033 primary), never the ``lightgbm`` rows; "(uncalibrated)" follows the null p while its
  row carries ``null_uncalibrated`` (ADR-0034);
- ``proteins_placed``: the distinct ``protein_id`` of ``mapping/protein_hits.tsv`` of an input
  of those plus the ``unmapped.tsv`` rows of ``kind = protein``; ``peptides_placed`` counts its
  input as the distinct hit peptides plus the ``unmapped.tsv`` rows of ``kind = peptide`` only;
  the card text says "of <input>".
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path

import pytest
from helpers.report_content import (
    det_rows,
    page,
    words,
    write_det_model,
    write_genes,
    write_islands,
    write_mapping,
    write_novel,
    write_operon_flags,
    write_tier_fdr,
    write_variants,
)
from helpers.report_fixtures import Page, cfg, run_out
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.io.run_view import RunView

ALL_CARDS = {"genes_detected", "peptides_placed", "proteins_placed", "islands", "operon_flags",
             "alt_calls", "novel_orfs", "model_r2"}  # fmt: skip


def _card(p: Page, name: str) -> str:
    cards = p.cards()
    assert name in cards, f"no {name} card; cards: {sorted(cards)}"
    return words(cards[name].text())


def _part(p: Page, name: str, cls: str) -> str:
    found = p.cards()[name].by_class(cls)
    assert found, f"card {name}: no .{cls}"
    return words(found[0].text())


def _full_out(tmp: Path) -> Path:
    """Every card's table present."""
    out = run_out(tmp)
    write_genes(out, 6, 3, n_2pep=4)
    write_tier_fdr(out, {"0.01": 0.0098, "0.05": 0.041})
    write_mapping(out, placed_q=[0.001, 0.005, 0.2], unmapped_peptides=2,
                  proteins=["P1", "P2"], unmapped_proteins=1,
                  strata={"7-8": 0.0072, "9-12": 0.0031})  # fmt: skip
    write_islands(out, [0.001, 0.04, 0.3])
    write_operon_flags(out, [0.02, 0.5])
    write_variants(out, [("V1", "high", "0.012"), ("V2", "candidate", "0.034")])
    write_novel(out, [("N1", True, "0.03"), ("N2", False, "0.5")])
    write_det_model(out, det_rows("elastic_net", 0.21, 0.03, 0.001, "null_uncalibrated"))
    return out


# ---------------------------------------------------------------- example 1, 2: genes detected


def test_example_1_decoy_picked_812_detected(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_genes(out, 812, 40, tier="decoy_picked", n_2pep=700)
    write_tier_fdr(out, {"0.01": 0.0098, "0.05": 0.0411})
    text = _card(page(out), "genes_detected")
    assert re.search(r"(?<![\d.])812 genes detected\b", text), text
    assert "picked target–decoy FDR 0.0098 at q ≤ 0.01" in text, text
    assert "0.0411" not in text, text


def test_example_1_tier_fdr_follows_alpha_identification(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_genes(out, 812, 40, tier="decoy_picked")
    write_tier_fdr(out, {"0.01": 0.0098, "0.05": 0.0411})
    text = _card(page(out, "stats.alpha.identification=0.05"), "genes_detected")
    assert "FDR 0.0411 at q ≤ 0.05" in text, text


def test_example_2_two_peptide_tier_unestimated(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_genes(out, 640, 10, tier="two_peptide", n_2pep=640)
    p = page(out)
    text = _card(p, "genes_detected")
    assert re.search(r"(?<![\d.])640 genes detected\b", text), text
    assert "two-peptide rule" in text and "FDR unestimated" in text, text
    assert "unestimated" in _part(p, "genes_detected", "estimate")


def test_two_peptide_card_ignores_any_tier_fdr_row(tmp_path: Path) -> None:
    """The two-peptide rule has no FDR (§7.2): a stray tier row never lends it one."""
    out = run_out(tmp_path)
    write_genes(out, 640, 10, tier="two_peptide")
    write_tier_fdr(out, {"0.01": 0.0098})
    text = _card(page(out), "genes_detected")
    assert "0.0098" not in text and "FDR unestimated" in text, text


def test_genes_card_shows_detected_2pep_count_beside(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_genes(out, 812, 40, n_2pep=537)
    write_tier_fdr(out, {"0.01": 0.0098})
    assert re.search(r"(?<![\d.])537(?![\d.])", _card(page(out), "genes_detected"))


def test_mayu_tier_prints_its_fdr(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_genes(out, 300, 10, tier="mayu")
    write_tier_fdr(out, {"0.01": 0.0123})
    text = _card(page(out), "genes_detected")
    assert "FDR 0.0123 at q ≤ 0.01" in text and "unestimated" not in text, text


# ---------------------------------------------------------------- example 3: model R²


def test_example_3_model_r2(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    rows = det_rows("elastic_net", 0.21, 0.03, 0.001, "null_uncalibrated")
    rows += det_rows("lightgbm", 0.35, 0.05, 0.002, "null_uncalibrated")
    write_det_model(out, rows)
    p = page(out)
    text = _card(p, "model_r2")
    assert "Model R² 0.21 ± 0.03 (CV); null p 0.001 (uncalibrated)" in text, text
    assert "0.35" not in text and "0.002" not in text, "LightGBM is secondary (ADR-0033)"


def test_model_r2_without_null_uncalibrated_drops_the_word(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_det_model(out, det_rows("elastic_net", 0.21, 0.03, 0.001))
    text = _card(page(out), "model_r2")
    assert "0.21 ± 0.03" in text and "null p 0.001" in text, text
    assert "(uncalibrated)" not in text, text


def test_model_r2_unestimated_when_elastic_net_has_no_value(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_det_model(out, det_rows("elastic_net", None, None, None, "too_few_rows"))
    p = page(out)
    assert "unestimated" in _part(p, "model_r2", "estimate"), _card(p, "model_r2")


# ---------------------------------------------------------------- mapping


def test_peptides_and_proteins_placed_of_input(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_mapping(out, placed_q=[0.001, 0.005, 0.2], unmapped_peptides=2,
                  proteins=["P1", "P2", "P2"], unmapped_proteins=3,
                  strata={"7-8": 0.0072, "9-12": 0.0031})  # fmt: skip
    p = page(out)
    pep = _card(p, "peptides_placed")
    assert re.fullmatch(r"\D*2\D*", _part(p, "peptides_placed", "value")), pep
    assert re.search(r"\bof 5\b", pep), pep  # 3 hit peptides + 2 unmapped peptides
    prot = _card(p, "proteins_placed")
    assert re.fullmatch(r"\D*2\D*", _part(p, "proteins_placed", "value")), prot
    assert re.search(r"\bof 5\b", prot), prot  # 2 distinct proteins + 3 unmapped proteins


def test_peptide_card_lists_placement_fdr_per_stratum(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_mapping(out, placed_q=[0.001], unmapped_peptides=0, proteins=["P1"],
                  unmapped_proteins=0, strata={"7-8": 0.0072, "9-12": 0.0031})  # fmt: skip
    text = _card(page(out), "peptides_placed")
    assert "0.0072" in text and "0.0031" in text, text
    assert "0.4321" not in text, "de novo precision rows are not placement strata"
    assert "0.01" in text, text


# ---------------------------------------------------------------- islands, operons, novel ORFs


def test_islands_card_fwer_at_alpha_context(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_islands(out, [0.001, 0.05, 0.3])
    p = page(out)
    assert re.fullmatch(r"\D*2\D*", _part(p, "islands", "value"))
    assert "FWER ≤ 0.05" in _part(p, "islands", "estimate")


def test_operon_card_fdr_at_alpha(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_operon_flags(out, [0.02, 0.1, 0.5])
    p = page(out)
    assert re.fullmatch(r"\D*2\D*", _part(p, "operon_flags", "value"))
    assert "FDR ≤ 0.1" in _part(p, "operon_flags", "estimate")


# ---------------------------------------------------------------- ALT calls


def test_alt_calls_by_confidence_with_error_rate_range(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_variants(out, [("V1", "high", "0.012"), ("V2", "high", "0.012"),
                         ("V3", "high", "0.012"), ("V4", "candidate", "0.034"),
                         ("V5", "candidate", "0.034")])  # fmt: skip
    p = page(out)
    text = _card(p, "alt_calls")
    assert re.search(r"(?<![\d.])3 high\b", text), text
    assert re.search(r"(?<![\d.])2 candidate\b", text), text
    assert "0.012–0.034" in _part(p, "alt_calls", "estimate"), text


def test_alt_calls_unestimated_error_rate(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_variants(out, [("V1", "candidate", "unestimated")])
    assert "unestimated" in _part(page(out), "alt_calls", "estimate")


# ---------------------------------------------------------------- degenerate: absent tables


@pytest.mark.parametrize(
    ("card", "rel"),
    [
        ("alt_calls", "tables/variants.tsv"),
        ("novel_orfs", "tables/novel_orfs.tsv"),
        ("model_r2", "tables/determinants_model.tsv"),
        ("islands", "tables/islands.tsv"),
        ("operon_flags", "tables/operon_flags.tsv"),
        ("genes_detected", "tables/genes.tsv"),
    ],
)
def test_card_absent_not_zero_when_its_table_is_absent(tmp_path: Path, card: str, rel: str) -> None:
    out = _full_out(tmp_path)
    (out / rel).unlink()
    cards = page(out).cards()
    assert card not in cards, f"{card} shown without {rel}"
    assert ALL_CARDS - {card} <= set(cards), sorted(cards)


def test_full_run_has_every_card(tmp_path: Path) -> None:
    assert set(page(_full_out(tmp_path)).cards()) == ALL_CARDS


# ---------------------------------------------------------------- decision 3: traceability


def _header(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").splitlines()[0].split("\t")


def test_every_card_names_existing_tables_and_columns(tmp_path: Path) -> None:
    out = _full_out(tmp_path)
    for name, node in page(out).cards().items():
        src = node.attrs.get("data-source", "")
        entries = [e for e in src.split(";") if e]
        assert entries, f"card {name}: no data-source"
        for entry in entries:
            rel, sep, col = entry.rpartition(":")
            assert sep and rel and col, (name, entry)
            assert not rel.startswith("/"), (name, entry)
            assert (out / rel).is_file(), (name, entry)
            assert col in _header(out / rel), (name, entry)


@pytest.mark.parametrize(
    ("card", "must"),
    [
        ("genes_detected", {"tables/genes.tsv:detected", "tables/genes.tsv:detection_tier",
                            "tables/genes.tsv:detected_2pep",
                            "tables/detection_summary.tsv:value"}),
        ("peptides_placed", {"mapping/peptide_hits.tsv:placement_q",
                             "mapping/mapping_fdr.tsv:fdr_stratum"}),
        ("islands", {"tables/islands.tsv:p_value"}),
        ("operon_flags", {"tables/operon_flags.tsv:q_value"}),
        ("alt_calls", {"tables/variants.tsv:confidence", "tables/variants.tsv:class_error_rate"}),
        ("novel_orfs", {"tables/novel_orfs.tsv:called", "tables/novel_orfs.tsv:class_fdr"}),
        ("model_r2", {"tables/determinants_model.tsv:value"}),
    ],
)  # fmt: skip
def test_card_sources_name_the_columns_it_reads(tmp_path: Path, card: str, must: set[str]) -> None:
    node = page(_full_out(tmp_path)).cards()[card]
    got = set(node.attrs.get("data-source", "").split(";"))
    assert must <= got, (card, sorted(got))


def test_cards_function_returns_sources(tmp_path: Path) -> None:
    from gippa.report.build import cards

    out = _full_out(tmp_path)
    found, notes = cards(RunView.from_outdir(out), cfg())
    assert not notes, notes
    for c in found:
        assert c.get("sources"), c
        assert all(":" in s for s in c["sources"]), c


@settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(conf=st.lists(st.sampled_from(["high", "candidate"]), min_size=1, max_size=12))
def test_property_alt_counts_are_the_table_counts(conf: list[str]) -> None:
    """The card's per-confidence counts are exactly the ``confidence`` counts of the table."""
    with tempfile.TemporaryDirectory() as d:
        out = run_out(Path(d))
        write_variants(out, [(f"V{i}", c, "0.02") for i, c in enumerate(conf)])
        text = _card(page(out), "alt_calls")
    for level in ("high", "candidate"):
        n = conf.count(level)
        if n:
            assert re.search(rf"(?<![\d.]){n} {level}\b", text), (conf, text)


@settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(det=st.lists(st.booleans(), min_size=1, max_size=30))
def test_property_genes_value_is_the_detected_count(det: list[bool]) -> None:
    with tempfile.TemporaryDirectory() as d:
        out = run_out(Path(d))
        write_genes(out, sum(det), len(det) - sum(det))
        p = page(out)
        value = _part(p, "genes_detected", "value")
    assert re.fullmatch(rf"\D*{sum(det)}\D*", value), (det, value)


# ---------------------------------------------------------------- checks.unlabelled_numbers


def test_unlabelled_numbers_flags_a_card_without_estimate() -> None:
    from gippa.report.checks import unlabelled_numbers

    cards = [
        {"name": "a", "value": "3", "estimate": "FWER ≤ 0.05"},
        {"name": "b", "value": "7", "estimate": ""},
        {"name": "c", "value": "2", "estimate": "unestimated"},
    ]
    assert unlabelled_numbers(cards) == ["b"]


def test_unlabelled_numbers_missing_estimate_key() -> None:
    from gippa.report.checks import unlabelled_numbers

    assert unlabelled_numbers([{"name": "x", "value": "1"}]) == ["x"]
    assert unlabelled_numbers([]) == []


def test_unlabelled_numbers_clean_on_the_full_run(tmp_path: Path) -> None:
    from gippa.report.build import cards
    from gippa.report.checks import unlabelled_numbers

    found, _ = cards(RunView.from_outdir(_full_out(tmp_path)), cfg())
    assert len(found) == len(ALL_CARDS)
    assert unlabelled_numbers(found) == []
