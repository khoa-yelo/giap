# ruff: noqa: RUF001, RUF002, RUF003  (the report's wording is Unicode: –, ℓ̂, −∞, σ70)
"""M9-16 behaviour 3 and 4: every mandatory caveat is present in its section when its condition
holds; the operon bin-count table; the ℓ̂ table and its degenerate rows
(docs/specs/M9-16-report-content.md; GIPPA_SPEC.md §7.4, §7.7.1, §7.8.2; ADR-0017, ADR-0021
Amendment 1, ADR-0025, ADR-0026, ADR-0027, ADR-0032, ADR-0033, ADR-0034; backlog M6-06).

Interface pinned here (M9-16 test author, beyond the spec text):
- conditions: operon note when ``pairs.tsv`` or ``operon_flags.tsv`` exists; §7.8.2 and
  ADR-0032/33/34 caveats when ``tables/determinants_model.tsv`` exists; GSEA caveats when
  ``tables/enrichment_gsea.tsv`` exists (ADR-0025 when a row is flagged ``genomic_cluster``);
  function coverage, Bakta-light and KEGG caveats when ``tables/category_summary.tsv`` exists;
  islands validity sentence when ``tables/islands.tsv`` exists;
- the KEGG caveat reads "KEGG retrieval date: <date>" (decision 2), the date the ``retrieved``
  value of the manifest's ``kegg*`` database entries; "KEGG retrieval date: not recorded in this
  run" without one; "KEGG not installed" (and no academic-use sentence) under
  ``layer_missing:kegg``;
- ``<table data-table="island_diagnostics">`` in ``genome-context``: one row per replicon, a
  header cell containing each of ℓ̂, ℓ′ and Î₁; ℓ̂ printed with
  ``gippa.context.diagnostics.render_l_hat`` (never ``0``);
- ``<table data-table="missed_operons_bins">`` in ``operons``: the ``bin_counts`` rows, bounds
  ``−∞`` (U+2212) / ``∞``, and π₀ printed beside it.
ADR sentences are matched case-insensitively on whitespace-normalised section text, with room
for the article and verb the prose needs.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from helpers.report_content import (
    DETERMINANTS_CAVEATS,
    ISLANDS_VALIDITY,
    KEGG_DATE,
    det_rows,
    manifest,
    page,
    section_text,
    write_category_summary,
    write_det_model,
    write_diagnostics,
    write_flags,
    write_gsea,
    write_islands,
    write_manifest,
    write_missed_event,
    write_operon_flags,
    write_pairs,
)
from helpers.report_fixtures import run_out, table_rows


def _has(pattern: str, text: str) -> bool:
    return re.search(pattern, text, re.IGNORECASE) is not None


# ---------------------------------------------------------------- operons (§7.4)


@pytest.mark.parametrize("which", ["pairs", "operon_flags"])
def test_operons_sigma70_intrinsic_note(tmp_path: Path, which: str) -> None:
    out = run_out(tmp_path)
    write_pairs(out) if which == "pairs" else write_operon_flags(out, [0.02])
    text = section_text(page(out), "operons")
    assert "σ70 promoters only" in text, text
    assert "intrinsic terminators only" in text, text
    assert _has(r"missing call (therefore )?never proves a missing element", text), text


# ---------------------------------------------------------------- determinants


def _det_out(tmp_path: Path) -> Path:
    out = run_out(tmp_path)
    write_det_model(out, det_rows("elastic_net", 0.21, 0.03, 0.001, "null_uncalibrated"))
    return out


def test_determinants_section_7_8_2_caveats_verbatim(tmp_path: Path) -> None:
    text = section_text(page(_det_out(tmp_path)), "determinants")
    for caveat in DETERMINANTS_CAVEATS:
        assert " ".join(caveat.split()) in text, caveat


def test_determinants_caveats_match_the_module_constant() -> None:
    """``gippa.determinants.CAVEATS`` (#178) is the source the report renders, verbatim §7.8.2."""
    import gippa.determinants as det

    if not hasattr(det, "CAVEATS"):
        pytest.skip("gippa.determinants.CAVEATS arrives with #178 (M8-10)")
    assert tuple(det.CAVEATS) == DETERMINANTS_CAVEATS  # type: ignore[attr-defined]


@pytest.mark.parametrize(
    ("adr", "pattern"),
    [
        ("ADR-0032", r"univariate\b.{0,40}\bCIs?\b.{0,40}\buncalibrated"),
        ("ADR-0033", r"elastic-net ΔR² is the primary (block )?attribution"),
        ("ADR-0033", r"LightGBM.{0,40}SHAP.{0,40}secondary"),
        ("ADR-0034", r"null p.{0,30}uncalibrated"),
    ],
)
def test_determinants_adr_caveats(tmp_path: Path, adr: str, pattern: str) -> None:
    text = section_text(page(_det_out(tmp_path)), "determinants")
    assert _has(pattern, text), f"{adr}: /{pattern}/ not in {text!r}"


def test_no_determinants_caveats_without_the_table(tmp_path: Path) -> None:
    text = section_text(page(run_out(tmp_path)), "determinants")
    assert "synthetic GFP library" not in text, text


# ---------------------------------------------------------------- GSEA (function)


@pytest.mark.parametrize(
    ("adr", "pattern"),
    [
        ("ADR-0026", r"conditional on detection,? (with )?no p_detect adjustment"),
        ("ADR-0027", r"NES\b.{0,30}\b(no|without (a|an))\b.{0,20}interval"),
    ],
)
def test_gsea_caveats(tmp_path: Path, adr: str, pattern: str) -> None:
    out = run_out(tmp_path)
    write_gsea(out, [None])
    text = section_text(page(out), "function")
    assert _has(pattern, text), f"{adr}: /{pattern}/ not in {text!r}"


def test_gsea_genomic_cluster_caveat(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_gsea(out, [None, "genomic_cluster"])
    text = section_text(page(out), "function")
    assert _has(r"no level claim.{0,60}genomic_cluster", text), text


# ---------------------------------------------------------------- function coverage, Bakta, KEGG


def _function_out(tmp_path: Path) -> Path:
    out = run_out(tmp_path)
    write_category_summary(out, {"cog": 0.62, "kegg": 0.41, "go_slim": 0.55})
    return out


def test_function_coverage_caveat_with_numbers(tmp_path: Path) -> None:
    text = section_text(page(_function_out(tmp_path)), "function")
    assert _has(r"annotation coverage leans toward well-studied genes", text), text
    for value in ("0.62", "0.41", "0.55"):
        assert value in text, (value, text)


def test_bakta_light_caveat_when_flagged(tmp_path: Path) -> None:
    out = _function_out(tmp_path)
    write_flags(out, "annotate", ["bakta_db_light"])
    assert "Bakta light DB" in section_text(page(out), "function")


def test_no_bakta_light_caveat_without_the_flag(tmp_path: Path) -> None:
    assert "Bakta light DB" not in section_text(page(_function_out(tmp_path)), "function")


def test_kegg_academic_use_and_retrieval_date_from_manifest(tmp_path: Path) -> None:
    out = _function_out(tmp_path)
    write_manifest(out)
    text = section_text(page(out), "function")
    assert "academic use only" in text and "rest.kegg.jp" in text, text
    assert f"KEGG retrieval date: {KEGG_DATE}" in text, text


@pytest.mark.parametrize("case", ["no_manifest", "no_kegg_entry", "empty_databases"])
def test_kegg_retrieval_date_not_recorded_never_guessed(tmp_path: Path, case: str) -> None:
    """Decision 2: no recorded date -> "not recorded in this run"; nothing date-like printed."""
    out = _function_out(tmp_path)
    if case == "no_kegg_entry":
        doc = manifest()
        doc["databases"] = [d for d in doc["databases"] if not d["name"].startswith("kegg")]
        write_manifest(out, doc)
    elif case == "empty_databases":
        write_manifest(out, manifest(databases=[]))
    text = section_text(page(out), "function")
    assert "KEGG retrieval date: not recorded in this run" in text, text
    m = re.search(r"KEGG retrieval date: not recorded in this run(.{0,30})", text)
    assert m is not None and not re.search(r"\d{4}-\d{2}-\d{2}", m.group(1)), text


def test_kegg_not_installed_replaces_the_caveat(tmp_path: Path) -> None:
    out = _function_out(tmp_path)
    write_flags(out, "function", ["layer_missing:kegg"])
    write_manifest(out)
    text = section_text(page(out), "function")
    assert "KEGG not installed" in text, text
    assert "academic use only" not in text, text


# ---------------------------------------------------------------- islands (ADR-0021)


def test_islands_validity_domain_sentence(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_islands(out, [0.001])
    assert ISLANDS_VALIDITY in section_text(page(out), "genome-context")


def _diag_rows() -> list[dict[str, object]]:
    return [
        {"contig": "chr1", "scanned": True, "block_rule": "auto", "l_hat": 4, "block_len": 9,
         "morans_i_lag1": 0.4213},
        {"contig": "pA", "scanned": True, "block_rule": "fixed", "l_hat": None, "block_len": 12,
         "morans_i_lag1": 0.1057},
        {"contig": "pB", "scanned": False, "block_rule": "not_scanned", "l_hat": None,
         "block_len": None, "morans_i_lag1": float("nan"), "flags": ["few_values"]},
    ]  # fmt: skip


def _col(header: list[str], symbol: str) -> int:
    hits = [i for i, h in enumerate(header) if symbol in h]
    assert len(hits) == 1, (symbol, header)
    return hits[0]


def test_island_diagnostics_table_per_replicon(tmp_path: Path) -> None:
    from gippa.context.diagnostics import render_l_hat

    out = run_out(tmp_path)
    write_islands(out, [0.001])
    write_diagnostics(out, _diag_rows())
    p = page(out)
    table = p.table("island_diagnostics")
    assert table is not None, "no <table data-table='island_diagnostics'>"
    assert table in p.section("genome-context").iter()
    header, body = table_rows(table)
    lh, lp, i1 = _col(header, "ℓ̂"), _col(header, "ℓ′"), _col(header, "Î₁")
    rows = {r[0]: r for r in body}
    assert sorted(rows) == ["chr1", "pA", "pB"], body
    assert rows["chr1"][lh] == "4" and rows["chr1"][lp] == "9", rows["chr1"]
    assert rows["chr1"][i1].startswith("0.421"), rows["chr1"]
    assert rows["pA"][lh] == render_l_hat(None, "fixed"), rows["pA"]
    assert rows["pA"][lp] == "12", rows["pA"]
    assert rows["pB"][lh] == render_l_hat(None, "not_scanned"), rows["pB"]
    for r in body:
        assert r[lh] != "0" and r[lp] != "0", r


def test_island_diagnostics_absent_gives_a_note(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_islands(out, [0.001])
    p = page(out)
    assert p.table("island_diagnostics") is None
    text = section_text(p, "genome-context")
    assert "island_diagnostics" in text, text


# ---------------------------------------------------------------- operons bin counts (M6-06)


def test_missed_operons_bin_counts_and_pi0(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_pairs(out)
    write_missed_event(out, pi0=0.5)
    p = page(out)
    table = p.table("missed_operons_bins")
    assert table is not None, "no <table data-table='missed_operons_bins'>"
    assert table in p.section("operons").iter()
    _, body = table_rows(table)
    assert len(body) == 3, body
    cells = [c for r in body for c in r]
    assert "−∞" in cells and "∞" in cells, body
    text = table.text()  # the page's own sort script (M9-14) mentions Infinity
    assert "Infinity" not in text, text
    for n in ("3", "7", "11", "1", "4", "6", "2", "5", "9"):
        assert n in cells, (n, body)
    text = section_text(p, "operons")
    assert re.search(r"(π₀|pi0)\D{0,5}0\.5\b", text), text


def test_missed_operons_table_omitted_without_event(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    write_pairs(out)
    assert page(out).table("missed_operons_bins") is None
