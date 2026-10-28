"""M8-06a item 5 and example 9: Frictionless schemas of ``category_summary``,
``enrichment_detection`` and ``enrichment_gsea`` (docs/specs/M8-06-function-step-and-run.md;
GIPPA_SPEC.md §9.3, §9.5; ADR-0029; M8-04 decisions 3 and 13).

Pinned here:

* field names == each module's ``COLUMNS`` (so a #179 column change fails loudly);
* ``collection`` enums (6 for the enrichment tables, {cog, kegg, go_slim} for the summary),
  ``direction`` {over, under}, ``method`` / ``null_model`` enums, ``n_perm`` = 0 in detection;
* the ``flags`` token patterns of item 5, plus ``tu_recal_skipped`` on ORA rows (lead decision 3 /
  M8-04 decision 13 put it on every ORA row of a skipped recalibration; this overrides item 5's
  "not in any table" for ``tu_recal_skipped`` only); the other run-scope tokens are rejected
  everywhere;
* counts are integers >= 0, ``p_*`` / ``q_value`` / ``coverage`` in [0, 1];
* ``missingValues ["NA", ""]``, the ``coordinates`` key as ``pairs.schema.json``, and "in TUs"
  in the descriptions of the TU counts.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

SCHEMA_DIR = Path(__file__).resolve().parents[3] / "src" / "gippa" / "schemas" / "tables"
ALL6 = ["cog", "kegg", "go_bp", "go_mf", "go_cc", "go_slim"]
RUN_ONLY = ["layer_missing:cog", "layer_missing:kegg", "layer_missing:go",
            "layer_missing:abundance", "bakta_db_light", "too_few_ranked_tus",
            "no_eligible_genes"]  # fmt: skip
FLAG_TOKENS = {
    "category_summary": ["ci_few_units"],
    "enrichment_detection": ["ci_unestimated", "unestimated", "odds_out_of_range", "elim_empty",
                             "null_degenerate", "p_detect_clipped", "ora_uncalibrated",
                             "tu_recal_skipped"],
    "enrichment_gsea": ["ci_unestimated", "unestimated", "size_filtered", "genomic_cluster",
                        "detection_confounded"],
}  # fmt: skip

SUMMARY_ROW = {"collection": "cog", "category": "J", "name": "Translation", "n_genes": "3",
               "n_detected": "2", "expected": "1.8", "o_minus_e": "0.2", "ci_low": "NA",
               "ci_high": "NA", "median_abundance_log2": "NA", "coverage": "0.75",
               "flags": "ci_few_units"}  # fmt: skip
DETECTION_ROW = {"collection": "go_bp", "term": "GO:0000101", "name": "toy", "direction": "over",
                 "n_universe": "40", "n_in_set": "12", "n_detected": "20",
                 "k_detected_in_set": "9", "n_in_set_genes": "15",
                 "k_detected_in_set_genes": "11", "odds_weight": "1.2", "p_classic": "0.01",
                 "p_noelim": "0.02", "parent_group": "NA", "coverage": "0.5", "stat": "9",
                 "effect": "2.5", "ci_low": "NA", "ci_high": "NA", "p_value": "0.01",
                 "q_value": "0.04", "method": "cpb_tu_ora",
                 "null_model": "conditional_poisson_binomial_tu", "n_perm": "0",
                 "flags": "ci_unestimated;ora_uncalibrated"}  # fmt: skip
GSEA_ROW = {"collection": "go_bp", "term": "GO:0000101", "name": "toy", "size_tus": "12",
            "es": "0.5", "nes": "1.3", "log2err": "NA", "leading_edge": "g1;g2",
            "redundant": "false", "parent_group": "GO:0000101", "coverage": "0.5", "stat": "0.5",
            "effect": "1.3", "ci_low": "NA", "ci_high": "NA", "p_value": "0.01",
            "q_value": "0.04", "method": "fgsea_simple", "null_model": "tu_set_permutation",
            "n_perm": "1000", "flags": "ci_unestimated"}  # fmt: skip
ROWS = {"category_summary": SUMMARY_ROW, "enrichment_detection": DETECTION_ROW,
        "enrichment_gsea": GSEA_ROW}  # fmt: skip
TABLES = list(ROWS)


def _schema(stem: str) -> dict[str, Any]:
    return json.loads((SCHEMA_DIR / f"{stem}.schema.json").read_text(encoding="utf-8"))  # type: ignore[no-any-return]


def _field(stem: str, name: str) -> dict[str, Any]:
    (f,) = [f for f in _schema(stem)["fields"] if f["name"] == name]
    return f  # type: ignore[no-any-return]


def _columns(stem: str) -> tuple[str, ...]:
    from gippa.function import gsea, ora, summary

    return {"category_summary": tuple(summary.COLUMNS),
            "enrichment_detection": tuple(ora.COLUMNS),
            "enrichment_gsea": tuple(gsea.COLUMNS)}[stem]  # fmt: skip


def _valid(stem: str, rows: list[dict[str, str]], tmp_path: Path) -> bool:
    frictionless = pytest.importorskip("frictionless")
    cols = _columns(stem)
    path = tmp_path / f"{stem}.tsv"
    lines = ["\t".join(cols)] + ["\t".join(r[c] for c in cols) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    d = {k: v for k, v in _schema(stem).items() if k != "coordinates"}
    res = frictionless.Resource(path=path.name, basepath=str(path.parent),
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    return bool(res.validate().valid)


def _with(stem: str, **cells: str) -> list[dict[str, str]]:
    return [ROWS[stem] | cells]


# ----------------------------------------------------------------------------- structure


@pytest.mark.parametrize("stem", TABLES)
def test_schema_fields_equal_module_columns(stem: str) -> None:
    assert tuple(f["name"] for f in _schema(stem)["fields"]) == _columns(stem)


@pytest.mark.parametrize("stem", TABLES)
def test_schema_metadata(stem: str) -> None:
    s = _schema(stem)
    assert s["missingValues"] == ["NA", ""]
    assert s["coordinates"] == _schema("pairs")["coordinates"]
    assert s.get("name") == stem
    assert all(f.get("description") for f in s["fields"]), "every field is described"


@pytest.mark.parametrize("stem", TABLES)
def test_base_row_is_valid(stem: str, tmp_path: Path) -> None:
    assert _valid(stem, [ROWS[stem]], tmp_path)


def test_tu_count_descriptions_state_the_unit() -> None:
    for name in ("n_universe", "n_in_set", "n_detected", "k_detected_in_set"):
        assert "in TUs" in _field("enrichment_detection", name)["description"], name
    assert "in TUs" in _field("enrichment_gsea", "size_tus")["description"]


# ----------------------------------------------------------------------------- enums


def test_collection_enums() -> None:
    for stem in ("enrichment_detection", "enrichment_gsea"):
        assert _field(stem, "collection")["constraints"]["enum"] == ALL6, stem
    want = ["cog", "kegg", "go_slim"]
    assert _field("category_summary", "collection")["constraints"]["enum"] == want


def test_method_and_null_model_enums() -> None:
    assert _field("enrichment_detection", "method")["constraints"]["enum"] == ["cpb_tu_ora"]
    assert _field("enrichment_gsea", "method")["constraints"]["enum"] == [
        "fgsea_simple",
        "fgsea_multilevel",
    ]
    assert _field("enrichment_detection", "null_model")["constraints"]["enum"] == [
        "conditional_poisson_binomial_tu"
    ]
    assert _field("enrichment_gsea", "null_model")["constraints"]["enum"] == ["tu_set_permutation"]
    assert _field("enrichment_detection", "direction")["constraints"]["enum"] == ["over", "under"]


@pytest.mark.parametrize(
    ("stem", "col", "value"),
    [
        ("category_summary", "collection", "go_bp"),
        ("category_summary", "collection", "pfam"),
        ("enrichment_detection", "collection", "pfam"),
        ("enrichment_gsea", "collection", "pfam"),
        ("enrichment_detection", "direction", "both"),
        ("enrichment_detection", "method", "wallenius_ora"),  # example 9
        ("enrichment_gsea", "method", "wallenius_ora"),  # example 9
        ("enrichment_detection", "null_model", "tu_set_permutation"),
        ("enrichment_gsea", "null_model", "conditional_poisson_binomial_tu"),
        ("enrichment_detection", "n_perm", "1000"),
    ],
)
def test_enum_violations_are_rejected(stem: str, col: str, value: str, tmp_path: Path) -> None:
    assert not _valid(stem, _with(stem, **{col: value}), tmp_path)


def test_fgsea_multilevel_is_already_accepted(tmp_path: Path) -> None:
    """Reserved for M8-11 so the switch needs no schema change."""
    assert _valid("enrichment_gsea", _with("enrichment_gsea", method="fgsea_multilevel"), tmp_path)


@pytest.mark.parametrize("coll", ALL6)
def test_every_collection_accepted_in_enrichment_tables(coll: str, tmp_path: Path) -> None:
    for stem in ("enrichment_detection", "enrichment_gsea"):
        assert _valid(stem, _with(stem, collection=coll), tmp_path), (stem, coll)


# ----------------------------------------------------------------------------- flags


@pytest.mark.parametrize("stem", TABLES)
def test_every_listed_flag_token_is_accepted(stem: str, tmp_path: Path) -> None:
    toks = FLAG_TOKENS[stem]
    rows = [ROWS[stem] | {"flags": t} for t in toks] + [ROWS[stem] | {"flags": ";".join(toks)}]
    rows.append(ROWS[stem] | {"flags": "NA"})
    assert _valid(stem, rows, tmp_path)


@pytest.mark.parametrize("stem", TABLES)
@pytest.mark.parametrize("bad", ["foo", *RUN_ONLY])
def test_unknown_and_run_scope_flags_are_rejected(stem: str, bad: str, tmp_path: Path) -> None:
    """Example 9 (flag ``foo``) and item 5 (run-only tokens go to the event, not to rows)."""
    assert not _valid(stem, _with(stem, flags=bad), tmp_path)
    joined = f"{FLAG_TOKENS[stem][0]};{bad}"
    assert not _valid(stem, _with(stem, flags=joined), tmp_path)


@pytest.mark.parametrize("stem", ["category_summary", "enrichment_gsea"])
def test_tu_recal_skipped_only_on_ora_rows(stem: str, tmp_path: Path) -> None:
    assert not _valid(stem, _with(stem, flags="tu_recal_skipped"), tmp_path)


@pytest.mark.parametrize("stem", TABLES)
def test_malformed_flag_lists_are_rejected(stem: str, tmp_path: Path) -> None:
    tok = FLAG_TOKENS[stem][0]
    for bad in (f"{tok};", f";{tok}", f"{tok};;{tok}", f"{tok} ", f"{tok},{tok}"):
        assert not _valid(stem, _with(stem, flags=bad), tmp_path), bad


# ----------------------------------------------------------------------------- numbers


@pytest.mark.parametrize(
    ("stem", "col"),
    [
        ("category_summary", "n_genes"), ("category_summary", "n_detected"),
        ("enrichment_detection", "n_universe"), ("enrichment_detection", "n_in_set"),
        ("enrichment_detection", "n_detected"), ("enrichment_detection", "k_detected_in_set"),
        ("enrichment_detection", "n_in_set_genes"),
        ("enrichment_detection", "k_detected_in_set_genes"),
        ("enrichment_gsea", "size_tus"), ("enrichment_gsea", "n_perm"),
    ],
)  # fmt: skip
def test_counts_are_nonnegative_integers(stem: str, col: str, tmp_path: Path) -> None:
    assert _field(stem, col)["type"] == "integer"
    assert not _valid(stem, _with(stem, **{col: "-1"}), tmp_path)
    assert not _valid(stem, _with(stem, **{col: "1.5"}), tmp_path)


@pytest.mark.parametrize(
    ("stem", "col"),
    [
        ("category_summary", "coverage"),
        ("enrichment_detection", "p_classic"), ("enrichment_detection", "p_noelim"),
        ("enrichment_detection", "p_value"), ("enrichment_detection", "q_value"),
        ("enrichment_detection", "coverage"),
        ("enrichment_gsea", "p_value"), ("enrichment_gsea", "q_value"),
        ("enrichment_gsea", "coverage"),
    ],
)  # fmt: skip
def test_probabilities_and_coverage_in_unit_interval(stem: str, col: str, tmp_path: Path) -> None:
    assert _field(stem, col)["type"] == "number"
    for ok in ("0", "1", "0.5", "1e-300"):
        assert _valid(stem, _with(stem, **{col: ok}), tmp_path), (col, ok)
    for bad in ("-0.01", "1.5"):
        assert not _valid(stem, _with(stem, **{col: bad}), tmp_path), (col, bad)


def test_gsea_redundant_is_boolean(tmp_path: Path) -> None:
    assert _field("enrichment_gsea", "redundant")["type"] == "boolean"
    for ok in ("true", "false", "NA"):
        assert _valid("enrichment_gsea", _with("enrichment_gsea", redundant=ok), tmp_path)
    assert not _valid("enrichment_gsea", _with("enrichment_gsea", redundant="maybe"), tmp_path)
