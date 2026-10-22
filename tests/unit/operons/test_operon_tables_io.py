"""M6-08 writers and schemas of ``tables/pairs.tsv``, ``tables/operons.tsv`` and
``tables/operon_flags.tsv`` (docs/specs/M6-08-operons-step-and-run.md, public interface items 3
and 4, example 6; GIAP_SPEC.md §6.5, §6.6, §9.2, §9.3, §9.4, §9.5).

Interface pinned here (test-author reading of the spec):

* ``giap.io.operon_tables.write_pairs_tsv(path, df)``, ``write_operons_tsv(path, df)`` and
  ``write_operon_flags_tsv(path, df)``. They take the internal frames of M6-04..M6-07
  (``merge_tests``, ``operons_table`` and ``flag_pairs``), which are 0-based half-open.
* Cells follow ``giap.io.tsv.format_value`` (M1-04): ``NA`` for missing, floats ``%.6g``,
  booleans ``true``/``false`` and ``;``-joined lists, with an empty list written as ``NA``.
* Coordinates are converted **only** in the writer (§6.5). In ``operons.tsv`` that means
  ``start + 1`` with ``end`` unchanged, and an origin-crossing TU keeps ``end > length``. The two
  pair tables carry no coordinates, so ``intergenic_bp`` is written as it is.
* Row order (§6.5): ``pairs.tsv`` and ``operon_flags.tsv`` are sorted by the number of
  ``PAIR_<n>``, because the pair-ID order is the contig/start order (M6-index decision 3).
  ``operons.tsv`` keeps the order ``operons_table`` gives (contig in C1 order, then start).
* ``pairs.tsv`` never carries the internal ``which`` column (lead decision 4).
* Schemas live in ``src/giap/schemas/tables/{pairs,operons,operon_flags}.schema.json``, the
  per-table resources the other tables already use. Each declares ``coordinates:
  1-based-closed``. The ``class``, ``test`` and ``flag`` enums are exactly those of §9.3.
* ``pairs.tsv`` ``flags`` accepts ``rna_between`` (§9.4 Operons) and exactly the M6-08 lead
  decision 7 tokens: ``abundance_imputed``, ``layer_missing:promoters``,
  ``layer_missing:terminators``, ``ci_unestimated``, ``unestimated``, ``too_few_op_pairs``,
  ``too_few_no_pairs``, ``no_detection_limit_model`` and ``relative_scale``. The flags field
  documents both ``abundance_imputed`` (pair scope) and ``imputed_abundance`` (gene scope). It
  rejects ``uncalibrated`` (index risk "uncalibrated must not come back"; ADR-0021 Amendment 1).
"""

from __future__ import annotations

import json
import math
import string
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.io.tsv import format_value
from giap.operons.flags import FLAG_COLUMNS, FLAGS
from giap.operons.pairs import PAIR_COLUMNS, TEST_COLUMNS
from giap.operons.tu_table import OPERON_COLUMNS

SCHEMA_DIR = Path(__file__).resolve().parents[3] / "src" / "giap" / "schemas" / "tables"

#: §9.3 enums, copied from GIAP_SPEC.md (not from the code under test)
SPEC_CLASS = {"OP", "NO", "CV", "DV", "SW"}
SPEC_TEST = {"both_detected_t", "one_detected_pmiss", "missed_operon_decoy", "none"}
SPEC_FLAG = {
    "internal_promoter_supported", "internal_terminator_supported", "unexplained_discordance",
    "conflict_terminator_in_operon", "novel_operon_candidate", "coabundant_despite_terminator",
    "other",
}  # fmt: skip
#: lead decision 7 (M6-08): the pair-scope tokens added through the ``pairs.tsv`` schema
DECISION_7 = (
    "abundance_imputed", "layer_missing:promoters", "layer_missing:terminators",
    "ci_unestimated", "unestimated", "too_few_op_pairs", "too_few_no_pairs",
    "no_detection_limit_model", "relative_scale",
)  # fmt: skip
#: every token ``pairs.tsv`` ``flags`` must accept: §9.4 Operons ``rna_between`` + decision 7
OPERON_TOKENS = ("rna_between", *DECISION_7)


def _writers() -> Any:
    from giap.io import operon_tables

    return operon_tables


def _lines(path: Path) -> list[list[str]]:
    text = path.read_text(encoding="utf-8")
    assert text.endswith("\n"), "TSV must end with a newline"
    return [line.split("\t") for line in text[:-1].split("\n")]


def _schema(stem: str) -> dict[str, Any]:
    path = SCHEMA_DIR / f"{stem}.schema.json"
    assert path.is_file(), f"missing Frictionless resource {path}"
    data: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return data


def _field(stem: str, name: str) -> dict[str, Any]:
    found = [f for f in _schema(stem)["fields"] if f["name"] == name]
    assert len(found) == 1, f"{stem}: field {name} declared {len(found)} times"
    return found[0]


def _validate(path: Path, stem: str) -> Any:
    """Frictionless report of ``path`` against its resource (``coordinates`` is a GIAP key)."""
    frictionless = pytest.importorskip("frictionless")
    d = {k: v for k, v in _schema(stem).items() if k != "coordinates"}
    schema = frictionless.Schema.from_descriptor(d)
    res = frictionless.Resource(path=path.name, basepath=str(path.parent), schema=schema)
    return res.validate()


# ----------------------------------------------------------------------------- frames


def pairs_frame() -> pd.DataFrame:
    """Three merged-pair rows, out of ``PAIR_<n>`` order, with the internal ``which`` column."""
    nan = math.nan
    rows: list[dict[str, Any]] = [
        {"pair_id": "PAIR_10", "class": "SW", "gene_a": "G9", "gene_b": "G11", "strand": "-",
         "intergenic_bp": 412, "uniop_p": nan, "promoter_between": 0.25,
         "tx_rate_between": nan, "terminator_between": 0.0, "detected_a": True,
         "detected_b": False, "p_detect_a": 0.75, "p_detect_b": 0.125,
         "abundance_log2_a": 11.5, "abundance_log2_b": nan, "delta_log2": nan, "test": "none",
         "flags": []},
        {"pair_id": "PAIR_2", "class": "OP", "gene_a": "G2", "gene_b": "G3", "strand": "+",
         "intergenic_bp": -4, "uniop_p": 0.9731234, "promoter_between": 0.1,
         "tx_rate_between": 123.456789, "terminator_between": 0.82, "detected_a": True,
         "detected_b": True, "p_detect_a": 0.5, "p_detect_b": 0.6, "abundance_log2_a": 10.0,
         "abundance_log2_b": 12.25, "delta_log2": 2.25, "test": "both_detected_t",
         "flags": ["rna_between", "ci_unestimated"]},
        {"pair_id": "PAIR_1", "class": "CV", "gene_a": "G1", "gene_b": "G2", "strand": ".",
         "intergenic_bp": 0, "uniop_p": nan, "promoter_between": nan, "tx_rate_between": nan,
         "terminator_between": nan, "detected_a": False, "detected_b": True,
         "p_detect_a": nan, "p_detect_b": 0.5, "abundance_log2_a": nan,
         "abundance_log2_b": 10.0, "delta_log2": nan, "test": "none",
         "flags": ["layer_missing:terminators"]},
    ]  # fmt: skip
    df = pd.DataFrame(rows)
    for col in TEST_COLUMNS:
        df[col] = None if col in ("method", "null_model") else nan
    df.loc[1, ["stat", "effect", "p_value", "q_value"]] = [2.25, 2.0, 0.0123456789, 0.05]
    df.loc[1, ["method", "null_model"]] = ["student_t", "within_operon_t"]
    df["n_perm"] = pd.array([pd.NA, 5, pd.NA], dtype="Int64")
    df["which"] = ["", "sandwich", ""]
    return df[[*PAIR_COLUMNS, "which"]]


def operons_frame() -> pd.DataFrame:
    """Two TUs (0-based half-open), the second crossing the origin of a 100-kb chr1."""
    return pd.DataFrame(
        {
            "tu_id": ["TU_1", "TU_2"],
            "contig": ["chr1", "chr1"],
            "start": np.array([0, 99_900], dtype=np.int64),
            "end": np.array([300, 100_200], dtype=np.int64),
            "strand": ["+", "-"],
            "genes": [["G1", "G2"], ["G7"]],
            "n_genes": np.array([2, 1], dtype=np.int64),
            "detected_frac": [0.5, math.nan],
            "median_abundance_log2": [11.0, math.nan],
            "max_abs_delta_log2": [0.333333333, math.nan],
            "n_flags": np.array([1, 0], dtype=np.int64),
        },
        columns=list(OPERON_COLUMNS),
    )


def flags_frame() -> pd.DataFrame:
    """Two ``flag_pairs`` rows, out of ``PAIR_<n>`` order."""
    return pd.DataFrame(
        {
            "pair_id": ["PAIR_12", "PAIR_3"],
            "gene_a": ["G12", "G3"],
            "gene_b": ["G13", "G4"],
            "uniop_class": ["NO", "OP"],
            "promoter_between_yes": [False, True],
            "terminator_between_yes": [True, False],
            "abundance_pattern": ["coabundant", "downstream_higher"],
            "flag": ["coabundant_despite_terminator", "internal_promoter_supported"],
            "q_value": [math.nan, 0.04],
            "lfdr": [math.nan, 0.2],
            "evidence_note": ["terminator 0.91 between", "promoter 0.80 between"],
        },
        columns=list(FLAG_COLUMNS),
    )


def _empty(cols: tuple[str, ...]) -> pd.DataFrame:
    return pd.DataFrame({c: pd.Series(dtype=object) for c in cols})


# ----------------------------------------------------------------------------- pairs.tsv


def test_pairs_header_is_the_9_3_dictionary_without_which(tmp_path: Path) -> None:
    out = tmp_path / "pairs.tsv"
    df = pairs_frame()
    try:
        _writers().write_pairs_tsv(out, df)
    except ValueError:  # refusing the internal column is also acceptable
        _writers().write_pairs_tsv(out, df.drop(columns="which"))
    header = _lines(out)[0]
    assert header == list(PAIR_COLUMNS)
    assert "which" not in header


def test_pairs_rows_sorted_by_pair_number(tmp_path: Path) -> None:
    out = tmp_path / "pairs.tsv"
    _writers().write_pairs_tsv(out, pairs_frame().drop(columns="which"))
    ids = [r[0] for r in _lines(out)[1:]]
    assert ids == ["PAIR_1", "PAIR_2", "PAIR_10"]


def test_pairs_cells_follow_the_tsv_conventions(tmp_path: Path) -> None:
    out = tmp_path / "pairs.tsv"
    df = pairs_frame().drop(columns="which")
    _writers().write_pairs_tsv(out, df)
    got = {r[0]: dict(zip(PAIR_COLUMNS, r, strict=True)) for r in _lines(out)[1:]}
    op = got["PAIR_2"]
    assert op["intergenic_bp"] == "-4"  # not a coordinate: never shifted
    assert op["detected_a"] == "true" and got["PAIR_1"]["detected_a"] == "false"
    assert op["uniop_p"] == format_value(0.9731234)
    assert op["tx_rate_between"] == format_value(123.456789)
    assert op["p_value"] == format_value(0.0123456789)
    assert op["n_perm"] == "5" and got["PAIR_1"]["n_perm"] == "NA"
    assert op["method"] == "student_t" and got["PAIR_1"]["method"] == "NA"
    assert op["flags"] == "rna_between;ci_unestimated"
    assert got["PAIR_1"]["flags"] == "layer_missing:terminators"
    assert got["PAIR_10"]["flags"] == "NA"  # empty list
    sw = got["PAIR_10"]
    for col in ("uniop_p", "tx_rate_between", "abundance_log2_b", "delta_log2", "ci_low", "lfdr"):
        assert sw[col] == "NA", col
    assert sw["terminator_between"] == "0"


def test_pairs_writer_does_not_modify_its_input(tmp_path: Path) -> None:
    df = pairs_frame().drop(columns="which")
    before = df.copy(deep=True)
    _writers().write_pairs_tsv(tmp_path / "pairs.tsv", df)
    pd.testing.assert_frame_equal(df, before)


def test_pairs_empty_frame_is_header_only(tmp_path: Path) -> None:
    out = tmp_path / "pairs.tsv"
    _writers().write_pairs_tsv(out, _empty(PAIR_COLUMNS))
    assert out.read_text(encoding="utf-8") == "\t".join(PAIR_COLUMNS) + "\n"


# ----------------------------------------------------------------------------- operons.tsv


def test_operons_header_and_one_based_closed_coordinates(tmp_path: Path) -> None:
    out = tmp_path / "operons.tsv"
    df = operons_frame()
    _writers().write_operons_tsv(out, df)
    lines = _lines(out)
    assert lines[0] == list(OPERON_COLUMNS)
    rows = [dict(zip(OPERON_COLUMNS, r, strict=True)) for r in lines[1:]]
    assert [r["tu_id"] for r in rows] == ["TU_1", "TU_2"]
    assert (rows[0]["start"], rows[0]["end"]) == ("1", "300")
    # origin-crossing TU: 1-based start, end unchanged and > replicon length (§6.5)
    assert (rows[1]["start"], rows[1]["end"]) == ("99901", "100200")
    assert rows[0]["genes"] == "G1;G2" and rows[1]["genes"] == "G7"
    assert rows[0]["n_genes"] == "2" and rows[1]["n_flags"] == "0"
    assert rows[0]["detected_frac"] == "0.5" and rows[1]["detected_frac"] == "NA"
    assert rows[0]["max_abs_delta_log2"] == format_value(0.333333333)


def test_operons_conversion_happens_only_in_the_writer(tmp_path: Path) -> None:
    df = operons_frame()
    _writers().write_operons_tsv(tmp_path / "operons.tsv", df)
    assert df["start"].tolist() == [0, 99_900]  # the internal frame stays 0-based half-open
    assert df["end"].tolist() == [300, 100_200]


def test_operons_empty_frame_is_header_only(tmp_path: Path) -> None:
    out = tmp_path / "operons.tsv"
    _writers().write_operons_tsv(out, _empty(OPERON_COLUMNS))
    assert out.read_text(encoding="utf-8") == "\t".join(OPERON_COLUMNS) + "\n"


@settings(
    max_examples=60, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(
    spans=st.lists(
        st.tuples(st.integers(0, 10_000_000), st.integers(1, 50_000)), min_size=1, max_size=8
    )
)
def test_property_operons_start_plus_one_end_unchanged(
    tmp_path: Path, spans: list[tuple[int, int]]
) -> None:
    n = len(spans)
    df = pd.DataFrame(
        {
            "tu_id": [f"TU_{i + 1}" for i in range(n)],
            "contig": ["chr1"] * n,
            "start": np.array([s for s, _ in spans], dtype=np.int64),
            "end": np.array([s + w for s, w in spans], dtype=np.int64),
            "strand": ["+"] * n,
            "genes": [[f"G{i}"] for i in range(n)],
            "n_genes": np.ones(n, dtype=np.int64),
            "detected_frac": [math.nan] * n,
            "median_abundance_log2": [math.nan] * n,
            "max_abs_delta_log2": [math.nan] * n,
            "n_flags": np.zeros(n, dtype=np.int64),
        },
        columns=list(OPERON_COLUMNS),
    )
    out = tmp_path / "operons_prop.tsv"
    _writers().write_operons_tsv(out, df)
    rows = _lines(out)[1:]
    i_s, i_e = OPERON_COLUMNS.index("start"), OPERON_COLUMNS.index("end")
    got = {r[0]: (int(r[i_s]), int(r[i_e])) for r in rows}  # the writer may sort by start (§6.5)
    assert got == {f"TU_{i + 1}": (s + 1, s + w) for i, (s, w) in enumerate(spans)}


# ----------------------------------------------------------------------------- operon_flags.tsv


def test_operon_flags_header_order_and_cells(tmp_path: Path) -> None:
    out = tmp_path / "operon_flags.tsv"
    _writers().write_operon_flags_tsv(out, flags_frame())
    lines = _lines(out)
    assert lines[0] == list(FLAG_COLUMNS)
    rows = [dict(zip(FLAG_COLUMNS, r, strict=True)) for r in lines[1:]]
    assert [r["pair_id"] for r in rows] == ["PAIR_3", "PAIR_12"]
    assert rows[0]["promoter_between_yes"] == "true"
    assert rows[0]["terminator_between_yes"] == "false"
    assert rows[0]["q_value"] == "0.04" and rows[1]["q_value"] == "NA"
    assert rows[1]["lfdr"] == "NA"
    assert rows[1]["flag"] == "coabundant_despite_terminator"


def test_operon_flags_empty_frame_is_header_only_and_valid(tmp_path: Path) -> None:
    """Degenerate case "empty operon_flags.tsv": header-only, valid against its schema."""
    out = tmp_path / "operon_flags.tsv"
    _writers().write_operon_flags_tsv(out, _empty(FLAG_COLUMNS))
    assert out.read_text(encoding="utf-8") == "\t".join(FLAG_COLUMNS) + "\n"
    report = _validate(out, "operon_flags")
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


# ----------------------------------------------------------------------------- schemas


@pytest.mark.parametrize(
    ("stem", "cols"),
    [("pairs", PAIR_COLUMNS), ("operons", OPERON_COLUMNS), ("operon_flags", FLAG_COLUMNS)],
)
def test_schema_fields_and_coordinates(stem: str, cols: tuple[str, ...]) -> None:
    d = _schema(stem)
    assert [f["name"] for f in d["fields"]] == list(cols)
    assert d.get("coordinates") == "1-based-closed"
    assert "NA" in d.get("missingValues", [])
    for f in d["fields"]:
        assert f.get("description"), f"{stem}.{f['name']}: no description"


def test_schema_class_and_test_enums_are_the_9_3_ones() -> None:
    assert set(_field("pairs", "class")["constraints"]["enum"]) == SPEC_CLASS
    assert set(_field("pairs", "test")["constraints"]["enum"]) == SPEC_TEST


def test_schema_flag_enum_is_the_9_3_one() -> None:
    enum = _field("operon_flags", "flag")["constraints"]["enum"]
    assert set(enum) == SPEC_FLAG
    assert set(FLAGS) == SPEC_FLAG  # the M6-07 code and §9.3 agree


@pytest.mark.parametrize(
    ("stem", "col"),
    [("pairs", "detected_a"), ("pairs", "detected_b"),
     ("operon_flags", "promoter_between_yes"), ("operon_flags", "terminator_between_yes")],
)  # fmt: skip
def test_schema_booleans(stem: str, col: str) -> None:
    assert _field(stem, col)["type"] == "boolean"


@pytest.mark.parametrize(
    ("stem", "col"),
    [("pairs", "intergenic_bp"), ("pairs", "n_perm"), ("operons", "start"), ("operons", "end"),
     ("operons", "n_genes"), ("operons", "n_flags")],
)  # fmt: skip
def test_schema_integers(stem: str, col: str) -> None:
    assert _field(stem, col)["type"] == "integer"


def test_example_6_written_tables_validate(tmp_path: Path) -> None:
    """Example 6: ``frictionless validate`` passes on the three resources."""
    w = _writers()
    cases = [
        (w.write_pairs_tsv, pairs_frame().drop(columns="which"), "pairs"),
        (w.write_operons_tsv, operons_frame(), "operons"),
        (w.write_operon_flags_tsv, flags_frame(), "operon_flags"),
    ]
    for write, df, stem in cases:
        out = tmp_path / f"{stem}.tsv"
        write(out, df)
        report = _validate(out, stem)
        assert report.valid, (stem, report.flatten(["rowNumber", "fieldName", "type", "note"]))


def _one_pair_with_flags(flags: list[str]) -> pd.DataFrame:
    df = pairs_frame().drop(columns="which").iloc[[1]].reset_index(drop=True)
    df["flags"] = pd.Series([flags], dtype=object)
    return df


@pytest.mark.parametrize("token", OPERON_TOKENS)
def test_pairs_flag_vocabulary_accepts_every_operon_token(tmp_path: Path, token: str) -> None:
    out = tmp_path / "pairs.tsv"
    _writers().write_pairs_tsv(out, _one_pair_with_flags([token, "rna_between"]))
    report = _validate(out, "pairs")
    assert report.valid, (token, report.flatten(["rowNumber", "fieldName", "type", "note"]))


@pytest.mark.parametrize("token", ["uncalibrated", "not_a_flag", "layer_missing:"])
def test_pairs_flag_vocabulary_rejects_unknown_tokens(tmp_path: Path, token: str) -> None:
    out = tmp_path / "pairs.tsv"
    try:
        _writers().write_pairs_tsv(out, _one_pair_with_flags(["rna_between", token]))
    except ValueError:
        return  # the writer refusing an unknown flag also keeps the vocabulary closed
    report = _validate(out, "pairs")
    assert not report.valid, f"pairs.tsv flags accepted {token!r}"


def test_decision_7_pairs_flags_documents_both_imputation_tokens() -> None:
    """Decision 7: ``abundance_imputed`` is the pair-scope token; the schema documents it and
    the distinct gene-scope ``imputed_abundance``."""
    text = json.dumps(_field("pairs", "flags"))
    assert "abundance_imputed" in text
    assert "imputed_abundance" in text


def test_decision_7_vocabulary_exact_in_one_cell(tmp_path: Path) -> None:
    """All decision 7 tokens together in one ``;``-joined cell validate."""
    out = tmp_path / "pairs.tsv"
    _writers().write_pairs_tsv(out, _one_pair_with_flags(list(OPERON_TOKENS)))
    report = _validate(out, "pairs")
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def test_schema_mentions_no_uncalibrated_anywhere() -> None:
    for stem in ("pairs", "operons", "operon_flags"):
        text = (SCHEMA_DIR / f"{stem}.schema.json").read_text(encoding="utf-8")
        assert "uncalibrated" not in text, stem


# ----------------------------------------------------------------------------- properties


_ids = st.lists(st.integers(1, 10_000), min_size=0, max_size=12, unique=True)


@settings(
    max_examples=60, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(nums=_ids, seed=st.integers(0, 2**32 - 1))
def test_property_flag_rows_sorted_by_pair_number(
    tmp_path: Path, nums: list[int], seed: int
) -> None:
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(nums)) if nums else np.array([], dtype=int)
    shuffled = [nums[i] for i in order]
    n = len(shuffled)
    df = pd.DataFrame(
        {
            "pair_id": [f"PAIR_{k}" for k in shuffled],
            "gene_a": [f"A{k}" for k in shuffled],
            "gene_b": [f"B{k}" for k in shuffled],
            "uniop_class": ["OP"] * n,
            "promoter_between_yes": [True] * n,
            "terminator_between_yes": [False] * n,
            "abundance_pattern": ["downstream_higher"] * n,
            "flag": ["internal_promoter_supported"] * n,
            "q_value": [math.nan] * n,
            "lfdr": [math.nan] * n,
            "evidence_note": ["promoter between"] * n,
        },
        columns=list(FLAG_COLUMNS),
    )
    out = tmp_path / "flags_prop.tsv"
    _writers().write_operon_flags_tsv(out, df)
    got = [r[0] for r in _lines(out)[1:]]
    assert got == [f"PAIR_{k}" for k in sorted(nums)]


_word = st.text(alphabet=string.ascii_lowercase + "_", min_size=1, max_size=8)


@settings(
    max_examples=60, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(genes=st.lists(st.lists(_word, min_size=1, max_size=5), min_size=1, max_size=5))
def test_property_operon_gene_lists_round_trip(tmp_path: Path, genes: list[list[str]]) -> None:
    n = len(genes)
    df = operons_frame().iloc[[0] * n].reset_index(drop=True)
    df["tu_id"] = [f"TU_{i + 1}" for i in range(n)]
    df["start"] = np.arange(n, dtype=np.int64) * 1000
    df["end"] = df["start"] + 500
    df["genes"] = pd.Series(genes, dtype=object)
    df["n_genes"] = np.array([len(g) for g in genes], dtype=np.int64)
    out = tmp_path / "operons_genes.tsv"
    _writers().write_operons_tsv(out, df)
    i_g = OPERON_COLUMNS.index("genes")
    assert [r[i_g].split(";") for r in _lines(out)[1:]] == genes
