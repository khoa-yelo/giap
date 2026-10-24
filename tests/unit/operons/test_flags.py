"""Evidence "yes", abundance patterns, the seven flags and ``operon_flags.tsv`` (issue M6-07;
GIAP_SPEC.md §7.4.5, §7.4.7, §9.3; docs/specs/M6-07-flags-and-operons-table.md).

Every row of the §7.4.5 flag table, spec examples 1-10 and 12, the full 2 x 2 x 5 x 2 truth table
and the degenerate cases, on hand-built merged pair frames (tests/helpers/flags_data.py), with
the lead decisions of 2026-10-21 (Q1: rows 3 and 5 still fire without a terminator layer, with
``layer_missing:terminators``; decision 5: ``ci_unestimated`` is a tested row, ``unestimated``
never makes evidence "yes"; decision 6: imputed abundances add the label ``abundance_imputed``).

Readings this file pins (the spec leaves them open):

- ``flag_for(cls, promoter_yes, terminator_yes, pattern) -> str | None`` is the public pure
  function behind the table (the property "the flag is a pure function of (class, promoter_yes,
  terminator_yes, pattern)" needs a handle), and ``PATTERNS`` lists the six pattern values.
- A pattern is one value per pair: ``downstream_higher`` / ``upstream_higher`` when the
  discordant OP pair has a direction, ``discordant`` only when it has none (Delta = 0). Rows 3
  and 7 ("discordant") accept all three; row 7 also takes ``coabundant``.
- The input of ``flag_pairs`` is the *merged* pair frame: ``merge_tests(pairs, discordance_rows,
  missed_rows, genes, *, missing_layers=())`` writes the M6-05/06 test columns into the M6-04
  frame (``which`` dropped, lead decision 4), adds ``abundance_imputed`` (decision 6) and, per
  missing layer, ``layer_missing:<layer>`` to every OP and NO row. Per-pair NA evidence with the
  layer present is "no call" (M2-15's ``terminator_between`` returns None for an empty window), so
  a missing layer must be named by the caller rather than guessed from NA.
- ``evidence_note`` names the numbers it used and says "missing" (never "no terminator") when
  the row carries ``layer_missing:<layer>``.
- ``operon_flags.tsv`` is sorted by the number in ``pair_id`` (``PAIR_2`` before ``PAIR_10``), the
  ``pairs.tsv`` order.
"""

from __future__ import annotations

import itertools
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from helpers.flags_data import (
    NAN,
    as_lists,
    cfg,
    frame,
    gene_table,
    row,
    rows_frame,
)
from helpers.pairs_data import PAIR_COLS

from giap.operons import flags as fl

ROOT = Path(__file__).resolve().parents[3]
FLAGS = (
    "internal_promoter_supported", "internal_terminator_supported", "unexplained_discordance",
    "conflict_terminator_in_operon", "novel_operon_candidate", "coabundant_despite_terminator",
    "other",
)  # fmt: skip
FLAG_COLUMNS = (
    "pair_id", "gene_a", "gene_b", "uniop_class", "promoter_between_yes",
    "terminator_between_yes", "abundance_pattern", "flag", "q_value", "lfdr", "evidence_note",
)  # fmt: skip
PATTERNS = (
    "downstream_higher",
    "upstream_higher",
    "discordant",
    "concordant",
    "coabundant",
    "none",
)
DISCORDANT = {"downstream_higher", "upstream_higher", "discordant"}


def _spec_line(prefix: str) -> str:
    for line in (ROOT / "GIAP_SPEC.md").read_text(encoding="utf-8").splitlines():
        if line.startswith(prefix):
            return line
    raise AssertionError(f"{prefix!r} not found in GIAP_SPEC.md")


def reference_flag(cls: str, p: bool, t: bool, pattern: str) -> str | None:
    """The §7.4.5 table as written in docs/specs/M6-07 step 3, first match wins."""
    table = [
        ("OP", True, False, {"downstream_higher"}, "internal_promoter_supported"),
        ("OP", False, True, {"upstream_higher"}, "internal_terminator_supported"),
        ("OP", False, False, DISCORDANT, "unexplained_discordance"),
        ("OP", None, True, {"concordant"}, "conflict_terminator_in_operon"),
        ("NO", False, False, {"coabundant"}, "novel_operon_candidate"),
        ("NO", None, True, {"coabundant"}, "coabundant_despite_terminator"),
    ]
    for c, pp, tt, pats, flag in table:
        if c == cls and (pp is None or pp == p) and tt == t and pattern in pats:
            return flag
    if cls in ("OP", "NO") and pattern in DISCORDANT | {"coabundant"}:
        return "other"
    return None


def flag_of(*rows: dict, c=None) -> list[str]:
    out = fl.flag_pairs(frame(*rows), c or cfg())
    return list(out["flag"])


# ---------------------------------------------------------------------- constants and schema


def test_flags_constant_is_exact() -> None:
    assert tuple(fl.FLAGS) == FLAGS


def test_flags_match_spec_9_3_enum() -> None:
    line = _spec_line("**`tables/operon_flags.tsv`:**")
    enum = re.search(r"`flag` enum\{([^}]*)\}", line)
    assert enum is not None
    assert tuple(s.strip() for s in enum.group(1).split(",")) == tuple(fl.FLAGS)


def test_flag_columns_match_spec_9_3() -> None:
    line = re.sub(r"enum\{[^}]*\}", "", _spec_line("**`tables/operon_flags.tsv`:**"))
    cols = tuple(re.findall(r"`([^`]+)`", line))[1:]
    assert cols == FLAG_COLUMNS
    assert tuple(fl.FLAG_COLUMNS) == FLAG_COLUMNS


def test_patterns_constant() -> None:
    assert tuple(fl.PATTERNS) == PATTERNS


# ------------------------------------------------------------------------------ truth table


def test_truth_table_2x2x5x2_matches_the_spec_table() -> None:
    pats = [p for p in PATTERNS if p != "none"]
    combos = list(itertools.product(("OP", "NO"), (False, True), (False, True), pats))
    assert len(combos) == 40
    for cls, p, t, pat in combos:
        assert fl.flag_for(cls, p, t, pat) == reference_flag(cls, p, t, pat), (cls, p, t, pat)


@pytest.mark.parametrize("cls", ["OP", "NO", "CV", "DV", "SW"])
def test_pattern_none_never_flags(cls: str) -> None:
    for p, t in itertools.product((False, True), repeat=2):
        assert fl.flag_for(cls, p, t, "none") is None


@pytest.mark.parametrize("cls", ["CV", "DV", "SW"])
def test_cv_dv_sw_never_flag_whatever_the_pattern(cls: str) -> None:
    for p, t, pat in itertools.product((False, True), (False, True), PATTERNS):
        assert fl.flag_for(cls, p, t, pat) is None


@pytest.mark.parametrize(
    ("args", "flag"),
    [
        (("OP", True, False, "downstream_higher"), "internal_promoter_supported"),
        (("OP", False, True, "upstream_higher"), "internal_terminator_supported"),
        (("OP", False, False, "discordant"), "unexplained_discordance"),
        (("OP", False, False, "downstream_higher"), "unexplained_discordance"),
        (("OP", False, False, "upstream_higher"), "unexplained_discordance"),
        (("OP", False, True, "concordant"), "conflict_terminator_in_operon"),
        (("OP", True, True, "concordant"), "conflict_terminator_in_operon"),
        (("NO", False, False, "coabundant"), "novel_operon_candidate"),
        (("NO", False, True, "coabundant"), "coabundant_despite_terminator"),
        (("NO", True, True, "coabundant"), "coabundant_despite_terminator"),
        (("OP", True, True, "downstream_higher"), "other"),
        (("OP", True, False, "upstream_higher"), "other"),
        (("NO", True, False, "coabundant"), "other"),
        (("OP", True, False, "discordant"), "other"),
        (("OP", True, False, "concordant"), None),
        (("OP", False, False, "concordant"), None),
    ],
)
def test_flag_table_row_by_row(args: tuple, flag: str | None) -> None:
    assert fl.flag_for(*args) == flag


# --------------------------------------------------------------------------- spec examples

EXAMPLES = {
    1: (row("OP", prom=0.81, term=0.10, delta=1.4, q=0.02), "internal_promoter_supported"),
    2: (row("OP", prom=0.10, term=0.80, delta=-1.4, q=0.02), "internal_terminator_supported"),
    3: (row("OP", prom=0.10, term=0.10, delta=1.4, q=0.02), "unexplained_discordance"),
    4: (row("OP", prom=0.10, term=0.80, delta=0.2, q=0.60), "conflict_terminator_in_operon"),
    5: (row("NO", prom=0.10, term=0.10, delta=0.01, q=0.04, lfdr=0.3), "novel_operon_candidate"),
    6: (row("NO", prom=0.10, term=0.80, delta=0.01, q=0.04, lfdr=0.3),
        "coabundant_despite_terminator"),
    7: (row("OP", prom=0.81, term=0.80, delta=1.4, q=0.02), "other"),
    8: (row("OP", prom=0.81, term=0.10, delta=-1.4, q=0.02), "other"),
    9: (row("NO", prom=0.81, term=0.10, delta=0.01, q=0.04, lfdr=0.3), "other"),
    12: (row("OP", prom=0.81, term=0.10, q=0.03, detected=(False, True)),
         "internal_promoter_supported"),
}  # fmt: skip


@pytest.mark.parametrize("n", sorted(EXAMPLES))
def test_spec_examples(n: int) -> None:
    r, flag = EXAMPLES[n]
    assert flag_of(r) == [flag]


def test_example_6_any_promoter() -> None:
    r = row("NO", prom=0.81, term=0.80, delta=0.01, q=0.04, lfdr=0.3)
    assert flag_of(r) == ["coabundant_despite_terminator"]


@pytest.mark.parametrize("term", [NAN, 0.10])
def test_example_10_concordant_without_terminator_gets_no_row(term: float) -> None:
    out = fl.flag_pairs(frame(row("OP", prom=0.10, term=term, delta=0.2, q=0.60)), cfg())
    assert len(out) == 0
    assert tuple(out.columns) == FLAG_COLUMNS


def test_example_12_pattern_is_downstream_higher() -> None:
    out = fl.flag_pairs(frame(EXAMPLES[12][0]), cfg())
    assert list(out["abundance_pattern"]) == ["downstream_higher"]


def test_one_detected_mirror_upstream_detected() -> None:
    r = row("OP", prom=0.10, term=0.80, q=0.03, detected=(True, False))
    out = fl.flag_pairs(frame(r), cfg())
    assert list(out["abundance_pattern"]) == ["upstream_higher"]
    assert list(out["flag"]) == ["internal_terminator_supported"]


def test_row_3_with_upstream_higher() -> None:
    assert flag_of(row("OP", prom=0.10, term=0.10, delta=-1.4, q=0.02)) == [
        "unexplained_discordance"
    ]


def test_one_detected_cannot_be_conflict_terminator_in_operon() -> None:
    for det in ((False, True), (True, False)):
        r = row("OP", prom=0.10, term=0.80, q=0.60, detected=det)
        assert flag_of(r) == []


# ------------------------------------------------------------------------------ evidence_yes


def test_evidence_yes_columns_dtype_and_index() -> None:
    df = frame(row("OP", prom=0.81, term=0.10), row("NO", prom=0.10, term=0.80))
    df.index = [7, 3]
    ev = fl.evidence_yes(df, cfg())
    assert list(ev.columns) == ["promoter_between_yes", "terminator_between_yes"]
    assert list(ev.index) == [7, 3]
    assert ev["promoter_between_yes"].dtype == bool
    assert ev["terminator_between_yes"].dtype == bool
    assert list(ev["promoter_between_yes"]) == [True, False]
    assert list(ev["terminator_between_yes"]) == [False, True]


def test_evidence_yes_cut_off_is_inclusive() -> None:
    df = frame(row("OP", prom=0.5, term=0.5), row("OP", prom=0.4999, term=0.4999))
    ev = fl.evidence_yes(df, cfg())
    assert list(ev["promoter_between_yes"]) == [True, False]
    assert list(ev["terminator_between_yes"]) == [True, False]


def test_evidence_yes_reads_config_cut_offs() -> None:
    df = frame(row("OP", prom=0.81, term=0.80))
    ev = fl.evidence_yes(df, cfg(promotech_yes=0.9, terminatornet_yes=0.85))
    assert not ev["promoter_between_yes"].iloc[0]
    assert not ev["terminator_between_yes"].iloc[0]
    ev = fl.evidence_yes(df, cfg(promotech_yes=0.8, terminatornet_yes=0.8))
    assert ev["promoter_between_yes"].iloc[0]
    assert ev["terminator_between_yes"].iloc[0]


def test_evidence_yes_na_is_false() -> None:
    ev = fl.evidence_yes(frame(row("OP"), row("CV"), row("DV")), cfg())
    assert not ev.to_numpy().any()


def test_unestimated_row_never_makes_evidence_yes() -> None:
    """Lead decision 5: a row flagged ``unestimated`` never makes evidence "yes"."""
    df = frame(
        row("OP", prom=0.9, term=0.9, flags=["no_detection_limit_model", "unestimated"]),
        row("OP", prom=0.9, term=0.9, flags=["too_few_op_pairs", "unestimated"]),
    )
    ev = fl.evidence_yes(df, cfg())
    assert not ev.to_numpy().any()


def test_ci_unestimated_is_a_tested_row() -> None:
    """``ci_unestimated`` is not ``unestimated`` (a substring match would be wrong)."""
    df = frame(row("OP", prom=0.9, term=0.9, delta=1.4, q=0.02, flags=["ci_unestimated"]))
    ev = fl.evidence_yes(df, cfg())
    assert ev.to_numpy().all()


def test_flags_as_semicolon_string_are_read() -> None:
    df = frame(row("OP", prom=0.9, term=0.9))
    df["flags"] = ["too_few_op_pairs;unestimated"]
    assert not fl.evidence_yes(df, cfg()).to_numpy().any()


# ------------------------------------------------------------------------- abundance_pattern


def _pat(*rows: dict, alpha: float = 0.1) -> list[str]:
    return list(fl.abundance_pattern(frame(*rows), alpha))


@pytest.mark.parametrize(
    ("r", "expected"),
    [
        (row("OP", delta=1.4, q=0.02), "downstream_higher"),
        (row("OP", delta=-1.4, q=0.02), "upstream_higher"),
        (row("OP", delta=0.0, q=0.02), "discordant"),
        (row("OP", q=0.03, detected=(False, True)), "downstream_higher"),
        (row("OP", q=0.03, detected=(True, False)), "upstream_higher"),
        (row("OP", delta=0.3, q=0.60), "concordant"),
        (row("OP", delta=1.4, q=0.10), "downstream_higher"),  # q <= alpha is inclusive
        (row("OP", delta=0.3, q=0.1000001), "concordant"),
        (row("OP", delta=0.3), "none"),  # untested: q NA is not concordant
        (row("OP", q=0.60, detected=(False, True)), "none"),  # concordant needs both detected
        (row("OP", detected=(False, False)), "none"),
        (row("NO", delta=0.01, q=0.04), "coabundant"),
        (row("NO", delta=0.01, q=0.10), "coabundant"),
        (row("NO", delta=0.01, q=0.30), "none"),
        (row("NO", delta=0.01), "none"),
        (row("CV", delta=0.01, q=0.01), "none"),
        (row("DV", delta=0.01, q=0.01), "none"),
        (row("SW", delta=0.01, q=0.01), "none"),
    ],
)
def test_abundance_pattern_rules(r: dict, expected: str) -> None:
    assert _pat(r) == [expected]


def test_abundance_pattern_reads_alpha() -> None:
    rs = (row("OP", delta=1.4, q=0.06), row("NO", delta=0.01, q=0.06))
    assert _pat(*rs, alpha=0.1) == ["downstream_higher", "coabundant"]
    assert _pat(*rs, alpha=0.05) == ["concordant", "none"]


def test_unestimated_row_has_pattern_none_even_with_a_q() -> None:
    r = row("OP", delta=1.4, q=0.02, flags=["too_few_op_pairs", "unestimated"])
    assert _pat(r) == ["none"]


def test_abundance_pattern_series_index_and_values() -> None:
    df = frame(row("OP", delta=1.4, q=0.02), row("CV"))
    df.index = [5, 1]
    s = fl.abundance_pattern(df, 0.1)
    assert isinstance(s, pd.Series)
    assert list(s.index) == [5, 1]
    assert set(s) <= set(PATTERNS)


# --------------------------------------------------------------------------------- flag_pairs


def test_flag_pairs_columns_and_values() -> None:
    df = frame(EXAMPLES[1][0], EXAMPLES[5][0], row("CV"))
    out = fl.flag_pairs(df, cfg())
    assert tuple(out.columns) == FLAG_COLUMNS
    assert list(out["pair_id"]) == ["PAIR_1", "PAIR_2"]
    assert list(out["gene_a"]) == ["G1a", "G2a"]
    assert list(out["gene_b"]) == ["G1b", "G2b"]
    assert list(out["uniop_class"]) == ["OP", "NO"]
    assert [bool(v) for v in out["promoter_between_yes"]] == [True, False]
    assert [bool(v) for v in out["terminator_between_yes"]] == [False, False]
    assert all(isinstance(v, bool | np.bool_) for v in out["promoter_between_yes"])
    assert list(out["abundance_pattern"]) == ["downstream_higher", "coabundant"]
    assert list(out["flag"]) == ["internal_promoter_supported", "novel_operon_candidate"]
    assert out["q_value"].tolist() == pytest.approx([0.02, 0.04])
    assert math.isnan(out["lfdr"].iloc[0])
    assert out["lfdr"].iloc[1] == pytest.approx(0.3)


def test_flag_pairs_sorted_by_pair_number() -> None:
    rs = [row("OP", prom=0.1, term=0.1, delta=1.4, q=0.02) for _ in range(12)]
    df = frame(*rs)
    shuffled = df.iloc[np.random.default_rng(3).permutation(len(df))].reset_index(drop=True)
    out = fl.flag_pairs(shuffled, cfg())
    assert list(out["pair_id"]) == [f"PAIR_{k}" for k in range(1, 13)]
    pd.testing.assert_frame_equal(
        out.reset_index(drop=True), fl.flag_pairs(df, cfg()).reset_index(drop=True)
    )


@pytest.mark.parametrize("cls", ["CV", "DV", "SW"])
def test_cv_dv_sw_rows_never_flagged(cls: str) -> None:
    for prom, term in itertools.product((0.1, 0.9), repeat=2):
        assert flag_of(row(cls, prom=prom, term=term, delta=0.01, q=0.01)) == []


def test_every_test_unestimated_gives_header_only() -> None:
    """sim100k: ``too_few_op_pairs`` / ``too_few_no_pairs`` everywhere -> header-only table."""
    rs = [
        row("OP", prom=0.9, term=0.9, delta=1.4, flags=["too_few_op_pairs", "unestimated"]),
        row("OP", prom=0.1, term=0.1, flags=["no_detection_limit_model", "unestimated"],
            detected=(False, True)),
        row("NO", prom=0.1, term=0.1, delta=0.01, flags=["too_few_no_pairs"]),
        row("CV"), row("SW"),
    ]  # fmt: skip
    out = fl.flag_pairs(frame(*rs), cfg())
    assert len(out) == 0
    assert tuple(out.columns) == FLAG_COLUMNS


def test_empty_pair_frame() -> None:
    out = fl.flag_pairs(frame(), cfg())
    assert len(out) == 0
    assert tuple(out.columns) == FLAG_COLUMNS
    assert fl.evidence_yes(frame(), cfg()).shape == (0, 2)
    assert len(fl.abundance_pattern(frame(), 0.1)) == 0


def test_abundance_imputed_label_does_not_change_the_flag() -> None:
    r = row("NO", prom=0.1, term=0.1, delta=0.01, q=0.04, flags=["abundance_imputed"])
    assert flag_of(r) == ["novel_operon_candidate"]


def test_flag_pairs_does_not_mutate_and_is_deterministic() -> None:
    df = frame(*(r for r, _ in EXAMPLES.values()), row("CV"), row("OP", flags=["unestimated"]))
    before = df.copy(deep=True)
    a = fl.flag_pairs(df, cfg())
    b = fl.flag_pairs(df, cfg())
    pd.testing.assert_frame_equal(df, before)
    pd.testing.assert_frame_equal(a, b)
    fl.evidence_yes(df, cfg())
    fl.abundance_pattern(df, 0.1)
    pd.testing.assert_frame_equal(df, before)


def test_flagged_rows_keep_pair_q_values() -> None:
    df = frame(*(r for r, _ in EXAMPLES.values()))
    out = fl.flag_pairs(df, cfg())
    by_id = df.set_index("pair_id")
    for pid, q in zip(out["pair_id"], out["q_value"], strict=True):
        assert by_id.loc[pid, "q_value"] == pytest.approx(q)
    assert "none" not in set(out["abundance_pattern"])
    assert set(out["flag"]) <= set(FLAGS)


# ------------------------------------------------------------------------------ evidence_note


def test_evidence_note_names_the_evidence_used() -> None:
    out = fl.flag_pairs(frame(EXAMPLES[1][0]), cfg())
    note = out["evidence_note"].iloc[0]
    assert isinstance(note, str) and note
    for token in ("0.81", "0.02", "1.4"):
        assert token in note, note


def test_evidence_note_is_a_nonempty_string_for_every_row() -> None:
    out = fl.flag_pairs(frame(*(r for r, _ in EXAMPLES.values())), cfg())
    assert len(out) == len(EXAMPLES)
    for note in out["evidence_note"]:
        assert isinstance(note, str) and note.strip()
        assert "absent" not in note.lower() and "absence" not in note.lower()


# --------------------------------------------------------------- layers missing (lead Q1)

TERM_MISSING = ["layer_missing:terminators"]
PROM_MISSING = ["layer_missing:promoters"]


def test_no_terminator_layer_rows_3_and_5_still_fire() -> None:
    rs = (
        row("OP", prom=0.1, delta=1.4, q=0.02, flags=["ci_unestimated", *TERM_MISSING]),
        row("NO", prom=0.1, delta=0.01, q=0.04, lfdr=0.3, flags=TERM_MISSING),
    )
    assert flag_of(*rs) == ["unexplained_discordance", "novel_operon_candidate"]


def test_no_terminator_layer_rows_2_4_6_impossible() -> None:
    rs = (
        row("OP", prom=0.1, delta=-1.4, q=0.02, flags=["ci_unestimated", *TERM_MISSING]),
        row("OP", prom=0.1, delta=0.2, q=0.6, flags=["ci_unestimated", *TERM_MISSING]),
        row("NO", prom=0.81, delta=0.01, q=0.04, flags=TERM_MISSING),
    )
    out = fl.flag_pairs(frame(*rs), cfg())
    assert list(out["pair_id"]) == ["PAIR_1", "PAIR_3"]
    assert list(out["flag"]) == ["unexplained_discordance", "other"]
    assert not out["terminator_between_yes"].astype(bool).any()


def test_no_terminator_layer_is_said_in_evidence_note() -> None:
    rs = (
        row("OP", prom=0.1, delta=1.4, q=0.02, flags=["ci_unestimated", *TERM_MISSING]),
        row("NO", prom=0.1, delta=0.01, q=0.04, flags=TERM_MISSING),
    )
    for note in fl.flag_pairs(frame(*rs), cfg())["evidence_note"]:
        low = note.lower()
        assert "terminator" in low and "missing" in low, note
        assert "no terminator" not in low, note


def test_no_promoter_layer_never_internal_promoter_supported() -> None:
    r = row("OP", term=0.1, delta=1.4, q=0.02, flags=["ci_unestimated", *PROM_MISSING])
    out = fl.flag_pairs(frame(r), cfg())
    assert list(out["flag"]) == ["unexplained_discordance"]
    note = out["evidence_note"].iloc[0].lower()
    assert "promoter" in note and "missing" in note and "no promoter" not in note


# -------------------------------------------------------------------------------- merge_tests


def _base_pairs() -> pd.DataFrame:
    return frame(
        row("OP", prom=0.81, term=0.1, delta=1.4, flags=[]),
        row("NO", prom=0.1, term=0.1, delta=0.01, flags=["rna_between"]),
        row("CV", delta=0.3),
        row("DV", delta=0.3),
        row("SW", prom=0.2, term=0.2, delta=0.3),
        row("OP", delta=0.2, flags=[]),
    )


def _disc_rows() -> pd.DataFrame:
    return rows_frame(
        [
            {"pair_id": "PAIR_1", "test": "both_detected_t", "stat": 1.4, "effect": 1.3,
             "p_value": 0.01, "q_value": 0.02, "method": "student_t",
             "null_model": "within_operon_t", "flags": ["ci_unestimated"]},
            {"pair_id": "PAIR_6", "test": "both_detected_t", "stat": 0.2, "effect": 0.1,
             "p_value": 0.5, "q_value": 0.6, "method": "student_t",
             "null_model": "within_operon_t", "flags": ["ci_unestimated"]},
        ]
    )  # fmt: skip


def _missed_rows() -> pd.DataFrame:
    return rows_frame(
        [
            {"pair_id": "PAIR_2", "test": "missed_operon_decoy", "stat": -0.01, "effect": 0.01,
             "q_value": 0.04, "lfdr": 0.3, "method": "decoy_fdr",
             "null_model": "convergent+sandwich_reweighted", "which": "sandwich",
             "flags": ["rna_between"]},
        ],
        missed=True,
    )  # fmt: skip


def _genes(pairs: pd.DataFrame, imputed: tuple[str, ...] = ()) -> pd.DataFrame:
    tags = list(dict.fromkeys([*pairs["gene_a"], *pairs["gene_b"]]))
    return gene_table({t: 10.0 for t in tags}, tags, imputed=imputed)


def test_merge_tests_fills_test_columns_and_drops_which() -> None:
    pairs = _base_pairs()
    out = fl.merge_tests(pairs, _disc_rows(), _missed_rows(), _genes(pairs))
    assert tuple(out.columns) == PAIR_COLS
    assert "which" not in out.columns
    assert list(out["pair_id"]) == list(pairs["pair_id"])
    by = out.set_index("pair_id")
    assert by.loc["PAIR_1", "test"] == "both_detected_t"
    assert by.loc["PAIR_1", "q_value"] == pytest.approx(0.02)
    assert by.loc["PAIR_1", "p_value"] == pytest.approx(0.01)
    assert by.loc["PAIR_1", "method"] == "student_t"
    assert by.loc["PAIR_6", "q_value"] == pytest.approx(0.6)
    assert by.loc["PAIR_2", "test"] == "missed_operon_decoy"
    assert by.loc["PAIR_2", "q_value"] == pytest.approx(0.04)
    assert by.loc["PAIR_2", "lfdr"] == pytest.approx(0.3)
    assert math.isnan(by.loc["PAIR_2", "p_value"])
    for pid in ("PAIR_3", "PAIR_4", "PAIR_5"):
        assert by.loc[pid, "test"] == "none"
        assert math.isnan(by.loc[pid, "q_value"])
    # the pair columns themselves are untouched
    for c in ("class", "gene_a", "gene_b", "promoter_between", "delta_log2"):
        pd.testing.assert_series_equal(out[c], pairs[c], check_names=False)


def test_merge_tests_flags_are_a_union_without_duplicates() -> None:
    pairs = _base_pairs()
    out = fl.merge_tests(pairs, _disc_rows(), _missed_rows(), _genes(pairs))
    flags = dict(zip(out["pair_id"], as_lists(out), strict=True))
    assert "ci_unestimated" in flags["PAIR_1"]
    assert flags["PAIR_2"].count("rna_between") == 1
    for f in flags.values():
        assert len(f) == len(set(f))
    assert flags["PAIR_3"] == []


def test_merge_tests_abundance_imputed_label() -> None:
    """Lead decision 6: a pair whose gene has ``abundance_imputed`` gains that flag."""
    pairs = _base_pairs()
    out = fl.merge_tests(pairs, _disc_rows(), _missed_rows(), _genes(pairs, ("G2b", "G1a")))
    flags = dict(zip(out["pair_id"], as_lists(out), strict=True))
    assert "abundance_imputed" in flags["PAIR_1"]
    assert "abundance_imputed" in flags["PAIR_2"]
    for pid in ("PAIR_3", "PAIR_4", "PAIR_5", "PAIR_6"):
        assert "abundance_imputed" not in flags[pid]


def test_merge_tests_without_imputed_column() -> None:
    pairs = _base_pairs()
    genes = _genes(pairs).drop(columns=["abundance_imputed"])
    out = fl.merge_tests(pairs, _disc_rows(), _missed_rows(), genes)
    assert not any("abundance_imputed" in f for f in as_lists(out))


def test_merge_tests_layer_missing_on_op_and_no_rows_only() -> None:
    pairs = _base_pairs()
    out = fl.merge_tests(
        pairs, _disc_rows(), _missed_rows(), _genes(pairs), missing_layers=("terminators",)
    )
    flags = dict(zip(out["pair_id"], as_lists(out), strict=True))
    for pid in ("PAIR_1", "PAIR_2", "PAIR_6"):
        assert "layer_missing:terminators" in flags[pid]
        assert "layer_missing:promoters" not in flags[pid]
    for pid in ("PAIR_3", "PAIR_4"):
        assert not any(f.startswith("layer_missing:") for f in flags[pid])


def test_merge_tests_no_layer_missing_by_default_even_with_na_evidence() -> None:
    """NA evidence with the layer present is "no call", not a missing layer."""
    pairs = _base_pairs()
    out = fl.merge_tests(pairs, _disc_rows(), _missed_rows(), _genes(pairs))
    assert not any(f.startswith("layer_missing:") for fs in as_lists(out) for f in fs)


def test_merge_tests_empty_test_rows() -> None:
    """No abundance: both families empty, every test stays ``none``."""
    pairs = _base_pairs()
    out = fl.merge_tests(pairs, rows_frame([]), rows_frame([], missed=True), _genes(pairs))
    assert tuple(out.columns) == PAIR_COLS
    assert set(out["test"]) == {"none"}
    assert out["q_value"].isna().all()
    assert len(fl.flag_pairs(out, cfg())) == 0


def test_merge_tests_does_not_mutate_and_is_deterministic() -> None:
    pairs, disc, missed = _base_pairs(), _disc_rows(), _missed_rows()
    genes = _genes(pairs, ("G1a",))
    copies = [x.copy(deep=True) for x in (pairs, disc, missed, genes)]
    a = fl.merge_tests(pairs, disc, missed, genes, missing_layers=("terminators",))
    b = fl.merge_tests(pairs, disc, missed, genes, missing_layers=("terminators",))
    for x, c in zip((pairs, disc, missed, genes), copies, strict=True):
        pd.testing.assert_frame_equal(x, c)
    pd.testing.assert_frame_equal(a, b)


def test_merge_then_flag_end_to_end() -> None:
    pairs = _base_pairs()
    merged = fl.merge_tests(pairs, _disc_rows(), _missed_rows(), _genes(pairs))
    out = fl.flag_pairs(merged, cfg())
    assert list(out["pair_id"]) == ["PAIR_1", "PAIR_2"]
    assert list(out["flag"]) == ["internal_promoter_supported", "novel_operon_candidate"]


def test_real_discordance_and_missed_rows_merge_and_flag() -> None:
    """The actual M6-05 / M6-06 frames feed M6-07 (lead decision 4)."""
    from helpers.missed_data import block, genes_for, stack
    from helpers.missed_data import cfg as mcfg

    from giap.operons.discordance import discordance
    from giap.operons.missed import missed_operons

    rng = np.random.default_rng(11)
    pairs = stack(
        block("OP", rng.normal(0.0, 0.4, 100)),
        block("OP", [6.0, -6.0]),
        block("NO", rng.normal(0.0, 1.0, 30)),
        block("CV", rng.normal(0.0, 1.0, 30)),
        block("SW", rng.normal(0.0, 1.0, 30)),
    )
    pairs["promoter_between"] = 0.1
    pairs["terminator_between"] = 0.1
    genes = genes_for(pairs)
    c = mcfg()
    d = discordance(pairs, genes, c, None, sigma_model="constant")  # the M6-05 null (M6-13)
    m = missed_operons(pairs, c, genes=genes)
    merged = fl.merge_tests(pairs, d.rows, m.rows, genes)
    assert tuple(merged.columns) == PAIR_COLS
    out = fl.flag_pairs(merged, c)
    by = merged.set_index("pair_id")
    assert "which" not in merged.columns
    discordant = set(merged.loc[(merged["class"] == "OP") & (merged["q_value"] <= 0.1), "pair_id"])
    assert {"PAIR_101", "PAIR_102"} <= discordant  # the planted |Delta| = 6 OP pairs
    op_out = out.loc[out["uniop_class"] == "OP"]
    assert set(op_out["pair_id"]) == discordant
    assert set(op_out["flag"]) == {"unexplained_discordance"}
    for pid, q in zip(out["pair_id"], out["q_value"], strict=True):
        assert by.loc[pid, "q_value"] == pytest.approx(q)
        assert by.loc[pid, "q_value"] <= 0.1
    assert set(out["uniop_class"]) <= {"OP", "NO"}
