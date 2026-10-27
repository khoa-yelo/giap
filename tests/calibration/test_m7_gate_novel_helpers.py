"""Matching helpers of the M7-13b novel-ORF gate, ``tests/calibration/_m7_gate.py`` (implementer;
docs/specs/M7-13-calibration-gate.md rules 5-6, examples 7-9, the novel rows of "Degenerate
cases"; M7-11b rule 11 side table; M7-11b statistician review, PR #256 comment 6365221592).
Unmarked: runs on every commit.

Interface pinned here (pure functions, pandas in, ints out). ``novel_orfs`` is ``novel_orfs.tsv``
as written (all strings, ``keep_default_na=False``); ``side`` is ``novel/novel_class_fdr.parquet``
as read (``q`` float64 with NaN, ``class_rule_met`` / ``called`` bool, ``reason`` str or None);
``truth`` is ``gippasim.truth.read_novel_orfs`` (``stop`` 1-based lowest base of the stop codon,
``true_start`` int, ``n_true_peptides`` int, ``rbs_planted`` bool).

* ``NOVEL_CLASSES = ("intergenic", "cds_outframe", "antisense", "nterm_extension")``;
  ``NOVEL_COUNT_KEYS = ("n_pred", "tp_p", "n_truth", "tp_r", "start_ok", "mismatch",
  "rbs_planted", "rbs_pass")``.
* ``novel_stop(strand, end, length=None) -> int``: the 1-based lowest base of the stop codon from
  the TSV ``end`` (last base of the stop in transcription direction): ``end - 2`` on ``+``,
  ``end`` on ``-``; with ``length`` wrapped into ``1..length`` (origin-crossing keys, rule 5).
* ``g3_predictions(novel_orfs, lengths=None) -> DataFrame`` with columns ``NOVEL_PRED_COLUMNS =
  ("orf_id", "class", "contig", "strand", "stop", "start_longest", "start_best_ostir",
  "rbs_supported")``: the ``called = true`` rows of the four classes; ``stop`` int; starts int or
  NA (``pd.isna``); ``rbs_supported`` bool (NA -> False).
* ``g4_predictions(novel_orfs, side, level, *, regime=None, lengths=None) -> DataFrame`` (same
  columns): side rows with ``class_rule_met`` and ``q <= level`` (NaN q never enters: pick losers
  and fallback calls are left out of the FDP), of ``regime`` when given, of the four classes by
  side ``class`` (the TSV class, so a class mismatch is counted in its own class); coordinates
  joined from ``novel_orfs`` on ``orf_id`` (a side row missing there raises ValueError).
  ``level=None``: no q filter (attack a: every candidate meeting the class rule, NaN q included).
* ``match_novel(pred, truth, lengths=None, min_peptides=2) -> {class: {NOVEL_COUNT_KEYS}}`` for
  every class of ``NOVEL_CLASSES`` (zeros when absent), rule 5: a prediction matches a truth row
  with equal ``(contig, strand, stop)`` and equal class (start ignored). ``tp_p`` = distinct
  predicted keys of the class found in that class's truth (two calls at one key: one TP);
  ``n_truth`` = truth rows with ``n_true_peptides >= min_peptides``; ``tp_r`` = those whose key is
  predicted in their class; ``mismatch`` = predictions at a key that is truth only under another
  class (an FP; that truth row stays unmatched, Q6); ``start_ok`` = TP keys whose truth
  ``true_start`` equals the first matching prediction's ``start_longest`` or
  ``start_best_ostir``; ``rbs_planted`` = TP keys whose truth row has ``rbs_planted``;
  ``rbs_pass`` = of those, the first matching prediction has ``rbs_supported``.
* ``list_targets(side, regime) -> {class: int}``: side rows of ``regime`` per ``fdr_class`` that
  entered the class list (finite ``q``, or ``reason = "lost_pick"``), rule 6 vacuity input.
* ``reason_counts(side) -> {class: {reason: int}}``: non-null ``reason`` values per side ``class``.
* ``fallback_calls(side) -> {class: int}``: ``called`` rows with a non-finite ``q`` (called by the
  >= 2-peptide fallback, not error-controlled; recorded, not gated).
"""

from __future__ import annotations

import math
from typing import Any

import _m7_gate as G
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.novel.outputs import NOVEL_ORFS_COLUMNS

TRUTH_COLS = ["orf_id", "class", "contig", "start", "end", "strand", "aa", "stop", "true_start",
              "n_true_peptides", "rbs_planted"]  # fmt: skip
SIDE_COLS = ["orf_id", "class", "fdr_class", "regime", "q", "class_fdr", "reason",
             "class_mismatch", "class_rule_met", "called"]  # fmt: skip
CLASSES = ("intergenic", "cds_outframe", "antisense", "nterm_extension")
ZERO = dict.fromkeys(("n_pred", "tp_p", "n_truth", "tp_r", "start_ok", "mismatch",
                      "rbs_planted", "rbs_pass"), 0)  # fmt: skip


def tsv(*rows: dict[str, Any]) -> pd.DataFrame:
    """A written ``novel_orfs.tsv`` (all strings) from partial rows; ``end`` and ``class`` given."""
    base = {c: "NA" for c in NOVEL_ORFS_COLUMNS}
    base.update(contig="c1", strand="+", called="true", class_fdr="0.001", rbs_supported="true",
                flags="")  # fmt: skip
    out = []
    for k, r in enumerate(rows):
        row = {**base, "orf_id": f"NOVEL_c1_{k + 1}_p", **{a: str(v) for a, v in r.items()}}
        out.append(row)
    return pd.DataFrame(out, columns=list(NOVEL_ORFS_COLUMNS), dtype=str)


def truth(*rows: tuple[str, str, int, int, int, bool]) -> pd.DataFrame:
    """Truth rows ``(class, strand, stop, true_start, n_true_peptides, rbs_planted)`` on c1."""
    data = [[f"NOVEL_{k:04d}", c, "c1", 1, 2, s, "M", stop, ts, n, rb]
            for k, (c, s, stop, ts, n, rb) in enumerate(rows)]  # fmt: skip
    return pd.DataFrame(data, columns=TRUTH_COLS)


def side(*rows: dict[str, Any]) -> pd.DataFrame:
    base = {"class": "intergenic", "regime": "picked_decoy", "q": 0.001, "reason": None,
            "class_mismatch": False, "class_rule_met": True, "called": True}  # fmt: skip
    out = []
    for k, r in enumerate(rows):
        row = {**base, "orf_id": f"NOVEL_c1_{k + 1}_p", **r}
        row.setdefault("fdr_class", row["class"])
        q = float(row["q"])
        row["class_fdr"] = f"{q:.6g}" if math.isfinite(q) else (
            "1" if row["reason"] == "lost_pick" else "unestimated")  # fmt: skip
        out.append(row)
    return pd.DataFrame(out, columns=SIDE_COLS).astype({"q": "float64", "class_mismatch": bool,
                                                        "class_rule_met": bool,
                                                        "called": bool})  # fmt: skip


def test_constants() -> None:
    assert (
        G.NOVEL_CLASSES == CLASSES == ("intergenic", "cds_outframe", "antisense", "nterm_extension")
    )
    assert tuple(ZERO) == G.NOVEL_COUNT_KEYS
    assert G.NOVEL_PRED_COLUMNS == (
        "orf_id",
        "class",
        "contig",
        "strand",
        "stop",
        "start_longest",
        "start_best_ostir",
        "rbs_supported",
    )


@pytest.mark.parametrize(("strand", "end", "length", "want"),
                         [("+", 5023, None, 5021), ("-", 5021, None, 5021), ("+", 2, 1000, 1000),
                          ("+", 1, 1000, 999), ("-", 7, 1000, 7)])  # fmt: skip
def test_novel_stop(strand: str, end: int, length: int | None, want: int) -> None:
    assert G.novel_stop(strand, end, length) == want


# --------------------------------------------------------------------------- G3 predictions


def test_g3_predictions_called_rows_of_the_four_classes() -> None:
    t = tsv({"class": "antisense", "strand": "-", "end": 5021, "start_longest": 5200,
             "start_best_ostir": "NA"},
            {"class": "intergenic", "end": 903, "called": "false"},
            {"class": "sorf", "end": 1203},
            {"class": "nterm_extension", "end": 2003, "rbs_supported": "NA"})  # fmt: skip
    p = G.g3_predictions(t)
    assert tuple(p.columns) == G.NOVEL_PRED_COLUMNS
    assert p["class"].tolist() == ["antisense", "nterm_extension"]
    assert p["stop"].tolist() == [5021, 2001]
    assert int(p["start_longest"].iloc[0]) == 5200 and pd.isna(p["start_best_ostir"].iloc[0])
    assert p["rbs_supported"].tolist() == [True, False]


# --------------------------------------------------------------------------- matching (rule 5)


def test_example7_antisense_tp_other_start() -> None:
    p = G.g3_predictions(tsv({"class": "antisense", "strand": "-", "end": 5021,
                              "start_longest": 5300, "start_best_ostir": 5240}))  # fmt: skip
    got = G.match_novel(p, truth(("antisense", "-", 5021, 5200, 3, False)))
    assert got["antisense"] == {**ZERO, "n_pred": 1, "tp_p": 1, "n_truth": 1, "tp_r": 1}
    assert set(got) == set(CLASSES)
    assert all(got[c] == ZERO for c in CLASSES if c != "antisense")


def test_start_accuracy_either_start_column() -> None:
    tr = truth(("intergenic", "+", 1201, 1000, 2, True))
    for longest, best in ((1000, 1012), (988, 1000)):
        p = G.g3_predictions(tsv({"class": "intergenic", "end": 1203, "start_longest": longest,
                                  "start_best_ostir": best}))  # fmt: skip
        got = G.match_novel(p, tr)["intergenic"]
        assert got["start_ok"] == 1 and got["rbs_planted"] == 1 and got["rbs_pass"] == 1


def test_example8_class_mismatch_is_fp_and_fn() -> None:
    p = G.g3_predictions(tsv({"class": "intergenic", "end": 3003}))
    got = G.match_novel(p, truth(("cds_outframe", "+", 3001, 2800, 4, False)))
    assert got["intergenic"] == {**ZERO, "n_pred": 1, "mismatch": 1}
    assert got["cds_outframe"] == {**ZERO, "n_truth": 1}


def test_two_calls_at_one_key_one_tp() -> None:
    p = G.g3_predictions(
        tsv(
            {"class": "cds_outframe", "end": 3003, "start_longest": 2800},
            {"class": "cds_outframe", "end": 3003, "start_longest": 2900},
        )
    )
    got = G.match_novel(p, truth(("cds_outframe", "+", 3001, 2800, 4, False)))["cds_outframe"]
    assert (got["n_pred"], got["tp_p"], got["tp_r"], got["start_ok"]) == (2, 1, 1, 1)


def test_recall_denominator_needs_two_peptides_precision_does_not() -> None:
    p = G.g3_predictions(tsv({"class": "intergenic", "end": 1203}))
    got = G.match_novel(
        p,
        truth(("intergenic", "+", 1201, 1000, 1, False), ("intergenic", "+", 4001, 3800, 2, False)),
    )["intergenic"]
    assert (got["n_pred"], got["tp_p"], got["n_truth"], got["tp_r"]) == (1, 1, 1, 0)


def test_no_prediction_contributes_truth() -> None:
    p = G.g3_predictions(tsv())
    got = G.match_novel(p, truth(("antisense", "-", 77, 300, 5, False)))
    assert got["antisense"] == {**ZERO, "n_truth": 1}


def test_minus_strand_and_origin_wrap() -> None:
    """Origin-crossing key compared modulo the replicon length (rule 5)."""
    p = G.g3_predictions(tsv({"class": "intergenic", "end": 2}), lengths={"c1": 1000})
    got = G.match_novel(p, truth(("intergenic", "+", 1000, 900, 2, False)), {"c1": 1000})
    assert got["intergenic"]["tp_p"] == 1 and got["intergenic"]["tp_r"] == 1
    other = G.match_novel(p, truth(("intergenic", "-", 1000, 1100, 2, False)), {"c1": 1000})
    assert other["intergenic"]["tp_p"] == 0  # strand differs


def test_nterm_extension_matches_by_the_extended_genes_stop() -> None:
    p = G.g3_predictions(tsv({"class": "nterm_extension", "strand": "-", "end": 8001,
                              "start_longest": 9100}))  # fmt: skip
    got = G.match_novel(p, truth(("nterm_extension", "-", 8001, 9100, 2, True)))
    assert got["nterm_extension"]["tp_p"] == 1 and got["nterm_extension"]["start_ok"] == 1


PRED = st.lists(st.tuples(st.sampled_from(CLASSES), st.sampled_from("+-"),
                          st.integers(1, 30)), max_size=12)  # fmt: skip
TRUTH = st.lists(st.tuples(st.sampled_from(CLASSES), st.sampled_from("+-"),
                           st.integers(1, 30), st.integers(0, 4)), max_size=12)  # fmt: skip


def _pred(rows: list[tuple[str, str, int]]) -> pd.DataFrame:
    return G.g3_predictions(tsv(*({"class": c, "strand": s, "end": 100 * k + (2 if s == "+"
                                                                            else 0)}
                                  for c, s, k in rows)))  # fmt: skip


@settings(max_examples=80, deadline=None)
@given(pr=PRED, tr=TRUTH, perm=st.randoms(use_true_random=False))
def test_match_counts_are_bounded_and_order_free(pr: Any, tr: Any, perm: Any) -> None:
    t = truth(*((c, s, 100 * k, 1, n, False) for c, s, k, n in tr))
    got = G.match_novel(_pred(pr), t)
    shuffled = list(pr)
    perm.shuffle(shuffled)
    again = G.match_novel(_pred(shuffled), t)
    for c in CLASSES:
        g = got[c]
        assert set(g) == set(G.NOVEL_COUNT_KEYS)
        assert g["tp_p"] + g["mismatch"] <= g["n_pred"] and g["tp_r"] <= g["n_truth"]
        assert g["n_pred"] == sum(1 for x in pr if x[0] == c)
        assert {k: g[k] for k in ("n_pred", "tp_p", "n_truth", "tp_r", "mismatch")} == {
            k: again[c][k] for k in ("n_pred", "tp_p", "n_truth", "tp_r", "mismatch")
        }


# --------------------------------------------------------------------------- G4 (rule 6)


def _g4_world() -> tuple[pd.DataFrame, pd.DataFrame]:
    t = tsv(
        *(
            {"class": c, "end": 1003 + 1000 * k}
            for k, c in enumerate(
                [
                    "intergenic",
                    "intergenic",
                    "intergenic",
                    "intergenic",
                    "cds_outframe",
                    "intergenic",
                ]
            )
        )
    )
    s = side({"q": 0.004},
             {"q": 0.02},
             {"q": math.nan, "reason": "lost_pick", "called": False},
             {"q": math.nan, "reason": "no_makedb_pair", "called": True},
             {"class": "cds_outframe", "q": 0.004, "class_rule_met": False, "called": False},
             {"q": 0.003, "regime": "denovo_subgroup"})  # fmt: skip
    return t, s


def test_g4_predictions_class_q_and_rule() -> None:
    t, s = _g4_world()
    p = G.g4_predictions(t, s, 0.01, regime="picked_decoy")
    assert tuple(p.columns) == G.NOVEL_PRED_COLUMNS
    assert p["orf_id"].tolist() == ["NOVEL_c1_1_p"] and p["stop"].tolist() == [1001]
    assert G.g4_predictions(t, s, 0.05, regime="picked_decoy")["orf_id"].tolist() == [
        "NOVEL_c1_1_p",
        "NOVEL_c1_2_p",
    ]
    both = G.g4_predictions(t, s, 0.01)["orf_id"].tolist()
    assert both == ["NOVEL_c1_1_p", "NOVEL_c1_6_p"]


def test_g4_predictions_without_q_filter_is_attack_a() -> None:
    t, s = _g4_world()
    got = G.g4_predictions(t, s, None, regime="picked_decoy")["orf_id"].tolist()
    assert got == ["NOVEL_c1_1_p", "NOVEL_c1_2_p", "NOVEL_c1_3_p", "NOVEL_c1_4_p"]


def test_g4_class_is_the_tsv_class_not_fdr_class() -> None:
    t = tsv({"class": "intergenic", "end": 3003})
    s = side({"class": "intergenic", "fdr_class": "cds_outframe", "class_mismatch": True})
    p = G.g4_predictions(t, s, 0.01)
    assert p["class"].tolist() == ["intergenic"]
    assert G.list_targets(s, "picked_decoy")["cds_outframe"] == 1


def test_g4_side_row_missing_from_tsv_fails_closed() -> None:
    with pytest.raises(ValueError):
        G.g4_predictions(tsv(), side({"q": 0.001}), 0.01)


def test_example9_fdp_counts_from_match() -> None:
    """Example 9: cds_outframe, 120 called, 4 not in truth -> n_false = n_pred - tp_p = 4."""
    t = tsv(*({"class": "cds_outframe", "end": 1003 + 10 * k} for k in range(120)))
    tr = truth(*(("cds_outframe", "+", 1001 + 10 * k, 1, 0, False) for k in range(116)))
    got = G.match_novel(
        G.g4_predictions(t, side(*({"class": "cds_outframe"} for _ in range(120))), 0.05), tr
    )
    c = got["cds_outframe"]
    assert c["n_pred"] - c["tp_p"] == 4
    assert (c["n_pred"] - c["tp_p"]) / max(c["n_pred"], 1) == pytest.approx(0.0333, abs=5e-5)


def test_list_targets_counts_winners_and_pick_losers() -> None:
    _, s = _g4_world()
    got = G.list_targets(s, "picked_decoy")
    assert got == {"intergenic": 3, "cds_outframe": 1, "antisense": 0, "nterm_extension": 0}
    assert G.list_targets(s, "denovo_subgroup")["intergenic"] == 1


def test_reason_and_fallback_counts() -> None:
    _, s = _g4_world()
    assert G.reason_counts(s)["intergenic"] == {"lost_pick": 1, "no_makedb_pair": 1}
    assert G.reason_counts(s)["antisense"] == {}
    assert G.fallback_calls(s) == {"intergenic": 1, "cds_outframe": 0, "antisense": 0,
                                   "nterm_extension": 0}  # fmt: skip
