"""Matching helpers of the M7-13a variants gate, ``tests/calibration/_m7_gate.py`` (implementer;
docs/specs/M7-13-calibration-gate.md rules 3, 4, 9, examples 1-6 and 10, the variant rows of
"Degenerate cases"; lead decisions 1-2). Unmarked: runs on every commit.

Interface pinned here (pure functions, pandas in, ints out):

* ``GATE_R`` = ``{"fast": 10, "full": 300}`` (lead decision 1, ADR-0048 pending);
  ``INDEL_MASK_NT`` = 3; ``ISOBARIC_ERRORS`` = ``{"block", "modification"}`` (giapsim's
  mass-equivalent error types; ``near_block`` and ``swap`` are not).
* ``fold_il(s)``: I, L and J become one letter, every other residue unchanged.
* ``g1_predictions(variants)``: ``variants.tsv`` as written (strings, ``codon_start`` 1-based) ->
  columns ``contig, codon_start`` (int), ``strand, alt_aa, mass_equivalent`` (bool): the rows with
  ``confidence = high`` and ``alt_type = substitution`` (a ``mixed`` site's ALT rows included),
  one per ``(contig, codon_start, strand, fold_il(alt_aa))``.
* ``match_variants(pred, alt_sites, indel_sites)`` -> ``{n_pred, tp_p, n_truth, tp_r, masked}``
  (rules 3-4). Truth frames as ``giapsim.truth.read_alt_sites`` / ``read_indel_sites`` give them.
  Masked: the codon ``[c, c+2]`` and a ``junction`` j on the same contig with
  ``c - 3 <= j <= c + 2 + 3``. A masked prediction is removed before matching (neither TP nor
  FP, so it cannot recall a truth row either) and counted in ``masked``.
* ``g2_violations(variants)``: ``high`` rows whose flags hold neither ``start_codon`` nor
  ``stop_codon``; ``start_codon_high(variants)``: ``high`` rows flagged ``start_codon``.
* ``a6a_unflagged(variants, site_peptides, denovo_truth)`` -> ``{n_rows, n_unflagged}``: rows
  (any confidence) with ``alt_type`` in {substitution, block} whose supporting de novo reads
  (``truth/denovo.tsv`` rows whose ``peptide`` equals one of the row's peptides after
  ``fold_il``) are non-empty and all have ``error_type`` in ``ISOBARIC_ERRORS``; ``n_unflagged``
  of them lack ``mass_equivalent``. ``site_peptides`` maps ``(contig, codon_start (int,
  1-based), strand, alt_type, alt_aa)`` to the row's ALT peptides.
* ``a6_high_swaps(variants, min_peptides=2)`` -> ``{folded, unfolded}``: swap rows grouped by
  ``(contig, codon_start, strand, fold_il(alt_aa))``; ``folded`` counts the groups with a
  ``high`` row or ``sum(n_alt) >= min_peptides`` (rule 9), ``unfolded`` the ``high`` rows.
"""

from __future__ import annotations

from typing import Any

import _criteria as C
import _m7_gate as G
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.variants.vcf_records import VARIANTS_COLUMNS

ALT_COLS = ["contig", "codon_start", "strand", "gene", "ref_aa", "alt_aa", "n_true_peptides",
            "frameshift_upstream"]  # fmt: skip


def vt(*rows: dict[str, Any]) -> pd.DataFrame:
    """A written ``variants.tsv`` (all strings) from partial rows."""
    base = {c: "NA" for c in VARIANTS_COLUMNS}
    base.update(contig="chr1", strand="+", ref_aa="K", alt_type="substitution", call="ALT",
                confidence="high", n_ref="0", n_alt="2", flags="")  # fmt: skip
    out = [{**base, **{k: str(v) for k, v in r.items()}} for r in rows]
    return pd.DataFrame(out, columns=list(VARIANTS_COLUMNS), dtype=str)


def alt_sites(*rows: tuple[int, str, str, str, int]) -> pd.DataFrame:
    """Truth rows ``(codon_start, strand, ref_aa, alt_aa, n_true_peptides)`` on chr1."""
    data = [["chr1", c, s, "G1", r, a, n, False] for c, s, r, a, n in rows]
    return pd.DataFrame(data, columns=ALT_COLS)


def indels(*junctions: int, contig: str = "chr1") -> pd.DataFrame:
    return pd.DataFrame({"contig": [contig] * len(junctions), "junction": list(junctions)})


def pred(*rows: tuple[int, str, str]) -> pd.DataFrame:
    """Predictions ``(codon_start, strand, alt_aa)`` on chr1 through ``g1_predictions``."""
    return G.g1_predictions(vt(*({"codon_start": c, "strand": s, "alt_aa": a} for c, s, a in rows)))


NONE = indels()


def test_constants() -> None:
    assert dict(G.GATE_R) == {"fast": 10, "full": 300}
    assert G.INDEL_MASK_NT == 3
    assert set(G.ISOBARIC_ERRORS) == {"block", "modification"}


def test_fold_il() -> None:
    assert G.fold_il("I") == G.fold_il("L") == G.fold_il("J")
    assert G.fold_il("IKLJ") == G.fold_il("LKIL") and G.fold_il("K") == "K"


def test_example1_counts_and_f1() -> None:
    recall = [(100 * k, "+", "K", "E", 2) for k in range(1, 10)]
    truth = alt_sites(*recall, (1000, "+", "K", "E", 1))
    p = pred(*((100 * k, "+", "E") for k in range(1, 8)), (1000, "+", "E"), (5000, "+", "E"),
             (6000, "+", "E"))  # fmt: skip
    m = G.match_variants(p, truth, NONE)
    assert m == {"n_pred": 10, "tp_p": 8, "n_truth": 9, "tp_r": 7, "masked": 0}
    assert C.pooled_f1([m])[2] == pytest.approx(0.789, abs=5e-4)


def test_example3_il_folded_on_minus_strand() -> None:
    truth = alt_sites((4201, "-", "K", "I", 2))
    for alt in ("L", "J", "I"):
        m = G.match_variants(pred((4201, "-", alt)), truth, NONE)
        assert (m["tp_p"], m["tp_r"]) == (1, 1), alt
    assert G.match_variants(pred((4201, "+", "L")), truth, NONE)["tp_p"] == 0  # strand differs


def test_example4_il_truth_leaves_recall_only() -> None:
    truth = alt_sites((300, "+", "L", "I", 5), (600, "+", "K", "E", 2))
    m = G.match_variants(pred((300, "+", "I")), truth, NONE)
    assert m == {"n_pred": 1, "tp_p": 1, "n_truth": 1, "tp_r": 0, "masked": 0}


@pytest.mark.parametrize(("junction", "masked"), [(104, 1), (105, 1), (106, 0), (97, 1), (96, 0),
                                                  (101, 1)])  # fmt: skip
def test_example5_indel_mask_boundaries(junction: int, masked: int) -> None:
    truth = alt_sites((100, "+", "N", "K", 2))
    m = G.match_variants(pred((100, "+", "K")), truth, indels(junction))
    assert m["masked"] == masked
    assert (m["n_pred"], m["tp_p"], m["tp_r"]) == ((0, 0, 0) if masked else (1, 1, 1))
    assert G.match_variants(pred((100, "+", "K")), truth, indels(104, contig="p1"))["masked"] == 0


def test_g1_predictions_filters_and_dedups() -> None:
    v = vt({"codon_start": 10, "alt_aa": "E", "flags": "mass_equivalent;snv"},
           {"codon_start": 10, "alt_aa": "D"},  # multi-allelic: one prediction per residue
           {"codon_start": 20, "alt_aa": "L"}, {"codon_start": 20, "alt_aa": "J"},  # folded
           {"codon_start": 30, "alt_aa": "E", "call": "mixed"},  # mixed: its high ALT counts
           {"codon_start": 40, "alt_aa": "E", "confidence": "candidate"},
           {"codon_start": 50, "alt_aa": "GG", "alt_type": "block"},
           {"codon_start": 60, "alt_aa": "KE", "alt_type": "swap"})  # fmt: skip
    p = G.g1_predictions(v)
    assert sorted(zip(p["codon_start"].astype(int), p["alt_aa"].map(G.fold_il), strict=True)) == [
        (10, "D"), (10, "E"), (20, G.fold_il("L")), (30, "E")]  # fmt: skip
    me = dict(zip(zip(p["codon_start"].astype(int), p["alt_aa"], strict=True),
                  p["mass_equivalent"], strict=True))  # fmt: skip
    assert me[(10, "E")] and not me[(10, "D")]


def test_no_prediction_and_no_truth() -> None:
    truth = alt_sites((100, "+", "K", "E", 2), (200, "+", "K", "E", 3))
    assert G.match_variants(pred(), truth, NONE) == {
        "n_pred": 0, "tp_p": 0, "n_truth": 2, "tp_r": 0, "masked": 0}  # fmt: skip
    assert G.match_variants(pred((100, "+", "E")), alt_sites(), NONE) == {
        "n_pred": 1, "tp_p": 0, "n_truth": 0, "tp_r": 0, "masked": 0}  # fmt: skip


def test_example6_g2() -> None:
    ok = vt({"codon_start": 1, "ref_aa": "V", "alt_aa": "M", "flags": "start_codon"},
            {"codon_start": 9, "alt_aa": "E", "confidence": "candidate"})  # fmt: skip
    assert G.g2_violations(ok) == 0 and G.start_codon_high(ok) == 1
    bad = vt({"codon_start": 300, "ref_aa": "N", "alt_aa": "D", "flags": "mass_equivalent"})
    assert G.g2_violations(bad) == 1 and G.start_codon_high(bad) == 0
    assert G.g2_violations(vt({"codon_start": 900, "alt_aa": "Q", "flags": "stop_codon"})) == 0


def test_a6a_rows_and_flags() -> None:
    v = vt(
        {"codon_start": 10, "ref_aa": "A", "alt_aa": "S", "flags": "mass_equivalent"},
        {
            "codon_start": 20,
            "ref_aa": "N",
            "alt_aa": "GG",
            "alt_type": "block",
            "confidence": "candidate",
        },  # unflagged, isobaric reads only
        {"codon_start": 30, "alt_aa": "E"},  # substitution error read: not isobaric
        {"codon_start": 40, "alt_aa": "KE", "alt_type": "swap"},  # swaps are A6b's
        {"codon_start": 50, "ref_aa": "N", "alt_aa": "D"},  # one read correct
        {"codon_start": 60, "ref_aa": "G", "alt_aa": "A"},
    )  # no read in the truth
    sp = {("chr1", c, "+", t, a): frozenset({p}) for c, t, a, p in [
        (10, "substitution", "S", "PEPSLDE"), (20, "block", "GG", "PGGTIDE"),
        (30, "substitution", "E", "PEPEKR"), (40, "swap", "KE", "PKEPR"),
        (50, "substitution", "D", "PDPR"), (60, "substitution", "A", "PAAPR")]}  # fmt: skip
    tr = pd.DataFrame(
        [["r", "1", "PEPSIDE", "PEPAIDE", "modification", "false"],  # I/L-folded match
         ["r", "2", "PGGTLDE", "PNTLDE", "block", "false"],
         ["r", "3", "PEPEKR", "PEPKKR", "substitution", "false"],
         ["r", "4", "PKEPR", "PEKPR", "swap", "false"],
         ["r", "5", "PDPR", "PNPR", "modification", "false"],
         ["r", "6", "PDPR", "PDPR", "none", "true"]],
        columns=["run", "spectrum_id", "peptide", "true_peptide", "error_type", "correct"],
    )  # fmt: skip
    assert G.a6a_unflagged(v, sp, tr) == {"n_rows": 2, "n_unflagged": 1}


def test_example10_a6b_swaps() -> None:
    v = vt({"codon_start": 10, "ref_aa": "IK", "alt_aa": "KL", "alt_type": "swap"},
           {"codon_start": 20, "ref_aa": "N", "alt_aa": "GG", "alt_type": "block",
            "flags": "mass_equivalent"})  # fmt: skip
    assert G.a6_high_swaps(v) == {"folded": 1, "unfolded": 1}
    split = vt(*({"codon_start": 30, "ref_aa": "LK", "alt_aa": a, "alt_type": "swap", "n_alt": 1,
                  "confidence": "candidate"} for a in ("KL", "KI")))  # fmt: skip
    assert G.a6_high_swaps(split) == {"folded": 1, "unfolded": 0}  # the I/L split cannot hide it
    assert G.a6_high_swaps(split, min_peptides=3) == {"folded": 0, "unfolded": 0}
    share = [1.0] + [0.0] * 19
    res = C.check_fwer(share, 0.05)
    assert res.ok and res.bound == pytest.approx(0.05 + 2 * (0.05 * 0.95 / 20) ** 0.5)


# --------------------------------------------------------------------------- properties

CODONS = st.lists(st.tuples(st.integers(1, 60).map(lambda k: 10 * k), st.sampled_from("+-"),
                            st.sampled_from("KEDLIN"), st.sampled_from("KEDLIN"),
                            st.integers(0, 4)), max_size=15)  # fmt: skip


@settings(max_examples=60, deadline=None)
@given(rows=CODONS, flip=st.booleans())
def test_truth_as_prediction_is_perfect_and_il_invariant(rows: list, flip: bool) -> None:
    rows = sorted({r for r in rows if r[2] != r[3]})
    truth = alt_sites(*rows)
    swap = {"I": "L", "L": "I"}
    p = pred(*{(c, s, swap.get(a, a) if flip else a) for c, s, _r, a, _n in rows})
    m = G.match_variants(p, truth, NONE)
    assert m["tp_p"] == m["n_pred"] and m["tp_r"] == m["n_truth"]
    assert m["n_truth"] == sum(n >= 2 and not {r, a} <= {"I", "L"} for _c, _s, r, a, n in
                               {(c, s, r, a, n) for c, s, r, a, n in rows})  # fmt: skip


@settings(max_examples=60, deadline=None)
@given(rows=CODONS, preds=st.lists(st.tuples(st.integers(1, 60).map(lambda k: 10 * k),
                                             st.sampled_from("+-"), st.sampled_from("KEDLIN")),
                                   max_size=15), junction=st.integers(0, 700))  # fmt: skip
def test_counts_are_consistent(rows: list, preds: list, junction: int) -> None:
    m = G.match_variants(pred(*preds), alt_sites(*rows), indels(junction))
    assert 0 <= m["tp_p"] <= m["n_pred"] and 0 <= m["tp_r"] <= m["n_truth"]
    assert m["n_pred"] + m["masked"] == len(pred(*preds))
