"""Chance shared peptides in G4's ``UNESTIMATED_BUGS`` check (M7-13b lead decision 11) and the
ADR-0050 strict xfails of G4-k3 and full-tier attack GC (lead decisions 10, 13). Unmarked: runs
on every commit.

Decision 11: ``no_makedb_pair`` is right when every peptide row of a GIAP ORF names another
planted ORF's makedb entry (a short Met-initial 7-mer recurring at another ATG; statistician
``local/stat/m713b-full/nopair.py``: replicate 35 ``MVTAAGR`` at ``NOVEL_chr1_463046_p``, whose
three rows name ``GIAPORF_cds_outframe_chr1_872115_872252_p``; replicate 36 ``MTPLWTK``). Such
ORFs are counted as ``shared_peptide_no_pair`` and exempted; any other ``no_makedb_pair`` still
fails. Interface (``helpers.m7_gate_novel_records``, test plumbing):

* ``parse_makedb_id(pid) -> (cls, contig, a, b, strand) | None``: ``GIAPORF_<cls>_<contig>_<a>_
  <b>_<p|m>`` (1-based closed plus-strand span without the stop codon; ``strand`` ``+``/``-``).
* ``planted_in_entry(entry, truth) -> set[(contig, strand, stop)]``: truth ORFs whose start
  codon (``start`` on ``+``, ``end - 2`` on ``-``) lies in frame inside the entry on its strand
  (a six-frame entry is the stop-to-stop ORF; an ``nterm_extension`` entry holds the upstream
  part, where the planted extension starts).
* ``shared_peptide_no_pair(peptides, rows, own_key, truth) -> bool``: the input rows whose
  peptide equals one of ``peptides`` (I/L-insensitive) are at least one, and every row names
  at least one ``GIAPORF_`` entry and each such entry holds a planted ORF other than
  ``own_key``. No rows, a gene locus, a decoy (``rev_``) or a non-planted entry: False.
* ``shared_no_pair_counts(novel_orfs, side, rows, truth, lengths) -> {class: int}`` over the
  four classes: side rows with ``reason == "no_makedb_pair"`` that pass the test above.
* ``unestimated_bugs(reasons, shared, bug_reasons) -> {reason: int}``: the bug reasons left
  after ``no_makedb_pair`` loses the ``shared`` count (``shared`` above it: ValueError).
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers import m7_gate_novel_records as N
from hypothesis import given, settings
from hypothesis import strategies as st

HERE = Path(__file__).parent
ENTRY_35 = "GIAPORF_cds_outframe_chr1_872115_872252_p"
BUGS = frozenset({"paired_share_low", "decoys_unscored", "no_makedb_pair", "decoys_not_makedb",
                  "denovo_inputs_missing"})  # fmt: skip


def truth(*rows: tuple[str, str, int, int, str]) -> pd.DataFrame:
    """``(class, contig, start, end, strand)`` -> a read_novel_orfs-shaped frame."""
    out = []
    for i, (cls, c, s, e, st_) in enumerate(rows):
        out.append({"orf_id": f"NOVEL_{i:04d}", "class": cls, "contig": c, "start": s, "end": e,
                    "strand": st_, "stop": e - 2 if st_ == "+" else s, "true_start": s,
                    "n_true_peptides": 3, "rbs_planted": False})  # fmt: skip
    return pd.DataFrame(out)


def rows(*pairs: tuple[str, str]) -> pd.DataFrame:
    return pd.DataFrame([{"peptide": p, "protein_ids": ids, "is_decoy": "true" if
                          ids.startswith("rev_") else "false"} for p, ids in pairs])  # fmt: skip


#: replicate 35 (nopair.out): the planted cds_outframe ORF of the named entry, start in frame
T35 = truth(("cds_outframe", "chr1", 872145, 872255, "+"),
            ("cds_outframe", "chr1", 462976, 463080, "+"))  # fmt: skip
OWN_35 = ("chr1", "+", 463082)  # NOVEL_chr1_463046_p, end 463084
R35 = rows(("MVTAAGR", ENTRY_35), ("MVTAAGR", ENTRY_35), ("MVTAAGR", ENTRY_35))


# --------------------------------------------------------------------------- parse / entry


def test_parse_makedb_id() -> None:
    assert N.parse_makedb_id(ENTRY_35) == ("cds_outframe", "chr1", 872115, 872252, "+")
    assert N.parse_makedb_id("GIAPORF_nterm_extension_chr_2_5_94_m") == (
        "nterm_extension", "chr_2", 5, 94, "-")  # fmt: skip
    for bad in ("GSIM_02040", "rev_" + ENTRY_35, "GIAPORF_other_chr1_1_9_p", "GIAPORF_x"):
        assert N.parse_makedb_id(bad) is None, bad


def test_planted_in_entry_plus_in_frame_only() -> None:
    e = N.parse_makedb_id(ENTRY_35)
    assert N.planted_in_entry(e, T35) == {("chr1", "+", 872253)}
    shifted = truth(("cds_outframe", "chr1", 872146, 872256, "+"))
    assert N.planted_in_entry(e, shifted) == set()
    other_strand = truth(("cds_outframe", "chr1", 872145, 872255, "-"))
    assert N.planted_in_entry(e, other_strand) == set()


def test_planted_in_entry_minus_strand() -> None:
    e = N.parse_makedb_id("GIAPORF_antisense_c1_961_990_m")
    t = truth(("antisense", "c1", 958, 990, "-"), ("antisense", "c1", 958, 989, "-"))
    assert N.planted_in_entry(e, t) == {("c1", "-", 958)}


def test_planted_in_entry_nterm_upstream_start() -> None:
    e = N.parse_makedb_id("GIAPORF_nterm_extension_c1_101_400_p")
    t = truth(("nterm_extension", "c1", 131, 2000, "+"), ("intergenic", "c1", 500, 600, "+"))
    assert N.planted_in_entry(e, t) == {("c1", "+", 1998)}


# --------------------------------------------------------------------------- the exemption


def test_replicate35_shape_is_shared() -> None:
    assert N.shared_peptide_no_pair(["MVTAAGR"], R35, OWN_35, T35)


def test_replicate36_shape_is_shared_across_classes() -> None:
    t = truth(("cds_outframe", "chr1", 309570, 309758, "+"))
    r = rows(*[("MTPLWTK", "GIAPORF_cds_outframe_chr1_309570_309755_p")] * 3)
    assert N.shared_peptide_no_pair(["MTPLWTK"], r, ("chr1", "+", 546230), t)


def test_il_insensitive_match() -> None:
    assert N.shared_peptide_no_pair(["MVTAAGR"], rows(("MVTAAGR", ENTRY_35)), OWN_35, T35)
    t = truth(("cds_outframe", "chr1", 872145, 872255, "+"))
    assert N.shared_peptide_no_pair(["MLNRLVR"], rows(("MINRIVR", ENTRY_35)), OWN_35, t)


def test_no_rows_is_not_shared() -> None:
    """Replicate 15's shape (nopair.out: no input row for the ORF's peptide) stays a bug."""
    assert not N.shared_peptide_no_pair(["MLNRLVR"], R35, ("chr1", "+", 282380), T35)
    assert not N.shared_peptide_no_pair([], R35, OWN_35, T35)


def test_non_planted_entry_is_not_shared() -> None:
    r = rows(("MVTAAGR", "GIAPORF_antisense_chr1_10_99_p"))
    assert not N.shared_peptide_no_pair(["MVTAAGR"], r, OWN_35, T35)


def test_one_row_off_planted_entry_breaks_the_exemption() -> None:
    for extra in ("GSIM_02040", "GIAPORF_antisense_chr1_10_99_p", "rev_" + ENTRY_35):
        r = pd.concat([R35, rows(("MVTAAGR", extra))], ignore_index=True)
        assert not N.shared_peptide_no_pair(["MVTAAGR"], r, OWN_35, T35), extra


def test_own_planted_orf_is_not_another() -> None:
    assert not N.shared_peptide_no_pair(["MVTAAGR"], R35, ("chr1", "+", 872253), T35)


def test_multi_id_row_needs_every_entry_planted() -> None:
    ok = rows(("MVTAAGR", f"{ENTRY_35};GIAPORF_cds_outframe_chr1_462976_463077_p"))
    assert N.shared_peptide_no_pair(["MVTAAGR"], ok, OWN_35, T35)
    bad = rows(("MVTAAGR", f"{ENTRY_35};GIAPORF_antisense_chr1_10_99_p"))
    assert not N.shared_peptide_no_pair(["MVTAAGR"], bad, OWN_35, T35)


def _novel(orf_id: str, cls: str, end: int, peptides: str) -> dict[str, str]:
    return {"orf_id": orf_id, "class": cls, "contig": "chr1", "strand": "+", "end": str(end),
            "peptides": peptides}  # fmt: skip


def test_shared_no_pair_counts_per_class() -> None:
    novel = pd.DataFrame([
        _novel("A", "cds_outframe", 463084, "MVTAAGR"),
        _novel("B", "antisense", 546232, "MTPLWTK"),
        _novel("C", "intergenic", 282382, "MVTAAGR"),
    ])  # fmt: skip
    side = pd.DataFrame({"orf_id": ["A", "B", "C"], "class": ["cds_outframe", "antisense",
                         "intergenic"], "reason": ["no_makedb_pair", "no_makedb_pair",
                         "lost_pick"]})  # fmt: skip
    r = pd.concat([R35, rows(("MTPLWTK", "GIAPORF_antisense_chr1_10_99_p"))], ignore_index=True)
    got = N.shared_no_pair_counts(novel, side, r, T35, {"chr1": 1_000_000})
    assert got == {"intergenic": 0, "cds_outframe": 1, "antisense": 0, "nterm_extension": 0}


# --------------------------------------------------------------------------- bug accounting


def test_unestimated_bugs_branches() -> None:
    reasons = {"no_makedb_pair": 3, "lost_pick": 5}
    assert N.unestimated_bugs(reasons, 3, BUGS) == {}
    assert N.unestimated_bugs(reasons, 2, BUGS) == {"no_makedb_pair": 1}
    assert N.unestimated_bugs({"decoys_unscored": 1, "no_makedb_pair": 1}, 1, BUGS) == {
        "decoys_unscored": 1}  # fmt: skip
    assert N.unestimated_bugs({"lost_pick": 2}, 0, BUGS) == {}
    with pytest.raises(ValueError):
        N.unestimated_bugs({"no_makedb_pair": 1}, 2, BUGS)


REASON = st.sampled_from(sorted(BUGS | {"lost_pick", "below_rule"}))


@settings(max_examples=100, deadline=None)
@given(reasons=st.dictionaries(REASON, st.integers(1, 50)), frac=st.floats(0, 1))
def test_unestimated_bugs_drops_exactly_the_shared(reasons: dict[str, int], frac: float) -> None:
    n = reasons.get("no_makedb_pair", 0)
    shared = int(frac * n)
    got = N.unestimated_bugs(reasons, shared, BUGS)
    assert sum(got.values()) == sum(v for k, v in reasons.items() if k in BUGS) - shared
    assert set(got) <= BUGS and all(v > 0 for v in got.values())
    assert got.get("no_makedb_pair", 0) == n - shared


@settings(max_examples=60, deadline=None)
@given(il=st.lists(st.booleans(), min_size=7, max_size=7), extra=st.sampled_from(
    ["GSIM_1", "GIAPORF_intergenic_chr1_10_99_p", "rev_" + ENTRY_35]))  # fmt: skip
def test_il_invariant_and_one_bad_row_vetoes(il: list[bool], extra: str) -> None:
    pep = "".join("I" if b else "L" for b in il)
    r = rows((pep, ENTRY_35))
    assert N.shared_peptide_no_pair([pep.replace("I", "L")], r, OWN_35, T35)
    r2 = pd.concat([r, rows((pep, extra))], ignore_index=True)
    assert not N.shared_peptide_no_pair([pep], r2, OWN_35, T35)


# --------------------------------------------------------------------------- ADR-0050 xfails


def _load(path: Path, name: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _xfails(param_sets: Any) -> dict[tuple[Any, ...], Any]:
    out = {}
    for p in param_sets:
        marks = [m for m in getattr(p, "marks", ()) if m.name == "xfail"]
        out[tuple(p.values) if hasattr(p, "values") else tuple(p)] = marks[0] if marks else None
    return out


def _check(mark: Any) -> None:
    assert mark.kwargs.get("strict") is True and mark.kwargs.get("raises") is AssertionError
    assert "ADR-0050" in mark.kwargs.get("reason", "")


def test_g4k3_xfail_cells_are_exactly_decision_10() -> None:
    T = _load(HERE / "test_novel_orfs_gate.py", "_novel_gate_for_marks")
    got = _xfails(T.K3_CASES)
    want = {("cds_outframe", 0.01), ("cds_outframe", 0.05), ("antisense", 0.01),
            ("antisense", 0.05), ("intergenic", 0.01)}  # fmt: skip
    assert {k for k, m in got.items() if m is not None} == want
    assert len(got) == 9
    for k in want:
        _check(got[k])


def test_gc_xfail_only_in_the_full_tier_at_001() -> None:
    """Lead decision 13: fast GC is a hard check at every level; the full-tier test (q 0.01
    only) carries the strict ADR-0050 xfail."""
    A = _load(HERE / "attacks" / "test_attack_m7_gate.py", "_attack_m7_for_marks")
    fast = A.test_attack_gc_poor_intergenic_fdp_at_level
    assert not [m for m in getattr(fast, "pytestmark", []) if m.name == "xfail"]
    par = [m for m in fast.pytestmark if m.name == "parametrize"]
    assert len(par) == 1
    assert all(not getattr(p, "marks", ()) for p in par[0].args[1])
    full = A.test_attack_gc_poor_intergenic_fdp_at_level_full
    names = [m.name for m in full.pytestmark]
    assert "calibration_full" in names and "parametrize" not in names
    xf = [m for m in full.pytestmark if m.name == "xfail"]
    assert len(xf) == 1
    _check(xf[0])
