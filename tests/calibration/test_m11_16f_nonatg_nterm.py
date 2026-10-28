"""Lead decision 11 amendment (M11-16f; statistician-approved oracle change on #481): a
``no_makedb_pair`` ORF whose every peptide row names a non-ATG annotated gene by locus tag, the
peptide being that gene's initiator-M N-terminal tryptic peptide, is a chance recurrence, not a
GIPPA bug. Unmarked: runs on every commit.

Diagnosis (code 25, ``picked_decoy``; replicates 14, 83 (two ORFs), 179, 213; seeds 581532778,
1832983731, 1835538576, 24629369): gippasim's peptide of a GTG-start gene keeps the initiator M
(``MLSVFLR``, ``protein_ids`` = the locus tag ``GSIM_01145``), six-frame translation reads the
gene's GTG as V (``VLSVFLR``), so the row misses its gene and lands on a chance ATG copy
(``NOVEL_chr1_407412_p``, ``start`` 407412, ``end`` 407507; ``novel_stop`` key ``("chr1", "+",
407505)``, the 1-based lowest stop base, i.e. 407504 0-based):
a novel ORF with no makedb entry. Lead decision 11's ``shared_peptide_no_pair`` is vetoed by the
locus tag, so 2 cds_outframe and 3 antisense ORFs were left as ``UNESTIMATED_BUGS``.

Interface (``helpers.m7_gate_novel_records``, test plumbing):

* ``gene_table(bakta_json) -> {locus: feature}``: the Bakta-like JSON's features by ``locus``.
* ``nterm_tryptic(peptide, aa) -> bool``: ``peptide`` is an N-terminal Trypsin/P product of
  ``aa`` (I/L-folded prefix ending after K/R, or the whole protein).
* ``nonatg_nterm_recurrence(peptides, rows, orf, genes) -> bool`` with ``orf = (contig, lo,
  hi)`` (1-based closed): the input rows carrying one of ``peptides``
  (I/L-folded) are at least one, and every name in every such row is a locus tag of ``genes``
  whose ``start_type`` is a codon other than ATG, whose ``aa`` begins with M, whose N-terminal
  tryptic peptide is the row's peptide, and which does not overlap the ORF (same contig, either
  strand; an origin-crossing gene or ORF counts as overlapping: fails closed). Anything else
  (no rows, a ``GIPPAORF_`` entry, a decoy, an unknown name, a missing field): False.
* ``orf_span(contig, strand, start_longest, end) -> (contig, lo, hi)``: the ORF's widest
  1-based closed span (on ``-`` the start is the upper base; origin-crossing: ``lo > hi``).
* ``nonatg_nterm_counts(novel_orfs, side, rows, genes) -> {class: int}`` over the four classes:
  side rows with ``reason == "no_makedb_pair"`` that pass the check above on ``orf_span``.
* ``_diag(run, rows, truth, lengths, genes)`` records it as ``NONATG_NTERM`` =
  ``"nonatg_nterm_recurrence"``, separately from ``shared_peptide_no_pair``; the G4 (and G4-k3)
  checks of arms 11, 4 and 25 subtract both from ``no_makedb_pair`` and log both.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import _m7_gate as G
import pandas as pd
from helpers import m7_gate_novel_records as M
from hypothesis import given, settings
from hypothesis import strategies as st

HERE = Path(__file__).parent
TAG = "GSIM_01145"
AA = "MLSVFLRTYSMGQVPELIGRRGLEYAGIEFATNNEWTGTWTYILLSEPTMDKVGQNLGAVNYIESPSRKPGHFHLTEIME"
PEP = "MLSVFLR"
ORF = ("chr1", 407412, 407507)  # NOVEL_chr1_407412_p (replicate 14, code 25)
OWN = ("chr1", "+", 407505)
ENTRY = "GIPPAORF_cds_outframe_chr1_872115_872252_p"
EMPTY_TRUTH = pd.DataFrame(columns=["orf_id", "class", "contig", "start", "end", "strand", "stop"])


def gene(**kw: Any) -> dict[str, Any]:
    """GSIM_01145 of replicate 14 (bakta_like.json), fields overridable."""
    g = {"locus": TAG, "type": "cds", "sequence": "chr1", "start": 248839, "stop": 249717,
         "strand": "+", "start_type": "GTG", "aa": AA}  # fmt: skip
    g.update(kw)
    return g


def genes(**kw: Any) -> dict[str, dict[str, Any]]:
    return {TAG: gene(**kw)}


def rows(*pairs: tuple[str, str]) -> pd.DataFrame:
    return pd.DataFrame([{"peptide": p, "protein_ids": ids, "is_decoy": "true" if
                          ids.startswith("rev_") else "false"} for p, ids in pairs])  # fmt: skip


R14 = rows((PEP, TAG), (PEP, TAG), (PEP, TAG))


# --------------------------------------------------------------------------- the repro


def test_repro_gtg_start_reads_v_at_code_25() -> None:
    from gippa.mapping.sixframe import translate_text

    assert translate_text("GTGTTATCTGTGTTCCTTCGC", 25) == "VLSVFLR"
    assert gene()["aa"].startswith(PEP)  # the annotation keeps the initiator M


def test_repro_orf_stop_key() -> None:
    """The own key the counts build (``novel_stop``: 1-based; 407504 is the 0-based base)."""
    contig, _, end = ORF
    assert (contig, "+", G.novel_stop("+", end, 1_000_000)) == OWN


def test_repro_is_a_nonatg_nterm_recurrence() -> None:
    assert M.nonatg_nterm_recurrence([PEP], R14, ORF, genes())


def test_repro_is_not_a_shared_peptide() -> None:
    """Lead decision 11 vetoes the locus tag: the amendment is a separate record."""
    assert not M.shared_peptide_no_pair([PEP], R14, OWN, EMPTY_TRUTH)
    assert M.NONATG_NTERM == "nonatg_nterm_recurrence" != M.SHARED_NO_PAIR


# --------------------------------------------------------------------------- each condition


def test_atg_start_gene_vetoes() -> None:
    assert not M.nonatg_nterm_recurrence([PEP], R14, ORF, genes(start_type="ATG"))
    assert not M.nonatg_nterm_recurrence([PEP], R14, ORF, genes(start_type="atg"))


def test_missing_start_type_vetoes() -> None:
    for bad in (None, "", "NA"):
        assert not M.nonatg_nterm_recurrence([PEP], R14, ORF, genes(start_type=bad)), bad
    g = genes()
    del g[TAG]["start_type"]
    assert not M.nonatg_nterm_recurrence([PEP], R14, ORF, g)


def test_ttg_start_also_counts() -> None:
    assert M.nonatg_nterm_recurrence([PEP], R14, ORF, genes(start_type="TTG"))


def test_aa_not_starting_with_peptide_vetoes() -> None:
    assert not M.nonatg_nterm_recurrence([PEP], R14, ORF, genes(aa="MKLSVFLRTYS"))
    assert not M.nonatg_nterm_recurrence([PEP], R14, ORF, genes(aa="TMLSVFLRTYS"))
    assert not M.nonatg_nterm_recurrence([PEP], R14, ORF, genes(aa=""))


def test_aa_not_beginning_with_m_vetoes() -> None:
    """A V-initial peptide at a V-initial aa is an exact match, not an initiator-M recurrence."""
    r = rows(("VLSVFLR", TAG))
    assert not M.nonatg_nterm_recurrence(["VLSVFLR"], r, ORF, genes(aa="VLSVFLRTYS"))


def test_peptide_must_end_at_a_tryptic_cut() -> None:
    assert not M.nonatg_nterm_recurrence(["MLSVFL"], rows(("MLSVFL", TAG)), ORF, genes())
    assert M.nterm_tryptic("MLSVFLR", AA)
    assert M.nterm_tryptic("MLSVFLRTYSMGQVPELIGR", AA)  # one missed cleavage, still N-terminal
    assert not M.nterm_tryptic("MLSVFL", AA)
    assert not M.nterm_tryptic("LSVFLR", AA)
    assert M.nterm_tryptic("MKP", "MKP")  # the whole protein
    assert not M.nterm_tryptic("", AA)


def test_il_folded() -> None:
    assert M.nonatg_nterm_recurrence(["MISVFIR"], R14, ORF, genes())
    assert M.nonatg_nterm_recurrence([PEP], rows(("MISVFLR", TAG)), ORF, genes())


# --------------------------------------------------------------------------- overlap


def test_overlapping_gene_vetoes() -> None:
    """The named gene covering the ORF (here on the other strand, as an antisense ORF's gene
    would): not a chance copy elsewhere."""
    covering = genes(start=407100, stop=408000, strand="-")
    assert not M.nonatg_nterm_recurrence([PEP], R14, ORF, covering)
    same_strand = genes(start=407100, stop=408000, strand="+")
    assert not M.nonatg_nterm_recurrence([PEP], R14, ORF, same_strand)


def test_overlap_boundaries_are_closed() -> None:
    """1-based closed spans: one shared base is an overlap, adjacency is not."""
    _, s, e = ORF
    assert not M.nonatg_nterm_recurrence([PEP], R14, ORF, genes(start=406000, stop=s))
    assert M.nonatg_nterm_recurrence([PEP], R14, ORF, genes(start=406000, stop=s - 1))
    assert not M.nonatg_nterm_recurrence([PEP], R14, ORF, genes(start=e, stop=409000))
    assert M.nonatg_nterm_recurrence([PEP], R14, ORF, genes(start=e + 1, stop=409000))


def test_other_contig_does_not_overlap() -> None:
    assert M.nonatg_nterm_recurrence([PEP], R14, ORF, genes(sequence="chr2", start=407100,
                                                            stop=408000))  # fmt: skip


def test_origin_crossing_fails_closed() -> None:
    assert not M.nonatg_nterm_recurrence([PEP], R14, ORF, genes(start=999_000, stop=200))
    assert not M.nonatg_nterm_recurrence([PEP], R14, ("chr1", 999_000, 200), genes())


def test_missing_coordinates_fail_closed() -> None:
    g = genes()
    del g[TAG]["start"]
    assert not M.nonatg_nterm_recurrence([PEP], R14, ORF, g)


# --------------------------------------------------------------------------- row vetoes


def test_row_vetoes() -> None:
    empty = pd.DataFrame(columns=["peptide", "protein_ids", "is_decoy"])
    assert not M.nonatg_nterm_recurrence([PEP], empty, ORF, genes())
    assert not M.nonatg_nterm_recurrence([], R14, ORF, genes())
    assert not M.nonatg_nterm_recurrence(["MTPLWTK"], R14, ORF, genes())  # no row carries it
    for ids in ("GSIM_09999", f"{TAG};{ENTRY}", f"{TAG};GSIM_09999", f"rev_{TAG}", "", ";"):
        r = pd.concat([R14, rows((PEP, ids))], ignore_index=True)
        assert not M.nonatg_nterm_recurrence([PEP], r, ORF, genes()), ids


def test_two_qualifying_genes_in_one_row() -> None:
    g = {**genes(), "GSIM_02000": gene(locus="GSIM_02000", start=600000, stop=600900,
                                       start_type="TTG")}  # fmt: skip
    assert M.nonatg_nterm_recurrence([PEP], rows((PEP, f"{TAG};GSIM_02000")), ORF, g)
    g["GSIM_02000"]["start_type"] = "ATG"
    assert not M.nonatg_nterm_recurrence([PEP], rows((PEP, f"{TAG},GSIM_02000")), ORF, g)


def test_every_peptide_row_must_qualify() -> None:
    """A second ORF peptide whose rows name another (ATG) gene vetoes the whole ORF."""
    g = {**genes(), "GSIM_00007": gene(locus="GSIM_00007", start=10, stop=900, start_type="ATG",
                                       aa="MTPLWTKAAA")}  # fmt: skip
    r = pd.concat([R14, rows(("MTPLWTK", "GSIM_00007"))], ignore_index=True)
    assert not M.nonatg_nterm_recurrence([PEP, "MTPLWTK"], r, ORF, g)
    assert M.nonatg_nterm_recurrence([PEP], r, ORF, g)


# --------------------------------------------------------------------------- disjoint records

NAMES = st.sampled_from([TAG, ENTRY, f"rev_{TAG}", "GSIM_09999", f"rev_{ENTRY}"])


@settings(max_examples=150, deadline=None)
@given(st.lists(st.lists(NAMES, min_size=1, max_size=3), min_size=1, max_size=4))
def test_exemptions_are_disjoint(name_lists: list[list[str]]) -> None:
    """An ORF is never in both records, so their sum never exceeds ``no_makedb_pair``; and the
    amendment holds only when every name is the qualifying locus tag."""
    truth = pd.DataFrame([{"orf_id": "N0", "class": "cds_outframe", "contig": "chr1",
                           "start": 872145, "end": 872255, "strand": "+",
                           "stop": 872253}])  # fmt: skip
    r = rows(*((PEP, ";".join(n)) for n in name_lists))
    nonatg = M.nonatg_nterm_recurrence([PEP], r, ORF, genes())
    shared = M.shared_peptide_no_pair([PEP], r, OWN, truth)
    assert not (nonatg and shared)
    assert nonatg == all(x == TAG for n in name_lists for x in n)


# --------------------------------------------------------------------------- counts and wiring


def _novel(oid: str, cls: str, start: int, end: int, peptides: str,
           strand: str = "+") -> dict[str, Any]:  # fmt: skip
    return {"orf_id": oid, "class": cls, "contig": "chr1", "strand": strand,
            "start_longest": str(start), "start_best_ostir": str(start), "end": str(end),
            "peptides": peptides}  # fmt: skip


def test_orf_span() -> None:
    assert M.orf_span("chr1", "+", "407412", "407507") == ORF
    # NOVEL_chr1_754806_m (replicate 83): start_longest 755033, end 754806
    assert M.orf_span("chr1", "-", "755033", "754806") == ("chr1", 754806, 755033)
    assert M.orf_span("chr1", "+", 999_000, 200) == ("chr1", 999_000, 200)  # origin-crossing


def test_minus_strand_counts_use_the_widest_span() -> None:
    """A minus-strand ORF: overlap is tested on ``end``..``start_longest``."""
    side = pd.DataFrame({"orf_id": ["M"], "class": ["antisense"], "reason": ["no_makedb_pair"]})
    clear = pd.DataFrame([_novel("M", "antisense", 755033, 754806, PEP, "-")])
    assert M.nonatg_nterm_counts(clear, side, R14, genes())["antisense"] == 1
    near = genes(start=755033, stop=756000)  # shares the ORF's upper (start) base
    assert M.nonatg_nterm_counts(clear, side, R14, near)["antisense"] == 0


def test_nonatg_nterm_counts_per_class() -> None:
    novel = pd.DataFrame([
        _novel("A", "antisense", 407412, 407507, PEP),
        _novel("B", "cds_outframe", 249000, 249095, PEP),  # inside GSIM_01145 itself
        _novel("C", "intergenic", 407412, 407507, PEP),  # lost_pick: not a no_makedb_pair
        _novel("D", "cds_outframe", 700000, 700095, "MTPLWTK"),
    ])  # fmt: skip
    side = pd.DataFrame({"orf_id": list("ABCD"),
                         "class": ["antisense", "cds_outframe", "intergenic", "cds_outframe"],
                         "reason": ["no_makedb_pair", "no_makedb_pair", "lost_pick",
                                    "no_makedb_pair"]})  # fmt: skip
    r = pd.concat([R14, rows(("MTPLWTK", ENTRY))], ignore_index=True)
    got = M.nonatg_nterm_counts(novel, side, r, genes())
    assert got == {"intergenic": 0, "cds_outframe": 0, "antisense": 1, "nterm_extension": 0}


def test_gene_table_reads_bakta_json(tmp_path: Path) -> None:
    p = tmp_path / "bakta_like.json"
    feats = [gene(), {"locus": None, "type": "crispr"}, gene(locus="GSIM_00002", start_type="ATG")]
    p.write_text(json.dumps({"features": feats, "genome": {}, "sequences": [], "version": {}}))
    t = M.gene_table(p)
    assert set(t) == {TAG, "GSIM_00002"}
    assert t[TAG]["start_type"] == "GTG" and t[TAG]["aa"] == AA
    assert M.nonatg_nterm_recurrence([PEP], R14, ORF, t)


def test_diag_records_both_keys() -> None:
    novel = pd.DataFrame([_novel("A", "antisense", 407412, 407507, PEP)])
    side = pd.DataFrame({"orf_id": ["A"], "class": ["antisense"], "fdr_class": ["antisense"],
                         "regime": ["picked_decoy"], "q": [float("nan")],
                         "class_fdr": ["unestimated"], "reason": ["no_makedb_pair"],
                         "class_mismatch": [False], "class_rule_met": [False],
                         "called": [False], "denovo_only_unsupported": [False]})  # fmt: skip
    run = M.NovelRun(Path("."), novel, side)
    d = M._diag(run, R14, EMPTY_TRUTH, {"chr1": 1_000_000}, genes())
    assert d[M.SHARED_NO_PAIR]["antisense"] == 0
    assert d[M.NONATG_NTERM] == {"intergenic": 0, "cds_outframe": 0, "antisense": 1,
                                 "nterm_extension": 0}  # fmt: skip
    exempt = d[M.SHARED_NO_PAIR]["antisense"] + d[M.NONATG_NTERM]["antisense"]
    assert M.unestimated_bugs({"no_makedb_pair": 1}, exempt, frozenset({"no_makedb_pair"})) == {}


def test_gate_checks_subtract_and_log_the_amendment() -> None:
    """Arms 11 (G4, G4-k3) and 4/25 (G4) read ``NONATG_NTERM`` from the records and log it."""
    for name in ("test_novel_orfs_gate.py", "test_m11_16f_novel_arms.py"):
        src = (HERE / name).read_text()
        assert src.count("M.NONATG_NTERM") >= (4 if name == "test_novel_orfs_gate.py" else 2), name
    helper = (HERE.parent / "helpers" / "m7_gate_novel_records.py").read_text()
    assert helper.count("gene_table(res.bakta_json)") >= 3
