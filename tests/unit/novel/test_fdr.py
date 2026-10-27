"""M7-11 class FDR for novel ORFs: examples 1-6 and every degenerate row of
docs/specs/M7-11-novel-class-fdr.md (GIPPA_SPEC.md §7.6.4, §8.1; ADR-0009, ADR-0010).

Lead decisions (2026-10-22) pinned here: per-ORF regime by peptide source, database first (1);
untagged classes ``pseudogene``/``sorf``/``cterm_extension`` are unestimated (2); de novo ORF q is
the minimum of its peptides' subgroup q (3); a class whose decoys are < 80 % paired falls back to
unestimated (4). Frame layouts and score encodings: ``tests/helpers/novel_fdr_fixture.py``.

Test-author resolutions of open points:

* A decoy row is *paired* when its stripped ID is a makedb entry (``id_map`` row) of its class;
  its class is the tag in the ID, so an ID absent from ``id_map`` is an unpaired decoy of that
  class. ADR-0009 rescaling as in M4-12: FDR(t) = (D(t)/f + 1)/T(t), f = paired share; unpaired
  decoys never enter the picked list.
* ``id_map=None``: entries are read from the class-tagged IDs themselves
  (``GIPPAORF_{class}_{contig}_{start1}_{end1}_{p|m}``).
* ``c`` (ADR-0010) is the largest π̂ over every searched de novo peptide (target or decoy rows),
  a missing π̂ counting as 1.
* ``regime(peptides, [orf])`` is the regime of that ORF (decision 1); ``regime`` on a whole run
  with both sources is not pinned.
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pandas as pd
import pytest
from helpers.novel_fdr_fixture import (
    ContigClassifier,
    PickedWorld,
    callable_at,
    hit_rows,
    make_orf,
    pep_name,
    pep_row,
    peptides_frame,
    q_of,
)

from gippa.mapping.loci import LocusClassifier
from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene, GeneModel


def _fdr():  # type: ignore[no-untyped-def]
    import gippa.novel.fdr as m  # fails (ModuleNotFoundError) until M7-11 lands

    return m


def picked_class_q(orfs, peptides, id_map):  # type: ignore[no-untyped-def]
    return _fdr().picked_class_q(orfs, peptides, id_map)


def denovo_class_q(orfs, target_hits, decoy_hits, classifier, pi):  # type: ignore[no-untyped-def]
    return _fdr().denovo_class_q(orfs, target_hits, decoy_hits, classifier, pi)


def regime(peptides, orfs):  # type: ignore[no-untyped-def]
    return _fdr().regime(peptides, orfs)


# ---------------------------------------------------------------- picked builders


def _class_block(
    w: PickedWorld,
    cls: str,
    contig: str,
    n: int,
    *,
    first_rank: int = 10,
    decoys: bool = True,
    strand: str = "+",
) -> list[str]:
    """``n`` target ORFs (ranks ``first_rank..``) each beating its own paired decoy (ranks 600+)."""
    ids = []
    for i in range(n):
        orf, did = w.target(cls, contig, first_rank + i, strand=strand)
        if decoys:
            w.decoy(did, 600 + i)
        ids.append(orf.orf_id)
    return ids


# ---------------------------------------------------------------- example 1


@pytest.mark.parametrize("prefix", ["rev_", "DECOY_"])
@pytest.mark.parametrize("with_map", [True, False])
def test_example1_150_targets_before_first_decoy_win(prefix: str, with_map: bool) -> None:
    w = PickedWorld()
    winners = []
    for i in range(150):
        orf, did = w.target("intergenic", "c1", 1 + i)
        w.decoy(did, 500 + i, prefix=prefix)
        winners.append(orf.orf_id)
    losers = []
    for i in range(5):  # the first decoy wins come after the 150
        orf, did = w.target("intergenic", "c1", 300 + i)
        w.decoy(did, 200 + i, prefix=prefix)
        losers.append(orf.orf_id)
    qs = picked_class_q(w.orfs, w.peptides, w.id_map if with_map else None)
    for o in winners:
        assert q_of(qs, o) == pytest.approx(1 / 150)
        assert round(q_of(qs, o), 4) == 0.0067
        assert callable_at(qs, o)
    for o in losers:
        assert not callable_at(qs, o, 1.0)  # lost the pick: no class q


# ---------------------------------------------------------------- example 2 / zero decoys


@pytest.mark.parametrize("decoys", [True, False], ids=["losing_decoys", "zero_decoys"])
def test_example2_50_targets_no_decoy_win(decoys: bool) -> None:
    """Example 2 and the degenerate row "a class with targets and zero decoys": q = 1/T."""
    w = PickedWorld()
    ids = _class_block(w, "intergenic", "c1", 50, decoys=decoys)
    other = _class_block(w, "cds_outframe", "c2", 3, first_rank=400)  # makedb decoys elsewhere
    qs = picked_class_q(w.orfs, w.peptides, w.id_map)
    for o in ids:
        assert q_of(qs, o) == pytest.approx(0.02)
        assert not callable_at(qs, o)
    for o in other:  # its own class list, never pooled with intergenic
        assert q_of(qs, o) == pytest.approx(1 / 3)


def test_small_class_cannot_pass_alpha() -> None:
    """B.6 logic: with the +1 rule a class needs >= 100 target wins before any ORF reaches 0.01."""
    w = PickedWorld()
    ids = _class_block(w, "antisense", "c3", 99)
    qs = picked_class_q(w.orfs, w.peptides, w.id_map)
    assert not any(callable_at(qs, o) for o in ids)
    w2 = PickedWorld()
    ids2 = _class_block(w2, "antisense", "c3", 100)
    qs2 = picked_class_q(w2.orfs, w2.peptides, w2.id_map)
    assert all(callable_at(qs2, o) for o in ids2)


# ---------------------------------------------------------------- example 3


def test_example3_decoy_outscores_target_counts_as_decoy_win() -> None:
    w = PickedWorld()
    entry = w.entry("cds_outframe", "c1", "+", start0=100, length=300)
    assert entry == ("GIPPAORF_cds_outframe_c1_101_400_p", 400)
    x, did = w.target("cds_outframe", "c1", 50, entry=entry)
    assert x.stop == 400
    assert w.decoy(did, 10)  # rev_GIPPAORF_cds_outframe_c1_101_400_p scores higher
    assert w.pep_rows[-1]["protein_ids"] == ["rev_GIPPAORF_cds_outframe_c1_101_400_p"]
    others = _class_block(w, "cds_outframe", "c1", 10, first_rank=100)
    qs = picked_class_q(w.orfs, w.peptides, w.id_map)
    assert not callable_at(qs, x.orf_id, 1.0)
    for o in others:  # one decoy win above all 10: (1 + 1) / 10
        assert q_of(qs, o) == pytest.approx(0.2)


# ---------------------------------------------------------------- pairing (rule 2)


@pytest.mark.parametrize(
    ("cls", "contig", "strand", "by_tag"),
    [("cds_outframe", "c2", "+", False), ("antisense", "c3", "-", False),
     ("intergenic", "c1", "-", False), ("nterm_extension", "c4", "+", True)],
)  # fmt: skip
def test_target_pairs_with_its_own_decoy(cls: str, contig: str, strand: str, by_tag: bool) -> None:
    """X beats its own decoy, which beats every other target of the class. Paired correctly the
    decoy leaves the list (q = 1/151 for all); unpaired it would be a decoy win (2/150).
    Plus strand: stop = entry end; minus: stop = entry start - 3; ``nterm_extension``: by locus
    tag (the entry's stop lies inside the gene, the ORF's stop is the gene's)."""
    w = PickedWorld()
    if by_tag:
        entry = w.entry(cls, contig, strand, locus_tag="g7", start0=500, length=189)
        x, did = w.target(cls, contig, 1, entry=entry, locus_tag="g7", stop=90_007)
        assert x.stop != entry[1]
    else:
        x, did = w.target(cls, contig, 1, strand=strand)
        if strand == "-":
            row = w.id_rows[-1]
            assert x.stop == row[5] - 1 - 3  # 1-based start -> 0-based, minus the stop codon
    w.decoy(did, 2)
    others = []
    for i in range(150):
        if by_tag:
            tag = f"g{100 + i}"
            e = w.entry(cls, contig, strand, locus_tag=tag)
            o, d = w.target(cls, contig, 10 + i, entry=e, locus_tag=tag, stop=95_000 + 600 * i)
        else:
            o, d = w.target(cls, contig, 10 + i, strand=strand)
        w.decoy(d, 600 + i)
        others.append(o.orf_id)
    qs = picked_class_q(w.orfs, w.peptides, w.id_map)
    for o in [x.orf_id, *others]:
        assert q_of(qs, o) == pytest.approx(1 / 151)


def test_target_score_uses_il_folded_peptides() -> None:
    """C9 keeps I/L as written; ``Orf.peptides`` is I/L-folded."""
    w = PickedWorld()
    ids = _class_block(w, "intergenic", "c1", 150)
    w.pep_rows[0]["peptide"] = "PEPTIDEKAAAA"
    w.orfs[0] = dataclasses.replace(w.orfs[0], peptides=("PEPTLDEKAAAA",))
    qs = picked_class_q(w.orfs, w.peptides, w.id_map)
    assert q_of(qs, ids[0]) == pytest.approx(1 / 150)


def test_orf_score_is_its_best_peptide() -> None:
    """Score of a target ORF = best (minimum) input q over its peptides: an ORF with peptides at
    ranks 1 and 700 beats a decoy at rank 5."""
    w = PickedWorld()
    x, did = w.target("intergenic", "c1", [700, 1])
    w.decoy(did, 5)
    others = _class_block(w, "intergenic", "c1", 150)
    qs = picked_class_q(w.orfs, w.peptides, w.id_map)
    for o in [x.orf_id, *others]:
        assert q_of(qs, o) == pytest.approx(1 / 151)


@pytest.mark.parametrize(
    ("t", "d", "decoy_wins"),
    [((1e-5, 999.0), (1e-5, 999.0), True),  # full tie goes to the decoy
     ((1e-5, 990.0), (1e-5, 999.0), True),  # equal q, decoy higher score
     ((1e-5, 999.0), (1e-5, 990.0), False),  # equal q, target higher score
     ((2e-5, 999.0), (1e-5, 990.0), True)],  # q first, then score
    ids=["tie", "decoy_score", "target_score", "q_first"],
)  # fmt: skip
def test_picked_ties_and_tie_break(
    t: tuple[float, float], d: tuple[float, float], decoy_wins: bool
) -> None:
    w = PickedWorld()
    x, did = w.target("intergenic", "c1", 1, q=t[0], score=t[1])
    w.decoy(did, 1, q=d[0], score=d[1])
    others = _class_block(w, "intergenic", "c1", 150)
    qs = picked_class_q(w.orfs, w.peptides, w.id_map)
    if decoy_wins:
        assert not callable_at(qs, x.orf_id, 1.0)
        assert all(q_of(qs, o) == pytest.approx(2 / 150) for o in others)
    else:
        assert all(q_of(qs, o) == pytest.approx(1 / 151) for o in [x.orf_id, *others])


# ---------------------------------------------------------------- ADR-0009 share and guard


def _share_world(n_wins: int, n_unpaired: int) -> tuple[PickedWorld, list[str], list[str]]:
    """Intergenic: ``n_wins`` paired decoys beat their targets at ranks 1..; 90 targets win
    (no decoy rows); ``n_unpaired`` intergenic-tagged decoys absent from ``id_map`` at the bottom.
    A cds_outframe class (150 winners, all decoys paired) rides along."""
    w = PickedWorld()
    for i in range(n_wins):
        _, did = w.target("intergenic", "c1", 300 + i)
        w.decoy(did, 1 + i)
    winners = _class_block(w, "intergenic", "c1", 90, first_rank=20, decoys=False)
    for _ in range(n_unpaired):
        did, _stop = w.entry("intergenic", "c9", in_map=False)
        w.decoy(did, 950)
    other = _class_block(w, "cds_outframe", "c2", 150, first_rank=100)
    return w, winners, other


@pytest.mark.parametrize(
    ("n_wins", "n_unpaired", "expected"),
    [(10, 0, 11 / 90), (10, 2, (10 / (10 / 12) + 1) / 90), (8, 2, (8 / 0.8 + 1) / 90)],
    ids=["f1", "f0.83", "f0.80_boundary"],
)
def test_decoy_wins_rescaled_by_paired_share(n_wins: int, n_unpaired: int, expected: float) -> None:
    """ADR-0009: FDR = (D/f + 1)/T; a share of exactly 80 % is still estimated."""
    w, winners, other = _share_world(n_wins, n_unpaired)
    qs = picked_class_q(w.orfs, w.peptides, w.id_map)
    for o in winners:
        assert q_of(qs, o) == pytest.approx(expected)
    for o in other:
        assert q_of(qs, o) == pytest.approx(1 / 150)


def test_paired_share_below_80_percent_falls_back_per_class() -> None:
    """Decision 4: f = 10/14 < 0.8 -> intergenic unestimated; cds_outframe keeps its q."""
    w, winners, other = _share_world(10, 4)
    qs = picked_class_q(w.orfs, w.peptides, w.id_map)
    assert all(math.isnan(q_of(qs, o)) for o in winners)
    assert all(q_of(qs, o) == pytest.approx(1 / 150) for o in other)


# ---------------------------------------------------------------- decision 2: untagged classes


@pytest.mark.parametrize("cls", ["pseudogene", "sorf", "cterm_extension"])
def test_untagged_classes_are_unestimated(cls: str) -> None:
    w = PickedWorld()
    _class_block(w, "intergenic", "c1", 150)
    rows = []
    orfs = []
    for i in range(150):
        tag = f"GSIM_{i:05d}"
        p = pep_name(10_000 + i)
        rows.append(pep_row(p, [tag], 1 + i, is_decoy=False))
        rows.append(pep_row(pep_name(20_000 + i, "DECY"), ["rev_" + tag], 700, is_decoy=True))
        orfs.append(make_orf(f"NOVEL_u_{i}", cls, "c5", "+", 1000 + 600 * i, [p],
                             overlapping_gene=tag))  # fmt: skip
    peptides = peptides_frame([*w.pep_rows, *rows])
    qs = picked_class_q([*w.orfs, *orfs], peptides, w.id_map)
    assert all(math.isnan(q_of(qs, o.orf_id)) for o in orfs)
    assert all(callable_at(qs, o.orf_id) for o in w.orfs)


# ---------------------------------------------------------------- example 6 and regimes


def _bakta_run(source: str = "database") -> tuple[pd.DataFrame, list]:  # type: ignore[type-arg]
    rows, orfs = [], []
    for i in range(150):
        tag = f"GSIM_{i:05d}"
        p = pep_name(i)
        rows.append(pep_row(p, [tag], 1 + i, is_decoy=False, source=source))
        rows.append(pep_row(pep_name(5000 + i, "DECY"), ["rev_" + tag], 600 + i, is_decoy=True))
        orfs.append(make_orf(f"NOVEL_c1_{i}", "intergenic", "c1", "+", 1300 + 600 * i, [p]))
    return peptides_frame(rows), orfs


def test_example6_bakta_decoys_are_unestimated() -> None:
    """Example 6 / degenerate row 1: ``rev_<locus_tag>`` decoys only -> unestimated."""
    peptides, orfs = _bakta_run()
    assert regime(peptides, orfs) == "unestimated"
    qs = picked_class_q(orfs, peptides, None)
    assert all(math.isnan(q_of(qs, o.orf_id)) for o in orfs)


def test_regime_picked_with_makedb_decoys() -> None:
    w = PickedWorld()
    _class_block(w, "intergenic", "c1", 5)
    assert regime(w.peptides, w.orfs) == "picked_decoy"


def test_regime_unestimated_without_decoy_rows() -> None:
    w = PickedWorld()
    _class_block(w, "intergenic", "c1", 5, decoys=False)
    assert regime(w.peptides, w.orfs) == "unestimated"


def test_regime_denovo_subgroup() -> None:
    rows = [pep_row(pep_name(i), [], 1 + i, is_decoy=False, source="denovo") for i in range(5)]
    orfs = [make_orf(f"N{i}", "intergenic", "c1", "+", 1300 + 600 * i, [pep_name(i)])
            for i in range(5)]  # fmt: skip
    assert regime(peptides_frame(rows), orfs) == "denovo_subgroup"


@pytest.mark.parametrize("makedb", [True, False])
def test_regime_per_orf_database_first(makedb: bool) -> None:
    """Decision 1: a mixed-source ORF takes the database regime when it is ``picked_decoy``,
    else de novo; a pure de novo ORF stays de novo; a database ORF follows the decoys."""
    if makedb:
        w = PickedWorld()
        _class_block(w, "intergenic", "c1", 5)
        rows, db_orf = list(w.pep_rows), w.orfs[0]
    else:
        peptides, orfs = _bakta_run()
        rows, db_orf = peptides.to_dict("records"), orfs[0]
    dn = [pep_name(90_000 + i, "DNVO") for i in range(2)]
    rows += [pep_row(p, [], 5, is_decoy=False, source="denovo") for p in dn]
    dn_orf = make_orf("N_dn", "intergenic", "c7", "+", 9000, dn)
    mixed = make_orf("N_mix", "intergenic", "c7", "+", 12000, [db_orf.peptides[0], dn[0]])
    peptides = peptides_frame(rows)
    assert regime(peptides, [db_orf]) == ("picked_decoy" if makedb else "unestimated")
    assert regime(peptides, [dn_orf]) == "denovo_subgroup"
    assert regime(peptides, [mixed]) == ("picked_decoy" if makedb else "denovo_subgroup")


# ---------------------------------------------------------------- de novo (rule 3)


@pytest.fixture(scope="module")
def tiny_gm() -> tuple[GeneModel, dict[str, str]]:
    """c1 (6 kb, linear), one plus CDS g1 at [3000, 3300): minus-strand intervals far from it
    are intergenic, plus intervals at 3001 + 3k cds_outframe, at 3003 cds_inframe."""
    rng = np.random.default_rng(1)
    seq = "".join(rng.choice(list("ACGT"), 6000))
    g = Gene("g1", "cds", Interval("c1", 3000, 3300, "+"), "M" * 99, seq[3000:3300], "ATG")
    gm = GeneModel(replicons={"c1": Replicon("c1", 6000, "linear")}, genes=(g,), features=())
    return gm, {"c1": seq}


def _dn_targets(
    cls_ivs: dict[str, list[Interval]], e0: float = 1e-6
) -> tuple[list, pd.DataFrame, list[str]]:  # type: ignore[type-arg]
    """One single-peptide de novo ORF per interval; E = e0 * (1, 2, ...) per class."""
    orfs, peps, ivs, es = [], [], [], []
    k = 0
    for cls, intervals in cls_ivs.items():
        for j, iv in enumerate(intervals):
            p = pep_name(k, "DNVO")
            k += 1
            orfs.append(make_orf(f"N_{cls}_{j}", cls, iv.contig, iv.strand,
                                 iv.end if iv.strand == "+" else iv.start - 3, [p]))  # fmt: skip
            peps.append(p)
            ivs.append(iv)
            es.append(e0 * (j + 1))
    return orfs, hit_rows(peps, ivs, es, search="target"), peps


def test_example4_decoy_counts_in_its_mirrored_class_only(
    tiny_gm: tuple[GeneModel, dict[str, str]],
) -> None:
    """The decoy row's peptide is a cds_outframe target, its mirrored interval intergenic: it
    counts in D_intergenic only. A cds_inframe decoy counts in no novel class. π̂ = 1."""
    gm, genome = tiny_gm
    clf = LocusClassifier(gm, genome)
    ig = [Interval("c1", 100 + 30 * k, 130 + 30 * k, "-") for k in range(50)]
    co = [Interval("c1", 3001 + 3 * k, 3031 + 3 * k, "+") for k in range(50)]
    assert clf.classify(ig[0])[0] == "intergenic" and clf.classify(co[-1])[0] == "cds_outframe"
    inframe = Interval("c1", 3003, 3033, "+")
    assert clf.classify(inframe)[0] == "cds_inframe"
    orfs, targets, peps = _dn_targets({"intergenic": ig, "cds_outframe": co})
    decoys = hit_rows([peps[60], peps[61]], [Interval("c1", 5000, 5030, "-"), inframe],
                      [1e-9, 1e-9], search="decoy")  # fmt: skip
    qs = denovo_class_q(orfs, targets, decoys, clf, dict.fromkeys(peps, 1.0))
    for o in orfs:
        expected = 2 / 50 if o.cls == "intergenic" else 1 / 50
        assert q_of(qs, o.orf_id) == pytest.approx(expected), o.orf_id


def test_example5_orf_q_is_min_of_its_peptides() -> None:
    """Intergenic list of 300 target peptides, 8 decoys (π̂ = 1) between E ranks 250 and 251:
    peptide q 1/250 = 0.004 for ranks <= 250, 9/300 = 0.03 after. The ORF holding ranks 1 and
    300 gets 0.004 (decision 3)."""
    clf = ContigClassifier({"ig": "intergenic"})
    peps = [pep_name(k, "DNVO") for k in range(300)]
    ivs = [Interval("ig", 100 + 30 * k, 130 + 30 * k, "-") for k in range(300)]
    e = [1e-6 * (k + 1) for k in range(300)]
    targets = hit_rows(peps, ivs, e, search="target")
    dpeps = [pep_name(k, "DCOY") for k in range(8)]
    decoys = hit_rows(dpeps, ivs[:8], [250.5e-6] * 8, search="decoy")
    combo = make_orf("N_combo", "intergenic", "ig", "-", 97, [peps[0], peps[299]])
    singles = [make_orf(f"N_{k}", "intergenic", "ig", "-", 1000 + k, [peps[k]])
               for k in range(1, 299)]  # fmt: skip
    pi = dict.fromkeys([*peps, *dpeps], 1.0)
    qs = denovo_class_q([combo, *singles], targets, decoys, clf, pi)
    assert q_of(qs, "N_combo") == pytest.approx(0.004)
    assert q_of(qs, "N_1") == pytest.approx(0.004)
    assert q_of(qs, "N_249") == pytest.approx(0.004)  # rank 250
    assert q_of(qs, "N_250") == pytest.approx(0.03)  # rank 251
    assert q_of(qs, "N_298") == pytest.approx(0.03)
    assert callable_at(qs, "N_combo") and not callable_at(qs, "N_250")


def test_denovo_class_without_decoy_hit_is_c_over_t() -> None:
    """Degenerate row: no decoy hit in a class -> q = c / T, c = max π̂ of the searched."""
    clf = ContigClassifier({"ig": "intergenic", "co": "cds_outframe"})
    orfs, targets, peps = _dn_targets({
        "intergenic": [Interval("ig", 100 + 30 * k, 130 + 30 * k, "-") for k in range(200)],
        "cds_outframe": [Interval("co", 100 + 30 * k, 130 + 30 * k, "+") for k in range(100)],
    })  # fmt: skip
    dpeps = [pep_name(k, "DCOY") for k in range(5)]
    decoys = hit_rows(dpeps, [Interval("co", 0, 30, "+")] * 5, [1e-9] * 5, search="decoy")
    pi = dict.fromkeys([*peps, *dpeps], 0.4)
    qs = denovo_class_q(orfs, targets, decoys, clf, pi)
    for o in orfs:
        expected = 0.4 / 200 if o.cls == "intergenic" else (5 * 0.4 + 0.4) / 100
        assert q_of(qs, o.orf_id) == pytest.approx(expected), o.orf_id


@pytest.mark.parametrize(("known", "expected"), [(False, 2 / 100), (True, 1 / 100)])
def test_denovo_missing_pi_counts_as_one(known: bool, expected: float) -> None:
    """Degenerate row: π̂ unavailable -> 1 (as a decoy weight and in c). 100 targets at π̂ 0.5,
    one decoy above them all: known π̂ 0.5 -> (0.5 + 0.5)/100; unknown -> (1 + 1)/100."""
    clf = ContigClassifier({"ig": "intergenic"})
    orfs, targets, peps = _dn_targets(
        {"intergenic": [Interval("ig", 100 + 30 * k, 130 + 30 * k, "-") for k in range(100)]}
    )
    decoys = hit_rows(["DCOYAAAAAAAA"], [Interval("ig", 0, 30, "-")], [1e-9], search="decoy")
    pi = dict.fromkeys(peps, 0.5)
    if known:
        pi["DCOYAAAAAAAA"] = 0.5
    qs = denovo_class_q(orfs, targets, decoys, clf, pi)
    assert all(q_of(qs, o.orf_id) == pytest.approx(expected) for o in orfs)


def test_denovo_weighted_decoys() -> None:
    """ADR-0010: W sums the decoys' own π̂; c = the largest searched π̂ (0.6 here)."""
    clf = ContigClassifier({"ig": "intergenic"})
    orfs, targets, peps = _dn_targets(
        {"intergenic": [Interval("ig", 100 + 30 * k, 130 + 30 * k, "-") for k in range(100)]}
    )
    dpeps = ["DCOYAAAAAAAA", "DCOYCCCCCCCC"]
    decoys = hit_rows(dpeps, [Interval("ig", 0, 30, "-")] * 2, [1e-9, 1e-9], search="decoy")
    pi = {**dict.fromkeys(peps, 0.2), dpeps[0]: 0.1, dpeps[1]: 0.6}
    qs = denovo_class_q(orfs, targets, decoys, clf, pi)
    assert all(q_of(qs, o.orf_id) == pytest.approx((0.1 + 0.6 + 0.6) / 100) for o in orfs)
