"""M7-11b class FDR wiring at the library level (docs/specs/M7-11b-class-fdr-wiring.md acceptance
criterion 1): ``gippa.novel.class_call`` (``class_fdr``, ``class_fdr_cell``, ``makedb_id_map``,
``denovo_inputs``), ``gippa.novel.fdr.denovo_class_fdr``, ``gippa.makedb.build.nterm_entries`` and
``gippa.novel.evidence.class_rules`` + the rule-6 call in ``apply_evidence``.

Examples 1-5 and 7-13 with their exact cells (``class_fdr`` strings, ``called``), and every
degenerate row that is not a step row. GIPPA_SPEC.md §7.6.4, §8.1, §9.3, §9.4; M7-11 decisions 5,
7-9, 11, 13, 14; lead decisions 2026-10-23 (the side-table bool ``class_mismatch`` is *not* a
``novel_orfs`` flag, decision 3).

Builders: ``helpers.class_call_fixture`` (world genome, ``assemble`` = the step body in memory)
and ``helpers.novel_fdr_fixture`` (``PickedWorld``: rank r -> ``q_value = r * 1e-5``,
``score = 1000 - r``).
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.class_call_fixture import (
    COST_COLUMNS,
    LENGTH,
    NTERM_GENE_AT,
    SUMMARY_FIELDS,
    Assembled,
    assemble,
    cc,
    config,
    costs_frame,
    decoy_hits,
    gene_model,
    nterm_world,
    orf_at,
    target_hits,
    world_genome,
)
from helpers.novel_fdr_fixture import (
    ContigClassifier,
    PickedWorld,
    db_id,
    pep_name,
    pep_row,
    peptides_frame,
)

from gippa.model.coords import Interval

GENOME = world_genome()
GM = gene_model(GENOME)
ID_MAP_COLUMNS = ("db_id", "kind", "locus_tag", "class", "contig", "start", "end", "strand")


def _asm(w: PickedWorld, alpha: float | None = None, **kw: Any) -> Assembled:
    return assemble(w.orfs, w.peptides, gm=GM, genome=GENOME, cfg=config(alpha), **kw)


def _decoy_only(w: PickedWorld, cls: str, contig: str, rank: float) -> str:
    """A makedb entry that only its decoy hit (no target ORF)."""
    did, _ = w.entry(cls, contig, in_map=False)
    w.decoy(did, rank)
    return did


def _assert_summary(s: Any, **want: Any) -> None:
    for k, v in want.items():
        got = getattr(s, k)
        if isinstance(v, float) and not isinstance(v, bool):
            assert got == pytest.approx(v), (k, got, v)
        else:
            assert got == v, (k, got, v)


# ---------------------------------------------------------------- interface shape


def test_public_names_and_dataclass_fields() -> None:
    import dataclasses

    m = cc()
    for name in ("OrfClassFdr", "ClassSummary", "DenovoInputs", "ClassFdrResult",
                 "makedb_id_map", "denovo_inputs", "class_fdr", "class_fdr_cell"):  # fmt: skip
        assert hasattr(m, name), name
    assert [f.name for f in dataclasses.fields(m.OrfClassFdr)] == [
        "orf_id", "regime", "fdr_class", "q", "reason", "class_mismatch"]  # fmt: skip
    assert tuple(f.name for f in dataclasses.fields(m.ClassSummary)) == SUMMARY_FIELDS
    assert [f.name for f in dataclasses.fields(m.ClassFdrResult)] == [
        "orfs", "classes", "n_decoy_rows", "decoys_not_makedb"]  # fmt: skip
    assert [f.name for f in dataclasses.fields(m.DenovoInputs)] == [
        "target_hits", "decoy_hits", "pi", "classifier"]  # fmt: skip


# ---------------------------------------------------------------- class_fdr_cell (rule 5)


@pytest.mark.parametrize(
    ("q", "reason", "want"),
    [(0.004, None, 0.004), (1.0, None, 1.0), (math.nan, "lost_pick", 1.0),
     *[(math.nan, r, "unestimated") for r in (
         "no_decoys", "decoys_not_makedb", "untagged_class", "no_makedb_pair",
         "paired_share_low", "decoys_unscored", "no_target_score", "no_target_hit",
         "denovo_inputs_missing")]],
)  # fmt: skip
def test_rule_5_cell_value(q: float, reason: str | None, want: float | str) -> None:
    m = cc()
    regime = "unestimated" if reason in ("no_decoys", "decoys_not_makedb") else "picked_decoy"
    r = m.OrfClassFdr("NOVEL_c1_1_p", regime, "intergenic", q, reason, False)
    got = m.class_fdr_cell(r)
    if isinstance(want, str):
        assert got == want
    else:
        assert isinstance(got, float) and not isinstance(got, bool)
        assert got == want


# ---------------------------------------------------------------- example 1


def _example1(alpha: float | None = None) -> tuple[Assembled, list[str]]:
    w = PickedWorld()
    ids = [w.target("intergenic", "c1", r)[0].orf_id for r in range(1, 151)]
    for j in range(30):
        _decoy_only(w, "intergenic", "c1", 151 + j)
    return _asm(w, alpha), ids


def test_example_1_150_paired_targets_before_30_decoy_only_entries() -> None:
    a, ids = _example1()
    for o in ids:
        assert a.cell(o) == "0.00666667", o
        assert a.called(o), o
        assert "class_fdr_unestimated" not in a.flags(o), o
        r = a.res.orfs[o]
        assert (r.regime, r.fdr_class, r.reason, r.class_mismatch) == (
            "picked_decoy", "intergenic", None, False)  # fmt: skip
        assert r.q == pytest.approx(1 / 150)
    s = a.summary("picked_decoy", "intergenic")
    _assert_summary(s, estimated=True, reason=None, n_targets=150, n_target_wins=150,
                    n_decoys=30, paired_share=1.0, c=None, q_floor=1 / 150,
                    floor_above_alpha=False)  # fmt: skip
    assert not a.res.decoys_not_makedb
    assert a.res.n_decoy_rows == 30
    assert set(a.res.orfs) == set(ids)


def test_result_classes_sorted_by_regime_then_class() -> None:
    w = PickedWorld()
    for i in range(3):
        w.target("intergenic", "c1", 10 + i)
        w.target("antisense", "c2", 20 + i)
        w.target("cds_outframe", "c3", 30 + i)
    _decoy_only(w, "antisense", "c2", 99)
    a = _asm(w)
    keys = [(s.regime, s.cls) for s in a.res.classes]
    assert keys == sorted(keys)
    assert {"antisense", "cds_outframe", "intergenic"} <= {c for _, c in keys}


# ---------------------------------------------------------------- examples 2, 3, 13 (q floor)


def _block(n: int) -> tuple[Assembled, list[str]]:
    """``n`` intergenic targets, no intergenic decoy; three antisense decoy-only entries make it
    a makedb run (the spec's "no decoys" is read as "none in this class": with no decoy row at
    all the run is ``no_decoys``)."""
    w = PickedWorld()
    ids = [w.target("intergenic", "c1", r)[0].orf_id for r in range(1, n + 1)]
    for j in range(3):
        _decoy_only(w, "antisense", "c3", 500 + j)
    return _asm(w), ids


def test_example_2_fifty_targets_floor_above_alpha() -> None:
    a, ids = _block(50)
    for o in ids:
        assert a.cell(o) == "0.02", o
        assert not a.called(o), o
        assert a.rules[o], o  # the class rule holds: q alone blocks the call
    s = a.summary("picked_decoy", "intergenic")
    _assert_summary(s, n_targets=50, n_target_wins=50, n_decoys=0, q_floor=0.02,
                    floor_above_alpha=True)  # fmt: skip


def test_example_2_decoy_only_class_has_no_floor() -> None:
    """Degenerate row "0 targets (or 0 target wins)": ``q_floor`` None, no warning."""
    a, _ = _block(50)
    s = a.summary("picked_decoy", "antisense")
    assert s.q_floor is None and s.floor_above_alpha is False


def test_example_3_hundred_targets_exactly_at_alpha() -> None:
    a, ids = _block(100)
    for o in ids:
        assert a.cell(o) == "0.01", o
        assert a.called(o), o  # one peptide each: no peptide minimum with a numeric q
    s = a.summary("picked_decoy", "intergenic")
    _assert_summary(s, n_target_wins=100, q_floor=0.01, floor_above_alpha=False)


def test_example_13_alpha_from_config() -> None:
    a, ids = _example1(alpha=0.005)
    for o in ids:
        assert a.cell(o) == "0.00666667", o
        assert not a.called(o), o
    s = a.summary("picked_decoy", "intergenic")
    _assert_summary(s, q_floor=1 / 150, floor_above_alpha=True)


# ---------------------------------------------------------------- example 4 (lost_pick)


@pytest.mark.parametrize("decoy_rank", [1, 2], ids=["decoy_better", "tie"])
def test_example_4_lost_pick_is_one_and_not_called(decoy_rank: int) -> None:
    w = PickedWorld()
    entry = w.entry("cds_outframe", "c1", "+", start0=100, length=300)
    assert entry == ("GIPPAORF_cds_outframe_c1_101_400_p", 400)
    x, did = w.target("cds_outframe", "c1", [2, 3, 4], entry=entry)
    assert x.stop == 400 and len(x.peptides) == 3
    w.decoy(did, decoy_rank)  # rank 2 = same q and score as the best target (decision 5)
    a = _asm(w)
    r = a.res.orfs[x.orf_id]
    assert r.reason == "lost_pick" and math.isnan(r.q)
    assert a.cell(x.orf_id) == "1"
    assert not a.called(x.orf_id)
    assert a.rules[x.orf_id]  # >= 2 peptides, but the 2-peptide fallback is not applied
    assert "class_fdr_unestimated" not in a.flags(x.orf_id)


def test_every_target_loses_the_pick() -> None:
    """Degenerate row: all ``1``, none called; ``n_target_wins`` 0 -> ``q_floor`` None."""
    w = PickedWorld()
    ids = []
    for i in range(3):
        o, did = w.target("cds_outframe", "c2", [10 + 2 * i, 11 + 2 * i])
        w.decoy(did, 1 + i)
        ids.append(o.orf_id)
    a = _asm(w)
    for o in ids:
        assert a.cell(o) == "1" and not a.called(o), o
        assert a.res.orfs[o].reason == "lost_pick"
    s = a.summary("picked_decoy", "cds_outframe")
    _assert_summary(s, estimated=True, n_targets=3, n_target_wins=0, n_decoys=3, q_floor=None,
                    floor_above_alpha=False)  # fmt: skip


# ---------------------------------------------------------------- example 5 (decoys_not_makedb)


def _bakta_peptides(n_decoys: int = 7) -> tuple[list[Any], pd.DataFrame]:
    a = orf_at("NOVEL_c1_1001_p", "cds_outframe", "c1", "+", 1300, ["PEPAAAAK", "PEPCCCCK"])
    b = orf_at("NOVEL_c2_2001_p", "intergenic", "c2", "+", 2300, ["PEPDDDDK"])
    rows = [pep_row(p, ["b0003"], k + 1, is_decoy=False)
            for k, p in enumerate([*a.peptides, *b.peptides])]  # fmt: skip
    rows += [pep_row(pep_name(k, "DECY"), [f"rev_b000{1 + (k >= 4)}"], 50 + k, is_decoy=True)
             for k in range(n_decoys)]  # fmt: skip
    return [a, b], peptides_frame(rows)


def test_example_5_bakta_decoys_not_makedb() -> None:
    orfs, peps = _bakta_peptides()
    a = assemble(orfs, peps, gm=GM, genome=GENOME)
    for o in orfs:
        r = a.res.orfs[o.orf_id]
        assert (r.regime, r.reason, r.fdr_class) == ("unestimated", "decoys_not_makedb", o.cls)
        assert math.isnan(r.q)
        assert a.cell(o.orf_id) == "unestimated"
        assert "class_fdr_unestimated" in a.flags(o.orf_id)
    assert a.called("NOVEL_c1_1001_p")  # cds_outframe, 2 peptides
    assert not a.called("NOVEL_c2_2001_p")  # intergenic, 1 peptide
    assert a.res.decoys_not_makedb is True
    assert a.res.n_decoy_rows == 7
    assert a.res.classes == ()


def test_no_decoy_rows_is_no_decoys() -> None:
    """Rule 4: an unestimated run without decoy rows says ``no_decoys`` (example 6's reason)."""
    orfs, peps = _bakta_peptides(n_decoys=0)
    a = assemble(orfs, peps, gm=GM, genome=GENOME)
    assert {r.reason for r in a.res.orfs.values()} == {"no_decoys"}
    assert {r.regime for r in a.res.orfs.values()} == {"unestimated"}
    assert a.res.decoys_not_makedb is False and a.res.n_decoy_rows == 0
    assert a.called("NOVEL_c1_1001_p") and not a.called("NOVEL_c2_2001_p")


def test_no_candidate_orf_still_reports_decoys_not_makedb() -> None:
    _, peps = _bakta_peptides()
    res = cc().class_fdr([], peps, gm=GM, genome=GENOME, cfg=config(), denovo=None)
    assert res.orfs == {} and res.classes == ()
    assert res.decoys_not_makedb is True and res.n_decoy_rows == 7


# ---------------------------------------------------------------- example 7 (untagged classes)


def test_example_7_untagged_classes_fall_back_to_two_peptides() -> None:
    w = PickedWorld()
    _decoy_only(w, "intergenic", "c3", 5)  # a makedb run
    s = orf_at("NOVEL_c1_5001_p", "sorf", "c1", "+", 5300, ["SORFAAAK", "SORFCCCK"])
    p = orf_at("NOVEL_c2_7001_p", "pseudogene", "c2", "+", 7300, ["PSEDAAAK"])
    rows = [pep_row(x, ["b0100"], 10 + k, is_decoy=False)
            for k, x in enumerate([*s.peptides, *p.peptides])]  # fmt: skip
    peps = peptides_frame([*w.pep_rows, *rows])
    a = assemble([s, p], peps, gm=GM, genome=GENOME)
    for o in (s, p):
        r = a.res.orfs[o.orf_id]
        assert (r.regime, r.reason, r.fdr_class) == ("picked_decoy", "untagged_class", o.cls)
        assert a.cell(o.orf_id) == "unestimated"
        assert "class_fdr_unestimated" in a.flags(o.orf_id)
    assert a.called(s.orf_id)
    assert not a.called(p.orf_id)


def test_cterm_extension_untagged_and_never_called() -> None:
    w = PickedWorld()
    _decoy_only(w, "intergenic", "c3", 5)
    t = orf_at("NOVEL_c1_5001_p", "cterm_extension", "c1", "+", 5300,
               ["CTRMAAAK", "CTRMCCCK", "CTRMDDDK"])  # fmt: skip
    rows = [pep_row(x, ["b0100"], 10 + k, is_decoy=False) for k, x in enumerate(t.peptides)]
    a = assemble([t], peptides_frame([*w.pep_rows, *rows]), gm=GM, genome=GENOME)
    assert a.res.orfs[t.orf_id].reason == "untagged_class"
    assert a.cell(t.orf_id) == "unestimated"
    assert a.rules[t.orf_id] is False
    assert not a.called(t.orf_id)


# ---------------------------------------------------------------- example 8 (class mismatch)


def test_example_8_class_mismatch_competes_in_the_entry_class() -> None:
    w = PickedWorld()
    entry = w.entry("cds_outframe", "c1", "+", start0=100, length=300)
    assert entry[0] == "GIPPAORF_cds_outframe_c1_101_400_p"
    x, _ = w.target("intergenic", "c1", 1, entry=entry)  # reconstructed intergenic, 1 peptide
    y, _ = w.target("cds_outframe", "c2", [2, 3])
    others = [w.target("intergenic", "c3", 10 + i)[0].orf_id for i in range(150)]
    _decoy_only(w, "antisense", "c4", 900)
    a = _asm(w)
    r = a.res.orfs[x.orf_id]
    assert (r.fdr_class, r.class_mismatch, r.reason) == ("cds_outframe", True, None)
    # the cds_outframe list holds X and Y: q = 1/2 (in the intergenic list it would be 1/151)
    assert r.q == pytest.approx(0.5)
    assert a.res.orfs[y.orf_id].q == pytest.approx(0.5)
    assert all(a.res.orfs[o].q == pytest.approx(1 / 150) for o in others)
    assert a.row(x.orf_id)["class"] == "intergenic"  # the TSV keeps the reconstructed class
    assert a.rules[x.orf_id] is True  # intergenic rule; the cds_outframe rule needs 2 peptides
    assert a.cell(x.orf_id) == "0.5" and not a.called(x.orf_id)
    # lead decision 3 (2026-10-23): a side-table bool only, no novel_orfs flag
    assert "class_mismatch" not in a.flags(x.orf_id)
    assert not any(a.res.orfs[o].class_mismatch for o in [y.orf_id, *others])


def test_flag_order_has_no_class_mismatch() -> None:
    """Lead decision 3: no ``class_mismatch`` flag (pending Khoa); ``FLAG_ORDER`` unchanged."""
    from gippa.novel.outputs import FLAG_ORDER

    assert "class_mismatch" not in FLAG_ORDER


# ---------------------------------------------------------------- examples 9, 10 (de novo)


def _denovo_world(
    n_orfs: int, per_orf: int, cls: str, pi: float, *, decoys: int = 0, decoy_cls: str = "c3"
) -> tuple[list[Any], pd.DataFrame, Any]:
    orfs, peps = [], []
    k = 0
    for i in range(n_orfs):
        ps = []
        for _ in range(per_orf):
            ps.append(pep_name(k, "DNVO"))
            k += 1
        stop = 1000 + 600 * i
        orfs.append(orf_at(f"NOVEL_c1_{stop - 299}_p", cls, "c1", "+", stop, ps))
        peps += ps
    dps = [pep_name(j, "DCOY") for j in range(decoys)]
    rows = [pep_row(p, [], 0, is_decoy=False, source="denovo", q=math.nan, score=50.0)
            for p in peps]  # fmt: skip
    costs = pd.concat([
        costs_frame(peps, e_target=1e-3, pi=pi),
        costs_frame(dps, e_target=math.nan, e_decoy=1e-2, pi=pi) if dps else None,
    ], ignore_index=True)  # fmt: skip
    hits = pd.concat([
        target_hits(peps, [Interval("c1", 100 + 30 * j, 130 + 30 * j, "+")
                           for j in range(len(peps))]),
        decoy_hits(dps, [Interval(decoy_cls, 100 + 60 * j, 130 + 60 * j, "-")
                         for j in range(len(dps))]),
    ], ignore_index=True)  # fmt: skip
    clf = ContigClassifier({"c1": cls, "c2": cls, "c3": "antisense", "c4": "cds_outframe"})
    d = cc().denovo_inputs(costs, hits, clf)
    return orfs, peptides_frame(rows), d


def test_example_9_denovo_c_over_t() -> None:
    orfs, peps, d = _denovo_world(120, 1, "intergenic", 0.5, decoys=5, decoy_cls="c3")
    a = assemble(orfs, peps, gm=GM, genome=GENOME, denovo=d)
    for o in orfs:
        r = a.res.orfs[o.orf_id]
        assert (r.regime, r.fdr_class, r.reason) == ("denovo_subgroup", "intergenic", None)
        assert a.cell(o.orf_id) == "0.00416667", o.orf_id
        assert a.called(o.orf_id), o.orf_id
    s = a.summary("denovo_subgroup", "intergenic")
    _assert_summary(s, estimated=True, reason=None, n_targets=120, n_target_wins=None,
                    n_decoys=0, paired_share=None, c=0.5, q_floor=0.5 / 120,
                    floor_above_alpha=False)  # fmt: skip


def test_example_10_denovo_floor_with_missing_pi() -> None:
    orfs, peps, d = _denovo_world(20, 2, "antisense", math.nan)
    assert dict(d.pi) == {}  # NaN dropped: missing = 1
    a = assemble(orfs, peps, gm=GM, genome=GENOME, denovo=d)
    s = a.summary("denovo_subgroup", "antisense")
    _assert_summary(s, n_targets=40, n_decoys=0, c=1.0, q_floor=0.025, floor_above_alpha=True)
    for o in orfs:
        assert a.cell(o.orf_id) == "0.025", o.orf_id
        assert a.rules[o.orf_id] and not a.called(o.orf_id), o.orf_id


def test_denovo_inputs_missing_reason() -> None:
    orfs, peps, _ = _denovo_world(3, 1, "intergenic", 0.5)
    a = assemble(orfs, peps, gm=GM, genome=GENOME, denovo=None)
    for o in orfs:
        r = a.res.orfs[o.orf_id]
        assert (r.regime, r.reason) == ("denovo_subgroup", "denovo_inputs_missing")
        assert a.cell(o.orf_id) == "unestimated"
    assert not any(s.regime == "denovo_subgroup" for s in a.res.classes)


def test_denovo_orf_without_finite_e_target_is_no_target_hit() -> None:
    orfs, peps, _ = _denovo_world(3, 1, "intergenic", 0.5)
    pep_all = [p for o in orfs for p in o.peptides]
    e = [1e-3, 1e-3, math.nan]
    costs = costs_frame(pep_all, e_target=e, pi=0.5)
    hits = target_hits(pep_all[:2], [Interval("c1", 100, 130, "+")] * 2)
    d = cc().denovo_inputs(costs, hits, ContigClassifier({"c1": "intergenic"}))
    a = assemble(orfs, peps, gm=GM, genome=GENOME, denovo=d)
    last = next(o for o in orfs if pep_all[2] in o.peptides)
    assert a.res.orfs[last.orf_id].reason == "no_target_hit"
    assert a.cell(last.orf_id) == "unestimated"
    for o in orfs:
        if o is not last:
            assert a.res.orfs[o.orf_id].reason is None


def test_nan_e_decoy_on_a_best_cost_decoy_row_raises() -> None:
    orfs, peps, _ = _denovo_world(2, 1, "intergenic", 0.5)
    pep_all = [p for o in orfs for p in o.peptides]
    costs = costs_frame(pep_all, e_target=1e-3, e_decoy=[math.nan, 1e-2], pi=0.5)
    costs["m_decoy"] = 0.0  # a decoy hit exists, its E is missing
    hits = decoy_hits(pep_all, [Interval("c1", 100, 130, "-")] * 2)
    with pytest.raises(ValueError):
        d = cc().denovo_inputs(costs, hits, ContigClassifier({"c1": "intergenic"}))
        cc().class_fdr(orfs, peps, gm=GM, genome=GENOME, cfg=config(), denovo=d)


# ---------------------------------------------------------------- denovo_inputs (rule 3, Q11)


def test_denovo_inputs_keeps_only_denovo_peptides_and_attaches_e() -> None:
    dn_t, dn_d, db_t = ["DNVOAAAK", "DNVOCCCK"], ["DNVODDDK"], ["DBPPAAAK"]
    costs = pd.concat([
        costs_frame(dn_t, e_target=[1e-3, 2e-3], e_decoy=[math.nan, 5e-3], pi=[0.3, math.nan]),
        costs_frame(dn_d, e_target=math.nan, e_decoy=7e-3, pi=0.6),
        costs_frame(db_t, e_target=1e-4, e_decoy=1e-4, pi=math.nan, source="database"),
    ], ignore_index=True).sort_values("peptide").reset_index(drop=True)  # fmt: skip
    costs.loc[costs["peptide"] == "DNVOCCCK", "m_decoy"] = 0.0
    hits = pd.concat([
        target_hits([*dn_t, *db_t], [Interval("c1", 100, 130, "+")] * 3),
        decoy_hits(["DNVOCCCK", "DNVODDDK", "DBPPAAAK"], [Interval("c3", 400, 430, "-")] * 3),
    ], ignore_index=True)  # fmt: skip
    clf = ContigClassifier({"c1": "intergenic", "c3": "antisense"})
    d = cc().denovo_inputs(costs, hits, clf)
    assert d.classifier is clf
    th = d.target_hits.set_index("peptide")["e_target"].to_dict()
    assert th == {"DNVOAAAK": pytest.approx(1e-3), "DNVOCCCK": pytest.approx(2e-3)}
    dh = d.decoy_hits
    assert set(dh["peptide"]) == {"DNVOCCCK", "DNVODDDK"}  # no database peptide's decoy (Q11)
    assert (dh["search"] == "decoy").all()
    assert {"contig", "start", "end", "strand", "n_diff"} <= set(dh.columns)
    assert dh.set_index("peptide")["e_decoy"].to_dict() == {
        "DNVOCCCK": pytest.approx(5e-3), "DNVODDDK": pytest.approx(7e-3)}  # fmt: skip
    assert dict(d.pi) == {"DNVOAAAK": pytest.approx(0.3), "DNVODDDK": pytest.approx(0.6)}


# ---------------------------------------------------------------- example 11 (nterm by locus tag)


def _example11(with_gene: bool) -> tuple[Assembled, str, str]:
    from gippa.makedb.build import _nterm_entries

    genome, gene = nterm_world()
    gm = gene_model(genome, [gene] if with_gene else [])
    full = gene_model(genome, [gene])
    recs = [e for e in _nterm_entries(full, genome) if e.locus_tag == "b0002"]
    assert len(recs) == 1
    entry = recs[0].db_id
    assert entry.startswith("GIPPAORF_nterm_extension_c1_") and entry.endswith("_p")
    p = NTERM_GENE_AT
    o = orf_at(f"NOVEL_c1_{p - 56}_p", "nterm_extension", "c1", "+", p + 297,
               ["NTRMAAAK", "NTRMCCCK"], start=p - 57, overlapping_gene="b0002")  # fmt: skip
    rows = [pep_row(x, [entry], 1 + k, is_decoy=False) for k, x in enumerate(o.peptides)]
    rows.append(pep_row("DECYAAAK", ["rev_" + db_id("intergenic", "c2", 1000, 1300, "+")], 9,
                        is_decoy=True))  # fmt: skip
    a = assemble([o], peptides_frame(rows), gm=gm, genome=genome)
    return a, o.orf_id, entry


def test_example_11_nterm_paired_by_locus_tag() -> None:
    a, oid, _ = _example11(with_gene=True)
    r = a.res.orfs[oid]
    assert (r.regime, r.fdr_class, r.reason) == ("picked_decoy", "nterm_extension", None)
    assert math.isfinite(r.q) and r.q == pytest.approx(1.0)  # one target win: (0 + 1) / 1
    assert a.cell(oid) == "1"
    assert "class_fdr_unestimated" not in a.flags(oid)


def test_example_11_gene_model_without_the_gene_is_no_makedb_pair() -> None:
    a, oid, _ = _example11(with_gene=False)
    r = a.res.orfs[oid]
    assert (r.regime, r.reason) == ("picked_decoy", "no_makedb_pair")
    assert a.cell(oid) == "unestimated"
    assert "class_fdr_unestimated" in a.flags(oid)


def test_nterm_entries_is_the_public_alias() -> None:
    from gippa.makedb.build import _nterm_entries, nterm_entries

    genome, gene = nterm_world()
    gm = gene_model(genome, [gene])
    assert list(nterm_entries(gm, genome)) == list(_nterm_entries(gm, genome))


# ---------------------------------------------------------------- makedb_id_map (rule 2)


def test_makedb_id_map_rows_from_tagged_ids() -> None:
    genome, gene = nterm_world()
    gm = gene_model(genome, [gene])
    from gippa.makedb.build import _nterm_entries

    nterm = next(e.db_id for e in _nterm_entries(gm, genome))
    stray = "GIPPAORF_nterm_extension_c2_101_400_p"  # no nterm_entries record
    rows = [
        pep_row("PEPAAAAK", ["GIPPAORF_cds_outframe_c1_101_400_p", "b0001"], 1, is_decoy=False),
        pep_row("PEPCCCCK", ["GIPPAORF_cds_outframe_c1_101_400_p"], 2, is_decoy=False),
        pep_row("PEPDDDDK", [nterm, stray], 3, is_decoy=False),
        pep_row("DECYAAAK", ["rev_GIPPAORF_intergenic_c2_1001_1300_m", "rev_b0002"], 4,
                is_decoy=True),
        pep_row("DECYCCCK", ["DECOY_GIPPAORF_antisense_c3_7_306_p"], 5, is_decoy=True),
    ]  # fmt: skip
    peps = peptides_frame(rows)
    m = cc().makedb_id_map(peps, gm, genome, ["rev_", "DECOY_"])
    assert tuple(m.columns) == ID_MAP_COLUMNS
    got = {r["db_id"]: r for r in m.to_dict("records")}
    assert set(got) == {
        "GIPPAORF_cds_outframe_c1_101_400_p",
        nterm,
        stray,
        "GIPPAORF_intergenic_c2_1001_1300_m",
        "GIPPAORF_antisense_c3_7_306_p",
    }
    assert len(m) == len(got)  # one row per distinct ID
    assert (m["kind"] == "orf").all()
    r = got["GIPPAORF_cds_outframe_c1_101_400_p"]
    assert (r["class"], r["contig"], int(r["start"]), int(r["end"]), r["strand"]) == (
        "cds_outframe", "c1", 101, 400, "+")  # fmt: skip
    assert pd.isna(r["locus_tag"])
    r = got["GIPPAORF_intergenic_c2_1001_1300_m"]
    assert (int(r["start"]), int(r["end"]), r["strand"]) == (1001, 1300, "-")
    assert got[nterm]["locus_tag"] == "b0002"
    assert pd.isna(got[stray]["locus_tag"])
    for k in ("GIPPAORF_intergenic_c2_1001_1300_m", "GIPPAORF_antisense_c3_7_306_p"):
        assert pd.isna(got[k]["locus_tag"])


# ---------------------------------------------------------------- example 12 (mixed run)


def test_example_12_mixed_run_lists_not_pooled() -> None:
    w = PickedWorld()
    a_orf, did = w.target("intergenic", "c1", [1, 2])
    w.decoy(did, 40)
    a_dn = "DNVOAAAK"  # A also carries a de novo peptide: A stays picked (database first)
    w.pep_rows.append(pep_row(a_dn, [], 0, is_decoy=False, source="denovo", q=math.nan))
    a_orf = orf_at(a_orf.orf_id, "intergenic", "c1", "+", a_orf.stop, [*a_orf.peptides, a_dn])
    b_pep = "DNVOCCCK"
    b = orf_at("NOVEL_c2_90001_p", "intergenic", "c2", "+", 90300, [b_pep])
    w.pep_rows.append(pep_row(b_pep, [], 0, is_decoy=False, source="denovo", q=math.nan))
    db_peps = list(w.orfs[0].peptides)
    costs = pd.concat([
        costs_frame([a_dn, b_pep], e_target=[1e-4, 1e-3], pi=[0.3, 0.5]),
        costs_frame(db_peps, e_target=1e-5, e_decoy=1e-5, source="database"),
    ], ignore_index=True)  # fmt: skip
    hits = pd.concat([
        target_hits([a_dn, b_pep, *db_peps], [Interval("c1", 100, 130, "+")] * 4),
        decoy_hits(db_peps[:1], [Interval("c2", 500, 530, "-")]),
    ], ignore_index=True)  # fmt: skip
    clf = ContigClassifier({"c1": "intergenic", "c2": "intergenic"})
    d = cc().denovo_inputs(costs, hits, clf)
    assert db_peps[0] not in set(d.decoy_hits["peptide"])
    a = assemble([a_orf, b], w.peptides, gm=GM, genome=GENOME, denovo=d)
    assert a.res.orfs[a_orf.orf_id].regime == "picked_decoy"
    assert a.res.orfs[b.orf_id].regime == "denovo_subgroup"
    pk = a.summary("picked_decoy", "intergenic")
    dn = a.summary("denovo_subgroup", "intergenic")
    assert pk.n_targets == 1
    assert dn.n_targets == 1  # B's peptide only: A's de novo peptide is not in the de novo list
    assert dn.c == pytest.approx(0.5)
    assert a.res.orfs[b.orf_id].q == pytest.approx(0.5)  # c / T = 0.5 / 1


# ---------------------------------------------------------------- picked degenerate rows


def test_unscored_decoys_make_the_class_unestimated() -> None:
    w = PickedWorld()
    ids = []
    for i in range(2):
        o, did = w.target("cds_outframe", "c2", [10 + 2 * i, 11 + 2 * i])
        w.decoy(did, 50 + i, q=math.nan, score=math.nan)
        ids.append(o.orf_id)
    a = _asm(w)
    for o in ids:
        assert a.res.orfs[o].reason == "decoys_unscored"
        assert a.cell(o) == "unestimated"
        assert a.called(o)  # the 2-peptide rule and the cds_outframe rule
    s = a.summary("picked_decoy", "cds_outframe")
    _assert_summary(s, estimated=False, reason="decoys_unscored", q_floor=None,
                    floor_above_alpha=False)  # fmt: skip


def test_decoys_without_q_rank_on_score() -> None:
    w = PickedWorld()
    ids = []
    for i in range(2):
        o, did = w.target("cds_outframe", "c2", [10 + 2 * i, 11 + 2 * i])
        w.decoy(did, 50 + i, q=math.nan)  # score 950, 949: below every target
        ids.append(o.orf_id)
    a = _asm(w)
    for o in ids:
        assert a.res.orfs[o].reason is None
        assert a.res.orfs[o].q == pytest.approx(0.5)
        assert a.cell(o) == "0.5"


def test_target_row_with_a_decoy_id_is_a_target() -> None:
    """An aggregated peptide that is target in one row and decoy in another (``is_decoy`` false
    after ``all``): its ``rev_`` ID is ignored, so X is not beaten by its own peptide."""
    w = PickedWorld()
    entry = w.entry("cds_outframe", "c1", "+", start0=100, length=300)
    x, did = w.target("cds_outframe", "c1", 3, entry=entry)
    w.pep_rows[-1]["protein_ids"] = [did, "rev_" + did]
    _decoy_only(w, "antisense", "c4", 50)
    a = _asm(w)
    r = a.res.orfs[x.orf_id]
    assert r.reason is None and r.q == pytest.approx(1.0)


def test_two_orfs_paired_to_one_entry_share_its_q() -> None:
    w = PickedWorld()
    entry = w.entry("cds_outframe", "c1", "+", start0=100, length=300)
    x1, _ = w.target("cds_outframe", "c1", 1, entry=entry)
    x2, _ = w.target("cds_outframe", "c1", 5, entry=entry)
    w.target("cds_outframe", "c2", 2)
    w.target("cds_outframe", "c2", 3)
    _decoy_only(w, "antisense", "c4", 50)
    a = _asm(w)
    q1, q2 = a.res.orfs[x1.orf_id].q, a.res.orfs[x2.orf_id].q
    assert q1 == pytest.approx(1 / 3) and q2 == pytest.approx(1 / 3)  # three list members
    assert a.summary("picked_decoy", "cds_outframe").n_target_wins == 3


def test_origin_crossing_orf_pairs_modulo_the_replicon_length() -> None:
    w = PickedWorld()
    did = db_id("intergenic", "c1", LENGTH - 150, LENGTH + 150, "+")
    assert did == f"GIPPAORF_intergenic_c1_{LENGTH - 149}_{LENGTH + 150}_p"
    x = orf_at(f"NOVEL_c1_{LENGTH - 149}_p", "intergenic", "c1", "+", 150, ["ORGNAAAK"],
               start=LENGTH - 150)  # fmt: skip
    w.pep_rows.append(pep_row("ORGNAAAK", [did], 1, is_decoy=False))
    _decoy_only(w, "antisense", "c4", 50)
    res = cc().class_fdr([x], w.peptides, gm=GM, genome=GENOME, cfg=config(), denovo=None)
    r = res.orfs[x.orf_id]
    assert r.reason is None and r.q == pytest.approx(1.0)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_open_3p_orf_on_a_linear_contig(strand: str) -> None:
    """Decisions 13-14: key ``stop +/- 3``, wrapped like the entry key."""
    genome = {"c5": world_genome(("c5",), 3000, seed=5)["c5"], **GENOME}
    gm = gene_model(genome, topology="linear")
    if strand == "+":
        x = orf_at("NOVEL_c5_2698_p", "intergenic", "c5", "+", 2997, ["OPENAAAK"], start=2697,
                   end_3p_open=True)  # fmt: skip
        did = db_id("intergenic", "c5", 2697, 3000, "+")
    else:
        x = orf_at("NOVEL_c5_301_m", "intergenic", "c5", "-", 0, ["OPENAAAK"], start=300,
                   end_3p_open=True)  # fmt: skip
        did = db_id("intergenic", "c5", 0, 303, "-")
    w = PickedWorld()
    w.pep_rows.append(pep_row("OPENAAAK", [did], 1, is_decoy=False))
    _decoy_only(w, "antisense", "c4", 50)
    res = cc().class_fdr([x], w.peptides, gm=gm, genome=genome, cfg=config(), denovo=None)
    r = res.orfs[x.orf_id]
    assert r.reason is None and r.q == pytest.approx(1.0), r


def test_minus_strand_orf_pairs_by_its_stop() -> None:
    w = PickedWorld()
    ids = [w.target("antisense", "c2", 1 + i, strand="-")[0].orf_id for i in range(4)]
    _decoy_only(w, "intergenic", "c4", 50)
    a = _asm(w)
    for o in ids:
        assert a.res.orfs[o].reason is None
        assert a.cell(o) == "0.25"


def test_no_cascade_from_picked_to_denovo() -> None:
    """Q6: a picked-regime ORF that is unestimated (``no_makedb_pair``) keeps its regime even
    though it also has a de novo peptide with a target hit."""
    w = PickedWorld()
    _decoy_only(w, "intergenic", "c3", 5)
    dn = "DNVOAAAK"
    x = orf_at("NOVEL_c1_50001_p", "cds_outframe", "c1", "+", 50300, ["DBPPAAAK", dn])
    w.pep_rows.append(pep_row("DBPPAAAK", ["b0007"], 2, is_decoy=False))
    w.pep_rows.append(pep_row(dn, [], 0, is_decoy=False, source="denovo", q=math.nan))
    costs = costs_frame([dn], e_target=1e-4, pi=0.01)
    d = cc().denovo_inputs(costs, target_hits([dn], [Interval("c1", 50100, 50130, "+")]),
                           ContigClassifier({"c1": "cds_outframe"}))  # fmt: skip
    res = cc().class_fdr([x], w.peptides, gm=GM, genome=GENOME, cfg=config(), denovo=d)
    r = res.orfs[x.orf_id]
    assert (r.regime, r.reason) == ("picked_decoy", "no_makedb_pair")
    assert math.isnan(r.q)


# ---------------------------------------------------------------- class_rules + call (rule 6)


def _rules_case(cls: str, n_peps: int, *, tir: float) -> tuple[bool, bool]:
    peps = [pep_name(k, "RULE") for k in range(n_peps)]
    o = orf_at("NOVEL_c1_1001_p", cls, "c1", "+", 1300, peps)
    rows = [pep_row(p, ["b0003"], 1 + k, is_decoy=False) for k, p in enumerate(peps)]
    a = assemble([o], peptides_frame(rows), gm=GM, genome=GENOME, tir=tir)
    return a.rules[o.orf_id], a.called(o.orf_id)


@pytest.mark.parametrize(
    ("cls", "n_peps", "tir", "rule"),
    [("cds_outframe", 1, 1e6, False), ("cds_outframe", 2, 1e6, True),
     ("antisense", 1, 1e6, False), ("antisense", 3, 1e6, True),
     ("intergenic", 2, 1e6, True), ("intergenic", 2, math.nan, False),
     ("sorf", 1, math.nan, True), ("pseudogene", 2, math.nan, True),
     ("cterm_extension", 3, 1e6, False)],
)  # fmt: skip
def test_class_rules_is_the_m7_09_rule_4(cls: str, n_peps: int, tir: float, rule: bool) -> None:
    got_rule, called = _rules_case(cls, n_peps, tir=tir)
    assert got_rule is rule
    # unestimated cell (no decoys): called = rule and n_peptides >= 2, as M7-09 rule 8
    assert called is (rule and n_peps >= 2)


def test_class_rules_covers_every_row() -> None:
    w = PickedWorld()
    for i in range(3):
        w.target("intergenic", "c1", 1 + i)
    a = _asm(w)
    assert set(a.rules) == set(a.pre["orf_id"])
    assert all(isinstance(v, bool) for v in a.rules.values())


def test_apply_evidence_keeps_the_cell_it_receives() -> None:
    a, _ = _block(50)
    assert set(a.pre["class_fdr"]) == {"0.02"}
    assert set(a.table["class_fdr"]) == {"0.02"}  # no longer overwritten with "unestimated"


def test_rank_covers_called_orfs_only() -> None:
    a, ids = _example1()
    ranks = sorted(int(a.row(o)["rank"]) for o in ids)
    assert ranks == list(range(1, 151))
    b, ids2 = _example1(alpha=0.005)
    assert all(pd.isna(b.row(o)["rank"]) for o in ids2)


# ---------------------------------------------------------------- denovo_class_fdr (fdr.py)


def test_denovo_class_fdr_per_class_stats_and_q() -> None:
    from helpers.novel_fdr_fixture import hit_rows

    from gippa.novel import fdr

    t = hit_rows(["TGTAAAAK", "TGTCCCCK"], [Interval("c1", 100, 130, "+")] * 2, [1e-4, 1e-3],
                 search="target")  # fmt: skip
    d = hit_rows(["DCYAAAAK", "DCYAAAAK", "DCYCCCCK", "DCYDDDDK"],
                 [Interval("c1", 400, 430, "-"), Interval("c3", 400, 430, "-"),
                  Interval("c1", 500, 530, "-"), Interval("c3", 600, 630, "-")],
                 [5e-4, 5e-4, 2e-3, 1e-4], search="decoy")  # fmt: skip
    d.loc[1, "n_diff"] = 1  # DCYAAAAK: only its best-cost row (intergenic) counts
    clf = ContigClassifier({"c1": "intergenic", "c3": "antisense"})
    pi = {"TGTAAAAK": 0.2, "DCYAAAAK": 0.4}  # TGTCCCCK, DCYCCCCK, DCYDDDDK missing = 1
    orfs = [
        orf_at("NOVEL_c1_1_p", "intergenic", "c1", "+", 300, ["TGTAAAAK"]),
        orf_at("NOVEL_c1_601_p", "intergenic", "c1", "+", 900, ["TGTCCCCK"]),
    ]
    out = fdr.denovo_class_fdr(orfs, t, d, clf, pi)
    assert out.q == fdr.denovo_class_q(orfs, t, d, clf, pi)
    st = out.classes["intergenic"]
    assert (st.n_targets, st.n_decoys) == (2, 2)
    assert st.w_decoys == pytest.approx(0.4 + 1.0)
    assert st.c == pytest.approx(1.0)
    assert set(COST_COLUMNS) >= {"e_target", "e_decoy", "pi_denovo"}  # layout guard
    assert np.isfinite(list(out.q.values())).all()
