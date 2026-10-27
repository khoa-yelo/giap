"""Pure TE2-TE5 functions of ``validation/pilot/te_real.py`` (M10-06 acceptance criteria 1-5).

Spec: docs/specs/M10-06-te-real-and-te7.md ("Definitions" TE2-TE5, examples 1-6, edge cases
"Empty input", "Circular and linear contigs", "Minus strand", "Ties and boundaries", "Duplicate
names"). Every share is compared with ``gippa.stats.cutoffs.wilson_ci`` to 1e-12 and with the digits
the spec prints. Inputs are hand-made frames, the committed public Nissle Bakta annotation
(criterion 2) and ``tests/data/contaminants_mini.faa`` (criterion 5).

Interface fixed here (test-author reading of the spec "Interface"; the spec names the functions
and their arguments, these tests fix the returned keys):

* every TE result is a dict with ``status`` in {pass, fail, not_evaluable} and, when not
  ``pass``, ``reason`` (a list of strings); a **share** is ``{"value", "ci_low", "ci_high", ...}``
  as in M10-05; ``warnings`` is a list of codes (strings or dicts with ``code``);
* ``te2(variants, genes, detection_summary)``: ``variants`` is a frame or None (file missing);
  returns ``counts`` (te2_result's five counts), ``violations`` (te2_result's first 50 examples),
  ``vacuous``, ``breakdown`` {``alt_type``, ``call``, ``flags``, ``n_alt``, ``replicon``,
  ``region``: {key: count}, ``aa_pos_1_unflagged``: int}, ``warnings``; no ``notes`` key;
  ``flags`` keys are the sorted other flags joined by ``;``; ``n_alt`` keys "1", "2" and one
  ``>= 3`` bin; ``region`` keys are ``region_ids`` prefixes without the ``_<n>`` suffix;
* ``ribosomal_set(genes)``: a frame with ``locus_tag`` or an iterable of locus tags;
* ``te3(genes)``: ``n_ribosomal``, ``detection`` {``share``, ``holds``}, ``abundance`` {``m``,
  ``q3``, ``n_r_measured``, ``n_measured``, ``percentile_rank``, ``holds``, ``abundance_type``},
  ``secondary`` {``without_type_b``, ``section_7_1_2``}: each {``n_ribosomal``, ``detection``,
  ``abundance``} as above; ``reason`` items ``detection``, ``abundance``,
  ``no_ribosomal_genes``;
* ``te4(genes)``: ``pairs`` {"rplL/rplJ", "atpA/atpD"}: each {``status``, ``reason``, ``low``,
  ``high``, ``log2_ratio``, ``top3_log2_ratio``}; ``abundance_type``; ``secondary``
  {``genes``: {name: {``n_peptides_unique``, ``n_observable_peptides``}}};
* ``contaminant_families(headers)``: header lines without ``>``; returns family -> set of IDs
  (``gippa.io.protein_fasta.resolve_protein_id`` of the first token);
* ``te5_arm(P, peptide_hits, unmapped, families, cfg)``: ``n_s``, ``families`` {family: n},
  ``counts`` {``genome_exact``, ``bucket``, ``match_other_entry``, ``placed_tolerant``,
  ``labelled_unmapped``, ``absent``}, ``share_bucket`` (a share), ``secondary``
  {``keratin_associated``: {``n_s``, ``counts``}}.
"""

from __future__ import annotations

import gzip
import json
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.pilot_metrics_m10 import load_te_proxies
from helpers.te_real_m10 import (
    ALBUMIN_PEPTIDE,
    EX1_GENES,
    EX6_FAMILIES,
    MINI_CONTAMINANTS,
    NISSLE_JSON,
    codes,
    ds_frame,
    ex1_variants,
    ex2_extra,
    ex6_inputs,
    gframe,
    grow,
    hits_frame,
    load_te_real,
    pframe,
    prow,
    tag_set,
    unmapped_frame,
    vframe,
    vrow,
)
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.config import resolve_config
from gippa.stats.cutoffs import wilson_ci

TOL = 1e-12


@pytest.fixture(scope="module")
def tr() -> Any:
    return load_te_real()


@pytest.fixture(scope="module")
def cfg() -> Any:
    return resolve_config()


def assert_wilson(s: Mapping[str, Any], k: int, n: int, printed: tuple[float, ...] = ()) -> None:
    wlo, whi = wilson_ci(k, n)
    assert s["value"] == pytest.approx(k / n, abs=TOL)
    assert s["ci_low"] == pytest.approx(wlo, abs=TOL)
    assert s["ci_high"] == pytest.approx(whi, abs=TOL)
    if printed:
        assert round(s["ci_low"], 3) == pytest.approx(printed[0])
        assert round(s["ci_high"], 3) == pytest.approx(printed[1])


# ------------------------------------------------------------------ criterion 1: TE2


def test_te2_example_1_pass(tr: Any) -> None:
    r = tr.te2(vframe(ex1_variants()), gframe(EX1_GENES), ds_frame(1))
    assert r["counts"] == {
        "alt_sites": 3,
        "alt_stop_codon": 1,
        "high": 2,
        "high_start_codon": 1,
        "high_stop_codon": 1,
        "violations": 0,
    }
    assert r["status"] == "pass"
    assert not r.get("vacuous")
    assert "te2_id_mixups_inconsistent" not in codes(r.get("warnings", []))
    assert "notes" not in r  # the proxy's in-silico notes are dropped


def test_te2_example_2_fail_and_breakdowns(tr: Any) -> None:
    v = vframe([*ex1_variants(), ex2_extra()])
    r = tr.te2(v, gframe(EX1_GENES), ds_frame(2))
    assert r["counts"]["violations"] == 1
    assert r["status"] == "fail" and r["reason"]
    b = r["breakdown"]
    assert b["alt_type"] == {"substitution": 1}
    assert b["call"] == {"ALT": 1}
    assert b["flags"] == {"database_derived;mass_equivalent": 1}
    assert b["n_alt"] == {"2": 1}
    assert b["replicon"] == {"plasmid": 1}
    assert b["region"] == {"prophage": 1}
    assert b["aa_pos_1_unflagged"] == 0
    assert [x["locus_tag"] for x in r["violations"]] == ["JMJJOG_04600"]
    assert r["violations"][0]["aa_pos"] == 57  # protein coordinate as written (minus strand)
    assert "te2_id_mixups_inconsistent" not in codes(r.get("warnings", []))


def test_te2_reuses_te2_result(tr: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """Criterion 1: ``te2`` calls ``te_proxies.te2_result`` (spy), no copy of its rule."""
    tp = load_te_proxies()
    real = tp.te2_result
    calls: list[int] = []

    def spy(sites: pd.DataFrame) -> Any:
        calls.append(len(sites))
        return real(sites)

    monkeypatch.setattr(tp, "te2_result", spy)
    if hasattr(tr, "te2_result"):
        monkeypatch.setattr(tr, "te2_result", spy)
    tr.te2(vframe(ex1_variants()), gframe(EX1_GENES), ds_frame(1))
    assert calls == [3]


def test_te2_cross_check_warns_when_violations_exceed_id_mixups(tr: Any) -> None:
    v = vframe([*ex1_variants(), ex2_extra()])
    r = tr.te2(v, gframe(EX1_GENES), ds_frame(0))
    assert "te2_id_mixups_inconsistent" in codes(r["warnings"])
    assert r["status"] == "fail"  # a warning, never the verdict


def test_te2_violations_equal_to_id_mixups_is_no_warning(tr: Any) -> None:
    """Surviving mutant 10: the cross-check is ``violations <= id_mixups``; equality is fine."""
    v = vframe([*ex1_variants(), ex2_extra()])
    r = tr.te2(v, gframe(EX1_GENES), ds_frame(1))
    assert r["counts"]["violations"] == 1
    assert "te2_id_mixups_inconsistent" not in codes(r["warnings"])


def test_te2_breakdown_disagreeing_with_te2_result_is_inconsistent(
    tr: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    """S5: the breakdown rows must be as many as te2_result's ``counts.violations``; otherwise
    an internal inconsistency (exit 1), never a silently drifting breakdown."""
    import dataclasses

    from validation.pilot.metrics import Inconsistent

    tp = load_te_proxies()
    real = tp.te2_result

    def drifted(sites: pd.DataFrame) -> Any:
        res = real(sites)
        counts = {**res.counts, "violations": res.counts["violations"] + 1}
        return dataclasses.replace(res, counts=counts)

    monkeypatch.setattr(tp, "te2_result", drifted)
    if hasattr(tr, "te2_result"):
        monkeypatch.setattr(tr, "te2_result", drifted)
    with pytest.raises(Inconsistent):
        tr.te2(vframe([*ex1_variants(), ex2_extra()]), gframe(EX1_GENES), ds_frame(2))


def test_te2_gene_not_in_genes_is_not_in_genes_for_replicon_and_region(tr: Any) -> None:
    """N1: a violation whose gene has no genes.tsv row gets the same token in both breakdowns."""
    v = vframe([vrow("v9", "JMJJOG_09999", 12, "N", "D", "high", "database_derived")])
    r = tr.te2(v, gframe(EX1_GENES), ds_frame(1))
    assert r["counts"]["violations"] == 1
    assert r["breakdown"]["replicon"] == {"not_in_genes": 1}
    assert r["breakdown"]["region"] == {"not_in_genes": 1}


def test_te2_id_mixups_na(tr: Any) -> None:
    r = tr.te2(vframe(ex1_variants()), gframe(EX1_GENES), ds_frame(None))
    assert "id_mixups_na" in codes(r["warnings"])
    assert r["status"] == "pass"


def test_te2_missing_variants_is_not_evaluable(tr: Any) -> None:
    r = tr.te2(None, gframe(EX1_GENES), ds_frame(None))
    assert r["status"] == "not_evaluable"
    assert "missing:tables/variants.tsv" in r["reason"]


def test_te2_empty_variants_is_a_vacuous_pass(tr: Any) -> None:
    r = tr.te2(vframe([]), gframe(EX1_GENES), ds_frame(0))
    assert r["status"] == "pass" and r["vacuous"] is True
    assert r["counts"]["alt_sites"] == 0
    assert "notes" not in r


def test_te2_no_high_call_but_rows_is_not_vacuous(tr: Any) -> None:
    v = vframe([vrow("v3", "JMJJOG_00030", 41, "N", "D", "candidate", "database_derived")])
    r = tr.te2(v, gframe(EX1_GENES), ds_frame(0))
    assert r["status"] == "pass" and not r.get("vacuous")


def test_te2_breakdown_counts_every_violation_not_only_50(tr: Any) -> None:
    genes = gframe([grow(f"G_{i:05d}") for i in range(60)])
    v = vframe(
        [vrow(f"v{i}", f"G_{i:05d}", 10, "N", "D", "high", "database_derived") for i in range(60)]
    )
    r = tr.te2(v, genes, ds_frame(60))
    assert r["counts"]["violations"] == 60
    assert len(r["violations"]) == 50
    assert sum(r["breakdown"]["alt_type"].values()) == 60
    assert r["breakdown"]["replicon"] == {"chromosome": 60}


def test_te2_n_alt_bins_and_start_codon_miss(tr: Any) -> None:
    genes = gframe([grow("G_1"), grow("G_2"), grow("G_3")])
    v = vframe(
        [
            vrow("a", "G_1", 1, "V", "M", "high", None, n_alt=3),  # aa_pos 1 without start_codon
            vrow("b", "G_2", 40, "N", "D", "high", "database_derived", n_alt=7),
            vrow("c", "G_3", 40, "N", "D", "high", "database_derived", n_alt=1),
        ]
    )
    r = tr.te2(v, genes, ds_frame(3))
    bins = r["breakdown"]["n_alt"]
    assert bins.get("1") == 1
    high = {k: n for k, n in bins.items() if k not in ("1", "2")}
    assert list(high.values()) == [2], bins  # 3 and 7 share the one >= 3 bin
    assert r["breakdown"]["aa_pos_1_unflagged"] == 1


def test_te2_circular_and_linear_replicons_counted_once(tr: Any) -> None:
    """Edge case: a violation on an origin-crossing gene of a circular plasmid and one on a
    linear contig; both counted once, the replicon breakdown is correct."""
    genes = gframe(
        [
            grow("P_1", contig="pMUT2", replicon_type="plasmid", crosses_origin=True),
            grow("L_1", contig="contig_7", replicon_type="unknown"),
        ]
    )
    v = vframe(
        [
            vrow("p", "P_1", 3, "N", "D", "high", "database_derived"),
            vrow("l", "L_1", 9, "Q", "K", "high", "database_derived"),
        ]
    )
    r = tr.te2(v, genes, ds_frame(2))
    assert r["counts"]["violations"] == 2
    assert r["breakdown"]["replicon"] == {"plasmid": 1, "unknown": 1}


# ------------------------------------------------------------------ criterion 2: set R

GAP_FIXTURE = [
    grow("R_01", gene="prmA", product="50S ribosomal protein L11 methyltransferase"),
    grow("R_02", gene="rimL", product="50S ribosomal protein L7/L12-serine acetyltransferase"),
    grow("R_03", gene="rimK", product="30S ribosomal protein S6--L-glutamate ligase"),
    grow("R_04", gene=None, product="30S ribosomal protein S12 RpsL"),
    grow("R_05", gene="rpmE2", product="type B 50S ribosomal protein L31"),
    grow("R_06", gene="rpmH", product="ribosomal protein, bL34 family"),
    grow("R_07", gene="rplA", product="50S ribosomal protein L1", feature_type="pseudogene"),
    grow("R_08", gene=None, product="50S ribosomal protein L7/L12"),
    grow("R_09", gene="rpsD", product="30S ribosomal protein S4"),
    grow("R_10", gene="rpsD", product="30S ribosomal protein S4"),
    grow("R_11", gene="thrA", product="aspartate kinase"),
    grow("R_12", gene="rimO", product="30S ribosomal protein S12 methylthiotransferase RimO"),
]


def test_ribosomal_set_on_the_gap_fixture(tr: Any) -> None:
    got = tag_set(tr.ribosomal_set(gframe(GAP_FIXTURE)))
    # r-proteins (two rpsD genes both in), never the enzymes or the pseudogene
    assert got == {"R_04", "R_05", "R_06", "R_08", "R_09", "R_10"}


def _nissle_genes() -> pd.DataFrame:
    """The committed public Nissle annotation as a genes.tsv-like frame (cds/sorf rows)."""
    doc = json.loads(gzip.decompress(NISSLE_JSON.read_bytes()))
    rows = [
        grow(f["locus"], gene=f.get("gene"), product=f.get("product"), feature_type=f["type"])
        for f in doc["features"]
        if f["type"] in ("cds", "sorf")
    ]
    return gframe(rows)


def _str_or_none(v: object) -> str | None:
    """A Bakta record carries None, not a pandas NaN, for a missing gene or product."""
    return v if isinstance(v, str) else None


def test_ribosomal_set_on_nissle_is_53_and_section_7_1_2_is_59(tr: Any) -> None:
    from types import SimpleNamespace

    from gippa.annotate.bakta_layers import is_ribosomal_protein

    genes = _nissle_genes()
    assert len(tag_set(tr.ribosomal_set(genes))) == 53
    n712 = sum(
        is_ribosomal_protein(SimpleNamespace(gene=_str_or_none(g), product=_str_or_none(p)))  # type: ignore[arg-type]
        for g, p in zip(genes["gene"], genes["product"], strict=True)
    )
    assert n712 == 59


def test_ribosomal_set_reuses_the_section_7_1_2_gene_rule(tr: Any) -> None:
    """``RIBOSOMAL_GENE`` is reused, not copied: patching it changes the set."""
    import gippa.annotate.bakta_layers as bl

    genes = gframe([grow("X_1", gene="rpmH", product="weird")])
    assert tag_set(tr.ribosomal_set(genes)) == {"X_1"}
    text = Path(tr.__file__).read_text(encoding="utf-8")
    assert "RIBOSOMAL_GENE" in text
    assert not re.search(r"re\.compile\(\s*r?[\"']\^\(rps", text), "reuse RIBOSOMAL_GENE, no copy"
    assert bl.RIBOSOMAL_GENE.pattern == r"^(rps|rpl|rpm)[A-Z]"


# ------------------------------------------------------------------ criterion 3: TE3


def _rp(i: int, **kw: Any) -> dict[str, Any]:
    return grow(f"RP_{i:03d}", product=f"50S ribosomal protein L{i}", **kw)


def ex3_genes(rp_values: list[float], other_values: list[float]) -> pd.DataFrame:
    rows = [_rp(i + 1, measured=v) for i, v in enumerate(rp_values)]
    rows.append(_rp(99))  # the tenth r-protein: not detected, not measured
    rows += [grow(f"O_{i:03d}", measured=v) for i, v in enumerate(other_values)]
    return gframe(rows)


def test_te3_example_3_pass(tr: Any) -> None:
    g = ex3_genes([float(v) for v in range(12, 21)], [float(v) for v in range(1, 12)])
    r = tr.te3(g)
    assert r["n_ribosomal"] == 10
    assert_wilson(r["detection"]["share"], 9, 10, (0.596, 0.982))
    assert r["detection"]["holds"] is True  # 0.90 inclusive
    a = r["abundance"]
    assert a["q3"] == pytest.approx(15.25, abs=TOL)
    assert a["m"] == pytest.approx(16.0, abs=TOL)
    assert a["percentile_rank"] == pytest.approx(0.80, abs=TOL)
    assert a["n_r_measured"] == 9 and a["n_measured"] == 20
    assert a["holds"] is True and a["abundance_type"] == "ibaq"
    assert r["status"] == "pass"


def test_te3_example_3_clause_2_fails(tr: Any) -> None:
    others = [1.0, 2.0, 3.0, *[float(v) for v in range(13, 21)]]
    g = ex3_genes([float(v) for v in range(4, 13)], others)
    r = tr.te3(g)
    assert r["abundance"]["m"] == pytest.approx(8.0, abs=TOL)
    assert r["abundance"]["q3"] == pytest.approx(15.25, abs=TOL)
    assert r["abundance"]["percentile_rank"] == pytest.approx(0.40, abs=TOL)
    assert r["status"] == "fail" and r["reason"] == ["abundance"]


def ex4_genes() -> pd.DataFrame:
    """|R| = 52, 46 detected; the 6 missed include rpmE2 (type B)."""
    rows = [grow("RP_B31", gene="rpmE2", product="type B 50S ribosomal protein L31")]
    rows += [_rp(i, measured=20.0) for i in range(1, 47)]  # 46 detected, measured high
    rows += [_rp(i) for i in range(47, 52)]  # 5 more missed
    rows += [grow(f"O_{i:03d}", measured=float(i % 10)) for i in range(100)]
    # two modification enzymes, caught by §7.1.2 only, not detected
    rows += [
        grow("E_1", gene="prmA", product="50S ribosomal protein L11 methyltransferase"),
        grow("E_2", gene="rimK", product="30S ribosomal protein S6--L-glutamate ligase"),
    ]
    return gframe(rows)


def test_te3_example_4_detection_boundary_and_secondaries(tr: Any) -> None:
    r = tr.te3(ex4_genes())
    assert r["n_ribosomal"] == 52
    assert_wilson(r["detection"]["share"], 46, 52, (0.770, 0.946))
    assert r["detection"]["holds"] is False
    assert r["status"] == "fail" and "detection" in r["reason"]
    tb = r["secondary"]["without_type_b"]
    assert tb["n_ribosomal"] == 51
    assert_wilson(tb["detection"]["share"], 46, 51, (0.790, 0.957))
    assert tb["detection"]["holds"] is True  # reported; the verdict is unchanged
    s712 = r["secondary"]["section_7_1_2"]
    assert s712["n_ribosomal"] == 54  # the 52 plus the two enzymes
    assert_wilson(s712["detection"]["share"], 46, 54)


def test_te3_secondaries_never_change_status(tr: Any) -> None:
    rows = [_rp(i, measured=20.0) for i in range(1, 10)] + [_rp(10)]
    rows += [grow(f"O_{i}", measured=float(i)) for i in range(10)]
    rows += [
        grow(f"E_{i}", gene="prmA", product="50S ribosomal protein L11 methyltransferase")
        for i in range(5)
    ]
    r = tr.te3(gframe(rows))
    assert r["secondary"]["section_7_1_2"]["detection"]["holds"] is False  # 9/15
    assert r["status"] == "pass"


def test_te3_boundaries_m_equals_q3_and_even_median(tr: Any) -> None:
    rows = [_rp(1, measured=4.0)] + [grow(f"O_{i}", measured=float(i)) for i in (1, 2, 3, 5)]
    r = tr.te3(gframe(rows))
    assert r["abundance"]["q3"] == pytest.approx(4.0) and r["abundance"]["m"] == pytest.approx(4.0)
    assert r["abundance"]["holds"] is True and r["status"] == "pass"
    rows = [_rp(1, measured=4.0), _rp(2, measured=5.0)]
    rows += [grow(f"O_{i}", measured=float(i)) for i in (1, 2, 3)]
    r = tr.te3(gframe(rows))
    assert r["abundance"]["m"] == pytest.approx(4.5, abs=TOL)  # NumPy mean of the middle two


def test_te3_empty_r_is_not_evaluable(tr: Any) -> None:
    r = tr.te3(gframe([grow("O_1", measured=3.0)]))
    assert r["status"] == "not_evaluable" and "no_ribosomal_genes" in r["reason"]


def test_te3_pseudogene_is_out(tr: Any) -> None:
    rows = [_rp(1, measured=9.0, feature_type="pseudogene"), grow("O_1", measured=3.0)]
    assert tr.te3(gframe(rows))["status"] == "not_evaluable"


def test_te3_abundance_none_or_unmeasured(tr: Any) -> None:
    # clause 1 holds, clause 2 not evaluable -> not_evaluable
    rows = [
        _rp(1, measured=9.0, abundance_type="none"),
        grow("O_1", measured=3.0, abundance_type="none"),
    ]
    assert tr.te3(gframe(rows))["status"] == "not_evaluable"
    rows = [_rp(1, detection_status="detected", detected=True), grow("O_1", measured=3.0)]
    assert tr.te3(gframe(rows))["status"] == "not_evaluable"  # no measured r-protein
    # clause 1 fails, clause 2 not evaluable -> fail
    rows = [_rp(1, detection_status="detected", detected=True), _rp(2), grow("O_1", measured=3.0)]
    r = tr.te3(gframe(rows))
    assert r["status"] == "fail" and "detection" in r["reason"]


def test_te3_imputed_values_are_not_measured(tr: Any) -> None:
    rows = [
        _rp(1, measured=9.0),
        grow("O_1", measured=3.0),
        grow("O_2", abundance_log2=50.0, abundance_imputed=True),
    ]
    a = tr.te3(gframe(rows))["abundance"]
    assert a["n_measured"] == 2
    assert a["q3"] == pytest.approx(float(np.quantile([9.0, 3.0], 0.75)))


@settings(max_examples=60, deadline=None)
@given(
    det=st.lists(st.booleans(), min_size=1, max_size=30),
    vals=st.lists(st.floats(-5, 30, allow_nan=False), min_size=1, max_size=40),
)
def test_te3_property_clauses_match_numpy(det: list[bool], vals: list[float]) -> None:
    """Invariant: clause 1 holds iff k/n >= 0.9; clause 2 iff median(R measured) >= q3(Q)."""
    tr = load_te_real()
    rows = [
        _rp(i + 1, measured=vals[i % len(vals)]) if d else _rp(i + 1) for i, d in enumerate(det)
    ]
    rows += [grow(f"O_{j}", measured=v) for j, v in enumerate(vals)]
    r = tr.te3(gframe(rows))
    k, n = sum(det), len(det)
    assert r["detection"]["holds"] is (k / n >= 0.90)
    if k:
        rv = [vals[i % len(vals)] for i, d in enumerate(det) if d]
        q = rv + list(vals)
        assert r["abundance"]["m"] == pytest.approx(float(np.median(rv)), abs=TOL)
        assert r["abundance"]["q3"] == pytest.approx(float(np.quantile(q, 0.75)), abs=TOL)
        holds = float(np.median(rv)) >= float(np.quantile(q, 0.75))
        expected = "pass" if (k / n >= 0.90 and holds) else "fail"
        assert r["status"] == expected


# ------------------------------------------------------------------ criterion 4: TE4


def te4_genes(**over: Any) -> pd.DataFrame:
    vals = {"rplL": 22.0, "rplJ": 20.1, "atpA": 19.0, "atpD": 19.6}
    top3 = {"rplL": 21.0, "rplJ": 19.5, "atpA": 18.0, "atpD": 18.25}
    rows = []
    for name, v in vals.items():
        kw: dict[str, Any] = {"measured": v, "top3_log2": top3[name], "n_observable_peptides": 9}
        kw.update(over.get(name, {}))
        rows.append(grow(f"T_{name}", gene=name, **kw))
    rows.append(grow("T_other", gene="thrA", measured=10.0))
    return gframe(rows)


def test_te4_example_5_pass(tr: Any) -> None:
    r = tr.te4(te4_genes())
    p1, p2 = r["pairs"]["rplL/rplJ"], r["pairs"]["atpA/atpD"]
    assert p1["log2_ratio"] == 22.0 - 20.1 == 1.8999999999999986  # float64, no rounding
    assert (p1["low"], p1["high"]) == (1, 3) and (p2["low"], p2["high"]) == (-1, 1)
    assert p2["log2_ratio"] == pytest.approx(-0.6, abs=TOL)
    assert p1["status"] == p2["status"] == "pass"
    assert r["status"] == "pass"
    assert r["abundance_type"] == "ibaq"
    assert p1["top3_log2_ratio"] == pytest.approx(1.5, abs=TOL)
    assert p2["top3_log2_ratio"] == pytest.approx(-0.25, abs=TOL)
    sec = r["secondary"]["genes"]
    assert set(sec) == {"rplL", "rplJ", "atpA", "atpD"}
    assert sec["rplL"]["n_observable_peptides"] == 9 and sec["rplL"]["n_peptides_unique"] == 3


@pytest.mark.parametrize(
    ("rplL", "expected"), [(23.0, "pass"), (23.25, "fail"), (21.0, "pass"), (20.99, "fail")]
)
def test_te4_closed_bounds(tr: Any, rplL: float, expected: str) -> None:
    rplJ = 20.0
    r = tr.te4(te4_genes(rplL={"measured": rplL}, rplJ={"measured": rplJ}))
    assert r["pairs"]["rplL/rplJ"]["status"] == expected
    assert r["status"] == expected


@pytest.mark.parametrize(("a", "b"), [(20.0, 21.0), (21.0, 20.0)])
def test_te4_atp_bounds_closed(tr: Any, a: float, b: float) -> None:
    r = tr.te4(te4_genes(atpA={"measured": a}, atpD={"measured": b}))
    assert r["pairs"]["atpA/atpD"]["status"] == "pass"


def test_te4_imputed_gene_makes_the_pair_not_evaluable(tr: Any) -> None:
    r = tr.te4(te4_genes(rplJ={"abundance_imputed": True}))
    assert r["pairs"]["rplL/rplJ"]["status"] == "not_evaluable"
    assert r["pairs"]["atpA/atpD"]["status"] == "pass"
    assert r["status"] == "not_evaluable"


def test_te4_fail_wins_over_not_evaluable(tr: Any) -> None:
    r = tr.te4(te4_genes(rplJ={"abundance_imputed": True}, atpA={"measured": 25.0}))
    assert r["status"] == "fail"


def test_te4_relative_scale(tr: Any) -> None:
    over = {
        n: {"abundance_scale": "relative", "abundance_type": "maxlfq"}
        for n in ("rplL", "rplJ", "atpA", "atpD")
    }
    r = tr.te4(te4_genes(**over))
    for p in r["pairs"].values():
        assert p["status"] == "not_evaluable" and "relative_scale" in p["reason"]
    assert r["status"] == "not_evaluable"


def test_te4_gene_missing_and_ambiguous(tr: Any) -> None:
    g = te4_genes()
    g = g[g["gene"] != "atpD"]
    r = tr.te4(g)
    assert "gene_missing:atpD" in r["pairs"]["atpA/atpD"]["reason"]
    assert r["status"] == "not_evaluable"
    g = te4_genes()
    dup = g[g["gene"] == "rplJ"].assign(locus_tag="T_rplJ_2")
    r = tr.te4(pd.concat([g, dup], ignore_index=True))
    assert r["pairs"]["rplL/rplJ"]["status"] == "not_evaluable"
    assert "gene_ambiguous:rplJ" in r["pairs"]["rplL/rplJ"]["reason"]


def test_te4_lookup_is_case_sensitive_and_cds_sorf_only(tr: Any) -> None:
    g = te4_genes()
    g.loc[g["gene"] == "rplL", "gene"] = "RplL"
    assert "gene_missing:rplL" in tr.te4(g)["pairs"]["rplL/rplJ"]["reason"]
    g = te4_genes()
    pseudo = g[g["gene"] == "atpA"].assign(locus_tag="T_atpA_ps", feature_type="pseudogene")
    r = tr.te4(pd.concat([g, pseudo], ignore_index=True))
    assert r["pairs"]["atpA/atpD"]["status"] == "pass"  # the pseudogene is not in G


@settings(max_examples=80, deadline=None)
@given(a=st.floats(10, 30, allow_nan=False), b=st.floats(10, 30, allow_nan=False))
def test_te4_property_pair_holds_iff_in_closed_interval(a: float, b: float) -> None:
    tr = load_te_real()
    r = tr.te4(te4_genes(rplL={"measured": a}, rplJ={"measured": b}))
    ratio = a - b
    assert r["pairs"]["rplL/rplJ"]["log2_ratio"] == ratio
    assert r["pairs"]["rplL/rplJ"]["status"] == ("pass" if 1 <= ratio <= 3 else "fail")


# ------------------------------------------------------------------ criterion 5: TE5

HAND_HEADERS = [
    "sp|Cont_Q29463|TRY1_BOVIN Anionic trypsin OS=Bos taurus OX=9913",
    "sp|Cont_P00760|TRY1_BOVIN Cationic trypsin OS=Bos taurus OX=9913",
    "sp|Cont_P00974|BPT1_BOVIN Pancreatic trypsin inhibitor OS=Bos taurus OX=9913",
    "sp|Cont_P00766|CTRA_BOVIN Chymotrypsinogen A OS=Bos taurus OX=9913",
    "sp|Cont_Q01726|TRYI_X Trypsin inhibitor OS=Some species",
    "sp|Cont_Q9BYQ7|KRA41_HUMAN Keratin-associated protein 4-1 OS=Homo sapiens OX=9606",
    "sp|Cont_P35527|K1C9_HUMAN Keratin, type I cytoskeletal 9 OS=Homo sapiens OX=9606",
    "sp|Cont_Q00001|TRYP_X Trypsin",
    "sp|Cont_P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913",
]


def test_contaminant_families_mini_and_hand_written(tr: Any) -> None:
    mini = [
        ln[1:].rstrip("\n")
        for ln in MINI_CONTAMINANTS.read_text(encoding="utf-8").splitlines()
        if ln.startswith(">")
    ]
    fam = tr.contaminant_families([*mini, *HAND_HEADERS])
    get = {k: set(v) for k, v in fam.items()}
    assert get.get("trypsin", set()) == {"P00761", "Cont_Q29463", "Cont_P00760", "Cont_Q00001"}
    assert get.get("keratin", set()) == {"P04264", "P13645", "Cont_P35527"}
    assert get.get("keratin_associated", set()) == {"Cont_Q9BYQ7"}
    every = set().union(*get.values())
    for out in ("Cont_P00974", "Cont_P00766", "Cont_Q01726", "Cont_P02769", "P02769"):
        assert out not in every


def test_te5_example_6_fail(tr: Any, cfg: Any) -> None:
    P, hits, unm = ex6_inputs()
    r = tr.te5_arm(P, hits, unm, EX6_FAMILIES, cfg)
    assert r["n_s"] == 4  # albumin is not a family
    assert {k: r["families"][k] for k in ("trypsin", "keratin")} == {"trypsin": 2, "keratin": 2}
    c = r["counts"]
    assert c["genome_exact"] == 1
    assert c["bucket"] == 2 and c["match_other_entry"] == 1
    assert c["placed_tolerant"] == 1
    assert c["labelled_unmapped"] == 0 and c["absent"] == 0
    assert_wilson(r["share_bucket"], 2, 3)
    assert r["status"] == "fail"


def test_te5_example_6_without_the_tolerant_peptide(tr: Any, cfg: Any) -> None:
    P, hits, unm = ex6_inputs(with_tolerant=False)
    r = tr.te5_arm(P, hits, unm, EX6_FAMILIES, cfg)
    assert r["n_s"] == 3
    assert_wilson(r["share_bucket"], 2, 2, (0.342, 1.0))
    assert r["status"] == "pass"


def test_te5_peptides_outside_select_queries_are_not_in_s(tr: Any, cfg: Any) -> None:
    P, hits, unm = ex6_inputs(with_tolerant=False)
    long41 = "A" * 40 + "K"
    extra = pframe(
        [
            prow(long41, "Cont_P00761"),
            prow("DECOYPEPK", "Cont_P00761", decoy=True),
            prow("HIGHQPEPK", "Cont_P13645", q=0.02),
        ]
    )
    P2 = pd.concat([P, extra], ignore_index=True)
    r = tr.te5_arm(P2, hits, unm, EX6_FAMILIES, cfg)
    assert r["n_s"] == 3 and r["status"] == "pass"  # none of the three counts as absent


def test_te5_decoy_row_with_a_family_id_does_not_put_a_target_peptide_in_s(
    tr: Any, cfg: Any
) -> None:
    """Surviving mutant 9: S takes family IDs from accepted *target* rows only. The peptide is in
    **U** through its target row (a genome protein); its accepted decoy row names trypsin."""
    from gippa.steps.mapping import select_queries

    P, hits, unm = ex6_inputs(with_tolerant=False)
    mixed = "MIXEDPEPTK"
    extra = pframe([prow(mixed, "JMJJOG_00099"), prow(mixed, "Cont_P00761", decoy=True)])
    P2 = pd.concat([P, extra], ignore_index=True)
    assert mixed in set(select_queries(P2, cfg)["peptide"])  # the fixture's premise
    r = tr.te5_arm(P2, hits, unm, EX6_FAMILIES, cfg)
    assert r["n_s"] == 3 and r["counts"]["absent"] == 0
    assert r["status"] == "pass"


@pytest.mark.parametrize("match", [None, ""], ids=["na", "empty"])
def test_te5_bucket_row_without_a_match_is_not_other_entry(
    tr: Any, cfg: Any, match: str | None
) -> None:
    """Surviving mutant 8: an empty ``contaminant_match`` names no other entry."""
    P = pframe([prow("VATVSLPR", "Cont_P00761")])
    unm = unmapped_frame([("VATVSLPR", "contaminant", match)])
    r = tr.te5_arm(P, hits_frame([]), unm, EX6_FAMILIES, cfg)
    assert r["counts"]["bucket"] == 1
    assert r["counts"]["match_other_entry"] == 0
    assert r["status"] == "pass"


def test_te5_labelled_unmapped_and_absent(tr: Any, cfg: Any) -> None:
    P = pframe(
        [
            prow("VATVSLPR", "Cont_P00761"),
            prow("NLDLDSLLK", "Cont_P04264"),
            prow("GSYGSGGSSYGSGGGK", "Cont_P13645"),
        ]
    )
    unm = unmapped_frame(
        [("VATVSLPR", "contaminant", "Cont_P00761"), ("NLDLDSLLK", "unmapped", None)]
    )
    r = tr.te5_arm(P, hits_frame([]), unm, EX6_FAMILIES, cfg)
    assert r["counts"]["labelled_unmapped"] == 1 and r["counts"]["absent"] == 1
    assert r["counts"]["bucket"] == 1 and r["status"] == "fail"


def test_te5_only_genome_exact_is_not_evaluable(tr: Any, cfg: Any) -> None:
    P = pframe([prow("LSSPATLNSR", "Cont_P00761;JMJJOG_01234")])
    r = tr.te5_arm(P, hits_frame([("LSSPATLNSR", "exact")]), unmapped_frame([]), EX6_FAMILIES, cfg)
    assert r["status"] == "not_evaluable"
    assert "no_trypsin_keratin_peptides" in r["reason"]
    assert r["counts"]["genome_exact"] == 1


def test_te5_peptide_in_both_families_counts_once_in_s(tr: Any, cfg: Any) -> None:
    P = pframe([prow("SHAREDPEPK", "Cont_P00761;Cont_P13645")])
    unm = unmapped_frame([("SHAREDPEPK", "contaminant", "Cont_P00761")])
    r = tr.te5_arm(P, hits_frame([]), unm, EX6_FAMILIES, cfg)
    assert r["n_s"] == 1
    assert {k: r["families"][k] for k in ("trypsin", "keratin")} == {"trypsin": 1, "keratin": 1}
    assert r["counts"]["bucket"] == 1 and r["status"] == "pass"


def test_te5_keratin_associated_is_secondary_only(tr: Any, cfg: Any) -> None:
    P = pframe([prow("VATVSLPR", "Cont_P00761"), prow("KAPPEPTIDEK", "Cont_Q9BYQ7")])
    unm = unmapped_frame([("VATVSLPR", "contaminant", "Cont_P00761")])
    hits = hits_frame([("KAPPEPTIDEK", "tolerant")])
    r = tr.te5_arm(P, hits, unm, EX6_FAMILIES, cfg)
    assert r["n_s"] == 1 and r["status"] == "pass"  # the tolerant KAP peptide is not judged
    sec = r["secondary"]["keratin_associated"]
    assert sec["n_s"] == 1 and sec["counts"]["placed_tolerant"] == 1


def test_te5_exact_beats_bucket_first_match(tr: Any, cfg: Any) -> None:
    P = pframe([prow("VATVSLPR", "Cont_P00761")])
    unm = unmapped_frame([("VATVSLPR", "contaminant", "Cont_P00761")])
    r = tr.te5_arm(
        P, hits_frame([("VATVSLPR", "exact"), ("VATVSLPR", "tolerant")]), unm, EX6_FAMILIES, cfg
    )
    assert r["counts"]["genome_exact"] == 1 and r["counts"]["bucket"] == 0


def test_te5_json_holds_no_peptide(tr: Any, cfg: Any) -> None:
    """te_real.json is aggregate-only: the arm result names no sequence. A private key (leading
    ``_``, review trim T1) may carry the per-peptide classes that ``compute`` pops before writing;
    the written file is checked in test_m10_06_cli.py."""
    P, hits, unm = ex6_inputs()
    r = tr.te5_arm(P, hits, unm, EX6_FAMILIES, cfg)
    text = json.dumps({k: v for k, v in r.items() if not str(k).startswith("_")}, default=str)
    for pep in ("VATVSLPR", "SLNNQFASFIDK", "LSSPATLNSR", "LASYLDK", ALBUMIN_PEPTIDE):
        assert pep not in text
    assert not re.search(r"[ACDEFGHIKLMNPQRSTVWY]{7,}", text)
