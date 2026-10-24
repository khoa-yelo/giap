"""Pure functions of ``validation/pilot/metrics.py`` (M10-05 acceptance criteria 1-4).

Spec: docs/specs/M10-05-pilot-metrics.md ("Metric definitions", "TE1 on real data", "Hypotheses
1 and 3", examples 1-14, edge cases). Every CI is compared with ``giap.stats.cutoffs.wilson_ci``
to 1e-12 and with the digits printed in the spec; the H1 p-values with ``scipy.stats.binomtest``
to 1e-12.

Interface fixed here (test-author reading of the spec "Interface"; the spec names the functions
and their arguments, these tests fix the returned keys):

* a **share** is a dict ``{"value", "ci_low", "ci_high", ...}``; with a zero denominator
  ``value`` is None and ``reason = "empty_denominator"``;
* a metric-1 **count** is ``{"all", "contaminant_only", "non_contaminant"}``;
* ``identifications(P, prefixes)`` -> ``{"union": {precursors, peptides, protein_groups,
  protein_groups_min2}, "per_run": {run: {...the same four...}}}``;
* ``mapping_shares(P, peptide_hits, unmapped, cfg)`` -> ``n_universe``, ``n_contaminant_only``,
  ``n_denominator``, ``counts`` {exact, tolerant, unmapped}, ``shares`` (same keys, shares),
  ``unmapped_reasons`` {no_hit, rejected_placement_fdr, below_identity, too_short, not_listed},
  ``contaminant_difference`` (|C| minus the distinct ``label = contaminant, kind = peptide``
  rows of ``unmapped.tsv``);
* ``genome_inconsistent(detection_summary)`` -> ``{"value", "n", "count"}``;
* ``strain_specific_set(nissle_proteins, k12_proteins)``: both ``{id: protein sequence}``;
  returns ``{locus_tag: "strain_specific" | "k12_shared" | "no_peptides"}``;
* ``coverage(genes, strain_specific)`` -> per subset (a dict keyed by subset name, or a list of
  entries with a ``subset`` key, in table order) ``{"genes", "detected", "share", "by_status"}``;
  per plasmid contig the subset ``plasmid:<contig>``; pseudogenes under ``pseudogene``;
* ``abundance_agreement`` / ``stability``: the compare rows copied (``stability`` keyed by
  family, ``novel_orfs`` and ``alt_sites`` under ``other_families``);
* ``runtime(manifest, trace_rows, sacct_rows)`` -> ``wall_time_s``, ``sum_task_wall_s``,
  ``max_step_peak_rss_bytes``, ``cached_steps``, ``runtime_complete``, ``steps``, ``hosts``;
* ``te1_real(P, peptide_hits, unmapped, genes, circular, non_atg, cfg, q_crossings=...)`` ->
  ``(te1, exceptions)``: ``te1`` holds ``n_peptides`` (|D2|), ``exact_inframe``,
  ``exact_inframe_masked``, ``share_exact_inframe``, ``ci_low``, ``ci_high``, ``passed``,
  ``ci_low_below_threshold``, ``target_q_na``, ``non_atg_nterm`` {peptides, misses},
  ``misses`` {non_atg_nterm, exact_other_gene, exact_other_class: {class: n}, tolerant_only,
  unplaced: {reason: n}} (first match; the counts sum to ``n_misses``), ``unknown_source``
  (review S1: a separate count, not a miss cause: peptides of D2 none of whose sources is a
  locus tag of ``genes``), ``share_excluding_non_atg`` (a share) and
  ``share_without_q_at_cutoff`` (a share plus ``removed``, ``q_crossings``,
  ``pass_becomes_fail``; ``value`` None with ``reason = "no_crossings"`` when q_crossings = 0;
  an empty superset denominator keeps ``pass_becomes_fail`` True with a ``reason``);
  ``exceptions`` is the uncapped miss list (rows ``kind``, ``locus_tag``, ``peptide``, ``hits``),
  kept out of ``te1`` (no sequence enters ``metrics.json``). ``q_crossings`` is a keyword: the
  spec's signature does not carry it, but example 12 needs it;
* ``hypothesis1(genes_a1, genes_a2, compare_detection, subsets)``; ``subsets`` maps
  ``plasmid``, ``prophage``, ``bgc``, ``strain_specific`` to locus-tag sets (the genes frames
  carry the same membership in ``replicon_type`` / ``region_ids``). Returns ``b``, ``c``,
  ``prop_a2_only``, ``coverage_diff``, ``b_tu``, ``c_tu``, ``mixed_tu``, ``p_tu``, ``p_gene``,
  ``unit_effect``, ``mostly``, ``mostly_base_rate``, ``strain_specific_share``, ``verdict``,
  ``failed_clauses``, ``warnings``;
* ``hypothesis3(m4, m5)`` -> ``verdict`` and ``criteria`` {rho, islands, operon_flags,
  enriched_categories}: a status string or a dict with ``status``;
* ``render(template, metrics, notes) -> str`` (an unresolved path renders ``**missing: <path>**``;
  a rejected notes text raises ``Inconsistent``, review B1: exit 1 means only that);
* review S2: ``mapping_shares``, ``te1_real`` and ``te1_diagnostic`` take the raw
  ``read_peptide_table`` frame and build **U** as GIAP's mapping step does,
  ``select_queries(validate_peptides(P, cfg.input)[0], cfg)``; ``identifications`` (metric 1)
  keeps the raw rows (its cross-check is against the converter's count on the same file);
* review S6: a missing compare row is ``reason = "missing_row"`` with ``flags`` None (no
  synthesized flag);
* review S9: shares (a share dict, and the scalar keys ``share*`` / ``ci_*``) render with 4
  significant figures and never as a bare "1" below 1; a dict without ``value`` renders as
  ``k: v`` pairs (``none`` when empty).
"""

from __future__ import annotations

import json
import math
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.compare_data import DETECTION_COLUMNS, gene, genes_frame
from helpers.pilot_metrics_m10 import GIB, load_metrics, load_te_proxies
from scipy.stats import binomtest

from giap.config import resolve_config
from giap.model.coords import Interval
from giap.model.genes import Gene
from giap.stats.cutoffs import wilson_ci
from giap.steps.mapping import select_queries

PREFIXES = ("Cont_", "CON__", "CONT_", "cRAP")
TOL = 1e-12


@pytest.fixture(scope="module")
def m() -> Any:
    return load_metrics()


@pytest.fixture(scope="module")
def cfg() -> Any:
    return resolve_config()


# ------------------------------------------------------------------ helpers


def _share(s: Mapping[str, Any]) -> tuple[Any, Any, Any]:
    return s["value"], s["ci_low"], s["ci_high"]


def assert_wilson(s: Mapping[str, Any], k: int, n: int, printed: tuple[float, ...] = ()) -> None:
    """``s`` is k/n with ``wilson_ci(k, n)`` to 1e-12, and matches the spec's printed digits."""
    value, lo, hi = _share(s)
    wlo, whi = wilson_ci(k, n)
    assert value == pytest.approx(k / n, abs=TOL)
    assert lo == pytest.approx(wlo, abs=TOL) and hi == pytest.approx(whi, abs=TOL)
    if printed:
        assert round(lo, 3) == pytest.approx(printed[0]), (lo, printed)  # the spec prints 3
        assert round(hi, 3) == pytest.approx(printed[1]), (hi, printed)


def _count(x: Any, cls: str = "all") -> Any:
    return x[cls] if isinstance(x, Mapping) else x


def _subset(cov: Any, name: str) -> Mapping[str, Any]:
    if isinstance(cov, Mapping) and "subsets" in cov:
        cov = cov["subsets"]
    if isinstance(cov, Mapping):
        return cov[name]  # type: ignore[no-any-return]
    for entry in cov:
        if entry.get("subset") == name:
            return entry  # type: ignore[no-any-return]
    raise KeyError(name)


def _status(c: Any) -> str:
    return str(c["status"] if isinstance(c, Mapping) else c)


def _peptides(rows: Iterable[Sequence[Any]]) -> pd.DataFrame:
    """Rows (run, peptide, modified, charge, q, is_decoy, protein_ids) as ``read_peptide_table``
    types them."""
    cols = ("run", "peptide", "modified_peptide", "charge", "q_value", "is_decoy", "protein_ids")
    df = pd.DataFrame([dict(zip(cols, r, strict=True)) for r in rows], columns=list(cols))
    return df.astype({"run": "string", "peptide": "string", "modified_peptide": "string",
                      "charge": "Int64", "q_value": float, "is_decoy": "boolean",
                      "protein_ids": "string"})  # fmt: skip


def _p(pep: str, ids: str, q: float | None = 0.001, *, decoy: bool = False, run: str = "WT_1",
       charge: int = 2) -> tuple[Any, ...]:  # fmt: skip
    return (run, pep, pep, charge, np.nan if q is None else q, decoy, ids)


HIT_COLS = ("peptide", "route", "locus_class", "gene", "contig", "start", "end", "strand")


def _hits(rows: Iterable[Sequence[Any]]) -> pd.DataFrame:
    """``peptide_hits.tsv`` rows (1-based closed ``start``/``end``, as written)."""
    return pd.DataFrame([dict(zip(HIT_COLS, r, strict=True)) for r in rows], columns=list(HIT_COLS))


def _unmapped(rows: Iterable[tuple[str, str, str]]) -> pd.DataFrame:
    """(sequence, label, reason) rows of ``unmapped.tsv``."""
    cols = ["sequence_id", "sequence", "kind", "label", "reason", "contaminant_match"]
    return pd.DataFrame([{"sequence_id": s, "sequence": s, "kind": "peptide", "label": lab,
                          "reason": r, "contaminant_match": None} for s, lab, r in rows],
                        columns=cols)  # fmt: skip


_ALPHA = "ACDEFGHMNQSTVWY"


def _pep(i: int, prefix: str = "PEP") -> str:
    """A distinct 8-mer (``prefix`` + 4 letters + K) per ``i``."""
    s = ""
    for _ in range(4):
        i, r = divmod(i, len(_ALPHA))
        s += _ALPHA[r]
    return prefix + s + "K"


def _gene(tag: str, contig: str, start: int, end: int, strand: str = "+") -> Gene:
    return Gene(tag, "cds", Interval(contig, start, end, strand), "M" + "A" * 20, None, "ATG")  # type: ignore[arg-type]


# ================================================================== metric 1


def test_example_1_identifications_union_and_per_run(m: Any) -> None:
    P = _peptides([
        _p("PEPTIDEK", "JMJJOG_00010", 0.001, run="WT_1", charge=2),
        _p("PEPTIDEK", "JMJJOG_00010", 0.004, run="WT_1", charge=3),
        _p("PEPTIDEK", "JMJJOG_00010", 0.002, run="WT_2", charge=2),
        _p("AAAAGGGK", "JMJJOG_00020", 0.02, run="WT_2"),
        _p("KEDITPEP", "rev_JMJJOG_00010", 0.005, run="WT_2", decoy=True),
        _p("LLLGGGR", "Cont_P00761", 0.01, run="WT_2"),  # q = 0.01 is accepted
    ])  # fmt: skip
    out = m.identifications(P, PREFIXES)
    u = out["union"]
    assert (u["precursors"]["all"], u["precursors"]["contaminant_only"],
            u["precursors"]["non_contaminant"]) == (3, 1, 2)  # fmt: skip
    assert (u["peptides"]["all"], u["peptides"]["contaminant_only"],
            u["peptides"]["non_contaminant"]) == (2, 1, 1)  # fmt: skip
    assert (u["protein_groups"]["all"], u["protein_groups"]["contaminant_only"],
            u["protein_groups"]["non_contaminant"]) == (2, 1, 1)  # fmt: skip
    assert _count(u["protein_groups_min2"]) == 0
    runs = out["per_run"]
    assert _count(runs["WT_1"]["precursors"]) == 2 and _count(runs["WT_1"]["peptides"]) == 1
    assert _count(runs["WT_2"]["precursors"]) == 2 and _count(runs["WT_2"]["peptides"]) == 2


def test_example_2_groups_sorted_and_multi_row_contaminant_rule(m: Any) -> None:
    P = _peptides([
        _p("AAAWWWK", "JMJJOG_00031;JMJJOG_00030"),
        _p("CCCWWWK", "JMJJOG_00030;JMJJOG_00031"),
        _p("SAMPLEPEPK", "Cont_P00722", run="WT_1"),
        _p("SAMPLEPEPK", "JMJJOG_00350", run="WT_2"),
    ])  # fmt: skip
    u = m.identifications(P, PREFIXES)["union"]
    # groups: {00030, 00031} (B;A == A;B), {Cont_P00722}, {00350}
    assert u["protein_groups"]["all"] == 3
    assert u["protein_groups"]["contaminant_only"] == 1  # Cont_P00722 classed by its own items
    assert _count(u["protein_groups_min2"]) == 1  # {00030, 00031}: 2 distinct peptides
    # SAMPLEPEPK has a strain item in one accepted row: not contaminant-only
    assert u["peptides"]["contaminant_only"] == 0 and u["peptides"]["non_contaminant"] == 3
    assert u["precursors"]["contaminant_only"] == 0


def test_metric_1_na_q_target_is_not_accepted(m: Any) -> None:
    P = _peptides([_p("PEPTIDEK", "G1", 0.001), _p("QQQNNNK", "G2", None)])
    assert m.identifications(P, PREFIXES)["union"]["peptides"]["all"] == 1


# ================================================================== metric 2


def _example_4() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    peps = [_pep(i) for i in range(10)]
    P = _peptides([_p(p, f"P{i:05d}") for i, p in enumerate(peps)]
                  + [_p("WWWGGGWWK", "Cont_P02769")])  # fmt: skip
    hits = _hits([(p, "exact" if i < 7 else "tolerant", "cds_inframe", "G1", "chr", 1, 24, "+")
                  for i, p in enumerate(peps[:9])])  # fmt: skip
    unm = _unmapped([(peps[9], "unmapped", "no_hit"),
                     ("WWWGGGWWK", "contaminant", "contaminant_id")])  # fmt: skip
    return P, hits, unm


def test_example_4_mapping_shares(m: Any, cfg: Any) -> None:
    P, hits, unm = _example_4()
    out = m.mapping_shares(P, hits, unm, cfg)
    assert (out["n_universe"], out["n_contaminant_only"], out["n_denominator"]) == (11, 1, 10)
    assert out["counts"] == {"exact": 7, "tolerant": 2, "unmapped": 1}
    assert_wilson(out["shares"]["exact"], 7, 10, (0.397, 0.892))
    assert_wilson(out["shares"]["tolerant"], 2, 10)
    assert_wilson(out["shares"]["unmapped"], 1, 10)
    assert sum(out["shares"][k]["value"] for k in ("exact", "tolerant", "unmapped")) == (
        pytest.approx(1.0, abs=TOL)
    )
    assert out["unmapped_reasons"]["no_hit"] == 1 and out["unmapped_reasons"]["not_listed"] == 0
    assert out["contaminant_difference"] == 0


def test_example_4_genome_inconsistent_row(m: Any, cfg: Any) -> None:
    ds = pd.DataFrame([
        {"section": "qc", "group_type": "run", "group": "all", "metric": "id_mixups",
         "value": 0.0, "ci_low": np.nan, "ci_high": np.nan, "n": 4, "flags": None},
        {"section": "qc", "group_type": "run", "group": "all",
         "metric": "genome_inconsistent_rate", "value": 0.3, "ci_low": np.nan, "ci_high": np.nan,
         "n": 10, "flags": None},
    ])  # fmt: skip
    gi = m.genome_inconsistent(ds)
    assert gi["value"] == pytest.approx(0.3) and gi["n"] == 10 and gi["count"] == 3
    P, hits, unm = _example_4()
    ms = m.mapping_shares(P, hits, unm, cfg)
    assert gi["count"] == ms["counts"]["tolerant"] + ms["counts"]["unmapped"]


def test_genome_inconsistent_na_row_for_a2(m: Any) -> None:
    ds = pd.DataFrame([{"section": "qc", "group_type": "run", "group": "all",
                        "metric": "genome_inconsistent_rate", "value": np.nan, "ci_low": np.nan,
                        "ci_high": np.nan, "n": pd.NA, "flags": None}])  # fmt: skip
    gi = m.genome_inconsistent(ds)
    assert gi["value"] is None and gi["count"] is None


def test_criterion_2_universe_is_select_queries(m: Any, cfg: Any) -> None:
    long41 = "A" * 40 + "K"
    P = _peptides([
        _p("PEPTIDEK", "G1"),
        _p("KEDITPEPK", "rev_G1", 0.001, decoy=True),  # decoy: out
        _p("AAAAGGGK", "G1", 0.02),  # q above 0.01: out
        _p(long41, "G1"),  # 41 aa: out (mapping.exact.max_len 40)
        _p("QQQNNNK", "G1", None),  # target with NA q: kept, as select_queries keeps it
        _p("LLLGGGR", "G1", 0.01),  # q = 0.01: kept
    ])  # fmt: skip
    hits = _hits([("PEPTIDEK", "exact", "cds_inframe", "G1", "chr", 1, 24, "+"),
                  ("QQQNNNK", "tolerant", "cds_inframe", "G1", "chr", 1, 21, "+")])  # fmt: skip
    unm = _unmapped([("LLLGGGR", "unmapped", "no_hit")])
    out = m.mapping_shares(P, hits, unm, cfg)
    want = set(select_queries(P, cfg)["peptide"].astype(str))
    assert want == {"PEPTIDEK", "QQQNNNK", "LLLGGGR"}
    assert out["n_universe"] == len(want) == out["n_denominator"]
    assert out["counts"] == {"exact": 1, "tolerant": 1, "unmapped": 1}
    total = sum(out["shares"][k]["value"] for k in ("exact", "tolerant", "unmapped"))
    assert total == pytest.approx(1.0, abs=TOL)


def test_unmapped_reason_first_in_order_and_not_listed(m: Any, cfg: Any) -> None:
    peps = [_pep(i) for i in range(3)]
    P = _peptides([_p(p, "G1") for p in peps])
    unm = _unmapped([(peps[0], "unmapped", "below_identity"), (peps[0], "unmapped", "no_hit"),
                     (peps[1], "unmapped", "too_short")])  # fmt: skip
    out = m.mapping_shares(P, _hits([]), unm, cfg)
    r = out["unmapped_reasons"]
    assert r["no_hit"] == 1  # no_hit comes before below_identity
    assert r["below_identity"] == 0 and r["too_short"] == 1 and r["not_listed"] == 1


def test_metric_2_empty_denominator_is_null(m: Any, cfg: Any) -> None:
    P = _peptides([_p("WWWGGGWWK", "Cont_P02769")])
    out = m.mapping_shares(P, _hits([]), _unmapped([]), cfg)
    assert out["n_denominator"] == 0
    for k in ("exact", "tolerant", "unmapped"):
        assert out["shares"][k]["value"] is None
        assert out["shares"][k]["reason"] == "empty_denominator"


# ================================================================== metric 3


def _coverage_genes(a2: bool) -> pd.DataFrame:
    """Example 5: 88 chromosome, 10 plasmid (6 pMUT1, 4 pMUT2), 2 unknown; one pseudogene."""
    rows: list[dict[str, object]] = []
    det_chr = 52 if a2 else 48
    det_pl = {"pMUT1": 5, "pMUT2": 3} if a2 else {"pMUT1": 1, "pMUT2": 1}
    for i in range(88):
        reg = "prophage_1" if i < 5 else ("bgc_1;prophage_2" if i < 8 else None)
        status = "detected" if i < det_chr else (
            "single_peptide_unestimated" if i < det_chr + 2 else "not_detected")  # fmt: skip
        rows.append(gene(f"C{i:03d}", start=i * 400, detected=status == "detected",
                         detection_status=status, replicon_type="chromosome",
                         region_ids=reg))  # fmt: skip
    for contig, n in (("pMUT1", 6), ("pMUT2", 4)):
        for i in range(n):
            d = i < det_pl[contig]
            rows.append(gene(f"{contig}_{i}", start=i * 400, contig=contig, detected=d,
                             replicon_type="plasmid"))  # fmt: skip
    for i in range(2):
        rows.append(gene(f"U{i}", start=i * 400, contig="ctg_x", detected=False,
                         replicon_type="unknown"))  # fmt: skip
    rows.append(gene("PS1", start=90000, detected=False, feature_type="pseudogene",
                     replicon_type="chromosome"))  # fmt: skip
    return genes_frame(rows)


def _classes() -> dict[str, str]:
    tags = [f"C{i:03d}" for i in range(88)] + [f"pMUT1_{i}" for i in range(6)]
    tags += [f"pMUT2_{i}" for i in range(4)] + ["U0", "U1"]
    out = {t: "k12_shared" for t in tags}
    for t in tags[:5]:
        out[t] = "no_peptides"
    for t in tags[88:98]:
        out[t] = "strain_specific"
    return out


def test_example_5_coverage(m: Any) -> None:
    a2 = m.coverage(_coverage_genes(True), _classes())
    a1 = m.coverage(_coverage_genes(False), _classes())
    all2, all1 = _subset(a2, "all"), _subset(a1, "all")
    assert (all2["genes"], all2["detected"]) == (100, 60)
    assert_wilson(all2["share"], 60, 100, (0.502, 0.691))
    assert_wilson(all1["share"], 50, 100, (0.404, 0.596))
    assert_wilson(_subset(a2, "plasmid")["share"], 8, 10, (0.490, 0.943))
    assert_wilson(_subset(a1, "plasmid")["share"], 2, 10, (0.057, 0.510))
    reps = [_subset(a2, k)["genes"] for k in ("chromosome", "plasmid", "unknown")]
    assert reps == [88, 10, 2] and sum(reps) == 100
    assert _subset(a2, "plasmid:pMUT1")["genes"] == 6
    assert _subset(a2, "plasmid:pMUT1")["detected"] == 5
    assert _subset(a2, "plasmid:pMUT2")["genes"] == 4
    # every status is inside the denominator, counted per status
    st = all2["by_status"]
    assert st["detected"] == 60 and st["single_peptide_unestimated"] == 2
    assert sum(st.values()) == 100


def test_coverage_regions_classes_and_pseudogene_row(m: Any) -> None:
    a2 = m.coverage(_coverage_genes(True), _classes())
    assert _subset(a2, "prophage")["genes"] == 8  # prophage_1 x5 and bgc_1;prophage_2 x3
    assert _subset(a2, "bgc")["genes"] == 3
    assert _subset(a2, "strain_specific")["genes"] == 10
    assert _subset(a2, "k12_shared")["genes"] == 85
    assert _subset(a2, "no_peptides")["genes"] == 5
    assert _subset(a2, "pseudogene")["genes"] == 1
    assert _subset(a2, "all")["genes"] == 100  # the pseudogene is not in G


def test_example_6_strain_specific_set(m: Any) -> None:
    nissle = {"X": "AAAGGGLKPEPTLDEK", "Y": "QQQWWWR", "Z": "MKAAK"}
    out = m.strain_specific_set(nissle, {"sp|P1|K12": "PEPTIDEK"})
    assert out == {"X": "k12_shared", "Y": "strain_specific", "Z": "no_peptides"}


# ================================================================== metrics 4, 5 and H3


def _abundance(**over: Any) -> pd.DataFrame:
    row: dict[str, Any] = {
        "run_a": "A1", "run_b": "A2", "n_shared": 45, "n_units": 40, "spearman_rho": 0.83,
        "ci_low": 0.76, "ci_high": 0.88, "ba_bias_log2": 0.1, "ba_loa_low": -0.9,
        "ba_loa_high": 1.1, "method": "spearman_tu_bootstrap", "n_boot": 1000, "flags": None,
    }  # fmt: skip
    row.update(over)
    other = {**row, "run_a": "A1", "run_b": "A3", "spearman_rho": 0.1, "ci_low": 0.0}
    return pd.DataFrame([other, row])


def _stab(rows: Sequence[tuple[str, Any, Any, Any, Any, Any]]) -> pd.DataFrame:
    cols = ("family", "n_a", "n_b", "n_shared", "jaccard", "flags")
    df = pd.DataFrame([dict(zip(cols, r, strict=True)) for r in rows])
    df.insert(1, "run_a", "A1")
    df.insert(2, "run_b", "A2")
    return df


EX10 = [
    ("islands", 10, 8, 6, 0.5, None),
    ("operon_flags", 0, 0, 0, None, "both_empty"),
    ("enriched_categories", None, 3, None, None, "missing_in:A1"),
    ("novel_orfs", 2, 4, 1, 0.2, None),
    ("alt_sites", 0, 0, 0, None, "both_empty"),
]


def test_example_10_metric_4_copied(m: Any) -> None:
    out = m.abundance_agreement(_abundance(), ("A1", "A2"))
    want = {"n_shared": 45, "n_units": 40, "spearman_rho": 0.83, "ci_low": 0.76,
            "ci_high": 0.88, "method": "spearman_tu_bootstrap", "n_boot": 1000,
            "ba_bias_log2": 0.1, "ba_loa_low": -0.9, "ba_loa_high": 1.1}  # fmt: skip
    for k, v in want.items():
        assert out[k] == (pytest.approx(v) if isinstance(v, float) else v), k


def test_metric_4_na_rho_copied_with_flags(m: Any) -> None:
    out = m.abundance_agreement(
        _abundance(spearman_rho=np.nan, ci_low=np.nan, ci_high=np.nan, flags="too_few_shared"),
        ("A1", "A2"),
    )
    assert out["spearman_rho"] is None and "too_few_shared" in str(out["flags"])


def test_example_10_metric_5_copied(m: Any) -> None:
    out = m.stability(_stab(EX10), ("A1", "A2"))
    assert out["islands"]["jaccard"] == pytest.approx(0.5)
    assert (out["islands"]["n_a"], out["islands"]["n_b"], out["islands"]["n_shared"]) == (10, 8, 6)
    assert out["operon_flags"]["jaccard"] is None
    assert "both_empty" in str(out["operon_flags"]["flags"])
    assert "missing_in:A1" in str(out["enriched_categories"]["flags"])
    assert set(out["other_families"]) == {"novel_orfs", "alt_sites"}
    assert out["other_families"]["novel_orfs"]["jaccard"] == pytest.approx(0.2)


def _h3(m: Any, abundance: pd.DataFrame, stab: Sequence[tuple[Any, ...]]) -> dict[str, Any]:
    m4 = m.abundance_agreement(abundance, ("A1", "A2"))
    m5 = m.stability(_stab(stab), ("A1", "A2"))
    out: dict[str, Any] = m.hypothesis3(m4, m5)
    return out


EX9 = [
    ("islands", 10, 8, 6, 0.5, None),
    ("operon_flags", 3, 1, 1, 1 / 3, None),
    ("enriched_categories", None, 3, None, None, "missing_in:A1"),
]


def test_example_9_h3_supported(m: Any) -> None:
    out = _h3(m, _abundance(), EX9)
    c = out["criteria"]
    assert _status(c["rho"]) == "holds"  # ci_low 0.76 >= 0.7
    assert _status(c["islands"]) == "holds"  # union 12, J 0.5
    assert _status(c["operon_flags"]) == "not_evaluable"  # union 3 < 5
    assert _status(c["enriched_categories"]) == "not_evaluable"  # missing_in
    assert out["verdict"] == "supported"


def test_example_9_variant_h3_not_supported(m: Any) -> None:
    stab = [("islands", 10, 8, 5, 5 / 13, None), *EX9[1:]]
    out = _h3(m, _abundance(), stab)
    assert _status(out["criteria"]["islands"]) == "fails"  # union 13, J 0.385
    assert out["verdict"] == "not_supported"


def test_h3_rho_fails_and_rho_not_evaluable(m: Any) -> None:
    out = _h3(m, _abundance(ci_low=0.65), EX9)
    assert _status(out["criteria"]["rho"]) == "fails" and out["verdict"] == "not_supported"
    out = _h3(m, _abundance(flags="too_few_shared"), EX9)
    assert _status(out["criteria"]["rho"]) == "not_evaluable"
    assert out["verdict"] == "not_evaluable"


def test_h3_every_family_not_evaluable(m: Any) -> None:
    stab = [("islands", 0, 0, 0, None, "both_empty"), ("operon_flags", 3, 1, 1, 1 / 3, None),
            ("enriched_categories", None, 3, None, None, "missing_in:A1")]  # fmt: skip
    out = _h3(m, _abundance(), stab)
    assert _status(out["criteria"]["rho"]) == "holds"
    assert out["verdict"] == "not_evaluable"


# ================================================================== metric 7


def test_example_11_runtime(m: Any) -> None:
    manifest = {
        "started_at": "2026-10-22T10:00:00+00:00",
        "finished_at": "2026-10-22T10:42:30+00:00",
        "exit_code": 0,
        "steps": [
            {"name": "mapping", "status": "completed", "n_tasks": 1, "n_retries": 0,
             "wall_time_s": 120.0, "cpu_time_s": 100.0, "peak_rss_bytes": 2 * GIB},
            {"name": "bakta", "status": "cached", "n_tasks": 1, "n_retries": 0},
            {"name": "gene_table", "status": "completed", "n_tasks": 1, "n_retries": 0,
             "wall_time_s": 30.0, "cpu_time_s": 25.0},
        ],
    }  # fmt: skip
    trace = [{"name": "A", "hostname": "n2"}, {"name": "B", "hostname": "n1"},
             {"name": "C", "hostname": "n2"}]  # fmt: skip
    out = m.runtime(manifest, trace, [])
    assert out["wall_time_s"] == pytest.approx(2550.0)
    assert out["sum_task_wall_s"] == pytest.approx(150.0)
    assert out["max_step_peak_rss_bytes"] == 2147483648
    assert out["cached_steps"] == ["bakta"] and out["runtime_complete"] is False
    assert [s["name"] for s in out["steps"]] == ["mapping", "bakta", "gene_table"]
    gt = out["steps"][2]
    assert gt["peak_rss_bytes"] is None and gt["wall_time_s"] == pytest.approx(30.0)
    assert sorted(out["hosts"]) == ["n1", "n2"]


def test_runtime_local_executor_and_sacct(m: Any) -> None:
    manifest = {
        "started_at": "2026-10-22T10:00:00+00:00", "finished_at": "2026-10-22T10:00:10+00:00",
        "exit_code": 0,
        "steps": [{"name": "mapping", "status": "completed", "n_tasks": 1, "n_retries": 0,
                   "wall_time_s": 5.0, "cpu_time_s": 4.0}],
    }  # fmt: skip
    out = m.runtime(manifest, [], [])
    assert out["hosts"] == [] and out["runtime_complete"] is True
    assert out["max_step_peak_rss_bytes"] is None
    sacct = [{"JobID": "123", "JobName": "giap-run-a2", "State": "COMPLETED",
              "Elapsed": "00:42:30", "MaxRSS": "2048K", "NodeList": "sh04-01n01",
              "AllocCPUS": "8"}]  # fmt: skip
    got = m.runtime(manifest, [], sacct)["sacct"]
    assert {str(v) for v in got.values()} >= {"00:42:30", "2048K", "sh04-01n01", "8", "COMPLETED"}


# ================================================================== TE1 on real data


G1 = _gene("G1", "chr", 0, 60000)
G2 = _gene("G2", "chr", 70000, 71000)
GM = _gene("GM", "chr", 100, 400, "-")
GO = _gene("GO", "c", 900, 1100, "+")  # on a 1000-bp circular contig, across the origin
JM = _gene("JMJJOG_00350", "chr", 80000, 81000)


def _te1_world(
    n_exact: int = 996,
    n_masked: int = 3,
    n_nonatg: int = 2,
    n_outframe: int = 1,
    n_nohit: int = 1,
    n_qcut: int = 0,
) -> dict[str, Any]:
    """|D2| = the sum of the counts; ``n_masked`` of the exact ones are ``rna_overlap`` in
    frame inside G1; ``n_qcut`` of the exact ones have every accepted row at q = 0.01."""
    rows, hits, unm = [], [], []
    non_atg: dict[str, str] = {}
    i = 0
    for j in range(n_exact):
        p = _pep(i)
        i += 1
        q = 0.01 if j >= n_exact - n_qcut else 0.001
        rows.append(_p(p, "G1", q))
        if j < n_masked:
            hits.append((p, "exact", "rna_overlap", "rRNA_1", "chr", 3 * j + 1, 3 * j + 24, "+"))
        else:
            hits.append((p, "exact", "cds_inframe", "G1", "chr", 3 * j + 1, 3 * j + 24, "+"))
    for _ in range(n_nonatg):
        p = _pep(i, "MEP")
        i += 1
        rows.append(_p(p, "G2"))
        non_atg[p] = "G2"
        hits.append((p, "tolerant", "cds_inframe", "G2", "chr", 70001, 70024, "+"))
    for _ in range(n_outframe):
        p = _pep(i)
        i += 1
        rows.append(_p(p, "G1"))
        hits.append((p, "exact", "cds_outframe", "G1", "chr", 2, 25, "+"))
    for _ in range(n_nohit):
        p = _pep(i)
        i += 1
        rows.append(_p(p, "G1"))
        unm.append((p, "unmapped", "no_hit"))
    return {"P": _peptides(rows), "hits": _hits(hits), "unm": _unmapped(unm),
            "genes": [G1, G2, GM, JM], "circular": {"chr": None}, "non_atg": non_atg}  # fmt: skip


def _te1(
    m: Any, cfg: Any, w: Mapping[str, Any], q_crossings: int = 0
) -> tuple[dict[str, Any], Any]:
    out = m.te1_real(w["P"], w["hits"], w["unm"], w["genes"], w["circular"], w["non_atg"], cfg,
                     q_crossings=q_crossings)  # fmt: skip
    te1, exceptions = out
    return te1, exceptions


def test_example_12_te1_pass(m: Any, cfg: Any) -> None:
    te1, exc = _te1(m, cfg, _te1_world())
    assert te1["n_peptides"] == 1000 and te1["exact_inframe"] == 996
    assert te1["share_exact_inframe"] == pytest.approx(0.996, abs=TOL)
    lo, hi = wilson_ci(996, 1000)
    assert te1["ci_low"] == pytest.approx(lo, abs=TOL) and te1["ci_high"] == pytest.approx(hi)
    assert (round(te1["ci_low"], 3), round(te1["ci_high"], 3)) == (0.990, 0.998)
    assert te1["passed"] is True
    assert te1["ci_low_below_threshold"] is True  # "CI lower end below 0.995"
    miss = te1["misses"]
    assert miss["non_atg_nterm"] == 2
    assert miss["exact_other_class"] == {"cds_outframe": 1}
    assert miss["unplaced"] == {"no_hit": 1} or (
        miss["unplaced"]["no_hit"] == 1 and sum(miss["unplaced"].values()) == 1
    )
    assert miss["exact_other_gene"] == 0 and miss["tolerant_only"] == 0
    assert te1["exact_inframe_masked"] == 3
    assert te1["non_atg_nterm"] == {"peptides": 2, "misses": 2}
    assert_wilson(te1["share_excluding_non_atg"], 996, 998, (0.993, 0.999))
    sq = te1["share_without_q_at_cutoff"]
    assert sq["value"] is None and sq["reason"] == "no_crossings"
    assert len(exc) == 4


def test_example_13_q_at_cutoff_superset(m: Any, cfg: Any) -> None:
    w = _te1_world(n_exact=995, n_masked=0, n_nonatg=0, n_outframe=0, n_nohit=5, n_qcut=5)
    te1, _ = _te1(m, cfg, w, q_crossings=3)
    assert te1["share_exact_inframe"] == pytest.approx(0.995, abs=TOL)
    assert (round(te1["ci_low"], 3), round(te1["ci_high"], 3)) == (0.988, 0.998)
    assert te1["passed"] is True  # the secondary numbers never change the verdict
    sq = te1["share_without_q_at_cutoff"]
    assert_wilson(sq, 990, 995)
    assert round(sq["value"], 5) == 0.99497
    assert sq["removed"] == 5 and sq["q_crossings"] == 3
    assert sq["pass_becomes_fail"] is True


def test_example_14_te1_fail_and_contaminants(m: Any, cfg: Any, monkeypatch: Any) -> None:
    w = _te1_world(n_exact=993, n_masked=0, n_nonatg=0, n_outframe=0, n_nohit=6)
    mixed = "MIXEDPEPK"
    extra = _peptides([_p(mixed, "Cont_P00722;JMJJOG_00350"), _p("WWWGGGWWK", "Cont_P00761")])
    w["P"] = pd.concat([w["P"], extra], ignore_index=True)
    hit = _hits([(mixed, "exact", "cds_inframe", "JMJJOG_00350", "chr", 80001, 80027, "+")])
    w["hits"] = pd.concat([w["hits"], hit], ignore_index=True)
    calls = _spy(m, monkeypatch)
    te1, _ = _te1(m, cfg, w)
    assert te1["n_peptides"] == 1000 and te1["exact_inframe"] == 994
    assert te1["share_exact_inframe"] == pytest.approx(0.994, abs=TOL)
    assert (round(te1["ci_low"], 3), round(te1["ci_high"], 3)) == (0.987, 0.997)
    assert te1["passed"] is False
    sources = max((c["sources"] for c in calls), key=len)
    assert set(sources[mixed]) == {"JMJJOG_00350"}
    assert "WWWGGGWWK" not in sources


def _spy(m: Any, monkeypatch: Any) -> list[dict[str, Any]]:
    """Wrap ``te_proxies.te1_result`` (in the proxy module and, if bound there, in metrics)."""
    te = load_te_proxies()
    real = te.te1_result
    calls: list[dict[str, Any]] = []

    def spy(*args: Any, **kwargs: Any) -> Any:
        names = ("sources", "excluded", "peptide_hits", "genes", "circular")
        bound = dict(zip(names, args, strict=False)) | kwargs
        calls.append(bound)
        return real(*args, **kwargs)

    monkeypatch.setattr(te, "te1_result", spy)
    if hasattr(m, "te1_result"):
        monkeypatch.setattr(m, "te1_result", spy)
    return calls


def test_criterion_3_te1_result_is_called_with_empty_excluded(
    m: Any, cfg: Any, monkeypatch: Any
) -> None:
    calls = _spy(m, monkeypatch)
    _te1(m, cfg, _te1_world())
    assert calls, "te1_result was not called"
    assert all(dict(c["excluded"]) == {} for c in calls)  # L5: nothing leaves the denominator
    assert max(len(c["sources"]) for c in calls) == 1000  # D2


def test_criterion_3_proxy_threshold_is_ignored(m: Any, cfg: Any) -> None:
    w = _te1_world(n_exact=997, n_masked=0, n_nonatg=0, n_outframe=0, n_nohit=3)
    te1, _ = _te1(m, cfg, w)
    assert te1["share_exact_inframe"] == pytest.approx(0.997)
    assert te1["passed"] is True  # the proxy's TE1_MIN 0.999 would fail it
    w = _te1_world(n_exact=994, n_masked=0, n_nonatg=0, n_outframe=0, n_nohit=6)
    assert _te1(m, cfg, w)[0]["passed"] is False


def test_criterion_3_na_q_target_in_d2(m: Any, cfg: Any) -> None:
    w = _te1_world(n_exact=10, n_masked=0, n_nonatg=0, n_outframe=0, n_nohit=0)
    w["P"] = pd.concat([w["P"], _peptides([_p("NAQPEPTK", "G1", None)])], ignore_index=True)
    te1, _ = _te1(m, cfg, w)
    assert te1["n_peptides"] == 11 and te1["target_q_na"] == 1


def test_te1_empty_denominator(m: Any, cfg: Any) -> None:
    w = _te1_world(n_exact=0, n_masked=0, n_nonatg=0, n_outframe=0, n_nohit=0)
    w["P"] = _peptides([_p("WWWGGGWWK", "Cont_P00761"), _p("AAAAGGGK", "G1", 0.02)])
    te1, exc = _te1(m, cfg, w)
    assert te1["n_peptides"] == 0 and te1["share_exact_inframe"] is None
    assert te1["passed"] is False and te1["reason"] == "empty_denominator"
    assert len(exc) == 0


def test_te1_miss_causes_first_match(m: Any, cfg: Any) -> None:
    peps = {k: _pep(i, "CAU") for i, k in enumerate(("other_gene", "tolerant", "nonatg_exact"))}
    P = _peptides([_p(peps["other_gene"], "G1"), _p(peps["tolerant"], "G1"),
                   _p(peps["nonatg_exact"], "G2")])  # fmt: skip
    hits = _hits([
        (peps["other_gene"], "exact", "cds_inframe", "G2", "chr", 70001, 70024, "+"),
        (peps["tolerant"], "tolerant", "cds_inframe", "G1", "chr", 1, 24, "+"),
        # a non-ATG product with an exact row at another gene: non_atg_nterm comes first
        (peps["nonatg_exact"], "exact", "cds_inframe", "G1", "chr", 1, 24, "+"),
    ])  # fmt: skip
    w = {"P": P, "hits": hits, "unm": _unmapped([]), "genes": [G1, G2], "circular": {"chr": None},
         "non_atg": {peps["nonatg_exact"]: "G2"}}  # fmt: skip
    te1, exc = _te1(m, cfg, w)
    miss = te1["misses"]
    assert miss["exact_other_gene"] == 1 and miss["tolerant_only"] == 1
    assert miss["non_atg_nterm"] == 1 and te1["exact_inframe"] == 0
    assert te1["non_atg_nterm"] == {"peptides": 1, "misses": 1}
    assert len(exc) == 3


def test_te1_exceptions_uncapped_sorted_and_kept_out_of_the_summary(m: Any, cfg: Any) -> None:
    w = _te1_world(n_exact=40, n_masked=0, n_nonatg=0, n_outframe=0, n_nohit=60)
    te1, exc = _te1(m, cfg, w)
    rows = list(exc.to_dict("records")) if isinstance(exc, pd.DataFrame) else list(exc)
    assert len(rows) == 60  # the proxy caps its examples at 50; the TSV does not
    keys = [(r["kind"], r["locus_tag"], r["peptide"]) for r in rows]
    assert keys == sorted(keys)
    assert {"kind", "locus_tag", "peptide", "hits"} <= set(rows[0])
    text = json.dumps(te1, default=str)
    assert not any(str(p) in text for p in w["P"]["peptide"])  # no sequence in metrics.json


# ---------------------------------------------------------------- criterion 4: coordinates


def _one(pep: str, gene_tag: str, hit: tuple[Any, ...], genes: list[Gene],
         circular: dict[str, int | None]) -> dict[str, Any]:  # fmt: skip
    return {"P": _peptides([_p(pep, gene_tag)]), "hits": _hits([hit]), "unm": _unmapped([]),
            "genes": genes, "circular": circular, "non_atg": {}}  # fmt: skip


def test_criterion_4_minus_strand_masked_hit_in_frame(m: Any, cfg: Any) -> None:
    pep = "MINUSPEK"  # 8 aa = 24 nt; 0-based [373, 397): (400 - 397) % 3 == 0
    w = _one(pep, "GM", (pep, "exact", "rna_overlap", "tRNA_1", "chr", 374, 397, "-"),
             [GM], {"chr": None})  # fmt: skip
    te1, _ = _te1(m, cfg, w)
    assert te1["exact_inframe"] == 1 and te1["exact_inframe_masked"] == 1
    w = _one(pep, "GM", (pep, "exact", "rna_overlap", "tRNA_1", "chr", 375, 398, "-"),
             [GM], {"chr": None})  # fmt: skip
    te1, _ = _te1(m, cfg, w)
    assert te1["exact_inframe"] == 0
    assert te1["misses"]["exact_other_class"] == {"rna_overlap": 1}


def test_criterion_4_origin_crossing_masked_hit_in_frame(m: Any, cfg: Any) -> None:
    pep = "ORIGINPK"
    # gene [900, 1100) on a 1000-bp circular contig; the hit is reported after the origin,
    # 0-based [2, 26): shifted by +L it is [1002, 1026), (1002 - 900) % 3 == 0
    hit = (pep, "exact", "rna_overlap", "rRNA_9", "c", 3, 26, "+")
    te1, _ = _te1(m, cfg, _one(pep, "GO", hit, [GO], {"c": 1000}))
    assert te1["exact_inframe"] == 1 and te1["exact_inframe_masked"] == 1
    te1, _ = _te1(m, cfg, _one(pep, "GO", hit, [GO], {"c": None}))  # linear: no shift
    assert te1["exact_inframe"] == 0
    shifted = (pep, "exact", "rna_overlap", "rRNA_9", "c", 4, 27, "+")
    te1, _ = _te1(m, cfg, _one(pep, "GO", shifted, [GO], {"c": 1000}))
    assert te1["exact_inframe"] == 0
    assert te1["misses"]["exact_other_class"] == {"rna_overlap": 1}


# ================================================================== hypothesis 1


def _h1_inputs(
    spec: Sequence[tuple[str, str | None, str, str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, set[str]]]:
    """Rows (locus, tu_id or None, A1 status, A2 status, region) with region in {"", "plasmid",
    "prophage", "bgc", "strain"}; returns genes_a1, genes_a2, compare_detection, subsets."""
    frames = []
    subsets: dict[str, set[str]] = {"plasmid": set(), "prophage": set(), "bgc": set(),
                                    "strain_specific": set()}  # fmt: skip
    for arm in (2, 3):
        rows = []
        for i, r in enumerate(spec):
            tag, tu, region = r[0], r[1], r[4]
            status = r[arm]
            row = gene(tag, start=i * 400, detected=status == "detected",
                       detection_status=status,
                       replicon_type="plasmid" if region == "plasmid" else "chromosome",
                       contig="pMUT1" if region == "plasmid" else "chr1",
                       region_ids=REGION_ID.get(region))  # fmt: skip
            row["tu_id"] = tu
            rows.append(row)
        frames.append(genes_frame(rows))
    for tag, _, _, _, region in spec:
        if region in ("plasmid", "prophage", "bgc"):
            subsets[region].add(tag)
        if region == "strain":
            subsets["strain_specific"].add(tag)
    det = []
    for tag, _, s1, s2, _ in spec:
        for label, s in (("A1", s1), ("A2", s2)):
            det.append({"section": "matrix", "locus_tag": tag, "run": label,
                        "detected": s == "detected", "detection_status": s,
                        "intersection": None, "n_genes": None})  # fmt: skip
    cmp = pd.DataFrame(det, columns=list(DETECTION_COLUMNS))
    return frames[0], frames[1], cmp, subsets


D, N = "detected", "not_detected"
REGION_ID = {"prophage": "prophage_1", "bgc": "bgc_1"}


def _ex7(n_plasmid: int = 7) -> list[tuple[str, str | None, str, str, str]]:
    spec: list[tuple[str, str | None, str, str, str]] = []
    for i in range(12):  # A2-only, empty tu_id
        region = "plasmid" if i < n_plasmid else ("strain" if i < 12 else "")
        spec.append((f"g{i:03d}", None, N, D, region))
    spec += [("g012", "TU_c1", D, N, ""), ("g013", "TU_c2", D, N, "")]
    spec += [(f"g{i:03d}", f"TU_{i}", D, D, "") for i in range(14, 60)]
    spec += [(f"g{i:03d}", f"TU_{i}", N, N, "") for i in range(60, 100)]
    return spec


def _h1(m: Any, spec: Sequence[tuple[str, str | None, str, str, str]]) -> dict[str, Any]:
    a1, a2, cmp, subsets = _h1_inputs(spec)
    out: dict[str, Any] = m.hypothesis1(a1, a2, cmp, subsets)
    return out


def test_example_7_h1_genes_are_tus(m: Any) -> None:
    out = _h1(m, _ex7())
    assert (out["b"], out["c"]) == (12, 2)
    assert_wilson(out["prop_a2_only"], 12, 14, (0.601, 0.960))
    assert out["coverage_diff"] == pytest.approx(0.10, abs=TOL)
    assert (out["b_tu"], out["c_tu"], out["mixed_tu"]) == (12, 2, 0)
    want = binomtest(12, 14, 0.5, alternative="greater").pvalue
    assert out["p_tu"] == pytest.approx(want, abs=TOL) == pytest.approx(106 / 16384, abs=TOL)
    assert out["p_gene"] == pytest.approx(want, abs=TOL)
    assert_wilson(out["mostly"], 7, 12, (0.320, 0.807))
    assert out["mostly_base_rate"]["value"] == pytest.approx(7 / 100, abs=TOL)
    assert out["strain_specific_share"]["value"] == pytest.approx(5 / 12, abs=TOL)
    assert out["verdict"] == "supported" and list(out["failed_clauses"]) == []


def _ex8() -> list[tuple[str, str | None, str, str, str]]:
    spec: list[tuple[str, str | None, str, str, str]] = []
    spec += [(f"g{i:03d}", "TU_A", N, D, "plasmid") for i in range(4)]
    spec += [(f"g{i:03d}", "TU_B", N, D, "plasmid") for i in range(4, 6)]
    spec += [("g006", "TU_C", D, N, ""), ("g007", "TU_D", D, N, "")]
    spec += [(f"g{i:03d}", f"TU_{i}", D, D, "") for i in range(8, 20)]
    return spec


def test_example_8_h1_tu_clustering(m: Any) -> None:
    out = _h1(m, _ex8())
    assert (out["b"], out["c"]) == (6, 2)
    assert_wilson(out["prop_a2_only"], 6, 8, (0.409, 0.929))
    assert (out["b_tu"], out["c_tu"]) == (2, 2)
    p_tu = binomtest(2, 4, 0.5, alternative="greater").pvalue
    assert out["p_tu"] == pytest.approx(p_tu, abs=TOL) == pytest.approx(11 / 16, abs=TOL)
    assert_wilson(out["unit_effect"], 2, 4, (0.150, 0.850))
    p_gene = binomtest(6, 8, 0.5, alternative="greater").pvalue
    assert out["p_gene"] == pytest.approx(p_gene, abs=TOL) == pytest.approx(37 / 256, abs=TOL)
    assert out["mostly"]["value"] == pytest.approx(1.0)
    # criterion 1: p_tu 0.6875 and not supported although p_gene is 0.145 and b > c
    assert out["verdict"] == "not_supported" and list(out["failed_clauses"]) == ["test"]


def test_criterion_2a_only_detected_counts(m: Any) -> None:
    spec = _ex7()
    spec[0] = ("g000", None, N, "single_peptide_unestimated", "plasmid")
    out = _h1(m, spec)
    assert out["b"] == 11 and out["b_tu"] == 11


def test_criterion_2a_mixed_tu(m: Any) -> None:
    spec = [*_ex7(), ("g100", "TU_M", N, D, ""), ("g101", "TU_M", D, N, "")]
    out = _h1(m, spec)
    assert (out["b"], out["c"]) == (13, 3)  # kept in b, c and p_gene
    assert (out["b_tu"], out["c_tu"], out["mixed_tu"]) == (12, 2, 1)
    want = binomtest(13, 16, 0.5, alternative="greater").pvalue
    assert out["p_gene"] == pytest.approx(want, abs=TOL)
    assert out["p_tu"] == pytest.approx(106 / 16384, abs=TOL)


def test_criterion_2a_direction_fails_alone(m: Any) -> None:
    spec: list[tuple[str, str | None, str, str, str]] = []
    spec += [(f"g{i:03d}", None, N, D, "plasmid") for i in range(8)]
    spec += [(f"g{i:03d}", "TU_C", D, N, "") for i in range(8, 18)]
    spec += [(f"g{i:03d}", f"TU_{i}", D, D, "") for i in range(18, 30)]
    out = _h1(m, spec)
    assert (out["b"], out["c"], out["b_tu"], out["c_tu"]) == (8, 10, 8, 1)
    assert out["p_tu"] == pytest.approx(binomtest(8, 9, 0.5, alternative="greater").pvalue)
    assert out["p_tu"] <= 0.05 and out["mostly"]["value"] == pytest.approx(1.0)
    assert out["verdict"] == "not_supported" and list(out["failed_clauses"]) == ["direction"]


def test_criterion_2a_mostly_fails_alone(m: Any) -> None:
    out = _h1(m, _ex7(n_plasmid=5))
    assert out["mostly"]["value"] == pytest.approx(5 / 12, abs=TOL)
    assert out["verdict"] == "not_supported" and list(out["failed_clauses"]) == ["mostly"]


def test_h1_not_evaluable_without_discordant_units(m: Any) -> None:
    spec = [("g0", "TU_M", N, D, "plasmid"), ("g1", "TU_M", D, N, ""),
            *[(f"g{i}", f"TU_{i}", D, D, "") for i in range(2, 10)]]  # fmt: skip
    out = _h1(m, spec)
    assert (out["b_tu"], out["c_tu"], out["mixed_tu"]) == (0, 0, 1)
    assert out["p_tu"] is None and out["verdict"] == "not_evaluable"


def test_h1_compare_matrix_cross_check_warns(m: Any) -> None:
    a1, a2, cmp, subsets = _h1_inputs(_ex7())
    assert not m.hypothesis1(a1, a2, cmp, subsets).get("warnings")
    sel = (cmp["locus_tag"] == "g000") & (cmp["run"] == "A1")
    cmp.loc[sel, "detected"] = True
    cmp.loc[sel, "detection_status"] = D
    out = m.hypothesis1(a1, a2, cmp, subsets)
    assert out["b"] == 12  # genes.tsv decides, the matrix only cross-checks
    assert out.get("warnings")


# ================================================================== render


METRICS_SMALL: dict[str, Any] = {
    "te1": {"share_exact_inframe": 0.996, "ci_low": 0.9897604437227377,
            "ci_high": 0.9984434118595916, "n_peptides": 1000},
    "h1": {"p_tu": 0.6875, "p_gene": 0.14453125},
    "m7": {"A1": {"wall_time_s": 2550.0, "max_step_peak_rss_bytes": 2147483648}},
}  # fmt: skip


def test_render_number_with_its_source(m: Any) -> None:
    out = m.render("TE1 {{te1.share_exact_inframe}} end", METRICS_SMALL, None)
    source = "(metrics.json \u203a te1)"  # spec: "0.996 [0.990, 0.998] (metrics.json > te1)"
    assert "0.996" in out and source in out
    assert out.index("0.996") < out.index(source)


@pytest.mark.parametrize(
    ("path", "text"),
    [
        ("h1.p_tu", "0.69"),  # p: 2 significant digits
        ("h1.p_gene", "0.14"),
        ("m7.A1.wall_time_s", "0:42:30"),  # h:mm:ss
        ("m7.A1.max_step_peak_rss_bytes", "2.0 GiB"),  # GiB, 1 decimal
    ],
)
def test_render_formats(m: Any, path: str, text: str) -> None:
    out = m.render("x {{" + path + "}} y", METRICS_SMALL, None)
    assert text in out, out


def test_render_missing_path_is_marked(m: Any) -> None:
    out = m.render("a {{te1.no_such_key}} b", METRICS_SMALL, None)
    assert "**missing: te1.no_such_key**" in out


@pytest.mark.parametrize(
    "notes",
    [
        "The misses come from PEPTIDEK and friends.",  # a peptide-like string
        "We saw 1205 peptides.",  # a 4-digit number not in metrics.json
    ],
)
def test_render_rejects_notes(m: Any, notes: str) -> None:
    with pytest.raises(m.Inconsistent):  # review B1: a rejected notes file is the exit-1 class
        m.render("# Report\n", METRICS_SMALL, notes)


def test_render_accepts_notes_with_numbers_from_metrics(m: Any) -> None:
    notes = "Of 1000 peptides the misses are non-ATG starts; the arm took 2550 s."
    out = m.render("# Report\n", METRICS_SMALL, notes)
    assert notes in out


def test_wilson_ci_reference_values_match_the_spec() -> None:
    """Self-check of the oracle: every printed CI of examples 4-14 is wilson_ci's."""
    printed = {(7, 10): (0.397, 0.892), (60, 100): (0.502, 0.691), (50, 100): (0.404, 0.596),
               (8, 10): (0.490, 0.943), (2, 10): (0.057, 0.510), (12, 14): (0.601, 0.960),
               (7, 12): (0.320, 0.807), (6, 8): (0.409, 0.929), (2, 4): (0.150, 0.850),
               (996, 1000): (0.990, 0.998), (996, 998): (0.993, 0.999),
               (995, 1000): (0.988, 0.998), (994, 1000): (0.987, 0.997)}  # fmt: skip
    for (k, n), (lo, hi) in printed.items():
        wlo, whi = wilson_ci(k, n)
        assert (round(wlo, 3), round(whi, 3)) == (lo, hi), (k, n)
    assert math.isclose(binomtest(12, 14, 0.5, alternative="greater").pvalue, 106 / 16384)


# ================================================================== review-round fixes
# Code review (rev-1005) and statistician review (stat-1005) of 83ba1a99; lead decisions of the
# review round. One block per finding; S10 rows name the mutant of rev-1005/mutants.log they kill.


def _miss_total(misses: Mapping[str, Any]) -> int:
    return sum(sum(v.values()) if isinstance(v, Mapping) else int(v) for v in misses.values())


def _add(
    w: dict[str, Any],
    rows: Sequence[tuple[Any, ...]] = (),
    hits: Sequence[tuple[Any, ...]] = (),
    unm: Sequence[tuple[str, str, str]] = (),
) -> None:
    """Append peptide rows, hit rows and unmapped rows to a ``_te1_world``."""
    w["P"] = pd.concat([w["P"], _peptides(rows)], ignore_index=True)
    if hits:
        w["hits"] = pd.concat([w["hits"], _hits(hits)], ignore_index=True)
    if unm:
        w["unm"] = pd.concat([w["unm"], _unmapped(unm)], ignore_index=True)


# ---------------------------------------------------------------- S1: unknown_source


def test_s1_unknown_source_is_a_separate_count(m: Any, cfg: Any) -> None:
    """A peptide whose sources hold no locus tag of ``genes`` (a UniProt-like ID, or no ID at
    all) is counted in ``te1.unknown_source``; it is still a miss under its first-match cause."""
    w = _te1_world(n_exact=4, n_masked=0, n_nonatg=0, n_outframe=0, n_nohit=0)
    uni, empty = _pep(900, "UNI"), _pep(901, "EMP")
    _add(w, [_p(uni, "P12345"), _p(empty, None)],  # type: ignore[arg-type]
         hits=[(uni, "exact", "cds_inframe", "G1", "chr", 1, 24, "+")],
         unm=[(empty, "unmapped", "no_hit")])  # fmt: skip
    te1, exc = _te1(m, cfg, w)
    assert te1["n_peptides"] == 6 and te1["exact_inframe"] == 4
    assert te1["unknown_source"] == 2
    assert "unknown_source" not in te1["misses"]  # a count beside the cause chain
    assert te1["n_misses"] == 2 == _miss_total(te1["misses"]) == len(exc)


def test_s1_unknown_source_zero_on_example_12(m: Any, cfg: Any) -> None:
    te1, _ = _te1(m, cfg, _te1_world())
    assert te1["unknown_source"] == 0
    assert _miss_total(te1["misses"]) == te1["n_misses"] == 4


# ---------------------------------------------------------------- S2: validate_peptides first


def test_s2_universe_goes_through_validate_peptides(m: Any, cfg: Any) -> None:
    from giap.io.peptides import validate_peptides

    P = _peptides([
        _p("PEPTIDEK", "G1"),
        _p("REVPEPTK", "rev_G1", decoy=None),  # type: ignore[arg-type]  # NA is_decoy: a decoy by prefix
        _p("PEPT1DEK", "G1"),  # not A-Z only: validate_peptides drops it
        _p("WWWGGGWWK", " Cont_P00761 ; Cont_P02769"),  # whitespace normalised: contaminant-only
    ])  # fmt: skip
    valid, _ = validate_peptides(P, cfg.input)
    oracle = set(select_queries(valid, cfg)["peptide"].astype(str))  # GIAP's mapping front end
    assert oracle == {"PEPTIDEK", "WWWGGGWWK"}
    hits = _hits([("PEPTIDEK", "exact", "cds_inframe", "G1", "chr", 1, 24, "+")])
    out = m.mapping_shares(P, hits, _unmapped([]), cfg)
    assert (out["n_universe"], out["n_contaminant_only"], out["n_denominator"]) == (2, 1, 1)
    te1, _ = m.te1_real(P, hits, _unmapped([]), [G1], {"chr": None}, {}, cfg)
    assert te1["n_peptides"] == 1 and te1["exact_inframe"] == 1
    assert m.te1_diagnostic(P, hits, cfg)["n_peptides"] == 1


# ---------------------------------------------------------------- S6: reason, not a made-up flag


def test_s6_missing_compare_rows_carry_a_reason_and_no_flag(m: Any) -> None:
    ab = _abundance()
    m4 = m.abundance_agreement(ab[ab["run_b"] != "A2"], ("A1", "A2"))  # only the A1/A3 row
    assert m4["reason"] == "missing_row" and m4["flags"] is None
    assert m4["spearman_rho"] is None and m4["ci_low"] is None
    m5 = m.stability(_stab([r for r in EX9 if r[0] != "islands"]), ("A1", "A2"))
    assert m5["islands"].get("reason") == "missing_row" and m5["islands"]["flags"] is None
    h3 = m.hypothesis3(m4, m5)
    assert _status(h3["criteria"]["rho"]) == "not_evaluable"
    assert _status(h3["criteria"]["islands"]) == "not_evaluable"
    assert h3["verdict"] == "not_evaluable"


# ---------------------------------------------------------------- H3 (should-fix 5, decision 6)


def test_h3_rho_present_but_ci_low_na_is_not_evaluable(m: Any) -> None:
    out = _h3(m, _abundance(ci_low=np.nan, ci_high=np.nan), EX9)
    assert _status(out["criteria"]["rho"]) == "not_evaluable"
    assert out["verdict"] == "not_evaluable"


def test_h3_criteria_carry_rho_ci_and_family_counts_even_when_not_evaluable(m: Any) -> None:
    """Lead decision 6: H3 is descriptive; rho with its CI and every family's n_a, n_b,
    n_shared, union and J are reported whatever the verdict."""
    out = _h3(m, _abundance(flags="too_few_shared"), EX10[:3])
    rho = out["criteria"]["rho"]
    assert _status(rho) == "not_evaluable" and out["verdict"] == "not_evaluable"
    assert rho["rho"] == pytest.approx(0.83)
    assert rho["ci_low"] == pytest.approx(0.76) and rho["ci_high"] == pytest.approx(0.88)
    isl = out["criteria"]["islands"]
    assert (isl["n_a"], isl["n_b"], isl["n_shared"], isl["union"]) == (10, 8, 6, 12)
    assert isl["jaccard"] == pytest.approx(0.5)
    op = out["criteria"]["operon_flags"]
    assert _status(op) == "not_evaluable"
    assert (op["n_a"], op["n_b"], op["n_shared"], op["union"]) == (0, 0, 0, 0)
    assert op["jaccard"] is None
    ec = out["criteria"]["enriched_categories"]
    assert ec["n_a"] is None and ec["n_b"] == 3 and ec["union"] is None


def test_s10_h3_union_exactly_5_is_evaluable(m: Any) -> None:  # mutant 52
    out = _h3(m, _abundance(), [("islands", 4, 4, 3, 0.6, None), *EX9[1:]])
    assert out["criteria"]["islands"]["union"] == 5
    assert _status(out["criteria"]["islands"]) == "holds"
    assert out["verdict"] == "supported"


# ---------------------------------------------------------------- TE1 (should-fix 6, S10)


def test_q_at_cutoff_empty_denominator_keeps_the_flag_with_a_reason(m: Any, cfg: Any) -> None:
    w = _te1_world(n_exact=5, n_masked=0, n_nonatg=0, n_outframe=0, n_nohit=0, n_qcut=5)
    te1, _ = _te1(m, cfg, w, q_crossings=2)
    assert te1["passed"] is True
    sq = te1["share_without_q_at_cutoff"]
    assert sq["value"] is None and sq["removed"] == 5
    assert sq["pass_becomes_fail"] is True and sq["reason"] == "empty_denominator"


def test_s10_q_at_cutoff_missed_and_mixed_rows(m: Any, cfg: Any) -> None:  # mutants 31, 32
    """Example 13 plus one missed peptide at q = 0.01 and one exact peptide with rows at 0.01
    and 0.001 (not all rows at the cut-off: stays in the superset's denominator)."""
    w = _te1_world(n_exact=995, n_masked=0, n_nonatg=0, n_outframe=0, n_nohit=5, n_qcut=5)
    cut, mix = _pep(950, "CUT"), _pep(951, "MIX")
    rows = [_p(cut, "G1", 0.01), _p(mix, "G1", 0.01, run="WT_1"), _p(mix, "G1", 0.001, run="WT_2")]
    _add(w, rows, hits=[(mix, "exact", "cds_inframe", "G1", "chr", 1, 24, "+")],
         unm=[(cut, "unmapped", "no_hit")])  # fmt: skip
    te1, _ = _te1(m, cfg, w, q_crossings=3)
    assert (te1["n_peptides"], te1["exact_inframe"]) == (1002, 996)
    sq = te1["share_without_q_at_cutoff"]
    assert sq["removed"] == 6  # 5 exact + 1 missed, all rows at 0.01; not the mixed one
    assert_wilson(sq, 991, 996)


def test_s10_non_atg_product_that_maps_leaves_both_sides(m: Any, cfg: Any) -> None:  # mutant 35
    w = _te1_world(n_exact=4, n_masked=0, n_nonatg=1, n_outframe=0, n_nohit=0)
    mapped = _pep(960, "MEP")
    _add(
        w,
        [_p(mapped, "G2")],
        hits=[(mapped, "exact", "cds_inframe", "G2", "chr", 70001, 70024, "+")],
    )
    w["non_atg"] = {**w["non_atg"], mapped: "G2"}
    te1, _ = _te1(m, cfg, w)
    assert (te1["n_peptides"], te1["exact_inframe"]) == (6, 5)
    assert te1["non_atg_nterm"] == {"peptides": 2, "misses": 1}
    assert_wilson(te1["share_excluding_non_atg"], 4, 4)


def test_s10_target_q_na_needs_every_row_na(m: Any, cfg: Any) -> None:  # mutant 36
    w = _te1_world(n_exact=3, n_masked=0, n_nonatg=0, n_outframe=0, n_nohit=0)
    half, alln = _pep(970, "HAL"), _pep(971, "ALL")
    _add(w, [_p(half, "G1", None, run="WT_1"), _p(half, "G1", 0.001, run="WT_2"),
             _p(alln, "G1", None)])  # fmt: skip
    te1, _ = _te1(m, cfg, w)
    assert te1["n_peptides"] == 5 and te1["target_q_na"] == 1


def test_s10_exceptions_sorted_by_locus_tag_before_peptide(m: Any, cfg: Any) -> None:  # mutant 73
    late, early = _pep(0, "CAA"), _pep(0, "DAA")  # CAA... < DAA..., but its locus tag is later
    P = _peptides([_p(late, "G2"), _p(early, "G1")])
    w = {"P": P, "hits": _hits([]), "unm": _unmapped([(late, "unmapped", "no_hit"),
                                                      (early, "unmapped", "no_hit")]),
         "genes": [G1, G2], "circular": {"chr": None}, "non_atg": {}}  # fmt: skip
    _, exc = _te1(m, cfg, w)
    rows = list(exc.to_dict("records")) if isinstance(exc, pd.DataFrame) else list(exc)
    assert [(r["locus_tag"], r["peptide"]) for r in rows] == [("G1", early), ("G2", late)]


def test_s10_a1_diagnostic_counts_in_frame_rows_only(m: Any, cfg: Any) -> None:  # mutant 38
    a, b = _pep(980, "DIA"), _pep(981, "DIA")
    P = _peptides([_p(a, "P00001"), _p(b, "P00002")])
    hits = _hits([(a, "exact", "cds_inframe", "G1", "chr", 1, 24, "+"),
                  (b, "exact", "cds_outframe", "G1", "chr", 2, 25, "+")])  # fmt: skip
    diag = m.te1_diagnostic(P, hits, cfg)
    assert diag["n_peptides"] == 2 and diag["exact_inframe_any_gene"] == 1
    assert_wilson(diag["share"], 1, 2)


# ---------------------------------------------------------------- metrics 1, 2 (should-fix 7, S10)


def test_metric_1_skips_rows_without_protein_ids_in_groups(m: Any) -> None:
    P = _peptides([_p("PEPTIDEK", "G1"), _p("QQQNNNK", None)])  # type: ignore[arg-type]
    u = m.identifications(P, PREFIXES)["union"]
    assert u["peptides"]["all"] == 2
    assert u["protein_groups"]["all"] == 1 and u["protein_groups"]["non_contaminant"] == 1


def test_s10_contaminant_row_last_does_not_win(m: Any) -> None:  # mutants 2, 3
    P = _peptides([_p("SAMPLEPEPK", "JMJJOG_00350", run="WT_1"),
                   _p("SAMPLEPEPK", "Cont_P00722", run="WT_2")])  # fmt: skip
    u = m.identifications(P, PREFIXES)["union"]
    assert u["peptides"]["contaminant_only"] == 0 and u["peptides"]["non_contaminant"] == 1
    assert u["precursors"]["contaminant_only"] == 0 and u["precursors"]["all"] == 1


def test_s10_d_keeps_a_peptide_with_one_strain_row(m: Any, cfg: Any) -> None:  # mutants 7, 10
    a, mixed = _pep(990, "STR"), _pep(991, "MXD")
    P = _peptides([_p(a, "G1"), _p(mixed, "Cont_P00722", run="WT_1"), _p(mixed, "G1", run="WT_2"),
                   _p("WWWGGGWWK", "Cont_P02769")])  # fmt: skip
    hits = _hits([(a, "exact", "cds_inframe", "G1", "chr", 1, 24, "+"),
                  (mixed, "exact", "cds_inframe", "G1", "chr", 1, 24, "+")])  # fmt: skip
    out = m.mapping_shares(P, hits, _unmapped([]), cfg)  # unmapped.tsv lists no contaminant
    assert (out["n_universe"], out["n_contaminant_only"], out["n_denominator"]) == (3, 1, 2)
    assert out["counts"]["exact"] == 2
    assert out["contaminant_difference"] == 1  # |C| = 1, none listed in unmapped.tsv


# ---------------------------------------------------------------- metric 7 (S10)


def test_s10_runtime_max_rss_sacct_head_and_trace_hosts(m: Any) -> None:  # mutants 18, 20, 21
    manifest = {
        "started_at": "2026-10-22T10:00:00+00:00", "finished_at": "2026-10-22T11:00:00+00:00",
        "exit_code": 0,
        "steps": [
            {"name": "a", "status": "completed", "wall_time_s": 5.0, "peak_rss_bytes": 2 * GIB},
            {"name": "b", "status": "completed", "wall_time_s": 6.0, "peak_rss_bytes": 3 * GIB},
        ],
    }  # fmt: skip
    trace = [{"hostname": "-"}, {"hostname": ""}, {"hostname": "n1"}, {"hostname": "None"}]
    sacct = [{"JobID": "123.batch", "State": "COMPLETED", "Elapsed": "00:59:00"},
             {"JobID": "123", "JobName": "giap-run-a2", "State": "COMPLETED",
              "Elapsed": "01:00:00", "NodeList": "sh04-01n01"},
             {"JobID": "123.extern", "State": "COMPLETED", "Elapsed": "01:00:01"}]  # fmt: skip
    out = m.runtime(manifest, trace, sacct)
    assert out["max_step_peak_rss_bytes"] == 3 * GIB  # the max of the steps, not their sum
    assert out["hosts"] == ["n1"]
    assert out["sacct"]["JobID"] == "123" and out["sacct"]["Elapsed"] == "01:00:00"


def test_s10_render_hours_minutes_seconds(m: Any) -> None:  # mutant 60
    out = m.render("x {{m7.A1.wall_time_s}} y", {"m7": {"A1": {"wall_time_s": 3725.0}}}, None)
    assert "1:02:05" in out


# ---------------------------------------------------------------- H1 (S10)


def test_s10_h1_mostly_exactly_half_fails(m: Any) -> None:  # mutant 44
    out = _h1(m, _ex7(n_plasmid=6))
    assert out["mostly"]["value"] == pytest.approx(0.5, abs=TOL)
    assert out["verdict"] == "not_supported" and list(out["failed_clauses"]) == ["mostly"]


def test_s10_h1_bgc_genes_count_for_mostly(m: Any) -> None:  # mutant 49
    spec = [(t, tu, s1, s2, "bgc" if i < 7 else r) for i, (t, tu, s1, s2, r) in
            enumerate(_ex7(n_plasmid=0))]  # fmt: skip
    out = _h1(m, spec)
    assert_wilson(out["mostly"], 7, 12)
    assert out["verdict"] == "supported"


def test_s10_h1_not_evaluable_has_no_failed_clauses(m: Any) -> None:  # mutant 46
    spec = [("g0", "TU_M", N, D, "plasmid"), ("g1", "TU_M", D, N, ""),
            *[(f"g{i}", f"TU_{i}", D, D, "") for i in range(2, 10)]]  # fmt: skip
    out = _h1(m, spec)
    assert out["verdict"] == "not_evaluable" and list(out["failed_clauses"]) == []


# ---------------------------------------------------------------- checklist (S10)


def test_s10_checklist_a3_failed_te1_no_data_and_m4_missing_row(m: Any) -> None:  # 56-58
    cl = m.checklist({"te1": {"n_peptides": 0, "passed": False}, "m4": {"reason": "missing_row"}},
                     {"A1": 0, "A2": 0, "A3": 1, "compare_ok": True})  # fmt: skip
    assert [_status(c) for c in cl[:3]] == ["fail", "fail (no data)", "fail"]
    cl = m.checklist({"te1": {"n_peptides": 10, "passed": False}, "m4": {"spearman_rho": 0.8}},
                     {"A1": 0, "A2": 0, "A3": None, "compare_ok": True})  # fmt: skip
    assert _status(cl[1]) == "fail" and _status(cl[2]) == "pass"
    assert _status(cl[0]).startswith("partial")


# ---------------------------------------------------------------- render (S9, S8, S10)


def test_s10_notes_reject_a_7_letter_peptide(m: Any) -> None:  # mutant 62
    with pytest.raises(m.Inconsistent):
        m.render("# Report\n", METRICS_SMALL, "We saw QQWWRAK here.")


def _first_token(text: str) -> str:
    return text.split()[1]  # "x <value> ..." -> <value>


@pytest.mark.parametrize(
    ("doc", "path"),
    [
        ({"te1": {"share_exact_inframe": 0.9996}}, "te1.share_exact_inframe"),
        ({"te1": {"ci_low": 0.9996}}, "te1.ci_low"),
        ({"m2": {"A2": {"shares": {"exact": {"value": 0.9996, "ci_low": 0.99, "ci_high": 0.9999,
                                              "k": 9996, "n": 10000}}}}}, "m2.A2.shares.exact"),
    ],
)  # fmt: skip
def test_s9_share_0_9996_renders_with_4_significant_figures(
    m: Any, doc: dict[str, Any], path: str
) -> None:
    out = m.render("x {{" + path + "}} y", doc, None)
    assert _first_token(out) == "0.9996", out


def test_s9_shares_near_one_never_render_as_one(m: Any) -> None:
    from hypothesis import given, settings
    from hypothesis import strategies as st

    @settings(max_examples=200, deadline=None)
    @given(st.floats(min_value=0.9995, max_value=0.99999))
    def prop(v: float) -> None:
        for doc, path in (
            ({"te1": {"share_exact_inframe": v}}, "te1.share_exact_inframe"),
            ({"x": {"s": {"value": v, "ci_low": None, "ci_high": None}}}, "x.s"),
        ):
            tok = _first_token(m.render("x {{" + path + "}} y", doc, None))
            assert tok not in ("1", "1.0"), (v, tok)
            assert abs(float(tok) - v) <= 5e-5, (v, tok)

    prop()


def test_render_rho_and_p_values_keep_their_formats(m: Any) -> None:
    out = m.render("x {{m4.spearman_rho}} y", {"m4": {"spearman_rho": 0.83}}, None)
    assert _first_token(out) == "0.83"


def test_s8_render_breakdown_dicts_as_k_v_pairs(m: Any) -> None:
    doc = {"te1": {"misses": {"unplaced": {"no_hit": 2, "too_short": 1}, "exact_other_class": {}}}}
    out = m.render("a {{te1.misses.unplaced}} b {{te1.misses.exact_other_class}} c", doc, None)
    assert "no_hit: 2" in out and "too_short: 1" in out
    assert "{'" not in out and "{}" not in out
    assert "b none" in out
