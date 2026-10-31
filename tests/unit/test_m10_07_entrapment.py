"""Universe, ``present``, precursors and discoveries of ``validation/pilot/entrapment.py`` (M10-07
acceptance criteria 1, 2, 4 and 6; examples 1, 1b, 1c, 2, 8, 9; the edge cases).

Spec: docs/specs/M10-07-free-entrapment.md ("The searchable space", "Definitions", "Examples",
"Edge cases") as amended by the binding "Statistician decisions" (S1: the assessment only at
precursor level, peptide and per-run blocks ``not_assessed``; S2: ``sensitivity.present_broad``
and example 1c; S3: the ``ci`` object; S4: ``sensitivity.r_prec_weighted`` and
``sensitivity.far``; S5: I/L collapse). Lead decisions L1-L7 take their defaults (L2: the
``strain_db`` substring rule is on). Inputs are the spec's hand-made sequences
(``tests/helpers/entrapment_m10_07.py``); nothing comes from a DIA-NN report.

Interface fixed here (the spec names the quantities, not the functions):

* ``Settings(enzyme="trypsin_p", missed=1, min_len=7, max_len=30, charges=(1, 4),
  mz=(300.0, 1800.0), var_mods=1)`` (a frozen dataclass; ``charges`` and ``mz`` are inclusive
  ``(min, max)`` as DIA-NN's ``--min/--max-pr-*`` flags) and ``DEFAULT_SETTINGS = Settings()``;
* ``norm(s)`` (I -> L); ``digest_set(proteins, settings=DEFAULT_SETTINGS) -> set[str]`` (D(F):
  normalised, standard-residue products);
* ``mz(peptide, charge, oxidised_at=None) -> float`` (carbamidomethyl on every C; ``oxidised_at``
  the 0-based index of the one oxidised M); ``n_prec(peptide, settings=DEFAULT_SETTINGS) -> int``;
* ``build_index(genome, replicons, strain_db=None, tolerant=None, *, annotation)`` (``genome``:
  contig -> nucleotides, ``replicons``: contig -> ``gippa.model.coords.Replicon``, ``strain_db``:
  protein sequences or None for "L2 answered no"; ``annotation``: the required
  ``gippa.model.genes.GeneModel``, ADR-0068), with the start cells in ``index.starts``
  (``gippa.mapping.initiator.StartCells``), built with ``resolve_config().annotation``'s
  ``start_codons`` and ``translation_table``;
* ``present(peptides, index) -> Mapping[str, Presence]`` keyed by the peptides passed, with
  ``Presence.rule`` in ``{"exact", "strain_db", "mass_equivalent", None}``; it is the one
  classification callable, looked up as a module global by ``evaluate``;
* ``evaluate(k12, contaminants, index, rows, settings=DEFAULT_SETTINGS) -> dict`` with the
  ``entrapment.json`` blocks ``universe``, ``peptide``, ``precursor``, ``per_run`` and
  ``warnings``; ``rows`` is a DataFrame of A1 rows (``run``, ``peptide``, ``modified_peptide``,
  ``charge``, ``q_value``, ``is_decoy``, ``protein_ids``), every row (``evaluate`` applies the
  acceptance rule);
* inside the peptide block: ``sensitivity.present_broad`` (``shared``, ``k12_only``, ``r``,
  ``N_T``, ``N_E``, ``fdp_combined``, ``fdp_lower``, ``ucb95``), ``sensitivity.r_prec_weighted``
  (``fdp_combined_rprec``, ``mean_nprec_shared``, ``mean_nprec_k12_only``) and
  ``sensitivity.far`` (``r_far``, ``N_T_10``, ``N_E_far``, ``fdp_combined_far``,
  ``fdp_lower_far``, ``ucb95_far``, ``breakdown.k12_only`` and ``breakdown.N_E``, each with
  ``lt10``, ``near_homolog``, ``far``);
* a null estimate carries its reason under a key containing ``reason`` in its block;
  ``above_one`` is ``block["above_one"] is True``, in ``block["flags"]`` or a warning code.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.entrapment_m10_07 import (
    CI_ASSUMES,
    CI_METHOD,
    CONTAMINANTS,
    EXAMPLE2_ROWS,
    EXAMPLE8_ROW,
    K12,
    K12_ONLY_1,
    N_PREC_1,
    NISSLE,
    PEPTIDE_NOT_ASSESSED,
    PER_RUN_NOT_ASSESSED,
    PLAIN,
    SHARED_1,
    SPACER,
    START_PEPTIDE,
    START_REST,
    back_translate,
    codes,
    coding_text,
    flagged,
    frame,
    load_entrapment,
    origin_split,
    read_json,
    reasons,
    replicons,
    revcomp,
    row,
    run_entrapment,
    start_contig,
    write_a1,
    write_annotation_json,
    write_faa,
    write_genome,
    write_out_a1,
    write_runs_manifest,
)
from helpers.initiator_toys import cds
from helpers.readthrough_toys import model
from helpers.readthrough_toys import write_annotation as write_annotation_parquet
from helpers.readthrough_toys import write_genome as write_genome_plain
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.model.coords import Interval

TOL = 1e-12


def E() -> Any:
    return load_entrapment()


def index_of(
    contigs: Mapping[str, str],
    circular: Iterable[str] = (),
    strain_db: Any = None,
    annotation: Any = None,
) -> Any:
    """``annotation`` None: a gene model without genes, so no start is annotated (no cell)."""
    reps = replicons(contigs, circular)
    gm = model(reps, []) if annotation is None else annotation
    return E().build_index(dict(contigs), reps, strain_db, annotation=gm)


def shared_index(*texts: str, strain_db: Sequence[str] | None = ()) -> Any:
    """The shared fixture's genome (one circular contig) holding ``texts`` (default NISSLE)."""
    contigs = {"nissle": coding_text(*(texts or NISSLE))}
    return index_of(contigs, ["nissle"], list(strain_db) if strain_db is not None else None)


def evaluate(index: Any, rows: Sequence[Mapping[str, Any]] = EXAMPLE2_ROWS,
             k12: Sequence[str] = K12, **kw: Any) -> dict[str, Any]:  # fmt: skip
    out: dict[str, Any] = E().evaluate(list(k12), list(CONTAMINANTS), index, frame(rows), **kw)
    return out


def rules(index: Any, peptides: Iterable[str]) -> dict[str, str | None]:
    peps = list(peptides)
    res = E().present(peps, index)
    return {p: res[p].rule for p in peps}


def approx(x: float) -> Any:
    return pytest.approx(x, abs=TOL)


@pytest.fixture(scope="module")
def ex1() -> Any:
    return shared_index()


@pytest.fixture(scope="module")
def ex2(ex1: Any) -> dict[str, Any]:
    return evaluate(ex1)


# ------------------------------------------------------------------ example 1: digest, universe


def test_example_1_digest_sets() -> None:
    assert E().digest_set(["MSAGLLKVEPTIDEKAAGGWWR"]) == {
        "MSAGLLK", "VEPTLDEK", "AAGGWWR", "MSAGLLKVEPTLDEK", "VEPTLDEKAAGGWWR"}  # fmt: skip
    assert E().digest_set(["AGNVLTEPWR"]) == {"AGNVLTEPWR"}
    assert E().digest_set(["QQWWLLDR"]) == {"QQWWLLDR"}
    # trypsin_p cuts R|P
    assert E().digest_set(list(CONTAMINANTS)) == {"QQWWLLDR", "PEPTLDES", "QQWWLLDRPEPTLDES"}


def test_example_1_present_rules(ex1: Any) -> None:
    got = rules(ex1, [*SHARED_1, *K12_ONLY_1])
    assert got == {
        "MSAGLLK": "exact", "VEPTLDEK": "exact", "MSAGLLKVEPTLDEK": "exact",
        "AGNVLTEPWR": "mass_equivalent",  # N <-> GG block, 10 aa, budget 1
        "AAGGWWR": None,  # 7 aa: no tolerant search
        "VEPTLDEKAAGGWWR": None,  # W -> Y substitution
    }  # fmt: skip


def test_example_1_universe(ex2: dict[str, Any]) -> None:
    u = ex2["universe"]
    assert u["k12_peptides"] == 7
    assert u["excluded_contaminant"] == 1  # QQWWLLDR
    assert u["non_standard"] == 0
    assert (u["shared"], u["k12_only"]) == (4, 2)
    assert u["shared"] + u["k12_only"] == u["k12_peptides"] - u["excluded_contaminant"]
    assert (u["shared_exact"], u["shared_mass_equivalent"], u["shared_strain_db_only"]) == (3, 1, 0)
    assert u["r"] == approx(0.5)
    assert (u["shared_precursors"], u["k12_only_precursors"]) == (18, 6)
    assert u["r_prec"] == approx(1 / 3)


def test_example_1b_short_variant_is_k12_only() -> None:
    """P2 = GNVTEPWR (8 aa), translation GGGVTEPWR: rule 2 needs >= 10 aa; r = 3/3."""
    idx = shared_index("MSAGLLKVEPTLDEKAAGGYWR", "GGGVTEPWR")
    assert rules(idx, ["GNVTEPWR"]) == {"GNVTEPWR": None}
    u = evaluate(idx, k12=(K12[0], "GNVTEPWR", K12[2]))["universe"]
    assert (u["shared"], u["k12_only"]) == (3, 3)
    assert u["r"] == approx(1.0)


# ------------------------------------------------------------------ criterion 4: searchable space


@pytest.mark.parametrize(
    ("pep", "charge", "ox", "want"),
    [("MSAGLLK", 1, None, 719.41201), ("MSAGLLK", 2, None, 360.20964),
     ("MSAGLLK", 1, 0, 735.40693), ("MSAGLLK", 2, 0, 368.20710),
     ("VEPTLDEK", 1, None, 930.47783)],
)  # fmt: skip
def test_example_1_mz(pep: str, charge: int, ox: int | None, want: float) -> None:
    assert E().mz(pep, charge, oxidised_at=ox) == pytest.approx(want, abs=1e-5)


@pytest.mark.parametrize("pep", sorted(N_PREC_1))
def test_example_1_n_prec(pep: str) -> None:
    assert E().n_prec(pep) == N_PREC_1[pep]


def test_several_m_each_position_is_a_form() -> None:
    assert E().n_prec("MMGGGGGK") == 6  # 3 forms x 2 charges


def test_literal_i_does_not_change_the_mass() -> None:
    assert E().mz("VEPTIDEK", 2) == E().mz("VEPTLDEK", 2)


def test_carbamidomethyl_on_every_c() -> None:
    from gippa.mapping.masses import MOD_DELTA, RESIDUE_MASS

    pep = "CAGCGGGK"
    mass = sum(RESIDUE_MASS[c] for c in pep) + 18.010565 + 2 * MOD_DELTA["Carbamidomethyl"]
    assert E().mz(pep, 1) == pytest.approx(mass + 1.007276, abs=1e-5)
    assert E().mz(pep, 2) == pytest.approx((mass + 2 * 1.007276) / 2, abs=1e-5)


@pytest.mark.parametrize("length", range(7, 31))
@pytest.mark.parametrize("residue", ["G", "W"])
def test_peptide_level_bounds_remove_nothing_at_the_extremes(length: int, residue: str) -> None:
    assert E().n_prec(residue * length) >= 1


@settings(max_examples=300, derandomize=True, deadline=None)
@given(st.text(alphabet="ACDEFGHKLMNPQRSTVWY", min_size=7, max_size=30))
def test_every_standard_peptide_has_an_in_bounds_precursor(pep: str) -> None:
    assert E().n_prec(pep) >= 1


def test_bounds_are_applied_not_hard_coded() -> None:
    """``--max-pr-charge 1 --max-pr-mz 600`` drops VEPTLDEK (z1 930.48)."""
    s = E().Settings(charges=(1, 1), mz=(300.0, 600.0))
    assert E().n_prec("VEPTLDEK", s) == 0
    assert E().n_prec("VEPTLDEK") == 3


def test_bounds_shrink_the_universe(ex1: Any) -> None:
    """z = 1 and m/z <= 900 keep only MSAGLLK (719.41, 735.41) and AAGGWWR (803.40)."""
    s = E().Settings(charges=(1, 1), mz=(300.0, 900.0))
    assert E().n_prec("MSAGLLK", s) == 2
    u = evaluate(ex1, settings=s)["universe"]
    assert (u["shared"], u["k12_only"]) == (1, 1)
    assert (u["shared_precursors"], u["k12_only_precursors"]) == (2, 1)


def test_default_settings_are_the_a1_search() -> None:
    s = E().DEFAULT_SETTINGS
    assert (s.enzyme, s.missed, s.min_len, s.max_len, s.var_mods) == ("trypsin_p", 1, 7, 30, 1)
    assert tuple(s.charges) == (1, 4) and tuple(s.mz) == (300.0, 1800.0)


# ------------------------------------------------------------------ example 2: discoveries (S1)


def test_example_2_peptide_block(ex2: dict[str, Any]) -> None:
    p = ex2["peptide"]
    assert (p["N_T"], p["N_E"], p["excluded_contaminant"], p["not_in_universe"]) == (3, 1, 1, 0)
    assert p["fdp_combined"] == approx(0.75)
    assert p["fdp_lower"] == approx(0.25)
    assert p["assessment"] == "not_assessed"  # S1 (the spec body's evidence_of_failure is void)
    assert p["assessment_reason"] == PEPTIDE_NOT_ASSESSED
    assert p["ci"] == {"ucb95": approx(1.0), "method": CI_METHOD, "assumes": CI_ASSUMES}


def test_example_2_precursor_block_carries_the_assessment(ex2: dict[str, Any]) -> None:
    p = ex2["precursor"]
    assert (p["N_T"], p["N_E"], p["excluded_contaminant"], p["not_in_universe"]) == (4, 1, 1, 0)
    assert p["fdp_combined"] == approx(0.8)
    assert p["fdp_lower"] == approx(0.2)
    assert p["assessment"] == "evidence_of_failure"
    assert p["ci"]["ucb95"] == approx(1.0)
    assert p["ci"]["method"] == CI_METHOD


def test_example_2_per_run_is_descriptive(ex2: dict[str, Any]) -> None:
    runs = ex2["per_run"]
    assert set(runs) == {"WT_1", "WT_2", "WT_3"}
    for name, b in runs.items():
        assert b["assessment"] == "not_assessed", name
        assert b["assessment_reason"] == PER_RUN_NOT_ASSESSED, name
    w3 = runs["WT_3"]
    assert (w3["N_T"], w3["N_E"]) == (0, 1)
    assert w3["fdp_combined"] == approx(3.0)  # unclipped
    assert w3["fdp_lower"] == approx(1.0)
    assert w3["ci"]["ucb95"] == approx(1.0)
    assert flagged(w3, "above_one", ex2["warnings"])
    w1 = runs["WT_1"]  # rows 1, 3, 7 accepted, row 5 (q 0.02) not
    assert (w1["N_T"], w1["N_E"]) == (3, 0)
    assert w1["fdp_combined"] == approx(0.0)
    w2 = runs["WT_2"]  # row 2 (oxidised MSAGLLK), row 6 contaminant, row 8 decoy
    assert (w2["N_T"], w2["N_E"], w2["excluded_contaminant"]) == (1, 0, 1)


def test_example_2_q_cutoff_inclusive_and_na_not_accepted(ex1: Any) -> None:
    rows = [*EXAMPLE2_ROWS[:3], row("WT_3", "AAGGWWR", 2, None), *EXAMPLE2_ROWS[4:]]
    p = evaluate(ex1, rows)["peptide"]
    assert (p["N_T"], p["N_E"]) == (3, 0)  # AGNVLTEPWR at q = 0.01 is in; NA q is out


def test_example_2_no_warning_but_above_one(ex2: dict[str, Any]) -> None:
    assert "digest_mismatch" not in codes(ex2["warnings"])
    assert "n_entrapment_zero" not in codes(ex2["warnings"])


def test_decoys_are_never_discoveries(ex1: Any) -> None:
    rows = [*EXAMPLE2_ROWS, row("WT_1", "AAGGWWR", 2, 0.001, True, ids="rev_x")]
    rows += [row("WT_1", "VEPTLDEKAAGGWWR", 2, 0.001, True, ids="rev_x")]
    p = evaluate(ex1, rows)["peptide"]
    assert (p["N_T"], p["N_E"], p["not_in_universe"]) == (3, 1, 0)


# ------------------------------------------------------------------ S4 sensitivities (example 2)


def test_example_2_r_prec_weighted(ex2: dict[str, Any]) -> None:
    s = ex2["peptide"]["sensitivity"]["r_prec_weighted"]
    assert s["fdp_combined_rprec"] == approx(1.0)  # 1 * (1 + 3) / 4
    assert s["mean_nprec_shared"] == approx(4.5)
    assert s["mean_nprec_k12_only"] == approx(3.0)
    assert "assessment" not in s


def test_example_2_far_stratum(ex2: dict[str, Any]) -> None:
    f = ex2["peptide"]["sensitivity"]["far"]
    assert f["breakdown"]["k12_only"] == {"lt10": 1, "near_homolog": 1, "far": 0}
    assert f["breakdown"]["N_E"] == {"lt10": 1, "near_homolog": 0, "far": 0}
    assert f["N_E_far"] == 0
    assert f["N_T_10"] == 1  # AGNVLTEPWR
    assert f["r_far"] == approx(0.0)  # |E_far| = 0 of |T_10| = 2
    assert f["fdp_combined_far"] is None
    assert "no_entrapment_space" in reasons(f)
    assert f["fdp_lower_far"] == approx(0.0)
    assert f["ucb95_far"] is None
    assert "assessment" not in f


def test_far_stratum_counts_a_k12_only_peptide_without_any_tolerant_row(ex1: Any) -> None:
    """A >= 10-aa K-12-only peptide with no Nissle neighbour is ``far``; discovered, it is
    N_E_far and makes r_far = 1/2."""
    far = "HHFFHHFFYYDR"  # 12 aa, nothing like it in the translation
    rows = [*EXAMPLE2_ROWS, row("WT_1", far, 2, 0.002)]
    out = evaluate(ex1, rows, k12=(*K12, far))
    f = out["peptide"]["sensitivity"]["far"]
    assert f["breakdown"]["k12_only"] == {"lt10": 1, "near_homolog": 1, "far": 1}
    assert f["breakdown"]["N_E"] == {"lt10": 1, "near_homolog": 0, "far": 1}
    assert f["N_E_far"] == 1 and f["N_T_10"] == 1
    assert f["r_far"] == approx(0.5)
    assert f["fdp_combined_far"] == approx(1 * 3 / 2)  # 1 (1 + 2) / (1 + 1)
    assert f["fdp_lower_far"] == approx(0.5)
    assert f["ucb95_far"] == approx(E().ucb95(1, 1, 1, 2))
    assert out["peptide"]["N_E"] == 2


def test_example_2_present_broad_equals_primary_on_the_shared_fixture(ex2: dict[str, Any]) -> None:
    b = ex2["peptide"]["sensitivity"]["present_broad"]
    assert (b["shared"], b["k12_only"], b["N_T"], b["N_E"]) == (4, 2, 3, 1)
    assert b["r"] == approx(0.5)
    assert b["fdp_combined"] == approx(0.75) and b["fdp_lower"] == approx(0.25)
    assert b["ucb95"] == approx(1.0)
    assert "assessment" not in b


def test_example_1c_swap_primary_and_present_broad() -> None:
    """S2 example 1c: translation AGNVLTEPRW. Primary: AGNVLTEPWR (swap) is k12_only, r 1.0;
    present_broad: shared, r 0.5."""
    idx = shared_index("MSAGLLKVEPTLDEKAAGGYWR", "AGNVLTEPRW")
    assert rules(idx, ["AGNVLTEPWR"]) == {"AGNVLTEPWR": None}
    out = evaluate(idx)
    u, p = out["universe"], out["peptide"]
    assert (u["shared"], u["k12_only"]) == (3, 3)
    assert u["r"] == approx(1.0)
    assert (p["N_T"], p["N_E"]) == (2, 2)
    assert p["fdp_combined"] == approx(1.0)  # 2 * 6 / (3 * 4)
    b = p["sensitivity"]["present_broad"]
    assert (b["shared"], b["k12_only"]) == (4, 2)
    assert b["r"] == approx(0.5)
    assert (b["N_T"], b["N_E"]) == (3, 1)
    assert b["fdp_combined"] == approx(0.75)
    f = p["sensitivity"]["far"]
    assert f["breakdown"]["k12_only"] == {"lt10": 1, "near_homolog": 2, "far": 0}


def test_near_block_is_present_only_in_present_broad() -> None:
    """AGWVLTEPWR vs the translation AGEGVLTEPWR: W <-> EG is a near_block (0.015 Da)."""
    idx = shared_index(*NISSLE, "AGEGVLTEPWR")
    assert rules(idx, ["AGWVLTEPWR"]) == {"AGWVLTEPWR": None}
    rows = [*EXAMPLE2_ROWS, row("WT_1", "AGWVLTEPWR", 2, 0.002)]
    p = evaluate(idx, rows, k12=(*K12, "AGWVLTEPWR"))["peptide"]
    assert (p["N_T"], p["N_E"]) == (3, 2)
    b = p["sensitivity"]["present_broad"]
    assert (b["shared"], b["k12_only"], b["N_T"], b["N_E"]) == (5, 2, 4, 1)


# ------------------------------------------------------------------ criterion 6: present rules


def test_block_row_makes_present_other_types_do_not() -> None:
    """block (N <-> GG) yes; substitution, swap, near_block and frameshift no."""
    fs = SPACER + back_translate("LLDEGT") + "A" + back_translate("SWAVPR") + SPACER
    contigs = {
        "c": coding_text("AGGGVLTEPWR", "SAGEGVLTEPWR", "AGNVLTEPRW", "VEPTLDEKAAGGYWR") + fs,
    }
    idx = index_of(contigs, ["c"], [])
    assert rules(idx, ["AGNVLTEPWR", "SAGWVLTEPWR", "VEPTLDEKAAGGWWR", "LLDEGTSWAVPR"]) == {
        "AGNVLTEPWR": "mass_equivalent",  # a block row (a swap row at another locus too)
        "SAGWVLTEPWR": None,  # near_block
        "VEPTLDEKAAGGWWR": None,  # substitution
        "LLDEGTSWAVPR": None,  # frameshift only
    }


def test_swap_only_is_not_present() -> None:
    idx = index_of({"c": coding_text("AGNVLTEPRW")}, ["c"], [])
    assert rules(idx, ["AGNVLTEPWR"]) == {"AGNVLTEPWR": None}


def test_origin_crossing_on_a_circular_contig_is_present() -> None:
    contigs = {"chr": coding_text("MSAGLLKVEPTLDEKAAGGYWR"), "p1": origin_split("AGGGVLTEPWR", 5)}
    circ = index_of(contigs, ["chr", "p1"], [])
    assert rules(circ, ["AGGGVLTEPWR", "AGNVLTEPWR"]) == {
        "AGGGVLTEPWR": "exact", "AGNVLTEPWR": "mass_equivalent"}  # fmt: skip
    assert evaluate(circ)["universe"]["shared"] == 4


def test_the_same_split_on_a_linear_contig_is_absent() -> None:
    contigs = {"chr": coding_text("MSAGLLKVEPTLDEKAAGGYWR"), "p1": origin_split("AGGGVLTEPWR", 5)}
    lin = index_of(contigs, ["chr"], [])
    assert rules(lin, ["AGGGVLTEPWR", "AGNVLTEPWR"]) == {"AGGGVLTEPWR": None, "AGNVLTEPWR": None}
    u = evaluate(lin)["universe"]
    assert (u["shared"], u["k12_only"]) == (3, 3)


def test_minus_strand_gives_the_same_classes(ex2: dict[str, Any]) -> None:
    idx = index_of({"nissle": revcomp(coding_text(*NISSLE))}, ["nissle"], [])
    assert rules(idx, [*SHARED_1, *K12_ONLY_1]) == rules(shared_index(), [*SHARED_1, *K12_ONLY_1])
    out = evaluate(idx)
    assert out["universe"] == ex2["universe"]
    assert out["peptide"] == ex2["peptide"]


def test_i_l_k12_peptide_with_l_in_the_genome_is_shared(ex1: Any) -> None:
    assert rules(ex1, ["VEPTIDEK"]) == {"VEPTIDEK": "exact"}


def test_present_rule_order_exact_before_strain_db() -> None:
    idx = shared_index(strain_db=["MSAGLLKVEPTIDEKAAGGYWR"])
    assert rules(idx, ["MSAGLLK"]) == {"MSAGLLK": "exact"}


# ------------------------------------------------------------------ example 8: digest mismatch


def test_example_8_digest_mismatch(ex1: Any, ex2: dict[str, Any]) -> None:
    out = evaluate(ex1, [*EXAMPLE2_ROWS, EXAMPLE8_ROW])
    p = out["peptide"]
    assert p["not_in_universe"] == 1
    assert (p["N_T"], p["N_E"]) == (ex2["peptide"]["N_T"], ex2["peptide"]["N_E"])
    assert p["fdp_combined"] == approx(0.75)
    assert "digest_mismatch" in codes(out["warnings"])
    w = next(
        w for w in out["warnings"] if not isinstance(w, str) and w.get("code") == "digest_mismatch"
    )
    assert 1 in w.values()  # the count


def test_contaminant_labelled_row_outside_d_c_is_not_in_universe(ex1: Any) -> None:
    rows = [*EXAMPLE2_ROWS, row("WT_1", "GGGGSGGGGK", 2, 0.001, ids="Cont_P99999")]
    out = evaluate(ex1, rows)
    assert out["peptide"]["not_in_universe"] == 1
    assert out["peptide"]["excluded_contaminant"] == 1
    assert "digest_mismatch" in codes(out["warnings"])


# ------------------------------------------------------------------ example 9: non-ATG start (L2)

P4 = "MTEKLLAGGPR"


def test_example_9_strain_db_rule() -> None:
    idx = shared_index(*NISSLE, "VTEKLLAGGPR", strain_db=["MTEKLLAGGPRAAGK"])
    assert rules(idx, ["LLAGGPR", P4]) == {"LLAGGPR": "exact", P4: "strain_db"}
    u = evaluate(idx, k12=(*K12, P4))["universe"]
    assert (u["shared"], u["k12_only"]) == (6, 2)
    assert (u["shared_exact"], u["shared_strain_db_only"], u["shared_mass_equivalent"]) == (4, 1, 1)


def test_example_9_without_strain_db_is_k12_only() -> None:
    """L2 answered "no": M <-> V is a substitution for rule 2."""
    idx = shared_index(*NISSLE, "VTEKLLAGGPR", strain_db=None)
    assert rules(idx, ["LLAGGPR", P4]) == {"LLAGGPR": "exact", P4: None}
    u = evaluate(idx, k12=(*K12, P4))["universe"]
    assert (u["shared"], u["k12_only"], u["shared_strain_db_only"]) == (5, 3, 0)


def test_strain_db_matches_after_i_to_l() -> None:
    idx = shared_index(*NISSLE, strain_db=["MTEKIIAGGPRAAGK"])
    assert rules(idx, [P4]) == {P4: "strain_db"}


# ------------------------------------------------------------------ edge cases


def test_two_k12_peptides_differing_by_i_l_are_one_element(ex2: dict[str, Any], ex1: Any) -> None:
    u = evaluate(ex1, k12=(*K12, "VEPTLDEK"))["universe"]
    assert u == ex2["universe"]
    assert (u["shared"], u["k12_only"]) == (4, 2)


def test_two_discoveries_differing_by_i_l_are_one(ex1: Any) -> None:
    rows = [*EXAMPLE2_ROWS, row("WT_2", "VEPTLDEK", 2, 0.002), row("WT_3", "VEPTIDEK", 3, 0.002)]
    out = evaluate(ex1, rows)
    assert out["peptide"]["N_T"] == 3
    assert out["precursor"]["N_T"] == 5  # (VEPTLDEK, 2) once, (VEPTLDEK, 3) new


def test_one_peptide_in_several_runs_is_one_peptide(ex1: Any) -> None:
    rows = [*EXAMPLE2_ROWS, row("WT_1", "AAGGWWR", 2, 0.001), row("WT_2", "AAGGWWR", 3, 0.004)]
    out = evaluate(ex1, rows)
    assert out["peptide"]["N_E"] == 1
    assert out["precursor"]["N_E"] == 2


def test_non_standard_products_are_dropped(ex2: dict[str, Any], ex1: Any) -> None:
    """GGUGGGGK and GGUGGGGKAAGGWWR carry U; GGXGGGGR carries X."""
    u = evaluate(ex1, k12=(*K12, "GGUGGGGKAAGGWWR", "GGXGGGGR"))["universe"]
    assert u["non_standard"] == 3
    assert (u["shared"], u["k12_only"]) == (4, 2)


def test_empty_discoveries(ex1: Any) -> None:
    rows = [r | {"q_value": 0.5} for r in EXAMPLE2_ROWS]
    out = evaluate(ex1, rows)
    for level in ("peptide", "precursor"):
        b = out[level]
        assert (b["N_T"], b["N_E"]) == (0, 0), level
        assert b["fdp_combined"] is None and b["fdp_lower"] is None, level
        assert b["ci"]["ucb95"] is None, level
        assert "empty_discoveries" in reasons(b), level
    assert out["precursor"]["assessment"] == "not_evaluable"
    assert out["peptide"]["assessment"] == "not_assessed"
    # rev-1007 S7 (mutant M21): no accepted row is not "no K-12-only peptide was accepted"
    assert "n_entrapment_zero" not in codes(out["warnings"])


def test_no_shared_peptides(ex1: Any) -> None:
    idx = index_of({"c": coding_text("HHHHHHHHHH")}, ["c"], [])
    out = evaluate(idx)
    assert out["universe"]["shared"] == 0
    assert out["universe"]["r"] is None
    for level in ("peptide", "precursor"):
        b = out[level]
        assert b["fdp_combined"] is None and b["fdp_lower"] is None, level
        assert b["ci"]["ucb95"] is None, level
        assert "no_shared_peptides" in reasons(b), level


def test_no_entrapment_space(ex1: Any) -> None:
    """K-12 = one protein whose products are all present: r = 0, 1/r undefined."""
    rows = [row("WT_1", "MSAGLLK", 2, 0.001), row("WT_2", "VEPTIDEK", 2, 0.001)]
    out = evaluate(ex1, rows, k12=("MSAGLLKVEPTIDEK",))
    assert out["universe"]["k12_only"] == 0
    assert out["universe"]["r"] == approx(0.0)
    p = out["peptide"]
    assert (p["N_T"], p["N_E"]) == (2, 0)
    assert p["fdp_combined"] is None and p["ci"]["ucb95"] is None
    assert "no_entrapment_space" in reasons(p)
    assert p["fdp_lower"] == approx(0.0)


def test_example_7_shape_n_entrapment_zero_warning(ex1: Any) -> None:
    rows = [r for r in EXAMPLE2_ROWS if r["peptide"] != "AAGGWWR"]
    out = evaluate(ex1, rows)
    p = out["precursor"]
    assert p["N_E"] == 0
    assert p["fdp_combined"] == approx(0.0) and p["fdp_lower"] == approx(0.0)
    assert p["assessment"] == "consistent_with_control"
    assert p["ci"]["ucb95"] == approx(E().ucb95(4, 0, 6, 18))
    assert "n_entrapment_zero" in codes(out["warnings"])


# ------------------------------------------------------------------ criterion 2: one function


def test_one_classification_function_for_universe_and_discoveries(
    ex1: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    mod = E()
    real = mod.present
    calls: list[tuple[list[str], Any]] = []

    def spy(peptides: Iterable[str], index: Any) -> Any:
        peps = list(peptides)
        calls.append((peps, index))
        return real(peps, index)

    monkeypatch.setattr(mod, "present", spy)
    out = evaluate(ex1)
    assert calls, "evaluate never called the module's present()"
    assert all(idx is ex1 for _, idx in calls), "present() was called with another index"
    seen = {mod.norm(p) for peps, _ in calls for p in peps}
    assert seen >= SHARED_1 | K12_ONLY_1  # the universe
    assert {"MSAGLLK", "VEPTLDEK", "AGNVLTEPWR", "AAGGWWR"} <= seen  # the discoveries in U
    assert out["peptide"]["N_T"] == 3


def test_discovery_classes_equal_universe_classes(ex1: Any) -> None:
    """``present`` on a discovery spelled as DIA-NN wrote it (literal I) gives the class of its
    universe element."""
    lit = rules(ex1, ["VEPTIDEK", "MSAGLLK", "AGNVLTEPWR", "AAGGWWR"])
    uni = rules(ex1, ["VEPTLDEK", "MSAGLLK", "AGNVLTEPWR", "AAGGWWR"])
    assert list(lit.values()) == list(uni.values())


# ------------------------------------------------------------------ rev-1007 S7/N6: mutation pins


def _tolerant_types(peptide: str, contigs: Mapping[str, str]) -> set[frozenset[str]]:
    """The ``diff_types`` of GIPPA's tolerant target rows for ``peptide`` (fixture guard)."""
    from gippa.config import resolve_config
    from gippa.mapping.sixframe import build_sixframe
    from gippa.mapping.tolerant import match_tolerant

    six = build_sixframe(dict(contigs), replicons(contigs, list(contigs)))
    cfg = resolve_config().mapping.tolerant
    res = match_tolerant([peptide], six, {"shuffle": ()}, {peptide: "shuffle"}, {peptide}, cfg)
    tgt = res[res["search"] == "target"]
    return {frozenset(x for x in str(d).split(";") if x) for d in tgt["diff_types"]}


def test_mixed_block_and_substitution_row_is_not_present() -> None:
    """rev-1007 S7 (M06): rule 2 needs ``diff_types`` non-empty and only ``block``. SAGNVLTEPWWR
    (12 aa, budget 2) against SAGGGVLTEPWYR is N <-> GG plus W -> Y: one ``block;substitution``
    row, so the peptide stays ``k12_only`` (also in ``present_broad``: substitution is not in
    {block, swap, near_block})."""
    pep, text = "SAGNVLTEPWWR", "SAGGGVLTEPWYR"
    contigs = {"nissle": coding_text(*NISSLE, text)}
    assert _tolerant_types(pep, contigs) == {frozenset({"block", "substitution"})}
    idx = index_of(contigs, ["nissle"], [])
    assert rules(idx, [pep]) == {pep: None}
    rows = [*EXAMPLE2_ROWS, row("WT_1", pep, 2, 0.002)]
    p = evaluate(idx, rows, k12=(*K12, pep))["peptide"]
    assert (p["N_T"], p["N_E"]) == (3, 2)
    b = p["sensitivity"]["present_broad"]
    assert (b["N_T"], b["N_E"]) == (3, 2)


def test_above_one_is_strictly_above_one(ex1: Any) -> None:
    """rev-1007 N6 (M15): WT_3 with N_T 2, N_E 1 and r 0.5 gives fdp_combined = 1 * 3 / 3 = 1.0
    exactly, which is not ``above_one``."""
    rows = [*EXAMPLE2_ROWS, row("WT_3", "MSAGLLK", 2, 0.001), row("WT_3", "VEPTIDEK", 2, 0.001)]
    out = evaluate(ex1, rows)
    w3 = out["per_run"]["WT_3"]
    assert (w3["N_T"], w3["N_E"]) == (2, 1)
    assert w3["fdp_combined"] == 1.0
    assert w3.get("above_one") is not True
    assert "above_one" not in (w3.get("flags") or [])


def test_na_is_decoy_is_a_target_row_as_in_m10_05(ex1: Any) -> None:
    """rev-1007 N6 (M16): M10-05's accepted-row rule reads a missing ``is_decoy`` as target
    (``metrics._target``: ``fillna(False)``); the entrapment uses the same rule."""
    rows = [r | {"is_decoy": None} if r["peptide"] == "AAGGWWR" else r for r in EXAMPLE2_ROWS]
    p = evaluate(ex1, rows)["peptide"]
    assert (p["N_T"], p["N_E"]) == (3, 1)


def test_a_run_with_no_accepted_row_has_an_empty_per_run_block(ex1: Any) -> None:
    """rev-1007 N6 (M28): every run of ``PEP_A1`` gets a block; WT_3's only row is rejected
    (q 0.5), so its block is the empty-discoveries one."""
    rows = [r | {"q_value": 0.5} if r["run"] == "WT_3" else r for r in EXAMPLE2_ROWS]
    out = evaluate(ex1, rows)
    assert set(out["per_run"]) == {"WT_1", "WT_2", "WT_3"}
    w3 = out["per_run"]["WT_3"]
    assert (w3["N_T"], w3["N_E"]) == (0, 0)
    assert w3["fdp_combined"] is None and w3["fdp_lower"] is None
    assert "empty_discoveries" in reasons(w3)
    assert w3["assessment"] == "not_assessed"


def test_no_shared_peptides_takes_precedence_over_empty_discoveries() -> None:
    """rev-1007 N6 (M32): with both |shared| = 0 and no accepted row, the reason is
    ``no_shared_peptides`` (r itself is undefined, whatever the discoveries)."""
    idx = index_of({"c": coding_text("HHHHHHHHHH")}, ["c"], [])
    rows = [r | {"q_value": 0.5} for r in EXAMPLE2_ROWS]
    out = evaluate(idx, rows)
    for level in ("peptide", "precursor"):
        b = out[level]
        assert b["fdp_combined"] is None and b["fdp_lower"] is None, level
        assert "no_shared_peptides" in reasons(b), level
        assert "empty_discoveries" not in reasons(b), level


# ------------------------------------------------------------------ ADR-0068: initiator M (M11-20)
# Since ADR-0068 GIPPA's exact route places a peptide whose first residue is an M read at an
# annotated non-ATG start (flag ``initiator_m``). ``present`` is the same matcher, so it must credit
# the same cells, and cross-check 1 stays as strict as before: a credited peptide agrees because it
# IS exact here, never because the guard looks away. All sequences are synthetic.

CODON_STRAND = [(c, s) for c in PLAIN for s in "+-"]
STARTS = ("ATG", "GTG", "TTG")
OTHER_REST = "AGGHHFFDK"
OTHER_PEPTIDE = "M" + OTHER_REST
#: ``cross_checks`` fixture rows: example 2 plus the start peptide's accepted precursor
START_ROWS = (*EXAMPLE2_ROWS, row("WT_1", START_PEPTIDE, 2, 0.002))


def _flip(gene: Any, n: int) -> Any:
    iv = gene.interval
    strand = "-" if iv.strand == "+" else "+"
    flipped = Interval(iv.contig, n - iv.end, n - iv.start, strand)  # type: ignore[arg-type]
    return dataclasses.replace(gene, interval=flipped)


def on_strand(text: str, genes: Sequence[Any], strand: str) -> tuple[dict[str, str], list[Any]]:
    """``{"nissle": text}`` and ``genes`` as given, or the reverse complement of both for ``-``."""
    if strand == "+":
        return {"nissle": text}, list(genes)
    return {"nissle": revcomp(text)}, [_flip(g, len(text)) for g in genes]


def start_world(
    codon: str, strand: str = "+", *, aa: str = START_PEPTIDE, feature_type: str = "cds"
) -> tuple[dict[str, str], list[Any]]:
    """The fixture's genes plus the start peptide's gene, annotated with a first codon ``codon``."""
    text, a, b = start_contig(codon)
    return on_strand(text, [cds("NISSLE_S", "nissle", a, b, "+", aa, codon, feature_type)], strand)


def gm_of(contigs: Mapping[str, str], genes: Sequence[Any]) -> Any:
    return model(replicons(contigs, ["nissle"]), genes)


def world_index(contigs: Mapping[str, str], genes: Sequence[Any], strain_db: Any = ()) -> Any:
    return index_of(contigs, ["nissle"], list(strain_db), gm_of(contigs, genes))


def presences(index: Any, peptides: Iterable[str]) -> dict[str, Any]:
    return dict(E().present(list(peptides), index))


def unannotated_world(case: str, codon: str, strand: str) -> tuple[dict[str, str], list[Any]]:
    """The start peptide's text at ``codon``, in a genome that does not make that codon a start."""
    if case == "internal_codon":  # in frame inside an annotated ATG gene
        text, a, b = start_contig(codon, upstream="MAA")
        aa = "MAA" + PLAIN[codon] + START_REST
        return on_strand(text, [cds("G", "nissle", a - 9, b, "+", aa, "ATG")], strand)
    text, a, b = start_contig(codon)
    if case == "no_gene":
        return on_strand(text, [], strand)
    if case == "pseudogene":
        return start_world(codon, strand, feature_type="pseudogene")
    if case == "protein_not_m":  # Bakta's protein starts with the plain residue (rule 1)
        return start_world(codon, strand, aa=PLAIN[codon] + START_REST)
    if case == "gene_on_the_other_strand":  # same interval: its first codon is the TAA's TTA
        return on_strand(text, [cds("G", "nissle", a, b, "-", START_PEPTIDE, "TTA")], strand)
    assert case == "other_start_only"  # the index has cells, none at this codon
    other = "GTG" + back_translate(OTHER_REST) + "TAA"
    gene = cds("G", "nissle", len(text), len(text) + len(other), "+", OTHER_PEPTIDE, "GTG")
    return on_strand(text + other + SPACER, [gene], strand)


UNANNOTATED = ("no_gene", "other_start_only", "internal_codon", "pseudogene", "protein_not_m",
               "gene_on_the_other_strand")  # fmt: skip


# item 1: annotated GTG and TTG starts, both strands


@pytest.mark.parametrize(("codon", "strand"), CODON_STRAND)
def test_m_at_an_annotated_non_atg_start_is_exact(codon: str, strand: str) -> None:
    idx = world_index(*start_world(codon, strand))
    got = presences(idx, [START_PEPTIDE])[START_PEPTIDE]
    assert got.rule == "exact"
    assert got.tolerant == ()  # rule 1 won: no tolerant search for it


def test_a_short_prefix_and_literal_i_are_credited_too() -> None:
    """Every length of window is credited at the cell, and I -> L still applies."""
    idx = world_index(*start_world("GTG"))
    assert rules(idx, ["MTEKLLA", "MTEKIIAGGPR"]) == {"MTEKLLA": "exact", "MTEKIIAGGPR": "exact"}


def test_exact_comes_before_strain_db_at_an_annotated_start() -> None:
    """Example 9's P4 was ``strain_db`` only; an annotated start now makes it ``exact``, which is
    what cross-check 1 needs (the A1 freeze run fired on exactly such peptides)."""
    contigs, genes = start_world("GTG")
    idx = world_index(contigs, genes, ["MTEKLLAGGPRAAGK"])
    assert rules(idx, ["LLAGGPR", P4]) == {"LLAGGPR": "exact", P4: "exact"}
    u = evaluate(idx, k12=(*K12, P4))["universe"]
    assert (u["shared"], u["k12_only"]) == (6, 2)
    assert (u["shared_exact"], u["shared_strain_db_only"], u["shared_mass_equivalent"]) == (5, 0, 1)


def test_an_annotated_start_makes_a_discovery_a_target_not_an_entrapment() -> None:
    """The estimate itself changes: example 9 without ``strain_db`` counts P4 as entrapment."""
    contigs, genes = start_world("GTG")
    plain = evaluate(world_index(contigs, []), START_ROWS, k12=(*K12, P4))
    credited = evaluate(world_index(contigs, genes), START_ROWS, k12=(*K12, P4))
    assert (plain["universe"]["shared"], plain["universe"]["k12_only"]) == (5, 3)
    assert (plain["peptide"]["N_T"], plain["peptide"]["N_E"]) == (3, 2)
    assert (credited["universe"]["shared"], credited["universe"]["k12_only"]) == (6, 2)
    assert (credited["peptide"]["N_T"], credited["peptide"]["N_E"]) == (4, 1)
    assert (P4, "target", "exact") in {tuple(d) for d in credited["discoveries"]}


@pytest.mark.parametrize(("codon", "strand"), CODON_STRAND)
def test_index_holds_the_start_cells_gippa_builds(codon: str, strand: str) -> None:
    from gippa.mapping.initiator import StartCells, annotated_starts, start_frame_cells

    contigs, genes = start_world(codon, strand)
    idx = world_index(contigs, genes)
    sites, _ = annotated_starts(gm_of(contigs, genes), contigs, STARTS, table=11)
    assert len(sites) == 1 and sites[0].codon == codon and sites[0].strand == strand
    assert isinstance(idx.starts, StartCells)
    assert idx.starts.sites == sites
    want = start_frame_cells(idx.sixframe, sites)
    assert [c.tolist() for c in idx.starts.target] == [c.tolist() for c in want]


def test_present_passes_the_index_start_cells(monkeypatch: pytest.MonkeyPatch) -> None:
    mod = E()
    idx = world_index(*start_world("GTG"))
    real, seen = mod.match_exact, []

    def spy(*args: Any, **kw: Any) -> Any:
        seen.append(kw.get("starts"))
        return real(*args, **kw)

    monkeypatch.setattr(mod, "match_exact", spy)
    mod.present([START_PEPTIDE], idx)
    assert seen and all(s is idx.starts for s in seen)


def test_origin_crossing_annotated_start_on_a_circular_contig() -> None:
    """The start codon just before the origin: the cell is in the origin-extended frame."""
    nt = "GTG" + back_translate(START_REST) + "TAA"
    pad = SPACER * 4
    text = nt[6:] + pad + nt[:6]
    contigs = {"nissle": text}
    a = len(text) - 6
    idx = world_index(contigs, [cds("G", "nissle", a, a + len(nt), "+", START_PEPTIDE, "GTG")])
    assert rules(idx, [START_PEPTIDE]) == {START_PEPTIDE: "exact"}


# item 2: the same sequence, unannotated


@pytest.mark.parametrize("case", UNANNOTATED)
@pytest.mark.parametrize(("codon", "strand"), CODON_STRAND)
def test_m_at_an_unannotated_non_atg_codon_is_not_exact(case: str, codon: str, strand: str) -> None:
    contigs, genes = unannotated_world(case, codon, strand)
    idx = world_index(contigs, genes)
    assert rules(idx, [START_PEPTIDE]) == {START_PEPTIDE: None}
    if case == "other_start_only":  # the fixture's cells are real: its own peptide is credited
        assert len(idx.starts.sites) == 1
        assert rules(idx, [OTHER_PEPTIDE]) == {OTHER_PEPTIDE: "exact"}
    else:
        assert not idx.starts.sites


def test_m_is_credited_only_as_the_first_residue() -> None:
    """A query with a residue before the M reads the codon as V at position 1: no credit."""
    text, a, b = start_contig("GTG", upstream="S")
    genes = [cds("G", "nissle", a, b, "+", START_PEPTIDE, "GTG")]
    idx = world_index({"nissle": text}, genes)
    assert rules(idx, ["S" + START_PEPTIDE, "S" + PLAIN["GTG"] + START_REST]) == {
        "S" + START_PEPTIDE: None,
        "S" + PLAIN["GTG"] + START_REST: "exact",
    }
    assert rules(idx, [START_PEPTIDE]) == {START_PEPTIDE: "exact"}


# item 3: the plain translation at that start behaves as today


@pytest.mark.parametrize(("codon", "strand"), CODON_STRAND)
def test_plain_translation_at_the_start_is_unchanged(codon: str, strand: str) -> None:
    contigs, genes = start_world(codon, strand)
    plain = PLAIN[codon] + START_REST
    other_residue = "A" + START_REST  # neither M nor what the codon reads: a substitution
    peptides = [plain, other_residue]
    with_cells = presences(world_index(contigs, genes), peptides)
    assert with_cells == presences(world_index(contigs, []), peptides)
    assert with_cells[plain].rule == "exact"
    assert with_cells[other_residue].rule is None
    assert with_cells[other_residue].tolerant  # its tolerant rows are still searched


def test_literal_i_for_the_ttg_plain_translation_is_unchanged() -> None:
    contigs, genes = start_world("TTG")
    assert rules(world_index(contigs, genes), ["ITEKLLAGGPR"]) == {"ITEKLLAGGPR": "exact"}


# item 4: nothing else changes when start cells exist and no peptide uses them


def test_unused_start_cells_change_no_block_of_the_output() -> None:
    contigs, genes = start_world("GTG")
    used = world_index(contigs, genes)
    assert used.starts.sites, "the fixture must have a start cell"
    got, want = evaluate(used), evaluate(world_index(contigs, []))
    assert got == want
    assert (got["universe"]["shared"], got["universe"]["k12_only"]) == (4, 2)
    assert (got["peptide"]["N_T"], got["peptide"]["N_E"]) == (3, 1)


def test_unused_start_cells_change_no_presence_of_the_universe() -> None:
    contigs, genes = start_world("TTG", "-")
    peptides = [*SHARED_1, *K12_ONLY_1, "VEPTIDEK", "GGGVTEPWR", "LLAGGPR"]
    assert presences(world_index(contigs, genes), peptides) == presences(
        world_index(contigs, []), peptides
    )


# the config values GIPPA's mapping step uses, read from the one place


def _config_with(**annotation: Any) -> Any:
    from gippa.config import resolve_config

    cfg = resolve_config()
    return cfg.model_copy(update={"annotation": cfg.annotation.model_copy(update=annotation)})


@pytest.mark.parametrize(
    ("codons", "want"),
    [(["ATG", "GTG"], {"GTG": "exact", "TTG": None}),
     (["ATG", "TTG"], {"GTG": None, "TTG": "exact"}),
     (["atg", "gtg", "ttg"], {"GTG": "exact", "TTG": "exact"})],
)  # fmt: skip
def test_start_codons_come_from_the_config(
    monkeypatch: pytest.MonkeyPatch, codons: list[str], want: dict[str, str | None]
) -> None:
    cfg = _config_with(start_codons=codons)
    monkeypatch.setattr(E(), "resolve_config", lambda *a, **k: cfg)
    for codon, rule in want.items():
        idx = world_index(*start_world(codon))
        assert rules(idx, [START_PEPTIDE]) == {START_PEPTIDE: rule}, codon


def test_translation_table_comes_from_the_config(monkeypatch: pytest.MonkeyPatch) -> None:
    """Table 4 reads TGA as W in the index the cells live in, as GIPPA's ``sixframe_index`` does."""
    text, a, b = start_contig("GTG")
    recoded = back_translate("AGGGK") + "TGA" + back_translate("HHDR") + SPACER
    contigs = {"nissle": text + recoded}
    genes = [cds("G", "nissle", a, b, "+", START_PEPTIDE, "GTG")]
    assert rules(world_index(contigs, genes), ["AGGGKWHHDR"]) == {"AGGGKWHHDR": None}  # table 11
    cfg = _config_with(translation_table=4)
    monkeypatch.setattr(E(), "resolve_config", lambda *a, **k: cfg)
    idx = world_index(contigs, genes)
    assert idx.sixframe.translation_table == 4
    assert rules(idx, ["AGGGKWHHDR", START_PEPTIDE]) == {
        "AGGGKWHHDR": "exact",
        START_PEPTIDE: "exact",
    }


# item 6: parity with GIPPA's own route


def _parity_world() -> tuple[dict[str, str], dict[str, Any], Any, set[str], set[str]]:
    """Contigs, replicons, gene model, the peptides exact and the peptides not exact, by hand.

    ``c1`` (circular): annotated GTG, annotated TTG, an ATG gene, an unannotated GTG and a TTG
    inside an ATG gene. ``c2`` (linear): annotated GTG and TTG on the minus strand. ``c3``
    (circular): an annotated GTG just before the origin.
    """
    parts = [("GTG", START_REST), ("TTG", "AGDEFHHNK"), ("ATG", "SAGLLKVEPTLDEK"),
             ("GTG", "WQQDEHHK"), ("ATG", "AATYYFFDEK")]  # fmt: skip
    c1, spans = SPACER, []
    for codon, rest in parts:
        nt = codon + back_translate(rest) + "TAA"
        spans.append((len(c1), len(c1) + len(nt)))
        c1 += nt + SPACER
    # the T of AATYYFFDEK's third residue is replaced by a TTG inside the ATG gene
    internal = spans[4][0] + 9
    c1 = c1[:internal] + "TTG" + c1[internal + 3 :]
    genes = [cds("A", "c1", *spans[0], "+", START_PEPTIDE, "GTG"),
             cds("B", "c1", *spans[1], "+", "MAGDEFHHNK", "TTG"),
             cds("C", "c1", *spans[2], "+", "MSAGLLKVEPTLDEK", "ATG"),
             cds("D", "c1", *spans[4], "+", "MAALYYFFDEK", "ATG")]  # fmt: skip
    minus, mspans = SPACER, []
    for codon, rest in (("GTG", "NNPPFFWR"), ("TTG", "HHDDYYQR")):
        nt = codon + back_translate(rest) + "TAA"
        mspans.append((len(minus), len(minus) + len(nt)))
        minus += nt + SPACER
    c2, n = revcomp(minus), len(minus)
    (a0, b0), (a1, b1) = mspans
    genes += [cds("E", "c2", n - b0, n - a0, "-", "MNNPPFFWR", "GTG"),
              cds("F", "c2", n - b1, n - a1, "-", "MHHDDYYQR", "TTG")]  # fmt: skip
    nt = "GTG" + back_translate("FAGHWLDK") + "TAA"
    c3 = nt[6:] + SPACER * 4 + nt[:6]
    a = len(c3) - 6
    genes.append(cds("H", "c3", a, a + len(nt), "+", "MFAGHWLDK", "GTG"))
    contigs = {"c1": c1, "c2": c2, "c3": c3}
    reps = replicons(contigs, ["c1", "c3"])
    yes = {START_PEPTIDE, "VTEKLLAGGPR", "MAGDEFHHNK", "LAGDEFHHNK", "MSAGLLKVEPTLDEK",
           "VWQQDEHHK", "MNNPPFFWR", "VNNPPFFWR", "MHHDDYYQR", "LHHDDYYQR", "MFAGHWLDK",
           "VFAGHWLDK", "MTEKLLA", "LYYFFDEK", "MAALYYFFDEK"}  # fmt: skip
    no = {"MWQQDEHHK", "MYYFFDEK", "ATEKLLAGGPR", "MHHHFFFWWK", "SMTEKLLAGGPR"}
    return contigs, reps, model(reps, genes), yes, no


def _step_exact(
    tmp: Path, contigs: Mapping[str, str], gm: Any, peptides: Sequence[str]
) -> set[str]:
    """Target peptides of GIPPA's registered ``sixframe_index`` and ``map_exact`` steps."""
    from gippa.config import resolve_config
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    tmp.mkdir(parents=True, exist_ok=True)
    ann = write_annotation_parquet(gm, tmp / "ann")
    peps = tmp / "peptides.parquet"
    pd.DataFrame(
        {"peptide": list(peptides), "source": "database", "q_value": 0.001, "pep": 0.01,
         "score": 0.9, "is_decoy": False, "search_db": "other"}
    ).to_parquet(peps, index=False)  # fmt: skip
    ins = {"genome": write_genome_plain(contigs, tmp / "genome.fa"), "peptides": peps,
           "replicons": ann["replicons"], "genes": ann["genes"]}  # fmt: skip
    cfg = resolve_config(overrides=["seed=7"])
    six, hx = tmp / "sixframe.npz", tmp / "hits_exact.parquet"
    reg.run_step("sixframe_index", cfg, {k: ins[k] for k in reg.get("sixframe_index").inputs},
                 {"sixframe": six})  # fmt: skip
    reg.run_step("map_exact", cfg, {"peptides": peps, "sixframe": six}, {"hits_exact": hx})
    hits = pd.read_parquet(hx)
    return {str(p) for p in hits.loc[hits["search"] == "target", "peptide"]}


def _route_exact(
    contigs: Mapping[str, str], reps: Any, gm: Any, peptides: Sequence[str]
) -> set[str]:
    """GIPPA's start-cell route as ``sixframe_index`` builds it, with the reverse decoy (cheap)."""
    from gippa.config import resolve_config
    from gippa.mapping.decoys import build_decoys
    from gippa.mapping.exact import match_exact
    from gippa.mapping.initiator import (
        StartCells,
        annotated_starts,
        mirror_start_cells,
        start_frame_cells,
    )
    from gippa.mapping.sixframe import build_sixframe
    from gippa.rng import for_step

    cfg = resolve_config()
    table = cfg.annotation.translation_table
    index = build_sixframe(dict(contigs), reps, table=table)
    decoys = {"reverse": build_decoys(index, "reverse", 50, for_step("sixframe_index", 7))}
    sites, _ = annotated_starts(gm, contigs, cfg.annotation.start_codons, table=table)
    cells = start_frame_cells(index, sites)
    dc = {"reverse": mirror_start_cells("reverse", decoys["reverse"], cells)}
    hits = match_exact(list(peptides), index, decoys, dict.fromkeys(peptides, "reverse"),
                       starts=StartCells(sites, cells, dc))  # fmt: skip
    return {str(p) for p in hits.loc[hits["search"] == "target", "peptide"]}


def test_parity_with_gippas_mapping_steps(tmp_path: Path) -> None:
    contigs, reps, gm, yes, no = _parity_world()
    peptides = sorted(yes | no)
    gippa = _step_exact(tmp_path, contigs, gm, peptides)
    assert gippa == yes  # the fixture says what it means to, in GIPPA's own words
    idx = E().build_index(contigs, reps, [], annotation=gm)
    got = E().present(peptides, idx)
    assert {p for p in peptides if got[p].rule == "exact"} == gippa


@settings(max_examples=40, derandomize=True, deadline=None)
@given(
    codon=st.sampled_from(sorted(PLAIN)),
    strand=st.sampled_from("+-"),
    annotated=st.booleans(),
    rest=st.text(alphabet="ACDEFGHKLNPQRSTVWY", min_size=6, max_size=18),
)
def test_property_present_agrees_with_gippas_route(
    codon: str, strand: str, annotated: bool, rest: str
) -> None:
    """For any N-terminal text: ``present``'s exact set is GIPPA's start-cell route's; the annotated
    M-initial peptide is exact; and non-M-initial peptides do not depend on the annotation."""
    text, a, b = start_contig(codon, rest)
    gene = cds("G", "nissle", a, b, "+", "M" + rest, codon)
    contigs, genes = on_strand(text, [gene] if annotated else [], strand)
    reps = replicons(contigs, ["nissle"])
    gm = model(reps, genes)
    plain = PLAIN[codon] + rest
    peptides = ["M" + rest, plain, "M" + rest[:6], "A" + rest]
    idx = E().build_index(contigs, reps, [], annotation=gm)
    got = E().present(peptides, idx)
    assert {p for p in peptides if got[p].rule == "exact"} == _route_exact(
        contigs, reps, gm, peptides
    )
    if annotated:
        assert got["M" + rest].rule == "exact"
    without = E().present([plain], E().build_index(contigs, reps, [], annotation=model(reps, [])))
    assert got[plain] == without[plain]


# items 1, 2, 4, 5 through the CLI: cross-check 1 and a missing annotation

STRAIN_DB = ("MSAGLLKVEPTLDEKAAGGYWR",)
START_PROTEIN = START_PEPTIDE + "AAGK"
EXIT1 = "1 discovered peptide(s) have a route = exact row"


def _cli_args(
    root: Path,
    contigs: Mapping[str, str],
    genes: Sequence[Any],
    *extra: str,
    rows: Sequence[Mapping[str, Any]] = START_ROWS,
    k12: Sequence[str] = (*K12, START_PEPTIDE),
    strain_db: Sequence[str] = (*STRAIN_DB, START_PROTEIN),
    annotation: Path | None = None,
) -> list[str]:
    ann = annotation or write_annotation_json(
        root / "nissle.bakta.json", contigs, ["nissle"], genes
    )
    flags = {
        "--pep-a1": write_a1(root / "a1.tsv.gz", list(rows)),
        "--k12-fasta": write_faa(root / "UP000000625.fasta", k12),
        "--contaminants": write_faa(root / "contaminants.fasta", CONTAMINANTS, "Cont_P"),
        "--genome": write_genome(root / "nissle.fna", contigs, ["nissle"]),
        "--strain-db": write_faa(root / "strain_db.faa", strain_db, "NISSLE_"),
        "--runs-manifest": write_runs_manifest(root / "files.tsv"),
        "--annotation": ann,
    }
    return ["compute", *(x for f, p in flags.items() for x in (f, str(p))), *extra]


def _show(res: Any) -> str:
    return f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


def _bucket_rule(out: Path, peptide: str) -> list[str]:
    lines = (out / "entrapment_discoveries.tsv").read_text(encoding="utf-8").splitlines()
    return next(ln.split("\t")[1:3] for ln in lines if ln.split("\t")[0] == peptide)


@pytest.mark.parametrize(("codon", "strand"), CODON_STRAND)
def test_guard_passes_for_an_annotated_start_with_a_route_exact_row(
    tmp_path: Path, codon: str, strand: str
) -> None:
    """The A1 freeze run's failure, reduced: an ``initiator_m`` row, present here only through
    ``strain_db`` before the fix."""
    contigs, genes = start_world(codon, strand)
    a1 = write_out_a1(tmp_path / "A1", ["MSAGLLK", "VEPTIDEK", START_PEPTIDE], 2, 5,
                      {START_PEPTIDE: "initiator_m"})  # fmt: skip
    out = tmp_path / "metrics"
    res = run_entrapment(*_cli_args(tmp_path / "in", contigs, genes), "--out", str(out),
                         "--out-a1", str(a1))  # fmt: skip
    assert res.returncode == 0, _show(res)
    assert _bucket_rule(out, START_PEPTIDE) == ["target", "exact"]
    run = read_json(out / "entrapment.json")["cross_checks"]["run"]
    assert run["exact_route_agrees"] is True
    assert run["n_exact_route"] == 3
    assert run["exact_rejected"] == 0


@pytest.mark.parametrize("case", ["no_gene", "other_start_only", "internal_codon"])
@pytest.mark.parametrize("strain_db", [STRAIN_DB, (*STRAIN_DB, START_PROTEIN)], ids=["no", "db"])
@pytest.mark.parametrize(("codon", "strand"), [("GTG", "+"), ("TTG", "-")])
def test_guard_still_fires_on_an_unannotated_codon(
    tmp_path: Path, case: str, strain_db: tuple[str, ...], codon: str, strand: str
) -> None:
    """Not weakened: an A1 ``route = exact`` row for an M at a codon this annotation does not make
    a start is a real disagreement, even when the row carries ``initiator_m`` (no exemption)."""
    contigs, genes = unannotated_world(case, codon, strand)
    a1 = write_out_a1(tmp_path / "A1", ["MSAGLLK", "VEPTIDEK", START_PEPTIDE], 2, 5,
                      {START_PEPTIDE: "initiator_m"})  # fmt: skip
    out = tmp_path / "metrics"
    args = _cli_args(tmp_path / "in", contigs, genes, strain_db=strain_db)
    res = run_entrapment(*args, "--out", str(out), "--out-a1", str(a1))
    assert res.returncode == 1, _show(res)
    assert "inconsistent" in res.stderr and EXIT1 in res.stderr, _show(res)
    assert not out.exists() or not any(out.iterdir())


def test_guard_counts_only_the_unannotated_peptide(tmp_path: Path) -> None:
    """Both M-initial, both ``initiator_m`` in A1, one start annotated: exactly one disagreement."""
    contigs, genes = unannotated_world("other_start_only", "GTG", "+")
    peps = ["MSAGLLK", "VEPTIDEK", START_PEPTIDE, OTHER_PEPTIDE]
    a1 = write_out_a1(tmp_path / "A1", peps, 2, 6, dict.fromkeys(peps[2:], "initiator_m"))
    rows = [*START_ROWS, row("WT_1", OTHER_PEPTIDE, 2, 0.002)]
    out = tmp_path / "metrics"
    args = _cli_args(tmp_path / "in", contigs, genes, rows=rows,
                     k12=(*K12, START_PEPTIDE, OTHER_PEPTIDE))  # fmt: skip
    res = run_entrapment(*args, "--out", str(out), "--out-a1", str(a1))
    assert res.returncode == 1, _show(res)
    assert EXIT1 in res.stderr, _show(res)


@pytest.mark.parametrize("codon", sorted(PLAIN))
def test_guard_passes_for_the_plain_translation_at_a_start(tmp_path: Path, codon: str) -> None:
    contigs, genes = start_world(codon)
    plain = PLAIN[codon] + START_REST
    a1 = write_out_a1(tmp_path / "A1", ["MSAGLLK", "VEPTIDEK", plain], 2, 5)
    out = tmp_path / "metrics"
    args = _cli_args(tmp_path / "in", contigs, genes, k12=(*K12, plain),
                     rows=[*EXAMPLE2_ROWS, row("WT_1", plain, 2, 0.002)])  # fmt: skip
    res = run_entrapment(*args, "--out", str(out), "--out-a1", str(a1))
    assert res.returncode == 0, _show(res)
    assert _bucket_rule(out, plain) == ["target", "exact"]


def test_unused_start_cells_change_nothing_but_the_annotation_input(tmp_path: Path) -> None:
    """Same genome, with and without the site, and no peptide using it: byte-identical TSV, and
    ``entrapment.json`` differs only in the record of the annotation file."""
    contigs, genes = start_world("GTG")
    docs, tsvs = [], []
    for name, gs in (("with", genes), ("without", [])):
        out = tmp_path / f"m_{name}"
        args = _cli_args(tmp_path / name, contigs, gs, rows=EXAMPLE2_ROWS, k12=K12,
                         strain_db=STRAIN_DB)  # fmt: skip
        res = run_entrapment(*args, "--out", str(out))
        assert res.returncode == 0, _show(res)
        docs.append(read_json(out / "entrapment.json"))
        tsvs.append((out / "entrapment_discoveries.tsv").read_bytes())
    assert tsvs[0] == tsvs[1]
    a, b = docs
    assert a["inputs"]["annotation"]["sha256"] != b["inputs"]["annotation"]["sha256"]
    assert {k: v for k, v in a.items() if k != "inputs"} == {
        k: v for k, v in b.items() if k != "inputs"
    }


def test_annotation_is_recorded_with_its_hash_and_a_gzipped_json_is_read(tmp_path: Path) -> None:
    """The recorded pilot annotation is ``.json.gz``; the credit works through it."""
    from helpers.entrapment_m10_07 import sha256

    contigs, genes = start_world("GTG")
    ann = write_annotation_json(tmp_path / "in" / "a.json.gz", contigs, ["nissle"], genes)
    out = tmp_path / "metrics"
    res = run_entrapment(*_cli_args(tmp_path / "in", contigs, genes, annotation=ann), "--out",
                         str(out))  # fmt: skip
    assert res.returncode == 0, _show(res)
    assert _bucket_rule(out, START_PEPTIDE) == ["target", "exact"]
    rec = read_json(out / "entrapment.json")["inputs"]["annotation"]
    assert rec["sha256"] == sha256(ann)
    assert Path(rec["path"]) == ann.resolve()


# item 5: a missing or unusable annotation fails clearly


def test_build_index_requires_the_annotation() -> None:
    contigs, _ = start_world("GTG")
    with pytest.raises(TypeError, match="annotation"):
        E().build_index(contigs, replicons(contigs, ["nissle"]), [])


def _exit2(res: Any, out: Path, *needles: str) -> None:
    """Exit 2 for this reason: argparse also exits 2 on a flag it does not know."""
    assert res.returncode == 2, _show(res)
    assert "unrecognized arguments" not in res.stderr, _show(res)
    assert all(n in res.stderr for n in needles), _show(res)
    assert not out.exists() or not any(out.iterdir())


def test_cli_without_annotation_exits_2_naming_the_flag(tmp_path: Path) -> None:
    contigs, genes = start_world("GTG")
    args = _cli_args(tmp_path / "in", contigs, genes)
    i = args.index("--annotation")
    out = tmp_path / "metrics"
    res = run_entrapment(*args[:i], *args[i + 2 :], "--out", str(out))
    _exit2(res, out, "required", "--annotation")


def test_cli_with_a_missing_annotation_file_exits_2_naming_it(tmp_path: Path) -> None:
    contigs, genes = start_world("GTG")
    missing = tmp_path / "no-such.bakta.json"
    args = _cli_args(tmp_path / "in", contigs, genes, annotation=missing)
    out = tmp_path / "metrics"
    _exit2(run_entrapment(*args, "--out", str(out)), out, "precondition", str(missing))


def _bad_annotation(kind: str, root: Path, contigs: dict[str, str], genes: list[Any]) -> Path:
    path = root / "bad.bakta.json"
    root.mkdir(parents=True, exist_ok=True)
    if kind == "not_json":
        path.write_text("{not json", encoding="utf-8")
    elif kind == "not_an_object":
        path.write_text("[]", encoding="utf-8")
    elif kind == "another_table":  # GIPPA's parse_bakta refuses this too (AnnotationTableMismatch)
        write_annotation_json(path, contigs, ["nissle"], genes, table=4)
    elif kind == "another_contig":
        other = {"plasmid9": contigs["nissle"]}
        gs = [dataclasses.replace(g, interval=dataclasses.replace(g.interval, contig="plasmid9"))
              for g in genes]  # fmt: skip
        write_annotation_json(path, other, ["plasmid9"], gs)
    else:
        assert kind == "another_length"
        write_annotation_json(path, {"nissle": contigs["nissle"] + "ACGTAC"}, ["nissle"], genes)
    return path


@pytest.mark.parametrize(
    "kind", ["not_json", "not_an_object", "another_table", "another_contig", "another_length"]
)
def test_cli_with_an_unusable_annotation_exits_2_naming_it(tmp_path: Path, kind: str) -> None:
    contigs, genes = start_world("GTG")
    bad = _bad_annotation(kind, tmp_path / "ann", contigs, genes)
    args = _cli_args(tmp_path / "in", contigs, genes, annotation=bad)
    out = tmp_path / "metrics"
    _exit2(run_entrapment(*args, "--out", str(out)), out, "precondition", str(bad))
