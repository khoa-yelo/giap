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
* ``build_index(genome, replicons, strain_db=None)`` (``genome``: contig -> nucleotides,
  ``replicons``: contig -> ``giap.model.coords.Replicon``, ``strain_db``: protein sequences or
  None for "L2 answered no"), an opaque index;
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

from collections.abc import Iterable, Mapping, Sequence
from typing import Any

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
    SHARED_1,
    SPACER,
    back_translate,
    codes,
    coding_text,
    flagged,
    frame,
    load_entrapment,
    origin_split,
    reasons,
    replicons,
    revcomp,
    row,
)
from hypothesis import given, settings
from hypothesis import strategies as st

TOL = 1e-12


def E() -> Any:
    return load_entrapment()


def index_of(
    contigs: Mapping[str, str], circular: Iterable[str] = (), strain_db: Any = None
) -> Any:
    return E().build_index(dict(contigs), replicons(contigs, circular), strain_db)


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
    from giap.mapping.masses import MOD_DELTA, RESIDUE_MASS

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
