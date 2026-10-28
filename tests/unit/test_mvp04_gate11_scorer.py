"""MVP-04 gate 11 scorer fixes (``scripts/dev/mvp04_summary.py`` ``gate11``).

The statistician's K-12 audit (docs/methods/variants-mvp04-k12-audit.md on stat/m7-k12-flags)
found four scorer/truth bugs; the lead fixed all four as plain bug fixes, thresholds unchanged
(precision >= 0.95, recall >= 0.90, n_high >= 30). Contract pinned here, on small hand-made
tables (no gippa, no gippasim):

1. **I/L folding.** A ``high`` row is true when (``contig``, ``codon_start``, folded
   ``alt_aa``) is a truth site, where ``I``, ``L`` and ``J`` fold to one token on both sides:
   ALT ``J`` matches truth ``I`` or ``L``; any other residue is compared literally.
2. **Initiator Met.** A ``high`` row is an initiator row when its ``flags`` hold
   ``start_codon``, or when ``alt_aa = M``, ``ref_codon`` is GTG/TTG/CTG and (``contig``,
   ``strand``, ``codon_start``) is the start codon of a ``truth/novel_orfs.tsv`` ORF (codon
   lowest base ``true_start`` on ``+``, ``true_start - 2`` on ``-``). Initiator rows leave the
   precision denominator and are reported as ``n_initiator_m_excluded`` (high rows only).
   ``value`` stays ``n_high`` (every high row).
3. **Novel-ORF SNVs.** The truth sites are ``truth/alt_sites.tsv`` plus the gene-less
   ``truth/novel_alt_sites.tsv`` (columns ``contig, codon_start, strand, orf_id, ref_aa,
   alt_aa, n_true_peptides``); an absent ``novel_alt_sites.tsv`` is an empty one.
4. **Recall denominator.** A truth site is a recall site when at least 2 distinct I/L-folded
   true peptides of ``truth/peptides.tsv`` (``origin != decoy``, ``is_wrong = false``,
   ``locus_tag`` = the site's ``gene`` / ``orf_id``, same ``contig``, ``start <= codon_start``
   and ``codon_start + 2 <= end``, 1-based closed) have length >= ``mapping.tolerant.min_len``
   of ``OUTDIR/resolved_config.yaml``. ``n_recall_sites`` counts them.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pytest
from helpers.mvp04_gate11_data import (
    CONTIG,
    FLAGS,
    GATE_TEXT,
    N_SITES,
    Case,
    base_case,
    codon,
    gene,
    mod,
)


@pytest.fixture
def case(tmp_path: Path) -> Case:
    return base_case(tmp_path)


# ---------------------------------------------------------------- the base case


def test_base_case_passes_with_no_initiator_rows(case: Case) -> None:
    res = case.gate()
    assert res["value"] == N_SITES
    assert res["precision"] == 1.0 and res["recall"] == 1.0
    assert res["n_recall_sites"] == N_SITES
    assert res["n_initiator_m_excluded"] == 0
    assert res["gate_ok"] is True


def test_gate_text_and_thresholds_are_unchanged(case: Case) -> None:
    """R3: the fix changes the scoring, not the gate (0.95 / 0.90 / 30)."""
    assert case.gate()["gate"] == GATE_TEXT


def test_precision_just_below_threshold_still_fails(case: Case) -> None:
    """4 false of 41 non-initiator high rows: 37/41 = 0.902 < 0.95, a fail as before."""
    for k in range(3):
        case.variants[k][8] = "W"
    case.var(codon(0) + 3, "W")
    res = case.gate()
    assert res["precision"] == pytest.approx(37 / 41)
    assert res["gate_ok"] is False and res["status"] == "fail"


# ---------------------------------------------------------------- bug 1: I/L folding


def test_alt_j_matches_truth_i_and_truth_l(case: Case) -> None:
    """Audit item 1: GIPPA writes the folded I/L ALT token ``J`` (M7-02 rule 5)."""
    for k in range(5):
        case.alt[k][5] = "I"
        case.variants[k][8] = "J"
    for k in range(5, 10):
        case.alt[k][5] = "L"
        case.variants[k][8] = "J"
    res = case.gate()
    assert res["precision"] == 1.0
    assert res["recall"] == 1.0


def test_alt_i_and_l_fold_against_each_other(case: Case) -> None:
    case.alt[0][5], case.variants[0][8] = "I", "L"
    case.alt[1][5], case.variants[1][8] = "L", "I"
    case.alt[2][5], case.variants[2][8] = "J", "I"
    res = case.gate()
    assert res["precision"] == 1.0 and res["recall"] == 1.0


def test_folding_does_not_touch_other_residues(case: Case) -> None:
    """Truth V, ALT J: still false (only I, L and J fold)."""
    case.variants[0][8] = "J"
    res = case.gate()
    assert res["precision"] == pytest.approx(39 / 40)
    assert res["recall"] == pytest.approx(39 / 40)


# ---------------------------------------------------------------- bug 2: initiator Met


def test_start_codon_flagged_row_leaves_the_precision_denominator(case: Case) -> None:
    """Audit item 2: ALT M at an annotated GTG start (``aa_pos`` 1, flag ``start_codon``) is no
    truth site (gippasim leaves initiators out of truth) and no longer counts as false."""
    case.var(50_000, "M", gene=gene(99), ref_codon="GTG", ref="V",
             flags=FLAGS + ";start_codon")  # fmt: skip
    res = case.gate()
    assert res["n_initiator_m_excluded"] == 1
    assert res["precision"] == 1.0
    assert res["value"] == N_SITES + 1  # value stays n_high


@pytest.mark.parametrize("ref_codon,ref", [("GTG", "V"), ("TTG", "L"), ("CTG", "L")])
def test_alt_m_at_a_plus_strand_novel_start_is_an_initiator(
    case: Case, ref_codon: str, ref: str
) -> None:
    """Audit item 2: ALT M over a non-ATG start of a planted novel ORF (gene-less, no flag)."""
    case.var(100_001, "M", ref_codon=ref_codon, ref=ref)
    res = case.gate()
    assert res["n_initiator_m_excluded"] == 1
    assert res["precision"] == 1.0


def test_alt_m_at_a_minus_strand_novel_start_is_an_initiator(case: Case) -> None:
    """On ``-`` the start codon's lowest base is ``true_start - 2`` (= ``end - 2``)."""
    case.var(200_298, "M", strand="-", ref_codon="TTG", ref="L")
    res = case.gate()
    assert res["n_initiator_m_excluded"] == 1
    assert res["precision"] == 1.0


def test_alt_m_off_a_novel_start_is_still_false(case: Case) -> None:
    """Controls: one codon downstream of the start, the other strand, and an ATG/CTC ref codon
    at the start are no initiator rows; each counts as a false high row."""
    case.var(100_004, "M", ref_codon="GTG", ref="V")  # next codon
    case.var(100_001, "M", strand="-", ref_codon="GTG", ref="V")  # other strand
    case.var(200_298, "M", strand="-", ref_codon="CTC", ref="L")  # not a start codon
    res = case.gate()
    assert res["n_initiator_m_excluded"] == 0
    assert res["precision"] == pytest.approx(40 / 43)


def test_alt_other_than_m_at_a_novel_start_is_not_an_initiator(case: Case) -> None:
    case.var(100_001, "V", ref_codon="TTG", ref="L")
    res = case.gate()
    assert res["n_initiator_m_excluded"] == 0
    assert res["precision"] == pytest.approx(40 / 41)


def test_initiator_count_covers_high_rows_only(case: Case) -> None:
    case.var(100_001, "M", conf="candidate", ref_codon="GTG", ref="V")
    case.var(50_000, "M", conf="candidate", gene=gene(99), ref_codon="GTG", ref="V",
             flags=FLAGS + ";start_codon")  # fmt: skip
    res = case.gate()
    assert res["n_initiator_m_excluded"] == 0
    assert res["value"] == N_SITES and res["precision"] == 1.0


def test_the_k12_audit_mix_scores_one(case: Case) -> None:
    """Audit shape in miniature: folded J rows, flagged and unflagged initiators and novel-ORF
    SNVs together give precision 1 (K-12: 55/55 once all three are fixed)."""
    case.alt[0][5], case.variants[0][8] = "I", "J"
    case.var(50_000, "M", gene=gene(99), ref_codon="GTG", ref="V",
             flags=FLAGS + ";start_codon")  # fmt: skip
    case.var(100_001, "M", ref_codon="TTG", ref="L")
    assert case.novel_alt is not None
    case.novel_alt.append([CONTIG, 100_031, "+", "NOVEL_0001", "P", "S", 2])
    case.true_peps("NOVEL_0001", 100_031)
    case.var(100_031, "S", ref_codon="CCC", ref="P")
    res = case.gate()
    assert res["precision"] == 1.0 and res["recall"] == 1.0
    assert res["n_initiator_m_excluded"] == 2
    assert res["value"] == N_SITES + 3


# ---------------------------------------------------------------- bug 3: novel-ORF SNVs


def test_a_gene_less_row_at_a_novel_orf_snv_is_true(case: Case) -> None:
    """Audit item 3: a real sample SNV inside a planted novel ORF (no gene) is in truth via
    ``truth/novel_alt_sites.tsv``."""
    assert case.novel_alt is not None
    case.novel_alt.append([CONTIG, 200_100, "-", "NOVEL_0002", "P", "S", 2])
    case.true_peps("NOVEL_0002", 200_100, strand="-")
    case.var(200_100, "S", strand="-", ref_codon="CCC", ref="P")
    res = case.gate()
    assert res["precision"] == 1.0
    assert res["recall"] == 1.0 and res["n_recall_sites"] == N_SITES + 1


def test_the_same_row_without_the_novel_truth_row_is_false(case: Case) -> None:
    case.var(200_100, "S", strand="-", ref_codon="CCC", ref="P")
    res = case.gate()
    assert res["precision"] == pytest.approx(40 / 41)


def test_novel_sites_enter_recall(case: Case) -> None:
    """An uncalled novel-ORF site with 2 long true peptides lowers recall."""
    assert case.novel_alt is not None
    case.novel_alt.append([CONTIG, 100_031, "+", "NOVEL_0001", "P", "S", 2])
    case.true_peps("NOVEL_0001", 100_031)
    res = case.gate()
    assert res["n_recall_sites"] == N_SITES + 1
    assert res["recall"] == pytest.approx(40 / 41)


def test_inputs_without_novel_alt_sites_still_score(case: Case) -> None:
    """Inputs written before the fix have no ``novel_alt_sites.tsv``: an empty table."""
    case.novel_alt = None
    res = case.gate()
    assert not (case.inp / "truth" / "novel_alt_sites.tsv").exists()
    assert res["precision"] == 1.0 and res["recall"] == 1.0


# ---------------------------------------------------------------- bug 4: recall denominator


def _short_site(case: Case, k: int, length: int) -> None:
    """Truth site ``k`` (>= N_SITES), uncalled, covered by 2 true peptides of ``length``."""
    case.alt.append([CONTIG, codon(k), "+", gene(k), "A", "V", 2, "false"])
    case.true_peps(gene(k), codon(k), length=length)


def test_sites_seen_only_by_short_peptides_leave_recall(case: Case) -> None:
    """Audit item 4: peptides shorter than ``mapping.tolerant.min_len`` (10) can never support
    an ALT, so a site covered only by 9-mers is no recall site."""
    _short_site(case, N_SITES, 9)
    res = case.gate()
    assert res["n_recall_sites"] == N_SITES
    assert res["recall"] == 1.0


def test_a_site_with_one_long_and_one_short_peptide_is_no_recall_site(case: Case) -> None:
    k = N_SITES
    case.alt.append([CONTIG, codon(k), "+", gene(k), "A", "V", 2, "false"])
    case.true_peps(gene(k), codon(k), n=1, length=10)
    case.true_peps(gene(k), codon(k), n=1, length=9)
    assert case.gate()["n_recall_sites"] == N_SITES


def test_peptides_at_min_len_count(case: Case) -> None:
    _short_site(case, N_SITES, 10)
    res = case.gate()
    assert res["n_recall_sites"] == N_SITES + 1
    assert res["recall"] == pytest.approx(40 / 41)


def test_min_len_is_read_from_the_resolved_config(case: Case) -> None:
    """Same tables, ``mapping.tolerant.min_len`` 8 then 13: not a hard-coded 10."""
    _short_site(case, N_SITES, 9)
    case.min_len = 8
    assert case.gate()["n_recall_sites"] == N_SITES + 1
    case.min_len = 13
    assert case.gate()["n_recall_sites"] == 0
    assert case.gate()["recall"] is None


def test_recall_ignores_the_truth_n_true_peptides_column(case: Case) -> None:
    """The count comes from ``truth/peptides.tsv``: a site whose column says 5 but has no
    covering peptide is no recall site; one whose column says 0 but has 2 long ones is."""
    case.alt[0][6] = 0
    case.alt.append([CONTIG, codon(N_SITES), "+", gene(N_SITES), "A", "V", 5, "false"])
    res = case.gate()
    assert res["n_recall_sites"] == N_SITES
    assert res["recall"] == 1.0


def test_wrong_decoy_duplicate_and_non_covering_peptides_do_not_count(case: Case) -> None:
    """Each extra site has 1 good long peptide plus one that must not count: a wrong one, a
    decoy, an I/L twin of the good one, one of another locus, and one ending 1 base short of
    the codon. None reaches 2."""
    kinds = ("wrong", "decoy", "twin", "locus", "short_cover")
    for j, kind in enumerate(kinds):
        k = N_SITES + j
        cs = codon(k)
        case.alt.append([CONTIG, cs, "+", gene(k), "A", "V", 2, "false"])
        case.true_peps(gene(k), cs, n=1)
        good = case.peptides[-1]
        if kind == "wrong":
            case.true_peps(gene(k), cs, n=1, is_wrong="true")
        elif kind == "decoy":
            case.true_peps(gene(k), cs, n=1, origin="decoy")
        elif kind == "twin":
            good[0] = "L" + good[0][1:]
            case.peptides.append(["I" + good[0][1:], *good[1:]])
        elif kind == "locus":
            case.true_peps(gene(k + 100), cs, n=1)
        else:
            case.true_peps(gene(k), cs, n=1, end=cs + 1)
    assert case.gate()["n_recall_sites"] == N_SITES


def test_two_folded_distinct_peptides_count_twice(case: Case) -> None:
    """Control for the twin case: two peptides that differ outside I/L are two."""
    _short_site(case, N_SITES, 11)
    assert case.gate()["n_recall_sites"] == N_SITES + 1


# ---------------------------------------------------------------- JSON output


def test_evidence_summary_json_reports_the_initiator_count(case: Case) -> None:
    """The gate's entry in ``evidence_summary.json`` carries ``n_initiator_m_excluded``."""
    case.var(100_001, "M", ref_codon="GTG", ref="V")
    case.var(50_000, "M", gene=gene(99), ref_codon="TTG", ref="L",
             flags=FLAGS + ";start_codon")  # fmt: skip
    case.write()
    rc = mod().main([str(case.out), str(case.inp)])
    assert rc == 0
    doc: Mapping[str, Any] = json.loads((case.out / "evidence_summary.json").read_text())
    g = doc["gates"]["11"]
    assert g["n_initiator_m_excluded"] == 2
    assert g["precision"] == 1.0 and g["recall"] == 1.0
    assert not any("sign_off" in k or "signoff" in k for k in g), "plain bug fix: no sign-off"


# ------------------------------------------------ min_len fails closed (#324 item 3)

_BAD_CONFIGS = {
    "missing_file": None,
    "invalid_yaml": "mapping: [unclosed\n  tolerant: {min_len: 10\n",
    "null_document": "",
    "mapping_null": "mapping: null\n",
    "missing_tolerant": "mapping:\n  exact:\n    min_len: 7\n",
    "missing_min_len": "mapping:\n  tolerant:\n    enabled: true\n",
    "non_numeric": "mapping:\n  tolerant:\n    enabled: true\n    min_len: ten\n",
    "list_value": "mapping:\n  tolerant:\n    min_len: [10]\n",
    "bool_true": "mapping:\n  tolerant:\n    enabled: true\n    min_len: true\n",
    "bool_false": "mapping:\n  tolerant:\n    enabled: true\n    min_len: false\n",
}


@pytest.mark.parametrize("kind", sorted(_BAD_CONFIGS))
def test_unreadable_min_len_gives_no_recall_and_a_fail(case: Case, kind: str) -> None:
    """No usable ``mapping.tolerant.min_len``: no recall sites, ``recall`` None, gate fails.
    Never a fallback (0 or 10): the base case otherwise passes with recall 1."""
    case.write()
    cfg = case.out / "resolved_config.yaml"
    text = _BAD_CONFIGS[kind]
    if text is None:
        cfg.unlink()
    else:
        cfg.write_text(text)
    res: dict[str, Any] = mod().gate11(case.out, case.inp)
    assert res["recall"] is None
    assert res["recall_min_len"] is None
    assert res["n_recall_sites"] == 0
    assert res["gate_ok"] is False and res["status"] == "fail"
    assert res["precision"] == 1.0  # precision does not depend on the config


def test_unreadable_min_len_is_a_fail_in_the_summary_json(case: Case) -> None:
    case.write()
    (case.out / "resolved_config.yaml").unlink()
    assert mod().main([str(case.out), str(case.inp)]) == 0
    doc: Mapping[str, Any] = json.loads((case.out / "evidence_summary.json").read_text())
    g = doc["gates"]["11"]
    assert g["recall"] is None and g["recall_min_len"] is None and g["status"] == "fail"


@pytest.mark.parametrize("min_len,want", [(7, 40 / 41), (12, 1.0)])
def test_the_configured_min_len_moves_recall(case: Case, min_len: int, want: float) -> None:
    """Same tables (an uncalled site seen by two 8-mers): min_len 7 counts it, 12 does not."""
    _short_site(case, N_SITES, 8)
    case.min_len = min_len
    res = case.gate()
    assert res["recall_min_len"] == min_len
    assert res["recall"] == pytest.approx(want)
    assert res["n_recall_sites"] == N_SITES + (min_len <= 8)


# ---------------------------------------------------------------- coverage boundaries (#324 item 4)


@pytest.mark.parametrize(
    "where,counted",
    [("start_at_codon", True), ("end_at_codon_plus_2", True),
     ("start_at_codon_plus_1", False), ("end_at_codon_plus_1", False)],
)  # fmt: skip
def test_coverage_is_1_based_closed_and_holds_the_whole_codon(
    case: Case, where: str, counted: bool
) -> None:
    """Two true 12-mers per extra (uncalled) site, both placed at the boundary under test: the
    codon ``cs..cs+2`` must lie in ``start..end`` (1-based closed)."""
    k = N_SITES
    cs = codon(k)
    case.alt.append([CONTIG, cs, "+", gene(k), "A", "V", 2, "false"])
    span = 3 * 12 - 1
    start, end = {
        "start_at_codon": (cs, cs + span),
        "end_at_codon_plus_2": (cs + 2 - span, cs + 2),
        "start_at_codon_plus_1": (cs + 1, cs + 1 + span),
        "end_at_codon_plus_1": (cs + 1 - span, cs + 1),
    }[where]
    case.true_peps(gene(k), cs, start=start, end=end)
    res = case.gate()
    assert res["n_recall_sites"] == N_SITES + counted
    assert res["recall"] == pytest.approx(40 / 41 if counted else 1.0)
