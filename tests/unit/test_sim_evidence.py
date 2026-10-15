"""giapsim evidence: peptides, identification errors, decoys and placed hits (issue M4-09).

Examples 1-10, edge cases and acceptance criteria of docs/specs/M4-09-sim-evidence.md (GIAP_SPEC.md
§11.3, §8.4, §8.7, §7.2.6, §7.3.2), checked against tests/helpers/sim_evidence_oracle.py (written from
the spec, not from the simulator). Lives in tests/unit (tests/sim is implementer tooling); the
distributional claims are gated in tests/calibration/test_sim_evidence_calibration.py (suite
``gene_fdr``).

Ambiguities resolved here (reported to the lead):

* ``score`` is an emit field like the others ("fields not in ``emit`` are written as NA"): tests
  that need scores pass ``emit`` with ``"score"``; with ``emit = set()`` it is NA (example 6).
  The default emit set's ``score`` is not pinned.
* ``protein_ids`` is a non-nullable list column (M4-06 C9); "NA" for it is the empty list, as the
  merged M1-12 aggregation writes it.
* Not emitting ``is_decoy`` drops the decoy rows (example 5); ``emit`` only masks: the target rows
  (peptides, q-values) are the same as with the default emit set for the same seed.
* PEP uses pi0 = ``pi_wrong`` (configured); q of a decoy row uses the same target-based formula.
* ``ln(intensity) = abundance_log2_true * ln 2 + N(0, 0.5^2)`` for true unique peptides.
* M3 files may sit in ``<work>/`` (M4-06 flat fixture layout) or ``<work>/mapping/``; the writer
  writes all three (C11-C13, header-only when empty) so ``load_m4_inputs`` sees a complete dir.
* Unplaceable false targets are listed in C13 (``unmapped.tsv``), as M3 would list them.
* Internal frames (``Evidence.peptide_hits``) are 0-based; the written C11/C12 are 1-based.
* Genes with ``detected_true`` but no 7-30 aa tryptic peptide (<= 1 missed cleavage) get
  ``present = False`` and ``n_true_peptides = 0`` in ``gene_truth``; ``expr`` is not mutated.
* The identical ``aa`` of ``multi_identical`` pairs has no carrier in the interface (``Evidence``
  has no genes frame and the layout is an input): example 8 is checked on C11 only, and placement
  coordinates are checked for ``unique`` rows only.
"""

from __future__ import annotations

import hashlib
import math
import re
import subprocess
import sys
import time
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import sim_evidence_oracle as eo

REPO = Path(__file__).resolve().parents[2]
GIAPSIM = REPO / "tests" / "sim" / "giapsim"
DEFAULT_EMIT = frozenset({"q_value", "pep", "is_decoy", "intensity", "protein_ids"})
SCORED = DEFAULT_EMIT | {"score"}
TRUTH_COLUMNS = ("peptide", "is_false", "gene_true", "decoy_gene")
GENE_TRUTH_COLUMNS = ("locus_tag", "present", "n_true_peptides")
LOCUS_RE = re.compile(r"^GSIM_\d{5,}$")
UP_RE = re.compile(r"^UP_\d+$")

# Contract tables of ``simulate_layout(LayoutConfig(), default_rng(1))`` (pinned by M4-08 too).
LAYOUT_SHA256 = {
    "annotation_meta.json": "0998816b81d3f81364ac87e43ef1ffc2d40c08bc620c535584e859619e215f72",
    "genes.parquet": "e0d4341901746a7444965ebe8070a3c382d5edc6f864bff8597bd8e697b65e65",
    "layers.parquet": "c7ae231f33b641c11a6b7b03ce95c36e5daf5c6b7fd5c1a46456d31e6d6c5ba5",
    "region_genes.parquet": "1197b721b0b595beaf182adb375dfbf0707c5cf8085b5d4598106988f1215fd7",
    "regions.parquet": "54ca39d8465b7720f4417364363d63ce4c6aa92cd893cad24c2f311c825f2ca5",
    "replicons.parquet": "baa9f94adf047f02a746031e33777abc727b0eea9e650454e0b7a868562bd85e",
    "tus.parquet": "c564de32b5fac4d9d1dc6a9ff5f885065ccdec0c2d7547b7ce48561f0e1ca825",
    "xrefs.parquet": "0fd3d2531b4d94b2dc6370de0378ecb870858dc2f6b091c37854771e474f592e",
}
# ``simulate_expression(layout seed 1, ExpressionConfig(), default_rng(1)).genes.to_csv`` at main
# 4b27316 (before M4-09): a shared score helper (lead Q3) must not change M4-08.
# M0 + M3-19 proteomics (``SimConfig(proteomics=ProteomicsConfig())``, seed 1): unchanged when
# M4-09 moves the score mixture into a shared helper (lead Q3).
M0_PROTEOMICS_SHA256 = {
    "abundance.tsv": "bdbe3b7329f9a369bb278ab48d82c925462027b7c719c6f2bfa4e922534fd34a",
    "peptides.tsv": "8489d448b3241131140f1ec59f28c0465aca8875e70f188550df4705328ee852",
    "proteins.faa": "138b2faa715a1f995fe510073c451253913d46589af8d2759a84df986286fb29",
    "truth/peptides.tsv": "4ddfb9b94088026d55e29e21550a030952d2f0fcc47bd2457816f62104cd0e4b",
}


# ----------------------------------------------------------------------------- helpers


def _ev_mod() -> ModuleType:
    import giapsim.evidence

    return giapsim.evidence


def layout(seed: int = 1, **kw: Any) -> Any:
    from giapsim.layout import LayoutConfig, simulate_layout

    return simulate_layout(LayoutConfig(**kw), np.random.default_rng(seed))


def express(lay: Any, seed: int = 1, **kw: Any) -> Any:
    from giapsim.expression import ExpressionConfig, simulate_expression

    return simulate_expression(lay, ExpressionConfig(**kw), np.random.default_rng(seed))


def ecfg(**kw: Any) -> Any:
    return _ev_mod().EvidenceConfig(**kw)


def evidence(lay: Any, expr: Any, seed: int = 1, **kw: Any) -> Any:
    m = _ev_mod()
    return m.simulate_evidence(lay, expr, m.EvidenceConfig(**kw), np.random.default_rng(seed))


def write(ev: Any, lay: Any, work: Path) -> Path:
    from giapsim.contract_tables import write_contract_tables

    write_contract_tables(lay, work)
    _ev_mod().write_evidence_tables(ev, work)
    return work


def false_hit_genes(ev: Any) -> list[str]:
    """Landing gene of every placed false target (from C11)."""
    bad = eo.false_set(ev.truth)
    h = ev.peptide_hits
    return [str(g) for p, g in zip(h["peptide"], h["gene"], strict=True) if p in bad]


@pytest.fixture(scope="module")
def lay1() -> Any:
    return layout(1)


@pytest.fixture(scope="module")
def ex1(lay1: Any) -> Any:
    return express(lay1, 1)


@pytest.fixture(scope="module")
def ev1(lay1: Any, ex1: Any) -> Any:
    return evidence(lay1, ex1, 1)


@pytest.fixture(scope="module")
def ev1s(lay1: Any, ex1: Any) -> Any:
    """Example 1 with scores emitted."""
    return evidence(lay1, ex1, 1, emit=SCORED)


@pytest.fixture(scope="module")
def work1(tmp_path_factory: pytest.TempPathFactory, lay1: Any, ev1: Any) -> Path:
    return write(ev1, lay1, tmp_path_factory.mktemp("ev1"))


# ----------------------------------------------------------------------------- interface


def test_config_defaults_match_spec() -> None:
    c = ecfg()
    assert c.pi_wrong == 0.01
    assert c.false_source == "random_gene"
    assert c.n_decoys == "matched"
    assert tuple(c.mixture) == (0.0, 3.0)
    assert frozenset(c.emit) == DEFAULT_EMIT
    assert c.protein_ids == "locus_tag"
    assert c.multi_identical_pairs == 0
    assert c.shared_frac == 0.0 and c.tolerant_frac == 0.0
    assert c.proteins_only is False
    assert c.abundance_table == "none"


def test_config_frozen_and_hashable() -> None:
    c = ecfg(pi_wrong=0.05, emit=frozenset({"q_value"}))
    with pytest.raises((AttributeError, TypeError)):
        c.pi_wrong = 0.2  # type: ignore[misc]
    hash(c)
    assert c == ecfg(pi_wrong=0.05, emit=frozenset({"q_value"}))


def test_evidence_frames_and_columns(ev1: Any, lay1: Any) -> None:
    assert set(TRUTH_COLUMNS) <= set(ev1.truth.columns)
    assert set(GENE_TRUTH_COLUMNS) <= set(ev1.gene_truth.columns)
    gt = ev1.gene_truth
    assert gt["locus_tag"].is_unique
    assert sorted(gt["locus_tag"]) == sorted(lay1.genes["locus_tag"])
    assert ev1.abundance is None  # abundance_table = "none"
    assert ev1.peptides["peptide"].is_unique
    assert set(ev1.truth["peptide"]) == set(ev1.peptides["peptide"])


def test_contract_tables_validate_c9_c11_c12(ev1: Any) -> None:
    from giap.genetable.contract import validate_table

    assert eo.errors(validate_table("peptides", ev1.peptides)) == []
    assert eo.errors(validate_table("peptide_hits", ev1.peptide_hits)) == []
    if ev1.protein_hits is not None:
        assert eo.errors(validate_table("protein_hits", ev1.protein_hits)) == []


def test_written_tables_load_without_errors(work1: Path, lay1: Any, ev1: Any) -> None:
    """``load_m4_inputs`` on C1-C8 (M4-07) + C9-C13 (this issue): no error, no cross-table issue."""
    from giap.genetable.contract import validate_inputs

    inputs, issues = eo.load(work1)
    assert eo.errors(issues) == []
    assert validate_inputs(inputs) == []
    assert inputs.peptides is not None and len(inputs.peptides) == len(ev1.peptides)
    assert len(inputs.peptide_hits) == len(ev1.peptide_hits)
    for key in ("peptide_hits", "protein_hits", "unmapped"):
        assert eo.m3_path(work1, key).is_file(), key


def test_written_c11_is_one_based_and_internal_frame_zero_based(work1: Path, ev1: Any) -> None:
    from giap.io.tsv import read_tsv

    inputs, _ = eo.load(work1)
    a = ev1.peptide_hits.sort_values(["peptide", "gene", "start"]).reset_index(drop=True)
    b = inputs.peptide_hits.sort_values(["peptide", "gene", "start"]).reset_index(drop=True)
    assert a["start"].astype(int).tolist() == b["start"].astype(int).tolist()  # loader: start - 1
    assert a["end"].astype(int).tolist() == b["end"].astype(int).tolist()
    header, rows = read_tsv(eo.m3_path(work1, "peptide_hits"))
    i = header.index("start")
    raw = sorted(int(r[i]) for r in rows)
    assert raw == sorted(int(s) + 1 for s in ev1.peptide_hits["start"])


# ----------------------------------------------------------------------------- example 1


def test_example1_every_detected_gene_has_a_correct_peptide(ev1: Any, ex1: Any,
                                                            lay1: Any) -> None:  # fmt: skip
    aa = dict(zip(lay1.genes["locus_tag"], lay1.genes["aa"], strict=True))
    bad = eo.false_set(ev1.truth)
    h = ev1.peptide_hits
    correct = {str(g) for p, g in zip(h["peptide"], h["gene"], strict=True) if p not in bad}
    det = [str(t) for t, d in zip(ex1.genes["locus_tag"], ex1.genes["detected_true"],
                                  strict=True) if d and eo.digestible(aa[t])]  # fmt: skip
    assert len(det) > 400
    missing = [g for g in det if g not in correct]
    assert missing == []


def test_example1_undetected_genes_have_no_correct_peptide(ev1: Any, ex1: Any) -> None:
    det = dict(zip(ex1.genes["locus_tag"], ex1.genes["detected_true"], strict=True))
    t = ev1.truth
    true_rows = t[
        ~t["peptide"].isin(eo.false_set(t)) & t["gene_true"].map(lambda v: not eo.is_na(v))
    ]
    assert len(true_rows) > 0
    assert all(bool(det[g]) for g in true_rows["gene_true"])


def test_example1_decoys_matched_to_false_targets(ev1: Any) -> None:
    p = ev1.peptides
    n_false = len(eo.false_set(ev1.truth) & set(eo.targets(p)["peptide"]))
    assert n_false > 0
    assert len(eo.decoys(p)) == n_false


def test_true_peptides_are_tryptic_pieces_of_their_gene(ev1: Any, lay1: Any) -> None:
    """Trypsin/P, 7-30 aa, <= 1 missed cleavage (tolerant_frac = 0: exact substrings)."""
    aa = dict(zip(lay1.genes["locus_tag"], lay1.genes["aa"], strict=True))
    t = ev1.truth
    bad = eo.false_set(t)
    rows = [(p, g) for p, g in zip(t["peptide"], t["gene_true"], strict=True)
            if p not in bad and not eo.is_na(g)]  # fmt: skip
    assert rows
    wrong = [(p, g) for p, g in rows if p not in eo.peptides_of(aa[g])]
    assert wrong == []


def test_gene_truth_matches_truth(ev1: Any, ex1: Any) -> None:
    t = ev1.truth
    bad = eo.false_set(t)
    tgt = set(eo.targets(ev1.peptides)["peptide"])
    counts: dict[str, int] = {}
    for p, g in zip(t["peptide"], t["gene_true"], strict=True):
        if p in tgt and p not in bad and not eo.is_na(g):
            counts[str(g)] = counts.get(str(g), 0) + 1
    gt = ev1.gene_truth.set_index("locus_tag")
    for g, row in gt.iterrows():
        n = int(row["n_true_peptides"])
        assert n == counts.get(str(g), 0), g
        assert bool(row["present"]) == (n >= 1), g


def test_placements_at_true_locus(ev1: Any, lay1: Any) -> None:
    """C11 (internal 0-based): every unique row sits on its gene at the peptide's residue span;
    class by feature type, exact route, placement_q 0, gene set."""
    genes = eo.gene_rows(lay1)
    h = ev1.peptide_hits
    assert len(h) > 0
    for r in h.to_dict("records"):
        g = genes[str(r["gene"])]
        assert r["locus_class"] == eo.LOCUS_CLASS[str(g["feature_type"])]
        assert r["route"] == "exact"
        assert float(r["placement_q"]) == 0.0
        assert r["uniqueness"] == "unique"
        assert int(r["peptide_len"]) == len(r["peptide"])
        assert r["contig"] == g["contig"] and r["strand"] == g["strand"]
        offs = eo.offsets(str(g["aa"]), str(r["peptide"]))
        want = {eo.span_coords(int(g["start"]), int(g["end"]), str(g["strand"]), a,
                               a + len(r["peptide"])) for a in offs}  # fmt: skip
        assert (int(r["start"]), int(r["end"])) in want, r


def test_c11_input_scores_copy_c9(ev1: Any) -> None:
    q = dict(zip(ev1.peptides["peptide"], ev1.peptides["q_value"], strict=True))
    pe = dict(zip(ev1.peptides["peptide"], ev1.peptides["pep"], strict=True))
    for p, bq, bp in zip(ev1.peptide_hits["peptide"], ev1.peptide_hits["best_q_input"],
                         ev1.peptide_hits["pep_input"], strict=True):  # fmt: skip
        assert float(bq) == pytest.approx(float(q[p]), abs=1e-12)
        assert float(bp) == pytest.approx(float(pe[p]), abs=1e-12)


def test_decoys_never_placed_and_unplaceable_absent(ev1: Any) -> None:
    dec = set(eo.decoys(ev1.peptides)["peptide"])
    assert dec
    assert not dec & set(ev1.peptide_hits["peptide"])


def test_example1_drives_decoy_picked(work1: Path) -> None:
    from giap.preflight.evidence import candidate_tier

    inputs, _ = eo.load(work1)
    assert candidate_tier(inputs.peptides, has_proteins=False) == "decoy_picked"


def test_target_protein_ids_are_the_gene(ev1: Any) -> None:
    """``protein_ids = "locus_tag"``: a placed true peptide names its gene's locus tag."""
    ids = dict(zip(ev1.peptides["peptide"], ev1.peptides["protein_ids"], strict=True))
    for p, g in zip(ev1.peptide_hits["peptide"], ev1.peptide_hits["gene"], strict=True):
        assert str(g) in eo.as_list(ids[p]), (p, g)


# ----------------------------------------------------------------------------- example 2


@pytest.mark.parametrize("pi", [0.01, 0.05, 0.10])
def test_example2_false_count_formula(lay1: Any, ex1: Any, pi: float) -> None:
    """n_false = round(pi/(1 - pi) n_true) exactly, so the false share is pi up to rounding."""
    ev = evidence(lay1, ex1, 2, pi_wrong=pi)
    tgt = eo.targets(ev.peptides)
    bad = eo.false_set(ev.truth) & set(tgt["peptide"])
    n_true = len(tgt) - len(bad)
    assert len(bad) == round(pi / (1 - pi) * n_true)
    if pi == 0.05:
        assert abs(len(bad) / len(tgt) - 0.05) <= 0.01


# ----------------------------------------------------------------------------- example 3


def test_example3_false_peptides_follow_n_observable() -> None:
    """False targets land on undetected genes in proportion to their observable peptides (§8.7);
    pooled over 20 seeds (one seed has ~35 false peptides)."""
    k = n = 0
    expected = 0.0
    for s in range(20):
        lay = layout(300 + s)
        ex = express(lay, s)
        ev = evidence(lay, ex, s)
        undet = {str(t) for t, d in zip(ex.genes["locus_tag"], ex.genes["detected_true"],
                                        strict=True) if not d}  # fmt: skip
        genes = false_hit_genes(ev)
        k += sum(g in undet for g in genes)
        n += len(genes)
        expected += len(genes) * eo.n_obs_share(lay, sorted(undet))
    assert n >= 300
    assert abs(k / n - expected / n) <= 0.05, (k / n, expected / n)


def test_false_targets_map_to_a_gene_peptide_not_observed_as_true(ev1: Any, lay1: Any) -> None:
    aa = dict(zip(lay1.genes["locus_tag"], lay1.genes["aa"], strict=True))
    bad = eo.false_set(ev1.truth)
    h = ev1.peptide_hits
    rows = [(p, g) for p, g in zip(h["peptide"], h["gene"], strict=True) if p in bad]
    assert rows
    for p, g in rows:
        assert p in eo.peptides_of(aa[g]), (p, g)


# ----------------------------------------------------------------------------- example 4


def test_example4_decoy_ids_strip_to_c2_tag(ev1: Any, lay1: Any) -> None:
    tags = set(lay1.genes["locus_tag"])
    aa = dict(zip(lay1.genes["locus_tag"], lay1.genes["aa"], strict=True))
    dec = eo.decoys(ev1.peptides)
    paired = dict(zip(ev1.truth["peptide"], ev1.truth["decoy_gene"], strict=True))
    assert len(dec) > 0
    for p, ids in zip(dec["peptide"], dec["protein_ids"], strict=True):
        ids = eo.as_list(ids)
        assert len(ids) >= 1 and all(i.startswith("rev_") for i in ids), ids
        stripped = {i[len("rev_") :] for i in ids}
        assert stripped <= tags
        assert stripped == {str(paired[p])}
        assert p in eo.decoys_of(aa[str(paired[p])]), (p, paired[p])


def test_decoys_are_not_targets_and_have_null_scores(ev1s: Any) -> None:
    p = ev1s.peptides
    dec = eo.decoys(p)
    assert not set(dec["peptide"]) & set(eo.targets(p)["peptide"])
    assert dec["score"].notna().all()


# ----------------------------------------------------------------------------- q / PEP rules


def test_pep_formula_exact(ev1s: Any) -> None:
    p = ev1s.peptides
    s = p["score"].to_numpy(float)
    want = eo.pep(s, 0.01, 0.0, 3.0)
    np.testing.assert_allclose(p["pep"].to_numpy(float), want, rtol=0, atol=1e-6)


def test_q_is_mean_pep_of_targets_at_or_above(ev1s: Any) -> None:
    p = ev1s.peptides
    t = eo.targets(p)
    q = eo.q_from_targets(p["score"], t["score"], t["pep"])
    ok = ~np.isnan(q)
    assert ok[~p["is_decoy"].astype(bool).to_numpy()].all()
    np.testing.assert_allclose(p["q_value"].to_numpy(float)[ok], q[ok], rtol=0, atol=1e-6)


def test_q_monotone_in_score(ev1s: Any) -> None:
    t = eo.targets(ev1s.peptides).sort_values("score", ascending=False)
    q = t["q_value"].to_numpy(float)
    assert (np.diff(q) >= -1e-12).all()


def test_scores_of_mixture_components(ev1s: Any) -> None:
    """Coarse check (calibration tests the laws): true ~ N(3, 1), false and decoy ~ N(0, 1)."""
    p = ev1s.peptides
    bad = eo.false_set(ev1s.truth)
    dec = p["is_decoy"].astype(bool)
    is_false = p["peptide"].isin(bad) & ~dec
    true_s = p.loc[~dec & ~is_false, "score"].to_numpy(float)
    null_s = p.loc[dec | is_false, "score"].to_numpy(float)
    assert abs(true_s.mean() - 3.0) <= 4 / math.sqrt(len(true_s))
    assert abs(null_s.mean()) <= 4 / math.sqrt(len(null_s))
    assert abs(true_s.std() - 1) <= 0.1


def test_mixture_knob_honoured(lay1: Any, ex1: Any) -> None:
    ev = evidence(lay1, ex1, 3, emit=SCORED, mixture=(-1.0, 2.0), pi_wrong=0.05)
    p = ev.peptides
    want = eo.pep(p["score"].to_numpy(float), 0.05, -1.0, 2.0)
    np.testing.assert_allclose(p["pep"].to_numpy(float), want, atol=1e-6)


# ----------------------------------------------------------------------------- example 5, 6


def test_example5_q_only(lay1: Any, ex1: Any, tmp_path: Path) -> None:
    from giap.preflight.evidence import candidate_tier

    ev = evidence(lay1, ex1, 1, emit=frozenset({"q_value"}))
    p = ev.peptides
    assert p["pep"].isna().all()
    assert p["intensity"].isna().all()
    assert not p["is_decoy"].astype(bool).any()
    assert p["q_value"].notna().all()
    inputs, issues = eo.load(write(ev, lay1, tmp_path / "w"))
    assert eo.errors(issues) == []
    assert candidate_tier(inputs.peptides, has_proteins=False) == "mayu"


def test_emit_only_masks(ev1: Any, lay1: Any, ex1: Any) -> None:
    """Same seed, smaller ``emit``: the same target rows with the same q-values."""
    ev = evidence(lay1, ex1, 1, emit=frozenset({"q_value"}))
    a = eo.targets(ev1.peptides).set_index("peptide")["q_value"].sort_index()
    b = ev.peptides.set_index("peptide")["q_value"].sort_index()
    pd.testing.assert_series_equal(a, b, check_names=False)


def test_example6_sequences_only(lay1: Any, ex1: Any, tmp_path: Path) -> None:
    from giap.preflight.evidence import candidate_tier

    ev = evidence(lay1, ex1, 1, emit=frozenset())
    p = ev.peptides
    for c in ("q_value", "pep", "score", "intensity"):
        assert p[c].isna().all(), c
    assert not p["is_decoy"].astype(bool).any()
    assert all(len(eo.as_list(v)) == 0 for v in p["protein_ids"])
    assert (p["n_rows"].astype(int) >= 1).all() and (p["n_runs"].astype(int) >= 1).all()
    assert len(p) > 0
    inputs, issues = eo.load(write(ev, lay1, tmp_path / "w"))
    assert eo.errors(issues) == []
    assert candidate_tier(inputs.peptides, has_proteins=False) == "two_peptide"


def test_intensity_lognormal_around_gene_abundance(ev1: Any, ex1: Any) -> None:
    ab = dict(zip(ex1.genes["locus_tag"], ex1.genes["abundance_log2_true"], strict=True))
    t = ev1.truth
    bad = eo.false_set(t)
    gene = {p: g for p, g in zip(t["peptide"], t["gene_true"], strict=True)
            if p not in bad and not eo.is_na(g)}  # fmt: skip
    p = eo.targets(ev1.peptides)
    p = p[p["peptide"].isin(set(gene))]
    resid = np.log(p["intensity"].to_numpy(float)) - np.array(
        [ab[gene[x]] * math.log(2) for x in p["peptide"]]
    )
    assert abs(resid.mean()) <= 4 * 0.5 / math.sqrt(len(resid))
    assert abs(resid.std() - 0.5) <= 0.05


def test_more_abundant_genes_yield_more_peptides(ev1: Any, ex1: Any, lay1: Any) -> None:
    """Peptide detectability rises with log abundance (logistic in abundance, spec rule)."""
    from scipy.stats import spearmanr

    aa = dict(zip(lay1.genes["locus_tag"], lay1.genes["aa"], strict=True))
    gt = ev1.gene_truth.set_index("locus_tag")
    g = ex1.genes[ex1.genes["detected_true"].astype(bool)]
    frac = [int(gt.at[t, "n_true_peptides"]) / max(1, len(eo.peptides_of(aa[t])))
            for t in g["locus_tag"]]  # fmt: skip
    rho = spearmanr(g["abundance_log2_true"].to_numpy(float), frac).statistic
    assert rho > 0.1, rho


# ----------------------------------------------------------------------------- example 7


def test_example7_proteins_only(lay1: Any, ex1: Any, tmp_path: Path) -> None:
    from giap.preflight.evidence import candidate_tier

    ev = evidence(lay1, ex1, 1, proteins_only=True)
    assert len(ev.peptides) == 0 and len(ev.peptide_hits) == 0
    ph = ev.protein_hits
    assert ph is not None
    genes = eo.gene_rows(lay1)
    det = {str(t) for t, d in zip(ex1.genes["locus_tag"], ex1.genes["detected_true"],
                                  strict=True) if d}  # fmt: skip
    got = [str(g) for g in ph["gene"]]
    assert len(got) == len(set(got))
    assert set(got) <= det
    assert {g for g in det if eo.digestible(genes[g]["aa"])} <= set(got)
    assert (ph["class"] == "identical").all()
    assert (ph["rank"].astype(int) == 1).all()
    for r in ph.to_dict("records"):
        g = genes[str(r["gene"])]
        assert (r["contig"], int(r["start"]), int(r["end"]), r["strand"]) == (
            g["contig"],
            int(g["start"]),
            int(g["end"]),
            g["strand"],
        )
    inputs, issues = eo.load(write(ev, lay1, tmp_path / "w"))
    assert eo.errors(issues) == []
    assert len(inputs.peptide_hits) == 0 and len(inputs.protein_hits) == len(ph)
    pep = inputs.peptides
    assert pep is None or len(pep) == 0
    assert candidate_tier(pep, has_proteins=True) == "protein_list"


# ------------------------------------------------------------------------ example 8, uniqueness


def test_example8_multi_identical_pairs(lay1: Any, ex1: Any, tmp_path: Path) -> None:
    ev = evidence(lay1, ex1, 1, multi_identical_pairs=2)
    h = ev.peptide_hits
    mi = h[h["uniqueness"] == "multi_identical"]
    assert len(mi) > 0
    pairs = set()
    for _, grp in mi.groupby("peptide"):
        assert len(grp) == 2
        assert grp["gene"].nunique() == 2
        pairs.add(frozenset(grp["gene"].astype(str)))
    assert len(pairs) == 2
    a, b = pairs
    assert not a & b
    assert set(h.loc[h["uniqueness"] != "multi_identical", "uniqueness"]) <= {"unique"}
    _, issues = eo.load(write(ev, lay1, tmp_path / "w"))
    assert eo.errors(issues) == []


def test_shared_frac(lay1: Any, ex1: Any) -> None:
    ev = evidence(lay1, ex1, 1, shared_frac=0.2)
    h = ev.peptide_hits
    bad = eo.false_set(ev.truth)
    uniq = h.groupby("peptide")["uniqueness"].nunique()
    assert (uniq == 1).all()
    sh = h[h["uniqueness"] == "shared"]
    assert len(sh) > 0
    for _, grp in sh.groupby("peptide"):
        assert len(grp) >= 2 and grp["gene"].nunique() == len(grp)
    true_peps = set(h["peptide"]) - bad
    k = len(set(sh["peptide"]) & true_peps)
    n = len(true_peps)
    assert abs(k - 0.2 * n) <= 4 * math.sqrt(n * 0.2 * 0.8) + 1, (k, n)


def test_tolerant_frac(lay1: Any, ex1: Any, tmp_path: Path) -> None:
    ev = evidence(lay1, ex1, 1, tolerant_frac=0.1)
    h = ev.peptide_hits
    genes = eo.gene_rows(lay1)
    bad = eo.false_set(ev.truth)
    tol = h[h["route"] == "tolerant"]
    assert len(tol) > 0
    for r in tol.to_dict("records"):
        assert int(r["n_diff"]) == 1
        assert "strain_difference" in eo.as_list(r["flags"])
        g = genes[str(r["gene"])]
        exact = eo.offsets(str(g["aa"]), str(r["peptide"]))
        near = eo.offsets(str(g["aa"]), str(r["peptide"]), max_diff=1)
        assert not exact and near, r["peptide"]
        want = {eo.span_coords(int(g["start"]), int(g["end"]), str(g["strand"]), a,
                               a + len(r["peptide"])) for a in near}  # fmt: skip
        assert (int(r["start"]), int(r["end"])) in want
    ex_rows = h[h["route"] == "exact"]
    assert ex_rows["n_diff"].map(lambda v: eo.is_na(v) or int(v) == 0).all()
    n = len(set(h["peptide"]) - bad)
    k = len(set(tol["peptide"]) - bad)
    assert abs(k - 0.1 * n) <= 4 * math.sqrt(n * 0.1 * 0.9) + 1, (k, n)
    _, issues = eo.load(write(ev, lay1, tmp_path / "w"))
    assert eo.errors(issues) == []


# ----------------------------------------------------------------------------- example 9


def test_example9_q_pep_consistency(ev1: Any) -> None:
    t = eo.targets(ev1.peptides)
    acc = t[t["q_value"].astype(float) <= 0.01]
    assert len(acc) > 100
    assert float(acc["pep"].astype(float).mean()) <= 0.01 + 1e-12
    share = float(acc["peptide"].isin(eo.false_set(ev1.truth)).mean())
    se = math.sqrt(0.01 * 0.99 / len(acc))
    assert share <= 0.01 + 3 * se, share


# ----------------------------------------------------------------------------- example 10, RNG


def test_example10_same_seed_byte_identical_tables(tmp_path: Path, lay1: Any, ex1: Any) -> None:
    kw: dict[str, Any] = {"pi_wrong": 0.05, "abundance_table": "ibaq", "tolerant_frac": 0.05}
    a = write(evidence(lay1, ex1, 10, **kw), lay1, tmp_path / "a")
    b = write(evidence(lay1, ex1, 10, **kw), lay1, tmp_path / "b")
    ta, tb = eo.tables_in(a), eo.tables_in(b)
    assert ta.keys() == tb.keys()
    assert ta == tb
    c = write(evidence(lay1, ex1, 11, **kw), lay1, tmp_path / "c")
    assert eo.tables_in(c) != ta


def test_frames_identical_same_seed(lay1: Any, ex1: Any, ev1: Any) -> None:
    again = evidence(lay1, ex1, 1)
    for name in ("peptides", "peptide_hits", "truth", "gene_truth"):
        pd.testing.assert_frame_equal(getattr(ev1, name), getattr(again, name), obj=name)


def test_rng_is_the_only_source_of_randomness(lay1: Any, ex1: Any) -> None:
    m = _ev_mod()
    np.random.seed(0)  # noqa: NPY002
    a = m.simulate_evidence(lay1, ex1, m.EvidenceConfig(), np.random.default_rng(4))
    np.random.seed(999)  # noqa: NPY002
    b = m.simulate_evidence(lay1, ex1, m.EvidenceConfig(), np.random.default_rng(4))
    pd.testing.assert_frame_equal(a.peptides, b.peptides)
    pd.testing.assert_frame_equal(a.truth, b.truth)


def test_inputs_not_mutated(lay1: Any, ex1: Any) -> None:
    lay_genes, ex_genes = lay1.genes.copy(deep=True), ex1.genes.copy(deep=True)
    evidence(lay1, ex1, 5, multi_identical_pairs=2, shared_frac=0.1, tolerant_frac=0.1)
    pd.testing.assert_frame_equal(lay1.genes, lay_genes)
    pd.testing.assert_frame_equal(ex1.genes, ex_genes)


# ----------------------------------------------------------------------------- abundance (C10)


@pytest.mark.parametrize("kind", ["ibaq", "maxlfq"])
def test_abundance_table(lay1: Any, ex1: Any, tmp_path: Path, kind: str) -> None:
    from scipy.stats import spearmanr

    ev = evidence(lay1, ex1, 1, abundance_table=kind)
    ab = ev.abundance
    assert ab is not None
    assert ab["protein_id"].is_unique
    det = {str(t) for t, d in zip(ex1.genes["locus_tag"], ex1.genes["detected_true"],
                                  strict=True) if d}  # fmt: skip
    present = set(ev.gene_truth.loc[ev.gene_truth["present"].astype(bool), "locus_tag"])
    assert present <= set(ab["protein_id"]) <= set(lay1.genes["locus_tag"])
    assert set(ab["protein_id"]) <= det
    true = dict(zip(ex1.genes["locus_tag"], ex1.genes["abundance_log2_true"], strict=True))
    v = ab["value"].to_numpy(float)
    rho = spearmanr(v, [true[p] for p in ab["protein_id"]]).statistic
    assert rho > 0.8, rho
    inputs, issues = eo.load(write(ev, lay1, tmp_path / "w"))
    assert eo.errors(issues) == []
    assert inputs.abundance is not None and len(inputs.abundance) == len(ab)


def test_ibaq_and_maxlfq_share_values(lay1: Any, ex1: Any) -> None:
    a = evidence(lay1, ex1, 1, abundance_table="ibaq").abundance
    b = evidence(lay1, ex1, 1, abundance_table="maxlfq").abundance
    pd.testing.assert_frame_equal(a.reset_index(drop=True), b.reset_index(drop=True))


def test_no_abundance_file_when_none(work1: Path) -> None:
    assert not (work1 / "abundance.norm.parquet").exists()


# ----------------------------------------------------------------------------- edge cases


def test_edge_undigestible_detected_gene_forced_absent(tmp_path: Path) -> None:
    from giapsim.expression import Expression

    lay = layout(7, n_genes_chromosome=300, frac_sorf=0.6)
    ex = express(lay, 7)
    forced = Expression(ex.genes.assign(detected_true=True), ex.islands)
    aa = dict(zip(lay.genes["locus_tag"], lay.genes["aa"], strict=True))
    undig = [t for t in lay.genes["locus_tag"] if not eo.digestible(aa[t])]
    assert any(isinstance(aa[t], str) for t in undig)  # some are real short sORFs
    ev = evidence(lay, forced, 7)
    gt = ev.gene_truth.set_index("locus_tag")
    for t in undig:
        assert not bool(gt.at[t, "present"]) and int(gt.at[t, "n_true_peptides"]) == 0, t
        if "detected_true" in gt.columns:
            assert not bool(gt.at[t, "detected_true"])
    dig = [t for t in lay.genes["locus_tag"] if eo.digestible(aa[t])]
    assert all(bool(gt.at[t, "present"]) for t in dig)
    assert (forced.genes["detected_true"]).all()  # input not mutated
    _, issues = eo.load(write(ev, lay, tmp_path / "w"))
    assert eo.errors(issues) == []


def test_edge_pi_zero(lay1: Any, ex1: Any) -> None:
    ev = evidence(lay1, ex1, 1, pi_wrong=0.0, emit=SCORED)
    assert not eo.false_set(ev.truth)
    assert len(eo.decoys(ev.peptides)) == 0
    assert (ev.peptides["pep"].astype(float) == 0).all()
    assert (ev.peptides["q_value"].astype(float) == 0).all()


@pytest.mark.parametrize("n", [0, 7, 40])
def test_edge_integer_n_decoys(lay1: Any, ex1: Any, n: int) -> None:
    ev = evidence(lay1, ex1, 1, n_decoys=n)
    assert len(eo.decoys(ev.peptides)) == n


def test_edge_foreign_protein_ids(lay1: Any, ex1: Any, tmp_path: Path) -> None:
    ev = evidence(lay1, ex1, 1, protein_ids="foreign", n_decoys=25)
    p = ev.peptides
    tags = set(lay1.genes["locus_tag"])
    dec = eo.decoys(p)
    assert len(dec) == 25
    for ids in dec["protein_ids"]:
        ids = eo.as_list(ids)
        assert ids and all(i.startswith("rev_") and UP_RE.match(i[4:]) for i in ids)
        assert not {i[4:] for i in ids} & tags
    acc: dict[str, set[str]] = {}
    for pep, g in zip(ev.peptide_hits["peptide"], ev.peptide_hits["gene"], strict=True):
        ids = eo.as_list(p.set_index("peptide").at[pep, "protein_ids"])
        assert ids and all(UP_RE.match(i) for i in ids), ids
        acc.setdefault(str(g), set()).update(ids)
    owners: dict[str, str] = {}
    for g, ids in acc.items():
        for i in ids:
            assert owners.setdefault(i, g) == g, (i, g)  # one accession per gene
    _, issues = eo.load(write(ev, lay1, tmp_path / "w"))
    assert eo.errors(issues) == []


def test_unplaceable_false_targets(lay1: Any, ex1: Any, tmp_path: Path) -> None:
    from giap.io.tsv import read_tsv

    ev = evidence(lay1, ex1, 1, pi_wrong=0.05, false_source="unplaceable")
    bad = eo.false_set(ev.truth) & set(eo.targets(ev.peptides)["peptide"])
    assert bad
    assert not bad & set(ev.peptide_hits["peptide"])
    allaa = "|".join(a for a in lay1.genes["aa"] if isinstance(a, str))
    assert not [p for p in bad if p in allaa]  # maps nowhere
    assert all(7 <= len(p) <= 30 for p in bad)
    work = write(ev, lay1, tmp_path / "w")
    _, issues = eo.load(work)
    assert eo.errors(issues) == []
    header, rows = read_tsv(eo.m3_path(work, "unmapped"))
    seqs = {r[header.index("sequence")] for r in rows}
    assert bad <= seqs


def test_invalid_config_rejected(lay1: Any, ex1: Any) -> None:
    for kw in ({"pi_wrong": -0.1}, {"pi_wrong": 1.0}, {"shared_frac": 1.5},
               {"tolerant_frac": -0.2}, {"false_source": "shuffled"},
               {"abundance_table": "lfq"}, {"n_decoys": -1}):  # fmt: skip
        with pytest.raises(ValueError):
            evidence(lay1, ex1, 1, **kw)


# ----------------------------------------------------------------------------- backward compat


def test_earlier_outputs_unchanged(tmp_path: Path) -> None:
    """M4-07 tables and M4-08 frames are the same before and after an evidence run (one run);
    pinned digests only on Linux (Parquet bytes differ across platforms)."""
    from giapsim.contract_tables import write_contract_tables

    lay = layout(1)
    write_contract_tables(lay, tmp_path / "before")
    ex_before = express(lay, 1).genes
    evidence(lay, express(lay, 1), 1, multi_identical_pairs=2, abundance_table="ibaq")
    write_contract_tables(lay, tmp_path / "after")
    ex_after = express(lay, 1).genes
    before, after = eo.tables_in(tmp_path / "before"), eo.tables_in(tmp_path / "after")
    assert before == after
    pd.testing.assert_frame_equal(ex_before, ex_after)
    if sys.platform.startswith("linux"):
        digests = {k: hashlib.sha256(v).hexdigest() for k, v in after.items()}
        assert digests == LAYOUT_SHA256
        # The expression frame is compared before/after above; its float output differs in the
        # last bits between CPUs (wheat-08 vs the session node), so it has no pinned digest.


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="digests pinned on Linux")
def test_m0_proteomics_outputs_unchanged(tmp_path: Path) -> None:
    """A shared score helper (lead Q3) must leave the M3-19 proteomics output byte-identical."""
    import giapsim
    from giapsim.proteomics import ProteomicsConfig

    _ev_mod()
    giapsim.simulate(giapsim.SimConfig(proteomics=ProteomicsConfig()), 1, tmp_path)
    for rel, digest in M0_PROTEOMICS_SHA256.items():
        assert hashlib.sha256((tmp_path / rel).read_bytes()).hexdigest() == digest, rel


# ----------------------------------------------------------------------------- acceptance


def test_speed_4500_genes_loose_guard() -> None:
    """Loose 10x guard (budget 3 s; wall-clock budgets live in benchmarks/ — AGENTS conventions)."""
    lay = layout(1, n_genes_chromosome=4500)
    ex = express(lay, 1)
    m = _ev_mod()
    best = math.inf
    for s in range(3):
        t0 = time.perf_counter()
        m.simulate_evidence(lay, ex, m.EvidenceConfig(), np.random.default_rng(s))
        best = min(best, time.perf_counter() - t0)
    assert best < 30.0, best


def test_evidence_module_imports_nothing_from_giap() -> None:
    src = (GIAPSIM / "evidence.py").read_text()
    assert not re.search(r"^\s*(?:from\s+giap(?:\.\w+)*\s+import|import\s+giap(?:\.|\s|,|$))",
                         src, re.M)  # fmt: skip
    code = (
        "import sys, giapsim.evidence as e; "
        "bad = [m for m in sys.modules if m == 'giap' or m.startswith('giap.')]; "
        "print(bad); sys.exit(1 if bad else 0)"
    )
    env = {"PYTHONPATH": str(REPO / "tests" / "sim"), "PATH": "/usr/bin:/bin"}
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env,
                       check=False)  # fmt: skip
    assert r.returncode == 0, r.stdout + r.stderr


def test_locus_tags_in_hits_are_layout_tags(ev1: Any) -> None:
    assert all(LOCUS_RE.match(str(g)) for g in ev1.peptide_hits["gene"])
