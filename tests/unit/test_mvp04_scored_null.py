"""giapsim on the real K-12 genome, part 2: ``decoys = "scored_null"``
(docs/specs/MVP-04-ecoli-evidence-run.md pipeline steps 5-6, examples 2-3, D6, D7, D11, Q2, lead
decisions 1 and 8; M7-12 decisions 21, 23 and 24 for the six-frame coin this mode extends).

On the test-time 200-kb excerpt (tests/helpers/mvp04_data.py, spec rule 3):

* rule 6: each null (wrong) PSM goes, on its own 1/2 coin, to the target or to the decoy of its
  entry (``rev_<locus_tag>`` for a protein entry, ``rev_GIAPORF_…`` for a six-frame entry), with
  the reversed peptide, the same q / PEP / score and intensity NA; the M3-19 ``match_wrong``
  decoys are off, so every decoy is a null PSM and carries finite scores;
* collisions (reversal = a target peptide, I = L) stay targets and are counted in
  ``provenance.json`` ``decoy_collisions``;
* Q2 / decision 1: null peptides split between protein and six-frame entries in proportion to the
  tryptic products of the simulated search space, 50/50 when that gives fewer than about 50
  gene-level null PSMs (the excerpt is in the fallback: ~0.14 x 18 x 3 ~ 8 PSMs);
* a fast fairness check of the coin (binomial, and a score-independence KS), pooled over three
  seeds (~160 null PSMs); the calibration file owns the powered version;
* GIAP's own ``pair_decoys`` + ``choose_tier`` give ``decoy_picked``, no warning, paired share
  >= 0.98 (``rev_GIAPORF_*`` leave pairing, M7-12 decision 23) and at least one decoy win over
  the three seeds (gate 5's shape; expected ~3 per seed);
* lead decision 9: protein-entry nulls are undetected real products of their locus, map to the
  genome there and give genes without a true peptide false target evidence, so gate 7's
  false-discovery share is estimable (non-zero over three seeds);
* the default ``decoys = "unscored"`` reproduces the part-1 tables byte for byte.

Assumed interface (spec "Interface"; names the implementer may refine within giapsim, see the
report on the commit):

* ``ProteomicsConfig.decoys: Literal["unscored", "scored_null"] = "unscored"``; another value is
  a ``ValueError``; ``scenarios.k12_evidence().proteomics.decoys == "scored_null"``.
* In ``scored_null`` mode the null peptides are drawn from both sources regardless of
  ``wrong_source``. Lead decision 9: a protein-entry null is a real Trypsin/P product (digest
  rules of the config) of an expressible recorded locus with no true peptide, the locus drawn
  uniformly over the loci with an undetected product; no shuffled peptides. Its truth row has
  ``is_wrong = true``, ``origin = undetected_product`` and is placed like a correct row
  (``locus_tag`` = the attributed locus, ``contig``, 1-based closed ``start`` / ``end``,
  ``strand`` on the sample); its target rows carry that locus tag and its decoy is
  ``rev_<locus>``. Six-frame nulls have ``origin = sixframe_noncoding`` and their target rows
  carry the makedb entry ID ``GIAPORF_…`` (step 5: "attributed to a ``GIAPORF_*`` entry").
* ``provenance.json``: ``decoys`` (``"scored_null"`` / ``"unscored"``), ``decoy_collisions``
  (int, collided null PSMs) and ``null_split`` = {``mode``: ``"search_space"`` |
  ``"fifty_fifty"``, ``search_space_protein_share``: float in (0, 1), ``n_protein``: int,
  ``n_sixframe``: int (distinct null peptides of each kind)}.
"""

from __future__ import annotations

import hashlib
import importlib
import math
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest
from helpers import mvp04_data as M
from helpers import mvp04_decoys as D
from scipy import stats

from giapsim.expression import IslandSpec

SEED = 20261023
SEEDS = (SEED, SEED + 1, SEED + 2)
N_NOVEL = 20
ISLANDS = ((2, 3.0), (18, -3.0), (34, 3.0), (50, -3.0), (66, 3.0), (82, -3.0))
ALPHA = 0.001

#: sha256 of the part-1 excerpt tables (head 693c7b2a, before part 2): seed 20261023, gap rule,
#: N_NOVEL = 20, ISLANDS, k12_evidence() otherwise. The ``unscored`` mode must reproduce them.
PART1_SHA = {
    "peptides.tsv": "0a9e7d75f3db2daaef7614e907ce9006223e5dbf8a6efbc04e78f315cba6d764",
    "abundance.tsv": "b2bbad9f8303aeaab2c27b3bd5b9ee1c65349678150468004168a26fb736f3b7",
    "truth/alt_sites.tsv": "a19a67ded082e119a3d3c0744f9995a8d40344806321a74ef337e7e3f3b1c31b",
    "truth/genes.tsv": "0fae7ed0e9d2c3cf0a586240a3a3ac2613d1ea230e03f08231d821a0fa5088d3",
    "truth/islands.tsv": "d0a489cfc63cc6128e3d46c88697449ccb6ea21700f0e9fbbc1c1eb43271279a",
    "truth/novel_orfs.tsv": "1bdede5b9236887381025ed93d5a830152875ca0eebd3757e2172ded5793b7ac",
    "truth/peptides.tsv": "1873828fe94b3164275bc602e414b0ebaee9b037ccffe9306a490283347f81d9",
}


def rg() -> Any:
    return importlib.import_module("giapsim.real_genome")


def excerpt_cfg(decoys: str | None = "scored_null", **prot: Any) -> Any:
    base = importlib.import_module("giapsim.scenarios").k12_evidence()
    expr = replace(
        base.expression, islands=tuple(IslandSpec(M.CONTIG, f, 8, e) for f, e in ISLANDS)
    )
    if decoys is not None:
        prot["decoys"] = decoys
    return replace(base, expression=expr, n_novel_per_class=N_NOVEL,
                   proteomics=replace(base.proteomics, **prot))  # fmt: skip


def run(ex: M.Excerpt, out: Path, cfg: Any, seed: int = SEED) -> Path:
    rg().simulate_on_genome(cfg, seed, out, genome_fasta=ex.fasta, bakta_json=ex.bakta_json,
                            tu_set=None, check_hashes=False)  # fmt: skip
    return out


@pytest.fixture(scope="module")
def ex(tmp_path_factory: pytest.TempPathFactory) -> M.Excerpt:
    return M.make_excerpt(tmp_path_factory.mktemp("k12_excerpt_p2"))


@pytest.fixture(scope="module")
def scored(ex: M.Excerpt, tmp_path_factory: pytest.TempPathFactory) -> dict[int, D.View]:
    """``scored_null`` excerpt runs at the three ``SEEDS``."""
    return {s: D.view(run(ex, tmp_path_factory.mktemp(f"scored_{s}"), excerpt_cfg(), s))
            for s in SEEDS}  # fmt: skip


@pytest.fixture(scope="module")
def v0(scored: dict[int, D.View]) -> D.View:
    return scored[SEED]


@pytest.fixture(scope="module")
def recorded(ex: M.Excerpt) -> dict[str, dict[str, Any]]:
    return D.annotation(ex.bakta_json)


@pytest.fixture(scope="module")
def space(ex: M.Excerpt, v0: D.View, tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """``giap makedb --six-frame-orfs`` on the excerpt sample (the simulated search, step 5)."""
    tmp = tmp_path_factory.mktemp("makedb_p2")
    fna = D.write_sample_fasta(v0.out, tmp / "sample.fna")
    return D.makedb_search_space(ex.bakta_json, fna, tmp)


def _mix() -> tuple[float, ...]:
    return tuple(float(x) for x in excerpt_cfg().proteomics.score_mixture)


# ---------------------------------------------------------------- config


def test_decoys_field_defaults_to_unscored() -> None:
    from giapsim.proteomics import ProteomicsConfig

    assert ProteomicsConfig().decoys == "unscored"
    assert ProteomicsConfig(decoys="scored_null").decoys == "scored_null"
    with pytest.raises(ValueError):
        ProteomicsConfig(decoys="scored")


def test_k12_evidence_uses_scored_null_decoys() -> None:
    """Spec "Interface": the evidence scenario's proteomics has ``decoys="scored_null"``."""
    cfg = importlib.import_module("giapsim.scenarios").k12_evidence()
    assert cfg.proteomics.decoys == "scored_null"


# ---------------------------------------------------------------- the default stays byte-identical


def test_unscored_reproduces_part1_byte_for_byte(ex: M.Excerpt, tmp_path: Path) -> None:
    """Every existing default stays byte-identical (spec "giapsim on a real genome")."""
    out = run(ex, tmp_path / "unscored", excerpt_cfg(decoys="unscored"))
    got = {k: hashlib.sha256((out / k).read_bytes()).hexdigest() for k in PART1_SHA}
    assert got == PART1_SHA
    prov = D.view(out).provenance
    assert prov["decoys"] == "unscored"
    assert prov["decoy_collisions"] == 0


# ---------------------------------------------------------------- rule 6: the coin


def test_every_null_psm_is_written_once(scored: dict[int, D.View]) -> None:
    """Rule 6: target xor decoy per (null peptide, run); "other runs draw their own coin"."""
    for seed, v in scored.items():
        assert len(v.runs) == 3, seed
        assert v.nulls, seed
        bad = D.accounting_violations(v)
        assert not bad, (seed, bad[:10])


def test_both_sides_of_the_coin_occur(scored: dict[int, D.View]) -> None:
    psms = [p for v in scored.values() for p in D.null_psms(v)]
    for kind in ("protein", "sixframe"):
        sides = {is_dec for n, _, _, is_dec in psms if n.kind == kind}
        assert sides == {False, True}, kind


def test_match_wrong_decoys_are_off(v0: D.View) -> None:
    """No reversed-random-protein decoys: the decoys are exactly the null PSMs that went to the
    decoy side, none without a score (rule 6, "No extra decoys")."""
    n_dec = sum(len(n.decoys) for n in v0.nulls.values())
    assert len(v0.decoy_rows) == n_dec
    n_null_psm = sum(len(n.decoys) + len(n.targets) for n in v0.nulls.values())
    assert n_null_psm == len(v0.nulls) * len(v0.runs)


def test_decoys_carry_finite_q_pep_score_and_na_intensity(scored: dict[int, D.View]) -> None:
    for seed, v in scored.items():
        assert v.decoy_rows, seed
        for r in v.decoy_rows:
            for c in ("q_value", "pep", "score"):
                assert D.finite(r[c]), (seed, c, r)
            assert r["intensity"] == "NA", (seed, r)
            assert r["source"] == "database", (seed, r)


def test_decoy_and_null_target_scores_come_from_one_mixture(v0: D.View) -> None:
    """Example 2: the decoy keeps q, PEP and score of its null draw, so q and PEP are the
    mixture's at that score and at the pre-coin wrong share ``pi_wrong_realised`` (rule 6, "Known
    bias"), for decoy and null-target rows alike."""
    pi = float(v0.meta["pi_wrong_realised"])
    mix = _mix()
    for n, run_, r, is_dec in D.null_psms(v0):
        s = float(r["score"])
        assert math.isclose(float(r["pep"]), D.mixture_pep(s, pi, mix), rel_tol=1e-6,
                            abs_tol=1e-12), (n.peptide, run_, is_dec, r)  # fmt: skip
        assert math.isclose(float(r["q_value"]), D.mixture_q(s, pi, mix), rel_tol=1e-6,
                            abs_tol=1e-12), (n.peptide, run_, is_dec, r)  # fmt: skip


def test_protein_entry_decoys_are_rev_locus_of_an_annotated_gene(
    scored: dict[int, D.View], recorded: dict[str, dict[str, Any]]
) -> None:
    """Example 2: ``rev_<locus_tag>`` of an expressible recorded protein gene (no ``rev_NOVEL_*``,
    M7-12 decision 23); the null target rows name the same locus."""
    n_seen = 0
    for seed, v in scored.items():
        for n in v.nulls.values():
            if n.kind != "protein":
                continue
            (entry,) = n.entries()
            assert entry in recorded and D.expressible(recorded[entry]), (seed, n.peptide, entry)
            for r in n.targets.values():
                assert r["protein_ids"] == entry, (seed, r)
            for r in n.decoys.values():
                assert r["protein_ids"] == f"rev_{entry}", (seed, r)
                assert r["peptide"] == n.peptide[::-1]
                n_seen += 1
    assert n_seen > 0


def test_sixframe_decoys_are_rev_of_their_makedb_entry(v0: D.View, space: dict[str, Any]) -> None:
    """Example 3: ``rev_GIAPORF_{class}_{contig}_{start1}_{end1}_{p|m}``; the entry is a
    ``giap makedb --six-frame-orfs`` entry of the sample whose protein holds the null peptide,
    and the null target rows carry the entry ID itself (step 5)."""
    entries = space["entries"]
    n_seen = 0
    for n in v0.nulls.values():
        if n.kind != "sixframe":
            continue
        (entry,) = n.entries()
        assert D.ORF_ID.match(entry), (n.peptide, entry)
        assert entry in entries, (n.peptide, entry)
        assert n.peptide in entries[entry], (n.peptide, entry)
        for r in n.targets.values():
            assert r["protein_ids"] == entry, r
        for r in n.decoys.values():
            assert r["protein_ids"] == f"rev_{entry}", r
            assert D.ORF_DECOY.match(r["protein_ids"]), r
            n_seen += 1
    assert n_seen > 0


def test_no_decoy_collides_with_a_target_and_collisions_are_counted(
    scored: dict[int, D.View],
) -> None:
    """Rule 6 "Collisions" / D6: a reversed peptide equal to a target (I = L) stays a target; the
    count goes to ``decoy_collisions``."""
    for seed, v in scored.items():
        targets = {D.il(r["peptide"]) for r in v.target_rows} | {D.il(p) for p in v.correct}
        for r in v.decoy_rows:
            assert D.il(r["peptide"]) not in targets, (seed, r)
        colliders = [n for n in v.nulls.values() if D.il(n.peptide[::-1]) in targets]
        k = v.provenance["decoy_collisions"]
        assert isinstance(k, int) and k >= 0, (seed, k)
        assert k <= len(colliders) * len(v.runs), (seed, k, len(colliders))


# ---------------------------------------------------------------- Q2: the split


def test_split_falls_back_to_fifty_fifty_on_the_excerpt(scored: dict[int, D.View]) -> None:
    """Q2 / decision 1: the search-space share (~0.14) gives ~8 gene-level null PSMs here, below
    about 50, so the split is 50/50 and the switch is recorded."""
    for seed, v in scored.items():
        split = v.provenance["null_split"]
        assert split["mode"] == "fifty_fifty", (seed, split)
        n_p = sum(n.kind == "protein" for n in v.nulls.values())
        n_s = sum(n.kind == "sixframe" for n in v.nulls.values())
        assert (split["n_protein"], split["n_sixframe"]) == (n_p, n_s), (seed, split)
        expected = split["search_space_protein_share"] * (n_p + n_s) * len(v.runs)
        assert expected < D.FALLBACK_PSMS, (seed, expected)
        assert n_p > 0 and n_s > 0, (seed, n_p, n_s)
    n_p = sum(v.provenance["null_split"]["n_protein"] for v in scored.values())
    n = n_p + sum(v.provenance["null_split"]["n_sixframe"] for v in scored.values())
    assert stats.binomtest(n_p, n, 0.5).pvalue >= ALPHA, (n_p, n)


def test_search_space_share_matches_makedb(v0: D.View, space: dict[str, Any]) -> None:
    """Q2: the recorded share is that of the tryptic products of protein entries among those of
    the simulated makedb search space (protein + six-frame entries)."""
    got = float(v0.provenance["null_split"]["search_space_protein_share"])
    assert 0.0 < got < 1.0
    assert abs(got - space["protein_share"]) <= 0.04, (got, space["protein_share"])


# ---------------------------------------------------------------- the coin is fair (fast check)


def test_coin_is_fair_and_independent_of_the_score(scored: dict[int, D.View]) -> None:
    """Pooled over three seeds (~160 null PSMs): D among D + T is Binomial(1/2) (two-sided exact
    test, p >= 0.001), and the decoy and null-target scores are one distribution (two-sample KS,
    p >= 0.001), so the coin does not look at the score."""
    psms = [p for v in scored.values() for p in D.null_psms(v)]
    d = [float(r["score"]) for _, _, r, is_dec in psms if is_dec]
    t = [float(r["score"]) for _, _, r, is_dec in psms if not is_dec]
    assert len(d) + len(t) >= 100, (len(d), len(t))
    assert stats.binomtest(len(d), len(d) + len(t), 0.5).pvalue >= ALPHA, (len(d), len(t))
    assert stats.ks_2samp(d, t).pvalue >= ALPHA


# ---------------------------------------------------------------- GIAP's tier on the output


def test_giap_tier_is_decoy_picked_with_full_pairing_and_decoy_wins(
    ex: M.Excerpt, scored: dict[int, D.View]
) -> None:
    """Acceptance 2, 3, 5 in shape (example 4 contrast): ``pair_decoys`` keeps the
    ``rev_<locus>`` decoys at q <= 0.01 and drops the ``rev_GIAPORF_*`` ones; ``choose_tier`` is
    ``decoy_picked`` with no warning and f >= 0.98; over the three seeds at least one paired
    decoy beats its gene's best target (not the D = 0 picture of 1869554)."""
    genes = ex.genes
    tags = [g.locus for g in genes]
    eligible = {g.locus for g in genes if g.eligible}
    proteins = {g.locus: g.aa for g in genes}
    wins = 0
    for seed, v in scored.items():
        t = D.giap_tier(v, tags)
        assert t.tier == "decoy_picked", (seed, t)
        assert t.warnings == (), (seed, t)
        assert t.n_pairing >= 1, seed
        assert t.paired_share >= 0.98, (seed, t.paired_share)
        prot = [r for n in v.nulls.values() if n.kind == "protein" for r in n.decoys.values()]
        n_prot = sum(float(r["q_value"]) <= D.Q_VALUE_MAX for r in prot)
        assert t.n_pairing == n_prot, (seed, t.n_pairing, n_prot)
        wins += sum(D.decoy_wins(v, proteins, eligible).values())
    assert wins >= 1


# ---------------------------------------------------------------- decision 9: protein nulls


def _protein_null_bad(ex: M.Excerpt, v: D.View) -> dict[str, list[str]]:
    recorded = D.annotation(ex.bakta_json)
    sample = D.read_seq(v.out / "truth" / "sample_genome.fna")
    proteins = D.locus_proteins(v.out, recorded, sample)
    return D.protein_null_violations(v, recorded, proteins, sample, ex.seq)


@pytest.fixture(scope="module")
def null_bad(ex: M.Excerpt, scored: dict[int, D.View]) -> dict[int, dict[str, list[str]]]:
    return {s: _protein_null_bad(ex, v) for s, v in scored.items()}


def test_protein_nulls_are_labelled_undetected_product(scored: dict[int, D.View]) -> None:
    """Decision 9: no shuffled peptides in this mode; the protein-entry origin label."""
    for seed, v in scored.items():
        origins = {n.origin for n in v.nulls.values()}
        assert origins == {D.PROTEIN_NULL_ORIGIN, "sixframe_noncoding"}, (seed, origins)


def test_protein_null_is_a_product_of_its_locus_on_the_sample(
    null_bad: dict[int, dict[str, list[str]]],
) -> None:
    """Decision 9: the target peptide is a Trypsin/P product (7-30 aa, <= 1 missed cleavage) of
    the sample protein of the locus it is attributed to, an expressible recorded locus."""
    for seed, bad in null_bad.items():
        assert not bad["product"], (seed, bad["product"][:5])


def test_protein_null_locus_has_no_true_peptide(
    null_bad: dict[int, dict[str, list[str]]],
) -> None:
    for seed, bad in null_bad.items():
        assert not bad["undetected"], (seed, bad["undetected"][:5])


def test_protein_null_target_maps_to_the_genome(
    null_bad: dict[int, dict[str, list[str]]],
) -> None:
    """Decision 9: its truth row places it on the attributed locus; the sample segment
    translates to it and the K-12 sequence GIAP is given differs by SNVs only."""
    for seed, bad in null_bad.items():
        assert not bad["maps"], (seed, bad["maps"][:5])


def test_gate7_false_discovery_share_is_estimable(ex: M.Excerpt, scored: dict[int, D.View]) -> None:
    """Decision 9 / acceptance 7: genes with no true peptide now get null target evidence that
    can win the picked competition (expected ~ 0.5 x 7/8 per protein null peptide, ~ 4 per
    seed here), so the false-discovery share among target-won genes is non-zero over three
    seeds, while true genes still dominate."""
    eligible = {g.locus for g in ex.genes if g.eligible}
    false_won = true_won = with_null = 0
    for v in scored.values():
        comp = D.gene_competition(v, eligible)
        false_won += sum(c.target_won and not c.has_true for c in comp)
        true_won += sum(c.target_won and c.has_true for c in comp)
        with_null += sum(not c.has_true for c in comp)
    assert with_null >= 5, with_null
    assert false_won >= 1, (false_won, with_null)
    assert 0 < false_won / (false_won + true_won) < 0.5, (false_won, true_won)


# ---------------------------------------------------------------- provenance and TE7


def test_provenance_records_scored_null(v0: D.View) -> None:
    prov = v0.provenance
    assert prov["decoys"] == "scored_null"
    assert isinstance(prov["decoy_collisions"], int)
    assert set(prov["null_split"]) >= {"mode", "search_space_protein_share", "n_protein",
                                       "n_sixframe"}  # fmt: skip


def test_scored_null_same_seed_is_byte_identical(ex: M.Excerpt, v0: D.View, tmp_path: Path) -> None:
    """Acceptance 9 (TE7 on the inputs) in the mode the evidence run uses."""
    again = run(ex, tmp_path / "again", excerpt_cfg())
    for name in ("peptides.tsv", "abundance.tsv", *PART1_SHA):
        assert (again / name).read_bytes() == (v0.out / name).read_bytes(), name
