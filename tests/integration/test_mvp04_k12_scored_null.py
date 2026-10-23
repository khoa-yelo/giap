"""MVP-04 part 2 on the full K-12 genome: ``decoys = "scored_null"`` at ``K12_EVIDENCE_SEED``
(docs/specs/MVP-04-ecoli-evidence-run.md steps 5-6, Q2, acceptance 2-5, lead decisions 1 and 8).

Slow: one whole-genome ``simulate_on_genome(k12_evidence(), K12_EVIDENCE_SEED)`` from the
committed FASTA and recorded Bakta JSON (hash-checked), outputs in pytest temporary directories.
RegulonDB TUs when ``$GIAP_DATA_DIR/regulondb/14.5.0/TUSet.txt`` exists, else the gap rule (the
decoy checks do not depend on the TUs). Interface: as tests/unit/test_mvp04_scored_null.py.

Claims (single fixed seed, as the acceptance gates):

* acceptance 4: eligible target-won genes T (truth-unique correct peptides at q <= 0.01 whose
  best key beats the gene's best paired decoy) is >= 110 (q_floor <= 0.00909) and inside the
  spec's expected 2,300-3,100 (truth-detected share 0.631 x 4,367 ~ 2,757);
* acceptance 2, 3: GIAP's ``pair_decoys`` + ``choose_tier`` give ``decoy_picked``, no warning,
  paired share f >= 0.98; acceptance 5: at least one decoy win;
* Q2 / decision 1: on the whole genome the search-space split gives well over 50 gene-level null
  PSMs (~0.14 x ~400 x 3), so the mode is ``search_space``; the protein share of the null
  peptides matches the recorded share (exact binomial, p >= 0.001), and the recorded share is
  the makedb search space's (+-0.04);
* lead decision 9: every protein-entry null is an undetected real product of its locus (no true
  peptide there) that maps to the genome, so genes without a true peptide win the picked
  competition (~ 0.5 x 7/8 x ~56 protein null peptides ~ 25 expected) and gate 7's
  false-discovery share among target-won genes is estimable and non-zero;
* rule 6 bookkeeping on every null PSM (once, target xor decoy; finite decoy scores).
"""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any

import pytest
from helpers import mvp04_data as M
from helpers import mvp04_decoys as D
from scipy import stats

pytestmark = pytest.mark.slow

T_RANGE = (2_300, 3_100)  # spec acceptance 4 "Expected"
T_MIN = 110  # spec acceptance 4 gate


@pytest.fixture(scope="module")
def k12(tmp_path_factory: pytest.TempPathFactory) -> D.View:
    sc = importlib.import_module("giapsim.scenarios")
    out = tmp_path_factory.mktemp("mvp04_k12_scored")
    tu = M.TUSET if M.TUSET.is_file() else None
    importlib.import_module("giapsim.real_genome").simulate_on_genome(
        sc.k12_evidence(), sc.K12_EVIDENCE_SEED, out, genome_fasta=M.K12_FASTA,
        bakta_json=M.K12_JSON, tu_set=tu,
    )  # fmt: skip
    return D.view(out)


@pytest.fixture(scope="module")
def genes() -> list[M.Gene]:
    return M.protein_genes(M.k12_json())


def test_provenance_and_bookkeeping(k12: D.View) -> None:
    assert k12.provenance["decoys"] == "scored_null"
    assert isinstance(k12.provenance["decoy_collisions"], int)
    bad = D.accounting_violations(k12)
    assert not bad, bad[:10]
    assert k12.decoy_rows
    for r in k12.decoy_rows:
        assert all(D.finite(r[c]) for c in ("q_value", "pep", "score")), r


def test_eligible_target_won_genes(k12: D.View, genes: list[M.Gene]) -> None:
    eligible = {g.locus for g in genes if g.eligible}
    t = len(D.target_won_genes(k12, eligible))
    assert t >= T_MIN, t
    assert T_RANGE[0] <= t <= T_RANGE[1], t


def test_giap_tier_pairing_and_decoy_wins(k12: D.View, genes: list[M.Gene]) -> None:
    tier = D.giap_tier(k12, [g.locus for g in genes])
    assert tier.tier == "decoy_picked", tier
    assert tier.warnings == (), tier
    assert tier.paired_share >= 0.98, tier.paired_share
    eligible = {g.locus for g in genes if g.eligible}
    wins = D.decoy_wins(k12, {g.locus: g.aa for g in genes}, eligible)
    assert sum(wins.values()) >= 1, len(wins)


def test_null_split_follows_the_search_space(
    k12: D.View, tmp_path_factory: pytest.TempPathFactory
) -> None:
    split = k12.provenance["null_split"]
    assert split["mode"] == "search_space", split
    share = float(split["search_space_protein_share"])
    n_p = sum(n.kind == "protein" for n in k12.nulls.values())
    n = len(k12.nulls)
    assert (split["n_protein"], split["n_sixframe"]) == (n_p, n - n_p), split
    assert n_p * len(k12.runs) >= D.FALLBACK_PSMS, (n_p, split)
    assert stats.binomtest(n_p, n, share).pvalue >= 0.001, (n_p, n, share)
    tmp = tmp_path_factory.mktemp("makedb_k12")
    fna = D.write_sample_fasta(k12.out, tmp / "sample.fna")
    space: dict[str, Any] = D.makedb_search_space(M.K12_JSON, fna, tmp)
    assert abs(share - space["protein_share"]) <= 0.04, (share, space["protein_share"])


def test_inputs_untouched(k12: D.View) -> None:
    assert M.sha256(M.K12_FASTA) == M.SHA_FASTA_FILE
    assert M.sha256(M.K12_JSON) == M.SHA_JSON_FILE
    assert Path(k12.out / "truth" / "sample_genome.fna").is_file()


def test_protein_nulls_are_undetected_products_that_map(k12: D.View) -> None:
    recorded = D.annotation(M.K12_JSON)
    sample = D.read_seq(k12.out / "truth" / "sample_genome.fna")
    proteins = D.locus_proteins(k12.out, recorded, sample)
    bad = D.protein_null_violations(k12, recorded, proteins, sample, M.k12_text()[1])
    assert not any(bad.values()), {k: x[:5] for k, x in bad.items() if x}


def test_gate7_false_discovery_share_is_estimable(k12: D.View, genes: list[M.Gene]) -> None:
    eligible = {g.locus for g in genes if g.eligible}
    sample = D.read_seq(k12.out / "truth" / "sample_genome.fna")
    proteins = D.locus_proteins(k12.out, D.annotation(M.K12_JSON), sample)
    comp = D.gene_competition(k12, eligible, proteins)
    false_won = sum(c.target_won and not c.has_true for c in comp)
    won = sum(c.target_won for c in comp)
    assert false_won >= 1, (false_won, won)
    assert 0 < false_won / won < 0.5, (false_won, won)
