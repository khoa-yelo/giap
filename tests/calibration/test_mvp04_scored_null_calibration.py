"""Calibration of the MVP-04 scored null decoys (docs/specs/MVP-04-ecoli-evidence-run.md pipeline
step 6, rule 6, Q2, lead decision 1; M7-12 decisions 21 and 24; GIAP_SPEC.md §7.3.2, §8.4 suite
"gene FDR"). The interface is the one in tests/unit/test_mvp04_scored_null.py, the oracle
tests/helpers/mvp04_decoys.py.

Each replicate is ``k12_evidence()`` with ``decoys = "scored_null"`` on the test-time 200-kb
excerpt (spec rule 3; gap-rule TUs, 20 novel ORFs per class), seed
``seeds.calibration_base + 8400 + r`` (fast) or ``+ 8500 + r`` (full). About 18 null peptides
x 3 runs = 54 null PSMs per replicate. Claims, each pooled over the replicates:

* **Coin (rule 6).** Each null PSM is a decoy with probability 1/2: exact two-sided binomial test
  of D among D + T, p >= 0.001; per entry kind (protein / six-frame), p >= 0.001 / 2
  (Bonferroni).
* **Coin, per row (rule 6, decision 24: "an independent 1/2 coin per row").** Per null peptide
  (reversal colliding with no target) the number of its 3 runs that went to the decoy side is
  Binomial(3, 1/2): chi-square goodness of fit on the 4 cells, p >= 0.001. One coin per peptide
  (all runs on one side) puts all the mass on 0 and 3 and fails.
* **Coin, dispersion.** sum_r (D_r - N_r / 2)^2 / (N_r / 4) is chi-square with R degrees of
  freedom; both tails >= 0.001 (a deterministic one-to-one split fails the lower tail).
* **Coin independent of the score.** The decoy and null-target scores are one distribution
  (two-sample KS, p >= ``ks_min_p``), and the decoy scores are N(mu0, sd0) of ``score_mixture``
  (one-sample KS, p >= ``ks_min_p``).
* **50/50 fallback (Q2, decision 1).** The excerpt is below the ~50 gene-level null PSM
  threshold, so the protein share of the null peptides is 1/2: exact binomial test, p >= 0.001.

* **Protein nulls drawn uniformly over loci (lead decision 9).** The locus of a protein-entry
  null is uniform over the expressible recorded loci with no true peptide and >= 1 undetected
  product, not over products: the mean number of undetected products of the drawn loci matches
  the pool mean (z test, |z| <= 3.29, two-sided 0.001). A draw over products is size-biased
  (mean of k^2 / mean of k) and fails.
* **The coin gives the picked gene FDR its null (rule 6, decision 9).** Among eligible genes with
  no true peptide and null evidence at q <= 0.01, target and decoy each win the picked
  competition with probability 1/2: exact binomial test, p >= 0.001.

Fast (``calibration``): 20 replicates, ~0.5 s each. Full (``calibration_full``):
``replicates_fast.gene_fdr`` replicates.
"""

from __future__ import annotations

import importlib
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from helpers import mvp04_data as M
from helpers import mvp04_decoys as D
from scipy import stats

from giapsim.expression import IslandSpec

pytestmark = pytest.mark.calsuite("gene_fdr")

ALPHA = 0.001
N_FAST = 20
ISLANDS = ((2, 3.0), (18, -3.0), (34, 3.0), (50, -3.0), (66, 3.0), (82, -3.0))


def _cfg() -> Any:
    base = importlib.import_module("giapsim.scenarios").k12_evidence()
    expr = replace(
        base.expression, islands=tuple(IslandSpec(M.CONTIG, f, 8, e) for f, e in ISLANDS)
    )
    return replace(base, expression=expr, n_novel_per_class=20,
                   proteomics=replace(base.proteomics, decoys="scored_null"))  # fmt: skip


def _rep(ex: M.Excerpt, out: Path, seed: int) -> D.View:
    importlib.import_module("giapsim.real_genome").simulate_on_genome(
        _cfg(), seed, out, genome_fasta=ex.fasta, bakta_json=ex.bakta_json, tu_set=None,
        check_hashes=False,
    )  # fmt: skip
    return D.view(out)


def _check(reps: list[D.View], ks_min_p: float, ex: M.Excerpt) -> None:
    _, _, mu0, sd0 = (float(x) for x in _cfg().proteomics.score_mixture)
    psms = [(i, p) for i, v in enumerate(reps) for p in D.null_psms(v)]
    ds = np.array([float(r["score"]) for _, (_, _, r, d) in psms if d])
    ts = np.array([float(r["score"]) for _, (_, _, r, d) in psms if not d])
    assert len(ds) >= 300 and len(ts) >= 300, (len(ds), len(ts))

    # coin, pooled and per entry kind
    n = len(ds) + len(ts)
    assert stats.binomtest(len(ds), n, 0.5).pvalue >= ALPHA, (len(ds), n)
    for kind in ("protein", "sixframe"):
        dk = sum(1 for _, (nl, _, _, d) in psms if d and nl.kind == kind)
        nk = sum(1 for _, (nl, _, _, _) in psms if nl.kind == kind)
        assert nk >= 50, (kind, nk)
        assert stats.binomtest(dk, nk, 0.5).pvalue >= ALPHA / 2, (kind, dk, nk)

    # coin per row: decoy runs per null peptide ~ Binomial(3, 1/2)
    cells = np.zeros(4)
    for v in reps:
        targets = {D.il(r["peptide"]) for r in v.target_rows} | {D.il(p) for p in v.correct}
        for nl in v.nulls.values():
            if D.il(nl.peptide[::-1]) not in targets:
                cells[len(nl.decoys)] += 1
    expect = cells.sum() * np.array([1, 3, 3, 1]) / 8
    assert expect.min() >= 5, cells
    assert stats.chisquare(cells, expect).pvalue >= ALPHA, cells

    # coin, dispersion across replicates
    chi = 0.0
    for i in range(len(reps)):
        d_r = sum(1 for j, (_, _, _, d) in psms if j == i and d)
        n_r = sum(1 for j, _ in psms if j == i)
        chi += (d_r - n_r / 2) ** 2 / (n_r / 4)
    df = len(reps)
    assert stats.chi2.cdf(chi, df) >= ALPHA, ("split looks deterministic", chi, df)
    assert stats.chi2.sf(chi, df) >= ALPHA, ("split overdispersed", chi, df)

    # the coin does not look at the score; decoys keep the null score
    assert stats.ks_2samp(ds, ts).pvalue >= ks_min_p
    assert stats.kstest(ds, stats.norm(mu0, sd0).cdf).pvalue >= ks_min_p

    # Q2 fallback: 50/50 between protein and six-frame null peptides
    n_p = sum(sum(nl.kind == "protein" for nl in v.nulls.values()) for v in reps)
    n_all = sum(len(v.nulls) for v in reps)
    assert all(v.provenance["null_split"]["mode"] == "fifty_fifty" for v in reps)
    assert stats.binomtest(n_p, n_all, 0.5).pvalue >= ALPHA, (n_p, n_all)

    # decision 9: loci uniform over the pool, not size-biased by their product counts
    recorded = D.annotation(ex.bakta_json)
    num = var = 0.0
    k_draws = 0
    for v in reps:
        sample = D.read_seq(v.out / "truth" / "sample_genome.fna")
        pool = D.candidate_loci(v, D.locus_proteins(v.out, recorded, sample))
        ks = np.array(list(pool.values()), dtype=float)
        for nl in v.nulls.values():
            if nl.kind != "protein":
                continue
            (tag,) = nl.entries()
            assert tag in pool, (tag, nl.peptide)
            num += pool[tag] - ks.mean()
            var += ks.var()
            k_draws += 1
    assert k_draws >= 100, k_draws
    z = num / np.sqrt(var)
    assert abs(z) <= 3.29, ("protein-null loci not uniform over the pool", z, k_draws)

    # the coin is the picked null: genes with only null evidence, target vs decoy wins
    eligible = {g.locus for g in ex.genes if g.eligible}
    comp = [c for v in reps for c in D.gene_competition(v, eligible) if not c.has_true]
    won = sum(c.target_won for c in comp)
    assert len(comp) >= 100, len(comp)
    assert stats.binomtest(won, len(comp), 0.5).pvalue >= ALPHA, (won, len(comp))


@pytest.fixture(scope="module")
def ex(tmp_path_factory: pytest.TempPathFactory) -> M.Excerpt:
    return M.make_excerpt(tmp_path_factory.mktemp("k12_excerpt_cal"))


@pytest.fixture(scope="module")
def base(tolerances: dict[str, Any]) -> int:
    return int(tolerances["seeds"]["calibration_base"])


@pytest.mark.calibration
def test_scored_null_coin_fast(
    ex: M.Excerpt, tmp_path_factory: pytest.TempPathFactory, base: int,
    tolerances: dict[str, Any],
) -> None:  # fmt: skip
    reps = [_rep(ex, tmp_path_factory.mktemp(f"mvp04c{r}"), base + 8400 + r)
            for r in range(N_FAST)]  # fmt: skip
    _check(reps, float(tolerances["ks_min_p"]), ex)


@pytest.mark.calibration_full
def test_scored_null_coin_full(
    ex: M.Excerpt, tmp_path_factory: pytest.TempPathFactory, base: int,
    tolerances: dict[str, Any],
) -> None:  # fmt: skip
    n = int(tolerances["replicates_fast"]["gene_fdr"])
    reps = [_rep(ex, tmp_path_factory.mktemp(f"mvp04cF{r}"), base + 8500 + r)
            for r in range(n)]  # fmt: skip
    _check(reps, float(tolerances["ks_min_p"]), ex)
