"""Attack nulls for the MVP-04 scored null decoys under the picked gene FDR (GIAP_SPEC.md §7.3.2,
§8.1 row "Gene detected", §8.4 "FDR", §8.5; docs/specs/MVP-04-ecoli-evidence-run.md pipeline
step 6, lead decision 9, acceptance 5-7; statistician review of PR #247).

The MVP-04 calibration file checks the ingredients of the null (the 1/2 coin per row, the
uniform locus draw, the fair picked label). These attacks check the claim the evidence run rests
on: that GIAP's picked gene FDR, run on ``scored_null`` output, controls the gene FDP, and that
gate 7 (realised gene FDP <= 0.02) could fail at all.

The gene competition mirrors ``giap.genetable.detection`` with an oracle that does not run the
mapping: target evidence of an eligible gene = target rows at q <= ``input.q_value_max`` of
correct peptides unique to it in the truth, plus protein-entry null target rows whose peptide
(I = L) lies in the gene's sample protein; a null peptide found in several loci gives evidence to
all of them only when their proteins are identical (GIAP accepts ``multi_identical``, not
``shared``). Decoy evidence = its ``rev_<gene>`` rows at q <= ``input.q_value_max``. Six-frame
nulls on a ``GIAPORF_nterm_extension_*`` entry also reach a gene: their target rows place in
the gene's ``nterm_extension`` locus class, which assigns (``contract.ASSIGNING_CLASSES``), and
their ``rev_GIAPORF_*`` decoys pair with that gene only when the entry was voted to it, i.e.
when the peptide has a target row at q <= ``input.q_value_max`` (``protein_gene_map`` step 3,
M7-12 decision 28). The oracle models both. Rows are ranked by the M4-24 shared key (q asc,
PEP asc, score desc; S = {q, pep, score} for scored nulls) and the best ranks go to
``giap.stats.tda.picked_fdr`` (ties to the decoy).

- **Dense null, FDR** (fast R = 50, full R = 300; ``calibrate``, ADR-0015): the 200-kb excerpt at
  pi_wrong = 0.05 with the fixed 50/50 split (lead decision 10), so about 47 protein-entry nulls
  meet about 100 truth-detected genes and the (D + 1) / T cut binds at q = 0.05 and 0.10. Mean
  FDP <= q + 2 SE at q in {0.01, 0.05, 0.10}.
- **Teeth**: the same replicates with the decoys ignored (a null-blind engine, every gene with
  target evidence detected) must be a confirmed failure at q = 0.01 and 0.05.
- **Gate 7 on the evidence run** (full, slow; lead decision 10): on the whole K-12 genome at
  ``K12_EVIDENCE_SEED`` the fixed 50/50 split (decision 16 pool) gives 80-89 protein-entry nulls
  and 36-54 decoy wins against T ~ 2,700 (Slurm 1872374, 10 seeds). A null-blind engine then has
  a gene FDP of 0.023-0.029, so gate 7 (<= 0.02) can fail, and the alpha = 0.01 cut binds
  (some target-won gene is left out; 10 of 10 seeds). Under the withdrawn search-space split
  both failed (FDP 0.0076; n_detected = T), which is why these were strict xfails before.
  "True" genes are the loci of truth rows with ``is_wrong = false`` (lead decision 11).
"""

from __future__ import annotations

import importlib
from dataclasses import replace
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest
from helpers import mvp04_data as M
from helpers import mvp04_decoys as D

from giap.stats.tda import picked_fdr
from giapsim.expression import IslandSpec

pytestmark = pytest.mark.calsuite("gene_fdr")

SUITE = "gene_fdr"
LEVELS = (0.01, 0.05, 0.10)
PI_DENSE = 0.05
GATE7 = 0.02  # spec acceptance 7 gate
ALPHA_ID = 0.01  # stats.alpha.identification
ISLANDS = ((2, 3.0), (18, -3.0), (34, 3.0), (50, -3.0), (66, 3.0), (82, -3.0))
STREAM = "gene_fdr.attack_mvp04_scored_null"

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]


# ---------------------------------------------------------------- the oracle gene competition


def _gene_keys(v: D.View, genes_: list[M.Gene], proteins: dict[str, list[str]]
               ) -> tuple[list[str], np.ndarray, np.ndarray]:  # fmt: skip
    """(genes, best target rank, best decoy rank) from ``D.gene_evidence`` (the oracle described
    in the module doc); rank higher = better, NaN = no evidence."""
    ev = D.gene_evidence(v, genes_, proteins)
    genes = sorted({e.gene for e in ev})
    if not ev:
        return genes, np.zeros(0), np.zeros(0)
    keys = np.array([e.key for e in ev], dtype=float)
    keys = np.where(np.isnan(keys), np.inf, keys)
    _uniq, inv = np.unique(keys, axis=0, return_inverse=True)
    rank = len(_uniq) - np.ravel(inv)  # dense, higher = better (as detection._key_rank)
    idx = {g: i for i, g in enumerate(genes)}
    best = np.full((2, len(genes)), np.nan)
    for e, rk in zip(ev, rank.tolist(), strict=True):
        j, k = int(e.is_decoy), idx[e.gene]
        best[j, k] = rk if np.isnan(best[j, k]) else max(best[j, k], rk)
    return genes, best[0], best[1]


def _fdp(genes: list[str], t: np.ndarray, d: np.ndarray, true: set[str], q: float) -> float:
    res = picked_fdr(t, d)
    det = [g for g, x in zip(genes, res.q.tolist(), strict=True) if np.isfinite(x) and x <= q]
    return sum(g not in true for g in det) / len(det) if det else 0.0


# ---------------------------------------------------------------- dense null on the excerpt


def _dense_cfg() -> Any:
    base = importlib.import_module("giapsim.scenarios").k12_evidence()
    expr = replace(
        base.expression, islands=tuple(IslandSpec(M.CONTIG, f, 8, e) for f, e in ISLANDS)
    )
    prot = replace(base.proteomics, decoys="scored_null", pi_wrong=PI_DENSE)
    return replace(base, expression=expr, n_novel_per_class=20, proteomics=prot)


@pytest.fixture(scope="module")
def ex(tmp_path_factory: pytest.TempPathFactory) -> M.Excerpt:
    return M.make_excerpt(tmp_path_factory.mktemp("k12_excerpt_attack"))


class Dense:
    """One simulated excerpt replicate per harness replicate (shared by every level)."""

    def __init__(self, ex: M.Excerpt, tmp: pytest.TempPathFactory) -> None:
        self.ex, self.tmp = ex, tmp
        self.memo: dict[int, tuple[list[str], np.ndarray, np.ndarray, set[str], str]] = {}
        self.eligible = {g.locus for g in ex.genes if g.eligible}

    def rep(self, rng: np.random.Generator) -> tuple[list[str], np.ndarray, np.ndarray, set[str]]:
        seed = int(rng.bit_generator.seed_seq.entropy) % 2**31  # type: ignore[attr-defined]
        if seed not in self.memo:
            out = self.tmp.mktemp(f"mvp04att{seed}")
            importlib.import_module("giapsim.real_genome").simulate_on_genome(
                _dense_cfg(), seed, out, genome_fasta=self.ex.fasta,
                bakta_json=self.ex.bakta_json, tu_set=None, check_hashes=False,
            )  # fmt: skip
            v = D.view(out)
            prots = D.locus_proteins(out, D.annotation(self.ex.bakta_json),
                                     D.read_seq(out / "truth" / "sample_genome.fna"))  # fmt: skip
            genes, t, d = _gene_keys(v, self.ex.genes, prots)
            mode = str(v.provenance["null_split"]["mode"])
            self.memo[seed] = (genes, t, d, D.true_loci(v), mode)
        genes, t, d, true, mode = self.memo[seed]
        assert mode == "fifty_fifty", mode  # the dense regime this attack is built on
        return genes, t, d, true


@pytest.fixture(scope="module")
def dense(ex: M.Excerpt, tmp_path_factory: pytest.TempPathFactory) -> Dense:
    return Dense(ex, tmp_path_factory)


@pytest.mark.parametrize("full", MODES)
def test_dense_scored_null_gene_fdr(dense: Dense, full: bool) -> None:
    R = C.replicates(SUITE, full)
    stream = STREAM + (".full" if full else "")
    statuses = {}
    for q in LEVELS:

        def fn(rng: np.random.Generator, q: float = q) -> float:
            genes, t, d, true = dense.rep(rng)
            return _fdp(genes, t, d, true, q)

        statuses[q] = C.calibrate(SUITE, "fdr", q, fn, R, stream=stream)["status"]
    assert "fail" not in statuses.values(), statuses
    # non-vacuous: the null is dense enough that the cut binds (false genes exist to be cut)
    reps = list(dense.memo.values())[:R]
    n_null = np.mean([sum(g not in tr for g in gs) for gs, _, _, tr, _ in reps])
    n_true = np.mean([sum(g in tr for g in gs) for gs, _, _, tr, _ in reps])
    assert n_null >= 0.15 * n_true, (n_null, n_true)


@pytest.mark.calibration
def test_teeth_null_blind_engine_fails(dense: Dense, teeth_log: Path) -> None:
    R = C.replicates(SUITE, False)
    statuses = {}
    for q in (0.01, 0.05):

        def fn(rng: np.random.Generator, q: float = q) -> float:
            genes, t, _d, true = dense.rep(rng)
            return _fdp(genes, t, np.full(t.shape, np.nan), true, q)

        statuses[q] = C.calibrate(SUITE, "fdr", q, fn, R, stream=STREAM, log_dir=teeth_log)[
            "status"
        ]
    assert set(statuses.values()) == {"fail"}, statuses


# ---------------------------------------------------------------- gate 7 on the evidence run


@pytest.fixture(scope="module")
def k12(tmp_path_factory: pytest.TempPathFactory) -> tuple[list[str], np.ndarray, np.ndarray,
                                                            set[str]]:  # fmt: skip
    sc = importlib.import_module("giapsim.scenarios")
    out = tmp_path_factory.mktemp("mvp04_k12_attack")
    tu = M.TUSET if M.TUSET.is_file() else None
    importlib.import_module("giapsim.real_genome").simulate_on_genome(
        sc.k12_evidence(), sc.K12_EVIDENCE_SEED, out, genome_fasta=M.K12_FASTA,
        bakta_json=M.K12_JSON, tu_set=tu,
    )  # fmt: skip
    v = D.view(out)
    prots = D.locus_proteins(out, D.annotation(M.K12_JSON),
                             D.read_seq(out / "truth" / "sample_genome.fna"))  # fmt: skip
    genes, t, d = _gene_keys(v, M.protein_genes(M.k12_json()), prots)
    return genes, t, d, D.true_loci(v)


@pytest.mark.slow
@pytest.mark.calibration_full
def test_gate7_can_fail_on_the_evidence_run(k12: Any) -> None:
    """A null-blind engine (decoys ignored) must break gate 7, else gate 7 cannot tell a correct
    picked FDR from none."""
    genes, t, _d, true = k12
    fdp = _fdp(genes, t, np.full(t.shape, np.nan), true, ALPHA_ID)
    assert fdp > GATE7, fdp


@pytest.mark.slow
@pytest.mark.calibration_full
def test_alpha_cut_binds_on_the_evidence_run(k12: Any) -> None:
    """At alpha = 0.01 some target-won gene is left out; otherwise every gene q is the floor of
    the list ((D + 1) / T at the loosest threshold) and the run is saturated in the sense of the
    spec's 1865149 contrast."""
    _genes, t, d, _true = k12
    res = picked_fdr(t, d)
    q = res.q[np.isfinite(res.q)]
    assert int((q <= ALPHA_ID).sum()) < q.size, (q.size, float(q.max()))
