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

**Two-block confirmation (ADR-0015, lead decision on Slurm 1871474).** Every p-value gate above
is a named check. All checks run on block 1 (replicates 0..R-1 of the seed offset); a check
that fails there is re-run on the independent block 2 (replicates R..2R-1, the next R seeds of
the same offset) and the test fails only if both blocks fail; a block-1 failure that block 2
does not confirm is ``unconfirmed`` (a warning, printed with both p values). Structural
preconditions (enough PSMs, the fallback mode, every drawn locus in the pool) are plain
asserts. Teeth: one coin per peptide (every run of a null peptide on one side) and a draw of
the locus weighted by its product count are each a *confirmed* failure of their check on the
real replicates of both blocks.

Fast (``calibration``): R = 20 replicates (~0.5 s each; block 2 is also built for the teeth
tests). Full (``calibration_full``): R = ``replicates_fast.gene_fdr``.
"""

from __future__ import annotations

import importlib
import warnings
from collections.abc import Callable
from dataclasses import dataclass, replace
from functools import cache
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


@dataclass(frozen=True)
class Rep:
    """What the checks read from one replicate."""

    psms: tuple[tuple[str, bool, float], ...]  # (entry kind, is_decoy, score) per null PSM
    cells: tuple[int, int, int, int]  # null peptides (no collision) by number of decoy runs
    n_protein: int
    n_nulls: int
    mode: str
    pool: dict[str, int]  # decision 9 pool: locus -> undetected products
    drawn: tuple[str, ...]  # attributed loci of the protein-entry nulls
    null_wins: tuple[bool, ...]  # target won, per eligible gene with only null evidence


def _summary(v: D.View, ex: M.Excerpt) -> Rep:
    recorded = _recorded(ex)
    targets = {D.il(r["peptide"]) for r in v.target_rows} | {D.il(p) for p in v.correct}
    cells = [0, 0, 0, 0]
    for nl in v.nulls.values():
        if D.il(nl.peptide[::-1]) not in targets:
            cells[len(nl.decoys)] += 1
    sample = D.read_seq(v.out / "truth" / "sample_genome.fna")
    prots = D.locus_proteins(v.out, recorded, sample)
    pool = D.candidate_loci(v, prots)
    drawn = []
    for nl in v.nulls.values():
        if nl.kind == "protein":
            (tag,) = nl.entries()
            assert tag in pool, (tag, nl.peptide)
            drawn.append(tag)
    eligible = {g.locus for g in ex.genes if g.eligible}
    wins = tuple(c.target_won for c in D.gene_competition(v, eligible, prots) if not c.has_true)
    psms = tuple((nl.kind, d, float(r["score"])) for nl, _, r, d in D.null_psms(v))
    return Rep(psms, (cells[0], cells[1], cells[2], cells[3]),
               sum(nl.kind == "protein" for nl in v.nulls.values()), len(v.nulls),
               str(v.provenance["null_split"]["mode"]), pool, tuple(drawn), wins)  # fmt: skip


@cache
def _annotation(path: Path) -> dict[str, dict[str, Any]]:
    return D.annotation(path)


def _recorded(ex: M.Excerpt) -> dict[str, dict[str, Any]]:
    return _annotation(ex.bakta_json)


def _preconditions(reps: list[Rep]) -> None:
    """Structural asserts (not p-value gates): sizes the tests need and the fallback mode."""
    ds = [s for r in reps for _, d, s in r.psms if d]
    ts = [s for r in reps for _, d, s in r.psms if not d]
    assert len(ds) >= 300 and len(ts) >= 300, (len(ds), len(ts))
    for kind in ("protein", "sixframe"):
        nk = sum(1 for r in reps for k, _, _ in r.psms if k == kind)
        assert nk >= 50, (kind, nk)
    cells = np.sum([r.cells for r in reps], axis=0)
    assert cells.sum() * 1 / 8 >= 5, cells
    assert sum(len(r.drawn) for r in reps) >= 100
    assert sum(len(r.null_wins) for r in reps) >= 100
    assert all(r.mode == "fifty_fifty" for r in reps), {r.mode for r in reps}


Check = Callable[[list[Rep]], float]  # the p value of one gate on one block


def _p_coin(kind: str | None) -> Check:
    def p(reps: list[Rep]) -> float:
        sel = [d for r in reps for k, d, _ in r.psms if kind is None or k == kind]
        return float(stats.binomtest(sum(sel), len(sel), 0.5).pvalue)

    return p


def _p_per_row(reps: list[Rep]) -> float:
    cells = np.sum([r.cells for r in reps], axis=0).astype(float)
    return float(stats.chisquare(cells, cells.sum() * np.array([1, 3, 3, 1]) / 8).pvalue)


def _chi_dispersion(reps: list[Rep]) -> float:
    chi = 0.0
    for r in reps:
        n_r = len(r.psms)
        chi += (sum(d for _, d, _ in r.psms) - n_r / 2) ** 2 / (n_r / 4)
    return chi


def _p_disp_low(reps: list[Rep]) -> float:
    return float(stats.chi2.cdf(_chi_dispersion(reps), len(reps)))


def _p_disp_high(reps: list[Rep]) -> float:
    return float(stats.chi2.sf(_chi_dispersion(reps), len(reps)))


def _p_ks2(reps: list[Rep]) -> float:
    ds = [s for r in reps for _, d, s in r.psms if d]
    ts = [s for r in reps for _, d, s in r.psms if not d]
    return float(stats.ks_2samp(ds, ts).pvalue)


def _p_ks_null(reps: list[Rep]) -> float:
    _, _, mu0, sd0 = (float(x) for x in _cfg().proteomics.score_mixture)
    ds = [s for r in reps for _, d, s in r.psms if d]
    return float(stats.kstest(ds, stats.norm(mu0, sd0).cdf).pvalue)


def _p_split(reps: list[Rep]) -> float:
    n_p, n = sum(r.n_protein for r in reps), sum(r.n_nulls for r in reps)
    return float(stats.binomtest(n_p, n, 0.5).pvalue)


def _p_locus_uniform(reps: list[Rep]) -> float:
    num = var = 0.0
    for r in reps:
        ks = np.array(list(r.pool.values()), dtype=float)
        for tag in r.drawn:
            num += r.pool[tag] - ks.mean()
            var += ks.var()
    return float(2 * stats.norm.sf(abs(num / np.sqrt(var))))


def _p_picked_null(reps: list[Rep]) -> float:
    w = [x for r in reps for x in r.null_wins]
    return float(stats.binomtest(sum(w), len(w), 0.5).pvalue)


def checks(ks_min_p: float) -> dict[str, tuple[Check, float]]:
    """Every p-value gate: name -> (p value on a block, the gate it must reach)."""
    return {
        "coin_pooled": (_p_coin(None), ALPHA),
        "coin_protein": (_p_coin("protein"), ALPHA / 2),
        "coin_sixframe": (_p_coin("sixframe"), ALPHA / 2),
        "coin_per_row": (_p_per_row, ALPHA),
        "dispersion_low": (_p_disp_low, ALPHA),
        "dispersion_high": (_p_disp_high, ALPHA),
        "ks_decoy_vs_target": (_p_ks2, ks_min_p),
        "ks_decoy_null": (_p_ks_null, ks_min_p),
        "split_fifty_fifty": (_p_split, ALPHA),
        "locus_uniform": (_p_locus_uniform, ALPHA),
        "picked_null": (_p_picked_null, ALPHA),
    }


def confirm(
    gates: dict[str, tuple[Check, float]], block1: list[Rep], block2: Callable[[], list[Rep]]
) -> dict[str, tuple[str, float, float | None]]:
    """ADR-0015: name -> (status pass | unconfirmed | fail, block-1 p, block-2 p or None)."""
    out: dict[str, tuple[str, float, float | None]] = {}
    for name, (p, gate) in gates.items():
        p1 = p(block1)
        if p1 >= gate:
            out[name] = ("pass", p1, None)
            continue
        p2 = p(block2())
        out[name] = ("fail" if p2 < gate else "unconfirmed", p1, p2)
    return out


def _verdict(res: dict[str, tuple[str, float, float | None]]) -> None:
    for name, (status, p1, p2) in res.items():
        print(f"mvp04 calibration {name}: {status} p1={p1:.4g} p2={p2}")
        if status == "unconfirmed":
            warnings.warn(f"{name}: block 1 p={p1:.4g} unconfirmed by block 2 p={p2}",
                          stacklevel=2)  # fmt: skip
    failed = {k: v for k, v in res.items() if v[0] == "fail"}
    assert not failed, failed


@pytest.fixture(scope="module")
def ex(tmp_path_factory: pytest.TempPathFactory) -> M.Excerpt:
    return M.make_excerpt(tmp_path_factory.mktemp("k12_excerpt_cal"))


@pytest.fixture(scope="module")
def base(tolerances: dict[str, Any]) -> int:
    return int(tolerances["seeds"]["calibration_base"])


class Blocks:
    """Block 1 = replicates 0..R-1, block 2 = R..2R-1 of seed ``base + offset + r``; each built
    once, block 2 only when asked for."""

    def __init__(self, ex: M.Excerpt, tmp: pytest.TempPathFactory, first: int, R: int) -> None:
        self.ex, self.tmp, self.first, self.R = ex, tmp, first, R
        self._b: dict[int, list[Rep]] = {}

    def block(self, k: int) -> list[Rep]:
        if k not in self._b:
            seeds = range(self.first + (k - 1) * self.R, self.first + k * self.R)
            self._b[k] = [_summary(_rep(self.ex, self.tmp.mktemp(f"mvp04c{s}"), s), self.ex)
                          for s in seeds]  # fmt: skip
        return self._b[k]


@pytest.fixture(scope="module")
def fast(ex: M.Excerpt, tmp_path_factory: pytest.TempPathFactory, base: int) -> Blocks:
    return Blocks(ex, tmp_path_factory, base + 8400, N_FAST)


@pytest.mark.calibration
def test_scored_null_coin_fast(fast: Blocks, tolerances: dict[str, Any]) -> None:
    b1 = fast.block(1)
    _preconditions(b1)
    _verdict(confirm(checks(float(tolerances["ks_min_p"])), b1, lambda: fast.block(2)))


def _one_coin_per_peptide(reps: list[Rep], seed: int) -> list[Rep]:
    """Teeth: every null peptide's runs all on one side (a coin per peptide, not per row)."""
    rng = np.random.default_rng(seed)
    out = []
    for r in reps:
        n = sum(r.cells)
        k = int(rng.binomial(n, 0.5))
        out.append(replace(r, cells=(n - k, 0, 0, k)))
    return out


def _weighted_draw(reps: list[Rep], seed: int) -> list[Rep]:
    """Teeth: the loci drawn with probability proportional to their product count."""
    rng = np.random.default_rng(seed)
    out = []
    for r in reps:
        tags = sorted(r.pool)
        w = np.array([r.pool[t] for t in tags], dtype=float)
        pick = rng.choice(len(tags), size=len(r.drawn), replace=False, p=w / w.sum())
        out.append(replace(r, drawn=tuple(tags[i] for i in pick.tolist())))
    return out


@pytest.mark.calibration
@pytest.mark.parametrize(
    ("name", "spoil"),
    [("coin_per_row", _one_coin_per_peptide), ("locus_uniform", _weighted_draw)],
)
def test_teeth_are_confirmed_failures(
    fast: Blocks, tolerances: dict[str, Any], name: str,
    spoil: Callable[[list[Rep], int], list[Rep]],
) -> None:  # fmt: skip
    gate = {name: checks(float(tolerances["ks_min_p"]))[name]}
    res = confirm(gate, spoil(fast.block(1), 1), lambda: spoil(fast.block(2), 2))
    assert res[name][0] == "fail", res


@pytest.mark.calibration_full
def test_scored_null_coin_full(
    ex: M.Excerpt, tmp_path_factory: pytest.TempPathFactory, base: int,
    tolerances: dict[str, Any],
) -> None:  # fmt: skip
    n = int(tolerances["replicates_fast"]["gene_fdr"])
    full = Blocks(ex, tmp_path_factory, base + 8500, n)
    b1 = full.block(1)
    _preconditions(b1)
    _verdict(confirm(checks(float(tolerances["ks_min_p"])), b1, lambda: full.block(2)))
