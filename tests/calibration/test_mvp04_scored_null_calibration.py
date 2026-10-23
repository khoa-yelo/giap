"""Calibration of the MVP-04 scored null decoys (docs/specs/MVP-04-ecoli-evidence-run.md pipeline
step 6, rule 6, Q2, lead decision 1; M7-12 decisions 21 and 24; GIAP_SPEC.md §7.3.2, §8.4 suite
"gene FDR"). The interface is the one in tests/unit/test_mvp04_scored_null.py, the oracle
tests/helpers/mvp04_decoys.py.

Each replicate is ``k12_evidence()`` with ``decoys = "scored_null"`` on the test-time 200-kb
excerpt (spec rule 3; gap-rule TUs, 20 novel ORFs per class). Replicate i of stream ``STREAM``
(``STREAM + ".full"`` in the full tier) has seed ``_criteria.suite_seed(stream, i) mod 2**31``,
the harness seed formula (statistician review of PR #247, item 7; lead decision 12). About 18
null peptides x 3 runs = 54 null PSMs per replicate. Claims, each pooled over the replicates:

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
* **Fixed 50/50 split (lead decision 10).** Per replicate ``mode = fifty_fifty`` and
  ``n_protein`` is half the null peptides (either rounding): an exact assert, not a binomial
  test, since the count is fixed by construction (review item 5).

* **Protein nulls drawn uniformly over loci (lead decision 9).** The locus of a protein-entry
  null is uniform over the expressible recorded loci with no true peptide and >= 1 undetected
  product, not over products: the mean number of undetected products of the drawn loci matches
  the pool mean (z test, |z| <= 3.29, two-sided 0.001). A draw over products is size-biased
  (mean of k^2 / mean of k) and fails.
* **The coin gives the picked gene FDR its null (rule 6, decision 9).** Among eligible genes with
  no true peptide (no truth row with ``is_wrong = false``, decision 11) and null evidence at
  q <= 0.01, target and decoy each win the picked competition with probability 1/2: exact
  binomial test, p >= 0.001. The oracle (``mvp04_decoys.gene_evidence``) includes the
  ``nterm_extension`` six-frame channel and ``multi_identical`` acceptance (decision 11).
  Power against P = 0.6 at the realised block size is printed and logged; the full tier
  (R = 300) must reach >= 0.99 (decision 12, review item 6).

**Two-block confirmation (ADR-0015) through the harness's seeds and log.** Every p-value gate
above is a named check. Block 1 is replicates 0..R-1, block 2 replicates R..2R-1 (the harness
blocks). A check that fails on block 1 is re-run on block 2 and the test fails only if both
blocks fail; a block-1 failure that block 2 does not confirm is ``unconfirmed`` (a warning).
Each check writes one record per block it ran to ``<GIAP_CALIBRATION_LOG_DIR>/gene_fdr.jsonl``
(``criterion = "pooled_p"``, ``check``, ``observed`` = p, ``bound`` = its gate, ``block``,
``status``, ``seed`` of the block's first replicate, the estimate and its 95 % CI), so
``unconfirmed`` reaches the nightly summary. The records carry no ``values``, so
``merge_shards`` leaves them alone; pooled checks need whole blocks, so they skip under
``--shard`` / ``--merge``. Each estimate is printed with its 95 % CI next to its p: decoy share
(Clopper-Pearson), share of peptides with all runs on one side (expected 1/4), dispersion
factor chi / R (chi-square CI), decoy-minus-target mean score, decoy mean minus mu0, mean-k
difference of the drawn loci, target-win share.

Structural preconditions (enough PSMs, the fixed split, every drawn locus in the pool) are
plain asserts. Teeth, each a *confirmed* failure of its check on the real replicates of both
blocks: one coin per peptide (``coin_per_row``); a draw of the locus weighted by its product
count (``locus_uniform``); a score-dependent coin, P(decoy) = 0.3 above mu0 and 0.7 below
(``ks_decoy_vs_target``, ``ks_decoy_null``, ``picked_null``; it keeps the pooled decoy share at
about 1/2); a p = 0.55 coin (``coin_pooled``).

Fast (``calibration``): R = ``replicates_fast.gene_fdr`` (50). Full (``calibration_full``):
R = ``replicates_full.gene_fdr`` (300; review item 2, lead decision 12).
"""

from __future__ import annotations

import importlib
import json
import time
import warnings
from collections.abc import Callable
from dataclasses import dataclass, replace
from functools import cache
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest
from helpers import mvp04_data as M
from helpers import mvp04_decoys as D
from scipy import stats

from giapsim.expression import IslandSpec

pytestmark = pytest.mark.calsuite("gene_fdr")

SUITE = "gene_fdr"
STREAM = "gene_fdr.mvp04_scored_null"
ALPHA = 0.001
P_ALT = 0.6  # picked_null power is reported against this target-win share
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
    null_rows: tuple[tuple[tuple[float, bool], ...], ...]  # (score, is_decoy) rows, same genes


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
    comp = [c for c in D.gene_competition(v, ex.genes, prots) if not c.has_true]
    wins = tuple(c.target_won for c in comp)
    by_gene: dict[str, list[tuple[float, bool]]] = {}
    for e in D.gene_evidence(v, ex.genes, prots):
        by_gene.setdefault(e.gene, []).append((-e.key[2], e.is_decoy))
    rows = tuple(tuple(sorted(by_gene[c.gene])) for c in comp)
    # the teeth re-draw sides on these rows: they must give the same picked labels (q and PEP
    # are monotone in the score, one mixture)
    assert tuple(_won(r) for r in rows) == wins
    psms = tuple((nl.kind, d, float(r["score"])) for nl, _, r, d in D.null_psms(v))
    return Rep(psms, (cells[0], cells[1], cells[2], cells[3]),
               sum(nl.kind == "protein" for nl in v.nulls.values()), len(v.nulls),
               str(v.provenance["null_split"]["mode"]), pool, tuple(drawn), wins,
               rows)  # fmt: skip


def _won(rows: tuple[tuple[float, bool], ...]) -> bool:
    """Picked label of one gene from its (score, is_decoy) rows: the best score wins, a tie goes
    to the decoy."""
    top = max(s for s, _ in rows)
    return not any(d for s, d in rows if s == top)


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
    for r in reps:  # decision 10, exact (review item 5): no binomial test of a fixed count
        assert r.mode == "fifty_fifty", r.mode
        n, cap = r.n_nulls, len(r.pool)
        assert r.n_protein in {min(n // 2, cap), min((n + 1) // 2, cap)}, (r.n_protein, n, cap)


@dataclass(frozen=True)
class Est:
    """One check on one block: its p value and the estimate it tests, with a 95 % CI."""

    p: float
    what: str
    est: float
    lo: float
    hi: float


Check = Callable[[list[Rep]], Est]


def _cp(k: int, n: int, what: str, p: float) -> Est:
    ci = stats.binomtest(k, n, 0.5).proportion_ci(0.95)  # Clopper-Pearson
    return Est(p, what, k / n, float(ci.low), float(ci.high))


def _mean_ci(x: np.ndarray, what: str, p: float, centre: float = 0.0) -> Est:
    m, se = float(x.mean()) - centre, float(x.std(ddof=1) / np.sqrt(x.size))
    return Est(p, what, m, m - 1.96 * se, m + 1.96 * se)


def _p_coin(kind: str | None) -> Check:
    def p(reps: list[Rep]) -> Est:
        sel = [d for r in reps for k, d, _ in r.psms if kind is None or k == kind]
        k, n = sum(sel), len(sel)
        return _cp(k, n, "decoy share", float(stats.binomtest(k, n, 0.5).pvalue))

    return p


def _p_per_row(reps: list[Rep]) -> Est:
    cells = np.sum([r.cells for r in reps], axis=0).astype(float)
    pv = float(stats.chisquare(cells, cells.sum() * np.array([1, 3, 3, 1]) / 8).pvalue)
    k, n = int(cells[0] + cells[3]), int(cells.sum())
    ci = stats.binomtest(k, n, 0.25).proportion_ci(0.95)
    return Est(pv, "share all-one-side (1/4)", k / n, float(ci.low), float(ci.high))


def _chi_dispersion(reps: list[Rep]) -> float:
    chi = 0.0
    for r in reps:
        n_r = len(r.psms)
        chi += (sum(d for _, d, _ in r.psms) - n_r / 2) ** 2 / (n_r / 4)
    return chi


def _disp(reps: list[Rep], p: float) -> Est:
    chi, df = _chi_dispersion(reps), len(reps)
    lo, hi = chi / stats.chi2.ppf(0.975, df), chi / stats.chi2.ppf(0.025, df)
    return Est(p, "dispersion chi/R (1)", chi / df, float(lo), float(hi))


def _p_disp_low(reps: list[Rep]) -> Est:
    return _disp(reps, float(stats.chi2.cdf(_chi_dispersion(reps), len(reps))))


def _p_disp_high(reps: list[Rep]) -> Est:
    return _disp(reps, float(stats.chi2.sf(_chi_dispersion(reps), len(reps))))


def _scores(reps: list[Rep]) -> tuple[np.ndarray, np.ndarray]:
    ds = np.array([s for r in reps for _, d, s in r.psms if d])
    ts = np.array([s for r in reps for _, d, s in r.psms if not d])
    return ds, ts


def _p_ks2(reps: list[Rep]) -> Est:
    ds, ts = _scores(reps)
    pv = float(stats.ks_2samp(ds, ts).pvalue)
    m = float(ds.mean() - ts.mean())
    se = float(np.sqrt(ds.var(ddof=1) / ds.size + ts.var(ddof=1) / ts.size))
    return Est(pv, "mean score decoy - target", m, m - 1.96 * se, m + 1.96 * se)


def _p_ks_null(reps: list[Rep]) -> Est:
    _, _, mu0, sd0 = (float(x) for x in _cfg().proteomics.score_mixture)
    ds, _ = _scores(reps)
    pv = float(stats.kstest(ds, stats.norm(mu0, sd0).cdf).pvalue)
    return _mean_ci(ds, "mean decoy score - mu0", pv, centre=mu0)


def _p_locus_uniform(reps: list[Rep]) -> Est:
    num = var = 0.0
    n = 0
    for r in reps:
        ks = np.array(list(r.pool.values()), dtype=float)
        for tag in r.drawn:
            num += r.pool[tag] - ks.mean()
            var += ks.var()
            n += 1
    pv = float(2 * stats.norm.sf(abs(num / np.sqrt(var))))
    half = 1.96 * np.sqrt(var) / n
    return Est(pv, "mean k drawn - pool mean k", num / n, num / n - half, num / n + half)


def _p_picked_null(reps: list[Rep]) -> Est:
    w = [x for r in reps for x in r.null_wins]
    k, n = sum(w), len(w)
    return _cp(k, n, "target-win share", float(stats.binomtest(k, n, 0.5).pvalue))


def picked_power(n: int, p_alt: float = P_ALT, alpha: float = ALPHA) -> float:
    """Power of the two-sided exact binomial test at 1/2 (level ``alpha``) on ``n`` null genes
    against a target-win share ``p_alt``."""
    k = np.arange(n + 1)
    pv = np.minimum(1.0, 2 * np.minimum(stats.binom.cdf(k, n, 0.5), stats.binom.sf(k - 1, n, 0.5)))
    return float(stats.binom.pmf(k[pv < alpha], n, p_alt).sum())


def checks(ks_min_p: float) -> dict[str, tuple[Check, float]]:
    """Every p-value gate: name -> (its estimate and p value on a block, the gate p must reach)."""
    return {
        "coin_pooled": (_p_coin(None), ALPHA),
        "coin_protein": (_p_coin("protein"), ALPHA / 2),
        "coin_sixframe": (_p_coin("sixframe"), ALPHA / 2),
        "coin_per_row": (_p_per_row, ALPHA),
        "dispersion_low": (_p_disp_low, ALPHA),
        "dispersion_high": (_p_disp_high, ALPHA),
        "ks_decoy_vs_target": (_p_ks2, ks_min_p),
        "ks_decoy_null": (_p_ks_null, ks_min_p),
        "locus_uniform": (_p_locus_uniform, ALPHA),
        "picked_null": (_p_picked_null, ALPHA),
    }


def _append(records: list[dict[str, Any]], log_dir: Path | None) -> None:
    directory = Path(log_dir) if log_dir is not None else C.log_dir()
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / f"{SUITE}.jsonl").open("a") as fh:  # one write: safe for xdist workers
        fh.write("".join(json.dumps(r) + "\n" for r in records))


def confirm(
    gates: dict[str, tuple[Check, float]], blocks: Blocks, *, log_dir: Path | None = None,
    spoil: Callable[[list[Rep], int], list[Rep]] | None = None,
) -> dict[str, dict[str, Any]]:  # fmt: skip
    """ADR-0015 per check on the harness blocks; one record per block run, appended to the suite
    log with the final status. Returns name -> the last record. ``spoil(reps, block)`` (teeth)
    alters the replicates of each block before the check."""
    out: dict[str, dict[str, Any]] = {}
    for name, (check, gate) in gates.items():
        recs = []
        for k in (1, 2):
            t0 = time.perf_counter()
            reps = blocks.block(k) if spoil is None else spoil(blocks.block(k), k)
            e = check(reps)
            recs.append({
                "suite": SUITE, "stream": blocks.stream, "criterion": "pooled_p", "check": name,
                "level": float(gate), "observed": e.p, "bound": float(gate), "R": blocks.R,
                "seed": C.suite_seed(blocks.stream, (k - 1) * blocks.R),
                "duration": time.perf_counter() - t0, "ok": e.p >= gate, "block": k,
                "estimate": e.est, "ci_low": e.lo, "ci_high": e.hi, "estimand": e.what,
            })  # fmt: skip
            if e.p >= gate:
                break
        status = "pass" if recs[0]["ok"] else ("unconfirmed" if recs[-1]["ok"] else "fail")
        for r in recs:
            r["status"] = status
        _append(recs, log_dir)
        out[name] = recs[-1]
    return out


def _verdict(res: dict[str, dict[str, Any]]) -> None:
    for name, r in res.items():
        print(f"mvp04 calibration {name}: {r['status']} block={r['block']} p={r['observed']:.4g}"
              f" {r['estimand']} = {r['estimate']:.4g} (95% CI {r['ci_low']:.4g},"
              f" {r['ci_high']:.4g})")  # fmt: skip
        if r["status"] == "unconfirmed":
            warnings.warn(f"{name}: block 1 failed, block 2 p={r['observed']:.4g} passes",
                          stacklevel=2)  # fmt: skip
    failed = {k: v for k, v in res.items() if v["status"] == "fail"}
    assert not failed, failed


@pytest.fixture(scope="module")
def ex(tmp_path_factory: pytest.TempPathFactory) -> M.Excerpt:
    return M.make_excerpt(tmp_path_factory.mktemp("k12_excerpt_cal"))


def _whole_blocks(R: int) -> None:
    if C.replicate_range(R) != range(R) or C.confirm_checks() is not None:
        pytest.skip("pooled checks need whole blocks; run this file unsharded")


class Blocks:
    """Harness blocks of stream ``stream``: block k = replicates (k-1)R..kR-1, replicate i with
    seed ``suite_seed(stream, i) mod 2**31``; each block built once, when first asked for."""

    def __init__(self, ex: M.Excerpt, tmp: pytest.TempPathFactory, stream: str, R: int) -> None:
        self.ex, self.tmp, self.stream, self.R = ex, tmp, stream, R
        self._b: dict[int, list[Rep]] = {}

    def seed(self, i: int) -> int:
        return C.suite_seed(self.stream, i) % 2**31

    def block(self, k: int) -> list[Rep]:
        if k not in self._b:
            seeds = [self.seed(i) for i in range((k - 1) * self.R, k * self.R)]
            self._b[k] = [_summary(_rep(self.ex, self.tmp.mktemp(f"mvp04c{s}"), s), self.ex)
                          for s in seeds]  # fmt: skip
        return self._b[k]


@pytest.fixture(scope="module")
def fast(ex: M.Excerpt, tmp_path_factory: pytest.TempPathFactory) -> Blocks:
    return Blocks(ex, tmp_path_factory, STREAM, C.replicates(SUITE, False))


def _power_line(b: list[Rep], R: int) -> float:
    n = sum(len(r.null_wins) for r in b)
    pw = picked_power(n)
    print(f"mvp04 calibration picked_null: {n} null genes per block (R = {R}); power against"
          f" P = {P_ALT} at {ALPHA}: {pw:.3f} per block, {pw**2:.3f} confirmed")  # fmt: skip
    return pw


@pytest.mark.calibration
def test_scored_null_coin_fast(fast: Blocks, tolerances: dict[str, Any]) -> None:
    _whole_blocks(fast.R)
    b1 = fast.block(1)
    _preconditions(b1)
    _power_line(b1, fast.R)
    _verdict(confirm(checks(float(tolerances["ks_min_p"])), fast))


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


def _recoin(reps: list[Rep], seed: int, p_decoy: Callable[[float], float]) -> list[Rep]:
    """Re-draw the side of every null PSM (and of every row of the null-only genes) with
    P(decoy) = ``p_decoy(score)``; the picked labels follow from the re-drawn rows."""
    rng = np.random.default_rng(seed)
    out = []
    for r in reps:
        psms = tuple((k, bool(rng.random() < p_decoy(s)), s) for k, _, s in r.psms)
        rows = tuple(tuple((s, bool(rng.random() < p_decoy(s))) for s, _ in g)
                     for g in r.null_rows)  # fmt: skip
        out.append(replace(r, psms=psms, null_rows=rows, null_wins=tuple(_won(g) for g in rows)))
    return out


def _score_coin(reps: list[Rep], seed: int) -> list[Rep]:
    """Teeth: a score-dependent coin, P(decoy) = 0.3 above mu0 and 0.7 below (pooled share
    stays ~1/2; accepted rows lie above mu0, so targets win ~70 % of the null genes)."""
    mu0 = float(_cfg().proteomics.score_mixture[2])
    return _recoin(reps, seed, lambda s: 0.3 if s > mu0 else 0.7)


def _biased_coin(reps: list[Rep], seed: int) -> list[Rep]:
    """Teeth: a p = 0.55 coin."""
    return _recoin(reps, seed, lambda s: 0.55)


@pytest.mark.calibration
@pytest.mark.parametrize(
    ("name", "spoil"),
    [
        ("coin_per_row", _one_coin_per_peptide),
        ("locus_uniform", _weighted_draw),
        ("ks_decoy_vs_target", _score_coin),
        ("ks_decoy_null", _score_coin),
        ("picked_null", _score_coin),
        ("coin_pooled", _biased_coin),
    ],
)
def test_teeth_are_confirmed_failures(
    fast: Blocks, tolerances: dict[str, Any], teeth_log: Path, name: str,
    spoil: Callable[[list[Rep], int], list[Rep]],
) -> None:  # fmt: skip
    gate = {name: checks(float(tolerances["ks_min_p"]))[name]}
    res = confirm(gate, fast, log_dir=teeth_log, spoil=spoil)
    print(res)
    assert res[name]["status"] == "fail", res


@pytest.mark.calibration
def test_score_coin_keeps_the_pooled_share(fast: Blocks, tolerances: dict[str, Any],
                                           teeth_log: Path) -> None:  # fmt: skip
    """The score-dependent coin is invisible to the pooled binomial: only the KS checks and
    ``picked_null`` catch it (why they carry their own teeth)."""
    gate = {"coin_pooled": checks(float(tolerances["ks_min_p"]))["coin_pooled"]}
    res = confirm(gate, fast, log_dir=teeth_log, spoil=_score_coin)
    assert res["coin_pooled"]["status"] != "fail", res


@pytest.mark.calibration_full
def test_scored_null_coin_full(
    ex: M.Excerpt, tmp_path_factory: pytest.TempPathFactory, tolerances: dict[str, Any],
) -> None:  # fmt: skip
    R = C.replicates(SUITE, True)
    _whole_blocks(R)
    full = Blocks(ex, tmp_path_factory, STREAM + ".full", R)
    b1 = full.block(1)
    _preconditions(b1)
    assert _power_line(b1, R) >= 0.99  # decision 12: picked_null powered against P = 0.6
    _verdict(confirm(checks(float(tolerances["ks_min_p"])), full))
