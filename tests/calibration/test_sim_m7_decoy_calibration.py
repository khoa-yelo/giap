"""Calibration of the M7-12 round-3 six-frame decoys and territories (lead decisions 21 and 22 on
the statistician's round-2 block on PR #225; GIAP_SPEC.md §7.6.6, §8.4 suite "Novel-ORF class
FDR"). The interface is the one in tests/unit/test_sim_m7_round3.py, and the oracle is
tests/helpers/m7_round3_oracle.py.

Claims checked, each pooled over replicates of ``novel_orf_power(seed, size)`` with seeds
``seeds.calibration_base + 5000 + r`` (fast) and ``+ 6000 + r`` (full; ``+ 6500 + r`` at
300 kb):

* **Decoy null (decision 21).** The six-frame decoy scores (``rev_GIAPORF_`` rows) are
  N(mu0, sd0) of ``score_mixture``: one-sample KS, p >= ``ks_min_p``. They are also one
  distribution with the wrong-target scores: two-sample KS, p >= ``ks_min_p``.
* **Coin (decision 21).** Each null PSM is a wrong target row or a six-frame decoy row with
  probability 1/2. Pooled exact binomial test of D among D + T, two-sided, p >= 0.001. The same
  test per makedb entry class (targets classed by the entry their locus lies in, decoys by their
  ID), p >= 0.001 / 4 (Bonferroni over the four classes).
* **Coin, dispersion.** Across replicates, sum_r (D_r - N_r / 2)^2 / (N_r / 4) is chi-square with
  R degrees of freedom. Both tails must be >= 0.001. A deterministic one-to-one split (statistic 0)
  fails the lower tail.
* **Territories (decision 22).** Pooled over replicates, each of the four classes holds at least
  0.25 - 3 SE of the six-frame wrong peptides. The SE is the binomial SE at 1/4: the draw takes
  the classes in turn, so the realised share is at least this.

Fast: 4 replicates at 100 kb (each about 2 s, plus 0.5 s for makedb). Full (``calibration_full``):
``replicates_fast.novel_orfs`` replicates at 100 kb plus 3 at 300 kb. Suite ``novel_orfs``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from helpers import m7_round2_oracle as R
from helpers import m7_round3_oracle as Q
from helpers import m7_sim_oracle as O
from scipy import stats

import giapsim

pytestmark = pytest.mark.calsuite("novel_orfs")

ALPHA = 0.001
N_FAST = 4


@dataclass
class Rep:
    t_scores: list[float] = field(default_factory=list)
    d_scores: list[float] = field(default_factory=list)
    t_by: dict[str, int] = field(default_factory=dict)
    d_by: dict[str, int] = field(default_factory=dict)
    territory: dict[str, int] = field(default_factory=dict)


def _rep(tmp: Path, seed: int, size: str) -> Rep:
    import giapsim.scenarios as sc

    out = tmp / "sim"
    giapsim.simulate(sc.novel_orf_power(seed, size), seed, out)
    entries = Q.makedb_entries(out, tmp / "db")
    index = Q.EntryIndex(entries, Q.lengths(out))
    loci = Q.wrong_loci(out)
    t, d6, _ = Q.split_rows(out)
    rep = Rep()
    for r in t:
        rep.t_scores.append(float(r["score"]))
        s = loci[r["peptide"]]
        hits = index.entry_of(s.contig, s.s1, s.e1, s.strand)
        assert len(hits) == 1, (r, hits)
        rep.t_by[hits[0].cls] = rep.t_by.get(hits[0].cls, 0) + 1
    for r in d6:
        rep.d_scores.append(float(r["score"]))
        m = Q.DECOY_ID.match(Q.ids(r["protein_ids"])[0])
        assert m is not None, r
        rep.d_by[m.group(1)] = rep.d_by.get(m.group(1), 0) + 1
    genome = O.receiving_genome(out)
    cds, genes = R.cds_spans(out), Q.true_starts(out)
    for span in loci.values():
        c = Q.wrong_class4(span, cds, genes, genome)
        rep.territory[c] = rep.territory.get(c, 0) + 1
    return rep


def _check(reps: list[Rep], ks_min_p: float, mix: tuple[float, ...]) -> None:
    _, _, mu0, sd0 = mix
    ds = np.concatenate([r.d_scores for r in reps])
    ts = np.concatenate([r.t_scores for r in reps])
    assert len(ds) >= 100 and len(ts) >= 100, (len(ds), len(ts))
    assert np.isfinite(ds).all()
    assert stats.kstest(ds, stats.norm(mu0, sd0).cdf).pvalue >= ks_min_p
    assert stats.ks_2samp(ds, ts).pvalue >= ks_min_p
    # coin, pooled and per class
    d, n = len(ds), len(ds) + len(ts)
    assert stats.binomtest(d, n, 0.5).pvalue >= ALPHA, (d, n)
    for c in Q.CLASSES4:
        dc = sum(r.d_by.get(c, 0) for r in reps)
        nc = dc + sum(r.t_by.get(c, 0) for r in reps)
        assert nc >= 20, (c, nc)
        assert stats.binomtest(dc, nc, 0.5).pvalue >= ALPHA / 4, (c, dc, nc)
    # coin, dispersion across replicates
    chi = sum((len(r.d_scores) - (len(r.d_scores) + len(r.t_scores)) / 2) ** 2
              / ((len(r.d_scores) + len(r.t_scores)) / 4) for r in reps)  # fmt: skip
    df = len(reps)
    assert stats.chi2.cdf(chi, df) >= ALPHA, ("split looks deterministic", chi, df)
    assert stats.chi2.sf(chi, df) >= ALPHA, ("split overdispersed", chi, df)
    # territories
    tot = {c: sum(r.territory.get(c, 0) for r in reps) for c in Q.CLASSES4}
    m = sum(tot.values())
    se = math.sqrt(0.25 * 0.75 / m)
    for c in Q.CLASSES4:
        assert tot[c] / m >= 0.25 - 3 * se, tot


def _mix() -> tuple[float, ...]:
    import giapsim.scenarios as sc

    return tuple(float(x) for x in sc.novel_orf_power(0).proteomics.score_mixture)


@pytest.fixture(scope="module")
def base(tolerances: dict[str, Any]) -> int:
    return int(tolerances["seeds"]["calibration_base"])


@pytest.mark.calibration
def test_sixframe_decoys_null_and_coin_fast(
    tmp_path_factory: pytest.TempPathFactory, base: int, tolerances: dict[str, Any]
) -> None:
    reps = [_rep(tmp_path_factory.mktemp(f"m7dec{r}"), base + 5000 + r, "100k")
            for r in range(N_FAST)]  # fmt: skip
    _check(reps, float(tolerances["ks_min_p"]), _mix())


@pytest.mark.calibration_full
def test_sixframe_decoys_null_and_coin_full(
    tmp_path_factory: pytest.TempPathFactory, base: int, tolerances: dict[str, Any]
) -> None:
    n = int(tolerances["replicates_fast"]["novel_orfs"])
    reps = [_rep(tmp_path_factory.mktemp(f"m7decF{r}"), base + 6000 + r, "100k")
            for r in range(n)]  # fmt: skip
    reps += [_rep(tmp_path_factory.mktemp(f"m7decM{r}"), base + 6500 + r, "300k")
             for r in range(3)]  # fmt: skip
    _check(reps, float(tolerances["ks_min_p"]), _mix())
