"""giapsim operon scenarios: the two lead decisions the PR #222 review found unpinned (issue M6-03,
docs/specs/M6-03-sim-operon-scenarios.md).

* **B1, decision 10.** In ``operon_attack_detection`` the detection logit is recomputed on the final
  (post-walk) abundance *and* ``detected_true`` is redrawn from it. The existing test pins only
  ``p_detect_true``. Here ``detected_true`` must be a Bernoulli draw from ``p_detect_true``: pooled
  over seeds 0-2 (5,400 genes), genes are binned by deciles of p, and the per-decile
  (observed - expected) detected counts, each scaled by its binomial SD ``sqrt(sum p(1-p))``, give
  a statistic that is chi-square with 10 df under the claim. The test rejects at chi2 > 29.6
  (alpha = 0.001). Measured on b2734b18 for the seed triples (0-2), (3-5), ..., (12-14): pooled
  chi2 of 4.1-8.6; keeping the pre-walk draw gives 544-1058.
* **B2, decision 5.** The planted co-abundant pairs of ``missed_operon_power`` have
  ``Delta = effect_sd * SD * z``, where SD is the background Delta SD of the both-detected NO
  pairs, so SD(planted Delta) / SD(unplanted both-detected true-NO Delta) is about ``effect_sd``.
  The test averages the ratio over seeds 0-2 at effect_sd = 0.1 and 0.5 and requires
  ratio / effect_sd to lie in [0.7, 1.4]. Measured on b2734b18 for the seed triples (0-2), (3-5),
  ..., (12-14): 0.82-1.06 at both sizes (seeds 0-2: 1.06). Without the ``* SD`` factor the effect
  is in log2 units and the ratio drops to about 0.33 x effect_sd.
"""

from __future__ import annotations

from functools import cache
from typing import Any

import numpy as np
import pandas as pd
import pytest

N = 1800
SEEDS = (0, 1, 2)
CHI2_10DF_CRIT = 29.588  # scipy.stats.chi2.isf(0.001, 10)


@cache
def _sim(name: str, args: tuple[Any, ...], seed: int) -> Any:
    import giapsim.scenarios as sc

    return sc.simulate_operons(getattr(sc, name)(*args), seed)


def _b(s: pd.Series) -> np.ndarray:
    return np.asarray(s.astype(bool).to_numpy(), dtype=bool)


def _true_delta(s: Any, p: pd.DataFrame) -> np.ndarray:
    g = s.expression.genes
    ab = dict(zip(g["locus_tag"].astype(str), g["abundance_log2_true"].astype(float), strict=True))
    return np.array(
        [ab[str(b)] - ab[str(a)] for a, b in zip(p["gene_a"], p["gene_b"], strict=True)]
    )


def test_b1_detection_attack_detected_true_is_drawn_from_p_detect_true() -> None:
    """Decision 10: pooled over seeds 0-2, the detected count in each decile of ``p_detect_true``
    matches sum(p) within binomial error (chi2 over the 10 deciles < 29.6, alpha = 0.001)."""
    p_all, d_all = [], []
    for seed in SEEDS:
        g = _sim("operon_attack_detection", (N, "prophage"), seed).expression.genes
        p_all.append(g["p_detect_true"].to_numpy(float))
        d_all.append(_b(g["detected_true"]))
    p, d = np.concatenate(p_all), np.concatenate(d_all)
    assert np.all((p >= 0) & (p <= 1))
    edges = np.quantile(p, np.linspace(0, 1, 11))
    k = np.clip(np.searchsorted(edges, p, side="right") - 1, 0, 9)
    chi2 = 0.0
    for j in range(10):
        m = k == j
        assert m.sum() > 100, (j, m.sum())
        var = float(np.sum(p[m] * (1 - p[m])))
        assert var > 0.0, (j, var)  # the top decile has p near 1, so var can be below 1
        chi2 += (float(d[m].sum()) - float(p[m].sum())) ** 2 / var
    assert chi2 < CHI2_10DF_CRIT, chi2


@pytest.mark.parametrize("effect_sd", [0.1, 0.5])
def test_b2_planted_effect_is_in_units_of_the_background_sd(effect_sd: float) -> None:
    """Decision 5: SD(planted Delta) / SD(unplanted both-detected true-NO Delta), averaged over
    seeds 0-2, is within [0.7, 1.4] x effect_sd."""
    ratios = []
    for seed in SEEDS:
        s = _sim("missed_operon_power", (N, 50, effect_sd), seed)
        p = s.pairs
        planted = p["pair_id"].isin(set(s.truth["missed_operons"]["pair_id"])).to_numpy()
        assert planted.sum() == 50
        both = _b(p["detected_a"]) & _b(p["detected_b"])
        bg = (p["class"] == "NO").to_numpy() & ~_b(p["true_operon"]) & both & ~planted
        assert bg.sum() >= 100, bg.sum()
        sd_planted = float(np.std(_true_delta(s, p[planted]), ddof=1))
        sd_bg = float(np.std(_true_delta(s, p[bg]), ddof=1))
        ratios.append(sd_planted / sd_bg)
    rel = float(np.mean(ratios)) / effect_sd
    assert 0.7 <= rel <= 1.4, (rel, ratios)
