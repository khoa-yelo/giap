"""Attack nulls for the part-B TerminatorNet composites (GIAP_SPEC.md §8.5; PR #264 statistician
review item 10). This ports ``local/logs/stat264_probe_tn.py``, cut down to the fast tier.

The composites are wired the way ``validation/calibrate_cutoffs.py`` wires them:
- c = the mean of one ``coincidence_rate(..., return_rates=True)`` call;
- SE = sd(ddof=1) / √R (the lead's decision of 2026-10-23);
- the window CI is the precision half of ``fixed_cutoff_ci`` with the call's window gene as
  the cluster;
- each bound is ``precision_lower_bound(<lower limit>, c)``.

World: a circular replicon with genes on both strands and windows of 150-3000 bp. There are
600 real terminators at |N(60, 120)| bp past the window's 3' side, and a share m of them is
listed (RegulonDB incompleteness). Calls are made at one confidence:
- True calls: clusters of 1-3 calls at 35 % of the real terminators. A call at an unlisted
  terminator does not match.
- False calls: clusters of 1-4 calls sharing an offset in a random window.
The true precision π is the share of true calls.

Claims (gated):
1. P(lb_window <= π) is within the one-sided binomial band of 0.95, that is
   >= 0.95 - k·√(0.95·0.05/n), pooled over m ∈ {1.0, 0.6}. The probe saw 1.000 over 200
   worlds per m at R = 200, B = 2000.
2. The window CI is wider than Wilson on average, because false calls cluster in windows.
3. The per-placement MC SE is within [0.6, 1.5] times the SD of c across 20 rng seeds on one world.

The point bound is a point estimate, not a confidence limit; its exceedance share is printed only.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene
from giap.rng import for_step
from giap.stats.cutoffs import (
    Window,
    coincidence_rate,
    downstream_windows,
    fixed_cutoff_ci,
    precision_lower_bound,
    terminator_match,
    wilson_ci,
    window_of,
)

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("regions")]

N_WORLDS = 12  # per completeness level; ~1.6 s each
B_BOOT = 200
R_PLACE = 200
LEVEL = 0.95


def world(
    g: np.random.Generator, m: float, n_genes: int = 1800
) -> tuple[Replicon, list[Window], list[Interval], list[Interval], np.ndarray]:
    gaps = g.integers(100, 3000, n_genes)
    lens = g.integers(300, 1500, n_genes)
    strands = g.choice(["+", "-"], n_genes)
    pos, genes = 0, []
    for i in range(n_genes):
        genes.append(
            Gene(
                locus_tag=f"g{i:04d}",
                feature_type="cds",
                interval=Interval("ring", int(pos), int(pos + lens[i]), str(strands[i])),  # type: ignore[arg-type]
                aa="M" * 99,
                nt=None,
                start_codon="ATG",
            )
        )
        pos += int(lens[i] + gaps[i])
    rep = Replicon("ring", int(pos), "circular")
    ws = [w for w in downstream_windows(genes, rep) if w.hi - w.lo >= 150]

    def place(w: Window, off: float, ln: int) -> Interval:
        o = int(min(max(off, 0), w.hi - w.lo - ln))
        s = (w.lo + o if w.strand == "+" else w.hi - o - ln) % rep.length
        return Interval("ring", s, s + ln, w.strand)  # type: ignore[arg-type]

    real_w = g.choice(len(ws), 600, replace=False)
    real = [place(ws[k], abs(g.normal(60, 120)), 25) for k in real_w]
    listed = g.random(600) < m
    truth = [t for t, keep in zip(real, listed, strict=True) if keep]
    calls: list[Interval] = []
    is_true: list[bool] = []
    for t in real:
        if g.random() < 0.35:
            for _ in range(int(g.integers(1, 4))):
                s = (t.start + int(g.integers(-15, 16))) % rep.length
                calls.append(Interval("ring", s, s + 30, t.strand))
                is_true.append(True)
    for _ in range(int(g.integers(230, 260))):
        w = ws[int(g.integers(0, len(ws)))]
        off = abs(g.normal(0, 300))
        for _ in range(int(g.integers(1, 5))):
            calls.append(place(w, off + int(g.integers(-20, 21)), 30))
            is_true.append(False)
    keep = [i for i, c in enumerate(calls) if window_of(c, ws, rep) is not None]
    return rep, ws, truth, [calls[i] for i in keep], np.array([is_true[i] for i in keep])


def composites(g: np.random.Generator, m: float) -> tuple[float, float, float, float, float, float]:
    """(π, lb_point, lb_wilson, lb_window, Wilson width, window width) for one world."""
    rep, ws, truth, calls, is_true = world(g, m)
    y = np.array([terminator_match(c, truth, rep) is not None for c in calls], dtype=int)
    gene = np.array([window_of(c, ws, rep) for c in calls])
    conf = np.full(len(calls), 0.9)
    rates = coincidence_rate(calls, conf, ws, truth, rep, 0.9, rng=g, R=R_PLACE, return_rates=True)
    c = float(np.mean(rates))
    win, _ = fixed_cutoff_ci(conf, y, gene, 0.9, rng=g, B=B_BOOT)
    wl = wilson_ci(int(y.sum()), len(y))
    return (
        float(is_true.mean()),
        precision_lower_bound(float(y.mean()), c),
        precision_lower_bound(wl[0], c),
        precision_lower_bound(win.ci_low, c),
        wl[1] - wl[0],
        win.ci_high - win.ci_low,
    )


def test_window_bound_is_a_lower_confidence_limit_under_clustering(
    tolerances: dict[str, object],
) -> None:
    k = float(tolerances["se_multiplier"])  # type: ignore[arg-type]
    seeds = tolerances["seeds"]
    assert isinstance(seeds, dict)
    g = for_step("calibration.m6_12.attack_driver", int(seeds["calibration_base"]))
    res = np.array([composites(g, m) for m in (1.0, 0.6) for _ in range(N_WORLDS)])
    pi, lbp, lbw, lbb, w_wil, w_win = res.T
    n = len(pi)
    cov_win, cov_wil, exceed = (
        float(np.mean(lbb <= pi)),
        float(np.mean(lbw <= pi)),
        float(np.mean(lbp > pi)),
    )
    ratio = float(np.mean(w_win / w_wil))
    print(
        f"M6-12 part-B attack (n={n}): P(lb_window<=pi)={cov_win:.3f} "
        f"P(lb_wilson<=pi)={cov_wil:.3f} P(lb_point>pi)={exceed:.3f} width ratio={ratio:.3f}"
    )
    floor = LEVEL - k * math.sqrt(LEVEL * (1 - LEVEL) / n)
    assert cov_win >= floor, f"window bound covers pi in {cov_win:.3f} < {floor:.3f}"
    assert ratio > 1.0, f"window CI not wider than Wilson under clustering ({ratio:.3f})"


def test_per_placement_mc_se_matches_the_seed_spread(tolerances: dict[str, object]) -> None:
    seeds = tolerances["seeds"]
    assert isinstance(seeds, dict)
    base = int(seeds["calibration_base"])
    rep, ws, truth, calls, _ = world(for_step("calibration.m6_12.attack_driver_se", base), 1.0)
    conf = np.full(len(calls), 0.9)
    cs, ses = [], []
    for s in range(20):
        rng = for_step("calibration.m6_12.attack_driver_se", base + 1 + s)
        rates = coincidence_rate(
            calls, conf, ws, truth, rep, 0.9, rng=rng, R=R_PLACE, return_rates=True
        )
        cs.append(float(np.mean(rates)))
        ses.append(float(np.std(rates, ddof=1) / math.sqrt(R_PLACE)))
    sd = float(np.std(cs, ddof=1))
    ratio = float(np.mean(ses)) / sd
    print(f"M6-12 part-B MC SE: mean SE {np.mean(ses):.5f}, sd(c) over 20 seeds {sd:.5f}")
    assert 0.6 <= ratio <= 1.5, f"MC SE / seed SD = {ratio:.3f}"
