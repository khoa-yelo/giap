"""Attack C through the step: the one-detected P(miss) arm of M6-05H under structured nulls
(GIAP_SPEC.md §7.4.3 "One detected, one not", §8.1 row "Operon discordance", §8.5; M6-09 attack C,
which ``test_attack_operons.py`` records with the *oracle* P(miss); M6-index lead decision 1).

Statistician's attack. ``test_attack_operons.py::test_attack_c_through_the_step`` only probes that
the step makes ``one_detected_pmiss`` rows; here the M6-02 model is *fitted* and P(miss) comes from
the step's own chain: ``giap.steps.operons.detection_limit`` (``fit_detection_limit`` with the
``operons.detection_limit`` RNG stream, ``p_miss_for_genes``) feeding
``giap.operons.discordance.discordance`` - the call ``run_operons`` makes. Only ``pair_frame`` is
bypassed: the simulated pair table is used directly, as everywhere in the M6-09 gate.

Data. Replicate ``i`` draws one 1,800-gene genome with ``giapsim.simulate_operons`` (TU blocks,
the ``OPERON_T`` within-TU walk, so the operon Delta are exactly the generating t) and redraws
detection at the *peptide* level so that M6-02's model family can hold:

- every observable peptide j of gene g (the spec item-1 digest of the simulated protein, 7-30 aa)
  is seen with logit = ``B0 + 0.5 a_g - 0.03 L_j + 0.8 GRAVY_j + o_g``; g is detected when >= 1
  peptide is seen (tier ``mayu``, k_min = 1); C11 holds the seen peptides; undetected genes have
  ``abundance_log2`` NA (as in a real gene table);
- arms (the offset o_g, which the model never sees):
  - ``null``: ``operon_null``, o_g = 0. M6-02 is correctly specified; the only departures from
    the ideal P(miss) are estimation and the ratio sample (Delta of both-detected pairs only);
  - ``gene_effect`` / ``gene_effect_05``: ``operon_null``, o_g ~ N(0, 1) / N(0, 0.5^2) per gene,
    shared by its peptides (protein-level detectability: ionisation, membrane, digestion - the §8.5
    "detectability" structure);
  - ``block``: ``operon_attack_detection(1800, "prophage")`` (ADR-0018 layout: a third of the
    chromosome with short genes), o_g = ``BLOCK_OFFSET`` inside the block: spatially clustered,
    TU-correlated detectability.

Quantities (pooled over replicates; alpha = the tolerances levels):

- *conditional share*: the share of ``one_detected_pmiss`` rows with P(miss) <= alpha (the M6-09
  attack C metric; under the complete null every such flag is false);
- *at-risk rate*: over every orientation (u, v) of an OP pair with v detected and u scorable, the
  share with u missed **and** P(miss) <= alpha. This is the level of P(miss) used as a discrete
  p-value of the event "u missed"; it is <= alpha when P(miss) is calibrated;
- *calibration*: logistic recalibration (miss ~ a + b logit P(miss)) of the step's P(miss) over
  the same orientations, the mean P(miss) against the miss rate, and the same for the oracle
  (true model, offset included, Delta from the generating t) and for the true model (offset
  included) with the step's ratio sample (which isolates the both-detected selection of the
  ratios);
- *mixture*: the per-pair rate of p <= alpha over every tested OP pair when the P(miss) rows join
  the t rows, and the share of replicates with any BH call at q <= 0.1 for the t family (shipped)
  and for the pooled family (the spec's literal "both kinds together"). Under this complete null
  the any-call share is the FWER and the FDR.

Recorded (printed and logged to ``operons_attack_c_step.jsonl``); the escalation thresholds fail
the test so a breach cannot pass silently: the conditional share above 2 alpha (the M6-09 C
threshold) on the ``null`` arm, the at-risk rate above alpha + 3 SE on the ``null`` arm, a step
P(miss) that differs from the vectorised recomputation, and the shipped t family's any-call share
above its own bound. The misspecified arms (``gene_effect*``, ``block``) are recorded only: they
measure how far P(miss) moves when detectability has structure the model does not see, which is
the evidence for keeping P(miss) out of the BH family (ADR-0056, proposed).
R = 20 fast / 300 full (module constants, as M6-09's ``R_C``).
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from functools import cache
from typing import Any

import _criteria as C
import _operon_gate as G
import numpy as np
import pandas as pd
import pytest
from scipy.special import log_expit, logit

from giapsim.scenarios import OPERON_T

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
DI = "discordance"
R_STEP = (20, 300)
ARMS = ("null", "gene_effect", "gene_effect_05", "block")
#: per-peptide truth; B0 puts gene-level detection near 0.75 (null arm) at the simulated abundances
BETA = {"intercept": -11.5, "abundance_log2": 0.5, "length": -0.03, "gravy": 0.8}
GENE_EFFECT_SD = {"gene_effect": 1.0, "gene_effect_05": 0.5}
BLOCK_OFFSET = -1.5
N_ORACLE_T = 2000
Q_FLAG = 0.1
HIT_COLS = ("peptide", "gene")


def _r(full: bool) -> int:
    return R_STEP[1] if full else R_STEP[0]


@cache
def _oracle_delta() -> np.ndarray:
    loc, scale, df = OPERON_T
    rng = C.replicate_rng(f"{DI}.attack_c.step.oracle_t", 0)
    return np.asarray(loc + scale * rng.standard_t(df, N_ORACLE_T), dtype=float)


def _observable(aa: Any, lo: int, hi: int, enzyme: str) -> list[str]:
    from giap.genetable.covariates import digest

    if not isinstance(aa, str) or not aa:
        return []
    return list(dict.fromkeys(p for p in digest(aa, enzyme) if lo <= len(p) <= hi))


def _gravy(pep: str) -> float:
    from giap.operons.detection_limit import _KD

    v = [_KD[c] for c in pep if c in _KD]
    return sum(v) / len(v) if v else 0.0


def _pmiss_vec(coef: dict[str, float], a_v: float, r: np.ndarray, x: np.ndarray,
               offset: float = 0.0) -> float:  # fmt: skip
    """mean_r prod_j (1 - p_j(a_v + r)): P(miss) at k_min = 1, vectorised."""
    eta = (coef["intercept"] + offset + coef["abundance_log2"] * (a_v + r)[:, None]
           + coef["length"] * x[None, :, 0] + coef["gravy"] * x[None, :, 1])  # fmt: skip
    return float(np.mean(np.exp(log_expit(-eta).sum(axis=1))))


@dataclass
class StepRep:
    pmiss: np.ndarray  # one_detected_pmiss p (the step's)
    t_p: np.ndarray  # both_detected_t p
    t_q: np.ndarray
    pi_step: np.ndarray  # per at-risk orientation
    pi_oracle: np.ndarray
    pi_true_r: np.ndarray
    miss: np.ndarray
    max_abs_diff: float  # step P(miss) vs the vectorised recomputation
    det_rate: float
    usable: bool
    auc: float


def _arm_cfgs(arm: str) -> Any:
    if arm == "block":
        return G.scenario("operon_attack_detection", G.N_GENES, "prophage")
    return G.scenario("operon_null", G.N_GENES)


def step_replicate(arm: str, i: int) -> StepRep:
    from giap.config import GiapConfig
    from giap.operons import detection_limit as dl
    from giap.operons.discordance import discordance
    from giap.steps.operons import detection_limit
    from giapsim.expression import detection_terms

    stream = f"{DI}.attack_c.step.{arm}"
    seed = C.suite_seed(stream, i)
    cfgs = _arm_cfgs(arm)
    sim = G.simulate(cfgs, seed)
    rng = C.replicate_rng(f"{stream}.peptides", i)
    cfg = GiapConfig(seed=int(seed % (2**31 - 1)))
    lo, hi = cfg.abundance.observable_len
    lg = sim.layout.genes.reset_index(drop=True)
    eg = sim.expression.genes.set_index("locus_tag").loc[lg["locus_tag"]]
    tags = lg["locus_tag"].astype(str).to_numpy()
    a_true = eg["abundance_log2_true"].to_numpy(float)
    n = len(tags)
    if arm in GENE_EFFECT_SD:
        off = rng.normal(0.0, GENE_EFFECT_SD[arm], n)
    elif arm == "block":
        _, shift, _, _ = detection_terms(sim.layout, cfgs[1])
        off = np.where(np.asarray(shift) != 0.0, BLOCK_OFFSET, 0.0)
    else:
        off = np.zeros(n)
    tu_of = {str(g): str(t) for t, gs in zip(sim.layout.tus["tu_id"], sim.layout.tus["genes"],
                                             strict=True) for g in gs}  # fmt: skip
    x_of: dict[str, np.ndarray] = {}
    det = np.zeros(n, dtype=bool)
    hits: list[tuple[str, str]] = []
    for k, (tag, aa) in enumerate(zip(tags, lg["aa"], strict=True)):
        seqs = _observable(aa, lo, hi, cfg.abundance.enzyme)
        x = np.array([(len(s), _gravy(s)) for s in seqs], dtype=float).reshape(-1, 2)
        x_of[tag] = x
        if not seqs:
            continue
        eta = (BETA["intercept"] + BETA["abundance_log2"] * a_true[k] + off[k]
               + BETA["length"] * x[:, 0] + BETA["gravy"] * x[:, 1])  # fmt: skip
        seen = rng.random(len(seqs)) < 1.0 / (1.0 + np.exp(-eta))
        det[k] = bool(seen.any())
        hits.extend((s, tag) for s, o in zip(seqs, seen, strict=True) if o)
    from helpers.detection_limit_data import hit_row, hits_frame

    c11 = hits_frame([hit_row(p, g) for p, g in hits])
    obs = np.where(det, a_true, np.nan)
    gt = pd.DataFrame({
        "locus_tag": pd.Series(tags, dtype=object),
        "tu_id": pd.Series([tu_of.get(t) for t in tags], dtype=object),
        "detected": det,
        "abundance_log2": obs,
        "detection_tier": "mayu",
    })  # fmt: skip
    fit, p_miss = detection_limit(gt, lg[["locus_tag", "aa"]], c11, cfg)
    pairs = sim.pairs.copy()
    pos = {t: k for k, t in enumerate(tags)}
    for side in ("a", "b"):
        ix = pairs[f"gene_{side}"].astype(str).map(pos).to_numpy(int)
        pairs[f"detected_{side}"] = det[ix]
        pairs[f"abundance_log2_{side}"] = obs[ix]
    pairs["delta_log2"] = pairs["abundance_log2_b"] - pairs["abundance_log2_a"]
    genes = pd.DataFrame({"locus_tag": tags})
    res = discordance(pairs, genes, cfg, p_miss)
    rows = res.rows
    tt = rows["test"].to_numpy()
    pm_rows = rows[tt == "one_detected_pmiss"]
    ratios = np.asarray(res.ratios, dtype=float)
    up = -ratios[::-1]
    od = _oracle_delta()
    m = fit.model
    pis: list[tuple[float, float, float, bool]] = []
    diffs = [0.0]
    step_pm = dict(zip(pm_rows["pair_id"].astype(str), pm_rows["p_value"].astype(float),
                       strict=True))  # fmt: skip
    op = pairs[pairs["class"] == "OP"]
    if p_miss is not None and ratios.size:
        for row in op.itertuples(index=False):
            for down in (True, False):  # u = gene_b downstream (v = gene_a), or u = gene_a
                v, u = (row.gene_a, row.gene_b) if down else (row.gene_b, row.gene_a)
                kv, ku = pos[str(v)], pos[str(u)]
                if not det[kv] or x_of[str(u)].shape[0] == 0:
                    continue
                coef = dict(dl.model_for_tu(m, tu_of.get(str(u), "")).coef)
                a_v = float(obs[kv])
                r = ratios if down else up
                xs = x_of[str(u)]
                pi = _pmiss_vec(coef, a_v, r, xs)
                pis.append((pi, _pmiss_vec(BETA, a_v, od if down else -od, xs, off[ku]),
                            _pmiss_vec(BETA, a_v, r, xs, off[ku]), not det[ku]))  # fmt: skip
                if not det[ku] and str(row.pair_id) in step_pm:
                    diffs.append(abs(step_pm[str(row.pair_id)] - pi))
    arr = np.array(pis, dtype=float).reshape(-1, 4)
    return StepRep(
        pmiss=pm_rows["p_value"].astype(float).to_numpy(),
        t_p=rows.loc[tt == "both_detected_t", "p_value"].astype(float).to_numpy(),
        t_q=rows.loc[tt == "both_detected_t", "q_value"].astype(float).to_numpy(),
        pi_step=arr[:, 0], pi_oracle=arr[:, 1], pi_true_r=arr[:, 2], miss=arr[:, 3].astype(bool),
        max_abs_diff=float(max(diffs)), det_rate=float(det.mean()), usable=bool(m.usable),
        auc=float(m.auc),
    )  # fmt: skip


def _recal(miss: np.ndarray, pi: np.ndarray) -> tuple[float, float]:
    """(intercept, slope) of the logistic recalibration miss ~ a + b logit(pi)."""
    import statsmodels.api as sm

    z = logit(np.clip(pi, 1e-9, 1 - 1e-9))
    res = sm.GLM(miss.astype(float), sm.add_constant(z), family=sm.families.Binomial()).fit()
    return float(res.params[0]), float(res.params[1])


def summarise(reps: list[StepRep], levels: tuple[float, ...]) -> dict[str, Any]:
    from giap.stats.multitest import bh

    lv = np.asarray(levels)
    pm = np.concatenate([r.pmiss for r in reps])
    pi = np.concatenate([r.pi_step for r in reps])
    po = np.concatenate([r.pi_oracle for r in reps])
    pt = np.concatenate([r.pi_true_r for r in reps])
    miss = np.concatenate([r.miss for r in reps])
    tp = np.concatenate([r.t_p for r in reps])
    n_risk = pi.size
    at_risk = np.array([float(((pi <= a) & miss).sum() / n_risk) for a in lv])
    at_risk_oracle = np.array([float(((po <= a) & miss).sum() / n_risk) for a in lv])
    cond = np.array([float((pm <= a).mean()) for a in lv])
    mix = np.array([float(((tp <= a).sum() + (pm <= a).sum()) / (tp.size + pm.size)) for a in lv])
    any_t = float(np.mean([bool((r.t_q <= Q_FLAG).any()) for r in reps]))
    pooled = [bh(np.concatenate([r.t_p, r.pmiss])) for r in reps]
    any_pool = float(np.mean([bool((q <= Q_FLAG).any()) for q in pooled]))
    n_pool = float(np.mean([int((q <= Q_FLAG).sum()) for q in pooled]))
    n_flag = [float(np.mean([int((r.pmiss <= a).sum()) for r in reps])) for a in lv]
    deciles = np.quantile(pi, np.linspace(0, 1, 11))
    bins = np.clip(np.searchsorted(deciles, pi, side="right") - 1, 0, 9)
    rel = [(float(pi[bins == b].mean()), float(miss[bins == b].mean())) for b in range(10)
           if (bins == b).any()]  # fmt: skip
    return {
        "R": len(reps),
        "usable": float(np.mean([r.usable for r in reps])),
        "auc_median": float(np.nanmedian([r.auc for r in reps])),
        "det_rate": float(np.mean([r.det_rate for r in reps])),
        "n_pmiss": int(pm.size),
        "n_at_risk": int(n_risk),
        "max_abs_diff": float(max(r.max_abs_diff for r in reps)),
        "levels": lv.tolist(),
        "conditional_share": cond.tolist(),
        "flags_per_genome": n_flag,
        "at_risk_rate": at_risk.tolist(),
        "at_risk_rate_oracle": at_risk_oracle.tolist(),
        "mean_pi": [float(pi.mean()), float(pt.mean()), float(po.mean())],
        "miss_rate": float(miss.mean()),
        "recal_step": _recal(miss, pi),
        "recal_true_r": _recal(miss, pt),
        "recal_oracle": _recal(miss, po),
        "reliability_deciles": rel,
        "pmiss_quantiles": (
            np.quantile(pm, [0.01, 0.05, 0.1, 0.25, 0.5]).tolist() if pm.size else []
        ),
        "mixture_rate": mix.tolist(),
        "any_call_t": any_t,
        "any_call_pooled": any_pool,
        "mean_calls_pooled": n_pool,
    }


@pytest.mark.calsuite(DI)
@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize("arm", ARMS)
def test_recorded_attack_c_through_the_step(arm: str, full: bool) -> None:
    if os.environ.get(C.CONFIRM_ENV) is not None:
        pytest.skip("recorded attack, not a calibrate() check; nothing to confirm")
    R = _r(full)
    if C.replicate_range(R).start != 0:
        pytest.skip("runs whole in the shard holding replicate 0")
    lv = G.levels()
    reps = [step_replicate(arm, i) for i in range(R)]
    s = summarise(reps, lv)
    s.update({"arm": arm, "tier": "full" if full else "fast"})
    print(f"C through the step [{arm}]: " + json.dumps(s))
    path = C.log_dir() / "operons_attack_c_step.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(s) + "\n")
    assert s["usable"] == 1.0 and s["n_pmiss"] > 0
    assert s["max_abs_diff"] < 1e-9, "step P(miss) differs from the k_min = 1 recomputation"
    # the shipped family is the t family; its any-call share stays within the M6-09 bound
    se = math.sqrt(Q_FLAG * (1 - Q_FLAG) / R)
    assert s["any_call_t"] <= Q_FLAG + 2 * se + 1 / R
    if arm == "null":
        lvl = np.asarray(lv)
        cond = np.asarray(s["conditional_share"])
        assert np.all(cond <= 2 * lvl), f"C above 2 alpha through the step: {cond}"
        se_r = np.sqrt(lvl * (1 - lvl) / s["n_at_risk"])
        risk = np.asarray(s["at_risk_rate"])
        assert np.all(risk <= lvl + 3 * se_r), f"P(miss) above level as a p-value: {risk}"
