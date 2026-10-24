"""Attack nulls of the M6 gate (issue M6-09; GIAP_SPEC.md §7.4.7 "Calibration" criterion 2, §8.5;
docs/specs/M6-09-calibration-gate.md, its lead defaults, and the attacks owed by the PR #147 / #148
statistician reviews in docs/methods/operons.md). Helpers: ``../_operon_gate.py``.

Gated (a confirmed ``fail`` fails the test; ADR-0015 two-block confirmation):

- **A** ``operon_attack_ar1(1800, 0.5)``: TU abundance AR(1) along the chromosome, no operon
  effect. Shipped configuration (expected: recorded vacuous pass) and the non-vacuous arm (100
  planted pairs at 10^-3 SD, convergent decoys only). R = ``missed_operons`` 50 / 300.
- **B** ``operon_attack_detection(1800, "prophage")``: detection and distance confounded
  (ADR-0018). Same two arms and R. M6-03's own measurement: B cannot move the FDR by more than
  ~0.005 at the FDR scale, so its pass shows level, not that the reweighting is needed.
- **M-2** (PR #148 review): Delta SD proportional to exp(-0.5 z) with z the standardised mean
  abundance of the pair, and the NO pairs' z shifted by +1 SD (NO pairs tighter than the decoys
  for a reason other than co-transcription). Non-vacuous arm only; R as A.

Gated, failing: strict xfails naming ADR-0055 (proposed, pending Khoa), the ADR-0049 pattern:

- **D-1** (PR #147 review): the discordance Delta SD proportional to exp(-0.5 z), z the
  standardised pair-mean abundance; the pooled level through the harness and the rate in the
  lowest-abundance quartile, in one test (``test_attack_d1_level``). Both fail (full tier, Slurm
  1884412: pooled 0.0521 / 0.1064 against
  0.0510 / 0.1013 at alpha 0.05 / 0.10; lowest quartile ~3x nominal). **The M6 discordance claim
  is not met** until Khoa decides ADR-0055. R = ``discordance`` 100 / 1,000. A sharded
  ``discordance`` task followed by ``--merge`` therefore exits 1 (a confirmed ``fail`` record;
  the ADR-0049 precedent of the variants task); see the test's docstring.

Recorded (numbers printed, logged and copied to docs/methods/operons.md "M6-09"; never gated, but
an escalation threshold fails the test so the breach cannot pass silently: A' and M-1 at mean FDP
> 0.15 at q 0.1, C at a rate > 2 alpha):

- **A'** ``operon_attack_ar1(1800, 0.9)``: both arms, R = 20 / 100 (module constants; the spec's
  counts, which no accepted ADR puts in tolerances.yaml). Escalation: mean FDP > 0.15 at q 0.1
  opens an ADR (lead default 1, the ADR-0021 islands precedent).
- **C** the one-detected arm on ``operon_null(1800)`` and on the detection attack: the share of
  one-detected OP pairs with P(miss) <= alpha. The step's arm (M6-05H) is not merged, so the
  oracle P(miss) (the generating detection model, the fitted within-operon Delta sample) is
  recorded; ``test_attack_c_through_the_step`` skips until the step makes ``one_detected_pmiss``
  rows. Escalation: > 2 alpha revisits M6-05 Q3. R = 50 / 300 (module constants).
- **M-1** co-imputation surplus in NO pairs (PR #148 review: Delta = 0 point masses). Recorded at
  1,800 genes (shipped and convergent-only decoys) and, full tier only, 6,000 genes; the PR #148
  oracle broke on its large genome, which the ADR on excluding imputed pairs addresses.
"""

from __future__ import annotations

import json
from functools import cache
from typing import Any

import _criteria as C
import _operon_gate as G
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels
from scipy.special import expit, logit

from giapsim.scenarios import OPERON_T

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
MO, DI = "missed_operons", "discordance"
#: A' and C replicate counts (fast, full): spec table, module constants (no ADR covers them)
R_A_PRIME = (20, 100)
R_C = (50, 300)
ESCALATE_FDR = 0.15
M_POOLED = 200
M1_LARGE_GENES = 6000
R_M1_LARGE = 50


def _r(pair: tuple[int, int], full: bool) -> int:
    return pair[1] if full else pair[0]


def _arms(name: str, *args: Any) -> list[tuple[str, Any, Any]]:
    """(label, scenario, giap cfg): the shipped configuration and the non-vacuous arm."""
    k, sd = G.NONVAC_PLANT
    return [
        ("shipped", G.scenario(name, *args), G.giap_cfg()),
        ("nonvacuous_cv_only", G.scenario(name, *args, n_planted=k, effect_sd=sd),
         G.giap_cfg(G.CV_ONLY)),
    ]  # fmt: skip


def _gate(stream: str, cfgs: Any, gcfg: Any, R: int, nonvac: bool, transform: Any = None) -> None:
    records = G.missed_replicate_fn(stream, cfgs, gcfg, transform)
    recs = G.missed_fdr_gate(MO, stream, records, R)
    G.assert_fdr_ok(recs, allow_vacuous=True, need_calls_at=0.1 if nonvac else None)
    print(f"{stream} at q 0.1: {G.fdp_summary(records, R, 0.1)}")


# --------------------------------------------------------------------------- A, A', B, M-2


@pytest.mark.calsuite(MO)
@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize("arm", [0, 1], ids=["shipped", "nonvacuous"])
def test_attack_a_ar1_rho05(arm: int, full: bool) -> None:
    label, cfgs, gcfg = _arms("operon_attack_ar1", G.N_GENES, 0.5)[arm]
    _gate(f"{MO}.attack_a.rho0.5.{label}", cfgs, gcfg, C.replicates(MO, full), arm == 1)


@pytest.mark.calsuite(MO)
@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize("arm", [0, 1], ids=["shipped", "nonvacuous"])
def test_attack_b_detection_distance(arm: int, full: bool) -> None:
    label, cfgs, gcfg = _arms("operon_attack_detection", G.N_GENES, "prophage")[arm]
    _gate(f"{MO}.attack_b.prophage.{label}", cfgs, gcfg, C.replicates(MO, full), arm == 1)


@pytest.mark.calsuite(MO)
@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize("arm", [0, 1], ids=["shipped", "nonvacuous"])
def test_recorded_attack_a_prime_rho09(arm: int, full: bool) -> None:
    label, cfgs, gcfg = _arms("operon_attack_ar1", G.N_GENES, 0.9)[arm]
    r = _r(R_A_PRIME, full)
    stream = f"{MO}.attack_a_prime.rho0.9.{label}"
    s = G.fdp_summary(G.missed_replicate_fn(stream, cfgs, gcfg), r, 0.1)
    print(f"recorded A' {label}: {json.dumps(s)}")
    _log("attack_a_prime", {"arm": label, **s})
    assert s["mean_fdp"] <= ESCALATE_FDR, f"A' mean FDP {s} > {ESCALATE_FDR}: escalate (ADR)"


def _noise_by_abundance(pairs: pd.DataFrame, classes: list[str], no_shift: float,
                        loc: float) -> pd.DataFrame:  # fmt: skip
    """Delta of the both-detected pairs of ``classes`` rescaled to ``loc + (Delta - loc) *
    exp(-0.5 z)``, z the standardised pair-mean abundance (+``no_shift`` for NO pairs)."""
    out = pairs.copy()
    a = out["abundance_log2_a"].astype(float).to_numpy()
    b = out["abundance_log2_b"].astype(float).to_numpy()
    keep = out["class"].isin(classes).to_numpy() & np.isfinite(a) & np.isfinite(b)
    m = (a + b) / 2
    z = np.zeros(len(out))
    z[keep] = (m[keep] - m[keep].mean()) / m[keep].std(ddof=1)
    z[keep & (out["class"] == "NO").to_numpy()] += no_shift
    d = b - a
    new = np.where(keep, loc + (d - loc) * np.exp(-0.5 * z), d)
    out["abundance_log2_b"] = np.where(keep, a + new, b)
    out["delta_log2"] = np.where(keep, new, out["delta_log2"].astype(float))
    return out


def _m2(sim: Any, rng: np.random.Generator) -> tuple[pd.DataFrame, None]:
    del rng
    return _noise_by_abundance(sim.pairs, ["NO", "CV", "SW"], 1.0, 0.0), None


@pytest.mark.calsuite(MO)
@pytest.mark.parametrize("full", MODES)
def test_attack_m2_abundance_mix(full: bool) -> None:
    k, sd = G.NONVAC_PLANT
    cfgs = G.scenario("missed_operon_power", G.N_GENES, k, sd)
    _gate(f"{MO}.attack_m2.no_shift1", cfgs, G.giap_cfg(G.CV_ONLY), C.replicates(MO, full), True,
          _m2)  # fmt: skip


# --------------------------------------------------------------------------- M-1


def _co_impute(sim: Any, rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame]:
    """5 % of detected genes imputed, plus both genes of 5 % of the both-detected NO pairs
    (co-imputation, rho = 1, in NO only); the fill is the median detected abundance (the default
    ``abundance.missing = median``), so a pair with both genes filled has Delta = 0 exactly."""
    g = sim.expression.genes
    tags = g["locus_tag"].astype(str).to_numpy()
    det = g["detected_true"].astype(bool).to_numpy()
    ab = g["abundance_log2_true"].astype(float).to_numpy()
    imp = det & (rng.random(len(g)) < 0.05)
    p = sim.pairs
    pos = {t: i for i, t in enumerate(tags)}
    ia = np.array([pos[str(t)] for t in p["gene_a"]])
    ib = np.array([pos[str(t)] for t in p["gene_b"]])
    no = (p["class"] == "NO").to_numpy() & det[ia] & det[ib]
    pick = no & (rng.random(len(p)) < 0.05)
    imp[ia[pick]] = imp[ib[pick]] = True
    obs = np.where(det, np.where(imp, float(np.median(ab[det])), ab), np.nan)
    out = p.copy()
    out["abundance_log2_a"], out["abundance_log2_b"] = obs[ia], obs[ib]
    out["delta_log2"] = obs[ib] - obs[ia]
    genes = pd.DataFrame({"locus_tag": tags, "abundance_imputed": imp})
    return out, genes


@pytest.mark.calsuite(MO)
@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize("decoys", ["shipped", "cv_only"])
def test_recorded_attack_m1_co_imputation(decoys: str, full: bool) -> None:
    gcfg = G.giap_cfg() if decoys == "shipped" else G.giap_cfg(G.CV_ONLY)
    cfgs = G.scenario("operon_null", G.N_GENES)
    stream = f"{MO}.attack_m1.co_impute.{decoys}"
    r = C.replicates(MO, full)
    s = G.fdp_summary(G.missed_replicate_fn(stream, cfgs, gcfg, _co_impute), r, 0.1)
    print(f"recorded M-1 {decoys}: {json.dumps(s)}")
    _log("attack_m1", {"decoys": decoys, **s})
    assert s["mean_fdp"] <= ESCALATE_FDR, f"M-1 {s} > {ESCALATE_FDR}: escalate (imputed-pair ADR)"


@pytest.mark.calsuite(MO)
@pytest.mark.calibration_full
def test_recorded_attack_m1_large_genome() -> None:
    """M-1 where the PR #148 oracle broke (~2,000 NO pairs): 6,000 genes, shipped decoys,
    R = ``R_M1_LARGE`` (module constant; full tier only, ~1 s per replicate)."""
    cfgs = G.scenario("operon_null", M1_LARGE_GENES)
    stream = f"{MO}.attack_m1.co_impute.large"
    s = G.fdp_summary(G.missed_replicate_fn(stream, cfgs, G.giap_cfg(), _co_impute),
                      R_M1_LARGE, 0.1)  # fmt: skip
    print(f"recorded M-1 large genome: {json.dumps(s)}")
    _log("attack_m1", {"decoys": "shipped", "n_genes": M1_LARGE_GENES, **s})
    assert s["mean_fdp"] <= ESCALATE_FDR, f"M-1 large {s} > {ESCALATE_FDR}: escalate (ADR)"


# --------------------------------------------------------------------------- D-1


def _d1(seed: int) -> G.DiscRun:
    """One D-1 genome. ``abundance`` is the pre-transform pair mean, the z the noise was planted
    on (the transform changes ``abundance_log2_b``), so the quartile is cut on the planted z."""
    sim = G.simulate(G.scenario("operon_null", G.N_GENES), seed)
    noisy = _noise_by_abundance(sim.pairs, ["OP"], 0.0, OPERON_T[0])
    return G.discordance_run(sim, pairs=noisy, abundance_pairs=sim.pairs)


_d1_cached = cache(_d1)
D1_STREAM = f"{DI}.attack_d1.noise_by_abundance"


#: ADR-0055 (proposed, pending Khoa): D-1 fails the discordance level; strict xfail, as the M7
#: gate's ADR-0049 xfails (test_variants_gate.py). Seeds are fixed, so the failure is
#: deterministic and an XPASS means the method or the generator changed. No bound is widened (R3).
ADR_0055 = "ADR-0055 pending: D-1 abundance-dependent noise; M6 discordance claim not met"
D1_XFAIL = pytest.mark.xfail(strict=True, raises=AssertionError, reason=ADR_0055)


@D1_XFAIL
@pytest.mark.calsuite(DI)
@pytest.mark.parametrize("full", MODES)
def test_attack_d1_level(full: bool) -> None:
    """Gated D-1, one test so both checks read the same genomes (one simulation per replicate):

    - the pooled level criterion through the harness (suite log ``discordance.jsonl``, ADR-0015):
      a confirmed ``fail`` at alpha 0.10 in the fast tier (0.110 vs 0.104) and at 0.05 and 0.10
      in the full tier (0.0521 / 0.1064 vs 0.0510 / 0.1013);
    - the lowest-abundance quartile (cut on the planted z): the rate over every tested pair of
      the quartile against ``level_bound(alpha, n)`` (n pairs treated as independent, the
      optimistic bound); about 3x nominal (full tier 0.0273 / 0.157 / 0.277).

    Both are recorded before the assertion. **Sharded runs:** a shard has no verdict and skips
    (after logging), but the records are in the suite log, so ``_criteria.py --shard ... --merge``
    of the ``discordance`` task merges them into a confirmed ``fail`` and exits 1 (MERGE_FAIL;
    ``_merge_main`` knows nothing of xfails). That is intended until ADR-0055 is decided: a
    failing attack blocks the gate, as the variants task does under ADR-0049. Do not move this
    log out of the suite to make the task green.
    """
    lv = G.levels()
    r = C.replicates(DI, full)

    def rep(rng: np.random.Generator) -> np.ndarray:
        run = _d1_cached(int(rng.bit_generator.seed_seq.entropy))
        p = rng.choice(run.p, M_POOLED, replace=False)
        return np.array([float((p <= a).mean()) for a in lv])

    gate = calibrate_levels(DI, "level", lv, rep, r, stream=D1_STREAM, m=M_POOLED)
    print("D-1 pooled:\n" + "\n".join(gate.detail()))
    G.skip_without_verdict(gate.records)
    hits, n = np.zeros(len(lv)), 0
    for i in range(r):
        run = _d1_cached(C.suite_seed(D1_STREAM, i))
        low = run.abundance <= np.quantile(run.abundance, 0.25)
        hits += [(run.p[low] <= a).sum() for a in lv]
        n += int(low.sum())
    rate = hits / n
    bounds = [C.level_bound(a, n) for a in lv]
    print(f"D-1 lowest quartile: R = {r}, n = {n}, rate {rate.round(4).tolist()}"
          f" (bounds {np.round(bounds, 4).tolist()})")  # fmt: skip
    _log("attack_d1_lowest_quartile", {"R": r, "n": n, "rate": rate.tolist(), "bounds": bounds})
    failed = gate.failed
    assert not failed and np.all(rate <= bounds), (gate.detail(), rate.tolist(), ADR_0055)


# --------------------------------------------------------------------------- C


def _oracle_pmiss(sim: Any, cfgs: Any) -> np.ndarray:
    """P(miss) of every eligible one-detected OP pair under the generating detection model:
    mean over the fitted Delta sample r of 1 - expit(L_u + b (a_v +- r)), with
    L_u = logit(p_detect_u) - b a_u (b = the abundance slope, 0 under covariates_only)."""
    from giap.operons.discordance import fit_rows

    ecfg = cfgs[1]
    b = float(ecfg.beta["log2_abundance"]) if ecfg.detection == "abundance" else 0.0
    g = sim.expression.genes
    ab = dict(zip(g["locus_tag"].astype(str), g["abundance_log2_true"].astype(float),
                  strict=True))  # fmt: skip
    p = sim.pairs
    ratios = fit_rows(p)
    op = p[p["class"] == "OP"]
    da, db = op["detected_a"].astype(bool).to_numpy(), op["detected_b"].astype(bool).to_numpy()
    out = []
    for row, a_det, b_det in zip(op.itertuples(index=False), da, db, strict=True):
        if a_det == b_det:
            continue
        if a_det:  # u = b, downstream: a_u = a_v + Delta
            a_v, u, pu, r = float(row.abundance_log2_a), str(row.gene_b), row.p_detect_b, ratios
        else:  # u = a, upstream: a_u = a_v - Delta
            a_v, u, pu, r = float(row.abundance_log2_b), str(row.gene_a), row.p_detect_a, -ratios
        lu = float(logit(np.clip(pu, 1e-12, 1 - 1e-12))) - b * ab[u]
        out.append(float(np.mean(1.0 - expit(lu + b * (a_v + r)))))
    return np.asarray(out)


@pytest.mark.calsuite(DI)
@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize("scen", ["operon_null", "operon_attack_detection"])
def test_recorded_attack_c_one_detected_oracle(scen: str, full: bool) -> None:
    lv = G.levels()
    args = (G.N_GENES,) if scen == "operon_null" else (G.N_GENES, "prophage")
    cfgs = G.scenario(scen, *args)
    r = _r(R_C, full)
    stream = f"{DI}.attack_c.oracle.{scen}"
    pm = np.concatenate([_oracle_pmiss(G.simulate(cfgs, C.suite_seed(stream, i)), cfgs)
                         for i in range(r)])  # fmt: skip
    rate = np.array([float((pm <= a).mean()) for a in lv])
    print(f"recorded C (oracle P(miss)) {scen}: R = {r}, one-detected pairs {pm.size},"
          f" rate {rate.round(4).tolist()} at {list(lv)}")  # fmt: skip
    _log("attack_c_oracle", {"scenario": scen, "R": r, "n": int(pm.size), "rate": rate.tolist()})
    assert pm.size > 0 and np.all((pm >= 0) & (pm <= 1))
    assert np.all(rate <= 2 * np.asarray(lv)), f"C above 2 alpha: revisit M6-05 Q3 ({rate})"


@pytest.mark.calsuite(DI)
@pytest.mark.calibration
def test_attack_c_through_the_step() -> None:
    """Probe of the step's one-detected arm (M6-05H): skips until the step makes
    ``one_detected_pmiss`` rows; then the oracle arm above moves here (statistician)."""
    cfgs = G.scenario("operon_null", G.N_GENES)
    sim = G.simulate(cfgs, C.suite_seed(f"{DI}.attack_c.step", 0))
    run = G.discordance_run(sim, p_miss=lambda u, a_v, ratios, v: 0.5)
    if run.pmiss.size == 0:
        n = run.n_one_detected
        pytest.skip(f"M6-05H not merged: no one_detected_pmiss rows ({n} one-detected OP pairs);"
                    " attack C is the oracle arm")  # fmt: skip
    print(f"C through the step (probe model P(miss) = 0.5): {run.pmiss.size} rows")
    assert np.all((run.pmiss >= 0) & (run.pmiss <= 1))


def _log(name: str, row: dict[str, Any]) -> None:
    path = C.log_dir() / f"operons_{name}.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(row) + "\n")
