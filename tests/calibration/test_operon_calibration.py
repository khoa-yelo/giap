"""The M6 gate, simulation half (issue M6-09; GIAP_SPEC.md §7.4.7 "Calibration" criteria 1 and 3,
§8.1 rows "Operon discordance" and "Missed operon", §8.4 suites ``missed_operons`` and
``discordance``; docs/specs/M6-09-calibration-gate.md with its lead defaults and Khoa's decision
of 2026-10-22). Criterion 2 (the attack null) is attacks/test_attack_operons.py.

R = ``replicates_fast/full.<suite>`` of the protected tolerances.yaml (missed_operons 50 / 300,
discordance 100 / 1,000); every check goes through ``_criteria`` (ADR-0015 two-block confirmation,
sharding, jsonl log). Helpers and the definition of a true call: ``_operon_gate``.

``missed_operons`` suite (§8.1: one decoy-FDR family per run, q <= 0.1, checked at every level):

- ``test_missed_operon_fdr_gate``: the spec generator ``missed_operon_power(1800, 50, 0.1)`` with
  the shipped configuration (CV + SW decoys, pi0 = 1). The method makes no call there (the
  identifiability bound, docs/methods/operons.md "Power"), so the expected verdict is a
  **recorded vacuous pass** (Khoa, 2026-10-22); any call is judged by the FDR criterion.
- ``test_missed_operon_fdr_nonvacuous``: the same check where calls exist: 100 planted pairs at
  10^-3 SD and the convergent decoy set alone (``operons.decoys = [convergent]``, a supported
  configuration). Must pass and must make calls at q = 0.1 (the §8.1 default; at 0.01 and 0.05
  the bound still leaves no call, recorded vacuous): the arm in which the gate can fail.
- ``test_teeth_pi0_misconfigured_fails``: the non-vacuous arm with ``operons.pi0 = 0.1`` (q scaled
  by 0.1 with no estimate) must be a confirmed ``fail``.
- ``test_missed_operon_power_grid`` (full tier only, R = 50 per cell): recorded table.

``discordance`` suite (§8.1: one BH family per run over the ``both_detected_t`` rows, q <= 0.1):

- ``test_discordance_pooled_level``: ``operon_null(1800)``, ``M_POOLED`` tested pairs drawn at
  random per replicate, pooled level at every alpha (the p of one run share the fitted t, so the
  pooled bound is the optimistic one; the per-replicate KS below is the independent check);
- ``test_discordance_p_uniform``: one tested pair at random per replicate, KS against U(0, 1);
- ``test_discordance_bh_family_fdr``: complete null, so FDP = 1{any BH call at q};
- ``test_teeth_normal_tail_fails``: the p computed from a normal with the fitted location and
  scale (the t tail ignored), pooled as above, must be a confirmed ``fail`` at alpha 0.01.
"""

from __future__ import annotations

import json
import time
from functools import cache

import _criteria as C
import _operon_gate as G
import numpy as np
import pytest
from _gate import calibrate_levels
from scipy import stats

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
MO, DI = "missed_operons", "discordance"
#: tested pairs pooled per replicate (~460 are tested on operon_null(1800))
M_POOLED = 200
#: the power grid (spec table: |Delta| in SD x planted count, R = 50 per cell, recorded)
POWER_EFFECTS = (0.25, 0.5, 1.0)
POWER_PLANTED = (20, 50, 100)
POWER_R = 50


# --------------------------------------------------------------------------- missed operons


@pytest.mark.calsuite(MO)
@pytest.mark.parametrize("full", MODES)
def test_missed_operon_fdr_gate(full: bool) -> None:
    stream = f"{MO}.gate.power_1800_50_0.1"
    cfgs = G.scenario("missed_operon_power", G.N_GENES, 50, 0.1)
    recs = G.missed_fdr_gate(
        MO, stream, G.missed_replicate_fn(stream, cfgs, G.giap_cfg()), C.replicates(MO, full)
    )
    G.assert_fdr_ok(recs, allow_vacuous=True)


@pytest.mark.calsuite(MO)
@pytest.mark.parametrize("full", MODES)
def test_missed_operon_fdr_nonvacuous(full: bool) -> None:
    k, sd = G.NONVAC_PLANT
    stream = f"{MO}.gate.nonvacuous_cv_only"
    cfgs = G.scenario("missed_operon_power", G.N_GENES, k, sd)
    records = G.missed_replicate_fn(stream, cfgs, G.giap_cfg(G.CV_ONLY))
    r = C.replicates(MO, full)
    recs = G.missed_fdr_gate(MO, stream, records, r)
    G.assert_fdr_ok(recs, allow_vacuous=True, need_calls_at=0.1)
    print(f"{stream} at q 0.1: {G.fdp_summary(records, r, 0.1)}")


@pytest.mark.calsuite(MO)
@pytest.mark.calibration
def test_teeth_pi0_misconfigured_fails(teeth_log: object) -> None:
    k, sd = G.NONVAC_PLANT
    stream = f"{MO}.teeth.pi0_0.1"
    cfgs = G.scenario("missed_operon_power", G.N_GENES, k, sd)
    records = G.missed_replicate_fn(stream, cfgs, G.giap_cfg(G.CV_ONLY, pi0=0.1))
    rec = C.calibrate_fdr_counts(
        MO, stream, "q0.1", 0.1, records, C.replicates(MO, False), log_dir=teeth_log
    )  # type: ignore[arg-type]
    print(f"teeth pi0 = 0.1: {rec['status']} meanFDP={rec['observed']} bound={rec['bound']}")
    assert rec["status"] == "fail", rec


@pytest.mark.calsuite(MO)
@pytest.mark.calibration_full
def test_missed_operon_power_grid() -> None:
    """Recorded, never gated: the family is near-powerless by construction (Khoa, 2026-10-22)."""
    t0 = time.perf_counter()
    rows = []
    for k in POWER_PLANTED:
        for e in POWER_EFFECTS:
            stream = f"{MO}.power.k{k}.sd{e:g}"
            cfgs = G.scenario("missed_operon_power", G.N_GENES, k, e)
            records = G.missed_replicate_fn(stream, cfgs, G.giap_cfg())
            rows.append({"n_planted": k, "effect_sd": e, "stream": stream,
                         **G.fdp_summary(records, POWER_R, 0.1)})  # fmt: skip
    for row in rows:
        print(json.dumps(row))
    path = C.log_dir() / "missed_operons_power.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps({"rows": rows, "seconds": time.perf_counter() - t0}) + "\n")
    assert len(rows) == 9 and all(np.isfinite(r["median_min_q"]) for r in rows)


# --------------------------------------------------------------------------- discordance

DISC_STREAM = f"{DI}.gate.operon_null_1800"


@cache
def _null_run(seed: int) -> G.DiscRun:
    return G.discordance_run(G.simulate(G.scenario("operon_null", G.N_GENES), seed))


def _seed(rng: np.random.Generator) -> int:
    return int(rng.bit_generator.seed_seq.entropy)  # = suite_seed(stream, i)


@pytest.mark.calsuite(DI)
@pytest.mark.parametrize("full", MODES)
def test_discordance_pooled_level(full: bool) -> None:
    lv = G.levels()

    def rep(rng: np.random.Generator) -> np.ndarray:
        run = _null_run(_seed(rng))
        assert run.p.size >= M_POOLED, run.p.size
        p = rng.choice(run.p, M_POOLED, replace=False)
        return np.array([float((p <= a).mean()) for a in lv])

    gate = calibrate_levels(DI, "level", lv, rep, C.replicates(DI, full), stream=DISC_STREAM,
                            m=M_POOLED)  # fmt: skip
    print("\n".join(gate.detail()))
    gate.assert_ok()


@pytest.mark.calsuite(DI)
@pytest.mark.parametrize("full", MODES)
def test_discordance_p_uniform(full: bool) -> None:
    def rep(rng: np.random.Generator) -> float:
        run = _null_run(_seed(rng))
        return float(rng.choice(run.p))

    rec = C.calibrate(DI, "uniform", 0.0, rep, C.replicates(DI, full), stream=DISC_STREAM)
    print(f"KS p = {rec['observed']} (R = {rec['R']}, status {rec['status']})")
    G.skip_without_verdict([rec])
    assert rec["status"] in {"pass", "unconfirmed"}, rec


@pytest.mark.calsuite(DI)
@pytest.mark.parametrize("full", MODES)
def test_discordance_bh_family_fdr(full: bool) -> None:
    lv = G.levels()

    def rep(rng: np.random.Generator) -> np.ndarray:
        q = _null_run(_seed(rng)).q
        return np.array([float((q <= a).any()) for a in lv])

    gate = calibrate_levels(DI, "fdr", lv, rep, C.replicates(DI, full), stream=DISC_STREAM)
    print("\n".join(gate.detail()))
    gate.assert_ok()


@pytest.mark.calsuite(DI)
@pytest.mark.calibration
def test_teeth_normal_tail_fails(teeth_log: object) -> None:
    """Reuses the gate's genomes (same stream, so the same seeds and, on one worker, the cached
    runs) and logs to the private teeth directory."""

    def rep(rng: np.random.Generator) -> float:
        run = _null_run(_seed(rng))
        x = rng.choice(run.ratios, M_POOLED, replace=False)
        p = 2 * stats.norm.sf(np.abs(x - run.fit.loc) / run.fit.scale)
        return float((p <= 0.01).mean())

    rec = C.calibrate(DI, "level", 0.01, rep, C.replicates(DI, False), m=M_POOLED,
                      stream=DISC_STREAM, log_dir=teeth_log)  # type: ignore[arg-type]  # fmt: skip
    print(f"teeth normal tail: {rec['status']} rate={rec['observed']} bound={rec['bound']}")
    assert rec["status"] == "fail", rec
