"""M6-13 attacks on the abundance-dependent discordance scale (ADR-0055 option (a);
docs/specs/M6-13-discordance-abundance-sigma.md "Calibration and the gate", the statistician's
measurement and decisions S-8, S-9, S-11; GIAP_SPEC.md §8.4 suite ``discordance``, §8.5).
Helpers: ``../_operon_sigma.py`` (genomes: the D-1 stream's, as the statistician measured them).

Gated (a confirmed ``fail`` fails the test; ADR-0015 two-block confirmation; R =
``discordance`` 100 / 1,000):

- **D-1, log-linear** (``test_attack_d1_loglinear_quartiles``): the lowest planted-z quartile
  against ``level_bound(alpha, n)`` (as ``test_attack_d1_level``) and every quartile recorded with
  the 2 alpha escalation (S-9). The pooled D-1 level is ``test_attack_d1_level`` itself once its
  strict xfail is removed (lead edit, same PR as the implementation).
- **teeth** (``test_teeth_d1_constant_fails``, fast tier): D-1 with ``sigma_model = "constant"``
  must still be a confirmed ``fail`` pooled (same genomes and subsample as
  ``test_attack_d1_level``; private log), and its lowest quartile above 2 alpha.
- **D-1b-floor** and **D-1c**: pooled level gated; every quartile recorded + 2 alpha escalation
  (D-1b on the planted z; D-1c on the upstream gene's true abundance).
- **D-1b-step**: the same criteria; measured to fail the pooled level at full R (0.0518 / 0.1024
  against 0.0510 / 0.1013; Q2 1.8x), so a **strict xfail** naming the pending ADR-0055 amendment
  in the full tier. The fast tier cannot see a 1.02x pooled excess (the oracle's first 100
  genomes: 0.0521 against 0.0531; Q2 1.82x < 2x), so there the xfail is not strict.

Recorded (with the 2 alpha escalation): the null quartile check on ``operon_null(1800)``, both
models, cut on the pair mean m and on the upstream abundance a (S-9: the m cut is not a valid
level check, m contains Delta / 2), and the gamma-hat records of every arm.
"""

from __future__ import annotations

from typing import Any

import _criteria as C
import _operon_gate as G
import _operon_sigma as S
import numpy as np
import pytest
from _gate import calibrate_levels

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
DI = S.DI
ADR_0055_AMENDMENT = (
    "ADR-0055 amendment pending: a step σ(abundance) is outside the smooth-monotone claim"  # noqa: RUF001
)
STEP_MODES = [
    pytest.param(False, marks=[pytest.mark.calibration,
                               pytest.mark.xfail(strict=False, raises=AssertionError,
                                                 reason=ADR_0055_AMENDMENT)], id="fast"),
    pytest.param(True, marks=[pytest.mark.calibration_full,
                              pytest.mark.xfail(strict=True, raises=AssertionError,
                                                reason=ADR_0055_AMENDMENT)], id="full"),
]  # fmt: skip


def _shape_runs(shape: str, model: str, R: int) -> list[S.SigmaRun]:
    return [S.shape_run(C.suite_seed(S.D1_STREAM, i), shape, model) for i in range(R)]


def _report(name: str, rec: dict[str, Any], lv: tuple[float, ...]) -> list[str]:
    bad = S.escalation(rec, lv)
    print(f"{name} quartiles on {rec['stratum']}: R = {rec['R']}, n = {rec['n']},"
          f" rate {np.round(rec['rate'], 4).tolist()}")  # fmt: skip
    S.log(name, rec)
    return bad


def _pooled(stream: str, get_run: Any, R: int) -> Any:
    lv = G.levels()
    gate = calibrate_levels(DI, "level", lv, S.pooled_rep(get_run, lv), R, stream=stream,
                            m=S.M_POOLED)  # fmt: skip
    print(f"{stream} pooled:\n" + "\n".join(gate.detail()))
    return gate


# --------------------------------------------------------------------------- D-1 and its teeth


@pytest.mark.calsuite(DI)
@pytest.mark.parametrize("full", MODES)
def test_attack_d1_loglinear_quartiles(full: bool) -> None:
    """Measured (R = 1,000): Q1 .0093 / .0478 / .0969 against .0106 / .0512 / .1017; worst
    quartile Q4 1.05x; gamma-hat -0.513 (planted -0.5)."""
    lv = G.levels()
    r = C.replicates(DI, full)
    runs = _shape_runs("d1", "loglinear", r)
    rec = S.quartile_rates(runs, "m", lv)
    bad = _report("d1_loglinear", rec, lv)
    g = S.gamma_record(runs)
    print(f"D-1 gamma-hat: {g}")
    S.log("d1_loglinear_gamma", g)
    low_ok = [rec["rate"][0][j] <= rec["bounds"][0][j] for j in range(len(lv))]
    assert all(low_ok), ("lowest quartile above level_bound", rec["rate"][0], rec["bounds"][0])
    assert not bad, bad
    assert -0.6 <= g["gamma_mean"] <= -0.4, g  # recovers the planted -0.5 through the step


@pytest.mark.calsuite(DI)
@pytest.mark.calibration
def test_teeth_d1_constant_fails(teeth_log: object) -> None:
    """The constant model under D-1 (merged step, Slurm 1884412: pooled 0.110 vs 0.104 at
    alpha 0.10 in the fast tier; lowest quartile ~3x) must still fail: the attack bites."""
    lv = G.levels()
    r = C.replicates(DI, False)

    def get(rng: np.random.Generator) -> S.SigmaRun:
        return S.shape_run(S.d1_seed(S.D1_STREAM, r, rng), "d1", "constant")

    gate = calibrate_levels(DI, "level", lv, S.pooled_rep(get, lv), r, stream=S.D1_STREAM,
                            m=S.M_POOLED, log_dir=teeth_log)  # type: ignore[arg-type]  # fmt: skip
    print("teeth D-1 constant:\n" + "\n".join(gate.detail()))
    assert any(rec["status"] == "fail" for rec in gate.records), gate.detail()
    rec = S.quartile_rates(_shape_runs("d1", "constant", r), "m", lv)
    print(f"teeth D-1 constant lowest quartile {np.round(rec['rate'][0], 4).tolist()}")
    assert all(rec["rate"][0][j] > S.ESCALATE_QUARTILE * a for j, a in enumerate(lv)), rec


# --------------------------------------------------------------------------- D-1b, D-1c


def _shape_gate(shape: str, full: bool) -> None:
    lv = G.levels()
    r = C.replicates(DI, full)
    stream = f"{DI}.attack_d1b.{shape}.loglinear"

    def get(rng: np.random.Generator) -> S.SigmaRun:
        return S.shape_run(S.d1_seed(stream, r, rng), shape, "loglinear")

    gate = _pooled(stream, get, r)
    G.skip_without_verdict(gate.records)
    runs = _shape_runs(shape, "loglinear", r)
    bad = _report(f"d1b_{shape}", S.quartile_rates(runs, "m", lv), lv)
    S.log(f"d1b_{shape}_gamma", S.gamma_record(runs))
    assert not gate.failed and not bad, (gate.detail(), bad)


@pytest.mark.calsuite(DI)
@pytest.mark.parametrize("full", MODES)
def test_attack_d1b_floor(full: bool) -> None:
    """Measured (R = 1,000): pooled .0100 / .0506 / .1004; worst quartile Q4 1.22x."""
    _shape_gate("floor", full)


@pytest.mark.calsuite(DI)
@pytest.mark.parametrize("full", STEP_MODES)
def test_attack_d1b_step(full: bool) -> None:
    """Pre-registered failure (statistician, Slurm 1885387): pooled .0095 / .0518 / .1024 against
    .0104 / .0510 / .1013; Q2 .0181 / .0910 / .1662 (1.8x). No smooth sigma(z) removes a jump;
    the claim is narrowed by the ADR-0055 amendment, never by widening a bound (R3). If the full
    tier XPASSes, the statistician re-reads the measurement before the xfail is removed."""
    _shape_gate("step", full)


@pytest.mark.calsuite(DI)
@pytest.mark.parametrize("full", MODES)
def test_attack_d1c_gene_noise(full: bool) -> None:
    """Per-gene measurement noise (errors in variables; S-11). Measured (R = 1,000): pooled
    .0098 / .0480 / .0973; worst quartile on the true upstream abundance Q4 1.68x at 0.01."""
    lv = G.levels()
    r = C.replicates(DI, full)
    stream = f"{DI}.attack_d1c.gene_noise.loglinear"

    def get(rng: np.random.Generator) -> S.SigmaRun:
        return S.d1c_run(S.d1_seed(stream, r, rng), "loglinear")

    gate = _pooled(stream, get, r)
    G.skip_without_verdict(gate.records)
    runs = [S.d1c_run(C.suite_seed(S.D1_STREAM, i), "loglinear") for i in range(r)]
    bad = _report("d1c", S.quartile_rates(runs, "a_true", lv), lv)
    S.log("d1c_gamma", S.gamma_record(runs))
    assert not gate.failed and not bad, (gate.detail(), bad)


# --------------------------------------------------------------------------- the null, recorded


@pytest.mark.calsuite(DI)
@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize("model", ["loglinear", "constant"])
def test_recorded_null_quartiles_on_m_and_a(model: str, full: bool) -> None:
    """``operon_null(1800)`` on the gate genomes, quartiles on the pair mean m and on the upstream
    abundance a. Recorded (S-9): on m even the constant method exceeds the independent-pairs bound
    in the outer quartiles (Q1 .0117 vs .0106 at 0.01) because m = a + Delta / 2; on a every
    quartile is at level. Escalation only at 2 alpha. The gamma-hat record is the S-10 claim:
    detection truncation is not measurable in the simulator (null gamma-hat 0.001 +- 0.041)."""
    lv = G.levels()
    r = C.replicates(DI, full)
    runs = [S.shape_run(C.suite_seed(S.DISC_STREAM, i), "none", model) for i in range(r)]
    bad = []
    for stratum in ("m", "a"):
        bad += _report(f"null_{model}", S.quartile_rates(runs, stratum, lv), lv)
    pooled = np.array([[float((x.p <= a).mean()) for a in lv] for x in runs]).mean(axis=0)
    print(f"null {model}: mean per-genome rate {np.round(pooled, 4).tolist()}")
    assert not bad, bad
    if model == "loglinear":
        g = S.gamma_record(runs)
        print(f"null gamma-hat: {g}")
        S.log("null_gamma", g)
        assert abs(g["gamma_mean"]) <= 0.02 + 3.0 * g["gamma_sd"] / np.sqrt(r), g
        assert g["models"].get("loglinear:None", 0) >= 0.95 * r, g  # fallbacks are rare (S-6)
