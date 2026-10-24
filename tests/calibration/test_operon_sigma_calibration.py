"""M6-13 power of the abundance-dependent discordance scale (docs/specs/M6-13-discordance-abundance-
sigma.md "Calibration and the gate" rows "power non-inferiority" and "power, BH recall"; the
statistician's decision S-8; GIAP_SPEC.md §8.4 suite ``discordance``). The level arms are in
attacks/test_attack_operons_sigma.py.

- ``test_power_non_inferior`` (gated): ``discordance_power(1800, 40, 2.0)``, both models on the
  same genomes. Per replicate d = (log-linear - constant) rejection rate of the planted pairs at
  alpha = 0.05; pass when mean(d) - 2 SE(d) >= -0.02 (an absolute margin, about 3 % of the
  0.647 power). Measured (R = 1,000, Slurm 1885388): d = -0.0060 (SE 0.0011). R =
  ``discordance`` 100 / 1,000. Not a level check, so not through ``calibrate``; the summary goes to
  ``operons_sigma_power.jsonl``.
- recorded: d at alpha 0.01 and 0.10, the BH recall at q 0.1 of both models (~0 under both;
  methods note, D-3), and effect 4.0 (``test_recorded_power_effect4``).
"""

from __future__ import annotations

from typing import Any

import _criteria as C
import _operon_gate as G
import _operon_sigma as S
import numpy as np
import pytest

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
DI = S.DI
N_PLANTED = 40
#: S-8: the non-inferiority margin (absolute, planted-pair rejection rate at alpha = 0.05)
NI_MARGIN = 0.02
NI_ALPHA = 0.05


def _power_rep(i: int, effect: float) -> dict[str, Any]:
    bs = S.base(C.suite_seed(f"{DI}.power.sigma.{effect:g}", i),
                ("discordance_power", G.N_GENES, N_PLANTED, effect))  # fmt: skip
    out: dict[str, Any] = {}
    for model in ("loglinear", "constant"):
        run = S.sigma_run(bs.pairs, bs.genes, model, {})
        planted = np.isin(run.pair_id, sorted(bs.truth_discordant))
        assert planted.sum() > 0, "no planted pair tested"
        out[model] = {
            "rej": [float((run.p[planted] <= a).mean()) for a in G.levels()],
            "bh_recall": float((run.q[planted] <= 0.1).mean()),
            "bh_calls": int((run.q <= 0.1).sum()),
            "k": int(planted.sum()),
        }
    return out


def _summary(effect: float, R: int) -> dict[str, Any]:
    reps = [_power_rep(i, effect) for i in range(R)]
    lv = G.levels()
    d = np.array([[r["loglinear"]["rej"][j] - r["constant"]["rej"][j] for j in range(len(lv))]
                  for r in reps])  # fmt: skip
    se = d.std(axis=0, ddof=1) / np.sqrt(R)
    return {
        "effect": effect, "R": R, "levels": list(lv),
        "d_mean": d.mean(axis=0).tolist(), "d_se": se.tolist(),
        "power": {m: np.mean([r[m]["rej"] for r in reps], axis=0).tolist()
                  for m in ("loglinear", "constant")},
        "bh_recall_q0.1": {m: float(np.mean([r[m]["bh_recall"] for r in reps]))
                           for m in ("loglinear", "constant")},
    }  # fmt: skip


@pytest.mark.calsuite(DI)
@pytest.mark.parametrize("full", MODES)
def test_power_non_inferior(full: bool) -> None:
    r = C.replicates(DI, full)
    s = _summary(2.0, r)
    print(f"power (effect 2.0): {s}")
    S.log("power", s)
    j = list(G.levels()).index(NI_ALPHA)
    lower = s["d_mean"][j] - 2.0 * s["d_se"][j]
    assert lower >= -NI_MARGIN, (f"log-linear power below constant by more than {NI_MARGIN}: "
                                 f"mean(d) - 2 SE = {lower:.4f}", s)  # fmt: skip


@pytest.mark.calsuite(DI)
@pytest.mark.parametrize("full", MODES)
def test_recorded_power_effect4(full: bool) -> None:
    """Recorded: measured d -0.0008 at alpha 0.05 and +0.0138 at 0.01 (R = 1,000)."""
    r = C.replicates(DI, full)
    s = _summary(4.0, r)
    print(f"power (effect 4.0): {s}")
    S.log("power", s)
    assert all(np.isfinite(s["d_mean"])) and all(0 <= x <= 1 for x in s["power"]["loglinear"])
