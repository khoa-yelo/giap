"""Attack nulls for the M8-05 TU-collapsed GSEA (GIPPA_SPEC.md §7.7.4, §8.2, §8.5; statistician
review of M8-05, 2026-10-21). Statistician-owned.

Each arm keeps "abundance independent of the set" at the gene level and adds one structure the
committed null arms (tests/calibration/test_gsea_tu_calibration.py) do not have. Probe evidence
(R = 2,000 single-set replicates through an independent implementation of spec rules 1-4 on
M8-02's ``fgsea_simple``; ``$GIPPA_ROOT/stat/m805/``), rejection rate at alpha 0.01 / 0.05 / 0.10:

* **(a) default median imputation** (gated): 20 % of genes imputed to the observed median, a
  contiguous 15-TU set: 0.0125 / 0.0475 / 0.090. The committed ties arm covers only ``min``.
* **(b) size-biased membership** (gated): a set of 40 genes drawn uniformly over *genes* in an
  E. coli-like operon-length world (geometric, ~half singletons) with weak within-TU correlation
  (rho 0.2), so the any-gene rule favours long TUs whose means have small variance:
  0.003 / 0.0225 / 0.054 (conservative: a power cost, not a level failure; M8-15 measures it).
* **(c) between-TU AR(1), scattered set** (gated): phi 0.5, 15 TUs drawn uniformly:
  0.0105 / 0.0555 / 0.109.
* **(d) between-TU AR(1), genomically clustered set** (strict xfail, ADR owed to Khoa): phi 0.5,
  one run of 15 adjacent TUs: 0.075 / 0.186 / 0.271 (5 runs of 3: 0.040 / 0.130 / 0.211;
  phi 0.2: 0.024 / 0.085 / 0.139). The TU-label permutation of fgsea treats TUs as exchangeable,
  which positional autocorrelation breaks for clustered sets (M8-02 review item (a)).
* **(e) detectability confounded with membership** (strict xfail, ADR owed to Khoa): detection
  ~ logistic(2 (a - 1.5 [in set] + 0.5)), only detected genes ranked (§7.7.4): the detected
  members are selected for high abundance: 0.515 / 0.786 / 0.868 (shift 1.0, slope 1:
  0.047 / 0.137 / 0.222). §8.2's p_detect adjustment covers detection tests only.

Not gated (borderline, recorded for M8-15): median imputation at 50 % inside the set vs 10 %
outside, 0.020 / 0.078 / 0.127; ``min`` imputation with the same confounding, 0.99 / 1.0 / 1.0
(by design of ``min``: unquantified counts as low; §8.6 reruns are where it is varied).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels
from helpers.gsea_data import Gene, frame, gene_sets

from gippa.config.model import GseaConfig
from gippa.rng import for_step

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("gsea")]

SUITE = "gsea"
CFG = GseaConfig()
SIZE = (10, 500)
ALPHA = 0.05
XFAIL_ADR = pytest.mark.xfail(
    strict=True, reason="ADR-0025/ADR-0026: no all-rows level claim under this structure"
)


def _p(genes: pd.DataFrame, members: list[str], rng: np.random.Generator) -> float:
    from gippa.function import gsea

    seed = int(rng.integers(0, 2**31 - 1))
    out = gsea.gsea_table(genes, [gene_sets("kegg", {"S": members})], SIZE, CFG, ALPHA,
                          for_step("function.gsea", seed))  # fmt: skip
    assert len(out) == 1, out
    return float(out["p_value"].iloc[0])


def _lengths(rng: np.random.Generator, n: int, kind: str) -> np.ndarray:
    if kind == "1-4":
        return rng.integers(1, 5, n)
    return np.minimum(rng.geometric(0.48, n), 14)  # E. coli-like: ~half singletons


def _world(rng: np.random.Generator, n: int, kind: str, rho: float, phi: float = 0.0
           ) -> tuple[np.ndarray, np.ndarray]:  # fmt: skip
    """Gene -> TU index and gene log2 abundance (TU effect sqrt(rho), gene noise sqrt(1-rho))."""
    e = rng.standard_normal(n)
    t = e.copy()
    for i in range(1, n):
        t[i] = phi * t[i - 1] + np.sqrt(1 - phi**2) * e[i]
    tu = np.repeat(np.arange(n), _lengths(rng, n, kind))
    x = np.sqrt(rho) * t[tu] + np.sqrt(1 - rho) * rng.standard_normal(tu.size) + 20.0
    return tu, x


def _frame(tu: np.ndarray, x: np.ndarray, det: np.ndarray | None = None,
           imp: np.ndarray | None = None) -> pd.DataFrame:  # fmt: skip
    tus: dict[str, list[Gene]] = {}
    for g, (t, v) in enumerate(zip(tu, x, strict=True)):
        d = True if det is None else bool(det[g])
        tus.setdefault(f"TU_{int(t) + 1}", []).append(
            Gene(f"g{g + 1}", float(v) if d else None, detected=d,
                 imputed=False if imp is None else bool(imp[g]))
        )  # fmt: skip
    return frame(tus)


def _first_genes(tu: np.ndarray, pick: np.ndarray) -> list[str]:
    first = {int(t): g for g, t in reversed(list(enumerate(tu)))}
    return [f"g{first[int(t)] + 1}" for t in pick]


# --------------------------------------------------------------------------- arms


def _median_imputed(rng: np.random.Generator) -> float:
    n = 600
    tu, x = _world(rng, n, "1-4", 0.8)
    imp = rng.random(tu.size) < 0.2
    x = x.copy()
    x[imp] = np.median(x[~imp])
    a = int(rng.integers(0, n - 15))
    return _p(_frame(tu, x, imp=imp), _first_genes(tu, np.arange(a, a + 15)), rng)


def _size_biased(rng: np.random.Generator) -> float:
    tu, x = _world(rng, 1500, "ecoli", 0.2)
    g = rng.choice(tu.size, 40, replace=False)
    return _p(_frame(tu, x), [f"g{int(i) + 1}" for i in g], rng)


def _ar1(shape: str) -> Any:
    def rep(rng: np.random.Generator) -> float:
        n = 1500
        tu, x = _world(rng, n, "ecoli", 0.8, phi=0.5)
        if shape == "run15":
            a = int(rng.integers(0, n - 15))
            pick = np.arange(a, a + 15)
        else:
            pick = rng.choice(n, 15, replace=False)
        return _p(_frame(tu, x), _first_genes(tu, pick), rng)

    return rep


def _detectability_confounded(rng: np.random.Generator) -> float:
    n = 2000
    tu, x = _world(rng, n, "ecoli", 0.8)
    pick = rng.choice(n, 60, replace=False)
    ins = np.isin(tu, pick)
    pdet = 1 / (1 + np.exp(-2.0 * (x - 20.0 - np.where(ins, 1.5, 0.0) + 0.5)))
    det = rng.random(tu.size) < pdet
    genes = _frame(tu, x, det=det).assign(p_detect=pdet)
    return _p(genes, [f"g{int(i) + 1}" for i in np.flatnonzero(ins)], rng)


GATED = {"median_imputed": _median_imputed, "size_biased": _size_biased,
         "ar1_scattered": _ar1("scattered")}  # fmt: skip


FAST_FULL = [
    pytest.param(False, id="fast"),
    pytest.param(True, id="full", marks=pytest.mark.calibration_full),
]


@pytest.mark.parametrize("full", FAST_FULL)
@pytest.mark.parametrize("arm", sorted(GATED))
def test_attack_gsea_structure_at_level(arm: str, full: bool) -> None:
    levels = C.load_tolerances().levels
    reps = C.replicates(SUITE, full=full)
    rec = calibrate_levels(SUITE, "level", levels, GATED[arm], reps, stream=f"{SUITE}.attack.{arm}")
    rec.assert_ok()


@XFAIL_ADR
@pytest.mark.calibration_full
@pytest.mark.parametrize("arm", ["ar1_run15", "detectability_confounded"])
def test_attack_gsea_known_limitation(arm: str, teeth_log: Path) -> None:
    fn = _ar1("run15") if arm == "ar1_run15" else _detectability_confounded
    rec = C.calibrate(SUITE, "level", ALPHA, fn, C.replicates(SUITE, full=True),
                      stream=f"{SUITE}.attack.{arm}", log_dir=teeth_log)  # fmt: skip
    assert rec["status"] != "fail", rec
