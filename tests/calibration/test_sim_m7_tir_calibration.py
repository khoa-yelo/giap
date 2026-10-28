"""Calibration of the M7-12 round-2 planted-TIR knob on the ``sim300k`` size (lead decisions 13 and
14 after the statistician's block on PR #225; GIPPA_SPEC.md §7.6.6, §8.4 suite "Novel-ORF class
FDR").

Claims checked (each pooled over replicates of ``novel_orf_power(seed, "300k")``, seeds
``seeds.calibration_base + 3000 + r`` (fast) / ``+ 4000 + r`` (full); interface as in
tests/unit/test_sim_m7_round2.py):

* **Distribution (decision 14).** The log ``tir_ostir`` of ``precomputed/rbs_novel.tsv`` at the
  planted starts (``rbs_planted = true`` rows of ``truth/novel_orfs.tsv`` at ``(contig, strand,
  true_start)``) is N(mu, sigma) of ``SimConfig.planted_tir``: one-sample KS, p >= ``ks_min_p``,
  for the default arm (8, 0.8) and the sanity arm (9, 0.5) (set by ``dataclasses.replace``).
* **Matched arm.** With the default (8, 0.8), planted log TIRs and the first-round log TIRs of the
  genes with a true planted RBS (``truth/rbs.tsv`` ``true_element = true``) are one distribution:
  two-sample KS, p >= ``ks_min_p``.
* **Estimable threshold (decision 13).** Every replicate has >= 200 annotated non-pseudo CDSs
  (decoy-picked q floor 1 / n < 0.01) and >= 20 detected genes with a TIR (M7-09 decision 1;
  detected = >= 1 correct peptide, the simulator's stand-in for GIPPA's rule).
* **Share calibration.** With threshold t_r = 10 % quantile (linear) of the detected genes'
  first-round TIRs, the pooled number of planted starts with TIR >= t_r is within 3 SE of its
  expectation sum_r n_r (1 - Phi((log t_r - mu) / sigma)) (SE from the binomial variances), in
  both arms. This is what lets M7-13 put an SE allowance on the matched arm.
* **Sanity bar (decision 8).** In the (9, 0.5) arm the pooled share is >= 0.90 (point estimate).

Fast: 3 replicates per arm at 300 kb. Full (``calibration_full``): ``replicates_fast.novel_orfs``
replicates per arm. Suite ``novel_orfs``.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from helpers import m7_sim_oracle as O
from scipy import stats

import gippasim

pytestmark = pytest.mark.calsuite("novel_orfs")

DEFAULT_ARM = (8.0, 0.8)
SANITY_ARM = (9.0, 0.5)
RBS_Q = 0.10
MIN_SHARE = 0.90
MIN_CDS = 200
MIN_GENES = 20
N_FAST = 3


@dataclasses.dataclass
class Rep:
    planted: list[float]  # TIR at planted starts
    genes_true: list[float]  # first-round TIR of genes with a true RBS
    thr: float
    n_cds: int
    n_detected: int


def _detected_tirs(out: Path) -> np.ndarray:
    data = json.loads((out / "bakta_like.json").read_text())
    lengths = {s["id"]: int(s["length"]) for s in data["sequences"]}
    ok = set()
    for f in data["features"]:
        if f["type"] not in ("cds", "sorf") or f.get("pseudo"):
            continue
        contig = f.get("contig", f.get("sequence"))
        if int(f["stop"]) < int(f["start"]) or int(f["stop"]) > lengths[contig]:
            continue
        ok.add(str(f["locus"]))
    detected = {r["locus_tag"] for r in O.correct_peptide_rows(out)} & ok
    tir: dict[str, float] = {}
    for r in O.rows(out / "precomputed" / "rbs.tsv"):
        if r["tir_ostir"] not in ("", "NA"):
            tir[r["locus_tag"]] = float(r["tir_ostir"])
    return np.array([tir[t] for t in sorted(detected) if t in tir], dtype=float)


def _rep(d: Path, seed: int, arm: tuple[float, float]) -> Rep:
    import gippasim.scenarios as sc

    cfg = sc.novel_orf_power(seed, "300k")
    if tuple(cfg.planted_tir) != arm:
        cfg = dataclasses.replace(cfg, planted_tir=arm)
    gippasim.simulate(cfg, seed, d)
    layer = {(r["contig"], r["strand"], int(r["start"])): float(r["tir_ostir"])
             for r in O.rows(d / "precomputed" / "rbs_novel.tsv")}  # fmt: skip
    rows = O.rows(d / "truth" / "novel_orfs.tsv")
    planted = [layer[(r["contig"], r["strand"], int(r["true_start"]))]
               for r in rows if O.truthy(r["rbs_planted"])]  # fmt: skip
    genes_true = [float(r["tir_ostir"]) for r in O.rows(d / "truth" / "rbs.tsv")
                  if O.truthy(r["true_element"]) and r["tir_ostir"] not in ("", "NA")]  # fmt: skip
    det = _detected_tirs(d)
    n_cds = sum(f["type"] == "cds" and not f.get("pseudo")
                for f in O.bakta_features(d / "bakta_like.json").values())  # fmt: skip
    thr = float(np.quantile(det, RBS_Q, method="linear")) if len(det) else float("nan")
    return Rep(planted, genes_true, thr, n_cds, len(det))


def _check(reps: list[Rep], arm: tuple[float, float], ks_min_p: float) -> float:
    mu, sigma = arm
    for r in reps:
        assert r.n_cds >= MIN_CDS and 1 / r.n_cds < 0.01, r.n_cds
        assert r.n_detected >= MIN_GENES, r.n_detected
        assert r.planted, "no planted start"
    logs = np.log(np.concatenate([r.planted for r in reps]))
    p1 = stats.kstest(logs, stats.norm(mu, sigma).cdf).pvalue
    assert p1 >= ks_min_p, (arm, p1, logs.mean(), logs.std())
    obs = sum(sum(t >= r.thr for t in r.planted) for r in reps)
    probs = np.concatenate([np.full(len(r.planted), stats.norm.sf((np.log(r.thr) - mu) / sigma))
                            for r in reps])  # fmt: skip
    exp, se = probs.sum(), float(np.sqrt((probs * (1 - probs)).sum()))
    assert abs(obs - exp) <= 3 * max(se, 0.5), (arm, obs, exp, se)
    if arm == DEFAULT_ARM:
        genes = np.log(np.concatenate([r.genes_true for r in reps]))
        p2 = stats.ks_2samp(logs, genes).pvalue
        assert p2 >= ks_min_p, (p2, logs.mean(), genes.mean())
    return obs / len(logs)


def _reps(tmp: pytest.TempPathFactory, seeds: list[int], arm: tuple[float, float]) -> list[Rep]:
    return [_rep(tmp.mktemp(f"m7tir{arm[0]:g}_{s}"), s, arm) for s in seeds]


@pytest.fixture(scope="module")
def base(tolerances: dict[str, Any]) -> int:
    return int(tolerances["seeds"]["calibration_base"])  # type: ignore[index]


@pytest.mark.calibration
@pytest.mark.parametrize("arm", [DEFAULT_ARM, SANITY_ARM], ids=["matched_8_0.8", "sanity_9_0.5"])
def test_planted_tir_calibration_300k_fast(
    tmp_path_factory: pytest.TempPathFactory,
    base: int,
    tolerances: dict[str, Any],
    arm: tuple[float, float],
) -> None:
    reps = _reps(tmp_path_factory, [base + 3000 + r for r in range(N_FAST)], arm)
    share = _check(reps, arm, float(tolerances["ks_min_p"]))  # type: ignore[arg-type]
    if arm == SANITY_ARM:
        assert share >= MIN_SHARE, share


@pytest.mark.calibration_full
@pytest.mark.parametrize("arm", [DEFAULT_ARM, SANITY_ARM], ids=["matched_8_0.8", "sanity_9_0.5"])
def test_planted_tir_calibration_300k_full(
    tmp_path_factory: pytest.TempPathFactory,
    base: int,
    tolerances: dict[str, Any],
    arm: tuple[float, float],
) -> None:
    n = int(tolerances["replicates_fast"]["novel_orfs"])  # type: ignore[index]
    reps = _reps(tmp_path_factory, [base + 4000 + r for r in range(n)], arm)
    share = _check(reps, arm, float(tolerances["ks_min_p"]))  # type: ignore[arg-type]
    if arm == SANITY_ARM:
        assert share >= MIN_SHARE, share
