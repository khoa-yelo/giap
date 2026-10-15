"""Calibration of giapsim's proteomics claims (issue M3-19; GIAP_SPEC.md §11.3, §8.4).

Claims checked over ``N_SEEDS`` = 20 replicate 100-kb genomes (seeds ``seeds.calibration_base + r``,
tests/calibration/tolerances.yaml):

* the simulated q-values are calibrated: mean FDP (vs truth ``is_wrong``) of the target rows with
  q <= q* is <= q* + ``se_multiplier`` SE for q* in ``levels`` (spec acceptance criterion);
* the wrong share of distinct target peptides is ``pi_wrong`` within 3 SE (pooled);
* PEPs are calibrated: the pooled sum of PEP equals the number of wrong rows within 4 SD;
* scores follow the two-component mixture: correct ~ N(3, 1), wrong ~ N(0, 1) (KS p >=
  ``ks_min_p``);
* detection is the logistic draw in length and log abundance: detected counts match the sum of
  the predicted probabilities over all tryptic products of the expressed proteins (within 4 SD,
  overall and in the low/high halves of the predicted probability);
* intensities are exp(log abundance + N(0, 0.5)): residuals have mean 0 (4 SE) and SD 0.5 ± 10 %.

Conventions as in tests/unit/test_sim_proteomics.py.
"""

from __future__ import annotations

import math
from typing import Any

import pytest
from helpers.peptide_oracle import Run, fnum, rows, truthy, trypsin_p
from scipy import stats

import giapsim

pytestmark = pytest.mark.calibration

N_SEEDS = 20
PI = 0.05
LOGIT = (-2.0, 0.08, 0.9)
FEATURES: dict[str, Any] = {"genome_length": 100_000, "n_sorfs": 3, "n_nterm_extensions": 2}


@pytest.fixture(scope="module")
def runs(tmp_path_factory: pytest.TempPathFactory, tolerances: dict[str, Any]) -> list[Run]:
    from giapsim.proteomics import ProteomicsConfig

    base = int(tolerances["seeds"]["calibration_base"])
    cfg = giapsim.SimConfig(proteomics=ProteomicsConfig(pi_wrong=PI), **FEATURES)
    out = []
    for r in range(N_SEEDS):
        d = tmp_path_factory.mktemp(f"cal{r}")
        giapsim.simulate(cfg, base + r, d)
        out.append(Run(d))
    return out


def labelled_targets(run: Run) -> list[tuple[float, float, float, bool]]:
    """(score, pep, q, is_wrong) of every target row."""
    return [
        (float(r["score"]), float(r["pep"]), float(r["q_value"]), run.is_wrong(r["peptide"]))
        for r in run.target_rows()
    ]


def test_q_values_calibrated(runs: list[Run], tolerances: dict[str, Any]) -> None:
    k = float(tolerances["se_multiplier"])
    for q in tolerances["levels"]:
        fdp = []
        for run in runs:
            acc = [w for _, _, qv, w in labelled_targets(run) if qv <= q]
            fdp.append(sum(acc) / len(acc) if acc else 0.0)
        mean = sum(fdp) / len(fdp)
        sd = math.sqrt(sum((x - mean) ** 2 for x in fdp) / (len(fdp) - 1))
        assert mean <= q + k * sd / math.sqrt(len(fdp)), (q, mean, sd)


def test_wrong_share(runs: list[Run]) -> None:
    n = w = 0
    for run in runs:
        labels = {r["peptide"]: truthy(r["is_wrong"]) for r in run.truth if r["origin"] != "decoy"}
        n += len(labels)
        w += sum(labels.values())
    se = math.sqrt(PI * (1 - PI) / n)
    assert abs(w / n - PI) <= 3 * se, (w / n, se)


def test_pep_calibrated(runs: list[Run]) -> None:
    peps = [(p, wr) for run in runs for _, p, _, wr in labelled_targets(run)]
    expected = sum(p for p, _ in peps)
    observed = sum(wr for _, wr in peps)
    sd = math.sqrt(sum(p * (1 - p) for p, _ in peps))
    assert abs(observed - expected) <= 4 * sd + 1, (observed, expected, sd)


def test_score_mixture(runs: list[Run], tolerances: dict[str, Any]) -> None:
    rows_ = [x for run in runs for x in labelled_targets(run)]
    correct = [s for s, _, _, w in rows_ if not w]
    wrong = [s for s, _, _, w in rows_ if w]
    assert len(wrong) > 100
    p_min = float(tolerances["ks_min_p"])
    assert stats.kstest(correct, stats.norm(3.0, 1.0).cdf).pvalue >= p_min
    assert stats.kstest(wrong, stats.norm(0.0, 1.0).cdf).pvalue >= p_min


def log_abundance(run: Run) -> dict[str, float]:
    out = {}
    for r in rows(run.out / "abundance.tsv"):
        if "log_abundance" in r:
            out[r["locus_tag"]] = float(r["log_abundance"])
        else:
            out[r["locus_tag"]] = math.log(float(r["abundance"]))
    return out


def test_detectability_logistic(runs: list[Run]) -> None:
    b0, b1, b2 = LOGIT
    halves: dict[bool, list[float]] = {False: [0.0, 0.0, 0.0], True: [0.0, 0.0, 0.0]}
    for run in runs:
        la = log_abundance(run)
        detected = {r["peptide"] for r in run.truth if r["origin"] == "gene"}
        seen: set[str] = set()
        for loc, p in run.proteins.items():
            assert loc in la, f"{loc} has no abundance"
            for a, b in trypsin_p(p.aa):
                pep = p.aa[a:b]
                if pep in seen:
                    continue
                seen.add(pep)
                prob = 1 / (1 + math.exp(-(b0 + b1 * len(pep) + b2 * la[loc])))
                h = halves[prob >= 0.5]
                h[0] += pep in detected
                h[1] += prob
                h[2] += prob * (1 - prob)
    tot = [halves[False][i] + halves[True][i] for i in range(3)]
    for obs, exp, var in (*halves.values(), tot):
        assert var > 0
        assert abs(obs - exp) <= 4 * math.sqrt(var) + 1, (obs, exp, var)


def test_intensity_lognormal(runs: list[Run]) -> None:
    res = []
    for run in runs:
        la = log_abundance(run)
        for r in run.target_rows():
            if run.is_wrong(r["peptide"]):
                continue
            locs = {t["locus_tag"] for t in run.truth_by_peptide[r["peptide"]]}
            if len(locs) != 1:
                continue
            x = fnum(r["intensity"])
            assert x is not None and x > 0
            res.append(math.log(x) - la[locs.pop()])
    n = len(res)
    mean = sum(res) / n
    sd = math.sqrt(sum((x - mean) ** 2 for x in res) / (n - 1))
    assert abs(mean) <= 4 * 0.5 / math.sqrt(n), mean
    assert 0.45 <= sd <= 0.55, sd


def test_charges_both_present(runs: list[Run]) -> None:
    charges = [r["charge"] for run in runs for r in run.peptides]
    assert set(charges) == {"2", "3"}
