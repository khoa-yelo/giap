"""Calibration of gippasim's de novo error claims (issue M3-20; GIPPA_SPEC.md §11.3, §8.4).

Claims checked over ``N_SEEDS`` replicate 60-kb genomes (seeds ``seeds.calibration_base + r``,
tests/calibration/tolerances.yaml):

* ``truth/denovo_precision.tsv`` matches the empirical precision per score decile within 3 SE
  (spec acceptance criterion), pooled over replicates, in every bin;
* the stated precision is the configured one (``precision_by_decile``) within 0.02;
* error types are drawn with the configured weights: among erroneous rows the share of each
  always-applicable type (substitution, swap, reversed_stretch, random) is its weight within
  3 SE (multinomial), with the other weights set to 0 so no fallback can distort the shares;
* with the default weights every type appears and ``block`` + ``near_block`` +
  ``modification`` together are within 3 SE of their weight or below it (fallbacks only move
  mass from site-dependent types to others, never the reverse).

Conventions as in tests/unit/test_sim_denovo.py. The M4-05 harness (``SUITE_REGISTRY`` in
tests/calibration/conftest.py) is not on main yet, so this suite is not registered there.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from helpers.denovo_oracle import decile_bins

import gippasim

pytestmark = pytest.mark.calibration

N_SEEDS = 8
SE_MULT = 3.0  # spec: "within 3 SE"
PRECISION = tuple(float(x) for x in np.round(np.linspace(0.3, 0.9, 10), 6))
ALWAYS = {"substitution": 0.4, "swap": 0.2, "reversed_stretch": 0.25, "random": 0.15}


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _sim(d: Path, seed: int, **dkw: Any) -> dict[str, list[dict[str, str]]]:
    from gippasim.denovo import DenovoConfig
    from gippasim.proteomics import ProteomicsConfig

    cfg = gippasim.SimConfig(
        genome_length=60_000,
        proteomics=ProteomicsConfig(pi_wrong=0.05),
        denovo=DenovoConfig(**dkw),
    )
    gippasim.simulate(cfg, seed, d)
    return {
        "denovo": _rows(d / "peptides_denovo.tsv"),
        "truth": _rows(d / "truth" / "denovo.tsv"),
        "precision": _rows(d / "truth" / "denovo_precision.tsv"),
    }


@pytest.fixture(scope="module")
def base(tolerances: dict[str, Any]) -> int:
    return int(tolerances["seeds"]["calibration_base"])


@pytest.fixture(scope="module")
def prec_runs(tmp_path_factory: pytest.TempPathFactory, base: int) -> list[Any]:
    return [
        _sim(tmp_path_factory.mktemp(f"p{r}"), base + r, precision_by_decile=PRECISION)
        for r in range(N_SEEDS)
    ]


def _bin_counts(run: dict[str, list[dict[str, str]]]) -> tuple[np.ndarray, np.ndarray]:
    truth = {(t["run"], t["spectrum_id"]): t["correct"].lower() == "true" for t in run["truth"]}
    s = np.array([float(r["score"]) for r in run["denovo"]])
    c = np.array([truth[(r["run"], r["spectrum_id"])] for r in run["denovo"]])
    b = decile_bins(s)
    return np.bincount(b, minlength=10), np.bincount(b, weights=c, minlength=10)


def test_stated_precision_is_configured(prec_runs: list[Any]) -> None:
    for run in prec_runs:
        stated = [float(r["precision_true"]) for r in run["precision"]]
        assert np.allclose(stated, PRECISION, atol=0.02), stated


def test_empirical_precision_per_bin_within_3se(prec_runs: list[Any]) -> None:
    n_tot, k_tot, p_exp = np.zeros(10), np.zeros(10), np.zeros(10)
    for run in prec_runs:
        n, k = _bin_counts(run)
        stated = np.array([float(r["precision_true"]) for r in run["precision"]])
        n_tot += n
        k_tot += k
        p_exp += n * stated
    p = p_exp / n_tot
    se = np.sqrt(p * (1 - p) / n_tot)
    emp = k_tot / n_tot
    bad = [(b, emp[b], p[b], se[b]) for b in range(10) if abs(emp[b] - p[b]) > SE_MULT * se[b]]
    assert not bad, bad


def test_error_type_shares_follow_weights(
    tmp_path_factory: pytest.TempPathFactory, base: int
) -> None:
    w = {
        "substitution": 0.0, "swap": 0.0, "block": 0.0, "near_block": 0.0,
        "reversed_stretch": 0.0, "modification": 0.0, "random": 0.0, **ALWAYS,
    }  # fmt: skip
    counts = dict.fromkeys(ALWAYS, 0)
    for r in range(N_SEEDS // 2):
        run = _sim(tmp_path_factory.mktemp(f"w{r}"), base + 100 + r, error_weights=w)
        for t in run["truth"]:
            if t["error_type"] != "none":
                counts[t["error_type"]] += 1
    n = sum(counts.values())
    assert n > 500, n
    for et, p in ALWAYS.items():
        se = math.sqrt(p * (1 - p) / n)
        assert abs(counts[et] / n - p) <= SE_MULT * se, (et, counts[et] / n, p, n)


def test_default_weights_site_dependent_share(
    tmp_path_factory: pytest.TempPathFactory, base: int
) -> None:
    from gippasim.denovo import DEFAULT_ERROR_WEIGHTS

    site = ("block", "near_block", "modification")
    n = k = 0
    seen: set[str] = set()
    for r in range(N_SEEDS // 2):
        run = _sim(tmp_path_factory.mktemp(f"d{r}"), base + 200 + r)
        for t in run["truth"]:
            if t["error_type"] != "none":
                n += 1
                k += t["error_type"] in site
                seen.add(t["error_type"])
    assert seen == set(DEFAULT_ERROR_WEIGHTS), seen
    p = sum(DEFAULT_ERROR_WEIGHTS[t] for t in site)
    assert k / n <= p + SE_MULT * math.sqrt(p * (1 - p) / n), (k / n, p)
    assert k / n >= 0.5 * p, (k / n, p)  # fallbacks are rare
