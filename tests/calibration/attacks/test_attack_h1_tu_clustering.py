"""Attack (GIAP_SPEC.md §8.5, statistician): TU-clustered discordance for the M10-05 H1 rule
(docs/specs/M10-05-pilot-metrics.md "Hypotheses 1 and 3", review S2; §8.2 "resample TUs, not
genes").

Null world (giapsim ``simulate_layout``, one layout per scenario): about 4,000 genes (3,700 on
``chr1``, plasmids of 120 and 180 genes, four 80-gene prophages and two 40-gene BGCs) in TUs of
geometric size with mean 2.5, so about 1,600 TUs. Half of the singleton TUs get an empty
``tu_id`` (the "a gene is its own unit" path). Per replicate, ``N_DISCORDANT`` TUs are discordant.
Each one has a fair-coin direction (A2-only or A1-only) that all its discordant genes share; each
gene of the TU is discordant with probability 0.8, and at least one is. The other TUs are
co-detected as a block (both arms or neither, with gene-level dropout that stays concordant). A
sprinkle of 1 % independent fair-coin gene discordance creates mixed TUs. The non-detected arm of
a discordant gene is ``not_detected`` or ``single_peptide_unestimated`` (neither counts). H0 is
true: the arms have equal marginal detection, and every discordant TU is a fair coin.

Scenarios:

- ``uniform``: the discordant TUs are drawn from all TUs.
- ``regional``: they are drawn from TUs with a plasmid, prophage or BGC gene, so the "mostly"
  clause holds almost always and only ``direction`` and ``test`` stand between the null and a
  ``supported`` verdict. This is the scenario a pilot with a divergent accessory genome looks like.

Gated (each scenario, R = 10 x ``replicates_fast.compare`` = 1,000; the full tier 10 x
``replicates_full.compare`` = 10,000, marked ``calibration_full``, Slurm):

- ``p_tu``'s rejection rate at alpha in ``levels`` (0.01, 0.05, 0.10) is at most
  ``level_bound(alpha, R)`` (§8.4);
- the H1 verdict ``supported`` occurs in at most ``level_bound(0.05, R)`` of replicates (the rule
  as a whole does not report support beyond its nominal rate).

Recorded (documentation, spec "Test plan"; ``compare_attack_h1_recorded.jsonl``): ``p_gene``'s
rejection rates, which the TU clustering pushes above level. One teeth check keeps the attack
honest: in the ``uniform`` scenario ``p_gene`` at 0.05 must exceed ``level_bound(0.05, R)``. If it
did not, the null would not be clustered enough to tell a TU-level test from a gene-level one.

Calls the implementation, ``validation.pilot.metrics.hypothesis1``, unchanged. Randomness only
through ``giap.rng.for_step`` seeded from ``tolerances.yaml`` (``_criteria.suite_seed``).
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest

from giap.rng import for_step
from giapsim import layout as la

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from helpers.pilot_metrics_m10 import load_metrics

SUITE = "compare"
pytestmark = [pytest.mark.calibration, pytest.mark.calsuite(SUITE)]
MODES = [pytest.param(False, id="fast"),
         pytest.param(True, marks=pytest.mark.calibration_full, id="full")]  # fmt: skip
SCENARIOS = ("uniform", "regional")
N_DISCORDANT = {"uniform": 300, "regional": 200}
P_GENE_DISCORDANT = 0.8  # share of a discordant TU's genes that are discordant
P_TU_DETECTED = 0.55  # concordant TUs detected in both arms (else in neither)
P_DROPOUT = 0.15  # gene of a co-detected TU detected in neither arm
P_SPRINKLE = 0.01  # independent fair-coin gene discordance (mixed TUs)
LAYOUT = la.LayoutConfig(
    n_genes_chromosome=3700,
    plasmids=(120, 180),
    mean_tu_size=2.5,
    regions=(*(la.RegionSpec("prophage", 80),) * 4, *(la.RegionSpec("bgc", 40),) * 2),
)
REGION_PREFIXES = ("prophage_", "bgc_")
NOT_DETECTED = np.array(["not_detected", "single_peptide_unestimated"], dtype=object)


@dataclass(frozen=True)
class World:
    """One giapsim layout as the H1 inputs need it."""

    tags: np.ndarray  # locus tags, layout order
    feature_type: np.ndarray
    tu_code: np.ndarray  # TU index per gene
    tu_id: np.ndarray  # written tu_id (empty for half of the singletons)
    subsets: dict[str, set[str]]
    pool: dict[str, np.ndarray]  # scenario -> TU indices eligible as discordant


def _world(scenario: str) -> World:
    rng = for_step(f"attack.h1_tu_clustering.layout.{scenario}", C.suite_seed(SUITE, 0))
    lay = la.simulate_layout(LAYOUT, rng)
    g = lay.genes.reset_index(drop=True)
    truth = lay.truth.set_index("locus_tag").loc[g["locus_tag"]].reset_index()
    codes, uniq = pd.factorize(truth["tu_id"])
    size = np.bincount(codes)
    blank = (size[codes] == 1) & (rng.random(len(codes)) < 0.5)
    tu_id = np.where(blank, "", truth["tu_id"].astype(str).to_numpy())
    tags = g["locus_tag"].astype(str).to_numpy()
    plasmid = g["contig"].astype(str).str.startswith("p").to_numpy()
    regions = [list(r) for r in truth["region_ids"]]
    pro = np.array([any(str(x).startswith("prophage_") for x in r) for r in regions])
    bgc = np.array([any(str(x).startswith("bgc_") for x in r) for r in regions])
    in_region = plasmid | pro | bgc
    regional = np.unique(codes[in_region])
    assert len(regional) >= 1.2 * N_DISCORDANT["regional"], len(regional)
    return World(
        tags=tags,
        feature_type=g["feature_type"].astype(str).to_numpy(),
        tu_code=codes,
        tu_id=tu_id,
        subsets={
            "plasmid": set(tags[plasmid]),
            "prophage": set(tags[pro]),
            "bgc": set(tags[bgc]),
            "strain_specific": set(),
        },
        pool={"uniform": np.arange(len(uniq)), "regional": regional},
    )


def _replicate(w: World, scenario: str, rng: np.random.Generator) -> tuple[Any, Any]:
    """Genes tables of A1 and A2 under H0 (TU-clustered, fair-coin discordance)."""
    n, n_tu = len(w.tags), int(w.tu_code.max()) + 1
    det_tu = rng.random(n_tu) < P_TU_DETECTED
    both = det_tu[w.tu_code] & (rng.random(n) >= P_DROPOUT)
    a1, a2 = both.copy(), both.copy()
    disc_tu = rng.choice(w.pool[scenario], size=N_DISCORDANT[scenario], replace=False)
    direction = np.zeros(n_tu, dtype=np.int8)  # 0 none, 1 A2-only, 2 A1-only
    direction[disc_tu] = np.where(rng.random(len(disc_tu)) < 0.5, 1, 2)
    disc = (direction[w.tu_code] > 0) & (rng.random(n) < P_GENE_DISCORDANT)
    for t in disc_tu:  # every discordant TU has at least one discordant gene
        members = np.flatnonzero(w.tu_code == t)
        if not disc[members].any():
            disc[rng.choice(members)] = True
    sprinkle = ~disc & (rng.random(n) < P_SPRINKLE)
    coin = rng.random(n) < 0.5
    to_a2 = (disc & (direction[w.tu_code] == 1)) | (sprinkle & coin)
    to_a1 = (disc & (direction[w.tu_code] == 2)) | (sprinkle & ~coin)
    a2[to_a2], a1[to_a2] = True, False
    a1[to_a1], a2[to_a1] = True, False

    def frame(det: np.ndarray) -> pd.DataFrame:
        status = np.where(det, "detected", NOT_DETECTED[rng.integers(0, 2, n)])
        return pd.DataFrame(
            {"locus_tag": w.tags, "feature_type": w.feature_type, "tu_id": w.tu_id,
             "detection_status": status}
        )  # fmt: skip

    return frame(a1), frame(a2)


def _run(scenario: str, R: int) -> dict[str, np.ndarray]:
    m = load_metrics()
    w = _world(scenario)
    p_tu, p_gene, supported, mixed = (np.empty(R) for _ in range(4))
    for i in range(R):
        rng = for_step(f"attack.h1_tu_clustering.{scenario}", C.suite_seed(SUITE, i))
        g1, g2 = _replicate(w, scenario, rng)
        h = m.hypothesis1(g1, g2, None, w.subsets)
        p_tu[i] = 1.0 if h["p_tu"] is None else h["p_tu"]
        p_gene[i] = 1.0 if h["p_gene"] is None else h["p_gene"]
        supported[i] = h["verdict"] == "supported"
        mixed[i] = h["mixed_tu"]
    return {"p_tu": p_tu, "p_gene": p_gene, "supported": supported, "mixed_tu": mixed}


def _log(record: dict[str, object]) -> None:
    """Recorded (not gated) numbers, kept out of ``compare.jsonl`` so the PR summary shows no
    FAIL row for a statistic that is expected above level."""
    path = C.log_dir() / f"{SUITE}_attack_h1_recorded.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(record) + "\n")


@pytest.mark.parametrize("scenario", SCENARIOS)
@pytest.mark.parametrize("full", MODES)
def test_h1_tu_level_and_verdict_under_tu_clustered_null(full: bool, scenario: str) -> None:
    R = 10 * C.replicates(SUITE, full=full)
    t0 = time.perf_counter()
    out = _run(scenario, R)
    dt = time.perf_counter() - t0
    seed = C.suite_seed(SUITE, 0)
    assert out["mixed_tu"].mean() > 0, "the sprinkle must create mixed TUs (tie path exercised)"
    failures = []
    for alpha in C.load_tolerances().levels:
        res = C.check_level(out["p_tu"], alpha)
        C.log_result(SUITE, f"h1_attack.{scenario}.p_tu", alpha, res, seed, dt)
        if not res.ok:
            failures.append(res.message)
        gene = C.check_level(out["p_gene"], alpha)  # recorded, expected above level
        _log({"scenario": scenario, "statistic": "p_gene", "alpha": alpha, "R": R,
              "rate": gene.observed, "bound": gene.bound, "full": full})  # fmt: skip
    rate = float(out["supported"].mean())
    bound = C.level_bound(0.05, R)
    verdict = C.CheckResult(
        rate <= bound, rate, bound, R, f"supported rate {rate:.4g} vs {bound:.4g}"
    )
    C.log_result(SUITE, f"h1_attack.{scenario}.supported", 0.05, verdict, seed, dt)
    if not verdict.ok:
        failures.append(verdict.message)
    assert not failures, failures


@pytest.mark.parametrize("full", MODES)
def test_teeth_gene_level_p_is_anticonservative(full: bool) -> None:
    R = 10 * C.replicates(SUITE, full=full)
    out = _run("uniform", R)
    rate = float(np.mean(out["p_gene"] <= 0.05))
    assert rate > C.level_bound(0.05, R), (
        f"p_gene rate {rate:.4g} within level: the null is not TU-clustered enough to bite"
    )
