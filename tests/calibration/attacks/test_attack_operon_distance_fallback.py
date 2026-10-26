"""Attack: the ADR-0022 distance-model operon fallback as the source of pair classes and TUs
(M11-14 14s1, docs/specs/M11-14-container-install.md "Statistical evidence that gates 14d2";
GIAP_SPEC.md §7.1.3, §8.1, §8.5; R3, R4). Statistician-owned.

When UniOP is absent or fails, ``operon_pairs`` scores every adjacent same-strand pair with
p_dist = 1 / (1 + exp((d - d0) / s)) (``calibration.yaml``: d0 = 102.69, s = 47.16) and links it
when p_dist >= ``annotation.uniop.pair_threshold`` (0.5, i.e. d <= 102 bp). The pair classes
(OP / NO) and the TUs, which are the exchangeability blocks of most §8 nulls, then come from
distance alone. No bound, alpha, R or tolerance differs from the gates the arms re-run (R3).

**Arm A, class re-derivation** (suites ``missed_operons`` and ``discordance``). giapsim
``operon_null``, ``operon_attack_ar1(rho 0.5)`` and ``operon_attack_detection(prophage)`` with
``uniop_p`` of every same-strand pair replaced by p_dist(``intergenic_bp``) (``source =
distance_model``) and the classes re-derived by the step's rule (``_operon_gate.reclass``). Then:

- the M6-09 missed-operon decoy-FDR gate (shipped decoys and the non-vacuous convergent-only arm
  with 100 planted pairs at 10^-3 SD; on ``operon_null``, which plants nothing, the non-vacuous
  arm is M6-09's ``missed_operon_power``), at every q, R = ``missed_operons``. A planted pair that the
  distance model links leaves the NO family and the truth list (it is then an OP pair);
- the M6-09 discordance checks (pooled level, KS uniformity of one p per replicate, BH family
  FDR), R = ``discordance``. The null is now *pair-level*: a re-derived OP pair that is not a true
  operon pair spans a TU boundary, so a promoter and a terminator lie between its genes
  (giapsim truth ``promoter_between_true``) and calling it discordant is a true discovery. Level
  and KS use the tested pairs with ``true_operon`` and no element between; the FDR counts a call
  as false on those pairs only. On giapsim layouts (intra-TU gaps -20...40 bp, inter-TU 50...300
  bp) the model never splits a TU, it only merges (inter-TU gaps <= 102 bp), so arm A attacks the
  merges; splits are arm B's.

**Arm B, K-12-like gaps** (suites ``gsea``, ``ora``, ``determinants``, ``autocorrelation``). Every
same-strand adjacent pair gets an intergenic distance drawn from the class-conditional ADR-0022 bin
table (operon pairs from the operon counts per bin, non-operon pairs from the rest; uniform
integer inside a bin, the open bins closed at -60 and 600 bp), so the distance model splits a true
operon pair with probability 0.084 and merges a same-strand TU boundary with probability 0.325,
the K-12 rates (``test_arm_b_plants_k12_rates``). Abundance, detection and set membership are
generated from the **true** TUs; the **distance-model TUs** are the analysis blocks (``tu_id``,
the TU collapse, the TU folds, the bootstrap blocks, the permutation blocks). In the abstract
worlds a TU boundary is same-strand with probability ``P_SAME`` (giapsim's default), otherwise it
is never linked. The arm feeds:

- GSEA: the gated arms of ``test_attack_gsea_structure`` (a) median imputation, (b) size-biased
  membership, (c) between-TU AR(1) with a scattered set, at level, R = ``gsea``. The contiguous
  set of (a) is cut at 15 *analysis* TUs (merges would otherwise drop it below the size filter);
- ORA: ``test_attack_ora_tusize`` 4b (TU-packed detection, c 1.5 and 2.5, ``large`` and
  ``single`` categories of whole true TUs, both tails), at level after the decision-12
  recalibration, with that test's ``ora_uncalibrated`` fallback rule, R = ``ora``;
- determinants: ``test_attack_determinants_tu_null::test_univariate_ci_level`` (95 % TU-block
  bootstrap CIs of Spearman rho and AUROC of a TU-constant feature; the feature is constant
  within the *analysis* TU, as the pipeline computes ``promotech_tu`` on the TUs it built), with
  that test's ``ci_unestimated`` rule, R = 5 x ``determinants``;
- autocorrelation: the TU-block null of ``context.autocorrelation`` on ``scenarios.tu_null``
  (tu_sd 2, decay 0.2, detection on abundance), both statistics at lags 1 and 20, R =
  ``autocorrelation``.

Each arm B check runs with ``source = distance_model`` (gated) and once with UniOP-like errors
(``uniop_like``: giapsim's Beta(8, 2) / Beta(2, 8) classes, recorded, not gated) on the same
worlds (the world stream does not depend on the source), so the two sources are compared on the
same data.

**Outcome rule** (spec 14s1): an arm that passes at the existing bounds clears the fallback for
that claim. A failing arm is a strict xfail naming the scope ADR: the affected rows must carry
``tu_source_distance_model`` with an R4 caveat and are scoped as uncalibrated under the fallback.
The bound never moves.
"""

from __future__ import annotations

import json
import sys
import warnings
from collections.abc import Callable
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import Any

import _criteria as C
import _operon_gate as G
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import autocorr_calib as acal
import autocorr_data as ad
import determinants_models_data as dm
from helpers.gsea_data import Gene, as_list, frame, gene_sets
from helpers.ora_calsim import K_FOLDS, collection, simulate, whole_tu_sets

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
NIGHTLY_MODES = [
    pytest.param(False, marks=[pytest.mark.calibration, pytest.mark.calibration_nightly],
                 id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]  # fmt: skip
MO, DI = "missed_operons", "discordance"
SOURCES = ("distance_model", "uniop_like")
GATED_SOURCE = "distance_model"
M_POOLED = 200
#: same-strand probability of a TU boundary in the abstract worlds (giapsim LayoutConfig default)
P_SAME = 0.6
#: giapsim's UniOP-like classes (LayoutConfig uniop_beta_in / uniop_beta_out)
BETA_IN, BETA_OUT = (8.0, 2.0), (2.0, 8.0)

#: ADR-0022 binned table (RegulonDB 14.5.0 K-12 labelled pairs): (lo, hi) bp (half-open; the open
#: end bins closed at -60 and 600), n pairs, observed operon share
K12_BINS: tuple[tuple[int, int, int, float], ...] = (
    (-60, -20, 47, 0.745), (-20, 0, 566, 0.979), (0, 20, 576, 0.965), (20, 40, 222, 0.914),
    (40, 60, 173, 0.549), (60, 80, 172, 0.442), (80, 100, 136, 0.449), (100, 150, 289, 0.273),
    (150, 200, 165, 0.164), (200, 300, 218, 0.119), (300, 600, 175, 0.103),
)  # fmt: skip


def _model() -> tuple[float, float, float]:
    from giap.config import GiapConfig
    from giap.data.calibration import load

    od = load()["operon_distance"]
    return float(od["d0"]), float(od["s"]), float(GiapConfig().annotation.uniop.pair_threshold)


def p_dist(d: np.ndarray) -> np.ndarray:
    """The step's ``distance_model_p`` at the bundled d0 and s, per pair."""
    from giap.annotate.operons import distance_model_p

    d0, s, _ = _model()
    return np.array([distance_model_p(int(v), d0, s) for v in np.asarray(d)], dtype=float)


def _log(name: str, row: dict[str, Any]) -> None:
    path = C.log_dir() / f"operon_fallback_{name}.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(row) + "\n")


# =========================================================================== arm B machinery


def _bin_counts() -> tuple[np.ndarray, np.ndarray]:
    n = np.array([b[2] for b in K12_BINS], dtype=float)
    op = np.round(n * np.array([b[3] for b in K12_BINS]))
    return op, n - op


def k12_gaps(rng: np.random.Generator, operon: np.ndarray) -> np.ndarray:
    """One intergenic distance per pair from its class's ADR-0022 bin counts."""
    op, no = _bin_counts()
    lo = np.array([b[0] for b in K12_BINS])
    hi = np.array([b[1] for b in K12_BINS])
    operon = np.asarray(operon, dtype=bool)
    k = np.where(operon, rng.choice(len(op), operon.size, p=op / op.sum()),
                 rng.choice(len(no), operon.size, p=no / no.sum()))  # fmt: skip
    return rng.integers(lo[k], hi[k])


def expected_link_rates() -> tuple[float, float]:
    """(split rate of operon pairs, merge rate of same-strand non-operon pairs) implied by the
    bins under the distance model (uniform integers inside a bin)."""
    op, no = _bin_counts()
    _, _, t = _model()
    share = []
    for lo, hi, _, _ in K12_BINS:
        d = np.arange(lo, hi)
        share.append(float((p_dist(d) >= t).mean()))
    w = np.array(share)
    return float(1 - (op * w).sum() / op.sum()), float((no * w).sum() / no.sum())


def link_scores(rng: np.random.Generator, operon: np.ndarray, source: str) -> np.ndarray:
    """Pair scores of ``source`` for same-strand pairs of the given true class."""
    operon = np.asarray(operon, dtype=bool)
    if source == "distance_model":
        return p_dist(k12_gaps(rng, operon))
    if source == "uniop_like":
        return np.where(operon, rng.beta(*BETA_IN, operon.size), rng.beta(*BETA_OUT, operon.size))
    raise ValueError(source)


def analysis_tus(true_tu: np.ndarray, same_next: np.ndarray, rng: np.random.Generator,
                 source: str) -> np.ndarray:  # fmt: skip
    """Analysis TU code per gene (genome order) from the true TU codes: gene j and j + 1 are
    linked iff they are on the same strand and the pair score is >= pair_threshold. Genes of a
    true TU are on one strand; ``same_next[j]`` is the strand relation at the boundaries."""
    _, _, t = _model()
    tu = np.asarray(true_tu)
    operon = tu[1:] == tu[:-1]
    cand = operon | np.asarray(same_next, dtype=bool)
    score = np.zeros(operon.size)
    score[cand] = link_scores(rng, operon[cand], source)
    link = cand & (score >= t)
    return np.concatenate([[0], np.cumsum(~link)])


def abstract_blocks(true_tu: np.ndarray, rng: np.random.Generator, source: str) -> np.ndarray:
    """``analysis_tus`` on an abstract world: a boundary is same-strand with ``P_SAME``."""
    tu = np.asarray(true_tu)
    same = (tu[1:] == tu[:-1]) | (rng.random(tu.size - 1) < P_SAME)
    return analysis_tus(tu, same, rng, source)


def _links_rng(stream: str, i: int, source: str) -> np.random.Generator:
    return C.replicate_rng(f"{stream}.links.{source}", i)


def _index(stream: str, R: int) -> dict[int, int]:
    return {C.suite_seed(stream, i): i for i in range(2 * R)}


def _rep_index(rng: np.random.Generator, index: dict[int, int]) -> int:
    return index[int(rng.bit_generator.seed_seq.entropy)]  # type: ignore[attr-defined]


def _assert_or_record(source: str, failed: list[Any], detail: list[str], what: str) -> None:
    """Gate the distance model; record the UniOP-like comparison (printed and logged)."""
    print(f"{what} [{source}]:\n" + "\n".join(detail))
    _log("arm_b", {"check": what, "source": source, "failed": bool(failed), "detail": detail})
    if source == GATED_SOURCE:
        assert not failed, detail


def test_arm_b_plants_k12_rates() -> None:
    """The arm really splits and merges at the K-12 rates (and UniOP-like at its Beta rates)."""
    split, merge = expected_link_rates()
    assert abs(split - 0.084) < 0.005 and abs(merge - 0.325) < 0.01, (split, merge)
    rng = np.random.default_rng(14)
    _, _, t = _model()
    n = 200_000
    for source, want in (("distance_model", (split, merge)), ("uniop_like", (None, None))):
        s_op = link_scores(rng, np.ones(n, dtype=bool), source)
        s_no = link_scores(rng, np.zeros(n, dtype=bool), source)
        got = (float((s_op < t).mean()), float((s_no >= t).mean()))
        if want[0] is None:
            from scipy.stats import beta

            want = (float(beta.cdf(t, *BETA_IN)), float(beta.sf(t, *BETA_OUT)))
        se = [np.sqrt(w * (1 - w) / n) for w in want]
        assert all(abs(g - w) < 5 * e for g, w, e in zip(got, want, se, strict=True)), (
            source,
            got,
            want,
        )
    tu = np.repeat(np.arange(5000), rng.geometric(0.48, 5000))
    blocks = abstract_blocks(tu, rng, "distance_model")
    assert blocks.size == tu.size and np.all(np.diff(blocks) >= 0)
    assert len(np.unique(blocks)) != len(np.unique(tu))


# =========================================================================== arm A


@dataclass(frozen=True)
class _SimA:
    layout: Any
    expression: Any
    pairs: pd.DataFrame
    truth: dict[str, pd.DataFrame]
    n_planted_linked: int


def simulate_a(cfgs: Any, seed: int) -> _SimA:
    """One giapsim operon replicate with distance-model pair scores and re-derived classes."""
    import giapsim.scenarios as sc

    s = sc.simulate_operons(cfgs, int(seed))
    pairs = s.pairs.copy()
    same = pairs["class"].isin(["OP", "NO"]).to_numpy()
    pairs.loc[same, "uniop_p"] = p_dist(pairs.loc[same, "intergenic_bp"].to_numpy())
    pairs["source"] = np.where(same, "distance_model", None)
    _, _, t = _model()
    pairs = G.reclass(pairs, t)
    truth = dict(s.truth)
    planted = truth["missed_operons"]
    keep = set(pairs.loc[pairs["class"] == "NO", "pair_id"].astype(str))
    still = planted["pair_id"].astype(str).isin(keep).to_numpy() if len(planted) else []
    truth["missed_operons"] = planted[still].reset_index(drop=True) if len(planted) else planted
    return _SimA(s.layout, s.expression, pairs, truth, int(len(planted) - int(np.sum(still))))


ARM_A_SCENARIOS: dict[str, tuple[str, tuple[Any, ...]]] = {
    "operon_null": ("operon_null", (G.N_GENES,)),
    "ar1_rho0.5": ("operon_attack_ar1", (G.N_GENES, 0.5)),
    "detection_prophage": ("operon_attack_detection", (G.N_GENES, "prophage")),
}


def _missed_records(stream: str, cfgs: Any, gcfg: Any) -> Callable[[int], dict[str, Any]]:
    @cache
    def records(i: int) -> dict[str, Any]:
        sim = simulate_a(cfgs, C.suite_seed(stream, i))
        out = G.missed_counts(sim, gcfg, G.levels())
        out["n_planted_linked"] = sim.n_planted_linked
        return out

    return records


@pytest.mark.calsuite(MO)
@pytest.mark.parametrize("full", NIGHTLY_MODES)
@pytest.mark.parametrize("arm", ["shipped", "nonvacuous"])
@pytest.mark.parametrize("scen", sorted(ARM_A_SCENARIOS))
def test_arm_a_missed_operon_fdr(scen: str, arm: str, full: bool) -> None:
    name, args = ARM_A_SCENARIOS[scen]
    if arm == "shipped":
        cfgs, gcfg = G.scenario(name, *args), G.giap_cfg()
    else:
        k, sd = G.NONVAC_PLANT
        # operon_null has no planting: its non-vacuous arm is M6-09's missed_operon_power
        cfgs = (G.scenario("missed_operon_power", G.N_GENES, k, sd) if scen == "operon_null"
                else G.scenario(name, *args, n_planted=k, effect_sd=sd))  # fmt: skip
        gcfg = G.giap_cfg(G.CV_ONLY)
    stream = f"{MO}.fallback_a.{scen}.{arm}"
    R = C.replicates(MO, full)
    records = _missed_records(stream, cfgs, gcfg)
    recs = G.missed_fdr_gate(MO, stream, records, R)
    if C.owns_direct_checks():
        s = G.fdp_summary(records, R, 0.1)
        s["planted_linked_mean"] = float(np.mean([records(i)["n_planted_linked"]
                                                  for i in range(R)]))  # fmt: skip
        print(f"{stream} at q 0.1: {s}")
        _log("arm_a_missed", {"stream": stream, **s})
    G.assert_fdr_ok(recs, allow_vacuous=True, need_calls_at=0.1 if arm == "nonvacuous" else None)


@dataclass(frozen=True)
class _DiscA:
    p: np.ndarray  # both_detected_t p of the null tested pairs
    q_null: np.ndarray  # their BH q (BH over the whole family)
    q_alt: np.ndarray  # BH q of the tested pairs spanning a TU boundary (true discoveries)


def discordance_a(sim: _SimA) -> _DiscA:
    from giap.config import GiapConfig
    from giap.operons.discordance import discordance

    res = discordance(sim.pairs, sim.expression.genes[["locus_tag"]], GiapConfig(), None)
    cols = ["pair_id", "true_operon", "promoter_between_true", "terminator_between_true"]
    rows = res.rows.merge(sim.pairs[cols], on="pair_id", how="left")
    t = (rows["test"] == "both_detected_t").to_numpy()
    null = (rows["true_operon"].astype(bool) & ~rows["promoter_between_true"].astype(bool)
            & ~rows["terminator_between_true"].astype(bool)).to_numpy()  # fmt: skip
    p = rows["p_value"].astype(float).to_numpy()
    q = rows["q_value"].astype(float).to_numpy()
    return _DiscA(p[t & null], q[t & null], q[t & ~null])


@cache
def _disc_run(scen: str, seed: int) -> _DiscA:
    name, args = ARM_A_SCENARIOS[scen]
    return discordance_a(simulate_a(G.scenario(name, *args), seed))


def _seed(rng: np.random.Generator) -> int:
    return int(rng.bit_generator.seed_seq.entropy)  # type: ignore[attr-defined]


@pytest.mark.calsuite(DI)
@pytest.mark.parametrize("full", NIGHTLY_MODES)
@pytest.mark.parametrize("scen", sorted(ARM_A_SCENARIOS))
def test_arm_a_discordance_pooled_level(scen: str, full: bool) -> None:
    lv = G.levels()
    stream = f"{DI}.fallback_a.{scen}"

    def rep(rng: np.random.Generator) -> np.ndarray:
        run = _disc_run(scen, _seed(rng))
        assert run.p.size >= M_POOLED, run.p.size
        p = rng.choice(run.p, M_POOLED, replace=False)
        return np.array([float((p <= a).mean()) for a in lv])

    gate = calibrate_levels(DI, "level", lv, rep, C.replicates(DI, full), stream=stream,
                            m=M_POOLED)  # fmt: skip
    print("\n".join(gate.detail()))
    gate.assert_ok()


@pytest.mark.calsuite(DI)
@pytest.mark.parametrize("full", NIGHTLY_MODES)
def test_arm_a_discordance_p_uniform(full: bool) -> None:
    stream = f"{DI}.fallback_a.operon_null"

    def rep(rng: np.random.Generator) -> float:
        return float(rng.choice(_disc_run("operon_null", _seed(rng)).p))

    rec = C.calibrate(DI, "uniform", 0.0, rep, C.replicates(DI, full), stream=stream)
    print(f"KS p = {rec['observed']} (R = {rec['R']}, status {rec['status']})")
    G.skip_without_verdict([rec])
    assert rec["status"] in {"pass", "unconfirmed"}, rec


@pytest.mark.calsuite(DI)
@pytest.mark.parametrize("full", NIGHTLY_MODES)
def test_arm_a_discordance_bh_family_fdr(full: bool) -> None:
    lv = G.levels()
    stream = f"{DI}.fallback_a.operon_null"
    R = C.replicates(DI, full)

    def rep(rng: np.random.Generator) -> np.ndarray:
        run = _disc_run("operon_null", _seed(rng))
        out = []
        for a in lv:
            v = int((run.q_null <= a).sum())
            out.append(v / max(v + int((run.q_alt <= a).sum()), 1))
        return np.array(out, dtype=float)

    gate = calibrate_levels(DI, "fdr", lv, rep, R, stream=stream)
    print("\n".join(gate.detail()))
    if C.owns_direct_checks():
        runs = [_disc_run("operon_null", C.suite_seed(stream, i)) for i in range(R)]
        alt = np.concatenate([r.q_alt for r in runs])
        print(f"boundary pairs tested per replicate {alt.size / R:.1f}, called at q 0.1:"
              f" {float((alt <= 0.1).mean()) if alt.size else float('nan'):.3f}")  # fmt: skip
    gate.assert_ok()


# =========================================================================== arm B: GSEA

GS = "gsea"


def _gsea_p(genes: pd.DataFrame, members: list[str], rng: np.random.Generator) -> float:
    from giap.config.model import GseaConfig
    from giap.function import gsea
    from giap.rng import for_step

    seed = int(rng.integers(0, 2**31 - 1))
    out = gsea.gsea_table(genes, [gene_sets("kegg", {"S": members})], (10, 500), GseaConfig(),
                          0.05, for_step("function.gsea", seed))  # fmt: skip
    assert len(out) == 1, out
    return float(out["p_value"].iloc[0])


def _gsea_world(rng: np.random.Generator, n: int, kind: str, rho: float, phi: float = 0.0
                ) -> tuple[np.ndarray, np.ndarray]:  # fmt: skip
    """``test_attack_gsea_structure._world``: true TU per gene and log2 abundance."""
    e = rng.standard_normal(n)
    t = e.copy()
    for i in range(1, n):
        t[i] = phi * t[i - 1] + np.sqrt(1 - phi**2) * e[i]
    sizes = rng.integers(1, 5, n) if kind == "1-4" else np.minimum(rng.geometric(0.48, n), 14)
    tu = np.repeat(np.arange(n), sizes)
    x = np.sqrt(rho) * t[tu] + np.sqrt(1 - rho) * rng.standard_normal(tu.size) + 20.0
    return tu, x


def _gsea_frame(blocks: np.ndarray, x: np.ndarray, imp: np.ndarray | None = None) -> pd.DataFrame:
    tus: dict[str, list[Gene]] = {}
    for g, (b, v) in enumerate(zip(blocks, x, strict=True)):
        tus.setdefault(f"TU_{int(b) + 1}", []).append(
            Gene(f"g{g + 1}", float(v), imputed=False if imp is None else bool(imp[g]))
        )
    return frame(tus)


def _first_genes(tu: np.ndarray, pick: np.ndarray) -> list[str]:
    first = {int(t): g for g, t in reversed(list(enumerate(tu)))}
    return [f"g{first[int(t)] + 1}" for t in pick]


Case = tuple[np.ndarray, np.ndarray, np.ndarray | None, list[str]]


def _gsea_case(arm: str, wrng: np.random.Generator) -> Case:
    """(true TU, abundance, imputed, set members) of one world of a gated GSEA arm."""
    if arm == "median_imputed":
        n = 600
        tu, x = _gsea_world(wrng, n, "1-4", 0.8)
        imp = wrng.random(tu.size) < 0.2
        x = x.copy()
        x[imp] = np.median(x[~imp])
        a = int(wrng.integers(0, n - 2 * SET_TUS))
        return tu, x, imp, _first_genes(tu, np.arange(a, n))  # cut by _contiguous
    if arm == "size_biased":
        tu, x = _gsea_world(wrng, 1500, "ecoli", 0.2)
        g = wrng.choice(tu.size, 40, replace=False)
        return tu, x, None, [f"g{int(i) + 1}" for i in g]
    if arm == "ar1_scattered":
        tu, x = _gsea_world(wrng, 1500, "ecoli", 0.8, phi=0.5)
        return tu, x, None, _first_genes(tu, wrng.choice(1500, 15, replace=False))
    raise ValueError(arm)


GSEA_ARMS = ("ar1_scattered", "median_imputed", "size_biased")
SET_TUS = 15


def _contiguous(members: list[str], blocks: np.ndarray, k: int = SET_TUS) -> list[str]:
    """The first genes of consecutive true TUs, in order, until they cover ``k`` analysis TUs:
    the contiguous set is "15 TUs" as the analysis sees them (merges would otherwise shrink it
    below the size filter of 10 and leave no row)."""
    seen: set[int] = set()
    out = []
    for m in members:
        b = int(blocks[int(m[1:]) - 1])
        if b not in seen and len(seen) == k:
            break
        seen.add(b)
        out.append(m)
    return out


@pytest.mark.calsuite(GS)
@pytest.mark.parametrize("full", NIGHTLY_MODES)
@pytest.mark.parametrize("source", SOURCES)
@pytest.mark.parametrize("arm", GSEA_ARMS)
def test_arm_b_gsea_at_level(arm: str, source: str, full: bool) -> None:
    R = C.replicates(GS, full)
    world = f"{GS}.fallback_b.{arm}"
    stream = f"{world}.{source}"
    index = _index(stream, R)

    def rep(rng: np.random.Generator) -> float:
        i = _rep_index(rng, index)
        tu, x, imp, members = _gsea_case(arm, C.replicate_rng(world, i))
        blocks = abstract_blocks(tu, _links_rng(world, i, source), source)
        if arm == "median_imputed":
            members = _contiguous(members, blocks)
        return _gsea_p(_gsea_frame(blocks, x, imp), members, rng)

    gate = calibrate_levels(GS, "level", C.load_tolerances().levels, rep, R, stream=stream)
    _assert_or_record(source, gate.failed, gate.detail(), f"gsea.{arm}")


# =========================================================================== arm B: ORA

ORA = "ora"
ORA_M = 10
ORA_SCENARIOS = {"4b_c1.5": 1.5, "4b_c2.5": 2.5}
ORA_KINDS = ("large", "single")
ORA_FLAG = "ora_uncalibrated"


def _reblock(sim: Any, blocks: np.ndarray, rng: np.random.Generator
             ) -> tuple[pd.DataFrame, dict[str, int]]:  # fmt: skip
    """The ORA gene frame with distance-model TUs (``tu_id``, ``tu_position``, folds)."""
    from giap.genetable.detectability import tu_folds

    genes = sim.genes.copy()
    ids = np.array([f"ATU_{int(b) + 1}" for b in blocks])
    genes["tu_id"] = ids
    genes["tu_position"] = genes.groupby("tu_id", sort=False).cumcount().to_numpy() + 1
    return genes, tu_folds(list(ids), K_FOLDS, rng)


@cache
def _ora_tables(scenario: str, source: str, i: int) -> dict[str, pd.DataFrame]:
    from giap.function import ora

    world = f"{ORA}.fallback_b.{scenario}"
    rng = C.replicate_rng(world, i)
    sim = simulate(rng, c=ORA_SCENARIOS[scenario])
    sets = {kind: whole_tu_sets(sim, rng, kind, ORA_M) for kind in ORA_KINDS}
    lrng = _links_rng(world, i, source)
    genes, folds = _reblock(sim, abstract_blocks(sim.tu_of_gene, lrng, source), lrng)
    out = {}
    for kind, s in sets.items():
        assert len(s) == ORA_M, (scenario, kind, len(s))
        out[kind] = ora.detection_enrichment(genes, [collection(s)], (5, 500), 0.05, folds=folds)
    return out


@pytest.mark.calsuite(ORA)
@pytest.mark.parametrize("full", NIGHTLY_MODES)
@pytest.mark.parametrize("source", SOURCES)
@pytest.mark.parametrize("scenario", sorted(ORA_SCENARIOS))
def test_arm_b_ora_tusize_at_level_or_flagged(scenario: str, source: str, full: bool) -> None:
    R = C.replicates(ORA, full)
    levels = list(C.load_tolerances().levels)
    failed: list[Any] = []
    detail: list[str] = []
    for kind in ORA_KINDS:
        for direction in ("over", "under"):
            stream = f"{ORA}.fallback_b.{scenario}.{source}.{kind}.{direction}"
            index = _index(stream, R)

            def rep(rng: np.random.Generator, k: str = kind, dr: str = direction,
                    ix: dict[int, int] = index) -> np.ndarray:  # fmt: skip
                tab = _ora_tables(scenario, source, _rep_index(rng, ix))[k]
                assert "tu_recal_skipped" not in tab.attrs["run_flags"]
                p = tab.loc[tab["direction"] == dr].sort_values("term")["p_value"].to_numpy(float)
                assert p.size == ORA_M and np.all(np.isfinite(p))
                return np.array([np.count_nonzero(p <= a) / ORA_M for a in levels])

            gate = calibrate_levels(ORA, "level", levels, rep, R, stream=stream, m=ORA_M)
            failed += gate.failed
            detail += gate.detail()
    flagged = all(ORA_FLAG in as_list(f) for t in _ora_tables(scenario, source, 0).values()
                  for f in t["flags"])  # fmt: skip
    _assert_or_record(source, failed and not flagged, detail, f"ora.{scenario}")


# =========================================================================== arm B: determinants

DET = "determinants"


def _det_world(rng: np.random.Generator, lrng: np.random.Generator, source: str, n: int = 1500
               ) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:  # fmt: skip
    """``test_attack_determinants_tu_null._world`` with analysis TUs: response and detection on
    the true TUs, the TU-constant feature on the analysis TUs (and, recorded, on the true TUs)."""
    sizes: list[int] = []
    while sum(sizes) < n:
        sizes.append(int(min(8, rng.geometric(0.45))))
    sizes[-1] -= sum(sizes) - n
    sizes = [x for x in sizes if x]
    k = len(sizes)
    true = np.repeat(np.arange(k), sizes)
    y = rng.normal(size=k)[true] + 0.7 * rng.normal(size=n)
    lat = rng.normal(size=k)[true] + 0.5 * rng.normal(size=n)
    x_true = rng.random(k)[true]
    blocks = abstract_blocks(true, lrng, source)
    x_an = lrng.random(int(blocks.max()) + 1)[blocks]
    tu = [f"ATU_{int(b):05d}" for b in blocks]
    return tu, x_an, x_true, y, lat > np.quantile(lat, 0.3)


@pytest.mark.calsuite(DET)
@pytest.mark.parametrize("full", NIGHTLY_MODES)
@pytest.mark.parametrize("source", SOURCES)
def test_arm_b_determinants_univariate_ci_level(source: str, full: bool) -> None:
    from giap.stats import bootstrap

    m = dm.mod()
    R = 5 * C.replicates(DET, full=full)
    world = f"{DET}.fallback_b.uni"
    miss: dict[str, list[bool]] = {"rho": [], "auroc": [], "rho_true_tu_x": [],
                                   "auroc_true_tu_x": []}  # fmt: skip
    flagged = True
    for i in range(R):
        rng = C.replicate_rng(world, i)
        tu, x_an, x_true, y, det = _det_world(rng, _links_rng(world, i, source), source)
        brng = C.replicate_rng(f"{world}.boot.{source}", i)
        f = dm.features({"promotech_tu": x_an.tolist()}, tu=tu, detected=det.tolist(),
                        abundance=np.where(det, y, np.nan).tolist())  # fmt: skip
        uni = m.univariate(f, 1000, brng)
        rho, auc = dm.one(uni, metric="spearman_rho"), dm.one(uni, metric="auroc")
        flagged &= all("ci_unestimated" in dm.flags_of(v) for v in uni["flags"])
        miss["rho"].append(not rho.ci_low <= 0.0 <= rho.ci_high)
        miss["auroc"].append(not auc.ci_low <= 0.5 <= auc.ci_high)
        tus = np.array(tu)
        gr = pd.DataFrame({"tu": tus[det], "x": x_true[det], "y": y[det]})
        ga = pd.DataFrame({"tu": tus, "x": x_true, "y": det})
        cr = bootstrap.tu_block_bootstrap(gr, "tu", m._rho, 1000, brng, columns=("x", "y"))
        ca = bootstrap.tu_block_bootstrap(ga, "tu", m._auroc, 1000, brng, columns=("x", "y"))
        miss["rho_true_tu_x"].append(not cr.ci_low <= 0.0 <= cr.ci_high)
        miss["auroc_true_tu_x"].append(not ca.ci_low <= 0.5 <= ca.ci_high)
    rates = {k: float(np.mean(v)) for k, v in miss.items()}
    bound = C.level_bound(0.05, R)
    _log("arm_b_determinants", {"source": source, "R": R, "rates": rates, "bound": bound,
                                "flagged": flagged})  # fmt: skip
    print(f"[{source}] CI miss rates, R={R}: {rates} (bound {bound:.3f}, flagged {flagged})")
    over = {k: v for k, v in rates.items() if k in ("rho", "auroc") and v > bound}
    if flagged:
        if over:
            warnings.warn(f"[{source}] flagged univariate CIs miss at {over}", stacklevel=1)
        return
    _assert_or_record(source, list(over), [f"miss {over} bound {bound:.4f}"], "determinants.uni")


# =========================================================================== arm B: autocorrelation

AC = "autocorrelation"
AC_GENES = (1000, 4500)
AC_LAGS = (1, 20)


def genome_blocks(genes: pd.DataFrame, rng: np.random.Generator, source: str) -> np.ndarray:
    """Distance-model TU ids of a giapsim gene table (per contig, genome order; the wrap pair of
    a circular replicon is left unlinked)."""
    out = np.empty(len(genes), dtype=object)
    for contig, sub in genes.groupby("contig", sort=False):
        sub = sub.sort_values("start", kind="stable")
        strand = sub["strand"].astype(str).to_numpy()
        true = pd.factorize(sub["tu_id"].astype(str))[0]
        b = analysis_tus(true, strand[1:] == strand[:-1], rng, source)
        out[sub.index.to_numpy()] = [f"{contig}_ATU_{int(v)}" for v in b]
    return out


@cache
def _ac_rows(n: int, source: str, i: int) -> pd.DataFrame:
    sc = acal.tu_null(n, name=f"fallback_b_tu_null_{n}")
    rng = C.replicate_rng(sc.stream, i)
    tables = acal.simulate(sc, rng)
    seed = int(rng.integers(0, 2**31 - 1))
    genes = tables.genes.reset_index(drop=True).copy()
    genes["tu_id"] = genome_blocks(genes, _links_rng(sc.stream, i, source), source)
    cfg = ad.config(acal.PERMUTATIONS, lags=(acal.LAGS[0], acal.LAGS[-1]))
    return ad.rows(genes, tables.replicons, cfg, seed)


@pytest.mark.calsuite(AC)
@pytest.mark.parametrize("full", NIGHTLY_MODES)
@pytest.mark.parametrize("source", SOURCES)
def test_arm_b_autocorrelation_tu_block_at_level(source: str, full: bool) -> None:
    R = C.replicates(AC, full)
    n = AC_GENES[1] if full else AC_GENES[0]
    levels = list(C.load_tolerances().levels)
    failed: list[Any] = []
    detail: list[str] = []
    for statistic in ad.STATISTICS:
        for lag in AC_LAGS:
            stream = f"{AC}.fallback_b.tu_null_{n}.{source}.{statistic}.tu_block.lag{lag}"
            index = _index(stream, R)

            def rep(rng: np.random.Generator, s: str = statistic, k: int = lag,
                    ix: dict[int, int] = index) -> float:  # fmt: skip
                return acal.lag_p(_ac_rows(n, source, _rep_index(rng, ix)), s, "tu_block", k)

            gate = calibrate_levels(AC, "level", levels, rep, R, stream=stream)
            failed += gate.failed
            detail += gate.detail()
    _assert_or_record(source, failed, detail, f"autocorrelation.tu_null_{n}")


def test_genome_blocks_split_and_merge_true_tus() -> None:
    """The autocorrelation arm's blocks differ from the true TUs in both directions."""
    sc = acal.tu_null(600, name="fallback_b_probe")
    tables = acal.simulate(sc, np.random.default_rng(3))
    genes = tables.genes.reset_index(drop=True)
    blocks = genome_blocks(genes, np.random.default_rng(4), "distance_model")
    pair = pd.DataFrame({"t": genes["tu_id"].astype(str), "b": blocks})
    assert (pair.groupby("t")["b"].nunique() > 1).any()  # a split
    assert (pair.groupby("b")["t"].nunique() > 1).any()  # a merge
