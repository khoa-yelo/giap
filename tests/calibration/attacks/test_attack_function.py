"""M8 gate, half 2: attacks 2 and 3 (docs/specs/M8-15-calibration-gate.md § Suites and criteria,
lead decisions 3-7; GIAP_SPEC.md §7.7.5, §8.4, §8.5 attacks 2 and 3; docs/methods/function.md
"Simulator worlds for the M8-15 gate"; ADR-0015, ADR-0029, ADR-0035, ADR-0038).
Statistician-owned nulls, written by the test-author as specified.

**Attack 2** (``ora`` suite, R = ``replicates_fast/full.ora``), two arms (lead decision 4):
``ora_attack_confounded(4000, -1.5)`` (separation) and ``ora_attack_confounded(4000, -0.75)``
(overlap). The confounded GO term (``sets.confounded``) is built through M8-03 from the world's
gene table and tested alone (its collection name kept; one term, so ``p_value = p_noelim``).

* *fitted* (gated): p_detect fitted by M4-13 with the default ``DetectabilityConfig()`` on
  covariates ``covariate_table`` recomputes from the sequences (never ``n_obs_true``), the M4-13
  folds to the TU recalibration; §8.4 level of the over **and** the under p at every alpha.
  Expected strongly conservative at -1.5 (decision 5: plug-in variance); a level gate only.
* *teeth* (gated, decision 3): the unweighted hypergeometric (omega forced to 1) on the same
  term, **gene-level** (M = eligible genes, n = members, N = detected, k = detected members):
  the under-representation p ``hypergeom(M, n, N).cdf(k)`` must be a confirmed ``fail`` through
  ``calibrate`` (both blocks) **and** reject in > 2 x alpha at alpha = 0.05 in each block. The
  -0.75 arm must also be a confirmed fail at R >= 100 (decision 4; statistician probe 0.95 —
  a gene-level number, as the M8-14 teeth ``test_example2_teeth_unweighted_test_rejects``).
  The TU-level omega = 1 variant (GIAP's TU units, every q equal) is **recorded** beside it: the
  test-author pre-run found it barely moves (replicate 0: p 0.021 at -1.5, 0.53 at -0.75),
  because members are random genes, so any-gene TU membership dilutes the confounding.
* *oracle* (recorded, decision 6): p_detect = p_true through the same ORA engine; the rate and
  the two-sided §8.4 band per alpha go to ``function_gate_diagnostics.jsonl``; below the lower
  bound (probes 0.025-0.067 at alpha 0.10) is reported with a warning, never gated.

**Attack 3** (``gsea`` suite; ADR-0038 amendment: R fgsea only, ``slow`` + ``requires_r``):
``gsea_attack_operon_packed(4000, 20)``; the 20 packed sets go through M8-03 as the world's
KEGG collection.

* *TU-collapsed* (gated): ``gsea_table(..., engine=RFgseaEngine)``; the 20 p of a replicate
  pooled (``check_level_pooled``, m = 20) at every alpha.
* *gene-level teeth* (gated): the same sets and genes ranked individually (each gene its own
  "TU", so the R engine permutes gene labels as if independent): a confirmed ``fail`` and a rate
  > 2 x alpha at 0.05. The spec table lists this arm as "recorded" at R 20 / 100; the lead's
  brief gates it as a confirmed fail, and R is the ``gsea`` suite's (R3: no R outside
  tolerances.yaml).

Nothing here needs M8-11 (#193): ADR-0038 routes every GSEA level claim to ``RFgseaEngine``
(``fgseaMultilevel`` in R), and ORA has no multilevel stage.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels
from helpers import function_gate_chain as F
from scipy.stats import hypergeom

ORA, GSEA = "ora", "gsea"
ORA_SIZE = (5, 500)
GSEA_SIZE = (10, 500)
N_GENES = 4000
ALPHA_TEETH = 0.05
SHIFTS = {"sep": -1.5, "overlap": -0.75}  # lead decision 4
N_PACKED = 20

NIGHTLY = pytest.mark.calibration_nightly  # fast tests over 60 s
MODES = [
    pytest.param(False, marks=[pytest.mark.calibration, NIGHTLY], id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
MODES_R = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
NO_VERDICT = frozenset({C.SHARD_STATUS, C.PENDING_STATUS, C.SKIPPED_STATUS})


def _levels() -> list[float]:
    return list(C.load_tolerances().levels)


def _sc() -> Any:
    import giapsim.scenarios as sc

    return sc


def _draw(pair: Any, rng: np.random.Generator) -> Any:
    import giapsim.function_sim as fs

    return fs.draw(pair, rng)


def _sharded() -> bool:
    import os

    return any(os.environ.get(v) for v in (C.SHARD_ENV, C.CONFIRM_ENV))


def _teeth_blocks(log: Path, suite: str, stream: str) -> list[dict[str, Any]]:
    recs = [json.loads(x) for x in (log / f"{suite}.jsonl").read_text().splitlines() if x]
    return [r for r in recs if r.get("stream") == stream]


# --------------------------------------------------------------------------- attack 2


def _sim2(tag: str) -> str:
    return f"ora.m815.attack2.{tag}.sim"


_A2: dict[tuple[str, int], dict[str, Any]] = {}


def _a2_world(tag: str, i: int) -> tuple[Any, pd.DataFrame, Any]:
    w = _draw(_sc().ora_attack_confounded(N_GENES, SHIFTS[tag]), C.replicate_rng(_sim2(tag), i))
    genes = F.gene_frame(w)
    one = F.single_term(F.collections(w, genes), F.confounded_term(w))
    return w, genes, one


def _a2_cheap(tag: str, i: int, w: Any, genes: pd.DataFrame, one: Any) -> dict[str, Any]:
    """Oracle p (p_true through the ORA engine) and the unweighted TU hypergeometric p."""
    from giap.function import ora
    from giap.genetable.detectability import tu_folds

    folds = tu_folds(list(genes["tu_id"]), 5, C.replicate_rng(_sim2(tag) + ".oracle", i))
    tab = ora.detection_enrichment(F.with_p_detect(genes, F.true_p_detect(w, genes)), [one],
                                   ORA_SIZE, 0.05, folds=folds)  # fmt: skip
    assert len(tab) == 2, tab
    by = tab.set_index("direction")
    mu, n, nd, k = (int(by.loc["under", c]) for c in
                    ("n_universe", "n_in_set", "n_detected", "k_detected_in_set"))  # fmt: skip
    el = ~genes["partial"].to_numpy(bool)
    det = genes["detected"].to_numpy(bool)[el]
    member = genes["locus_tag"].astype(str).isin(next(iter(one.sets.values()))).to_numpy()[el]
    gm, gn, gnd, gk = int(el.sum()), int(member.sum()), int(det.sum()), int((det & member).sum())
    return {
        "oracle_over": float(by.loc["over", "p_value"]),
        "oracle_under": float(by.loc["under", "p_value"]),
        "teeth_under": float(min(1.0, hypergeom(gm, gn, gnd).cdf(gk))),
        "teeth_over": float(min(1.0, hypergeom(gm, gn, gnd).sf(gk - 1))),
        "teeth_tu_under": float(min(1.0, hypergeom(mu, n, nd).cdf(k))),
        "n_in_set_tus": n, "k_tus": k, "n_in_set_genes": gn, "k_genes": gk,
        "expected_genes_unweighted": gn * gnd / gm,
    }  # fmt: skip


def _a2(tag: str, i: int, fitted: bool) -> dict[str, Any]:
    key = (tag, i)
    rec = _A2.get(key)
    if rec is not None and (not fitted or "fitted_over" in rec):
        return rec
    from giap.function import ora

    w, genes, one = _a2_world(tag, i)
    rec = rec or _a2_cheap(tag, i, w, genes, one)
    if fitted:
        p, folds = F.fitted_p_detect(w, genes, C.replicate_rng(_sim2(tag) + ".fit", i))
        tab = ora.detection_enrichment(F.with_p_detect(genes, p), [one], ORA_SIZE, 0.05,
                                       folds=folds)  # fmt: skip
        assert "tu_recal_skipped" not in tab.attrs["run_flags"], tab.attrs["run_flags"]
        by = tab.set_index("direction")
        rec = {**rec, "fitted_over": float(by.loc["over", "p_value"]),
               "fitted_under": float(by.loc["under", "p_value"]),
               "odds_weight": float(by.loc["under", "odds_weight"])}  # fmt: skip
    _A2[key] = rec
    return rec


@pytest.mark.calsuite(ORA)
@pytest.mark.parametrize("tag", sorted(SHIFTS))
@pytest.mark.parametrize("full", MODES)
def test_attack2_fitted_p_detect_confounded_term_at_level(full: bool, tag: str) -> None:
    """Gate: the weighted ORA with the fitted p_detect is at level on the confounded term, over
    and under, at every alpha (decision 4; conservative expected, decision 5)."""
    R = C.replicates(ORA, full)
    levels = _levels()
    ran: set[int] = set()
    gates = []
    for d in ("over", "under"):
        stream = f"ora.m815.attack2.{tag}.fitted.{d}{'.full' if full else ''}"

        def one(i: int, d: str = d) -> float:
            ran.add(i)
            return _a2(tag, i, True)[f"fitted_{d}"]

        gates.append(calibrate_levels(ORA, "level", levels, F.by_index(stream, R, one), R,
                                      stream=stream))  # fmt: skip
    if ran:
        recs = [_A2[(tag, i)] for i in sorted(ran)]
        F.diagnostics({
            "check": f"attack2.{tag}.fitted", "shift": SHIFTS[tag], "full": full,
            "R_run": len(recs),
            **{d: {"rates": F.rate_record([r[f"fitted_{d}"] for r in recs], levels),
                   "ks": F.ks_record([r[f"fitted_{d}"] for r in recs])} for d in ("over", "under")},
            "mean_odds_weight": float(np.mean([r["odds_weight"] for r in recs])),
        })  # fmt: skip
    for g in gates:
        g.assert_ok()


@pytest.mark.calibration
@NIGHTLY
@pytest.mark.calsuite(ORA)
@pytest.mark.parametrize("tag", sorted(SHIFTS))
def test_teeth_attack2_unweighted_under_rejects_beyond_twice_alpha(
    tag: str, teeth_log: Path
) -> None:
    """Gated teeth (decision 3): omega forced to 1 is a confirmed fail, > 2 x alpha per block."""
    R = C.replicates(ORA, False)
    assert R >= 100  # decision 4: the overlap arm must fail at R >= 100
    stream = f"ora.m815.attack2.{tag}.teeth_unweighted_under"
    rec = C.calibrate(ORA, "level", ALPHA_TEETH,
                      F.by_index(stream, R, lambda i: _a2(tag, i, False)["teeth_under"]), R,
                      stream=stream, log_dir=teeth_log)  # fmt: skip
    blocks = _teeth_blocks(teeth_log, ORA, stream)
    tu_p = np.array([_A2[(tag, i)]["teeth_tu_under"] for i in range(2 * R) if (tag, i) in _A2])
    F.diagnostics({"check": f"attack2.{tag}.teeth", "shift": SHIFTS[tag], "R": R,
                   "status": rec["status"],
                   "observed_by_block": {r["block"]: r["observed"] for r in blocks},
                   "tu_level_omega1_under": F.rate_record(tu_p, [ALPHA_TEETH])})  # fmt: skip
    assert rec["status"] == "fail", f"the unweighted arm is not a confirmed fail: {rec}"
    assert len(blocks) == 2 and all(r["observed"] > 2 * ALPHA_TEETH for r in blocks), blocks


@pytest.mark.calibration
@NIGHTLY
@pytest.mark.calsuite(ORA)
@pytest.mark.parametrize("tag", sorted(SHIFTS))
def test_attack2_oracle_arm_level_recorded(tag: str) -> None:
    """Recorded (decision 6): the oracle arm's rate and §8.4 band per alpha; below the lower
    bound is reported (warning), never gated."""
    if _sharded():
        pytest.skip("recorded arm: not run under --shard / confirm runs")
    R = C.replicates(ORA, False)
    levels = _levels()
    recs = [_a2(tag, i, False) for i in range(R)]
    out: dict[str, Any] = {"check": f"attack2.{tag}.oracle", "shift": SHIFTS[tag], "R": R}
    for d in ("over", "under"):
        p = np.array([r[f"oracle_{d}"] for r in recs])
        assert np.all(np.isfinite(p)) and np.all((p >= 0) & (p <= 1)), p
        out[d] = {"rates": F.rate_record(p, levels), "ks": F.ks_record(p)}
        for a, cell in out[d]["rates"].items():
            if cell["below_lower"]:
                warnings.warn(f"attack 2 {tag} oracle {d}: rate {cell['rate']:.3f} below the "
                              f"§8.4 lower bound {cell['lower']:.3f} at alpha {a} (recorded, "
                              "decision 6)", stacklevel=1)  # fmt: skip
    F.diagnostics(out)


# --------------------------------------------------------------------------- attack 3


@pytest.fixture(scope="module")
def r_engine(tmp_path_factory: pytest.TempPathFactory) -> Any:
    return F.r_engine(Path(tmp_path_factory.mktemp("m815_fgsea_a3")))


SIM3 = "gsea.m815.attack3.sim"
_A3: dict[int, dict[str, np.ndarray]] = {}


def _a3(engine: Any, i: int) -> dict[str, np.ndarray]:
    """TU-collapsed and gene-level R-engine p of the 20 packed sets of replicate ``i``."""
    if i in _A3:
        return _A3[i]
    from giap.config.model import GseaConfig
    from giap.function import gsea
    from giap.rng import for_step

    w = _draw(_sc().gsea_attack_operon_packed(N_GENES, N_PACKED), C.replicate_rng(SIM3, i))
    genes = F.gene_frame(w)
    genes = F.with_p_detect(genes, F.true_p_detect(w, genes))
    kegg = [c for c in F.collections(w, genes) if c.collection == "kegg"]
    assert len(kegg) == 1 and len(kegg[0].sets) == N_PACKED, [len(c.sets) for c in kegg]
    rng = C.replicate_rng(SIM3 + ".engine", i)
    out: dict[str, np.ndarray] = {}
    gene_level = genes.assign(tu_id=genes["locus_tag"], tu_position=1)
    for name, frame in (("tu", genes), ("gene", gene_level)):
        tab = gsea.gsea_table(frame, kegg, GSEA_SIZE, GseaConfig(), 0.05,
                              for_step("function.gsea", int(rng.integers(0, 2**31 - 1))),
                              engine=engine)  # fmt: skip
        assert set(tab["method"]) == {"fgsea_multilevel"}, set(tab["method"])
        p = tab["p_value"].to_numpy(float)
        assert p.size == N_PACKED and np.all(np.isfinite(p)), (name, tab)
        out[name] = p
    _A3[i] = out
    return out


@pytest.mark.calsuite(GSEA)
@pytest.mark.slow
@pytest.mark.requires_r
@pytest.mark.parametrize("full", MODES_R)
def test_attack3_tu_collapsed_r_engine_pooled_level(full: bool, r_engine: Any) -> None:
    R = C.replicates(GSEA, full)
    levels = _levels()
    stream = f"gsea.m815.attack3.tu{'.full' if full else ''}"
    ran: set[int] = set()

    def one(i: int) -> np.ndarray:
        ran.add(i)
        return F.shares(_a3(r_engine, i)["tu"], levels)

    gate = calibrate_levels(GSEA, "level", levels, F.by_index(stream, R, one), R, stream=stream,
                            m=N_PACKED)  # fmt: skip
    if ran:
        p = np.concatenate([_A3[i]["tu"] for i in sorted(ran)])
        F.diagnostics({"check": "attack3.tu.r", "full": full, "R_run": len(ran), "m": N_PACKED,
                       "rates": F.rate_record(p, levels), "ks": F.ks_record(p)})  # fmt: skip
    gate.assert_ok()


@pytest.mark.calibration
@pytest.mark.calsuite(GSEA)
@pytest.mark.slow
@pytest.mark.requires_r
def test_teeth_attack3_gene_level_r_engine_confirmed_fail(r_engine: Any, teeth_log: Path) -> None:
    R = C.replicates(GSEA, False)
    stream = "gsea.m815.attack3.teeth_gene_level"
    fn = F.by_index(stream, R, lambda i: F.shares(_a3(r_engine, i)["gene"], [ALPHA_TEETH])[0])
    rec = C.calibrate(GSEA, "level", ALPHA_TEETH, fn, R, stream=stream, log_dir=teeth_log,
                      m=N_PACKED)  # fmt: skip
    blocks = _teeth_blocks(teeth_log, GSEA, stream)
    F.diagnostics({"check": "attack3.gene.teeth.r", "R": R, "status": rec["status"],
                   "observed_by_block": {r["block"]: r["observed"] for r in blocks}})  # fmt: skip
    assert rec["status"] == "fail", f"the gene-level arm is not a confirmed fail: {rec}"
    assert len(blocks) == 2 and all(r["observed"] > 2 * ALPHA_TEETH for r in blocks), blocks
