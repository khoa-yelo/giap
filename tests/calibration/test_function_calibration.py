"""M8 gate, half 2: the ``ora`` and ``gsea`` nulls at level (docs/specs/M8-15-calibration-gate.md
§ Suites and criteria, lead decisions 1-7; GIAP_SPEC.md §7.7.5, §8.4 rows ORA and GSEA, §13 row
M8; ADR-0015, ADR-0029, ADR-0035, ADR-0038). Attacks 2 and 3: ``attacks/test_attack_function.py``.

Criterion: the §8.4 **level** check at every alpha of ``levels`` (lead decision 2: Wallenius /
CPB p are discrete, fgsea p are permutation-based, so "uniform" is gated as level); the KS
statistic is computed and **recorded** (``function_gate_diagnostics.jsonl``), never gated.
Through ``_gate.calibrate_levels`` -> ``_criteria.calibrate`` (two-block confirmation, shards,
jsonl). R = ``replicates_fast/full.<suite>`` (protected tolerances.yaml: ora 100 / 1,000, gsea
50 / 300).

* ``ora`` (gate): ``ora_null(4000)``; the gene table of the world (``helpers.function_gate_chain``),
  the six M8-03 collections of the synthetic ontology (``go_slim`` is a synthetic slim, M8-14
  gap 4), p_detect **fitted** by M4-13 with the default ``DetectabilityConfig()`` on covariates
  recomputed from the sequences by ``covariate_table``, the M4-13 folds passed to the TU
  recalibration. ``detection_enrichment`` over all six collections; per replicate **one** tested
  term per collection (m = ``M_ORA`` = 6; lead decision 9 B2), drawn at random (own stream) in
  the fixed order ``ORA_ORDER``, ``go_slim`` excluding term IDs already picked; their
  ``p_value`` pooled (``check_level_pooled``), **per direction** (over, under). Recorded: the
  oracle arm (p_true) on the same terms, ``p_noelim``, per-collection rates over every tested
  term by TU-size band (``SIZE_BANDS``), KS.
* ``gsea`` (gate; ADR-0038 amendment: Slurm only, ``slow`` + ``requires_r``): ``gsea_null(4000)``,
  the six collections, ``gsea_table(..., engine=RFgseaEngine)``; ``M_GSEA`` = 100 distinct
  tested terms per replicate pooled (spec: ">= 100 sets per replicate"). KS recorded.

Fast-tier timings (test-author pre-run, one core, iron-01): ORA ~5.5 s per replicate (fit 3.5 s),
so ~9 min for R = 100 (``calibration_nightly``: > 60 s); GSEA ~6 s per replicate with the R call.
"""

from __future__ import annotations

import os
from collections.abc import Callable
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels
from helpers import function_gate_chain as F

ORA, GSEA = "ora", "gsea"
ORA_SIZE = (5, 500)  # function.ora_size
GSEA_SIZE = (10, 500)  # function.gsea_size
ORA_ORDER = ("cog", "kegg", "go_bp", "go_mf", "go_cc", "go_slim")  # B2: fixed pick order
M_ORA = len(ORA_ORDER)
SIZE_BANDS = ((5, 9), (10, 49), (50, 500))  # n_in_set (TUs), as statistician probe 5
M_GSEA = 100
N_GENES = 4000
DIRECTIONS = ("over", "under")
SIM_ORA = "ora.m815.null.sim"
SIM_GSEA = "gsea.m815.null.sim"

MODES_ORA = [
    pytest.param(False, marks=[pytest.mark.calibration, pytest.mark.calibration_nightly],
                 id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]  # fmt: skip
MODES_GSEA = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]


def _levels() -> list[float]:
    return list(C.load_tolerances().levels)


def _sc() -> Any:
    import giapsim.scenarios as sc

    return sc


def _draw(pair: Any, rng: np.random.Generator) -> Any:
    import giapsim.function_sim as fs

    return fs.draw(pair, rng)


def by_index(stream: str, R: int, fn: Callable[[int], np.ndarray]) -> Callable[[Any], np.ndarray]:
    """``calibrate``'s ``fn(rng)`` for a function of the replicate index (one shared world per
    index feeds every check; the stream only names the check)."""
    index = {C.suite_seed(stream, i): i for i in range(2 * R)}
    return lambda rng: fn(index[int(rng.bit_generator.seed_seq.entropy)])


def shares(p: np.ndarray, levels: list[float]) -> np.ndarray:
    assert p.ndim == 1 and p.size > 0 and np.all(np.isfinite(p)), p
    return np.array([np.count_nonzero(p <= a) / p.size for a in levels])


def _distinct_terms(tab: pd.DataFrame) -> list[str]:
    """Tested terms with a finite p, one per term ID (``go_slim`` repeats GO terms)."""
    ok = tab[np.isfinite(tab["p_value"].to_numpy(float))]
    return sorted(set(ok["term"].astype(str)))


# --------------------------------------------------------------------------- ora


_ORA: dict[int, dict[Any, Any]] = {}


def _ora_null(i: int) -> dict[Any, Any]:
    """Replicate ``i`` (memoised): fitted and oracle ORA p of ``M_ORA`` random tested terms, per
    direction; the world, fit, oracle folds and pick each draw from their own stream."""
    if i not in _ORA:
        _ORA[i] = _ora_null_run(i)
    return _ORA[i]


def _ora_null_run(i: int) -> dict[Any, Any]:
    from giap.function import ora
    from giap.genetable.detectability import tu_folds

    w = _draw(_sc().ora_null(N_GENES), C.replicate_rng(SIM_ORA, i))
    genes = F.gene_frame(w)
    p_fit, folds = F.fitted_p_detect(w, genes, C.replicate_rng(SIM_ORA + ".fit", i))
    colls = F.collections(w, genes)
    fitted = ora.detection_enrichment(F.with_p_detect(genes, p_fit), colls, ORA_SIZE, 0.05,
                                      folds=folds)  # fmt: skip
    assert "tu_recal_skipped" not in fitted.attrs["run_flags"], fitted.attrs["run_flags"]
    o_rng = C.replicate_rng(SIM_ORA + ".oracle", i)
    oracle = ora.detection_enrichment(
        F.with_p_detect(genes, F.true_p_detect(w, genes)), colls, ORA_SIZE, 0.05,
        folds=tu_folds(list(genes["tu_id"]), 5, o_rng),
    )  # fmt: skip
    pick = pick_one_per_collection(fitted, C.replicate_rng(SIM_ORA + ".pick", i))
    out: dict[Any, Any] = {"collections": sorted({c.collection for c in colls}), "pick": pick}
    for d in DIRECTIONS:
        for name, tab in (("fitted", fitted), ("oracle", oracle)):
            sub = tab[tab["direction"] == d].set_index(["collection", "term"])
            out[(name, d)] = sub.loc[pick, "p_value"].to_numpy(float)
            out[(name, d, "noelim")] = sub.loc[pick, "p_noelim"].to_numpy(float)
        sub = fitted[fitted["direction"] == d]
        out[("by_collection", d)] = {
            str(c): (g["p_value"].to_numpy(float), g["n_in_set"].to_numpy(int))
            for c, g in sub.groupby("collection")
        }
    return out


def pick_one_per_collection(tab: pd.DataFrame, rng: np.random.Generator) -> list[tuple[str, str]]:
    """B2: one tested term (finite p) per collection, in ``ORA_ORDER``; a term ID picked in an
    earlier collection is excluded later (``go_slim`` repeats GO IDs)."""
    ok = tab[np.isfinite(tab["p_value"].to_numpy(float))]
    picked: list[tuple[str, str]] = []
    seen: set[str] = set()
    for c in ORA_ORDER:
        cand = sorted(set(ok.loc[ok["collection"] == c, "term"].astype(str)) - seen)
        assert cand, f"no tested term left in {c}"
        t = str(cand[int(rng.integers(0, len(cand)))])
        picked.append((c, t))
        seen.add(t)
    return picked


def _banded(p: np.ndarray, n: np.ndarray, levels: list[float]) -> dict[str, Any]:
    out: dict[str, Any] = {"all": F.rate_record(p, levels)}
    for lo, hi in SIZE_BANDS:
        sel = (n >= lo) & (n <= hi)
        out[f"{lo}-{hi}"] = {"n": int(sel.sum()), **F.rate_record(p[sel], levels)}
    return out


@pytest.mark.calsuite(ORA)
@pytest.mark.parametrize("full", MODES_ORA)
def test_ora_null_fitted_p_detect_pooled_level_per_direction(full: bool) -> None:
    R = C.replicates(ORA, full)
    levels = _levels()
    gates = []
    ran: set[int] = set()
    F.prefetch(_ORA, _ora_null_run, F.block1_indices(R))  # wall time only (pure in i)

    def one(i: int, d: str) -> np.ndarray:
        ran.add(i)
        return shares(_ora_null(i)[("fitted", d)], levels)

    for d in DIRECTIONS:
        stream = f"ora.m815.null.fitted.{d}{'.full' if full else ''}"
        fn = by_index(stream, R, lambda i, d=d: one(i, d))
        gates.append(calibrate_levels(ORA, "level", levels, fn, R, stream=stream, m=M_ORA))
    if ran:
        recs = [_ora_null(i) for i in sorted(ran)]
        rec: dict[str, Any] = {"check": "ora.null", "full": full, "R_run": len(ran), "m": M_ORA,
                               "collections": recs[0]["collections"]}  # fmt: skip
        for d in DIRECTIONS:
            fit = np.concatenate([r[("fitted", d)] for r in recs])
            orc = np.concatenate([r[("oracle", d)] for r in recs])
            rec[d] = {
                "fitted": F.rate_record(fit, levels), "fitted_ks": F.ks_record(fit),
                "fitted_noelim": F.rate_record(
                    np.concatenate([r[("fitted", d, "noelim")] for r in recs]), levels),
                "oracle": F.rate_record(orc, levels), "oracle_ks": F.ks_record(orc),
                "by_collection": {
                    c: _banded(
                        np.concatenate([r[("by_collection", d)].get(c, ([], []))[0]
                                        for r in recs]).astype(float),
                        np.concatenate([r[("by_collection", d)].get(c, ([], []))[1]
                                        for r in recs]).astype(int),
                        levels)
                    for c in recs[0]["collections"]
                },
            }  # fmt: skip
        F.diagnostics(rec)
    for g in gates:
        g.assert_ok()


# --------------------------------------------------------------------------- gsea


@pytest.fixture(scope="module")
def r_engine(tmp_path_factory: pytest.TempPathFactory) -> Any:
    from giap.io.fgsea_r import r_fgsea_engine, resolve_launcher

    return r_fgsea_engine(
        resolve_launcher(os.environ), tmp_root=Path(tmp_path_factory.mktemp("m815_fgsea"))
    )


def _gsea_null(engine: Any, i: int) -> np.ndarray:
    from giap.config.model import GseaConfig
    from giap.function import gsea
    from giap.rng import for_step

    w = _draw(_sc().gsea_null(N_GENES), C.replicate_rng(SIM_GSEA, i))
    genes = F.gene_frame(w)
    genes = F.with_p_detect(genes, F.true_p_detect(w, genes))
    colls = F.collections(w, genes)
    rng = C.replicate_rng(SIM_GSEA + ".engine", i)
    tab = gsea.gsea_table(genes, colls, GSEA_SIZE, GseaConfig(), 0.05,
                          for_step("function.gsea", int(rng.integers(0, 2**31 - 1))),
                          engine=engine)  # fmt: skip
    assert set(tab["method"]) == {"fgsea_multilevel"}, set(tab["method"])
    terms = _distinct_terms(tab)
    assert len(terms) >= M_GSEA, f"only {len(terms)} tested sets (spec: >= 100 per replicate)"
    pick = sorted(C.replicate_rng(SIM_GSEA + ".pick", i).choice(terms, M_GSEA, replace=False))
    p = tab.drop_duplicates("term").set_index("term").loc[pick, "p_value"].to_numpy(float)
    assert np.all((p > 0) & (p <= 1)), p
    return p


@pytest.mark.calsuite(GSEA)
@pytest.mark.slow
@pytest.mark.requires_r
@pytest.mark.parametrize("full", MODES_GSEA)
def test_gsea_null_r_engine_pooled_level(full: bool, r_engine: Any) -> None:
    R = C.replicates(GSEA, full)
    levels = _levels()
    stream = f"gsea.m815.null.r{'.full' if full else ''}"
    seen: dict[int, np.ndarray] = {}

    def one(i: int) -> np.ndarray:
        if i not in seen:
            seen[i] = _gsea_null(r_engine, i)
        return shares(seen[i], levels)

    gate = calibrate_levels(GSEA, "level", levels, by_index(stream, R, one), R, stream=stream,
                            m=M_GSEA)  # fmt: skip
    if seen:
        p = np.concatenate(list(seen.values()))
        F.diagnostics({"check": "gsea.null.r", "full": full, "R_run": len(seen), "m": M_GSEA,
                       "rates": F.rate_record(p, levels), "ks": F.ks_record(p),
                       "block_seconds": gate.block_seconds})  # fmt: skip
    gate.assert_ok()
