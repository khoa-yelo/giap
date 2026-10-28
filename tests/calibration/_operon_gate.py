"""Shared helpers of the M6-09 operon gate (docs/specs/M6-09-calibration-gate.md; GIPPA_SPEC.md
§7.4.7 "Calibration", §8.1 rows "Operon discordance" and "Missed operon", §8.4, §8.5).

Replicate ``i`` of stream ``s`` draws one genome with ``gippasim.simulate_operons(cfgs,
suite_seed(s, i))`` and runs the merged steps on its pair table:

- missed operons: ``gippa.operons.missed.missed_operons(pairs, cfg)``; a call is a scored NO row
  with ``q_value <= q``. A call is **true** when the pair is a planted co-abundant pair
  (``truth["missed_operons"]``) or a UniOP false negative (``true_operon``: the two genes share a
  TU); every other call is false. The §8.1 family is "all scored NO pairs of the run" at q, one
  family per run, as the step forms it.
- discordance: ``gippa.operons.discordance.discordance(pairs, genes, cfg, p_miss)``; the §8.1
  family is every OP row with ``test = both_detected_t`` (one BH family per run; P(miss) stays
  out of it, M6-index decision 1).

Pair classes are re-derived from ``uniop_p`` with the calibrated UniOP threshold in
``gippa/data/calibration.yaml`` (ADR-0022, M6-12c; 0.5 by Khoa's decision 2026-10-24), so the gate
uses the cut-off the pipeline ships. Pure apart from the module-level replicate caches.
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd

#: genome size of every M6-09 suite (M6-03: >= 1,000 genes; ~490 OP, 300 NO, 150 CV, 80 SW pairs
#: with both genes detected)
N_GENES = 1800
#: the shipped decoy sets, and the convergent-only configuration of the non-vacuous arms
DEFAULT_DECOYS: tuple[str, ...] = ("convergent", "sandwich")
CV_ONLY: tuple[str, ...] = ("convergent",)
#: planted pairs of the non-vacuous FDR arms: 100 co-abundant NO pairs at 10^-3 background SD.
#: The identifiability bound (docs/methods/operons.md "Power") allows calls only for a planted
#: spread <~ 10^-3 SD, and only with the convergent set: SW decoys are absent from the < 150 bp
#: bins of these layouts, so with both sets every q is >= the target share of those bins.
NONVAC_PLANT: tuple[int, float] = (100, 0.001)

Cfgs = tuple[Any, Any]
Transform = Callable[[Any, np.random.Generator], tuple[pd.DataFrame, pd.DataFrame | None]]


def levels() -> tuple[float, ...]:
    return tuple(C.load_tolerances().levels)


def uniop_threshold() -> float:
    from gippa.data.calibration import load

    return float(load()["uniop_pair_threshold"]["value"])


def gippa_cfg(decoys: Sequence[str] = DEFAULT_DECOYS, pi0: float = 1.0) -> Any:
    from gippa.config import GippaConfig

    return GippaConfig(operons={"decoys": list(decoys), "pi0": pi0})


def scenario(name: str, *args: Any, **kw: Any) -> Cfgs:
    import gippasim.scenarios as sc

    cfgs: Cfgs = getattr(sc, name)(*args, **kw)
    return cfgs


def reclass(pairs: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Same-strand pairs re-classed OP / NO at ``threshold`` by the step's own rule,
    ``gippa.operons.pairs.classify`` (OP iff ``uniop_p >= threshold``); CV, DV and SW unchanged.
    ``classify`` reads only the strand and the locus tag of the two genes, so they are passed as
    stand-ins."""
    from types import SimpleNamespace

    from gippa.operons.pairs import classify

    def gene(tag: object, strand: object) -> Any:
        return SimpleNamespace(locus_tag=str(tag), interval=SimpleNamespace(strand=str(strand)))

    out = pairs.copy()
    same = out["class"].isin(["OP", "NO"]).to_numpy()
    sub = out.loc[same]
    out.loc[same, "class"] = [
        classify(gene(a, st), gene(b, st), float(u), threshold)
        for a, b, st, u in zip(sub["gene_a"], sub["gene_b"], sub["strand"], sub["uniop_p"],
                               strict=True)
    ]  # fmt: skip
    return out


def simulate(cfgs: Cfgs, seed: int) -> Any:
    """One replicate, its pair classes at the calibrated UniOP threshold. The planted
    co-abundant pairs must still be NO (gippasim draws them with ``uniop_p < 0.5``): a threshold
    move that turns them OP would take them out of the missed-operon family."""
    import gippasim.scenarios as sc

    s = sc.simulate_operons(cfgs, int(seed))
    pairs = reclass(s.pairs, uniop_threshold())
    planted = set(s.truth["missed_operons"]["pair_id"].astype(str))
    cls = pairs.loc[pairs["pair_id"].astype(str).isin(planted), "class"]
    if not (cls == "NO").all():
        raise ValueError(f"planted missed operons re-classed {sorted(set(cls))} at the UniOP "
                         f"threshold {uniop_threshold()}")  # fmt: skip
    return _Sim(s.layout, s.expression, pairs, s.truth)


@dataclass(frozen=True)
class _Sim:
    layout: Any
    expression: Any
    pairs: pd.DataFrame
    truth: Mapping[str, pd.DataFrame]


def _quiet() -> None:
    # per-replicate bin warnings of the step would flood the logs (one per sparse bin per rep)
    logging.getLogger("gippa.operons.missed").setLevel(logging.ERROR)


# --------------------------------------------------------------------------- missed operons


def missed_counts(
    sim: Any,
    gcfg: Any,
    qs: Sequence[float],
    *,
    pairs: pd.DataFrame | None = None,
    genes: pd.DataFrame | None = None,
) -> dict[str, Any]:
    """Count record per level (``f"q{q:g}"``: n_called, n_false, n_planted_called, n_planted)
    plus ``min_q`` and ``n_targets`` of one run of the missed-operon step."""
    from gippa.operons.missed import missed_operons

    _quiet()
    p = sim.pairs if pairs is None else pairs
    res = missed_operons(p, gcfg, genes=genes)
    rows = res.rows.merge(p[["pair_id", "true_operon"]], on="pair_id", how="left")
    planted = set(sim.truth["missed_operons"]["pair_id"].astype(str))
    is_planted = rows["pair_id"].astype(str).isin(planted).to_numpy()
    true = is_planted | rows["true_operon"].fillna(False).astype(bool).to_numpy()
    q = rows["q_value"].astype(float).to_numpy()
    out: dict[str, Any] = {
        "min_q": float(np.nanmin(q)) if np.isfinite(q).any() else float("nan"),
        "n_targets": int(res.n_targets),
    }
    for a in qs:
        called = np.isfinite(q) & (q <= a)
        out[f"q{a:g}"] = {
            "n_called": int(called.sum()),
            "n_false": int((called & ~true).sum()),
            "n_planted_called": int((called & is_planted).sum()),
            "n_planted": len(planted),
        }
    return out


def missed_replicate_fn(
    stream: str, cfgs: Cfgs, gcfg: Any, transform: Transform | None = None
) -> Callable[[int], dict[str, Any]]:
    """``records(i)`` for ``calibrate_fdr_counts``: replicate ``i`` of ``stream`` (memoised, so
    every level reads one simulation). ``transform(sim, rng)`` returns (pairs, genes) to feed the
    step instead of the simulated pair table (attacks M-1, M-2)."""

    @cache
    def records(i: int) -> dict[str, Any]:
        seed = C.suite_seed(stream, i)
        sim = simulate(cfgs, seed)
        if transform is None:
            return missed_counts(sim, gcfg, levels())
        rng = np.random.default_rng(np.random.SeedSequence([seed, 1]))
        pairs, genes = transform(sim, rng)
        return missed_counts(sim, gcfg, levels(), pairs=pairs, genes=genes)

    return records


def missed_fdr_gate(
    suite: str,
    stream: str,
    records: Callable[[int], dict[str, Any]],
    R: int,
    *,
    log_dir: Path | None = None,
) -> list[dict[str, Any]]:
    """The FDR criterion at every q through ``calibrate_fdr_counts`` (ADR-0015). A level with no
    call in block 1 is ``vacuous`` (Khoa's decision 2026-10-22: recorded as a vacuous pass, never
    read as power evidence): the floor 0.5 / R on the mean call count is crossed by one call."""
    return [
        C.calibrate_fdr_counts(
            suite,
            stream,
            f"q{a:g}",
            a,
            records,
            R,
            vacuous={"n_called": 0.5 / R},
            log_dir=log_dir,
        )
        for a in levels()
    ]


VACUOUS_REASON = (
    "vacuous pass: no missed-operon call in block 1. The family is near-powerless by "
    "construction (a call needs near-identical abundance of the two genes; docs/methods/"
    "operons.md 'Power'); the absence of calls is not evidence against missed operons"
)


#: statuses with no verdict: a ``--shard`` range (``shard``), a merged check awaiting block 2
#: (``pending``) and a check not listed in a ``--merge`` confirm run (``skipped``); as
#: ``_gate.NOT_FAILED``, they are never a failure
NO_VERDICT = frozenset({C.SHARD_STATUS, C.PENDING_STATUS, C.SKIPPED_STATUS})
NO_VERDICT_SKIP = "shard or confirm run: the verdict comes from --merge"


def skip_without_verdict(recs: Sequence[Mapping[str, Any]]) -> None:
    """``pytest.skip`` when any record has no verdict. Called after ``calibrate*`` has logged, so
    the shard records still reach ``--merge``."""
    import pytest

    if any(r["status"] in NO_VERDICT for r in recs):
        pytest.skip(NO_VERDICT_SKIP)


def assert_fdr_ok(
    recs: Sequence[Mapping[str, Any]], *, allow_vacuous: bool, need_calls_at: float | None = None
) -> list[str]:
    """Skip a shard or confirm run (no verdict, ``skip_without_verdict``); otherwise fail on a
    confirmed ``fail``, accept ``vacuous`` only when ``allow_vacuous`` and never at the level
    ``need_calls_at``. The vacuity check runs after the skip, so a shard record (whose status is
    ``shard``, never ``vacuous``) can neither trip nor satisfy it. Returns one summary line per
    level (printed by the caller)."""
    lines = [
        f"{r['stream']} q={r['level']:g} status={r['status']} meanFDP={r['observed']}"
        f" bound={r['bound']} R={r['R']} block={r['block']}"
        + (f" [{VACUOUS_REASON}]" if r["status"] == "vacuous" else "")
        for r in recs
    ]
    print("\n".join(lines))
    skip_without_verdict(recs)
    ok = {"pass", "unconfirmed", *(["vacuous"] if allow_vacuous else [])}
    if need_calls_at is not None:
        hit = [r for r in recs if r["level"] == need_calls_at]
        assert hit and hit[0]["status"] != "vacuous", f"no call at q = {need_calls_at}: {hit}"
    bad = [line for r, line in zip(recs, lines, strict=True) if r["status"] not in ok]
    assert not bad, bad
    return lines


def fdp_summary(records: Callable[[int], dict[str, Any]], R: int, q: float) -> dict[str, float]:
    """Mean FDP, its SE, mean calls and planted recall at ``q`` over replicates 0..R-1 (from the
    memoised records; for the recorded arms and the methods note)."""
    vals = [records(i) for i in range(R)]
    c = [v[f"q{q:g}"] for v in vals]
    fdp = np.array([x["n_false"] / max(x["n_called"], 1) for x in c])
    planted = np.array([x["n_planted"] for x in c], dtype=float)
    rec = np.array([x["n_planted_called"] for x in c], dtype=float)
    return {
        "R": float(R),
        "mean_fdp": float(fdp.mean()),
        "se_fdp": float(fdp.std(ddof=1) / np.sqrt(R)) if R > 1 else float("nan"),
        "mean_calls": float(np.mean([x["n_called"] for x in c])),
        "planted_recall": float(rec.sum() / planted.sum()) if planted.sum() else float("nan"),
        "median_min_q": float(np.nanmedian([v["min_q"] for v in vals])),
    }


# --------------------------------------------------------------------------- discordance


@dataclass(frozen=True)
class DiscRun:
    p: np.ndarray  # both_detected_t p-values, pair order
    q: np.ndarray  # their BH q
    abundance: np.ndarray  # mean log2 abundance of each tested pair (of ``abundance_pairs``)
    n_one_detected: int
    pmiss: np.ndarray  # one_detected_pmiss p-values (empty until M6-05H)
    ratios: np.ndarray  # the Delta sample of the t fit (``Discordance.ratios``)
    fit: Any  # the fitted ``StudentT``


def discordance_run(
    sim: Any,
    *,
    pairs: pd.DataFrame | None = None,
    p_miss: Any = None,
    abundance_pairs: pd.DataFrame | None = None,
) -> DiscRun:
    """The discordance step on ``pairs`` (default: the simulated table). ``abundance`` is the
    pair mean read from ``abundance_pairs`` (default ``pairs``), so an attack can stratify on the
    abundances it planted from rather than on its transformed ones."""
    from gippa.config import GippaConfig
    from gippa.operons.discordance import discordance

    p = sim.pairs if pairs is None else pairs
    genes = sim.expression.genes[["locus_tag"]]
    res = discordance(p, genes, GippaConfig(), p_miss)
    ap = p if abundance_pairs is None else abundance_pairs
    rows = res.rows.merge(ap[["pair_id", "abundance_log2_a", "abundance_log2_b"]], on="pair_id")
    t = rows["test"].to_numpy() == "both_detected_t"
    m = rows["test"].to_numpy() == "one_detected_pmiss"
    op = p[p["class"] == "OP"]
    one = op["detected_a"].astype(bool).to_numpy() ^ op["detected_b"].astype(bool).to_numpy()
    ab = (rows["abundance_log2_a"].astype(float) + rows["abundance_log2_b"].astype(float)) / 2
    return DiscRun(
        p=rows.loc[t, "p_value"].astype(float).to_numpy(),
        q=rows.loc[t, "q_value"].astype(float).to_numpy(),
        abundance=ab[t].to_numpy(float),
        n_one_detected=int(one.sum()),
        pmiss=rows.loc[m, "p_value"].astype(float).to_numpy(),
        ratios=np.asarray(res.ratios, dtype=float),
        fit=res.fit,
    )
