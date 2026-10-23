"""giapsim function worlds -> GIAP inputs for the M8-15 gate (docs/specs/M8-15-calibration-gate.md;
M8-14 scenarios; GIAP_SPEC.md §7.3, §7.7, §7.8, §8.4, §8.5 attacks 2 and 3).

Test-side glue only: every statistic comes from GIAP (``covariate_table``, ``fit_detectability``,
``build_collections``, ``detection_enrichment``, ``gsea_table``, ``cross_validate``, ``null_p``,
``shap_blocks``, ``attribution_table``). The simulator's truth columns are used only where the
spec names them (the oracle arm, the attack-2 confounded set, the recovery truth).

- ``gene_frame(w)``: one row per simulated gene with the gene-table columns the M8 steps read
  (``locus_tag``, ``tu_id``, ``tu_position``, ``partial``, ``detected``, ``abundance_log2`` (NaN
  when undetected), ``abundance_imputed``, ``abundance_scale``, ``cog_category``, ``kegg_ko``,
  ``go_terms``; the last three from the layout's C3 xrefs through ``function_columns``).
- ``fitted_p_detect(w, genes, rng)``: the M4-13 cross-fitted ``p_detect`` with the **default**
  ``DetectabilityConfig()`` on ``covariate_table(layout.genes, GiapConfig())`` (covariates
  recomputed from the sequences; never the ``n_obs_true`` column; lead decision 4) and the M4-13
  ``fold_of_tu`` (passed to the ORA TU recalibration, M8-04 decision 13).
- ``resources(w, tmp)`` / ``collections(w, genes, tmp)``: the M8-03 collections of the synthetic
  ontology, KEGG maps and COG letters, **six** collections: ``go_slim`` is a synthetic slim, every
  third non-root term of the OBO by ID (M8-14 gap 4: "a synthetic go_slim goes to the backlog for
  M8-15"; test-side, so no simulator change).
- ``confounded_term(w)``: the GO term ``ora_attack_confounded`` flags (``sets.confounded``).
"""

from __future__ import annotations

import re
import tempfile
from collections.abc import Sequence
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

SLIM_EVERY = 3  # synthetic slim: every third non-root synthetic term (by ID)
_TERM = re.compile(r"\[Term\]\n(.*?)(?=\n\[|\Z)", re.S)


def _fs() -> Any:
    import giapsim.function_sim as fs

    return fs


def tu_maps(w: Any) -> tuple[dict[str, str], dict[str, int]]:
    """locus_tag -> tu_id and locus_tag -> 1-based position in the TU (layout order)."""
    tu: dict[str, str] = {}
    pos: dict[str, int] = {}
    for t, gs in zip(w.layout.tus["tu_id"], w.layout.tus["genes"], strict=True):
        for k, g in enumerate(gs, start=1):
            tu[str(g)], pos[str(g)] = str(t), k
    return tu, pos


def gene_frame(w: Any) -> pd.DataFrame:
    """The gene-table columns of the M8 steps, one row per gene in layout order."""
    from giap.genetable.assemble import function_columns

    lg = w.layout.genes
    tags = lg["locus_tag"].astype(str)
    e = w.expression.genes.set_index(w.expression.genes["locus_tag"].astype(str)).loc[tags]
    tu, pos = tu_maps(w)
    det = e["detected_true"].to_numpy(bool)
    fc = function_columns(w.layout.xrefs, w.cog_def).set_index("locus_tag").reindex(tags)
    return pd.DataFrame({
        "locus_tag": tags.to_numpy(), "tu_id": [tu[t] for t in tags],
        "tu_position": [pos[t] for t in tags],
        "partial": lg["partial"].astype(bool).to_numpy(), "detected": det,
        "abundance_log2": np.where(det, e["abundance_log2_true"].to_numpy(float), np.nan),
        "abundance_imputed": False, "abundance_scale": "molar",
        "cog_category": fc["cog_category"].to_numpy(object),
        "kegg_ko": [v if isinstance(v, list) else [] for v in fc["kegg_ko"]],
        "go_terms": [v if isinstance(v, list) else [] for v in fc["go_terms"]],
    })  # fmt: skip


def recomputed_covariates(w: Any) -> pd.DataFrame:
    """M4-10 ``covariate_table`` on the simulated sequences (``aa``, else translated ``nt``)."""
    from giap.config.model import GiapConfig
    from giap.genetable.covariates import covariate_table

    return covariate_table(w.layout.genes, GiapConfig())


def fitted_p_detect(
    w: Any, genes: pd.DataFrame, rng: np.random.Generator
) -> tuple[pd.Series, dict[str, int]]:
    """Cross-fitted p_detect (default ``DetectabilityConfig()``) and the M4-13 TU folds."""
    from giap.config.model import DetectabilityConfig
    from giap.genetable.detectability import fit_detectability

    tags = genes["locus_tag"].astype(str)
    det = pd.Series(genes["detected"].to_numpy(bool), index=tags.to_numpy(), dtype="boolean")
    det[genes["partial"].to_numpy(bool)] = pd.NA
    tus = pd.Series(genes["tu_id"].to_numpy(), index=tags.to_numpy())
    fit = fit_detectability(det, recomputed_covariates(w), tus, DetectabilityConfig(), rng)
    assert not fit.constant, "detectability_constant on a giapsim world"
    return fit.p_detect.reindex(tags.to_numpy()), dict(fit.fold_of_tu)


def with_p_detect(genes: pd.DataFrame, p: Sequence[float] | pd.Series) -> pd.DataFrame:
    return genes.assign(p_detect=np.asarray(p, dtype=float))


def true_p_detect(w: Any, genes: pd.DataFrame) -> np.ndarray:
    e = w.expression.genes.set_index(w.expression.genes["locus_tag"].astype(str))
    return e.loc[genes["locus_tag"].astype(str), "p_detect_true"].to_numpy(float)


def slim_obo(obo: str) -> str:
    """A synthetic GO slim: every ``SLIM_EVERY``-th non-root term of ``obo`` (id, name,
    namespace only), in the reserved range like its source."""
    keep: list[str] = []
    stanzas = _TERM.findall(obo)
    non_root = [s for s in stanzas if "\nis_a:" in "\n" + s or "relationship:" in s]
    ids = sorted(re.search(r"^id: (\S+)", s, re.M).group(1) for s in non_root)  # type: ignore[union-attr]
    chosen = set(ids[::SLIM_EVERY])
    for s in stanzas:
        tid = re.search(r"^id: (\S+)", s, re.M).group(1)  # type: ignore[union-attr]
        if tid in chosen:
            lines = [
                ln for ln in s.splitlines() if ln.split(":", 1)[0] in {"id", "name", "namespace"}
            ]
            keep.append("[Term]\n" + "\n".join(lines) + "\n")
    head = "format-version: 1.2\ndata-version: giapsim/synthetic-slim\nontology: go/subsets/slim\n"
    return head + "\n" + "\n".join(keep)


def resources(w: Any, tmp: Path) -> Any:
    from giap.function.collections import Resources
    from giap.io.function_resources import read_go

    go, slim = tmp / "go.obo", tmp / "slim.obo"
    go.write_text(w.obo)
    slim.write_text(slim_obo(w.obo))
    return Resources(go=read_go(go), slim=read_go(slim), kegg_links=w.kegg_links,
                     kegg_names=w.kegg_names, cog_names=w.cog_names)  # fmt: skip


def collections(w: Any, genes: pd.DataFrame) -> list[Any]:
    """The six M8-03 collections (cog, kegg, go_bp, go_mf, go_cc, go_slim) of the world."""
    from giap.function.collections import build_collections

    with tempfile.TemporaryDirectory(prefix="giap-m815-") as d:
        res = resources(w, Path(d))
    el = ~genes["partial"].to_numpy(bool)
    colls, flags = build_collections(genes, el, res)
    assert not [f for f in flags if f.startswith("layer_missing")], flags
    return colls


def confounded_term(w: Any) -> str:
    s = w.sets
    ids = sorted(set(s.loc[s["confounded"].astype(bool), "set_id"].astype(str)))
    assert len(ids) == 1, ids
    return ids[0]


def single_term(colls: Sequence[Any], term: str) -> Any:
    """The GeneSets of ``term`` alone, collection name kept (no elim partner: p = p_noelim)."""
    from giap.function.collections import GeneSets

    for c in colls:
        if term in c.sets:
            return GeneSets(c.collection, {term: c.sets[term]}, {term: c.names.get(term, term)},
                            {}, c.coverage, c.n_unknown_ids)  # fmt: skip
    raise KeyError(f"{term} in no collection")


@cache
def fs_constants() -> tuple[str, float]:
    fs = _fs()
    return str(fs.CONFOUND_COVARIATE), float(fs.SHIFT_TOL)


# --------------------------------------------------------------------------- records

DIAGNOSTICS = "function_gate_diagnostics.jsonl"


def diagnostics(record: dict[str, Any]) -> None:
    """Append one recorded-arm line to ``<log_dir>/function_gate_diagnostics.jsonl`` (never the
    suite log: a recorded arm carries no verdict, so ``--merge`` must not see it)."""
    import json

    import _criteria as C

    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / DIAGNOSTICS).open("a") as fh:
        fh.write(json.dumps(record, default=float) + "\n")


def ks_record(p: Sequence[float] | np.ndarray) -> dict[str, float]:
    """KS of pooled p against U(0, 1), recorded only (lead decision 2): statistic, p and the
    tie share (discrete p: Wallenius / CPB, permutation and multilevel p)."""
    from scipy import stats

    a = np.asarray(p, dtype=float).ravel()
    a = a[np.isfinite(a)]
    if a.size == 0:
        return {"n": 0, "ks_stat": float("nan"), "ks_p": float("nan"), "tie_share": float("nan")}
    ks = stats.kstest(a, "uniform")
    return {"n": int(a.size), "ks_stat": float(ks.statistic), "ks_p": float(ks.pvalue),
            "tie_share": float(1.0 - np.unique(a).size / a.size)}  # fmt: skip


def rate_record(p: Sequence[float] | np.ndarray, levels: Sequence[float]) -> dict[str, Any]:
    """Rejection rate per alpha with the §8.4 two-sided band alpha -+ k SE (recorded arms)."""
    import _criteria as C

    a = np.asarray(p, dtype=float).ravel()
    a = a[np.isfinite(a)]
    k = C.load_tolerances().se_multiplier
    out: dict[str, Any] = {}
    for al in levels:
        se = float(np.sqrt(al * (1 - al) / max(a.size, 1)))
        rate = float(np.mean(a <= al)) if a.size else float("nan")
        out[f"{al:g}"] = {"rate": rate, "upper": al + k * se, "lower": al - k * se,
                          "above_upper": rate > al + k * se,
                          "below_lower": rate < al - k * se}  # fmt: skip
    return out


def by_index(stream: str, R: int, fn: Any) -> Any:
    """``calibrate``'s ``fn(rng)`` for a function of the replicate index ``i`` (blocks 1 and 2):
    one shared world per index feeds every check; the stream only names the check."""
    import _criteria as C

    index = {C.suite_seed(stream, i): i for i in range(2 * R)}
    return lambda rng: fn(index[int(rng.bit_generator.seed_seq.entropy)])


def shares(p: np.ndarray, levels: Sequence[float]) -> np.ndarray:
    """Per-alpha rejection share of one replicate's pooled p (``check_level_pooled`` values)."""
    a = np.asarray(p, dtype=float)
    assert a.ndim == 1 and a.size > 0 and np.all(np.isfinite(a)), a
    return np.array([np.count_nonzero(a <= al) / a.size for al in levels])


def r_engine(tmp: Path) -> Any:
    import os

    from giap.io.fgsea_r import r_fgsea_engine, resolve_launcher

    return r_fgsea_engine(resolve_launcher(os.environ), tmp_root=tmp)


def cpu_share() -> int:
    """This process's share of the job's CPUs: ``SLURM_CPUS_PER_TASK`` (default 1) divided by
    the xdist worker count (``PYTEST_XDIST_WORKER_COUNT``, default 1); at least 1."""
    import os

    cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    n = max(1, int(os.environ.get("PYTEST_XDIST_WORKER_COUNT", "1")))
    return max(1, cpus // n)


def workers() -> int:
    """Processes for ``prefetch``: ``GIAP_M815_WORKERS`` (default 4), at most ``cpu_share()``
    (review nit 1: the same per-xdist-worker share as the determinants ``_threads``)."""
    import os

    return max(1, min(int(os.environ.get("GIAP_M815_WORKERS", "4")), cpu_share()))


def prefetch(cache: dict[Any, Any], run: Any, keys: Sequence[Any], done: Any = None) -> None:
    """Fill ``cache[k] = run(k)`` for the missing ``keys`` in a forked process pool (results are
    pure functions of ``k``, so the values equal a sequential run; only wall time changes).
    ``done(value)`` False marks a cached value as incomplete (recomputed)."""
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context

    missing = [k for k in keys if k not in cache or (done is not None and not done(cache[k]))]
    n = workers()
    if not missing:
        return
    if n == 1 or len(missing) == 1:
        for k in missing:
            cache[k] = run(k)
        return
    with ProcessPoolExecutor(n, mp_context=get_context("fork")) as ex:
        for k, v in zip(missing, ex.map(run, missing), strict=True):
            cache[k] = v


def prefetch_indices(R: int) -> list[int]:
    """Absolute replicate indices the first ``calibrate`` pass of this process runs: unsharded,
    block 1 (``0..R-1``); under ``--shard K/N`` the shard's range of block
    ``GIAP_CALIBRATION_BLOCK``, offset by ``(block - 1) * R`` (review B2); none in a confirm run
    (``GIAP_CALIBRATION_CONFIRM``: its block-2 checks are only known to ``calibrate``)."""
    import os

    import _criteria as C

    if os.environ.get(C.CONFIRM_ENV):
        return []
    local = C.replicate_range(R)
    if local == range(R):
        return list(range(R))
    offset = (C.replicate_block() - 1) * R
    return [offset + i for i in local]


def sharded() -> bool:
    """A ``--shard`` or ``--merge`` confirm run: recorded-only arms (no ``calibrate`` verdict)
    skip there, so they run once per suite run, not once per shard (review B1)."""
    import os

    import _criteria as C

    return any(os.environ.get(v) for v in (C.SHARD_ENV, C.CONFIRM_ENV))
