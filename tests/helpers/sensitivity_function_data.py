"""Inputs for the M8-16 tests: the ``enriched_categories`` and ``determinants`` sensitivity families
(docs/specs/M8-16-sensitivity-function-determinants.md, lead decisions of 2026-10-23; GIAP_SPEC.md
§8.6, §6.2 ``sensitivity`` row).

Interface the tests assume (test author; each point is reported to the lead):

* ``giap.function.sensitivity.enriched_categories_family(genes, cfg, rng, *, collections, tus,
  engine=...)`` returns a set of ``"{collection}:{term}:{kind}"`` IDs, kind in ``gsea``,
  ``ora_over``, ``ora_under``, with q <= ``stats.alpha.enrichment``. ``engine`` is a keyword-only
  ``GseaEngine`` (the spec's amendment: "its tests inject a fake or ``SIMPLE_ENGINE``"); the step
  binds R fgsea. GSEA runs at ``n_perm_simple = 199`` and ``sample_size = 101`` whatever
  ``function.gsea`` says; the ORA is ``detection_enrichment`` (deterministic).
* ``giap.determinants.sensitivity.determinants_family(genes, cfg, rng, *, features)`` returns a
  ``DeterminantsSensitivity`` with ``flagged`` (frozenset of the locus tags flagged
  ``abundance_not_explained_by_sequence``), ``dominant`` (top-decile locus tag -> dominant block,
  ``None`` for a top-decile gene without one) and ``cv_r2`` (mean LightGBM CV R² over
  ``determinants.repeats``); below the 100-row minimum, or without abundance, it returns ``None``
  or a result with ``cv_r2`` NaN and empty ``flagged`` / ``dominant``.
* ``giap.context.sensitivity.sensitivity_rows`` accepts ``DeterminantsSensitivity`` values for a
  family and fills ``agreement`` for it; ``n_*`` / ``jaccard`` use ``flagged`` (lead decision 2).

The gene tables here are hand-built §9.3 frames (``function_step_data.full_frame``): one gene per
TU, a circular ``chr1``. Nothing here imports the modules under test at import time.
"""

from __future__ import annotations

import inspect
import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from helpers import determinants_models_data as dm
from helpers.function_step_data import CONTIG, full_frame

SENS_B = 199  # §8.6 "Sensitivity reruns use B = 199"
SAMPLE_SIZE = 101  # spec item 1: the multilevel stage at sampleSize = 101
SEED = 20261023
EC = "enriched_categories"
DET = "determinants"


def config(**over: Any) -> Any:
    """``GiapConfig`` with a fast determinants block (2 repeats x 3 folds, 30 trees, one thread)
    and nested overrides."""
    from giap.config.model import GiapConfig

    raw: dict[str, Any] = {
        "seed": SEED,
        "threads": 1,
        "stats": {"bootstrap": 30},
        "determinants": {"repeats": 2, "outer_folds": 3, "lightgbm": {"n_estimators": 30},
                         "elastic_net": {"inner_folds": 2}},
    }  # fmt: skip
    for k, v in over.items():
        if isinstance(v, Mapping) and isinstance(raw.get(k), dict):
            raw[k] = {**raw[k], **v}
        else:
            raw[k] = v
    return GiapConfig.model_validate(raw)


def base_config(**over: Any) -> Any:
    """What the step hands every family: ``config`` at ``stats.permutations = 199``."""
    stats = {"bootstrap": 30, "permutations": SENS_B, **dict(over.pop("stats", {}))}
    return config(stats=stats, **over)


# --------------------------------------------------------------------------- gene tables


def gene_table(
    abundance: Sequence[float],
    detected: Sequence[bool],
    p_detect: Sequence[float],
    *,
    imputed: Sequence[bool] | None = None,
    status: Sequence[str] | None = None,
    cog: Sequence[str | None] | None = None,
    tier: str = "one_percent_fdr",
) -> pd.DataFrame:
    """A §9.3 gene table, one gene per TU (``g1``... in ``TU_1``...); NaN abundance where not
    given. ``status`` defaults to detected / not_detected."""
    n = len(abundance)
    det = [bool(d) for d in detected]
    st = list(status) if status is not None else ["detected" if d else "not_detected" for d in det]
    ab = [float(a) for a in abundance]
    cols: dict[str, Sequence[Any]] = {
        "locus_tag": [f"g{i + 1}" for i in range(n)],
        "feature_type": ["cds"] * n,
        "contig": [CONTIG] * n,
        "start": [1000 * i for i in range(n)],
        "end": [1000 * i + 900 for i in range(n)],
        "strand": ["+"] * n,
        "partial": [False] * n,
        "crosses_origin": [False] * n,
        "hypothetical": [False] * n,
        "replicon_type": ["chromosome"] * n,
        "tu_id": [f"TU_{i + 1}" for i in range(n)],
        "tu_position": [1] * n,
        "tu_size": [1] * n,
        "detected": det,
        "detection_tier": [tier] * n,
        "detection_status": st,
        "p_detect": [float(p) for p in p_detect],
        "abundance_log2": ab,
        "abundance_raw": [2.0**a if math.isfinite(a) else math.nan for a in ab],
        "abundance_imputed": [False] * n if imputed is None else [bool(x) for x in imputed],
        "imputation_policy": ["median"] * n,
        "abundance_type": ["ibaq"] * n,
        "abundance_scale": ["absolute"] * n,
        "top3_log2": [math.nan] * n,
        "cog_category": [None] * n if cog is None else list(cog),
        "annotation_source": ["bakta_full"] * n,
        "flags": [["imputed_abundance"] if (imputed is not None and imputed[i]) else []
                  for i in range(n)],
    }  # fmt: skip
    df = full_frame(cols)
    for c in ("detected", "partial", "abundance_imputed"):
        df[c] = df[c].astype(bool)
    for c in ("p_detect", "abundance_log2", "abundance_raw", "top3_log2"):
        df[c] = df[c].astype("float64")
    return df


def tags(idx: Sequence[int] | range) -> list[str]:
    return [f"g{i + 1}" for i in idx]


def gene_sets(collection: str, sets: Mapping[str, Sequence[str]]) -> Any:
    """``GeneSets`` of ``collection`` (names = term IDs; GO terms without parents)."""
    from giap.function.collections import GeneSets

    s = {t: frozenset(m) for t, m in sets.items()}
    parents = {t: frozenset() for t in s} if collection.startswith("go_") else {}
    return GeneSets(collection, s, {t: t for t in s}, parents, math.nan, 0)  # type: ignore[arg-type]


def tus_frame(genes: pd.DataFrame) -> pd.DataFrame:
    """C5 ``tus`` of a ``gene_table`` (one TU per gene)."""
    return pd.DataFrame({
        "tu_id": genes["tu_id"].astype(str), "contig": genes["contig"].astype(str),
        "start": genes["start"].astype("int64"), "end": genes["end"].astype("int64"),
        "strand": genes["strand"].astype(str), "genes": [[t] for t in genes["locus_tag"]],
    })  # fmt: skip


N_EX1 = 200
K_MEMBERS = range(0, 12)
X_MEMBERS = range(12, 24)
MIN_GENE = 199  # the lowest observed abundance (a background gene)


def example1_genes() -> pd.DataFrame:
    """Spec example 1's table (primary ``median``): COG K = g1-g12, detected at p_detect 0.99,
    observed log2 5; GO BP X = g13-g24, detected at p_detect 0.99 but imputed (filled with the
    observed median); background g25-g200 alternate detected (log2 in [-1, 1]) / not detected at
    p_detect 0.5, with g200 the lowest observed (-4)."""
    n = N_EX1
    det = [True] * 24 + [i % 2 == 0 for i in range(24, n)]
    det[MIN_GENE] = True
    p = [0.99] * 24 + [0.5] * (n - 24)
    ab = [math.nan] * n
    for i in K_MEMBERS:
        ab[i] = 5.0
    bg = [i for i in range(24, n) if det[i] and i != MIN_GENE]
    for j, i in enumerate(bg):
        ab[i] = -1.0 + 2.0 * j / max(len(bg) - 1, 1)
    ab[MIN_GENE] = -4.0
    observed = [a for a in ab if math.isfinite(a)]
    med = float(np.median(observed))
    imp = [False] * n
    for i in X_MEMBERS:
        ab[i], imp[i] = med, True
    return gene_table(ab, det, p, imputed=imp)


def example1_collections() -> list[Any]:
    return [gene_sets("cog", {"K": tags(K_MEMBERS)}),
            gene_sets("go_bp", {"X": tags(X_MEMBERS)})]  # fmt: skip


def example2_genes() -> pd.DataFrame:
    """Spec example 2 (``two_peptide``): COG J = g1-g20, every one single-peptide (not detected
    in the primary, p_detect 0.5); COG C = g21-g40, detected and observed; background g41-g200
    alternate detected / not detected at p_detect 0.5."""
    n = 200
    det = [False] * 20 + [True] * 20 + [i % 2 == 0 for i in range(40, n)]
    status = (["single_peptide_unestimated"] * 20
              + ["detected" if d else "not_detected" for d in det[20:]])  # fmt: skip
    ab = [math.nan if not d else 0.05 * (i % 40) for i, d in enumerate(det)]
    return gene_table(ab, det, [0.5] * n, status=status, tier="two_peptide")


def example2_collections() -> list[Any]:
    return [gene_sets("cog", {"J": tags(range(20)), "C": tags(range(20, 40))})]


def null_genes(n: int = 240) -> pd.DataFrame:
    """Every gene detected (p_detect 0.9), log2 abundance spread over [-2, 2]: the ORA has no
    signal (N = every TU), so only an injected GSEA engine decides the flagged set."""
    ab = [-2.0 + 4.0 * ((7 * i) % n) / (n - 1) for i in range(n)]
    return gene_table(ab, [True] * n, [0.9] * n)


def null_collections(collection: str, n_terms: int, size: int = 12) -> list[Any]:
    terms = {f"T{j + 1}": tags(range(j * size, (j + 1) * size)) for j in range(n_terms)}
    return [gene_sets(collection, terms)]


# --------------------------------------------------------------------------- fake GSEA engines


@dataclass
class FakeEngine:
    """A ``GseaEngine`` whose p per term is ``rule(job_id, term, member_scores)``; records every
    call. ``es`` is +1 (sign irrelevant to the flagged set)."""

    rule: Callable[[str, str, np.ndarray], float]
    calls: list[Any] = field(default_factory=list)
    method: str = "fake_engine"

    def __call__(self, call: Any) -> Mapping[str, tuple[Any, ...]]:
        from giap.function.gsea import EngineRow

        self.calls.append(call)
        score = dict(zip(call.tu_ids, np.asarray(call.scores, dtype=float), strict=True))
        out: dict[str, tuple[Any, ...]] = {}
        for job in call.jobs:
            rows = []
            for term, members in job.sets.items():
                s = np.array([score[t] for t in members], dtype=float)
                p = float(self.rule(job.job_id, term, s))
                rows.append(EngineRow(term, len(set(members)), 1.0, 1.0, p, 0.1, tuple(members)))
            out[job.job_id] = tuple(rows)
        return out


def constant_engine(p: float) -> FakeEngine:
    return FakeEngine(lambda _j, _t, _s: p)


def table_engine(p: Mapping[str, float], default: float = 0.9) -> FakeEngine:
    """p by term ID (``default`` otherwise)."""
    return FakeEngine(lambda _j, t, _s: float(p.get(t, default)))


def example1_engine() -> FakeEngine:
    """p = 0.001 when no member TU scores below -2 (median-centred), else 0.9."""
    return FakeEngine(lambda _j, _t, s: 0.001 if s.size and float(s.min()) > -2.0 else 0.9)


def ec_family() -> Callable[..., Any]:
    from giap.function.sensitivity import enriched_categories_family

    return enriched_categories_family


def call_ec(genes: pd.DataFrame, cfg: Any, rng: np.random.Generator, colls: Sequence[Any],
            engine: Any) -> frozenset[str]:  # fmt: skip
    """The family as the step calls it (``collections``, ``tus`` bound), ``engine`` injected."""
    res = ec_family()(genes, cfg, rng, collections=list(colls), tus=tus_frame(genes), engine=engine)
    assert isinstance(res, set | frozenset), type(res)
    return frozenset(res)


def with_engine(family: Callable[..., Any], engine: Any) -> Callable[..., Any]:
    """``family`` with ``engine`` forced, its signature (minus ``engine``) kept so the step's
    ``_bind`` binds the same tables by name."""
    sig = inspect.signature(family)

    def wrapped(genes: pd.DataFrame, cfg: Any, rng: np.random.Generator, **kw: Any) -> Any:
        kw["engine"] = engine
        return family(genes, cfg, rng, **kw)

    wrapped.__signature__ = sig.replace(  # type: ignore[attr-defined]
        parameters=[p for p in sig.parameters.values() if p.name != "engine"]
    )
    return wrapped


# --------------------------------------------------------------------------- oracles


def ora_ids(genes: pd.DataFrame, colls: Sequence[Any], cfg: Any) -> frozenset[str]:
    """``{collection}:{term}:ora_<direction>`` of the M8-04 ORA rows with q <= alpha (no folds)."""
    from giap.function.ora import detection_enrichment

    lo, hi = cfg.function.ora_size
    alpha = cfg.stats.alpha.enrichment
    res = detection_enrichment(genes, list(colls), (lo, hi), alpha, folds=None)
    q = res["q_value"].astype(float)
    keep = res[q.notna() & (q <= alpha)]
    cols = zip(keep["collection"], keep["term"], keep["direction"], strict=True)
    return frozenset(f"{c}:{t}:ora_{d}" for c, t, d in cols)


def gsea_ids(p_by_term: Mapping[str, float], collection: str, alpha: float) -> frozenset[str]:
    """BH within one collection over the engine p-values."""
    from giap.stats.multitest import bh

    terms = list(p_by_term)
    q = bh(np.array([p_by_term[t] for t in terms], dtype=float))
    return frozenset(f"{collection}:{t}:gsea" for t, qv in zip(terms, q, strict=True)
                     if qv <= alpha)  # fmt: skip


# --------------------------------------------------------------------------- determinants


def planted_world(
    n: int = 180, seed: int = 5, signal: float = 2.0
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(genes, features): ``cai`` ~ N(0, 1) drives log2 abundance (``signal`` * cai + N(0, 1));
    noise columns in the other four blocks; TUs of 1-3 genes (>= 25 TUs); every gene detected,
    observed, not ribosomal. The features frame is the M8-07 TSV frame of the same genes."""
    r = np.random.default_rng(seed)
    sizes = dm.mixed_sizes(r, n, 1, 3)
    tu = dm.tu_ids(sizes)
    tg = [f"g{i + 1}" for i in range(n)]
    cols = {c: r.normal(0.0, 1.0, n).tolist() for c in dm.PLANTED_BLOCKS}
    y = signal * np.asarray(cols["cai"]) + r.normal(0.0, 1.0, n)
    feats = dm.features(cols, tu=tu, tags=tg, abundance=y.tolist())
    genes = gene_table(y.tolist(), [True] * n, [0.5] * n)
    genes["tu_id"] = tu
    return genes, feats


def merged_features(genes: pd.DataFrame, features: pd.DataFrame) -> pd.DataFrame:
    """``features`` with ``detected`` / ``abundance_log2`` / ``abundance_imputed`` taken from
    ``genes`` (by locus tag): what a variant's gene table changes."""
    g = genes.set_index(genes["locus_tag"].astype(str))
    f = features.copy()
    tags_ = f["locus_tag"].astype(str)
    f["detected"] = pd.array(g.loc[tags_, "detected"].astype(bool).to_numpy(), dtype="boolean")
    f["abundance_log2"] = g.loc[tags_, "abundance_log2"].astype("float64").to_numpy()
    f["abundance_imputed"] = pd.array(g.loc[tags_, "abundance_imputed"].astype(bool).to_numpy(),
                                      dtype="boolean")  # fmt: skip
    return f


@dataclass(frozen=True)
class DetOracle:
    flagged: frozenset[str]
    dominant: dict[str, str]  # top-decile genes with a dominant block
    top: frozenset[str]
    cv_r2: float


def det_oracle(genes: pd.DataFrame, features: pd.DataFrame, cfg: Any,
               rng: np.random.Generator) -> DetOracle:  # fmt: skip
    """The M8-10 chain without its null: model rows, ``model_matrix``, folds from ``rng``
    (``tu_folds``), ``cross_validate`` (seed ``cfg.seed``), SHAP and ``attribution_table``."""
    from giap.determinants.attribution import UNEXPLAINED, attribution_table, shap_blocks
    from giap.determinants.features import model_matrix
    from giap.determinants.models import cross_validate, model_rows, tu_folds

    f = merged_features(genes, features)
    rows = model_rows(f)
    X_all, blocks = model_matrix(f)
    X = X_all[rows]
    y = pd.to_numeric(f.loc[rows, "abundance_log2"]).to_numpy(float)
    tu = f.loc[rows, "tu_id"].astype(str).tolist()
    d = cfg.determinants
    folds = tu_folds(tu, d.outer_folds, d.repeats, rng)
    cv = cross_validate(X, y, folds, blocks, d, cfg.seed, 1, tu=tu)
    shap = shap_blocks(cv["lightgbm"], X, blocks)
    shap.index = pd.Index(f.loc[shap.index, "locus_tag"].astype(str))
    att = attribution_table(f, rows, cv["lightgbm"], shap)
    flagged = frozenset(t for t, fl in zip(att["locus_tag"], att["flags"], strict=True)
                        if UNEXPLAINED in str(fl).split(";"))  # fmt: skip
    top = att[att["top_decile"].astype(bool)]
    dom = {str(t): str(b) for t, b in zip(top["locus_tag"], top["dominant_block"], strict=True)
           if not pd.isna(b)}  # fmt: skip
    return DetOracle(flagged, dom, frozenset(top["locus_tag"].astype(str)),
                     float(np.mean(cv["lightgbm"].r2)))  # fmt: skip


def det_family() -> Callable[..., Any]:
    from giap.determinants.sensitivity import determinants_family

    return determinants_family


def det_result(
    genes: pd.DataFrame, features: pd.DataFrame, cfg: Any, rng: np.random.Generator
) -> Any:
    return det_family()(genes, cfg, rng, features=features)


def non_none(dominant: Mapping[str, str | None]) -> dict[str, str]:
    return {str(k): str(v) for k, v in dominant.items() if v is not None and not pd.isna(v)}


def det_sens(flagged: Sequence[str], dominant: Mapping[str, str | None], cv_r2: float) -> Any:
    from giap.determinants.sensitivity import DeterminantsSensitivity

    return DeterminantsSensitivity(flagged=frozenset(flagged), dominant=dict(dominant),
                                   cv_r2=cv_r2)  # fmt: skip
