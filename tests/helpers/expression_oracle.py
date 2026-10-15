"""Independent oracle for the giapsim expression model (docs/specs/M4-08-sim-expression.md).

Written from the spec text, not from ``giapsim.expression``: Trypsin/P digestion and distinct
7-30 aa peptides (via ``covariates_oracle``, the M4-10 definition), Kyte-Doolittle GRAVY as a plain
mean, the detection logit of the spec, TU order per replicon from the ``TU_<n>`` numbering (contig
order, then start; spec §6.5), the AR(1) whitening transform, lag-1 autocorrelation and a small
Newton-Raphson logistic regression (numpy only).
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd
from covariates_oracle import KD, n_observable

_TU_NUM = re.compile(r"^TU_(\d+)$")


def gravy(aa: str) -> float:
    """Mean Kyte-Doolittle hydropathy over the residues (unknown letters 0)."""
    return float(np.mean([KD.get(c, 0.0) for c in aa])) if aa else math.nan


def n_obs(aa: str) -> int:
    """Distinct Trypsin/P peptides (0 missed cleavages) of 7-30 aa."""
    return n_observable(aa, "trypsin_p", (7, 30))


def sigmoid(x: Any) -> Any:
    return 1.0 / (1.0 + np.exp(-np.asarray(x, dtype=float)))


def logit(p: float) -> float:
    return math.log(p / (1.0 - p))


def is_na(v: object) -> bool:
    if v is None:
        return True
    try:
        return bool(pd.isna(v))
    except (TypeError, ValueError):
        return False


def as_list(v: object) -> list[Any]:
    if v is None:
        return []
    if isinstance(v, str):
        return [v]
    return list(v)  # type: ignore[call-overload]


def tu_number(tu_id: str) -> int:
    m = _TU_NUM.match(str(tu_id))
    assert m, tu_id
    return int(m.group(1))


def tu_order(layout: Any) -> dict[str, list[str]]:
    """TU ids per contig in replicon order (``TU_<n>`` numbering = contig order, then start)."""
    tus = layout.tus.assign(_n=[tu_number(t) for t in layout.tus["tu_id"]]).sort_values("_n")
    out: dict[str, list[str]] = {}
    for tid, contig in zip(tus["tu_id"], tus["contig"], strict=True):
        out.setdefault(str(contig), []).append(str(tid))
    return out


def tu_members(layout: Any) -> dict[str, list[str]]:
    """TU id -> member locus tags in transcription order (C5 ``genes``)."""
    return {str(t): [str(g) for g in as_list(gs)] for t, gs in
            zip(layout.tus["tu_id"], layout.tus["genes"], strict=True)}  # fmt: skip


def gene_tu_position(layout: Any) -> dict[str, tuple[str, int]]:
    """locus tag -> (TU id, 0-based position in transcription order)."""
    out: dict[str, tuple[str, int]] = {}
    for t, gs in tu_members(layout).items():
        for k, g in enumerate(gs):
            out[g] = (t, k)
    return out


def region_types_of_genes(layout: Any) -> dict[str, set[str]]:
    """locus tag -> set of region types it belongs to (C7 + C8)."""
    typ = dict(zip(layout.regions["region_id"], layout.regions["region_type"], strict=True))
    out: dict[str, set[str]] = {}
    for rid, g in zip(layout.region_genes["region_id"], layout.region_genes["locus_tag"],
                      strict=True):  # fmt: skip
        out.setdefault(str(g), set()).add(str(typ[rid]))
    return out


def tu_region_status(layout: Any, region_type: str) -> dict[str, str]:
    """TU id -> ``"all"`` (every member in a region of the type), ``"none"`` or ``"some"``."""
    rt = region_types_of_genes(layout)
    out = {}
    for t, gs in tu_members(layout).items():
        k = sum(region_type in rt.get(g, set()) for g in gs)
        out[t] = "all" if k == len(gs) else ("none" if k == 0 else "some")
    return out


def merged(layout: Any, expr_genes: pd.DataFrame) -> pd.DataFrame:
    """Layout genes (C2) joined 1:1 with the expression gene frame on ``locus_tag``."""
    g = layout.genes[["locus_tag", "contig", "aa", "feature_type", "partial"]]
    m = g.merge(expr_genes, on="locus_tag", how="inner", validate="one_to_one",
                suffixes=("", "_expr"))  # fmt: skip
    assert len(m) == len(g) == len(expr_genes)
    return m


def expected_logit(
    m: pd.DataFrame,
    beta: Mapping[str, float],
    detection: str,
    region_shift: Mapping[str, float],
    region_types: Mapping[str, set[str]],
) -> np.ndarray:
    """Spec logit for genes with an amino-acid sequence (rows of ``m`` with non-NA ``aa``)."""
    out = []
    for g, aa, ab in zip(m["locus_tag"], m["aa"], m["abundance_log2_true"], strict=True):
        x = (
            beta["intercept"]
            + beta["log_length"] * math.log(len(aa))
            + beta["log1p_n_obs"] * math.log1p(n_obs(aa))
            + beta["gravy"] * gravy(aa)
        )
        if detection == "abundance":
            x += beta["log2_abundance"] * float(ab)
        x += sum(region_shift.get(t, 0.0) for t in region_types.get(g, set()))
        out.append(x)
    return np.array(out)


def tu_series(layout: Any, expr_genes: pd.DataFrame) -> dict[str, np.ndarray]:
    """Per contig: the TU-level ``tu_log2`` in replicon TU order (one value per TU)."""
    tl = dict(zip(expr_genes["locus_tag"], expr_genes["tu_log2"], strict=True))
    mem = tu_members(layout)
    return {c: np.array([float(tl[mem[t][0]]) for t in ts]) for c, ts in tu_order(layout).items()}


def lag1_autocorr(x: np.ndarray) -> float:
    """Sample lag-1 autocorrelation (Box-Jenkins estimator)."""
    d = x - x.mean()
    return float(np.sum(d[1:] * d[:-1]) / np.sum(d * d))


def ar1_whiten(x: np.ndarray, rho: float, sigma: float) -> np.ndarray:
    """z_0 = x_0 sqrt(1 - rho^2) / sigma, z_t = (x_t - rho x_{t-1}) / sigma: iid N(0, 1) under the
    spec's stationary AR(1) with innovation SD ``sigma``."""
    z = np.empty_like(x, dtype=float)
    z[0] = x[0] * math.sqrt(1.0 - rho * rho) / sigma
    z[1:] = (x[1:] - rho * x[:-1]) / sigma
    return z


def logistic_fit(X: np.ndarray, y: np.ndarray, iters: int = 50) -> tuple[np.ndarray, np.ndarray]:
    """Maximum-likelihood logistic regression (Newton-Raphson); returns (coef, Wald SE)."""
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        p = sigmoid(X @ b)
        w = p * (1 - p)
        H = X.T @ (X * w[:, None])
        step = np.linalg.solve(H, X.T @ (y - p))
        b = b + step
        if np.max(np.abs(step)) < 1e-10:
            break
    p = sigmoid(X @ b)
    H = X.T @ (X * (p * (1 - p))[:, None])
    return b, np.sqrt(np.diag(np.linalg.inv(H)))


def covariate_matrix(m: pd.DataFrame, extra: Sequence[str] = ()) -> np.ndarray:
    """[1, log length, log1p n_obs, GRAVY, *extra] for rows with an amino-acid sequence."""
    cols = [
        np.ones(len(m)),
        np.log([len(a) for a in m["aa"]]),
        np.log1p([n_obs(a) for a in m["aa"]]),
        np.array([gravy(a) for a in m["aa"]]),
    ]
    cols += [m[c].to_numpy(dtype=float) for c in extra]
    return np.column_stack(cols)
