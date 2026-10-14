"""De novo precision per score bin, estimated from the genome (GIAP_SPEC.md §7.2.5, §8.1, §9.3).

For score bin b and length class l: m_bl = share of de novo peptides matching the target
translation exactly, r_bl = share matching the shuffled decoy exactly, c_l = share of reference
peptides (q <= 0.01, targets) matching the target exactly.

    prec_b = sum_l w_bl (m_bl - r_bl) / (c_l - r_bl),  clipped to [0, 1]

w_bl is the share of bin-b peptides in class l. Classes with c_l - r_bl < 0.05 (or without c_l)
are dropped from the bin and the remaining weights renormalised; a bin with no class left is
unreported (precision NA, π̂ = 1). The 95 % CI is the percentile interval of ``n_boot`` resamples of
the bin's peptides with replacement, c_l held fixed. The de novo π̂_wrong handed to M3-10 is the
conservative bound 1 - ci_low (ADR-0010). Method: docs/methods/denovo-precision.md.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

C_DEFAULT = 0.99
MIN_DENOMINATOR = 0.05
REF_MAX_Q = 0.01
_EPS = 1e-12

__all__ = [
    "BinPrecision",
    "c_by_class",
    "class_label",
    "denovo_precision",
    "pi_wrong_denovo",
    "score_bins",
]


@dataclass(frozen=True, slots=True)
class BinPrecision:
    """De novo precision of one score bin (``mapping_fdr.tsv``, section ``denovo_precision``).

    ``c_value`` is Σ_l w_bl c_l over the classes used for the estimate (NaN when no class has c_l).
    """

    score_bin: int
    n: int
    precision: float | None
    ci_low: float | None
    ci_high: float | None
    c_value: float
    flags: tuple[str, ...]


def class_label(bounds: Sequence[int]) -> str:
    """Label of a length class ``[lo, hi]``: ``"8"`` or ``"10-11"`` (M3-10 labels)."""
    lo, hi = int(bounds[0]), int(bounds[1])
    return str(lo) if lo == hi else f"{lo}-{hi}"


def score_bins(scores: np.ndarray, n_bins: int = 10) -> np.ndarray:
    """Quantile bin of each score: edges ``np.quantile(scores, [1/n, ..., (n-1)/n])``, right-closed.

    Coinciding edges merge bins (fewer than ``n_bins``; tied scores always share a bin). Scores that
    are NaN get bin -1 (no bin).
    """
    s = np.asarray(scores, dtype=float)
    out = np.full(s.shape, -1, dtype=np.int64)
    ok = ~np.isnan(s)
    if not ok.any():
        return out
    edges = np.unique(np.quantile(s[ok], np.arange(1, n_bins) / n_bins))
    out[ok] = np.searchsorted(edges, s[ok], side="left")
    return out


def c_by_class(
    ref_exact: pd.DataFrame | None, classes: Sequence[Sequence[int]]
) -> tuple[dict[str, float], bool]:
    """c_l per length class from the reference peptides, and whether c was assumed.

    Reads ``peptide``, ``q_value``, ``is_decoy`` and ``target_exact``; only targets with
    q <= 0.01 count, and the class comes from the peptide length. Without usable reference
    peptides every c_l = 0.99 and the flag ``c_assumed`` applies (returned True). Otherwise a class
    without reference peptides gets c_l = NaN (its bins are not reportable).
    """
    labels = [class_label(b) for b in classes]
    if ref_exact is not None and len(ref_exact) > 0:
        keep = (ref_exact["q_value"].to_numpy(dtype=float) <= REF_MAX_Q) & ~ref_exact[
            "is_decoy"
        ].to_numpy(dtype=bool)
        ref = ref_exact.loc[keep]
        if len(ref) > 0:
            lengths = ref["peptide"].astype(str).str.len().to_numpy()
            exact = ref["target_exact"].to_numpy(dtype=bool)
            c: dict[str, float] = {}
            for lab, b in zip(labels, classes, strict=True):
                sel = (lengths >= int(b[0])) & (lengths <= int(b[1]))
                c[lab] = float(exact[sel].mean()) if sel.any() else math.nan
            return c, False
    return dict.fromkeys(labels, C_DEFAULT), True


def _estimate(n_cls: np.ndarray, m_cnt: np.ndarray, r_cnt: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Estimator over the last axis (classes); leading axes are resamples. NaN = not reportable.

    ``n_cls``, ``m_cnt``, ``r_cnt`` are per-class counts of peptides, target and decoy matches.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        m = m_cnt / n_cls
        r = r_cnt / n_cls
        denom = c - r
        use = (n_cls > 0) & ~np.isnan(c) & (denom >= MIN_DENOMINATOR - _EPS)
        w = np.where(use, n_cls, 0.0)
        tot = w.sum(axis=-1)
        term = np.where(use, (m - r) / np.where(use, denom, 1.0), 0.0)
        raw = (w * term).sum(axis=-1) / tot
    out: np.ndarray = np.where(tot > 0, np.clip(raw, 0.0, 1.0), np.nan)
    return out


def _c_value(n_cls: np.ndarray, m_cnt: np.ndarray, r_cnt: np.ndarray, c: np.ndarray) -> float:
    """Σ_l w_bl c_l over the classes used in the estimate (all classes with c_l if none is)."""
    r = r_cnt / n_cls
    known = ~np.isnan(c)
    use = known & (c - r >= MIN_DENOMINATOR - _EPS)
    sel = use if use.any() else known
    if not sel.any():
        return math.nan
    return float(np.sum(n_cls[sel] * c[sel]) / np.sum(n_cls[sel]))


def _bin_row(
    sub: pd.DataFrame,
    score_bin: int,
    c: Mapping[str, float],
    flags0: tuple[str, ...],
    rng: np.random.Generator,
    n_boot: int,
) -> BinPrecision:
    """Point estimate and percentile bootstrap CI for one bin."""
    cls = sub["length_class"].astype(str).to_numpy()
    tgt = sub["target_exact"].to_numpy(dtype=bool)
    dec = sub["decoy_exact"].to_numpy(dtype=bool)
    labels = sorted(set(cls.tolist()))
    k = len(labels)
    cvec = np.array([float(c.get(lab, math.nan)) for lab in labels], dtype=float)
    # cell counts: class x (target, decoy) -> 4k categories, independent of row order
    idx = np.searchsorted(np.array(labels), cls) * 4 + tgt.astype(int) * 2 + dec.astype(int)
    cells = np.bincount(idx, minlength=4 * k).astype(float).reshape(k, 4)

    def summarise(cl: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return cl.sum(axis=-1), cl[..., 2] + cl[..., 3], cl[..., 1] + cl[..., 3]

    n_cls, m_cnt, r_cnt = summarise(cells)
    point = float(_estimate(n_cls, m_cnt, r_cnt, cvec))
    c_value = _c_value(n_cls, m_cnt, r_cnt, cvec)
    n = len(sub)
    if math.isnan(point):
        return BinPrecision(score_bin, n, None, None, None, c_value, (*flags0, "unreported"))
    # resampling n peptides with replacement == multinomial draw over the cells
    boot_cells = rng.multinomial(n, cells.ravel() / n, size=n_boot).astype(float)
    boot = _estimate(*summarise(boot_cells.reshape(n_boot, k, 4)), cvec)
    boot = boot[~np.isnan(boot)]
    if boot.size == 0:
        lo = hi = point
    else:
        lo, hi = (float(v) for v in np.quantile(boot, [0.025, 0.975]))
    return BinPrecision(score_bin, n, point, lo, hi, c_value, flags0)


def denovo_precision(
    df: pd.DataFrame,
    c: Mapping[str, float],
    c_assumed: bool,
    rng: np.random.Generator,
    n_boot: int = 1000,
) -> list[BinPrecision]:
    """Precision per score bin (deciles of ``score``), in bin order.

    ``df`` columns: ``peptide``, ``score``, ``length_class``, ``target_exact``, ``decoy_exact``
    (exact matches only, I = L). ``c`` maps length-class label to c_l (from :func:`c_by_class`);
    ``rng`` should be ``giap.rng.for_step("placement", seed)``. Unscored peptides form no bin.
    """
    if len(df) == 0:
        return []
    bins = score_bins(df["score"].to_numpy(dtype=float))
    flags0: tuple[str, ...] = ("c_assumed",) if c_assumed else ()
    rows: list[BinPrecision] = []
    for b in np.unique(bins[bins >= 0]):
        sub = df.loc[bins == b]
        rows.append(_bin_row(sub, int(b), c, flags0, rng, n_boot))
    return rows


def pi_wrong_denovo(bins: Sequence[BinPrecision], bin_of_peptide: np.ndarray) -> np.ndarray:
    """π̂_wrong per de novo peptide: 1 - ci_low of its bin (ADR-0010); 1 when not estimable."""
    lookup = {
        b.score_bin: min(1.0, max(0.0, 1.0 - b.ci_low))
        for b in bins
        if b.precision is not None and b.ci_low is not None
    }
    idx = np.asarray(bin_of_peptide, dtype=np.int64)
    return np.array([lookup.get(int(i), 1.0) for i in idx], dtype=float)
