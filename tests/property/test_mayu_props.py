"""Hypothesis properties of the MAYU closed form (issue M4-02; GIAP_SPEC.md §7.3.2, §8.3).

- the curve equals an independent reference of §7.3.2 given the function's own bins;
- FDR and fdr_2pep in [0, 1]; fdr_2pep <= fdr; gene_q in [0, 1] or NaN exactly where best_q is NaN;
- gene_q non-decreasing in best_q and equal to monotone_q of the curve;
- adding false-peptide mass (PEP) at or below every threshold never lowers FDR;
- bins are equal-count rank bins (counts differ by <= 1, ordered by a_g, ties by input order).
"""

from __future__ import annotations

import math

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.stats.mayu import MayuResult, mayu_gene_fdr
from giap.stats.tda import monotone_q

FAST = settings(max_examples=150, deadline=None)
NAN = float("nan")

Q = st.sampled_from([0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05])  # few values: frequent ties
Q_OR_NAN = st.one_of(Q, st.just(NAN))


@st.composite
def inputs(draw: st.DrawFn, with_pep: bool | None = None) -> dict[str, object]:
    n = draw(st.integers(1, 40))
    best_q = draw(st.lists(Q_OR_NAN, min_size=n, max_size=n))
    sizes = draw(st.lists(st.integers(0, 60), min_size=n, max_size=n))
    m = draw(st.integers(0, 60))
    peptide_q = draw(st.lists(Q, min_size=m, max_size=m))
    use_pep = draw(st.booleans()) if with_pep is None else with_pep
    pep = draw(st.lists(st.floats(0.0, 1.0), min_size=m, max_size=m)) if use_pep else None
    n_bins = draw(st.integers(1, 12))
    return {
        "best_q": np.array(best_q),
        "peptide_q": np.array(peptide_q, dtype=float),
        "peptide_pep": None if pep is None else np.array(pep),
        "sizes": np.array(sizes, dtype=float),
        "n_bins": n_bins,
    }


def _run(x: dict[str, object]) -> MayuResult:
    return mayu_gene_fdr(x["best_q"], x["peptide_q"], x["peptide_pep"], x["sizes"], x["n_bins"])  # type: ignore[arg-type]


def _reference(x: dict[str, object], bins: np.ndarray) -> list[tuple[float, float, float]]:
    """(t, fdr, fdr_2pep) straight from §7.3.2 using the given bin labels."""
    best_q = np.asarray(x["best_q"])
    pq = np.asarray(x["peptide_q"])
    pep = x["peptide_pep"]
    a = np.maximum(np.asarray(x["sizes"], dtype=float), 1.0)
    labels = np.unique(bins)
    w = {b: a[bins == b].mean() / a.sum() for b in labels}
    out = []
    for t in np.unique(best_q[np.isfinite(best_q)]):
        inp = pq <= t
        F = float(np.asarray(pep)[inp].sum()) if pep is not None else float(t * inp.sum())
        fp = fp2 = 0.0
        n_tot = 0
        for b in labels:
            Nb = int((bins == b).sum())
            nb = int(((bins == b) & (best_q <= t)).sum())
            n_tot += nb
            p = 1.0 - (1.0 - w[b]) ** F
            lam = F * w[b]
            p2 = 1.0 - math.exp(-lam) * (1.0 + lam)
            fp += (Nb - nb) if p >= 1.0 else (Nb - nb) * p / (1.0 - p)
            fp2 += (Nb - nb) if p2 >= 1.0 else (Nb - nb) * p2 / (1.0 - p2)
        out.append((float(t), min(1.0, fp / n_tot), min(1.0, fp2 / n_tot)))
    return out


@FAST
@given(inputs())
def test_curve_matches_reference(x: dict[str, object]) -> None:
    res = _run(x)
    ref = _reference(x, res.bin_of_gene)
    curve = res.curve
    assert len(curve) == len(ref)
    np.testing.assert_allclose(
        np.asarray(curve["t"], dtype=float), [r[0] for r in ref], rtol=0, atol=0
    )
    np.testing.assert_allclose(
        np.asarray(curve["fdr"], dtype=float), [r[1] for r in ref], rtol=1e-7, atol=1e-12
    )
    np.testing.assert_allclose(
        np.asarray(curve["fdr_2pep"], dtype=float), [r[2] for r in ref], rtol=1e-6, atol=1e-12
    )


@FAST
@given(inputs())
def test_estimates_in_unit_interval(x: dict[str, object]) -> None:
    res = _run(x)
    fdr = res.curve["fdr"].to_numpy()
    fdr2 = res.curve["fdr_2pep"].to_numpy()
    assert np.all((fdr >= 0) & (fdr <= 1))
    assert np.all((fdr2 >= 0) & (fdr2 <= 1))
    assert np.all(fdr2 <= fdr + 1e-12)
    assert np.all(res.curve["fp"].to_numpy() >= 0)
    best_q = np.asarray(x["best_q"])
    assert np.array_equal(np.isnan(res.gene_q), np.isnan(best_q))
    fin = res.gene_q[~np.isnan(res.gene_q)]
    assert np.all((fin >= 0) & (fin <= 1))


@FAST
@given(inputs())
def test_gene_q_monotone_in_best_q(x: dict[str, object]) -> None:
    res = _run(x)
    best_q = np.asarray(x["best_q"])
    ok = np.isfinite(best_q)
    order = np.argsort(best_q[ok], kind="stable")
    assert np.all(np.diff(res.gene_q[ok][order]) >= -1e-15)
    # equal best_q -> equal gene_q
    for t in np.unique(best_q[ok]):
        assert np.unique(res.gene_q[best_q == t]).size == 1


@FAST
@given(inputs())
def test_gene_q_is_monotone_q_of_curve(x: dict[str, object]) -> None:
    res = _run(x)
    curve = res.curve
    q_at = dict(zip(curve["t"].tolist(), monotone_q(curve["fdr"].to_numpy()).tolist(), strict=True))
    best_q = np.asarray(x["best_q"])
    for g in np.flatnonzero(np.isfinite(best_q)):
        assert res.gene_q[g] == q_at[float(best_q[g])]


@FAST
@given(inputs(with_pep=True), st.lists(st.floats(0.0, 1.0), min_size=1, max_size=20))
def test_more_false_mass_never_lowers_fdr(x: dict[str, object], extra: list[float]) -> None:
    base = _run(x)
    y = dict(x)
    y["peptide_q"] = np.r_[np.asarray(x["peptide_q"]), np.full(len(extra), 0.0001)]
    y["peptide_pep"] = np.r_[np.asarray(x["peptide_pep"]), extra]
    more = _run(y)
    assert np.all(more.curve["F"].to_numpy() >= base.curve["F"].to_numpy() - 1e-12)
    assert np.all(more.curve["fdr"].to_numpy() >= base.curve["fdr"].to_numpy() - 1e-12)
    assert np.all(
        more.gene_q[np.isfinite(base.gene_q)] >= base.gene_q[np.isfinite(base.gene_q)] - 1e-12
    )


@FAST
@given(inputs())
def test_equal_count_rank_bins(x: dict[str, object]) -> None:
    res = _run(x)
    b = res.bin_of_gene
    n = len(np.asarray(x["best_q"]))
    k = min(int(x["n_bins"]), n)  # type: ignore[call-overload]
    assert b.dtype == np.int64 and b.shape == (n,)
    counts = np.bincount(b, minlength=k)
    assert len(counts) == k and counts.min() >= 1 and counts.max() - counts.min() <= 1
    a = np.maximum(np.asarray(x["sizes"]), 1.0)
    order = np.argsort(a, kind="stable")
    assert np.all(np.diff(b[order]) >= 0)


@FAST
@given(inputs())
def test_deterministic(x: dict[str, object]) -> None:
    r1, r2 = _run(x), _run(x)
    assert r1.curve.equals(r2.curve)
    np.testing.assert_array_equal(r1.gene_q, r2.gene_q)
    np.testing.assert_array_equal(r1.bin_of_gene, r2.bin_of_gene)
