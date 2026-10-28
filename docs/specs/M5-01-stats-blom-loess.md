# M5-01: Stats engine: Blom normal scores with left-censoring, LOESS detrending (circular padding)

Milestone M5 · track S (engine) · depends on: none · estimate ≈ 200 code lines (+ ≈ 250 test
lines) · **statistical** (statistician reviews; label `stats`) · **can start now** · no tools, no
Slurm

## Spec sections implemented

- §8.3 "`blom_scores(x, censored_mask)` | ties → average ranks".
- §7.3.5(d) step 1 quoted: "Within each replicon, rank eligible genes by log2 abundance, with
  undetected genes tied at the lowest rank (left-censored). Convert to Blom normal scores
  Φ⁻¹((r − 3/8)/(n + 1/4)), using average ranks for ties."
- §7.3.5(c) "Values. x = log2-abundance residuals of detected genes, after removing the replicon mean
  and a wide LOESS trend along gene index (span = 10 % of genes, circular padding on circular
  replicons)" and §7.3.5(d) step 3 "Detrend with LOESS along TU order (span 10 % of TUs, circular
  padding on circular replicons)". The spec names no engine function for LOESS; this issue adds one
  shared helper so M5-05 and M5-06 detrend identically.
- §7.3.6 unit test "Blom scores with ties".

## Public interface

`src/gippa/stats/scores.py` (pure)
```python
def blom_scores(x: ArrayLike, censored_mask: ArrayLike) -> npt.NDArray[np.float64]: ...
```
- n = len(x). Censored entries (mask True) all get the average of ranks 1 … n₀ (n₀ = number censored),
  i.e. (n₀ + 1)/2; uncensored entries get n₀ + their average rank among the uncensored values.
  Score = Φ⁻¹((r − 3/8)/(n + 1/4)) (`scipy.special.ndtri`). Output in input order.
- The value of a censored entry is ignored (it may be NaN). An uncensored entry that is NaN or ±inf
  raises `ValueError` (the caller decides who is scored; see M5-06).
- Lengths of `x` and `censored_mask` differ → `ValueError`. n = 0 → empty array.

`src/gippa/stats/smooth.py` (pure)
```python
@dataclass(frozen=True, slots=True)
class LoessFit:
    fit: npt.NDArray[np.float64]       # fitted trend at each input point
    resid: npt.NDArray[np.float64]     # y - fit
    k: int                             # neighbours used per local fit (0 when n = 0)

def loess_span_points(n: int, frac: float, min_points: int = 3) -> int: ...
def loess_detrend(y: ArrayLike, pos: ArrayLike, frac: float, *, circular: bool,
                  period: float | None = None, min_points: int = 3) -> LoessFit: ...
```
- **Span:** k = min(n, max(`min_points`, ⌈frac·n − 1e-9⌉)) (the 1e-9 stops 0.1·30 from rounding up
  to 4).
- **Fit:** local linear regression with tricube weights over the k nearest points by `pos` and no
  robustness iterations: numerically equal (≤ 1e-10 absolute) to
  `statsmodels.nonparametric.smoothers_lowess.lowess(y_pad, pos_pad, frac=k/n_pad, it=0, delta=0,
  return_sorted=False)` restricted to the original points. Calling statsmodels (BSD-3) is allowed;
  a reimplementation must match it.
- **Circular padding** (`circular=True`, `period` required): prepend the last min(k, n) points with
  `pos − period` and append the first min(k, n) points with `pos + period`; linear: no padding.
- **Degenerate inputs (the estimator does not exist):** n = 0 → empty fit, k = 0; 1 ≤ n < `min_points`
  → `fit` = mean(y) everywhere (residuals centred), k = n; all y equal → residuals 0 (the local fits
  exist). `pos` not strictly increasing, any non-finite `y`/`pos`, `frac` ∉ (0, 1], or (circular and
  `pos` outside [0, `period`)) → `ValueError`.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `blom_scores([5, nan, 3, nan], [F, T, F, T])` | ranks [4, 1.5, 3, 1.5]; scores [1.0491, −0.6289, 0.2993, −0.6289] (4 dp) |
| 2 | `blom_scores([2, 2, 1], [F, F, F])` | ranks [2.5, 2.5, 1]; scores [0.3957, 0.3957, −0.8694] |
| 3 | all four censored | all scores exactly 0.0 (r = 2.5 = (n + 1)/2 maps to Φ⁻¹(0.5)) |
| 4 | n = 1, uncensored | [0.0] |
| 5 | n = 0 | empty array |
| 6 | uncensored entry NaN | `ValueError` |
| 7 | `loess_span_points(30, 0.10)` / `(25, 0.10)` / `(10, 0.10)` | 3 / 3 / 3 (min_points) ; `(1000, 0.10)` → 100 |
| 8 | y = 2·pos + 1, pos = 0 … 99, frac 0.10, linear | max abs residual < 1e-9 |
| 9 | y = cos(2π·pos/100), pos = 0 … 99, frac 0.10 | circular: max abs residual 0.00711 ± 1e-5; linear: 0.00977 ± 1e-5 (largest at the ends; padding removes the end effect) |
| 10 | n = 2, y = [1, 3] | fit [2, 2], resid [−1, 1], k = 2 |
| 11 | y constant 5.0, n = 50 | resid all 0 (≤ 1e-12) |

## Edge cases

- Blom: ties among detected values and the censored block are averaged separately; a detected gene
  whose abundance is lower than nothing else still ranks above every censored gene (censoring is "at
  the lowest rank", not "at its own value").
- Blom scores are symmetric: with no ties and no censoring, scores of ranks r and n + 1 − r sum to 0.
- LOESS on a circular replicon with k = n (tiny replicon): padding duplicates every point once on each
  side; the result must still equal statsmodels on the padded series.
- `pos` are gene indices (M5-05) or TU ordinals (M5-06), never base pairs; gaps in `pos` (undetected
  genes) are allowed and are respected by the nearest-neighbour search.

## Acceptance criteria

- [ ] `tests/unit/stats/test_scores.py`: examples 1–6; symmetry property (Hypothesis, no ties).
- [ ] `tests/unit/stats/test_smooth.py`: examples 7–11; equality with statsmodels on random inputs
  (linear and circular, n ∈ [3, 300], frac ∈ [0.05, 0.5]).
- [ ] Property (`tests/property/test_smooth_props.py`): **rotation equivariance** on circular input:
  rotating (y, pos) by any offset (pos taken mod period, re-sorted) rotates the residuals (≤ 1e-9);
  the linear fit reproduces any affine y exactly (≤ 1e-9).
- [ ] mypy and ruff clean; no randomness in either function.

## Methods note (statistician sign-off)

Sections in `docs/methods/islands.md` ("Gene scores: left-censored Blom") and
`docs/methods/autocorrelation.md` ("Detrending"): why left-censored ranks (undetected genes are
below the detection limit, so their order is unknown but their position is not), the choice it = 0
(no robustness iterations; see Q2), span in points, circular padding, the n < 3 rule.

## Questions for the lead

1. Is a separate engine helper for LOESS acceptable (§8.3 does not list one)? **Proposed default:** yes,
   `gippa.stats.smooth.loess_detrend`; it is an internal helper, not a new interface of the CLI or
   tables, so no ADR.
2. Robust LOESS (statsmodels default `it=3`) or plain local linear (`it=0`)? Robustness iterations
   would treat islands as outliers and keep them out of the trend (more power), but they make the
   detrend data-adaptive in a way the calibration must cover. **Proposed default:** `it=0` (what
   "LOESS" means in Cleveland's original definition without the robust step), revisit only if the M5-07
   power curve is poor.
3. Minimum span of 3 points: on small replicons (a 60-gene plasmid has ≈ 33 TUs → 10 % = 3 TUs) the
   LOESS absorbs almost all structure, so islands cannot be found there. **Proposed default:** keep the
   spec's 10 % with the 3-point floor (valid, low power on small replicons); record power by replicon
   size in M5-07.

## Lead decisions (2026-10-20, after the statistician's attack check)

1. Q1: yes — `gippa.stats.smooth.loess_detrend` as an internal engine helper, no ADR. It accepts positions
   0…n−1 (the ordinal among genes with values), so the smoother is one matrix for every permutation
   (ADR-0020 needs this; TU-block permutations move the valued positions).
2. Q2: plain local-linear LOESS (`it = 0`); robust iterations would make the trend depend on the island.
3. Q3: keep the 10 % span with the 3-point floor; record power by replicon size (the span costs about half
   the power at 1,000 genes; at 4,500 genes it is ≈ 255 TUs, far longer than an island). A plasmid floor
   is a later ADR candidate.

## Lead decisions after the test-author's report (2026-10-20)
4. The 3-point floor stays (statsmodels equality, examples 7 and 10). For any strictly increasing positions a k = 3
   window gives weight 0 to both neighbours, so the fit equals y and every residual is 0: the detrended
   series does not exist. Methods note says so. Consumers decide: M5-05 and M5-06 report NA with flag
   `few_values` for any replicon whose span gives k < 4 (≈ n ≤ 30 ordinals at the 10 % span), instead of
   computing statistics over zero residuals. (Recorded in M5-05 decision 21 and M5-06 decision 40.)
5. `loess_detrend` raises `ValueError` when `y` and `pos` differ in length.
6. Add `loess_matrix(pos, frac, *, circular, period=None, min_points=3) -> NDArray[float64]` (n × n,
   `loess_detrend(y, …).fit == loess_matrix(…) @ y`), so ADR-0020 applies one matrix to B permuted
   replicates; `loess_detrend` stays 1-D.
7. No calibration test for M5-01 (deterministic helpers, no statistical claim).
8. (test-author follow-up) `loess_matrix` validates exactly as `loess_detrend` (ValueError on tied or
   unsorted pos, NaN, circular without period, pos outside [0, period), frac outside (0, 1]); Hypothesis
   tolerance 1e-12 × max(1, max|y|) accepted; 2-D `y` to `loess_detrend` is unspecified (no test).
