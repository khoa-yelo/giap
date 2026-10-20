# Autocorrelation of abundance along the replicon

Issue M5-01 adds the shared detrending helper; the autocorrelation statistics are M5-05
(GIAP_SPEC.md §7.3.5(c), ADR-0020). Statistician sign-off for the statistics goes at the end of
this note when M5-05 lands.

## Detrending

`giap.stats.smooth.loess_detrend(y, pos, frac, *, circular, period=None, min_points=3)` and
`loess_matrix(...)` (the n × n smoother H with `loess_detrend(y, …).fit == H @ y`). M5-05 detrends
log2-abundance residuals along gene index and M5-06 detrends Blom scores along TU order, both with
span 10 % and circular padding on circular replicons, so both use this one helper.

- **Estimator.** Plain local-linear LOESS: at each point a weighted least-squares line over its k
  nearest neighbours by `pos`, tricube weights (1 − (d/radius)³)³ with radius the distance to the
  farthest neighbour, and **no robustness iterations** (`it = 0`, lead decision 2). Robust
  reweighting would treat an island as an outlier and keep it out of the trend, so the detrend
  would depend on the signal being tested and the calibration would have to cover that
  data-adaptive step. Plain LOESS is linear in y: one matrix per set of positions, so ADR-0020 can
  apply it to B permuted replicates. The implementation reproduces
  `statsmodels.nonparametric.smoothers_lowess.lowess(y, pos, frac=k/N, it=0, delta=0)` operation
  for operation (≤ 1e-10; neighbourhood ties keep the window on the left, fewer than two weights
  above 1e-12 give fit = y, weighted variance of x floored at 1e-12).
- **Span in points.** k = min(n, max(3, ⌈frac·n − 1e-9⌉)) (the 1e-9 keeps 0.1·30 at 3). `pos` are
  gene indices or TU ordinals, never base pairs; gaps are allowed and respected by the neighbour
  search. ADR-0020 passes ordinals 0 … n − 1 among the genes with values.
- **Circular padding.** On a circular replicon the last k points are copied in front at
  `pos − period` and the first k behind at `pos + period`; the fit at the original points is kept
  and the padded copies feed back into their original columns of H. This removes the end effect
  (a cosine over 100 points, span 10: max residual 0.0071 circular vs 0.0098 linear, the latter at
  the ends) and makes the fit equivariant under rotation of the origin.
- **Small n.** n = 0 gives an empty fit (k = 0). For 1 ≤ n < 3 no local line exists and the fit is
  mean(y) (residuals centred, k = n).
- **The 3-point floor is degenerate** (lead decision 4). With k = 3 the farthest neighbour sits at
  the radius and gets weight 0, and the line through the remaining two points (one is the point
  itself) returns y: fit = y and every residual is 0, for any strictly increasing positions. At the
  10 % span this covers n ≤ 30 values. The detrended series does not exist there, so M5-05 and
  M5-06 report NA with flag `few_values` for any replicon whose k < 4 instead of computing
  statistics over zero residuals. Power by replicon size is recorded in M5-07.
