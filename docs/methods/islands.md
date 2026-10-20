# Islands: expression islands along the replicon

Issue M5-01 adds the engine helpers; the island scan itself is M5-06 (GIAP_SPEC.md §7.3.5(d)).
Statistician sign-off for the scan goes at the end of this note when M5-06 lands.

## Gene scores: left-censored Blom

`giap.stats.scores.blom_scores(x, censored_mask)` (§7.3.5(d) step 1, §8.3). Within a replicon the
eligible genes are ranked by log2 abundance. Undetected genes are **left-censored**: they lie below
the detection limit, so their order among themselves is unknown, but their position (below every
detected gene) is known. All n₀ of them therefore share the average of ranks 1 … n₀,
r = (n₀ + 1)/2, and a detected gene gets n₀ + its average rank among the detected values; a
detected gene with a very low abundance still ranks above every undetected one (censoring is "at
the lowest rank", not "at its own value"). Ties among detected values and the censored block are
averaged separately. The score is the Blom normal score Φ⁻¹((r − 3/8)/(n + 1/4)) (`scipy.special.ndtri`).

Properties: scores are invariant under strictly increasing transforms of the values and
equivariant under permutation; without ties or censoring the scores of ranks r and n + 1 − r sum
to 0; a fully censored replicon scores exactly 0 everywhere (r = (n + 1)/2 maps to Φ⁻¹(0.5)). An
uncensored NaN or ±inf raises `ValueError`: the caller decides who is scored.

The scores are detrended with the shared LOESS helper (docs/methods/autocorrelation.md,
"Detrending") along TU order, span 10 % of TUs.
