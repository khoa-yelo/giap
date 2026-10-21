# Islands: expression islands along the replicon

Issues M5-01 (scores, detrending) and M5-03 (scan engine) add the engine helpers; the island step
itself is M5-06 (GIAP_SPEC.md §7.3.5(d)). Each section records its own statistician sign-off; the
block-bootstrap null and its calibration are signed off with M5-06/M5-07.

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

## Scan and maximum statistic (M5-03)

`giap.stats.scan` (§7.3.5(d) steps 4–6, §8.3). **Windows** are lengths w = wmin … min(wmax, cap)
with cap = n − 1 on circular replicons (no window covers the whole circle; after detrending such a
window has S ≈ 0) and n on linear ones. Circular windows start at every position and wrap. They are
ordered by length, then start. S_w = Σ_{t∈w} z_t / √n_w comes from prefix sums. The observed series
and the null rows use the same kernel, so a null row equal to the observed series reproduces its
maximum exactly.

**Family and FWER.** One family is all windows of one replicon, both signs. The statistic is M =
max_w |S_w|, and p_w = (1 + #{M_b ≥ |S_w|})/(1 + B) (single-step maxT, Westfall–Young), so high and
low islands share one α. Ties within 1e-12·max(1, |S|) count as exceedances (conservative). If the
null rows are exchangeable with the observed series under the complete null, P(min_w p_w ≤ α) ≤ α
exactly. The excess over α is discreteness only, at most 1/(B + 1). With true islands present, the
signal is spread into the null rows (permutation or block bootstrap), which raises M_b. The
procedure is then conservative: in simulation with planted islands of 5–60 TUs the false-island rate
was 0–0.033 at α = 0.05 (PR #118 review). Being single-step, it is also conservative for secondary
islands. The engine is only as valid as the exchangeability of its null rows: iid permutations are
not valid under autocorrelation. The production null is the circular block bootstrap with ℓ′
(ADR-0021), calibrated in M5-06/M5-07. FWER is per replicon; across k replicons the genome-wide FWER
is up to 1 − (1 − α)^k. A replicon with no window (n < wmin) has `max_abs` NaN and p NaN, which is
never significant. A non-finite |S| gives NaN; a negative finite |S| is a caller error and
raises (decision 22).

**Why S/√n_w is not standardised for TU size.** z_t is the mean of the gene scores of TU t, so Var
z_t falls with the TU's gene count. S_w/√n_w has unit variance only when every TU has unit variance.
Windows of single-gene TUs have the largest variance and dominate M_b. Windows over large operons
then need a larger true shift to exceed the same threshold, so power is lower in operon-dense
regions. This is a power cost, not a size cost, as long as the null rows keep the positional mix of
TU sizes (the block bootstrap does; iid permutation does not). TU-size-standardised scores were
evaluated and rejected (ADR-0021: FWER 0.18 at ρ = 0).

**Merging and the representative.** Significant windows of the same direction that share ≥ 1 TU
(cyclically on circular replicons) are linked. Islands are the connected components; touching
windows stay separate. A window with S = 0 has no direction and is dropped (it cannot be
significant, because M_b ≥ 0 gives p = 1). The best window has the largest |S|, then fewer TUs, then
the smaller start. Island-level guarantee: a false island (no window overlapping a true island)
contains a falsely rejected window, so P(any false island) ≤ the window FWER. Not controlled: the
island's extent. Boundary-straddling windows are non-null and significant, so the union of a cluster
usually extends past the true island (27–100 % of replicates in simulation). The reported span is
therefore the best window's span (≤ wmax = 20 TUs), not the union. An island longer than 20 TUs is
represented by a best window of ≤ 20 TUs and can split into two islands; M5-10's Jaccard inherits
both effects. S and the effect of the best window are selected (winner's curse; see ADR-0021 on the
effect CI).

Signed off by the statistician (PR #118 review).
