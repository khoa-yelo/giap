# Islands: expression islands along the replicon

Issues M5-01 (scores, detrending) and M5-03 (scan engine) add the engine helpers; the island step
itself is M5-06 (GIAP_SPEC.md §7.3.5(d)). Each section records its own statistician sign-off; the
block-bootstrap null and its calibration are signed off with M5-06/M5-07.

## Report text (M9-16)

The block below is copied into the run report's methods section by
`scripts/dev/sync_report_snippets.py` (plain paragraphs and lists only).

<!-- report:begin -->
Islands are runs of adjacent transcription units whose detrended abundance scores are jointly high
or low along a replicon. Each window sum is compared with the maximum of a TU-block bootstrap of
the same scores, taken over every scanned replicon, so one family-wise error rate covers the whole
run; a window is reported when its p-value is at most stats.alpha.context.

The FWER is claimed for a dependence regime, not for a replicon size. The calibration supports
it for TU-level autocorrelation up to rho = 0.5 at every replicon size measured, and up to
rho = 0.9 at 4,500 genes; there is no evidence beyond that (rho above 0.95, or dependence longer
than the 20-TU bootstrap block cap). It is not claimed where strong TU-level autocorrelation is
resolved by Moran's I only to a short lag (block length 7 or less): a high Moran's I with
block_len 7 or less is outside the claimed domain. No per-replicon statistic certifies the domain, so nothing is flagged; instead the report
lists, per replicon:

- the estimated correlation length (l-hat), or why it was not estimated;
- the bootstrap block length used (block_len);
- the lag-1 Moran's I of the detrended TU scores.
<!-- report:end -->

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
procedure is then conservative: in simulation with planted high islands of 5–60 TUs the
false-island rate was 0–0.033 at α = 0.05 (PR #118 review). A planted low (silenced) island, which
the step cannot detect, is not conservative in this sense: it is a low-variance stretch at the
censored atom, not extra variability spread into the null rows. It gives 0.034–0.036 at R = 500,
against 0.020 on the complete null ("Calibration and power" below), which is within α. Being single-step, it is also conservative for secondary
islands. The engine is only as valid as the exchangeability of its null rows: iid permutations are
not valid under autocorrelation. The production null is the circular block bootstrap with ℓ′
(ADR-0021), calibrated in M5-06/M5-07. **Which family applies where** (M5-06 decision 46): the
engine used alone, one call per replicon against that replicon's own M_b, controls the FWER per
replicon, and across k replicons the genome-wide rate is then up to 1 − (1 − α)^k. The M5-06 island
step does **not** use it that way — it pools one genome-wide M_b over every scanned replicon, so one
FWER covers the run ("The null and the genome-wide family" below). A replicon with no window
(n < wmin) has `max_abs` NaN and p NaN, which is never significant. A non-finite |S| gives NaN; a
negative finite |S| is a caller error and raises (decision 22).

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

## The island step (M5-06)

`giap.context.islands` (§7.3.5(d) steps 1–10, §8.1 row "Islands", §9.3 `tables/islands.tsv`).
This section is the claim the island table makes; the two sections above are the pieces it is
built from.

### Scores, TU series, detrending

Per replicon: the eligible genes (not partial, `p_detect` present) are Blom-scored with the
undetected ones left-censored ("Gene scores" above); z_t is the mean score of the scored genes of
TU t; the series is detrended by LOESS over **TU ordinals** (not base pairs) with span
k = ⌈0.1·T⌉ points and circular padding on circular replicons.

- **Missing-value policy** (decision 22). A detected gene whose abundance is NaN — possible only
  under `abundance.missing = exclude` — is *not scored*: it is detected, so censoring it at the
  bottom would be a lie, and its rank is unknown. It therefore leaves the ranking entirely, which
  changes every other gene's Blom score slightly; a TU with no scored gene left is dropped from the
  series. Under `abundance.missing = min` the same gene is scored at the imputed minimum. The two
  policies are the §8.6 `basic` sensitivity axis, so the island *membership* must not move with the
  policy for its own reasons — see "What is reported, and what the gene list is" below.
- **The detrend attenuates the signal it is scanning for.** A boxcar island of w TUs keeps about
  1 − 1.7·w/k of its amplitude while w ≲ k/3 and far less beyond. Measured retained amplitude
  (circular LOESS, span 10 %):

  | T (scored TUs) | k | w = 2 | w = 5 | w = 10 | w = 20 |
  |---|---|---|---|---|---|
  | 2,500 | 250 | 0.99 | 0.97 | 0.93 | 0.86 |
  | 1,000 | 100 | 0.97 | 0.91 | 0.83 | 0.66 |
  | 400 | 40 | 0.91 | 0.79 | 0.58 | 0.32 |
  | 200 | 20 | 0.83 | 0.58 | 0.32 | 0.16 |

  Since w ≤ `max_tus` = 20, this is a few per cent on a chromosome and a near-total loss of the long
  windows on a replicon of a few hundred TUs: on small replicons the step effectively scans only its
  short windows. It is a power cost, not a level cost — the null rows are detrended by the same
  linear smoother.

### The null and the genome-wide family

Each scanned replicon is resampled by `circular_block_bootstrap(detrended, ℓ′, B, rng)` with
ℓ′ = min(20, 2ℓ̂ + 1, T) per replicon (ADR-0021 decision 1; the rule, the evidence for ℓ′ over ℓ̂
and the join bias live in docs/methods/bootstrap.md "Block length"). The null keeps what the data
have: whole TU scores move together, so the TU-size mix and the detectability pattern travel with
the scores, and blocks of ℓ′ consecutive TUs keep the local autocorrelation. It is a *stationary*
null: it does not keep variation of the local scale along the replicon, which is the source of both
breaches recorded below.

**One family, genome-wide** (decisions 21 and 46). For replicate b, M_b is the maximum of |S_w|
over the windows of **every** scanned replicon, taken element-wise across replicons at the same b;
the replicons are resampled independently, which is the complete null (replicons are independent
under it), and the pairing at index b is arbitrary but valid. p_w = (1 + #{M_b ≥ |S_w|})/(1 + B)
and a window is reported when p_w ≤ `stats.alpha.context`. So **one FWER of 0.05 covers the whole
run**, matching the single §8.1 row, instead of 1 − (1 − α)^k over k replicons. The price is paid by
the smaller replicons: a 60-TU plasmid raises the common threshold a little and loses power to the
chromosome's maximum, which is the correct direction for a run-level claim. A replicon that is not
scanned (below) contributes nothing to M_b, consumes no draws and cannot carry an island.

### Validity domain (ADR-0021 Amendment 1)

**The FWER is claimed for a dependence regime, not for a replicon size.** It is claimed for giapsim
`ar1_background` at ρ ≤ 0.5 at every size measured (557 / 1,001 / 2,501 scored TUs: 0.020 / 0.031 /
0.037 at α = 0.05, R = 1,000, bound 0.0638), for the 4,500-gene chromosome at ρ = 0.9 (0.024, where
ℓ̂ ≈ 8 makes ℓ′ = 17 span the 20-TU windows) and for the §8.5 attacks that stay in that regime.

It is **not** claimed where strong TU-level autocorrelation is resolved by Moran's I only to a short
lag (ℓ′ ≤ 7): ρ = 0.9 gives 0.113 at 557 TUs and 0.095 at 1,001 TUs, and 1,000 single-gene TUs give
0.068 / 0.110 / 0.208 / 0.327 at ρ = 0.7 / 0.8 / 0.9 / 0.95. The scored-TU count decides nothing on
its own (400 TUs at ρ 0.5: 0.005; 1,000 TUs at ρ 0.9: 0.095–0.208; 2,500 TUs at ρ 0.9: 0.024); it
enters only through the power of the ℓ̂ test, together with ρ, the TU-size distribution and the
censored share — the same ρ = 0.9 replicon gives ℓ̂ = 9 and 0.007 uncensored, ℓ̂ = 2 and 0.208 at
55 % detected.

Two mechanisms, both measured at 1,000 TUs and ρ = 0.9 with the block forced:

1. **Autocovariance beyond the block.** Uncensored, Î₁ = 0.78: FWER 0.576 at block 5 against 0.016
   at block 19. The bootstrap breaks the dependence at every join, so a block shorter than the
   correlation length makes the null rows too tame and M_b too small.
2. **Local-scale heterogeneity at matched variance.** With 45 % censored the overall variance ratio
   is ≈ 1.0 and the FWER is still 0.158 at block 5; the block-SD ratio q90/q10 is 1.73 at ρ = 0.9
   against 1.24 at ρ = 0. A stationary null gives every window the same scale, so the windows that
   sit in a locally noisy stretch are tested against too narrow a null. This is the same mechanism
   as the ori–ter gradient attack (decision 49): left censoring compresses the low half of the Blom
   scores against the censored atom while the high half spreads, so a gradient makes the detrended
   series heteroscedastic (local-SD ratio 1.74 against 1.35 without it) and removing censoring
   removes the inflation (0.066 → 0.011).

A long block fixes both and costs all the power, which ADR-0021 decision 1 rejected. Local
standardisation (dividing by a LOESS fit of |residual|) was measured and rejected: it reduces the
inflation (ρ 0.8: 0.106 → 0.032; ρ 0.9: 0.188 → 0.084, still a breach) and absorbs the signal
(power on a 20-TU × +2 SD island at ρ 0.5: 0.103 → 0.003). A scale-aware null is on the backlog.

**No flag, and why not** (decision 48). Nothing is flagged `uncalibrated`: no per-replicon statistic
certifies the domain, and selecting on one is *adverse*. ℓ̂, ℓ′, the lag-1 Moran's I, a Bartlett
variance ratio of the window sums and a local-scale ratio all fail as keep/flag rules — at ρ = 0.9
and 557 TUs the overall rate is 0.113 but **0.220** among the replicates that look calmest (ℓ̂ ≤ 1,
flat local scale), because a small ℓ̂ is exactly what makes the block too short, and the same rule
flags 94 % of the perfectly calibrated 2,501-TU ρ = 0.5 cell. Instead the §9.4 run report carries
**ℓ̂, ℓ′ (`block_len`) and the lag-1 Moran's I of the detrended TU scores per replicon**, so a
reader can see where their data sit against the paragraphs above. Reading them: ℓ̂ is reported as 0
when it was not estimated — the replicon was not scanned, or `context.islands.block_len` fixed it —
and `block_len` is 0 only for a replicon that was not scanned. ℓ̂ = 0 is a sentinel, never a
measurement: the estimator is bounded to [1, 20], so a report must render it as "not estimated"
and never as the number 0, which would read as "no autocorrelation" — the false reassurance the
amendment refuses. The lag-1 Moran's I is computed from the data either way, so a fixed block does
not blind the reader. With `block_len` fixed the calibration above does not apply at all: it is
evidence about ℓ′, and the FWER claim is conditional on the default rule.

### The low direction is a structural zero at normal detection rates

Left censoring is not symmetric (decisions 47 and 50). Every undetected gene shares the censored
Blom score c = Φ⁻¹(q/2) with q the undetected share, so a *fully silenced* window of w TUs cannot
go below −√w·|c|: about 3.4 at q = 0.45 and w = 20, against a genome-wide critical value of about
4.1–4.9 across the settings measured. `low_direction_ceiling(q, w)` returns this bound. Measured on
a ρ = 0.5 null with 10 / 20 / 40
consecutive TUs silenced entirely: ceiling 3.42 / 3.38 / 3.35, M_b 95th percentile 4.92 / 4.57 /
4.49, best low p = 1.000, no low island at any α. The high direction has no such cap — the top Blom
score is about +3.3 at these sizes — which is why planted low islands are missed on giapsim at the
default detection rate while high ones are found.

The bound is exact for the raw TU scores; the detrended statistic differs from it by √w times the
window's mean trend, so it can sit slightly below the ceiling when the window lies on a locally high
trend (measured −2.77 detrended against −2.42 raw at 10 silenced TUs) — far below the critical value
either way. The procedure is two-sided and costs the high direction almost nothing, but **a null
result in the low direction is not evidence of no silenced island** unless the censored share is
small (roughly below a quarter undetected). M5-07 records low-island power separately and labels it
a structural zero with the ceiling and the critical value printed beside it; the run report prints
the ceiling against the 95th percentile of M_b per replicon.

### What is reported, and what the gene list is

Each cluster of merged significant same-direction windows is reported by its **best window**, never
by the union (decision 42, "Merging and the representative" above): `tu_first`/`tu_last`, `n_tus`,
unwrapped `start`/`end` (`end > length` when the window crosses the origin), `scan_s`, `block_len`,
`p_value` (FWER), `q_value` NA, `method = scan_max_stat`, `null_model = circular_block_bootstrap`.

`genes`, `n_genes` and `detected_frac` cover **all eligible genes of the window's TUs**, including a
detected gene that was dropped from the scores for want of an abundance. The asymmetry with the
scoring is deliberate: the statistic is an inference about the scored series, while the gene list is
a description of the locus, and it is the gene list that the §8.6 sensitivity Jaccard compares
across `abundance.missing ∈ {min, exclude}`. If membership followed the scored set, that comparison
would mostly measure the missing-value policy instead of the island's extent, and `detected_frac`
would stop being comparable to the replicon-level detected fraction and to the §7.3.5(c) region
tests, which also count detection over eligible genes. Two consequences to keep in mind:
`detected_frac`·`n_genes` can exceed the number of genes behind `mean_abundance_log2` (which
averages only detected genes that have a value), and a TU with no scored gene at all is absent from
the series, so the span `start … end` can contain eligible genes that `genes` does not list.

### The effect and its CI: a selected interval

`effect` is the mean detrended score of the best window's TUs — identically `scan_s`/√`n_tus`, the
test statistic on a per-TU scale — with a 95 % percentile `tu_block_bootstrap` CI over that window's
TUs (decision 26). Read it as a descriptive size, not as a test. Four caveats, all in the
anti-conservative direction:

- **Winner's curse.** The window was chosen as the genome-wide maximum of |S|, so both the point
  estimate and the interval are conditional on that selection: `effect` is biased away from zero,
  the CI under-covers the underlying island mean, and it essentially never contains zero — so it
  carries no independent evidence that the island is real. The `p_value` is the only test. On a
  600-gene null with the threshold opened to α = 0.9, the three reported windows have effects
  +2.35, +2.33, +1.67 with CIs [1.88, 2.83], [1.70, 2.96], [1.44, 1.99] — every one excluding zero
  — at p = 0.70, 0.86, 0.13.
- **Within-window dependence is ignored.** The CI resamples TUs independently inside the window;
  with positive autocorrelation that understates the variance and the interval is too narrow.
- **Scale.** The units are detrended Blom scores, not log2 fold change, and the detrend has already
  attenuated the island by the factor tabulated above.
- **Few units.** Windows of 2–9 TUs carry the M4-16 flag `ci_few_units`: a percentile bootstrap over
  so few units is coarse (at 2 TUs the interval is essentially the two values).

The realised coverage of this CI is not yet measured (backlog).

### Degenerate cases

| Situation | Behaviour |
|---|---|
| LOESS span gives k < 4 points (≈ T ≤ 30 at the 10 % span), or T < 2·`min_tus` | not scanned: run flag `few_values`, no table row, no draws, out of both maxima (decisions 40–41) |
| no replicon scanned | empty table, every `island_id` NA, `null_max` empty |
| all genes undetected | every Blom score 0 → S = 0 → p = 1 → no island |
| no valid M_b | p NaN; `p <= alpha` is false, so NaN never reports an island |
| low island with no detected gene | `mean_abundance_log2` NA, `detected_frac` 0, still reported |
| `context.islands.block_len` > T | bounded to T; ℓ̂ reported as 0 (not estimated) |
| `max_tus` ≥ T on a circular replicon | windows capped at T − 1 by the engine |

### Calibration evidence

Fast tier (`tests/calibration/test_islands_calibration.py`, R = 50, B = 199, bounds 0.0381 / 0.1116
/ 0.1849 at α = 0.01 / 0.05 / 0.10): 1,000-gene null at ρ = 0 → 0.00 / 0.00 / 0.04, at ρ = 0.5 →
0.00 / 0.02 / 0.06; chromosome + 60-gene plasmid in one family → 0.00 / 0.02 / 0.04. Teeth: the
bare ℓ̂ as the block length is a confirmed fail at ρ = 0.9 (0.62 and 0.72 in the two ADR-0015 blocks
against 0.1116). Attacks (`tests/calibration/attacks/test_attack_islands.py`, §8.5): the ori–ter
gradient at a realistic 2 log2 (0.003 / 0.042 / 0.094) and the same with detection driven by
abundance (0.000 / 0.017 / 0.057) are gated, as are large TUs, a draft assembly of 133 short
contigs and a chromosome with three plasmids; the 4-log2 gradient (0.009 / 0.062 / 0.136) and
ρ = 0.9 at three sizes are recorded, not gated, as the known breaches of the domain above. M5-07
carries the full-size gate and the power table.

Signed off by the statistician (PR #131 review).

### Calibration and power (M5-07)

**The M5 gate** (§8.4 row "Islands"; `tests/calibration/test_island_fwer.py`). Each replicate draws
giapsim `ar1_background` with no island planted; the indicator is "at least one island reported".
Settings: `stats.bootstrap = stats.permutations = 1000`, default `context.islands`, criterion `fwer`
with the ADR-0015 two-block confirmation. Levels and bounds come from
`tests/calibration/tolerances.yaml` (0.05 + 2 SE: 0.0189 / 0.0695 / 0.1268 at α = 0.01 / 0.05 /
0.10, R = 500). The layout is a circular 4,500-gene chromosome (median 2,501 scored TUs) plus a
circular 60-gene plasmid, in one genome-wide family. Branch run: Slurm 1868005 (commit fdbfc77).
Main run: Slurm 1868687 (commit ed129c4, barley-01, 1 h 10 min), rates identical to the branch run. Every cell passed on block 1. Rates carry Wilson 95 % intervals.

| Cell (R = 500) | α = 0.01 | α = 0.05 | α = 0.10 |
|---|---|---|---|
| ρ = 0 | 0.004 | 0.022 [0.012, 0.039] | 0.050 [0.034, 0.073] |
| ρ = 0.5 | 0.006 | 0.020 [0.011, 0.036] | 0.068 [0.049, 0.094] |
| ρ = 0.9 | 0.000 | 0.024 [0.014, 0.041] | 0.096 [0.073, 0.125] |
| attack (a): ρ = 0.95 | 0.002 | 0.026 [0.015, 0.044] | 0.094 [0.071, 0.123] |
| attack (b): ori–ter cosine 2 log2, ρ = 0.5 | 0.006 | 0.048 [0.032, 0.070] | 0.114 [0.089, 0.145] |
| attack (c): `tu_null`, mean TU size 3 | 0.006 | 0.022 [0.012, 0.039] | 0.030 [0.018, 0.049] |
| attack (d): ρ = 0.5, detection from abundance | 0.000 | 0.010 [0.004, 0.023] | 0.070 [0.051, 0.096] |
| attack (e): `tu_null(1800)`, short genes in the prophage | 0.000 | 0.010 [0.004, 0.023] | 0.022 [0.012, 0.039] |

Attacks (a)–(d) use 4,500 genes and no plasmid (`tests/calibration/attacks/test_attack_islands_scenarios.py`,
§8.5 attack 1, gated). Arm (b) has the least margin at both levels. It is the local-scale mechanism
of "Validity domain" item 2, whose 4-log2 version is the recorded breach. Recorded, not gated
(ADR-0021 Amendment 1): ρ = 0.9 gives 0.096 [0.073, 0.125] at 1,000 genes (557 TUs) and 0.130
[0.103, 0.162] at 1,800 genes (999 TUs), at α = 0.05 against 0.0695. Both are breaches outside the
claimed domain. Fast tier: 1,000 genes, ρ = 0.5, R = 50, 13.6 s (`benchmarks/baseline.json`).
Teeth: the iid bootstrap (`block_len = 1`) on the fast cell is a confirmed fail (0.23 at R = 100 in
the statistician's probe, against 0.1116). A fixed block of 3 on the full ρ = 0.9 cell gives 0.765
(R = 200).

**Block lengths observed** (chromosome, R = 500 per cell):

| Cell | ℓ̂ median [10 %, 90 %] | ℓ′ counts |
|---|---|---|
| ρ 0, 4,500 genes | 1 [1, 1] | 3: 459, 5: 40, 7: 1 |
| ρ 0.5, 4,500 | 2 [2, 3] | 3: 26, 5: 309, 7: 147, 9: 18 |
| ρ 0.9, 4,500 | 8 [6, 11] | 7: 2, 9: 12, 11: 27, 13: 65, 15: 87, 17: 89, 19: 85, 20: 133 |
| ρ 0.95, 4,500 (attack a) | 11 | 5: 2, 7: 2, 9: 8, 11: 13, 13: 24, 15: 31, 17: 36, 19: 64, 20: 320 |
| ρ 0.9, 1,800 (recorded) | 3 [1, 5] | 3: 96, 5: 135, 7: 127, 9: 77, 11: 40, 13: 14, 15: 9, 17: 2 |
| ρ 0.9, 1,000 (recorded) | 1 [1, 2] | 3: 352, 5: 111, 7: 28, 9: 7, 11: 2 |

ℓ′ = min(20, 2ℓ̂ + 1, T) reaches the cap whenever ℓ̂ ≥ 10: in 133 of 500 replicates at ρ 0.9 and
320 of 500 at ρ 0.95. The cap does not show in the rate. Capped replicates: 0/133 and 2/320 at
α = 0.05 (9/133 and 26/320 at 0.10); the rest: 12/367 and 11/180. The excess sits where ℓ̂ is
small for the dependence (uncapped ρ 0.95: 0.061 at α = 0.05, R = 180), the adverse selection of
decision 48, not where the cap binds. This comes from the statistician's probe for PR #188, which
reran ℓ′ per replicate and joined it to the gate indicators. Limitation: the blocks underrepresent
dependence at the joins (docs/methods/bootstrap.md, "Block length"). At 4,500 genes ℓ′ spans the
dependence well enough for 20-TU windows; the evidence does not extend to correlation lengths
beyond the cap (ρ > 0.95, or longer-range structure).

**The plasmid is not always scanned.** The 60-gene plasmid has 21–46 scored TUs (median 34). With
T ≤ 30 the LOESS span ⌈0.1·T⌉ is below 4 points, so the replicon is not scanned (`few_values`,
decision 40). That happened in 93 / 89 / 98 of 500 replicates at ρ = 0 / 0.5 / 0.9, about 19 %. The
two-replicon family was therefore exercised in 402–411 replicates per cell, with FWER 7/407, 9/411
and 11/402 at α = 0.05. An unscanned plasmid cannot carry an island and adds nothing to M_b; the
genome-wide-family claim rests on the 81 % of replicates where it was scanned.

**Power** (recorded, not gated; §7.3.6). `scripts/dev/island_power.py`, Slurm 1868246 (commit
20e5273), R = 100 per cell. Layout: the 4,500-gene `island_power` scenario at ρ = 0.5, island
starting at chromosome TU 625, `stats.bootstrap = stats.permutations = 1000`, α = 0.05. The layout
carries no regions: the scenario has no region effect and the island step reads no region, so this
changes the random stream, not the distribution of anything the step sees. *Power* is the share of
replicates with a same-direction island sharing at least one TU with the planted island; *J* is the
TU Jaccard of the best such island's best window with the planted TUs, averaged over hits (unstable
below 5 hits); *false* is the share of replicates with an island that shares no TU with the planted
one (an island counts as false when its best window misses the planted TUs, even if its cluster
touches them, which biases this rate upward). Null FWER of this layout: 0.020 [0.011, 0.036] (gate
cell ρ = 0.5).

High direction: power [Wilson 95 %] / J / false.

| TUs × SD | covariates_only | detection from abundance |
|---|---|---|
| 5 × 1 | 0.00 [0.00, 0.04] / – / 0.02 | 0.00 [0.00, 0.04] / – / 0.04 |
| 5 × 1.5 | 0.01 [0.00, 0.05] / 0.62 / 0.02 | 0.00 [0.00, 0.04] / – / 0.01 |
| 5 × 2 | 0.00 [0.00, 0.04] / – / 0.01 | 0.01 [0.00, 0.05] / 0.28 / 0.01 |
| 10 × 1 | 0.01 [0.00, 0.05] / 0.67 / 0.03 | 0.00 [0.00, 0.04] / – / 0.03 |
| 10 × 1.5 | 0.02 [0.01, 0.07] / 0.65 / 0.01 | 0.05 [0.02, 0.11] / 0.76 / 0.02 |
| 10 × 2 | 0.05 [0.02, 0.11] / 0.84 / 0.03 | 0.11 [0.06, 0.19] / 0.83 / 0.00 |
| 20 × 1 | 0.02 [0.01, 0.07] / 0.55 / 0.01 | 0.04 [0.02, 0.10] / 0.83 / 0.02 |
| 20 × 1.5 | 0.10 [0.06, 0.17] / 0.63 / 0.04 | 0.30 [0.22, 0.40] / 0.78 / 0.03 |
| 20 × 2 | 0.20 [0.13, 0.29] / 0.59 / 0.01 | 0.65 [0.55, 0.74] / 0.86 / 0.00 |

The low direction is a **structural zero**: power is 0/100 in all 18 low cells in both detection
modes (upper Wilson bound 0.04). The ceiling √w·|Φ⁻¹(q/2)| at q ≈ 0.45 is 1.68 / 2.38 / 3.37 for
w = 5 / 10 / 20, while the genome-wide critical value (the 95th percentile of M_b) has median
4.40–4.58 and minimum 3.71; across 1,200 probe replicates the ceiling was never at or above it.
False-island rates in the low cells (R = 100), covariates_only, 5 / 10 / 20 TUs × 1 / 1.5 / 2 SD:
0.04, 0.05, 0.05 / 0.03, 0.00, 0.01 / 0.06, 0.05, 0.07; detection from abundance: 0.02, 0.05, 0.05
/ 0.04, 0.03, 0.00 / 0.00, 0.03, 0.02.

The planning target of ≥ 0.8 is not met under giapsim defaults: 0.20 at 20 TUs × 2 SD. When
detection follows abundance (0.65), a high island also raises detection, and the left-censored
Blom score turns that into larger scores. Pooled false-island rates with a planted island: high
cells 18/900 = 0.020 (covariates_only) and 16/900 = 0.018 (abundance); low cells 36/900 = 0.040 and
24/900 = 0.027. At R = 500 the low 20-TU cells give 0.034 (−1 SD) and 0.036 (−2 SD): within α but
above the complete-null 0.020. All of these are high-direction islands, mostly more than 300 TUs
from the planted one (see "Scan and maximum statistic").

Signed off by the statistician (PR #188 review).
