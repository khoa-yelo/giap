# Origin-to-terminus gradient

Issue M5-08 (GIPPA_SPEC.md §7.3.5(e), §9.3), module `gippa.context.oriter`. The step plumbing that
writes the rows into `autocorrelation.tsv` and the two `genes.tsv` columns is M5-09.

## Why the regression is descriptive

Replication from a single origin makes gene dosage decrease from oriC to the terminus during
exponential growth, and highly expressed transcription and translation genes sit near oriC
*because* of it. The ori–ter gradient of protein abundance is therefore an evolved property of the
gene order, not a deviation from an exchangeable null: no permutation of the genes (or of the
values) reproduces "the same genome without a gradient", because the placement of the genes is the
effect we would be testing. GIPPA consequently reports the slope and an interval, never a p-value:
every row carries `p_value = NA`, `q_value = NA` and the flag `descriptive` (§7.3.5(e)). The
interval answers "how precisely is this genome's gradient measured", not "could it have arisen by
chance".

The rows therefore enter no §8.1 multiple-testing family. §8.1 lists no ori–ter claim; the
autocorrelation BH family is keyed on (contig, statistic, null_model) over lags with a p-value, and
these rows carry `statistic = oriter_slope`, `null_model = bootstrap` and a NA p, so they can only
form a family of their own with nothing in it. M5-09 must keep it that way: `p_value` NA ⇒
`q_value` NA, and no ori–ter row is counted in any BH denominator.

## Replichore geometry

The oriC midpoint o comes from Bakta (C15 `features`, `type = oric`); the terminus is taken as
o + L/2. For a point at coordinate m on a circular replicon of length L, with
d = (m − o) mod L:

| | condition | label | position |
|---|---|---|---|
| right replichore | d ≤ L/2 | `right` | d / (L/2) |
| left replichore | d > L/2 | `left` | (L − d) / (L/2) |

Positions lie in [0, 1]: 0 at the origin, 1 at the terminus. `right` is the arm of increasing
coordinate from the oriC (lead decision 32); the labels are directions on the circle, not
Bakta-annotated leading/lagging strands. The geometry is intrinsic to the circle: rotating every
coordinate *and* the oriC by the same offset leaves both columns unchanged, mirror images across
the oriC share a position on opposite arms, and the strand never enters. Midpoints are taken
modulo L, so an origin-crossing interval written unwrapped (`end > L`, §6.5) lands where it
belongs.

Only **circular chromosomes** are analysed (lead decision 33): the geometry above assumes a circle
with one origin, so linear replicons, `unknown`-kind contigs and plasmids get no row and NA gene
columns. Plasmid `oriV`/`oriT` calls are not origins of chromosomal replication and are ignored.

## Several oriC calls

Bakta can call more than one oriC (the recorded K-12 MG1655 annotation has two, 42 kb apart). The
calls are summarised by their **circular mean** (the angle of the mean of the unit vectors
2π·mᵢ/L), which handles calls either side of coordinate 0, but only when every pair lies within
0.05·L circular distance; farther apart the module refuses to guess and the replicon is `no_oric`
(lead decision 31). Requiring the *maximum pairwise* circular distance to be ≤ 0.05·L (not the
spread about the mean) is what makes the mean well defined: with 0.05·L < L/3 the calls then
necessarily lie inside one arc of that length, so the circular mean is inside the arc and no
antipodal cancellation is possible. Both branches log the event
`context_stats.oriter.multiple_oric` with the midpoints, the spread and the decision, so a
questionable origin is visible in `logs/events.jsonl`. For K-12 the two calls differ by 0.009 of a
replichore, far below the resolution of a descriptive slope. Without any call — or without a
`features` table at all — the replicon yields one row with `value` NA and the flags
`descriptive;no_oric`, and its genes get NA ori–ter columns.

An oriC that is wrong by δ bp shifts every x by 2δ/L and swaps arms for the genes in between. The
slope is not invariant to that (only the geometry is), so a genome whose oriC call is doubtful is
better read from the `multiple_oric` event than from the slope.

## Points and fits

The regression is at **TU level**, so that an operon counts once rather than once per gene:
y_t is the mean `abundance_log2` of the TU's detected eligible genes with a value (censoring is
ignored — the quantity regressed is the observed abundance of what was detected, lead decision
34), and x_t is the replichore position of the TU midpoint. TUs without a value are dropped.
Three OLS fits are reported per replicon: the TUs of the `right` arm, of the `left` arm and
`pooled`. The slope b is in log2 per replichore, i.e. the fitted terminus minus the fitted origin,
and is written to both `value` and `effect`. A fit with fewer than three TUs, or with no variance
in x, does not exist: its row is NA with `descriptive;not_estimable`.

Each arm is a subset of the pooled points, so a non-estimable pooled fit forces both arms to be
non-estimable: with fewer than three pooled TUs neither arm can have three, and if the pooled x has
no variance then neither subset has any. That is why ℓ̂ is only ever estimated when the pooled fit
exists, and the arms can still use it. The invariant is about subsets, not about the circle, so a
future change to the per-fit minimum would have to re-check it; `block_len_prime` raises on
ℓ̂ < 1, so a violation fails loudly rather than quietly resampling with block 1.

The pooled fit regresses on |position|, not on a signed coordinate, so it is the gradient averaged
over the two arms, weighted by their TU counts. Read the two arm rows before the pooled one when
the arms disagree; a single pooled slope hides an asymmetric origin call or a genuine asymmetry.

## The interval and its block length

Neighbouring TUs are correlated (operons, local expression domains), so an iid bootstrap of the
TUs would be too narrow — measurably so: resampling single TUs instead of blocks on the calibration
cell below gives mean width 0.78 against 1.11 and coverage 0.785 against 0.943 (R = 200/300). The
CI comes from the **circular block bootstrap** (docs/methods/bootstrap.md) of each fit's own TU
sequence — the pooled fit in genome order, an arm ordered from oriC to terminus — resampling index
blocks of (x, y) pairs and refitting; the reported interval is the 95 % percentile interval of the
`stats.bootstrap` refitted slopes. Replicates whose resampled x has no variance give a NaN slope
and are dropped; `n_perm` counts the valid ones, and fewer than B/2 of them leaves the CI NA with
`not_estimable` added (B ≥ 1 is assumed, as everywhere else in `gippa.stats.bootstrap`).

Resampling *pairs* in blocks, rather than residuals, is deliberate. A TU's y is the mean of however
many of its genes were detected, so Var y_t falls with that count and unweighted OLS is
heteroscedastic; the pairs bootstrap carries each TU's own noise level with it, so the interval
does not assume equal variances even though the point estimate ignores them. The resample also
keeps the realised detection pattern: it never re-draws which genes were detected, so the interval
is conditional on the detected set (see the caveats).

The block length follows ADR-0021 decision 1, shared with the island null: ℓ′ = min(20, 2ℓ̂ + 1, T),
with ℓ̂ the first lag at which Moran's I of the **pooled fit's residuals** is no longer significant
(`choose_block_len`, `stats.permutations` gene-label permutations). ℓ̂ is estimated once per
replicon; the clamp to T is applied **per fit** with that fit's own TU count, because an arm may
hold far fewer TUs than the pooled sequence the length was estimated from. An integer
`context.islands.block_len` overrides the rule (clamped to T as well). Those residuals still carry
whatever the straight line failed to absorb of the true (cosine-shaped) dosage curve, but that does
not drive ℓ̂ in practice: on the calibration cell ℓ̂ averages 3.91 with a 2-log2 cosine planted and
3.87 with no gradient at all, so the AR(1) background, not the leftover curvature, sets the block.

On an arm the circular wrap joins the TU nearest the terminus to the TU nearest the origin, which
are not neighbours on the chromosome; ℓ′/T of the blocks straddle that join (≈ 1 % at 20 TUs and
2,000). The pooled sequence wraps at coordinate 0 between genuine neighbours.

## What is measured

Coverage of the 95 % pooled interval for the OLS slope of the planted curve, on
`scenarios.oriter_gradient` (4,500 genes, ≈ 1,720 TUs with a value, B = 199, the real
`gippa.context.oriter.oriter_rows`). The estimand is the OLS slope of the planted cosine over the
TUs that enter the fit, which the closed form puts at −12a/π² = −2.432 for a = 2 log2 under a
uniform design; the estimator is unbiased for it (mean −2.424 against mean estimand −2.432 at
R = 300). Statistician re-measurement for this review, Slurm 1863166
(`local/scratch/stat-132-work/`); the implementer's recorded run (Slurm 1863081) is the first 50
replicates of the first row and reproduces exactly.

| cell (ρ, amplitude) | ℓ̂ | ℓ′ | R | pooled | right | left | width |
|---|---|---|---|---|---|---|---|
| ρ 0.5, 2 log2 (the gate cell) | 3.9 | 8.8 | 300 | **0.943** (SE 0.013) | 0.920 | 0.930 | 1.11 |
| ρ 0, 2 log2 | 1.1 | 3.2 | 200 | 0.925 (0.019) | 0.930 | 0.945 | 0.69 |
| ρ 0.9, 2 log2 | 17.9 | 20 | 200 | 0.900 (0.021) | 0.890 | 0.875 | 4.41 |
| ρ 0.5, **no gradient** | 3.9 | 8.7 | 300 | 0.917 (0.016) | 0.923 | 0.937 | 1.09 |
| ρ 0.5, 2 log2, abundance-driven detection | 4.0 | 9.0 | 200 | 0.940 (0.017) | 0.945 | 0.920 | 1.11 |
| ρ 0.5, 2 log2, 1,000 genes (381 TUs) | 3.1 | 7.1 | 200 | 0.905 (0.021) | 0.905 | 0.905 | 2.29 |
| ρ 0.5, 2 log2, blocks off (`block_len = 1`) | — | 1 | 200 | 0.785 (0.029) | 0.775 | 0.815 | 0.78 |

**Validity domain.** The interval is nominal to within a few points over the whole range measured:
0.90–0.94 against a nominal 0.95, never grossly wrong, and never conservative. What is claimed is
**coverage ≥ 0.90 everywhere in the table above**, and ≈ 0.94 for gippasim `ar1_background` at
ρ ≤ 0.5 on a ≈ 4,500-gene chromosome (the gate cell); the gate's numeric bound below is the
finite-R form of the first claim, not a second one. Unlike the island FWER
(ADR-0021 Amendment 1) there is no regime here in which the procedure breaks down — at ρ 0.9 the
cap ℓ′ = 20 binds and the interval becomes nearly vacuous instead (width 4.4 for an effect of
−2.4), which is a loss of precision, not of validity. No ori–ter row carries `uncalibrated`
(Amendment 1 dropped the flag).

The shortfall is systematic, not noise: it is the percentile interval's own skewness and the
blocks' under-representation of dependence across joins. ℓ′ is shared with the island scan only
because ADR-0021 made it so, and the island scan's power trade-off (which is what kept ℓ′ short)
does not apply to a slope CI; a longer block, a studentised interval or BCa would likely recover
the missing 1–3 points. That is hardening-pass work (backlog), not a merge condition for a
descriptive statistic.

## Caveats a reader must apply

The island effect CI is conditional on a *selected* window (winner's curse). Nothing is selected
here — three fixed fits per replicon, no scan, no maximum — so that caveat does not carry over.
These do:

1. **The interval is not a test.** With no gradient planted at all, the 95 % interval misses the
   true slope 0 in 8.3 % of replicates (R = 300, above). Reading "the CI excludes 0" as a 5 % test
   of "this genome has no ori–ter gradient" is wrong twice over: the error rate is not 5 %, and
   §7.3.5(e) refuses that test on purpose, because gene placement is the effect, not the null.
2. **Conditional on what was detected.** y is observed abundance among detected proteins. The
   bootstrap never re-draws detection, so the interval is for the gradient of the *observed*
   abundances of *this* detected set, not for the complete-data gradient. Under gippasim's
   abundance-driven detection the coverage of the planted curve is unchanged (0.940, above), but
   that is a statement about the simulator's selection strength, not a licence: if detectability
   covaries with position — ribosomal and transcription genes near oriC are short, abundant and
   easy — the slope is confounded with it and no interval repairs that.
3. **Conditional on the design.** The estimand is the OLS slope over the realised TU positions.
   Slopes are comparable across genomes only to the extent that the oriC call, the TU density and
   the detected share are; the dosage curve is not a straight line (the planted one is a cosine),
   so b is a linear summary whose value depends on how the TUs are spread over [0, 1].
4. **No simultaneity.** The three rows of a replicon, and the rows of different replicons, are each
   95 % intervals, not a 95 % simultaneous set; the pooled row and its two arms use the same TUs
   and are strongly dependent.
5. **The arms are the weaker rows.** They carry half the TUs, 40 % wider intervals and 1–2 points
   less coverage than the pooled row at every cell measured (0.875–0.945). Only the pooled row is
   gated.

## Statistician sign-off

Signed off by the statistician (PR #132 review). The review covered the estimand and its closed
form, the per-fit block clamp, the percentile interval, the degenerate rules and the measurements
above, and settled M5-08 lead decision 8:

**Example 7 is promoted to a gate, on the full (nightly) arm only.** The claim is the non-coverage
rate of the pooled interval on the gate cell, R = `replicates_full.islands` = 500, bounded by the
house rule `level_bound(0.07, R)` = 0.0929, i.e. **coverage ≥ 0.907**, with a non-estimable CI
counted as a miss (fail closed) and ADR-0015 two-block confirmation. At the measured 0.943 that
gate fails 0.1 % of blocks and ~1e-6 of confirmed nightlies, and catches a regression to 0.88 —
what the spec's ℓ̂ produced before ADR-0021 — with probability 0.98. The tolerated non-coverage is
0.07 and not the nominal 0.05 because the measured coverage is 0.943 ± 0.013, not 0.95: a 0.05
tolerance (bound 0.9305) would fail 12–20 % of blocks and 1.6–3.8 % of confirmed nightlies. The
fast arm (R = 50) stays **recorded**: at that R no bound has both a confirmed false-failure rate
below 1 % and any power (bound 0.90: 0.4 % false failures, power 0.32 against 0.88). The two arm
rows and every other cell of the table stay recorded.
