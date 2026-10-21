# Origin-to-terminus gradient

Issue M5-08 (GIAP_SPEC.md §7.3.5(e), §9.3), module `giap.context.oriter`. The step plumbing that
writes the rows into `autocorrelation.tsv` and the two `genes.tsv` columns is M5-09.

## Why the regression is descriptive

Replication from a single origin makes gene dosage decrease from oriC to the terminus during
exponential growth, and highly expressed transcription and translation genes sit near oriC
*because* of it. The ori–ter gradient of protein abundance is therefore an evolved property of the
gene order, not a deviation from an exchangeable null: no permutation of the genes (or of the
values) reproduces "the same genome without a gradient", because the placement of the genes is the
effect we would be testing. GIAP consequently reports the slope and an interval, never a p-value:
every row carries `p_value = NA`, `q_value = NA` and the flag `descriptive` (§7.3.5(e)). The
interval answers "how precisely is this genome's gradient measured", not "could it have arisen by
chance".

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
(lead decision 31). Both branches log the event `context_stats.oriter.multiple_oric` with the
midpoints, the spread and the decision, so a questionable origin is visible in
`logs/events.jsonl`. For K-12 the two calls differ by 0.009 of a replichore, far below the
resolution of a descriptive slope. Without any call — or without a `features` table at all — the
replicon yields one row with `value` NA and the flags `descriptive;no_oric`, and its genes get NA
ori–ter columns.

## Points and fits

The regression is at **TU level**, so that an operon counts once rather than once per gene:
y_t is the mean `abundance_log2` of the TU's detected eligible genes with a value (censoring is
ignored — the quantity regressed is the observed abundance of what was detected, lead decision
34), and x_t is the replichore position of the TU midpoint. TUs without a value are dropped.
Three OLS fits are reported per replicon: the TUs of the `right` arm, of the `left` arm and
`pooled`. The slope b is in log2 per replichore, i.e. the fitted terminus minus the fitted origin,
and is written to both `value` and `effect`. A fit with fewer than three TUs, or with no variance
in x, does not exist: its row is NA with `descriptive;not_estimable`.

## The interval and its block length

Neighbouring TUs are correlated (operons, local expression domains), so an iid bootstrap of the
TUs would be too narrow. The CI comes from the **circular block bootstrap** (docs/methods/
bootstrap.md) of each fit's own TU sequence — the pooled fit in genome order, an arm ordered from
oriC to terminus — resampling index blocks and refitting; the reported interval is the 95 %
percentile interval of the `stats.bootstrap` refitted slopes. Replicates whose resampled x has no
variance give a NaN slope and are dropped; `n_perm` counts the valid ones, and below B/2 of them
the CI is NA with `not_estimable` added.

The block length follows ADR-0021 decision 1, shared with the island null: ℓ′ = min(20, 2ℓ̂ + 1, T),
with ℓ̂ the first lag at which Moran's I of the **pooled fit's residuals** is no longer significant
(`choose_block_len`, `stats.permutations` gene-label permutations). ℓ̂ is estimated once per
replicon; the clamp to T is applied **per fit** with that fit's own TU count, because an arm may
hold far fewer TUs than the pooled sequence the length was estimated from. An integer
`context.islands.block_len` overrides the rule (clamped to T as well).

Residual dependence inside a block is what the interval is calibrated for. ADR-0021 Amendment 1
states the validity domain as a dependence regime — giapsim `ar1_background` at ρ ≤ 0.5 — rather
than a replicon size, and drops the `uncalibrated` flag: no ori–ter row carries it. The coverage
of the pooled interval against the planted gradient is recorded (not gated, the analysis being
descriptive) in `tests/calibration/test_oriter_calibration.py`.
