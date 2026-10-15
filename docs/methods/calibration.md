# Calibration harness: criteria, error rates, seeds

Implements `tests/calibration/_criteria.py` and the nightly job array (GIAP_SPEC.md §8.4, §8.5,
§11.2, §12; issue M4-05). Constants come only from the protected `tests/calibration/tolerances.yaml`
(`levels` = {0.01, 0.05, 0.10}, `se_multiplier` k = 2, `ks_min_p` = 0.001, `seeds.calibration_base`,
`replicates_fast`, `replicates_full`). Statistician sign-off: see the end of this note.

## What R means

R is the number of **independent** null (or truth) replicates: replicate i of a stream is a fresh
simulated data set with its own generator, and the method is run once on it. Each criterion is a
statement about R i.i.d. Bernoulli or bounded variables, so the SE formulas below are valid only if
replicates do not share draws (see Seeds). Units *within* a replicate (genes, regions) are not
replicates: a suite that pools p-values from many regions of one genome must either keep one p per
replicate or use a criterion that accounts for within-replicate dependence.

## The four criteria (all one-sided: conservative is allowed, anti-conservative is not)

| Criterion | Statistic | Passes when | SE |
|---|---|---|---|
| Level | r̂ = #{p ≤ α}/R | r̂ ≤ α + k·√(α(1−α)/R) | nominal binomial SE at the null rate α |
| Uniformity | KS p-value of the R p-values vs U(0, 1) | KS p > 0.001 | none (exact KS) |
| FDR | mean FDP over R | mean ≤ q + k·sd(FDP; ddof 1)/√R | sample SE of the mean |
| FWER | share of replicates with ≥ 1 false discovery | share ≤ α + k·√(α(1−α)/R), α = 0.05 | nominal binomial SE |

- **Level** uses p ≤ α (inclusive), which is the definition of a valid p-value (P(p ≤ α) ≤ α). The
  SE is evaluated at the nominal α, not at r̂, so a method that never rejects is not rewarded with a
  tighter bound.
- **Uniformity** is for continuous p only. Discrete p (permutation p = (1 + #{T* ≥ T})/(B + 1),
  Besag–Clifford p) have atoms, and KS against U(0, 1) rejects them for large R even when they are
  exactly valid. The harness refuses `check_uniform` when more than 5 % of values are repeated
  (1 − #unique/R > 0.05) and such p-values use the level check only (§8.4). Permutation p with
  α(B + 1) an integer (B = 199: 2, 10, 20 for the three levels) hits α exactly, so their level checks
  sit **on** the boundary and have the full false-alarm rate below.
- **FDR** uses the sample SD because FDP has no known variance. When every FDP is 0, SE = 0 and the
  bound is q itself (the run passes, correctly). FDP is skewed (an atom at 0 at small q), so the
  normal approximation behind "q + 2 SE" is rough at small R·q; the Bernoulli-FDP worst case is
  computed exactly below.
- **FWER** reuses the level bound with α = 0.05 (the spec's "0.05 + 2 SE"; the nominal binomial SE).

**Inputs (fail closed).** Empty inputs raise. `check_level`, `check_uniform`, `check_fwer` and
`check_fdr` raise ValueError on any non-finite value or value outside [0, 1]. `check_fwer` accepts
only 0/1 outcomes, and `check_fdr` requires R ≥ 2, since SE(FDP) is undefined for R = 1.

## False-alarm rate of each check (procedure exactly at level)

Exact binomial P(fail) for a level/FWER check when the true rejection rate equals α, and power when
it is 1.5 α (k = 2, bound floor-ed to an integer count):

| R | α = 0.01: false alarm / power | α = 0.05 | α = 0.10 |
|---|---|---|---|
| 20 | 0.017 / 0.04 | 0.076 / 0.19 | 0.043 / 0.17 |
| 50 | 0.089 / 0.17 | 0.038 / 0.17 | 0.025 / 0.21 |
| 100 | 0.079 / 0.19 | 0.028 / 0.22 | 0.021 / 0.33 |
| 300 | 0.033 / 0.17 | 0.029 / 0.49 | 0.025 / 0.76 |
| 500 | 0.031 / 0.22 | 0.030 / 0.69 | 0.025 / 0.93 |
| 1,000 | 0.026 / 0.34 | 0.028 / 0.92 | 0.028 / 0.998 |

The nominal one-sided normal rate of "mean > α + 2 SE" is 0.023; discreteness pushes it to 0.026–0.09.
The worst cells are α = 0.01 at R = 50/100 (fail iff ≥ 2 resp. ≥ 3 rejections; 8–9 % false alarms,
< 20 % power): the fast α = 0.01 check is mostly noise and should be read as a smoke test only.

FDR check, exact for Bernoulli FDP (global null, FDP ∈ {0, 1}, the most skewed case):
false alarm 0.0001 / 0.003 / 0.009 at R = 50 and q = 0.01 / 0.05 / 0.10, and 0.004 / 0.009 / 0.011
at R = 300. For realistic FDP (BH, m = 1,000, π₀ = 0.9, Beta(0.1, 1) alternatives, oracle level
q/π₀ so FDR = q exactly) simulation gives 0.01–0.07 per q at R = 50–300.

Validity of the seeding (below) was checked empirically: 400 independent streams of R = 1,000
exact U(0, 1) p-values failed the level check at rates 0.018 / 0.013 / 0.025 (theory 0.026 /
0.028 / 0.028), and the pooled 400,000 draws pass KS (p = 0.83).

### Family-wise false-alarm rate per nightly run

A nightly run has roughly 46 checks once all suites exist (6 level suites × 3 α, regions and
autocorrelation with several test types, 3 engine level streams, 5 FDR streams × 3 q, 1 FWER). If
every method sat exactly at level and checks were independent, P(≥ 1 red) ≈ 1 − ∏(1 − fᵢ) ≈ **0.75**.
Checks of one stream at three α are positively correlated and most methods are conservative (BH has
FDR = π₀q), so the real rate is lower, but permutation and shift tests are exact-level and alone
contribute ≈ 20–40 %.

A single-block gate with these numbers would not be acceptable. Seeds are fixed, so an unchanged
tree gives the same verdict every night. The false alarm appears when a PR changes how a method
consumes randomness and a correct method lands in the 2–3 % tail. The attack stream
`regions.attack_harness_valid` did exactly that: 122/1,000 U(0, 1) p ≤ 0.10, P = 0.013.
Seed-shopping to make such a check green is forbidden.

### Two-block confirmation (ADR-0015, accepted)

`calibrate()` runs a check on replicates 0 … R − 1 (block 1). If the check fails there, it is re-run
on the disjoint, deterministic block R … 2R − 1 of the same stream (block 2). The gate fails only if
both blocks fail; a block-1 failure that block 2 does not repeat is logged `unconfirmed` (shown in
the summary, not red). Records carry `block` and `status` (pass | unconfirmed | fail). The two
blocks are independent, so:

- per-check confirmed false alarm = f² ≤ 0.0081 (worst cell α = 0.01, R = 50), ≈ 0.0008 at R = 1,000;
- **confirmed family-wise false-alarm rate per nightly run ≈ 0.04** (46 checks all exactly at level,
  independence; lower in practice), against ≈ 0.75 for a single block;
- power of a confirmed fail = (single-block power)²: against 1.5 α at R = 1,000 it is 0.85 (α = 0.05)
  and 0.996 (α = 0.10); against 1.5 q FDR with realistic BH FDP it is 0.36 (q = 0.01) / 1.0 at fast
  R = 50 and 1.0 at full R = 300;
- `unconfirmed` occurs at rate ≈ f (≈ 3 % of checks, about one or two per nightly run) and is
  informative: a method with persistent `unconfirmed` checks on the same stream deserves a look.

Empirical check (attack `test_two_block_confirmed_false_alarm_matches_theory`): 200 streams × 3 α
of exact U(0, 1) p at R = 1,000 gave 0 confirmed `fail` in 600 checks (expected 0.5) and 24
`unconfirmed` (5 / 8 / 11 at α = 0.01 / 0.05 / 0.10; expected ≈ 5–6 each, within the 99.9 %
binomial band). Oracle BH at FDR = 1.5 q is a confirmed `fail` at every q, fast and full.

**Sharded runs.** `--shard K/N` sets `GIAP_CALIBRATION_REPLICATES`; `calibrate()` then runs only
that range of block 1 and logs a `shard` record with the per-replicate values. `--merge` puts the
full R back together (refusing gaps or conflicting duplicates), evaluates block 1, and runs block 2
itself for every failing check. It never reports `fail` from block 1 alone. The merged verdict
equals the unsharded one (tested).

**Scope.** ADR-0015 and sharding apply only to checks that go through `calibrate()`. Suites that
still call `run_replicates` / `check_*` directly (today all suites, including the M4-05 engine
tests) keep single-block semantics until they are migrated (M4-05 lead decision 7, backlog). Until
then the nightly family-wise false-alarm rate of those checks is the single-block figure above.

## Power to detect a 1.5× FDR inflation

For the FDR check, power ≈ Φ(0.5·q·√R / sd(FDP) − 2). With FDP sd from simulation (BH, m = 1,000,
π₀ = 0.9: sd 0.014 / 0.029 / 0.038 at q = 0.01 / 0.05 / 0.10):

| R | q = 0.01 | q = 0.05 | q = 0.10 |
|---|---|---|---|
| fast, 50 | 0.6 | 1.0 | 1.0 |
| full, 300 | 1.0 | 1.0 | 1.0 |

For the Bernoulli-FDP worst case (global null) power is far lower: 0.001 / 0.03 / 0.12 at R = 50 and
0.04 / 0.32 / 0.65 at R = 300. A suite whose null replicates make few discoveries (FDP mostly 0 or 1)
cannot see a 1.5× inflation at q = 0.01 even at full R; such suites should plant enough true
discoveries that each replicate calls tens of features, or report the per-replicate FDP histogram.

The attack `tests/calibration/attacks/test_attack_calibration_harness.py` demonstrates both sides:
oracle BH at FDR = q exactly passes, and at FDR = 1.5 q is flagged at every q, fast and full; the
level check flags p = U/1.5 at α = 0.10 at full R; and exact-level p fail across 200 streams at the
binomial false-alarm rate, not more.

## Seeds and independence

`suite_seed(stream, rep) = SeedSequence([calibration_base, crc32(stream), rep]).generate_state(1,
uint64)`, and replicate `rep` draws from `default_rng(suite_seed(stream, rep))`.

- **Independence.** SeedSequence hashes its whole entropy list, so distinct (stream, rep) triples
  give statistically independent streams (NumPy's guarantee for distinct entropy). Folding to one
  uint64 and re-seeding loses nothing that matters: the chance that any two of N replicate seeds
  collide is ≈ N²/2⁶⁵ (< 10⁻¹² for N = 10⁴).
- **crc32.** A 32-bit hash of the stream name; two names collide with probability 2⁻³², and a
  collision would make two streams draw *identical* data (not merely correlated). With n names the
  chance of any collision is ≈ n²/2³³ (< 10⁻⁶ for n = 100). The 23 names in use today (13 suites
  plus `gene_fdr.*`, `regions.*` streams) have distinct crc32. Test streams must have distinct names;
  two tests sharing a stream name share draws.
- **Sharding.** The seed depends only on (stream, rep), never on the shard, so replicate i is the
  same data whether run alone or in a block, and `shard_replicates(R, n, k)` partitions 0 … R − 1
  into contiguous blocks (tested). Criteria are only valid on the merged R replicates (`--merge`, above).
- **Fixed seeds** (§8.4): `calibration_base` lives in the protected file; tests are deterministic.

## Suites and `tolerances.yaml`

Array task i runs the i-th key of `replicates_full` (13 suites → `--array=0-12`); `--fast` uses
`replicates_fast`. Engine tests reuse `gene_fdr` (BH) and `regions` (permutation, circular shift,
Besag–Clifford), per M4-index lead decision 6. Each calibration test belongs to exactly one suite.
The suite comes from a `@pytest.mark.calsuite("<suite>")` marker or from `SUITE_REGISTRY` in
`tests/calibration/conftest.py` (keyed by function, then file). A test then gets the keyword
`calsuite_<suite>`, and a unit test checks that every calibration test lands in exactly one array
task.

## Review findings (M4-05, statistician)

First review (PR #54, BLOCK), all resolved at head 6fba049:

1. `check_level` failed open on NaN and p > 1. **Fixed:** all checks fail closed. The attack tests
   `test_level_check_fails_closed_on_nan` and `..._out_of_range_p` pass.
2. `--shard` was inert. **Fixed:** shards are honoured; `--merge` evaluates the union and runs the
   confirmation block.
3. `check_fdr` accepted R = 1 and FDP outside [0, 1]. **Fixed.**
4. Suite tagging by source regex. **Fixed:** explicit marker or registry.
5. Family-wise false alarms ≈ 0.75. **Resolved by ADR-0015** (≈ 0.04 confirmed).

Open follow-up (not blocking this PR): migrate the existing suites to `calibrate()` (lead decision
7). Until then ADR-0015 protects only new suites.

## Sign-off

Statistician, 2026-10-14 (M4-05 re-review at 6fba049): **signed.** The criteria, SE formulas,
one-sidedness, k, discrete-p rule, fail-closed inputs, seed derivation, shard partition and merge,
and ADR-0015 confirmation match §8.4, ADR-0015 and the M4-05 spec with lead decisions 1–7. The
false-alarm and power numbers above are exact binomial or simulated as stated. The confirmed nightly
family-wise false-alarm rate is ≈ 0.04 once suites run through `calibrate()`.
