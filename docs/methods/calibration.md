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

**Inputs.** Empty inputs raise. `check_level` and `check_fwer` must fail closed on non-finite or
out-of-range inputs (see Review findings); `check_fdr` with R = 1 has SE 0 and is not a test, so
suites must use R ≥ 2.

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

**This is not acceptable as a hard gate without a confirmation rule.** Seeds are fixed, so an
unchanged tree gives the same verdict every night; the false alarm appears when a PR changes how a
method consumes randomness and a correct method lands in the 2–3 % tail (the attack
`regions.attack_harness_valid` stream did exactly that: 122/1,000 U(0, 1) p ≤ 0.10, P = 0.013).
Seed-shopping to make such a check green is forbidden. Recommended rule (needs a lead decision or
ADR, because it changes the gate semantics of §8.4): a red check is re-run on the next disjoint,
deterministic replicate block (replicates R … 2R − 1 of the same stream) and the failure stands only
if that block fails too. Per-check false alarm becomes fᵢ² (≤ 0.008), the nightly family rate
≈ 0.04, and power against 1.5 α at R = 1,000 drops only from 0.92 to 0.85 (α = 0.05) and
0.998 to 0.996 (α = 0.10).

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
  into contiguous blocks (tested). Criteria are only valid on the merged R replicates; see Review
  findings for the state of `--shard`.
- **Fixed seeds** (§8.4): `calibration_base` lives in the protected file; tests are deterministic.

## Suites and `tolerances.yaml`

Array task i runs the i-th key of `replicates_full` (13 suites → `--array=0-12`); `--fast` uses
`replicates_fast`. Engine tests reuse `gene_fdr` (BH) and `regions` (permutation, circular shift,
Besag–Clifford), per M4-index lead decision 6. A test joins a suite through the keyword
`calsuite_<suite>`, assigned when the suite name appears quoted in the test function's source, else
anywhere in its module (so a module naming several suites runs in several tasks).

## Review findings (M4-05, statistician)

1. **Blocking.** `check_level` fails open: NaN p and p > 1 count as "not rejected", so 100 NaN
   p-values pass at α = 0.05 (observed 0, ok = True). A method that returns NaN for every region
   would pass its level gate. `check_level` and `check_fwer` must raise on non-finite or
   out-of-[0, 1] input. Attack tests `test_level_check_fails_closed_on_nan` and
   `..._out_of_range_p` fail until fixed.
2. **Blocking.** `--shard K/N` exports `GIAP_CALIBRATION_REPLICATES` but no test reads it: every
   shard runs all R replicates and logs a full-R verdict (N× compute, duplicated records). Either
   honour the range in `run_replicates`/tests **and** add a merge step that evaluates the criteria on
   the union, or remove `--shard` until then.
3. Minor: `check_fdr` accepts R = 1 (SE 0, trivially passes) and FDP outside [0, 1]; require R ≥ 2
   and 0 ≤ FDP ≤ 1.
4. Minor: suite tagging by source regex is fragile (a helper naming two suites puts a test in two
   array tasks); an explicit `@pytest.mark.calsuite("regions")` marker would be exact.
5. Policy (lead/ADR): adopt the two-block confirmation rule above, or accept a ≈ 50–75 % nightly
   family false-alarm rate when all methods are at level.

## Sign-off

Statistician, 2026-10-14 (M4-05 review). Criteria, SE formulas, one-sidedness, k, discrete-p rule
and seed derivation match §8.4 and the M4-05 spec. The false-alarm and power numbers above are
exact binomial or simulated as stated. Sign-off on the harness is conditional on findings 1–2;
finding 5 is a recommendation to the lead.
