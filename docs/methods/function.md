# Function module: detection enrichment primitives

Implements `giap.stats.wallenius` and `giap.function.elim` (GIAP_SPEC.md §7.7.3, §7.7.5, §8.1 row
"Category detection", §8.2, §8.3; issue M8-01, docs/specs/M8-01-stats-wallenius-elim.md, lead
decisions 1–4; M8-index lead decision 5). The families, the effect size and `parent_group` are
M8-04's ("Detection enrichment", still to be written). The `ora` suite and attack 2 are M8-15's.

**Status:** statistician review of #153 (2026-10-21). The primitives are **approved as
implemented**. The **method is not yet approved as a calibrated test**: the mean-odds Wallenius
null fails the §8.4 ORA null and attack 2 (see "Calibration evidence"). See Sign-off.

## Wallenius ORA

### Quantities

Per category C over the eligible universe U that the caller passes:

- M = |U|, n = |C ∩ U|, N = the number detected, k = the number detected in C;
- the odds ω = mean(p_detect in C) / mean(p_detect outside C). This is goseq's weighting with
  p_detect as the PWF. Algebraically it equals `avg_w·(M − n)/(Σp − n·avg_w)`, and a unit test pins
  the identity.

ω is NaN when C or its complement is empty, or when the outside mean is 0. Those categories are
not tested; M8-04 writes them `unestimated`. ω = 0 is allowed and is exact: X = lo almost surely.
ω = 1 uses the central hypergeometric pmf.

### Tails: direct sums of the pmf

p_over = Σ_{j=k}^{hi} pmf(j) and p_under = Σ_{j=lo}^{k} pmf(j), with lo = max(0, N − (M − n)) and
hi = min(n, N). Both are `math.fsum` sums (correctly rounded) of scipy's
`nchypergeom_wallenius.pmf`, clipped to [0, 1].

These are the same quantities as the spec's `.sf(k − 1)` / `.cdf(k)`. The sums are used because
scipy computes `sf` as 1 − cdf, which cancels to 0.0 below about 1e-16 (M8-01 example 1: true
3.2e-14, `sf` returns 0). p_over and p_under are each a valid one-sided p under the urn model:
they include the observed k, so p_over + p_under − pmf(k) = 1.

### Fail-closed range (lead decision 4)

Outside ω ∈ [1e-6, 1e6] scipy's pmf turns NaN, and near 1e±300 it aborts the process. A pmf that
is non-finite or does not sum to 1 ± 1e-6 also raises `ValueError`. M8-04 reports such a category
`unestimated`, flag `odds_out_of_range`.

**Is the range statistically acceptable for pass 1? Yes.**
- **Upper bound.** ω ≤ 1/mean(p outside C). With a size cap of 500, the outside mean is close to
  the run's detection rate (about 0.3–0.7), so ω > 1e6 cannot occur. With fewer than 50 detected
  genes, `detectability_constant` sets ω = 1.
- **Lower bound.** ω < 1e-6 needs every member of C to have p_detect below 1e-6 × the outside mean.
- **Simulation.** I used 300 giapsim genomes (M ≈ 4,450, M4-13 cross-fitted p̂). The smallest p̂
  was 2.7e-4. The worst possible size-5 category (the five smallest p̂) had ω ≥ 1.1e-3, with a
  median of 1.1e-2.
- **Accuracy inside the range.** Against an exact forward urn DP, the tails at ω ∈ {1e-6, 1e-4,
  1e-2, 1e2, 1e4, 1e6} are within 5e-11 relative error, for M up to 4,500 and tails ≥ 1e-8.
- **No spurious failures.** Over 400 random (M ≤ 12,000, n ≤ 500, ω ∈ [1e-3, 1e3]) cases, the
  pmf-sum check never raised.

So the fail-closed path is a safety net, not an expected outcome. It is R4-honest because the
category is labelled `unestimated`, never given a p.
- **Elim caveat.** elim recomputes ω on a reduced set, which can be a single gene; a single gene
  with p̂ < 1e-6 × the outside mean would also fail closed. This is not seen in the simulations.

### Far-tail precision of scipy's pmf (lead decision 2)

Compared against the exact urn DP over 60 random (M ≤ 1,500) urns and every k (15,585 tails),
the maximum relative error per tail magnitude was:

| true tail | max relative error |
|---|---|
| ≥ 1e-4 | 5e-10 |
| 1e-8 – 1e-4 | 3e-7 |
| 1e-11 – 1e-8 | 2e-4 |
| 1e-13 – 1e-11 | 0.04 |
| < ~1e-14 | up to 1 (scipy's pmf returns exactly 0.0, so p = 0) |

**Decisions at BH q ≤ 0.1 and at the elim cut are unaffected.**
- **BH.** BH q is a continuous, coordinatewise-monotone function of the p vector. A relative
  error ε in each p moves each q by at most ε relative. A q near 0.1 in a family of m ≤ ~3,000
  terms needs p ≥ 0.1/m ≈ 3e-5, where ε ≤ 5e-10.
- **elim.** The elim cut at 0.01 sees ε ≤ 1e-11.
- **Deep tails.** Tails below 1e-11 have q below about 1e-7, far below any threshold. A deep
  tail returned as 0 is still below every threshold, and ties at 0 do not change which terms
  pass.

What it does affect:
- **Reported precision.** `p_value = 0` can appear in `enrichment_detection.tsv` for true values
  around 1e-14 and below. That precision is not real.
- **M8-13 golden.** The relative-1e-6 comparison with goseq cannot hold for p below about 1e-9
  unless both sides are exact there.

Both are pass-2 items: a log-space or exact pmf, an ADR if needed, and never a tolerance change.

### The approximation and the attack it implies

**Mean-odds Wallenius is exact only for its own model.** That model is sequential weighted
sampling of N genes without replacement, with weight ∝ p_detect and a single weight inside C and
another outside. The fast calibration test (`tests/calibration/test_wallenius_calibration.py`)
draws exactly that urn, so it checks the code, not the method.

**The §8.4 ORA null is a different model.** There, detection depends only on p_detect: each gene
is detected independently with probability p_i. Given N, that law is the multivariate Fisher
non-central hypergeometric distribution with odds p_i/(1 − p_i), not a Wallenius urn with weights
p_i. Mean-odds Wallenius departs from it in two ways:
1. **Saturation.** With N/M ≈ 0.5, as in proteomics, the urn's mean of k is pulled toward the
   central value. Example: p = 0.9 in C and 0.5 outside, n = 60, M = 4,000, ω = 1.8. Wallenius
   gives E[k | N] ≈ 43, but the truth is 54. So |log ω| is effectively too small. A low-p_detect
   category then looks under-detected, and a high-p_detect one over-detected. The two
   approximations agree only when p_detect ≪ 1, which is goseq's DE setting.
2. **Heterogeneity within C (Jensen).** A category that mixes p = 0.02 and p = 0.9 genes is not
   one urn colour. A hazard-weight variant (ω from mean −log(1 − p)) fixes 1 but not 2: it stays
   at level for homogeneous categories, yet rejects 52 % at α = 0.01 for the mixed one.

**The exact test under the §8.4 null** is the conditional Poisson-binomial ("conditional
Bernoulli"):

  P(k | N) = PB_C(k) · PB_{U∖C}(N − k) / PB_U(N),

where PB_A is the Poisson-binomial pmf of the p_i in A. The DP costs O(n² + (M − n)·N) per
category: tens of milliseconds vectorised at M = 4,500, and the `poisson_binomial_sf` primitive
of §8.3 already exists. Multiplying every p_i/(1 − p_i) by a constant does not change it, so a
global intercept miscalibration of p_detect cancels. This is the **Poisson-binomial fallback**
the issue asks the note to name. Swapping it in is an ADR for Khoa (see Sign-off).

**TU dependence (new ORA attack, §8.5 "new attacks").**
- **What it is.** Detection really depends on abundance, which TUs share. p_detect uses
  physicochemical covariates only. So detection is positively correlated within a TU beyond
  p_detect, and a category packed into few operons has an over-dispersed k.
- **Who fails.** Every gene-independent null fails this attack, the exact conditional
  Poisson-binomial included. It is the ORA analogue of attack 3, and §8.2 "resample TUs, not
  genes" applies.
- **Candidate fixes.** A TU-level count (TUs as units, as in the GSEA collapse), or a design-effect
  variance from the within-TU correlation of d − p_detect. M8-15 must calibrate either.

### Calibration evidence (statistician simulations, 2026-10-21)

**Where to reproduce.** Scripts are in `$GIAP_ROOT/stat/m8-01/` (`part_a.py`, `part_b.py`,
`part_c.py`, `elim_sim.py`, `acc.py`). Slurm jobs were 1865116 (A and B), 1865118 (elim) and
1865148 (C).

**Level bounds.** The §8.4 bound is α + 2√(α(1 − α)/R): 0.0145 / 0.0597 at R = 2,000, 0.0215 /
0.0752 at R = 300, and 0.0211 / 0.0744 at R = 320 (α = 0.01 / 0.05).

**Methods compared.** Every method was evaluated on the same d and p:
- W = mean-odds Wallenius (`wallenius_ora`);
- CPB = the conditional Poisson-binomial;
- H = unweighted hypergeometric.

**A. Independent Bernoulli(p_i) detection, true p, M = 4,000, R = 2,000.** Each cell is the
rejection rate at α = 0.01, with α = 0.05 in brackets. Rates above the bound are failures.

| scenario | n | W over | W under | CPB over | CPB under | H (the failing side) |
|---|---|---|---|---|---|---|
| p = 0.5 everywhere | 300 | 0.011 (0.049) | 0.007 (0.047) | same as W | same as W | same as W |
| C 0.9, outside 0.5 | 60 | **0.871 (0.976)** | 0 (0) | 0.002 (0.025) | 0.005 (0.038) | over 1.000 |
| C 0.9, outside 0.5 | 300 | **1.000 (1.000)** | 0 (0) | 0.006 (0.038) | 0.009 (0.050) | over 1.000 |
| C 0.2, outside 0.6 | 60 | 0.001 (0.003) | **0.074 (0.231)** | 0.006 (0.049) | 0.007 (0.038) | under 1.000 |
| C 0.2, outside 0.6 | 300 | 0 (0) | **0.555 (0.791)** | 0.009 (0.047) | 0.006 (0.045) | under 1.000 |
| C 0.02, outside 0.05 (rare detection) | 300 | 0.004 (0.030) | 0.004 (0.025) | 0.004 (0.033) | 0.003 (0.021) | under 0.564 |
| Beta(2, 2), C random | 300 | 0.003 (0.032) | 0.005 (0.033) | 0.005 (0.049) | 0.009 (0.051) | at level |
| attack 2: C ~ Beta(1, 4), outside Beta(2, 2) | 60 | 0 (0.003) | **0.032 (0.132)** | 0.006 (0.036) | 0.009 (0.037) | under 0.994 |
| attack 2: same | 300 | 0 (0) | **0.218 (0.490)** | 0.011 (0.043) | 0.005 (0.038) | under 1.000 |
| attack 2: bimodal background, C 80 % low | 300 | 0 (0) | **0.219 (0.558)** | 0.008 (0.041) | 0.006 (0.042) | under 1.000 |
| attack 2: C half 0.02, half 0.9 | 300 | 0 (0) | 0 (0.006) | 0.011 (0.044) | 0.006 (0.040) | under 0.025 (0.237) |

n = 20 behaves like n = 60, with smaller excesses. Weighting helps compared with H, but W stays
anti-conservative wherever ω ≠ 1 and detection is not rare.

**B. giapsim layout with `short_genes_in = "prophage"` and `covariates_only` detection.** This
uses the M4-13 cross-fitted p̂ (W/CPB "hat") and the true p ("true"), M ≈ 4,450, N/M ≈ 0.55,
R = 300. The table shows the under tail; the over tails are at level unless stated.

| category | n | ω̂ | W hat | W true | CPB hat | CPB true | H |
|---|---|---|---|---|---|---|---|
| prophage (attack 2, §8.5) | ~40 | 0.28 | **0.027 (0.137)** | **0.030 (0.123)** | 0.003 (0.030) | 0.003 (0.030) | 1.000 (1.000) |
| 60 lowest-p genes | 60 | 0.05 | **0.027 (0.157)** | 0 (0) | 0 (0.040) | 0 (0) | 1.000 (1.000) |
| tm_segments ≥ 2, over tail | ~1,090 | 1.21 | **0.463 (0.913)** | **0.483 (0.760)** | 0 (0) | 0.013 (0.040) | 1.000 |
| random genes, 10 / 50 / 200 | | ≈ 1 | at level | at level | at level | at level | at level |

The tm ≥ 2 category lies above the size cap of 500, but categories near the cap with ω ≈ 1.2
behave the same way. **The §7.7.5 attack ("stays at level with weights and fails without them")
therefore fails for mean-odds Wallenius and passes for CPB.**

**C. Detection driven by TU-shared abundance.** giapsim `detection = "abundance"` with the
log2-abundance logit coefficient b₄ ∈ {0.7, 1.5}; p̂ from physicochemical covariates; N/M ≈ 0.53;
R = 320. Categories are random whole TUs, so the category-level null holds.

| category | b₄ | W over | W under | CPB over | CPB under |
|---|---|---|---|---|---|
| ~200 genes in whole TUs | 0.7 | **0.025** (0.066) | 0.019 **(0.094)** | **0.031 (0.087)** | **0.025 (0.100)** |
| ~200 genes in whole TUs | 1.5 | **0.037 (0.116)** | **0.050 (0.091)** | **0.037 (0.119)** | **0.050 (0.091)** |
| ~50 genes in whole TUs | 1.5 | **0.025 (0.100)** | **0.028** (0.075) | **0.025 (0.106)** | **0.028** (0.075) |
| 200 random genes | 1.5 | 0.009 (0.037) | 0.016 (0.050) | same | same |

Both gene-level nulls are inflated by TU packing, and the inflation grows with how strongly
abundance drives detection. Real detection depends on abundance more strongly than b₄ = 1.5, and
bacterial GO and COG categories are operon-packed.

### Effect size

`OraResult.odds` is the weighting ω, not an effect. §8.2 requires an effect with a 95 % CI on
every p or q. M8-04 reports O − E = k − Σ_C p_detect with the CI NA in pass 1; that CI is owed.
A test-consistent effect would be the conditional MLE of a log-odds shift δ for members of C (the
odds p_i/(1 − p_i) multiplied by e^δ), with an exact CI by inverting the CPB test. The TU-block
bootstrap CI of O − E (M4-16 construction) is the TU-aware alternative.

## GO elim

### Algorithm

Terms are visited deepest first: level = the longest path to a root, ties broken by term ID. Each
term is tested on its genes minus those removed for it. A term whose elim p is strictly below
0.01 (topGO's `< cutOff`) removes the genes it was tested on from every ancestor (the transitive
closure). M and N stay fixed. M8-04 recomputes ω on the reduced set.

`p_classic` is the test on the full set. An empty remainder gives NaN, and M8-04 writes it
`unestimated`. Cycles raise before any test runs. M8-04 runs elim for the `over` direction of
go_bp, go_mf and go_cc only; every other family uses `p_value = p_classic`.

### Validity of BH on elim p (§8.1 family = (collection, direction))

The elim p is conditional on the removal steps, which depend on the same data. Alexa et al.
(2006) and the topGO documentation present elim p as conditional scores and do not derive
multiplicity control for them. **BH on elim p is not covered by the BH theorem**: its
dependence is neither independent nor PRDS by any argument I know. In the direction it is run,
the selection bias is conservative:
- removal takes out genes of a significantly **over**-detected descendant;
- this leaves the ancestor's remainder depleted, while N still counts the removed genes;
- conditioning on non-significant descendants also biases their genes downwards.

Applying the `over` removals to an `under` test would be anti-conservative. M8-04 avoids this:
`under` is never elim'd.

**Simulation.** I built a synthetic GO-like DAG: 2,000 genes; 200 leaves of 5–30 genes; 60
mid-level terms, each the union of 2–6 leaves plus 5 own genes; 15 top terms; the true-path rule;
multiple parents; size filter [5, 500]. p = 0.5 everywhere, so the central hypergeometric is the
exact conditional null. "Planted" leaves get p_plant. BH ran at 0.05 within the one family, with
R = 1,000 per row.

| setting | FWER, complete-null terms (classic / elim) | mean FDP, terms with no planted gene (classic / elim) | FDP counting inherited ancestors as false (classic / elim) |
|---|---|---|---|
| complete null | 0.016 / 0.012 | 0.016 / 0.012 | same |
| 5 planted, p 0.95 | 0.006 / 0.005 | 0.0015 / 0.0015 | 0.155 / 0.013 |
| 5 planted, p 0.75 | 0.004 / 0.002 | 0.0023 / 0.0013 | 0.056 / 0.018 |
| 20 planted, p 0.95 | 0.004 / 0.002 | 0.0002 / 0.0001 | 0.309 / 0.020 |
| 20 planted, p 0.75 | 0 / 0 | 0 / 0 | 0.173 / 0.023 |
| 20 planted, p 0.65 | 0.001 / 0.001 | 0.001 / 0.001 | 0.065 / 0.024 |

Under a correct per-term null, BH on elim p is conservative for both readings of a false
discovery:
- the strict one: a term with no signal genes at all;
- the elim one: an ancestor whose only signal is inherited from planted descendants.

elim keeps 81–98 % of the classic true positives (lowest with 20 weak plants, p 0.75).

**Residual risk.** A truly enriched child that misses p < 0.01 leaves its signal in its
ancestors, which can then be called. The rate here is ≤ 0.024 FDP at q = 0.05.

**Verdict.** The elim p are acceptable for the (collection, direction) BH families under R4 as
long as the report and this note say:
- the q-values are BH on elim p;
- the claim is "over-detected after removing the genes of its elim-significant descendants";
- FDR control for that claim is empirical (simulation), not a theorem.

M8-04's labels (`p_classic` kept, `p_value` = elim) meet this. The per-term null must itself be
valid, and with mean-odds Wallenius it is not (above). elim does not repair that.

## Sign-off

**Statistician (#153 review, 2026-10-21).**

**Approved: the M8-01 primitives** (`category_odds`, `wallenius_tails`, `wallenius_ora`, `elim`).
They compute what the spec defines:
- the direct tail sums are exact to 5e-10 relative for every decision-relevant p;
- the fail-closed range is never reached by realistic p̂;
- elim matches topGO's elim with a strict cut, and BH on elim p is empirically conservative.

**Not approved: mean-odds Wallenius as the calibrated ORA claim of §8.1.**
- It fails the §8.4 ORA null (independent detection with p_detect) and the §7.7.5/§8.5 attack 2,
  both in synthetic urns and in giapsim with cross-fitted p̂.
- Every gene-level null, the exact one included, fails the new TU-packing attack.

**Forwarded:**
- M8-15:
  - the Bernoulli generator for the `ora` gate (the urn-model test is a code check only);
  - attack 2 on the prophage short-gene layout;
  - the TU-packed attack with abundance-driven detection.
- M8-04: the O − E CI, or the CPB log-odds CI.
- M8-13: the far-tail precision and `p = 0` below ~1e-14.
- Khoa, by ADR: replace the ORA p with the conditional Poisson-binomial; make ORA TU-aware.
