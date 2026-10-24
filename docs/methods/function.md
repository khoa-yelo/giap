# Function module: detection enrichment primitives

Implements `giap.stats.wallenius` and `giap.function.elim` (GIAP_SPEC.md §7.7.3, §7.7.5, §8.1 row
"Category detection", §8.2, §8.3; issue M8-01, docs/specs/M8-01-stats-wallenius-elim.md, lead
decisions 1–4; M8-index lead decision 5). The families, the effect size and `parent_group` are
M8-04's: see "Detection enrichment" below, which replaces the Wallenius p with the TU-counted
conditional Poisson-binomial (ADR-0029). The `ora` suite and attack 2 are M8-15's.

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

## Detection enrichment

Implements `giap.function.ora`, `giap.function.summary` and `tilt_pair(comp=)` in
`giap.stats.poisson_binomial` (GIAP_SPEC.md §7.7.2, §7.7.3, §8.1 row "Category detection", §8.2;
issue M8-04, docs/specs/M8-04-category-summary-ora.md, lead decisions 1–13; ADR-0029, ADR-0035).
Outputs: the rows of `enrichment_detection.tsv` and `category_summary.tsv` (§9.3).

**Status:** PR #179. The statistician's review blocked it (B1–B6); decision 13 answers each point,
and the re-check lifted the block. Statistician sign-off: given (#179 re-check of 2026-10-22, after
edits M1–M6). Every ORA row carries `ora_uncalibrated` (see "Pass-1 limit").

### The unit: TUs
- **Eligible genes** are the non-partial genes. Only they count, for q and for detection.
- **Universe U:** the TUs with ≥ 1 eligible gene, in natural TU-ID order. Results do not depend on
  the input row order.
- **Detection:** TU t is detected (d_t = 1) when any of its eligible genes is detected.
- **Membership rule (decision 1):** TU t is in category C when any of its eligible genes is in C.
  The rule runs over all eligible genes, never over the detected ones, so membership cannot depend
  on the response.
- **Partial membership dilutes.** A TU with one member gene counts as one unit of C, and it may be
  detected through a gene outside C. This costs power, not level. Its power (spec arm 7) is owed.
- **Size filter (decision 6):** n = |C| counted in TUs must lie in `function.ora_size` ([5, 500],
  both ends included). A term outside it writes no row and does not enter elim. Under elim the
  filter applies to the full TU set; a term reduced below 5 TUs is still tested.
- **Count columns are in TUs:** `n_universe`, `n_in_set`, `n_detected`, `k_detected_in_set`.
  `n_in_set_genes` and `k_detected_in_set_genes` are the gene counts behind `p_classic`.
- **Fail closed (decision 13):** an eligible gene with a NaN `p_detect`, a `p_detect` outside
  [0, 1], or no `tu_id` raises `ValueError` naming the locus tags. So does a TU missing from a
  non-empty `folds` mapping.

### The statistic and the null
- **Statistic:** k = the detected TUs in C (`stat`).
- **Null:** the TUs are independent, and TU t is detected with probability q_t. The test is the
  law of K = Σ_{t∈C} d_t given S = Σ_{t∈U} d_t = N (the conditional Bernoulli law):

  P(K = j | N) = PB_C(j) · PB_Ū(N − j) / Σ_i PB_C(i) · PB_Ū(N − i),

  where PB_A is the Poisson-binomial pmf of the q_t in A and Ū = U ∖ C.
- **Tails:** p_over = Σ_{j ≥ k} P(K = j | N) and p_under = Σ_{j ≤ k} P(K = j | N). Both include k
  and both are sums of non-negative terms, never 1 − cdf.

**Computation (exact; no approximation, no Monte Carlo).**
1. **Truncated DP.** The convex-combination recurrence pmf[j] ← pmf[j](1 − q) + pmf[j − 1] q,
   numba-compiled and stopped at the highest index read. PB_C runs to min(|C|, N); PB_Ū is
   computed only on the band that the support reads.
2. **Tilt, once per run.** `tilt_pair(q, N, comp=)` multiplies every TU odds q/(1 − q) by one
   factor e^θ, with θ solving Σ_t expit(logit q_t + θ) = N. The solver is Newton with a bisection
   guard; it targets |Σ − N| ≤ 1e-13 N, inside the spec contract of 1e-9 N.
   - The logit is log q − log comp, and the tilted complement is expit(−(logit q + θ)). Neither
     is ever re-derived as 1 − q, so a TU with q that rounds to 1.0 keeps its complement
     (decision 11).
   - q = 0 and comp = 0 stay exact. Every open entry stays ≥ the smallest normal double, so an
     open TU is never treated as certain.
   - **The tilt is numerics only.** It puts the mass of PB_Ū near the indices that are read, so no
     weight underflows. The conditional law is invariant to a common odds factor, so the tilt does
     not change the p.
   - **Conditioning on N cancels only a common TU-level odds factor** (ADR-0035). A mis-scaling of
     q that varies with TU size, with q itself or with the category does not cancel. That is why q
     is recalibrated (next subsection).
3. **Support-sum normalisation.** The weights PB_C(j) · PB_Ū(N − j) are normalised by their sum
   over the support, not by a separate PB_U(N). This saves a DP and makes the tails sum to 1 by
   construction. A zero or non-finite sum raises `ValueError`.
4. **p domain (decision 13, review B1).** The numerator and the support sum can round apart by one
   ulp, which gave p = 1 + 2⁻⁵² and aborted BH: 663 of 4,000 strongly under-detected categories,
   and 120 of 3,000 replicates of 4b at c = 3.5. Both tails are now clamped to [0, 1], in the
   primitive and again before BH (test: tests/unit/function/test_ora_pvalue_domain.py).

**Precision contract** (tests/unit/stats/test_conditional_poisson_binomial.py: enumeration for
T ≤ 14 units, a `Fraction` DP for T ≤ 60):

| input | tails ≥ 1e-300, and the mean | smaller tails |
|---|---|---|
| untilted q | relative 1e-12 of the exact rational value | the float underflow: ≥ 0, never NaN |
| `tilt_pair` output with its complements | relative 1e-11 of the untilted exact law, including q within 1e-15 of 1 and θ > 0 | the same |

A tail of 0.0 is reachable only below about 1e-300, far below the ~1e-14 floor of scipy's
Wallenius pmf ("Far-tail precision" above).

### TU detection probability
**Clipping (decision 12).** For the null only, p̂ is clipped to [1e-12, 1 − 1e-12]. Exact 0 and 1
from M4-13 are therefore never treated as certain.
- The run, and every row whose full TU set holds a TU with a moved gene, gets `p_detect_clipped`.
- `expected` in `category_summary` and the Wallenius `odds_weight` / `p_classic` use the unclipped
  p̂. The difference is ≤ 1e-12 per gene (statistician, #179 review).

**Raw q.** Per TU, log_comp = Σ log1p(−p̂) over its eligible genes and q_raw = −expm1(log_comp)
= 1 − ∏(1 − p̂) (ADR-0029).
- q_raw assumes that the genes of a TU are detected independently given p̂.
- Shared abundance makes co-detection positive, so q_raw overstates q for multi-gene TUs, and more
  so for larger TUs.
- The ADR-0035 attack gave p̂ as the exact gene marginal over a shared TU effect of strength c.
  Before recalibration, the rejection rates were (3,000 p per cell; bounds 0.0136 / 0.058; Slurm
  1867453, 1867454):

| scenario | categories | tail | α 0.01 | α 0.05 |
|---|---|---|---|---|
| c 1.5 | large TUs | under | 0.107 | 0.248 |
| c 1.5 | single-gene TUs | over | 0.019 | 0.079 |
| c 2.5 | large TUs | under | 0.278 | 0.490 |
| c 2.5 | single-gene TUs | over | 0.044 | 0.143 |
| c 0, gene-level p̂ intercept +0.5 | large TUs | under | 0.0137 | 0.067 |

**Recalibration (ADR-0035, decision 12).**
- **Model:** d_t ~ Bernoulli(expit(β₀ + β₁ logit q_raw,t + β₂ log n_t)), where n_t is the number
  of eligible genes in TU t. The log n term is dropped when every training TU has the same size.
- **Cross-fitting:** for each of M4-13's TU folds (`fold_of_tu`), the model is fitted on the TUs
  of the other folds and predicts the held-out fold. A TU's own detection never enters its q.
- **Output:** η = x·β̂, then q = expit(η) and comp = expit(−η), each computed directly. These
  pairs are what the tilt and the conditional law receive.
- **Leakage left:** cross-fitting does not remove own-category leakage when a category is a large
  share of a covariate stratum. In the statistician's worlds this was only conservative (a power
  loss).

**Fit (decision 13, review B4).**
- An unpenalised logistic MLE by Newton, started at (0, 1, 0), the identity on logit q.
- **Damped:** each step is halved until the log-likelihood does not decrease beyond rounding
  noise (relative 1e3 · ε), up to 60 halvings.
- **Converged** when the full Newton step is below 1e-10 in every coordinate, within 50
  iterations.
- **The noise slack cannot fake convergence.** Convergence is tested on the full Newton step
  H⁻¹g at the current iterate, which is a stationarity test. A concave log-likelihood has exactly
  one stationary point, the MLE. The slack can only let a near-optimal step through, or use up
  iterations and end in a skipped fit; it cannot stop the fit away from the MLE.
- **Checked against scipy** (statistician re-check at c1e3ca0, `$GIAP_ROOT/stat/m804b/recheck/`
  `newton_check.py`, `newton_summary.txt`): 9 worlds × 40 seeds × 5 folds = 1,800 fits against
  scipy `trust-exact` plus Newton polishing. 0 fits failed. Max |Δβ| was 5.3e-15, and the max
  gradient at the returned β was 9.7e-13. The slack branch was taken in 74 folds, and every one
  ended at the reference MLE. Complete separation, quasi-separation and a collinear design return
  no fit, so the fit is skipped (`newton_adversarial.py`).
- Before decision 13 the full-step Newton diverged where the MLE exists: 2,898 of 3,000 replicates
  of the 4b + 4c world (c = 2.5, a = +0.5, b = 0.8). Each fell back to the raw q, where large-TU
  under rejected 0.49 at α 0.01.
- **Tests:** tests/unit/function/test_ora_recal_convergence.py (the divergent case) and
  test_ora_recal_oracle.py (an independent scipy fit per fold, to 1e-8 relative, including a case
  where the size term matters).

**Skipped fit (`tu_recal_skipped`).** The fit is skipped, and q = q_raw with comp = exp(log_comp),
when any of these holds:
- `folds` is None or empty. `folds` is a required keyword, so this is the caller's explicit choice
  (it is the path of M4-13's constant model);
- a fold's training TUs have fewer than 10 detected or fewer than 10 undetected TUs;
- the Hessian is singular or a step is non-finite;
- 60 halvings find no increase, or 50 iterations do not converge.

The run flag and every row get `tu_recal_skipped` (review B3).
- **How often:** at c1e3ca0, 0 of 2,000 replicates skip in any tested world, deep coverage
  included (Slurm 1867941, 1867942).
- **Why the flag matters:** a raw-q fallback over-rejects. At 3f82e0f the undamped Newton diverged
  in 417 of 3,000 deep-coverage replicates (eta0 = +1.5, about 75 % of genes detected) and fell
  back to the raw q. There, large-TU under rejected 0.45 at α 0.01, and the ~200-TU categories
  0.98 (raw q, 3f82e0f). The statistician's first review attributed those skips to `MIN_CLASS`;
  the re-check corrected this.
- **A skipped run keeps `ora_uncalibrated` even after the pass-1 flag is lifted.**

### Tests per term: directions, elim and BH
Every tested term gives one row per direction, `over` and `under`.

| column | what it is | unit | elim | enters BH |
|---|---|---|---|---|
| `p_value` | CPB tail of the tested TU set (elim-reduced for GO `over`, else the full set) | TUs | GO `over` only | yes |
| `p_noelim` | CPB tail of the full TU set | TUs | no | no |
| `p_classic` | M8-01 mean-odds Wallenius tail on the full gene set | genes | no | no |

**elim over TU sets (decision 2).**
- **Scope:** `go_bp`, `go_mf` and `go_cc`, direction `over` only. Everywhere else (`under`, `cog`,
  `kegg`, `go_slim`) `p_value = p_noelim`.
- **Walk:** M8-01's walker, fed TU index sets. A term is tested on its remaining TUs against
  U ∖ remaining, with the run's N and tilt.
- **Removal:** a term whose elim p_over is strictly below 0.01 removes the TUs it was tested on
  from every ancestor that passed the size filter.
- **Why TUs, not genes:** removing genes would leave a shared TU in the ancestor and count the
  same detection event twice (spec example 7).
- **Empty remainder:** `p_value` NA, flags `unestimated;elim_empty`; `p_noelim` and `p_classic`
  are kept.
- **Evidence** (statistician, `$GIAP_ROOT/stat/m804b/probe_elim.py` and `elim.log`, 12,000
  parents, R = 600): the elim test, with N fixed, is conservative relative to the exact law that
  also conditions on the removed child's count. Rates at α 0.01 / 0.05 / 0.10:

| world | elim `p_value` | `p_noelim` | bound |
|---|---|---|---|
| global null | 0.0048 / 0.031 / 0.070 | 0.0065 / 0.036 / 0.076 | 0.0118 / 0.054 / 0.1055 |
| one planted child per parent (logit +3, child power 0.67) | 0.0008 / 0.0073 / 0.019 | 0.0093 / 0.056 / 0.117 | same |

- **R4 label:** "BH on elim p; conditional on removing elim-significant descendants; FDR control is
  empirical."

**BH (decision 5).** One family per (collection, direction), on `p_value`, with `giap.stats.bh`.
- A NaN p stays out of its family and the row gets `unestimated`.
- `p_classic` and `p_noelim` are never adjusted.

**`parent_group` (GO only).** For a row with `q_value ≤ stats.alpha.enrichment` it is the most
general significant ancestor in the same direction and namespace (by depth, ties by ID), or the
term's own ID if it has none. Otherwise it is NA.

**Why Wallenius stays only as `p_classic` (ADR-0029).**
- Mean-odds Wallenius fails the §8.4 ORA null and attack 2 ("The approximation and the attack it
  implies" above).
- It is kept, on genes, for goseq parity (M8-13).
- If it fails closed, `p_classic` is NA with `odds_out_of_range`; the CPB columns stand.

### Effect size
- `effect` = k − E[K | N] for the tested TU set (the full set, or the elim-reduced one), in TUs.
  It is the mean of the same conditional law as the p.
- A positive effect means more detected TUs than expected given N.
- `ci_low` and `ci_high` are NA and every row carries `ci_unestimated` (decision 6, the accepted
  Q6 pass-1 waiver).
- It is not the gene-level O − E of `category_summary`.

### Degenerate cases

| situation | behaviour |
|---|---|
| C covers every TU of U | not a test: `p_value`, `p_noelim` NA, `unestimated`, out of BH |
| N = 0, or N = \|U\| | one-point support; every p = 1 |
| all q equal | the central hypergeometric on TUs |
| `detectability_constant` from M4-13 | ω = 1 for `p_classic`, but q still grows with TU size, so the TU null is not hypergeometric |
| an undetected TU with comp = 0, or a detected TU with q = 0 | the observation has probability 0 (`null_degenerate`): every row of the run has `p_value` NA with `unestimated;null_degenerate`, plus the run flag and one `UserWarning`. Reachable when comp underflows (see "Hardening owed") |
| a tail of 1 + ulp | clamped to [0, 1] |
| `p_detect` exactly 0 or 1 | clipped; `p_detect_clipped` |
| NaN `p_detect`, `p_detect` outside [0, 1], no `tu_id`, a TU without a fold | `ValueError` naming them |
| ω undefined or out of range | `odds_weight`, `p_classic` NA; `odds_out_of_range` |
| elim empties a term | `elim_empty` (above) |

### Calibration evidence
Bounds are the §8.4 pooled bound α + 2√(α(1 − α)/(R·m)), where m is the number of categories per
replicate. The categories are negatively dependent given N, so pooling is conservative.

**A. Committed arms, fast tier (R = 100), after decision 13.** These are the last records of each
stream in `local/logs/calibration/ora.jsonl` of the m8-ora worktree. The records carry no commit
or job ID. Linking them to Slurm 1867928 at 48019fc (370 passed) is an inference: the file's mtime
(09:52) lies inside that job's window (09:43–09:58). They are smoke tests, not level claims. The
fast tier does not resolve 4b or 4f: 4f passes here (0.011 against 0.0163 at α 0.01) but fails at
R = 2,000 (below). Rates are at α 0.01 / 0.05 / 0.10.

| arm | cell | over | under | bound |
|---|---|---|---|---|
| CPB pair, base (m = 20) | q ~ Beta(1.2, 2.5), q-rank blocks | 0.0065 / 0.0315 / 0.0615 | 0.0030 / 0.0155 / 0.0425 | 0.0144 / 0.0597 / 0.1134 |
| CPB pair, tilt stress (m = 20) | odds × e⁻⁸, 5 % of units at q = 1.0 in float | 0.0035 / 0.027 / 0.065 | 0.0045 / 0.020 / 0.047 | same |
| TU-ORA null (m = 10) | 15 whole TUs | 0.004 / 0.029 / 0.060 | 0.003 / 0.023 / 0.056 | 0.0163 / 0.0638 / 0.119 |
| TU-ORA null | 40 random genes | 0.003 / 0.032 / 0.062 | 0.007 / 0.030 / 0.065 | same |
| 4b, c = 1.5 | large TUs | 0 / 0.017 / 0.035 | 0.006 / 0.044 / 0.094 | same |
| 4b, c = 1.5 | single-gene TUs | 0.002 / 0.027 / 0.049 | 0.005 / 0.029 / 0.071 | same |
| 4b, c = 2.5 | large TUs | 0 / 0.014 / 0.032 | 0.015 / 0.057 / 0.104 | same |
| 4b, c = 2.5 | single-gene TUs | 0.006 / 0.023 / 0.045 | 0.007 / 0.037 / 0.080 | same |
| 4c, a = ±0.5 | large, single, random | every cell passes; the highest at α 0.01 is 0.009 (a = −0.5, single, under) | | same |
| 4f, c = 2.5, deep coverage | large TUs | 0 / 0 / 0 | 0.011 / 0.060 / 0.104 | same |

The same run fails 4d and 4e in both blocks of 100:
- **4d** (high-q ≥ 4-gene TUs, c = 1.5, m = 5), under: 0.032 and 0.028 at α 0.01 (bound 0.0189),
  0.120 and 0.102 at α 0.05 (0.0695), 0.218 and 0.180 at α 0.10 (0.1268).
- **4e** (all ≥ 4-gene TUs in two halves, c = 1.5, m = 2), under: 0.040 and 0.035 (0.0241), 0.150
  and 0.100 (0.0808), 0.280 and 0.160 (0.1424). The ≥ 2-gene thirds (m = 3), over, at α 0.10:
  0.153 and 0.127 against 0.1346, unconfirmed.

These tests pass only because every row carries `ora_uncalibrated` (the ADR-0034 fallback).

The teeth passed in the same job: the TU hypergeometric on gene categories, the raw q on 4b large
under, and the unweighted CPB on q-rank blocks. Their rates are not logged; owed.

The implementer's first log (Slurm 1867826 at dcce9e7, R = 100, 15 of 15 passed) had 4b c = 2.5
large under at 0.015 against 0.0163. That R does not resolve the cell.

**B. Resolving runs (statistician, #179 review).** Files are in `$GIAP_ROOT/stat/m804b/`
(`probe.py`, `repo_arm.py`, `summ.py`, `res_*.jsonl`). The PR model was at 3f82e0f, before the
damped Newton, in Slurm 1867861, 1867862 and 1867888 (R = 3,000). Every cell also carries the
oracle q (the exact TU marginal, same detections), which is at level in every cell below. The
conditional law is sound; the failure is the misfit of the 3-parameter recalibration.

4b, c = 2.5, large TUs, under, on the repo's own arm and seeds (`repo_arm.py`):

| block | R | α 0.01 | α 0.05 | α 0.10 | bound |
|---|---|---|---|---|---|
| 1 (replicates 0–999) | 1,000 | 0.0140 | 0.0555 | 0.1082 | 0.0120 / 0.0544 / 0.1060 |
| 2 (replicates 1,000–1,999) | 1,000 | 0.0130 | 0.0587 | 0.1098 | same |
| probe | 3,000 | 0.0130 | 0.0576 | 0.1104 | 0.0111 / 0.0525 / 0.1035 |
| probe, c = 3.5 | 3,000 | 0.0113 | 0.0552 | 0.1086 | same |

This is a confirmed fail in both ADR-0015 blocks.

The other recalibrated cells at R = 3,000, at α 0.01 / 0.05 / 0.10:

| arm | cell | rate | bound |
|---|---|---|---|
| 4d | high-q ≥ 4-gene TUs, c = 1.5, under | 0.030 / 0.108 / 0.193 | 0.0117 / 0.054 / 0.105 |
| 4e | all ≥ 4-gene TUs in halves (~115 TUs), c = 1.5, under | 0.023 / 0.114 / 0.212 | 0.0126 / 0.056 / 0.108 |
| 4e | the same, c = 2.5, under | 0.030 / 0.143 / 0.256 | same |
| 4e | all ≥ 2-gene TUs in thirds (~200 TUs), c = 1.5, over | 0.015 / 0.069 / 0.135 | 0.0121 / 0.055 / 0.106 |
| 4e | the same, c = 2.5, over | 0.018 / 0.086 / 0.167 | same |
| 4f | deep coverage, c = 2.5, large, under (at c1e3ca0 every replicate is fitted; R = 2,000, re-check below) | 0.0177 / 0.0692 / 0.1270 | 0.0120 / 0.0544 / 0.1060 per block |
| — | co-detection by TU size (c = 3 for ≥ 4 genes, 0.5 otherwise), large, under | 0.075 / 0.212 / 0.321 | 0.0111 / 0.0525 / 0.1035 |
| 4c | a = +0.5, large, under | 0.0057 / 0.0336 / 0.0656 | 0.0111 / 0.0525 / 0.1035 |
| 4c | a = −0.5, large, under | 0.0028 / 0.0186 / 0.0439 | same |

- 4c was at level in every recalibrated cell of both signs (`res_4c_a±0.5.jsonl`).
- The 4e c = 1.5 rows and the 4c rows are from `summ.py` on `res_4b_c1.5.jsonl` and
  `res_4c_a±0.5.jsonl`; the others are quoted in the review.
- In these worlds the old fit never diverged (0 skipped), so the damped Newton reaches the same
  MLE and the rates carry over. The re-check confirms this: its 4b, 4d, 4e and csize rates match.
- The 4b + 4c combination world skipped 2,898 of 3,000 at 3f82e0f. With the damped fit it is
  measured in the re-check below, and it fails.

**B′. Re-check at c1e3ca0 (statistician; `$GIAP_ROOT/stat/m804b/recheck/summ_blocks.txt`, Slurm
1867941, 1867942).** R = 2,000 per world, split into the ADR-0015 blocks 0–999 and 1,000–1,999.
0 replicates skipped in any world. These are the fallback-free p, with the flag ignored. A cell
fails when it is above the bound in both blocks. The oracle q is at level in every cell of all 13
worlds.

| world | cell | block 1 | block 2 | bound (per block) |
|---|---|---|---|---|
| 4b c 2.5 (repo arm, its own seeds) | large, under | 0.0140 / 0.0555 / 0.1082 | 0.0130 / 0.0587 / 0.1098 | 0.0120 / 0.0544 / 0.1060 |
| 4b c 2.5 (probe) | large, under | 0.0135 / 0.0586 / 0.1103 | 0.0129 / 0.0595 / 0.1118 | same |
| 4b c 3.5 | large, under | 0.0120 / 0.0525 / 0.1060 | 0.0108 / 0.0591 / 0.1113 | same (borderline) |
| 4d (4b c 1.5) | large_hiq, under | 0.0299 / 0.1075 / 0.1956 | 0.0301 / 0.1118 / 0.1979 | 0.0129 / 0.0564 / 0.1089 |
| 4b c 2.5 | large_hiq, under | 0.0430 / 0.1428 / 0.2271 | 0.0409 / 0.1523 / 0.2403 | same |
| 4e (4b c 1.5) | big_large, under | 0.0235 / 0.1215 / 0.2195 | 0.0260 / 0.1175 / 0.2135 | 0.0144 / 0.0597 / 0.1134 |
| 4e (4b c 1.5) | big_multi, over | 0.0163 / 0.0723 / 0.1340 | 0.0107 / 0.0610 / 0.1293 | 0.0136 / 0.0580 / 0.1110 (fails at 0.05 and 0.10) |
| 4b c 2.5 | big_large, under | 0.0280 / 0.1365 / 0.2510 | 0.0360 / 0.1405 / 0.2550 | 0.0144 / 0.0597 / 0.1134 |
| 4f deep (eta0 +1.5) | large, under | 0.0184 / 0.0700 / 0.1278 | 0.0171 / 0.0684 / 0.1263 | 0.0120 / 0.0544 / 0.1060 |
| 4f deep | loq_any, under | 0.0250 / 0.1000 / 0.1700 | 0.0224 / 0.0974 / 0.1792 | 0.0128 / 0.0562 / 0.1085 |
| 4f deep | big_multi, over | 0.0170 / 0.0827 / 0.1633 | 0.0157 / 0.0877 / 0.1677 | 0.0136 / 0.0580 / 0.1110 |
| combo (4b + 4c), damped fit | large, under | 0.0168 / 0.0750 / 0.1336 | 0.0206 / 0.0776 / 0.1381 | 0.0120 / 0.0544 / 0.1060 |
| combo | large_hiq, under | 0.0807 / 0.2202 / 0.3219 | 0.0810 / 0.2165 / 0.3207 | 0.0129 / 0.0564 / 0.1089 |
| csize (c 3 for ≥ 4 genes, 0.5 otherwise) | large, under | 0.0743 / 0.2116 / 0.3181 | 0.0788 / 0.2144 / 0.3249 | 0.0120 / 0.0544 / 0.1060 |
| small genome (800 genes) | hiq_any, under | 0.0140 / 0.0840 / 0.1475 | 0.0175 / 0.0795 / 0.1530 | 0.0144 / 0.0597 / 0.1134 (fails at 0.05 and 0.10) |

At level in both blocks: the null; 4c at a = ±0.5 (every cell); slope b = 1.4; the random
whole-TU and single-gene categories in every world.

**C. Flexible comparator (statistician; `probe_flex.py`, Slurm 1867892, 1867900, R = 3,000).**
- **Model:** a cross-fitted logistic on [1, logit q_raw, 4 quantile hinges of logit q_raw, a size
  one-hot {2, 3, 4, 5–6, 7+}, logit q_raw × log n], with ridge 1e-4 on the non-intercept terms.
- **Result:** at level in every cell of all 8 worlds: 4b at c = 1.5 and 2.5, the combination, the
  size-dependent c, deep coverage (including the 417 replicates in which the old undamped fit diverged at 3f82e0f), low detection and
  a small genome.
- **Examples:** 4b c = 2.5 large under 0.0053 / 0.0303 / 0.0623; the high-q cells 0.003 / 0.023 /
  0.052. The large categories are conservative.

### Pass-1 limit (decision 13; ADR-0035 decision 3)
- **Every ORA row carries `ora_uncalibrated`** (`CALIBRATION_FLAGS`).
- **Why:** arm 4b (c = 2.5, large TUs, under) fails at R = 1,000 in both blocks (0.0140 and
  0.0130 against 0.0120 at α 0.01). Its close variants 4d–4f fail by more, and so does the
  size-dependent co-detection world. These are routine category kinds within [5, 500]: high-q
  operons (ribosome-like), and categories of 100–200 multi-gene TUs.
- **What holds:** the conditional law is exact for the q it is given (oracle-q cells at level;
  precision contract). 4c and the random whole-TU and random-gene categories are at level.
- **What fails:** q itself, through the misfit of the 3-parameter recalibration (re-check table
  B′). It is anti-conservative for:
  - under-detection of large or high-q TU sets;
  - under-detection of low-q TU sets at deep coverage (`loq_any`, 0.024 at α 0.01);
  - over-detection of large multi-gene categories (`big_multi` over, 0.016 to 0.026 at α 0.01);
  - high-q categories in a small genome (800 genes).

  Any category selected on TU size or on q is at risk. Only the random whole-TU and single-gene
  categories were at level in every world. The conditional law itself is still sound: the oracle
  q is at level in every cell.
- **Consequence:** the p and q are not discovery claims. `ora.py` writes no "calibrated" claim
  (ADR-0029). A skipped fit (`tu_recal_skipped`) is a second, independent reason for the flag.

### Hardening owed (pass 2)
- **Flexible recalibration.** Adopt the statistician's comparator model (C above), or another
  that passes. It needs an ADR-0035 amendment from Khoa.
- **Lifting the flag.** Only on a pass at R = 1,000 in both blocks of every cell of 4b–4f, with
  rates logged; the fast tier is a smoke test. Skipped runs keep the flag regardless. The
  fallback in the attack files cannot trigger on its own, because the flag is a constant (backlog
  N10).
- **Comp underflow (backlog N1).** comp = exp(log_comp) underflows to 0 once log_comp < about
  −745: about 27 genes at the clip 1 − 1e-12 in one TU, or a very large TU. On the fitted path,
  expit(−η) underflows once η > about 745. An undetected TU with comp = 0 makes the whole run
  `null_degenerate`. Fix: carry log_comp into the tilt instead of comp.
- **Clip semantics (backlog N3).** State which rows count as affected. The clip also moves p̂ in
  (0, 1e-12) and (1 − 1e-12, 1).
- **The CI on `effect`.** Candidates: a TU-block bootstrap of k − E, or the conditional MLE of a
  log-odds shift δ for C's TUs with an exact CI by inverting the CPB test.
- **The `category_summary` CI for small categories** (below).
- **Power of partial membership** (spec arm 7): owed.
- **The 4b + 4c combination** fails with the damped fit (large under 0.0168 / 0.0750 / 0.1336 and
  0.0206 / 0.0776 / 0.1381 in the two blocks; large_hiq under 0.081 / 0.218 / 0.321 pooled). The
  pass-2 model must pass it, together with every other failing cell of B′, in both blocks at
  R ≥ 1,000.

### `category_summary`
- **Scope:** collections `cog`, `kegg` and `go_slim`. It is descriptive: no test, and every
  category is reported whatever its size. It is gene-level, as §9.3 names it.
- **Columns, per category:**
  - `n_genes` (eligible members) and `n_detected`;
  - `expected` = Σ p̂ over members (unclipped);
  - `o_minus_e` = `n_detected` − `expected`, with a 95 % CI;
  - `median_abundance_log2` over detected members with a finite abundance (imputed values
    included);
  - the collection's `coverage`.
- **CI:** a percentile TU block bootstrap (`giap.stats.tu_block_bootstrap`, M4-16 `basic_rows`
  construction). Per-TU sums of d − p̂ over the category's members are resampled as whole TUs and
  summed. B is the caller's (`stats.bootstrap`), and the draws come from the step's rng in row
  order.
- **Few units:** below 10 TUs (`CI_MIN_UNITS`) the CI is NA with `ci_few_units`.
- **Coverage is low for small categories** (statistician, `$GIAP_ROOT/stat/m804b/probe_boot.py`
  and `boot.log`; B = 999, 1,200 CIs per cell, true O − E = 0):

| category | coverage, c = 0 | coverage, c = 2.5 |
|---|---|---|
| 10 TUs | 0.886 | 0.887 |
| 15 TUs | 0.908 | 0.898 |
| 30 TUs | 0.922 | 0.933 |
| 100 TUs | 0.963 | 0.932 |
| 40 random genes | 0.933 | 0.942 |

  - Nominal 95 % CIs covered 0.89 to 0.93 at 10 to 30 TUs (Monte Carlo SE about 0.009).
  - Under TU-packed detection (c = 2.5) they covered 0.93 even at 100 TUs.
  - The misses are mostly intervals wholly below the truth (hi < 0), because the sum is skewed:
    8.5 % at 10 TUs against 2.5 % nominal (0.085 at c = 0, 0.087 at c = 2.5; lo > 0: 0.029 and
    0.025).
  - So an apparent under-detection (a CI entirely below 0) in a category of 30 TUs or fewer is
    weaker evidence than its nominal level.

  Raising `CI_MIN_UNITS` to about 30, or a BCa or studentized bootstrap, is owed (backlog).

## Simulator worlds for the M8-15 gate (giapsim, M8-14)

`tests/sim/giapsim/function_sim.py` and the `giapsim.scenarios` function pairs (issue M8-14, PR
#251; GIAP_SPEC.md §8.4, §8.5 attacks 2 and 3, §11.3). Statistician review of PR #251
(2026-10-23); the attack arms are in `tests/calibration/attacks/test_attack_sim_function.py`.

**Annotations.** A random GO-shaped DAG in the reserved range `GO:9000000`–`GO:9999998` (three
namespaces, `is_a` and `part_of`), 35–55 % of genes annotated directly, COG letters on synthetic
`COG<n>` IDs, 2–20 synthetic KEGG maps (`K9xxxx`, `map9xxxx`). Annotations are drawn from their
own stream, independently of every covariate, of abundance and of detection. `world.sets` keeps
sets of 10–500 TUs (the GSEA filter), so the 5–9-TU sets the ORA filter (`function.ora_size`
[5, 500]) would test are absent there; an ORA gate should build its sets through M8-03 from
`world.obo` and `world.layout.xrefs`.

**ORA null (`ora_null`).** Detection is redrawn as Bernoulli(p_true), with p_true the M4-09
covariate logistic (length, observable peptides, GRAVY). Membership is independent of detection
given p_true, so the conditional Poisson-binomial law holds exactly with q from p_true. With the
fitted M4-13 p_detect the null is the M8-04 claim itself (arm `ora_null.fitted`).

**Attack 2 (`ora_attack_confounded`).** One random GO set covering 8–20 % of genes; the
confounding is carried by `n_obs` (`CONFOUND_COVARIATE`, M8-14 lead decision 6). Each member with
an `aa` and `n_obs > 0` has its K/R residues moved, with their codons, so that its recomputed
observable-peptide count gives the logit shift nearest `member_shift`. Composition is preserved, so
GRAVY and length are unchanged. Short sORFs that cannot lose a peptide keep shift 0. The mean shift
over movable members is −1.44 to −1.50 at −1.5; about [−2, +0.5] is reachable, and a draw whose
mean misses `member_shift` by more than `SHIFT_TOL` = 0.15 raises (M8-14 decision 9).
- *Visibility:* over 200 worlds, GIAP's `covariates()` recomputed `n_observable_peptides` equal to
  `n_obs_true` for every gene, and GRAVY was unchanged exactly. A side effect: 66 % of members
  change `tm_segments` (+1.56 on average) and pI moves by up to 0.22. These are not in the truth
  logit.
- *Results at −1.5:* the oracle arm is at level and the unweighted test rejects at 1.0. The fitted
  M4-13 arm is at level and strongly conservative: 0.00 / 0.01 / 0.01 at R = 100 (Slurm 1873380).
  The default five-covariate fit, recomputed from `aa`, gives 0 / 0 / 0.005 at R = 200.
- *Why it is conservative:* 99 % of members sit at or below the non-members' 5th `n_obs`
  percentile, so the cross-fit learns the members' own rate. Their O − Σp̂ has SD 0.30 in units of
  √Σp̂(1−p̂). This comes from the plug-in variance (q̂ treated as known) and shrinks with overlap:
  SD 0.67 at −0.5, 0.55 at −0.75, 0.43 at −1.0.
- *Power:* the same plug-in costs power against a real category effect beyond the covariates. For
  an extra −0.5 logit, the fitted arm detects it in 0.19 of worlds and the oracle in 0.68 at shift
  −0.5; at −1.5 the figures are 0.01 and 0.33 (statistician probes, Slurm 1874044 and 1874094).
- *M8-15:* gate −1.5 as a separation check and −0.75 as an overlap arm (unweighted 0.95 at 0.05).

**GSEA null (`gsea_null`).** The M4-08 model (AR(1) TU levels with ρ 0.5 and SD 2, stoichiometry
0.7, gene noise 0.5) with `position_decay_log2 = 0`. Gene-level abundance is independent of
membership for any decay. The TU score (mean of the scoring genes) falls with TU size under decay,
and any-gene membership favours large TUs, so decay is a potential TU-level confounder. It is not
one at the M4-08 default: with decay 0.2, the pooled rate over every in-range set (about 268 per
world, R = 80, `SIMPLE_ENGINE`) was 0.0084 / 0.047 / 0.096, against 0.0087 / 0.045 / 0.094 at
decay 0; at 0.5 it was 0.011 / 0.060 / 0.115 (pooled over dependent sets, so indicative only).
The arm `gsea_null.decay0.2` gates one random set per world at the default decay.

**Attack 3 (`gsea_attack_operon_packed`).** Each set is 1–3 packed chr1 operons of 10–30 genes
plus `PACKED_SINGLES` = 40 single-gene TUs (at least 9 on small layouts), whole TUs only, so every
set is inside `function.gsea_size` [10, 500]. Packed TUs draw their level like every other TU; the
shared TU level, moving 10–30 genes together, is the attack.
- The packed TUs are exchangeable with the other chr1 TUs: Fligner test 0 / 0.02 / 0.08.
- The TU-collapsed GSEA is at level: 0 / 0.04 / 0.10 at R = 50.
- Gene-level ranking with gene-label permutation is a confirmed fail (Slurm 1872913, 1873380).
