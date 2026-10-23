# Operons: the within-operon distribution and the missed-operon decoy FDR

Issue M6-01 (GIAP_SPEC.md §7.4.2–§7.4.4, §8.1–§8.3), modules `giap.stats.tdist` and
`giap.stats.decoy_fdr`, with the lead decisions of `docs/specs/M6-01-lead-decisions.md`. The pair
table is M6-04, the discordance tests M6-05, the missed-operon step M6-06 and the operon-level
calibration gate M6-09. This note covers the engine and the discordance step (M6-05); the other step sections come
later.

## The within-operon distribution

**Model.** Δ = log2 abundance(downstream) − log2 abundance(upstream) over OP pairs whose two genes
are both detected (imputed abundances are excluded by the caller) is modelled as a Student-t with
free location μ, scale σ and degrees of freedom ν ∈ [2, 50] (`fit_student_t`, default
`df_bounds = (2, 50)`). μ is free because operon polarity makes the downstream gene systematically
lower; a heavy tail is allowed because operon partners legitimately differ 1–10-fold (§7.4).

**Fit.** Maximum likelihood, Nelder–Mead on (μ, log σ, ν) from μ₀ = median, σ₀ = 1.4826 · MAD
(floored at 10⁻⁶), ν₀ = 8, `maxiter` 2,000. The ν bound is passed to scipy's bounded Nelder–Mead
rather than through the logit reparameterisation the issue text names; the two have the same
optimum. A simplex pressed against a bound stops just short of it (ν = 49.9997), so after the fit
both bounds are scored at the fitted (μ, σ) and the best of the three ν is kept
(`_snap_to_bound`). That step can only lower the negative log-likelihood, so it cannot move the fit
away from the optimum; a ν on a bound is a legitimate optimum and is returned with
`converged = True`. `converged` is otherwise scipy's own flag (False when the optimiser stops on
iterations); M6-05 marks the tests `unestimated` when it is False.

**Two-sided p.** p = 2 · min(F, 1 − F), F = F_t((Δ − μ)/σ; ν) (`student_t_two_sided`): exactly 1 at
Δ = μ, asymmetric about 0 whenever μ ≠ 0, NaN in → NaN out, all-NaN under an unestimated fit.

**What the p means.** M6-05 fits the t on the same OP pairs it then tests, so its p is a
**plug-in** p: the parameters absorb part of each pair's own deviation. With a fit on 400 draws the
plug-in rejection rate is 0.0102 / 0.0506 / 0.1007 at α = 0.01 / 0.05 / 0.10 (R = 1,000, recorded
in `local/logs/calibration/discordance_plugin.jsonl`), so the plug-in effect is negligible at that
size; it grows as the number of OP pairs falls and is gated at the operon level in M6-09. The fit
is the null model for OP pairs only: it assumes the OP pairs are exchangeable draws from one t, so
a genome whose OP pairs mix two very different populations (e.g. highly polar long operons and
short ones) violates it and the BH family of M6-05 inherits that.

**Degenerate cases.** Fewer than 5 finite values → `StudentT(nan, nan, nan, n, False)`. When more
than ν_min/(ν_min + 1) = 2/3 of the values are tied (all equal, or a dominant exact tie, which quantised
abundances can produce), the likelihood is unbounded as σ → 0. The issue requires
`converged = False` there. The review of PR #140 found the fit returning σ ≈ 10⁻³⁰⁵ with
`converged = True` (50 equal values; 40 tied of 50), which would have given p = 0 for every untied
OP pair. Fixed in 7ad17fe (M6-01 decision 13): a scale below 10⁻⁸ · max(1, range) is a collapse and
the fit is returned unestimated, so M6-05 tests nothing.

**Not yet provided (pass 2).** A CI for (μ, σ, ν) (profile likelihood), so the reported fit
carries its uncertainty. The per-pair effect M6-05 reports is described under *Discordance*
below.

## Discordance (M6-05, pass 1)

Module `giap.operons.discordance` (GIAP_SPEC.md §7.4.3, §8.1 row "Operon discordance"; lead
decisions in `docs/specs/M6-05-discordance.md` and M6-index decision 1).

**Test.** The t of the previous section is fitted to Δ over the *usable* OP pairs: both genes
eligible (not `partial`), both detected, both log2 abundances finite, neither imputed. Fewer than
20 usable pairs, or a fit with `converged = False`, leaves every OP pair `test = none` with
`too_few_op_pairs` and `unestimated`, and no family is formed. Otherwise every eligible
both-detected OP pair with a finite Δ is tested: `p = 2 · min(F, 1 − F)` under the fit. Pairs whose
abundance is imputed are tested but not fitted (spec §7.4.2, example 5). The fit and the test use
the same pairs, so p is a plug-in p (below).

**Effect and CI.** `effect = Δ − μ̂` in log2 units: how far the pair sits from the typical
within-operon ratio. No CI is reported for it (`ci_low`, `ci_high` NA): a single pair has no
sampling replicate, and the t prediction interval is a statement about the null, not about the
effect (M6-05 lead decision 2). The claim GIAP makes is "discordant at q ≤ 0.1", and that claim
carries its q.

**Family.** One BH family per run: every OP pair with `test = both_detected_t`. In pass 1 the
one-detected branch is `unestimated (no_detection_limit_model)`, so this is literally §7.4.3's
"all tested OP pairs". When M6-05H turns on P(miss), M6-index decision 1 keeps P(miss) out of the
family and without a q. That departs from §7.4.3's "both kinds together", so it needs an ADR at
that point. P(miss) is not uniform under any null, so BH over the mixture has no justification.

**Relative scale.** With `abundance_scale = relative` (MaxLFQ, summed intensity, unknown;
§7.3.4) Δ is not a ratio of copy numbers, because each protein's normalisation constant enters it.
A pair that is "unusual" on that scale is not a stoichiometry statement. Every abundance test is
therefore `none` with `relative_scale` and `unestimated`, as the spec requires. This is the right
behaviour, not merely a conservative one.

### Plug-in level, contamination and power (statistician review of PR #147)

Simulation (2026-10-21, Slurm 1864589): Δ = −0.3 + 0.5 · t_ν, ν ∈ {5, 30}; a share π₁ of pairs
shifted by ±s log2. The shifted pairs stay in the fit, as they would in a real run. R = 2,000
replicates at n = 20 and 50, and 400 at n = 500. The convergence rate was 1.000 everywhere.
*In-sample* is the M6-05 p (fit and test on the same pairs), measured on the null pairs. *Held-out*
is the p of 50 fresh null draws under the same fit: the honest out-of-sample p, which is also what
leave-one-out would give.

| n | ν | π₁, s | in-sample rate at α = .01 / .05 / .10 | held-out rate | ν̂ median (share at 50 / at 2) |
|---|---|---|---|---|---|
| 20 | 5 | 0 | 0.0017 / 0.0540 / 0.1050 | 0.029 / 0.080 / 0.137 | 19.2 (0.47 / 0.09) |
| 20 | 30 | 0 | 0.0004 / 0.0455 / 0.1016 | 0.021 / 0.075 / 0.138 | 50 (0.70 / 0.02) |
| 50 | 5 | 0 | 0.0093 / 0.0527 / 0.1022 | 0.018 / 0.062 / 0.115 | 6.4 (0.22 / 0.02) |
| 50 | 30 | 0 | 0.0062 / 0.0474 / 0.1013 | 0.013 / 0.059 / 0.114 | 50 (0.62 / 0.00) |
| 500 | 5 | 0 | 0.0104 / 0.0507 / 0.0996 | 0.011 / 0.052 / 0.103 | 5.1 (0 / 0) |
| 500 | 30 | 0 | 0.0096 / 0.0495 / 0.1002 | 0.010 / 0.050 / 0.101 | 33.3 (0.41 / 0) |
| 20 | 5 | 0.10, 2 | 0.0003 / 0.0130 / 0.0418 | 0.003 / 0.023 / 0.062 | 2.6 |
| 50 | 5 | 0.10, 2 | 0.0004 / 0.0119 / 0.0400 | 0.001 / 0.016 / 0.046 | 2.6 |
| 500 | 5 | 0.10, 2 | 0.0006 / 0.0119 / 0.0396 | 0.001 / 0.012 / 0.040 | 2.6 |
| 500 | 5 | 0.10, 4 | 0.0002 / 0.0063 / 0.0280 | 0.000 / 0.007 / 0.028 | 2.0 (bound) |

Monte Carlo SE of the in-sample rates is ≤ 0.0009. What the table shows:

- **Double use is not anti-conservative where it matters.** The largest excess is +0.005 at
  α = 0.10 (n = 20, ν = 5), and the rate is *below* nominal at α = 0.01 at every n. The small-p
  tail is the one BH uses. The in-sample point pulls the fit toward itself, and that cancels the
  downward bias of the MLE scale at small n. The held-out p is the anti-conservative one: 0.029 at
  α = 0.01 for n = 20. So **leave-one-out would make the level worse**, not better, and it is not
  recommended.
- **True discordant pairs in the fit make the null pairs conservative.** They lower ν̂ toward the
  bound 2 and widen the tails, so contamination never inflates the level.
- **Family-wise under the complete null.** The share of replicates with any BH call at q ≤ 0.1 is
  0.0015 at n = 20, 0.0065 at n = 50 and 0.0725 at n = 500 (ν = 5); the bound is 0.1 + 2 SE.
- **TU-chained Δ.** Consecutive pairs of one TU share a gene, so their Δ are negatively
  correlated. With TUs of 2–6 genes the in-sample rates are 0.0103 / 0.0495 / 0.0993 at n = 500 and
  the any-call share is 0.0725. That passes.
- **Abundance-dependent noise** (Δ SD × exp(−0.5 · z_abundance)) passes the *pooled* criterion
  (0.0070 / 0.0522 / 0.1076 at n = 500) but **not within the lowest-abundance quartile**
  (0.024 / 0.152 / 0.271), at 2.4–3× nominal. One t for all pairs puts null calls on noisy,
  low-abundance pairs. This is attack D-1 for M6-09 (below).
- **Power is close to zero at realistic effects.** With ν free, the fit absorbs discordant pairs
  into a heavy tail. BH power at q ≤ 0.1 by n and π₁:
  - s = 2 log2 (4-fold): 0.000–0.003 at every n and every π₁ ∈ {0.05, 0.10, 0.20}.
  - s = 4 log2 (16-fold, 8 null scale units), π₁ = 0.10: 0.12 at n = 20, 0.005 at n = 50, 0.000 at
    n = 500.
  - s = 4 log2, π₁ = 0.01 (5 of 500 pairs): still 0.000, because ν̂ falls to 3.3.

  Power *falls* as n grows, because BH needs ever smaller p while ν̂ stays near the bound. This is
  a property of the specified method (§7.4.2 fits ν by ML on all OP pairs, discordant ones
  included), not of the code.

**What M6-05 can claim.** Discordance calls control the per-pair level and the family-wise false
call rate on null OP pairs. Their FDR is controlled at q ≤ 0.1, conservatively. The module will
rarely make a call: a near-empty `operon_flags.tsv` is the expected output, not evidence of
concordance. Like the missed-operon family (M6-01 decision 12), M6-09 gates level and FDR and
**claims no power**. A method change that would restore power needs an ADR draft for Khoa. Options
include: ν fixed (e.g. 4) or bounded below at 3–4; a robust or trimmed fit of the bulk; a
two-component null+outlier mixture; σ modelled as a function of mean abundance, which also
answers attack D-1.

**Minimum of 20 usable pairs.** At n = 20, ν is not identified: ν̂ sits on a bound in 56 % (ν = 5)
and 72 % (ν = 30) of replicates. The p is still at level, for the cancellation reason above, and is
conservative in the BH tail (rate 0.0017 at α = 0.01). At n = 20 a single BH call needs
p ≤ 0.005. Sizes below 20 were not measured. So 20 is defensible as a floor for *level*. It is not
a floor for useful *power*, and that is stated rather than gated.

**Attacks owed to M6-09** (`tests/calibration/attacks/`, §8.5):

- **D-1.** Δ variance that depends on abundance, with level checked per abundance quartile as
  well as pooled. The pooled check passes; the per-quartile check fails, as measured above. A
  failing attack blocks the gate, so M6-09 needs either the ADR fix or an ADR that scopes the
  claim to the pooled level.
- **D-2.** TU-chained Δ (shared genes), level and any-call share. Passes as measured above.
- **D-3.** Contamination masking: report power at π₁ ∈ {0.01, 0.05, 0.10} and s ∈ {2, 4}. This is
  a power record only; it is not gated.

## Decoy FDR for missed operons

**Targets and decoys.** Targets are NO pairs (same strand, UniOP p below threshold) with both genes
detected, scored s = −|Δ| (higher = more co-abundant). Decoys are CV (convergent) and SW
(sandwich, → ← →) pairs with both genes detected, which cannot share a transcript. Each pair
carries a distance bin b from `bin_index` (half-open, `operons.distance_bins_bp`, default
[−∞, 0, 50, 150, 300, ∞) → 5 bins; NaN distance → excluded everywhere).

**Estimator** (`reweighted_fdr`):

    FDR(t) = π₀ · Σ_b N_NO,b · (D_b(t) + 1)/(N_decoy,b + 1) / T(t)

with N_NO,b **all** targets in bin b, D_b(t) the decoys of bin b with s ≥ t and T(t) the targets
with s ≥ t. (D_b + 1)/(N_decoy,b + 1) is the per-bin null tail probability with the +1 correction,
and N_NO,b times it is the expected number of null targets of bin b above t if every target were
null; the sum is therefore a π₀ = 1 (Storey-type) decoy estimate of the false count, reweighted to
the targets' distance mix. A bin with targets and no decoys contributes N_NO,b, the maximally
conservative value. π₀ is `operons.pi0` = 1 (π̂₀ estimation is pass-2 work, needs an ADR).

**Which decoy set.** The FDR is computed with each set and the **element-wise maximum** is taken
first, then made monotone with `giap.stats.tda.monotone_q` (q at a score = the smallest FDR over
thresholds loose enough to accept it), then clipped to [0, 1] (`max_decoy_fdr`). Maximising first
gives q ≥ each set's own monotone q, so q is at least as conservative as the more conservative set
at every threshold; the reverse order can report a q below one of the two decoy FDRs (lead
decision 1). `which` names the set that gave the larger raw FDR (first key on a tie, which is
forced at the loosest threshold where every set gives π₀); it is "" for an excluded target. A set
with no usable decoy at all is skipped; if none survives, q is NaN.

**Exchangeability assumption.** Within a distance bin, |Δ| of a null NO pair and |Δ| of a decoy
pair have the same distribution. Distance is binned because co-abundance falls with distance
(shared local chromatin/supercoiling effects, read-through) and CV/SW distances differ from NO
distances (an SW distance is the sum of two gaps). The assumption is untestable on one genome; CV
pairs can look *less* alike than same-strand neighbours (transcriptional interference), which is
conservative and is why the larger FDR is reported; anything that makes decoys look *more* alike
than null NO pairs (e.g. an abundance or detection mix that differs between the classes when the
noise of Δ depends on abundance) is anti-conservative and is an M6-09 attack.

**Local fdr** (`local_fdr`, `DecoyFdr.lfdr`), for ranking only, never for calls:
lfdr(s) = min(1, π₀ f₀(|Δ|)/f_NO(|Δ|)), f₀ the Gaussian KDE of the decoy |Δ| with per-decoy weights
N_NO,b/N_decoy,b (normalised), f_NO the unweighted KDE of the target |Δ|, both reflected at 0, each
with scipy's `gaussian_kde(bw_method="silverman")` bandwidth, h = (3n/4)^(−1/5) · sd, computed from
its own **unweighted** sample (so the weights move the density, not the smoothing). Over decoy sets
the larger lfdr wins (`np.fmax`, so a set whose KDE is undefined does not mask a defined one); the
set giving the larger lfdr need not be the one in `which`. Fewer than 5 values in either sample →
NaN.

*Numerical evaluation.* The KDE sums are evaluated by linear binning on an h/40 grid, one
convolution with the kernel truncated at 10 h, and linear interpolation. Measured against the exact
reflected kernel sum (statistician review, n = 50 and 2,000, half-normal, half-t₂ and an outlier
at 40): maximum relative error 1.1 × 10⁻³ – 2.9 × 10⁻³ where the density exceeds 10⁻⁶ of its
maximum, ∫₀^∞ f = 1.000000, slope at 0 ≈ 0 (the reflection is correct). The KDE's own sampling
error is O(n^(−2/5)), i.e. several per cent even at n = 2,000, so the approximation is immaterial
for an lfdr used as a ranking. Beyond 10 h from every decoy, f₀ is exactly 0 and lfdr = 0, where
the exact KDE gives a value of order e^(−50); also immaterial.

*Known asymmetry.* A bin with targets but no decoys contributes the conservative N_NO,b to the FDR
but weight 0 to f₀, whose mass is then spread over the other bins; if that bin's null is narrower
than the others, its targets get a too-small lfdr. Harmless while lfdr only ranks; must be fixed
(lfdr NaN, or 1, for those targets) before lfdr is ever thresholded.

## Power: what this family can and cannot find

Measured by the implementer (PR #140) and re-derived and re-measured independently in the
statistician review. **Conclusion: at q ≤ 0.1 this family has no power for co-abundance effects in
the realistic range. The cause is the score, not the +1, the bins or π₀.** The calls it makes are
correct (FDR control holds); it rarely makes any.

**The identifiability bound.** Let the target |Δ| have density π₀ f₀ + (1 − π₀) f₁ (f₁ = missed
operons) and the decoys f₀. In the large-sample limit FDR(τ) at the threshold |Δ| ≤ τ is
π₀,used · F₀(τ)/F(τ). When f₁/f₀ decreases in |Δ| (true when both are centred scale families and f₁
is the narrower), F₀/F increases in τ, so its infimum is the limit τ → 0:

    min over τ of FDR = π₀,used / (π₀ + (1 − π₀) · L),    L = f₁(0)/f₀(0)   (= σ₀/σ₁ for normals)

and the true FDP of even the very best-scoring calls is at least π₀ / (π₀ + (1 − π₀) L). A call at
q ≤ α is possible only if L ≥ (π₀,used/α − π₀)/(1 − π₀); with π₀,used = 1 and α = 0.1 that is
L ≥ 46 for 20 % missed operons and L ≥ 23.5 for 40 %. Even an oracle with the true π₀ needs
L ≥ 9π₀/(1 − π₀) = 36 and 13.5. In words: a missed operon has to be ~25–50 times more concentrated
near Δ = 0 than an unrelated neighbour pair before |Δ| alone can single it out at 10 % FDR, because
the null |Δ| also has its mode at 0. The M6-09 power grid's planted effects of 0.25–1.0 SD
correspond to L ≤ 4 and are predicted to give zero calls whatever the decoy counts, bins, +1 or π₀.

**The +1 per bin (secondary).** At the top of the list, before any decoy passes, FDR ≥
π₀ Σ_b N_NO,b/(N_decoy,b + 1)/k for k calls: one phantom decoy per bin, ≈ B · N_target/N_decoy in
total for balanced bins, versus N_target/N_decoy for single-bin target–decoy. With B = 5 and
equal-size sets about 50 calls must sit above every decoy before q can reach 0.1. The sum is
dominated by the bin with the fewest decoys per target; on a real genome the overlapping-gap bin
(< 0 bp) is common for NO pairs and rare for SW pairs, so one sparse bin can hold q near 1 on its
own. This floor binds only when the bound above allows calls at all (planted spread ≲ 0.05 null SD).

**Numbers** (statistician review; 1,200 targets, 1,200 N(0, 1) decoys per set, two sets, null
targets N(0, 1), planted N(0, σ₁), 20 replicates; "calls" = mean targets at q ≤ 0.1, "FDP" = mean
FDP of those calls, "bound" = min FDR above with π₀,used = 1):

| planted share, σ₁ | bound | engine, 3 bins | 5 bins | 1 bin | 3 bins, no +1 | 3 bins, true π₀ | 3 bins, 10× decoys |
|---|---|---|---|---|---|---|---|
| 20 %, 0.02 | 0.093 | 52 (FDP 0.03) | 15 | 82 | 121 | 165 | 91 |
| 20 %, 0.05 | 0.208 | 0 | 0 | 0 | 7.7 (FDP 0.14) | 0 | 0 |
| 20 %, 0.10 | 0.357 | 0 (min q 0.36) | 0 | 0 | 1.9 (FDP 0.22) | 0 | 0 |
| 20 %, 0.25 | 0.625 | 0 (min q 0.57) | 0 | 0 | 0.8 (FDP 0.34) | 0 | 0 |
| 20 %, 0.50 | 0.833 | 0 (min q 0.73) | 0 | 0 | 0.8 (FDP 0.35) | 0 | 0 |
| 40 %, 0.02 | 0.049 | 484 (FDP 0.05) | 472 | 498 | 505 | 533 | 508 |
| 40 %, 0.05 | 0.116 | 8.4 (FDP 0.01) | 0 | 24 | 80 | 440 (FDP 0.09) | 10 |
| 40 %, 0.10 | 0.217 | 0 (min q 0.22) | 0 | 0 | 6.8 (FDP 0.14) | 0 | 0 |
| 40 %, 0.25 | 0.455 | 0 (min q 0.42) | 0 | 0 | 1.4 (FDP 0.20) | 0 | 0 |
| 40 %, 0.50 | 0.714 | 0 (min q 0.65) | 0 | 0 | 0.9 (FDP 0.23) | 0 | 0 |

The engine's mean minimum q tracks the bound (0.36 vs 0.357, 0.22 vs 0.217, slightly below it at
larger σ₁ because the minimum over a noisy curve is selected). Dropping the +1 is the only variant
that produces calls in the realistic range, and those calls are **false at 14–35 %**, i.e. the +1
is what keeps the procedure valid at small k. The M6-06 issue's example 6 (100 planted pairs with
|Δ| ≈ 0 among 2,000 NO pairs, 5 bins, ≥ 50 calls expected) is feasible only if "≈ 0" means a
planted spread ≤ about 10⁻³ null SD: at 10⁻³ the engine calls 82 of 100 (≥ 50 in 83 % of 30
replicates); at 5 × 10⁻³ it calls none.

**What may be claimed.** FDR control of the missed-operon list (gated). Not power: the operons
report and M6-09 must state that a missed-operon call requires near-identical abundance of the two
genes and that absence of calls is not evidence against missed operons. The lfdr column ranks
candidates and is the useful output of this family as specified.

**Fix space (each is a method change → ADR, Khoa decides; tolerances do not change).**

1. *Keep the method, label it* (no ADR): report calls as above, keep lfdr as the ranking, say
   "near-powerless by construction" in the report and in M6-09's power table.
2. *A stronger score* — the only option that addresses the bound. Use evidence the spec already
   computes: the likelihood ratio of Δ under the fitted within-operon t versus the decoy density
   (signed, so the operon polarity μ ≠ 0 helps), combined with UniOP p, distance and the
   promoter/terminator evidence of §7.4.5 into one score, with the same decoy FDR on top (decoys get
   the same score). The bound then applies to the combined likelihood ratio, which can be large.
   Risk: decoys must mimic targets on every covariate that enters the score.
3. *π̂₀ estimated from the decoys* (pass-2 item already listed, ADR-0019 is the precedent): gains at
   most 1/π₀ (≈ 1.25× for 20 % missed operons); does not cross the bound for L < 36.
4. *One phantom in total instead of one per bin* (e.g. the +1 carried by the bin with the largest
   weight N_NO,b/N_decoy,b): removes the B× floor at small k; needs a finite-sample argument (the
   weighted knockoff+ literature) and its own calibration including sparse-bin attacks. Gains only
   where planted spread ≲ 0.05 SD.
5. *Fewer or coarser bins* (config `operons.distance_bins_bp`, but a default change still needs an
   ADR): shrinks the floor by B; weakens the distance-confounding protection that attack A tests.
6. *A higher level for this family* (e.g. 0.25, exploratory): cannot cross the bound either
   (0.36 at 20 %, σ₁ = 0.1).
7. *More decoys*: CV and SW counts are fixed by the genome; adding other decoy classes (e.g. more
   distant opposite-strand pairs) gives more decoys but they mimic NO pairs less well. Shrinks only
   the +1 floor.

## The missed-operon step (M6-06)

`giap.operons.missed.missed_operons(pairs, cfg, *, genes=None)` applies the engine above to the
M6-04 pair table.

- **Targets.** NO pairs whose two genes are eligible (not `partial`), both detected, with finite
  `abundance_log2` and Δ. Imputed abundances and `rna_between` pairs are kept (lead decision 2),
  since §7.4.4 excludes imputation only from the t fit.
- **Decoys.** CV and SW pairs with both genes detected and finite abundances. The code applies
  `partial` to targets only. In a production gene table a partial gene is never detected
  (`detection_status = partial_excluded`), so both sides are filtered the same way in practice.
- **Configuration.** `operons.decoys` selects the sets; an empty list is a config error. A
  configured set with no usable pair is dropped with a log line.
- **Binning.** Both sides are binned by `intergenic_bp` over `operons.distance_bins_bp`, and
  `operons.pi0` is passed through. A bin with targets and no decoys of a surviving set is named in
  a warning: it contributes its whole target count, which can push every q to 1 (spec example 7).

**Output rows.** One row per NO pair. A scored target has:

- `test = missed_operon_decoy`, `stat = −|Δ|`, `effect = |Δ|`;
- `q_value` from `max_decoy_fdr`, and `lfdr` (ranking only);
- `which` (the deciding decoy set), `method = decoy_fdr`, `null_model = <sets>_reweighted`;
- `p_value`, `ci_low`, `ci_high` and `n_perm` NA (lead decision 1: a decoy FDR has no p, and a
  pseudo-p would invite a second BH).

`bin_counts` (N_NO, N_cv, N_sw per bin) goes to the run report so a reader can see the
reweighting. Nothing is scored, and every NO row is flagged, in three cases:

- the scale is relative (`relative_scale`);
- there are fewer than 20 scored NO pairs (`too_few_no_pairs`);
- no decoy set survives (`too_few_no_pairs`).

**The claim and its error estimate.** The claim is "co-abundant at q ≤ 0.1", and its error
estimate is the decoy q (§8.1). There is no CI for |Δ| of a single pair, so `ci_low` and
`ci_high` are NA; the row must say so with `ci_unestimated`, as in M6-05. The lfdr is not a claim.
In a bin with no decoys it is too small (see *Known asymmetry*), so it must not be shown as an
error estimate.

**π₀.** The default `operons.pi0 = 1` is the conservative choice. The NO class contains the missed
operons, so the true null share is below 1. A configured π₀ < 1 multiplies every q and lfdr by π₀
with no estimate behind it. Until the manifest records the π₀ used (M6-08, backlog), such a q
cannot be read correctly from `pairs.tsv` alone.

**Power.** As stated in the power section above, the step controls the FDR and will make almost no
calls at realistic effects. Its unit example 6 makes calls only because the planted pairs have
|Δ| ≤ 10⁻³ null SD. On a real genome an empty `q_value ≤ 0.1` list is the expected result. It is
not evidence against missed operons, and the lfdr ranking is the useful output.

### Level under structured nulls, and the imputation point mass (statistician review of PR #148)

**Simulation setup.**

- *Classes and bins.* Every pair is null, so any call is false. The five distance bins have
  NO/CV/SW mixes that deliberately differ: NO 0.15/0.35/0.30/0.15/0.05, CV 0.05/0.15/0.30/0.30/0.20,
  SW 0.05/0.10/0.25/0.30/0.30.
- *Sizes.* An E. coli-like genome (700 NO, 450 CV, 300 SW) and a large one (2,000 of each).
  R = 500 replicates.
- *Δ.* Δ ~ N(0, 1).
- *Imputation.* Each gene is filled at rate r, as under the default `abundance.missing = median`.
  A pair with both genes filled has Δ = 0 exactly, because the fill is one shared value. A pair
  with one gene filled has Δ = observed − fill.
- *Co-imputation.* "Co-imputation ρ" makes the partner of a filled NO gene also filled with extra
  probability ρ. This models same-strand neighbours that lose quantification together, such as
  tandem paralogs that share their peptides.
- *Script.* `local/scratch/stat-pr148/missed_sim.py` (unversioned; M6-09 turns it into attack M-1).

The share of replicates with any call at q ≤ 0.1 equals the FDR here, because every call is false.
The gate bound is 0.1 + 2 SE = 0.127. The table gives that share, with the median of the smallest
q in parentheses.

| Null | E. coli-like | large |
|---|---|---|
| no imputation | 0 (0.97) | 0 (0.96) |
| r = 5 %, the same in every class | 0 (0.96) | 0 (0.96) |
| r = 15 %, the same in every class | 0 (0.97) | 0 (0.96) |
| r = 10 % in NO, 5 % in decoys | 0 (0.97) | 0 (0.66) |
| r = 5 %, co-imputation ρ = 0.3 | 0 (0.93) | 0 (0.40) |
| r = 5 %, co-imputation ρ = 0.6 | 0 (0.78) | 0 (0.23) |
| r = 10 % in NO, 5 % in decoys, ρ = 1 | 0 (0.29) | **0.96** (0.07), ≈ 198 false calls |
| the last row with imputed pairs excluded on both sides | 0 (0.97) | 0 (0.95) |
| Δ SD ∝ exp(−0.5 a), NO abundance a shifted by +0.5 SD | 0 (0.81) | 0 (0.76) |
| the same, shifted by +1.0 SD | 0 (0.65) | 0 (0.60) |

**What the table shows.**

- **Continuous nulls stay at level.** That includes NO and decoy pairs whose abundance mix, and so
  their Δ noise, differ. This is the identifiability bound at work: with a continuous null, a false
  call needs the null NO pairs to be about 10 times tighter than the decoys.
- **Imputation removes that protection.** Pairs with both genes filled sit at Δ = 0 exactly, a point
  mass with L = f₁(0)/f₀(0) = ∞. Ties do go to the decoys, so only the *excess* of tied targets
  over the reweighted tied decoys counts. Even so, once the number of tied targets T₀ exceeds about
  (1/α) · [reweighted tied decoys + Σ_b N_NO,b/(N_decoy,b + 1)], every one of them is called. In
  the ρ = 1 row that is ≈ 200 false calls, and the per-bin +1 floor is the only margin. That
  floor is about 7–8 for the large bins and 11–19 for the E. coli-like bins, so the risk grows with
  genome size.
- **When the attack is harmless.** With equal per-gene imputation rates in every class it is
  harmless. It needs filled genes to pair up preferentially in NO pairs.
- **Excluding imputed pairs fixes it.** Excluding every pair with an imputed gene, from targets and
  decoys alike, restores the level and costs nothing measurable, because the family has no power
  at realistic effects anyway.

**Decision needed (ADR, Khoa).** Should a pair with an imputed abundance be excluded from both
targets and decoys, as the t fit already excludes it? §7.4.4 says only "both genes detected", so
exclusion is a method change.

Until that is decided, two things hold:

- the `abundance_imputed` state must reach every pair row and every call (lead decision 2 assumes
  it does);
- the §8.6 `exclude` sensitivity rerun is the built-in check: a call that vanishes under `exclude`
  rested on an imputed value.

**Attacks owed to M6-09** (`tests/calibration/attacks/`):

- **M-1.** Co-imputation excess in NO pairs, the rows above. Gated on FDR.
- **M-2.** Abundance and detection mix differing between NO and decoys, the last two rows. Gated.
- **M-3.** The M6-01 attacks: AR(1) abundance (attack A), and detection confounded with distance
  (attack B).

## Simulations (M6-03)

The operon generators live in `tests/sim/giapsim/scenarios.py` (configs) and
`tests/sim/giapsim/operons.py` (`simulate_operons(cfgs, seed)`: layout, expression truth, the §9.3
pair table the simulator can know and the truth files `truth/pairs.tsv`, `missed_operons.tsv`,
`discordant_pairs.tsv`, `operon_t.tsv`, declared in `giapsim.contract_tables.OPERON_TRUTH_TYPES`).
Common to all: half the TU boundaries flip strand (`p_same_strand = 0.5`, so CV, DV and SW decoys
are frequent), no position decay, ~80 % of genes detected (`covariates_only`, intercept −7.3; the
detection attack is the exception), and
the within-TU increment Δ = −0.3 + 0.5 · t₅ (the methods-note simulation, Slurm 1864589): each TU's
first gene keeps the merged TU + stoichiometry + noise model and each later gene adds one t draw.
A promoter sits upstream of every TU and a terminator after it; nothing inside a TU unless planted.
The M6-09 suites need ≥ 1,000 genes (1,800 give ≈ 490 OP, 300 NO, 150 CV and 80 SW pairs with both
genes detected).

| Scenario | What it plants | What it is built to break |
|---|---|---|
| `operon_null(n, rho=0)` | nothing; pure UniOP classes (true OP p ≥ 0.5, true NO p < 0.5) | level: the discordance t fit and the missed-operon decoy FDR must hold their α |
| `missed_operon_power(n, k, s)` | k true NO pairs (UniOP false negatives) with Δ ~ s · SD(NO Δ) · N(0, 1); UniOP keeps its own Beta errors; convergent TU boundaries −20–300 bp | power of the missed-operon step |
| `discordance_power(n, k, e)` | ceil(k/2) OP pairs +e over an internal promoter, floor(k/2) −e over an internal terminator | power of discordance and the §7.4.5 flag rows 1–2 |
| `operon_attack_ar1(n, rho, *, n_planted=0, effect_sd=0)` | = `operon_null(n, rho=rho)`: TU abundance AR(1) along the chromosome; optionally planted pairs as `missed_operon_power` | decoys that ignore position: neighbours look alike without co-transcription |
| `operon_attack_detection(n, type, *, n_planted=0, effect_sd=0)` | the `type` block (a third of the chromosome) has short genes (length × 0.3) and intra-TU gaps between all of them; detection depends on the final (post-walk) abundance, logit = covariates + 2 · (log2 abundance − 20) + 7 (≈ 94 % detected outside the block) and −6 on the block (≈ 49 % inside, ≈ 30 % of block pairs both detected); convergent TU boundaries −20–300 bp; optionally planted pairs as `missed_operon_power` | detection and distance confounded (ADR-0018): the both-detected filter keeps only the block's high-abundance genes, so its NO pairs, the targets of the < 50 bp bins, have smaller \|Δ\| than the CV decoys from outside the block that share those bins |

Decoy geometry: with the merged gaps no decoy pair is closer than 50 bp and SW decoys are ≥ 100 bp
apart. The opt-in `LayoutConfig.short_convergent_gaps` (on in `missed_operon_power` and the
detection attack; `convergent_gap_bp = (−20, 300)`, lead decision 23) draws every convergent
TU-boundary gap uniformly from −20–300 bp; negative gaps are 3′ overlaps of the two convergent
genes, as real genomes have. CV decoys then fill every bin up to 300 bp; off, the layouts are
byte-identical to before. The earlier default (10, 150) left the 150–300 bp bin with 1.4 CV decoys
per replicate against 188 NO targets (statistician, PR #222 round 2). Pairs per replicate by
distance bin at 1,800 genes with both genes detected (R = 20, seeds of the fast tier, commit
9873d53, Slurm 1869401, `local/logs/m6_03_cvbins.py`); a future gap change that empties a bin
shows here:

| scenario | class | < 0 bp | 0–50 | 50–150 | 150–300 | ≥ 300 |
|---|---|---|---|---|---|---|
| `missed_operon_power(1800, 50, 0.1)` | CV | 9.5 | 23.3 | 48.6 | 73.0 | 0.6 |
| | SW | 0.0 | 0.2 | 7.3 | 32.9 | 47.9 |
| | NO targets | 3.5 | 7.6 | 122.5 | 187.5 | 1.0 |
| `operon_attack_detection(1800, "prophage")` | CV | 14.8 | 33.8 | 48.6 | 69.0 | 0.4 |
| | SW | 2.4 | 6.6 | 9.8 | 31.9 | 46.0 |
| | NO targets | 16.2 | 34.8 | 116.4 | 173.7 | 1.6 |

Generator checks (tests/calibration/test_sim_operon_calibration.py and
test_sim_operon_attack_calibration.py, R = 20, seed streams `m6_03_<scenario>`; commit 9873d53,
`local/logs/m6_03_calnums_r2.py`, Slurm 1869401):
null NO-vs-decoy small-|Δ| excess 0.003 (z = 0.58); power excess 0.108 (z = 15.3); null OP tail
beyond 2 scale 0.104 vs the t tail 0.102 (z = 0.81); discordance power tail 0.168 (z = 21.9). The
convergent gaps move no abundance and no class, so these four do not depend on the gap range.
Detection attack (lead decision 12: share of |Δ| below the null's 25th percentile c among true NO
targets vs CV/SW decoys in the same distance bin, CMH one-sided): z = 1.56 (c = 0.98), with 29 % of
the block's adjacent pairs both detected (decision 10) and 48 % vs 94 % of genes detected inside vs
outside the block. At R = 300 (Slurm 1869401; the full tier of `missed_operons`, which is where
decision 12's z ≥ 2 is asserted, Slurm 1869399 passed): z = 4.08. At R = 1,000 (Slurm 1869400,
commit 9873d53; null excess and power excess at R = 300): null excess 0.001 (z = 0.86); power
excess 0.116 (z = 62.4); null OP tail 0.1011 vs 0.1019 (z = −2.05); discordance tail 0.170 (z =
158.9); detection attack z = 8.49 (c = 1.00), block both-detected 30.2 %, detected 48.5 % vs
94.1 %. Per bin at R = 1,000 (targets vs decoys below c): < 0 bp 0.362 vs 0.321, 0–50 bp 0.365 vs
0.314, ≥ 50 bp flat (0.266 vs 0.268, 0.269 vs 0.265, 0.274 vs 0.267). The R = 1,000 z implies an
expected z of about 1.2 at R = 20 and 4.7 at R = 300, so the fast tier only prints z and checks
that a sub-50-bp bin holds targets and decoys (lead decision 21).

What the attack does not show (statistician, PR #222 round 2, `bias.py`, R = 150, gaps −20–300
with UniOP errors): the generator-level excess barely reaches the decoy-FDR tail. The ratio of
actual false NO targets to the decoy-estimated count below the null's 1–25 % |Δ| quantiles is
≤ ≈ 1.05 with bin-matched CV decoys and ≈ 1.00 with one pooled bin, against ≈ 1.00 for the null,
so attack B cannot make M6-09's FDR gate fail at the FDR scale (0.1 → ≈ 0.105 at worst). A pass of
M6-09's gate B is therefore not evidence that the distance reweighting is needed or that it
defends against a confounding that breaks a simpler method.

Property of the CV decoy, kept on purpose (statistician, PR #222): CV decoys compare two TU-final
genes, while NO and SW pairs compare a TU-final gene with a TU-first gene, so the CV Δ carries one
extra within-TU walk. Real convergent pairs are TU-final on both sides too, so the method must cope
with it. Measured SD of Δ: NO 3.106, CV 3.126, SW 3.066. That is under 1 % in SD, in the
anti-conservative direction, and below the gate's resolution.

## K-12 cut-off calibration (M6-12)

Issue M6-12 (GIAP_SPEC.md §15.4 steps 3–5, §7.1.3, §7.1.5, §7.4.6), spec
`docs/specs/M6-12-uniop-terminator-cutoffs.md` with its lead decisions of 2026-10-23. Module
`giap.stats.cutoffs` holds the pure functions (`auroc`, `threshold_curve`, `choose_cutoff`,
`fit_distance_model`). `validation/calibrate_cutoffs.py` is the driver: it does the re-map,
directons, windows, matching, the bootstrap and the writers. The truth is RegulonDB 14.5.0 on
K-12 MG1655 (`U00096.3`, circular, 4,641,652 bp).

**There is no hypothesis test in M6-12.** Every output is an estimate (AUROC, precision, recall,
d₀, s) with a 95 % CI, or a cut-off chosen by a fixed rule. No p or q value is produced, so there
is no §8.1 family and nothing to correct. The M6 gate clause is "UniOP-alone AUROC **recorded**"
(§13 row M6). It has no threshold, so no test is run on it.

Notation: L is the labelled pairs (label ≠ ambiguous) that have a UniOP p. y = 1 for positive and
0 for negative. P = Σy and N = |L| − P. An item is **called at threshold t iff score ≥ t**; this
is the pipeline's comparison (`classify`, `build_tus`, the flag code).

### Directons: the resampling unit

A directon is a maximal run of consecutive genes on one strand, in genome order. It is computed on
the **whole** gene model (CDS, sORF, pseudogene; all 4,373 genes), *before* ambiguous pairs are
dropped. On a circular replicon, the last and first runs are one directon when they share a
strand. Each same-strand adjacent pair belongs to exactly one directon, the one holding both of
its genes. `n_directons` counts the directons that hold at least one pair of L, because those are
the resampling units.

Why directons, and not pairs or TUs:
- A RegulonDB TU is a run of consecutive same-strand genes, so every TU lies inside one directon.
  The dependence between pairs comes from shared TUs (one operon gives several positive pairs,
  and its boundary gives a negative). That dependence is therefore within a directon, never
  across directons.
- TUs themselves cannot be the unit. A negative pair joins two TUs, and overlapping TUs
  (alternative promoters) give pairs several TU IDs, so a pair has no unique TU.
- Directons are the finest unit that contains all the TU dependence. Coarser units, such as
  genomic windows, would add little and lose clusters.
- Across directons, the remaining dependence (UniOP's genome-wide features; RegulonDB's study
  bias towards well-studied regions) is not modelled. It is small next to the within-TU
  dependence. The directon CI is therefore the right CI, but it is not a guaranteed upper bound
  on the variance.

The resampler is `giap.stats.bootstrap.tu_block_bootstrap(df, tu_col="directon_id", stat_fn,
B=2000, rng, level=0.95)`, without strata. Each replicate draws `n_directons` directons with
replacement and keeps all their L rows. Replicates have variable size and a variable P/N mix.
That is intended: the class mix is part of the sampling variability of a single genome's truth
set.

The **single-pair comparison line** (lead decision 2) is the same function with
`tu_col = "row_id"`, a column of unique pair IDs, so each pair is its own cluster. It is
reported as `*_ci_iid`, for comparison only.

### UniOP AUROC

- **Estimate.** The Mann–Whitney AUROC = (1/(P·N)) Σ_{i∈pos} Σ_{j∈neg} ([p_i > p_j] + ½[p_i = p_j]).
- **Computation.** Use average ranks: AUROC = (R₊ − P(P+1)/2)/(P·N), where R₊ is the sum of the
  average ranks (`scipy.stats.rankdata`, method `average`) of the positives over all of L. Average
  ranks give each tie exactly ½. UniOP's output is heavily tied (153 distinct values over 2,739
  pairs), so this matters.
- **Undefined case.** P = 0 or N = 0 makes the AUROC undefined. `auroc` returns NaN.
  - In a bootstrap replicate, that NaN is dropped. With fewer than B/2 valid replicates the CI is
    (NaN, NaN), written as `null`.
  - On the full data, P = 0 or N = 0 is a `ValueError` (the "< 500 labelled pairs" guard).
- **CI.** The 95 % percentile CI of the directon bootstrap: linear quantiles at 0.025 and 0.975
  of the valid replicates (`percentile_ci`).
- **Reported.** The estimate, the CI, P, N and `n_directons`.

**Distance AUROC.** The distance-model AUROC is the same statistic on −d, with the same CI. It
does not depend on the fitted d₀ and s, because the model is monotone in d. It is the quality of
the fallback that runs when UniOP fails.

Recommended (R4): report the **paired** difference AUROC(UniOP) − AUROC(−d), with its directon
CI. It is one more `stat_fn` that computes both AUROCs on the same replicate. On K-12 the two
AUROCs are close (see the design check), and only the paired CI can say whether UniOP adds
anything over distance.

### The cut-off rule (lead decision 1: the stable reading)

**Candidates and counts.** The candidates are the distinct observed scores, sorted in descending
order: u₁ > u₂ > … > u_m. For each u_k:
- TP_k = #{y = 1, score ≥ u_k} and FP_k = #{y = 0, score ≥ u_k};
- precision_k = TP_k/(TP_k + FP_k) and recall_k = TP_k/P;
- F1_k = 2TP_k/(2TP_k + FP_k + (P − TP_k)) = 2TP_k/(TP_k + FP_k + P).

These are cumulative counts. All items tied at u_k enter together, at the exact value u_k.

**The precision test is in integers.** With the target 0.90 = 9/10, "precision_k ≥ 0.90" is
`10·TP_k ≥ 9·(TP_k + FP_k)`. For a general target a/b in lowest terms, it is
`b·TP_k ≥ a·(TP_k + FP_k)`. A precision of exactly 0.9 counts as reached. No floating-point
division is used.

**Rule (stable reading).**
- Let k* be the largest k such that every j ≤ k passes the test. This is the last candidate
  before the first failing one in descending order.
- If k* exists, t* = u_{k*} and `rule = precision_0.90`.
- If u₁ already fails, no t qualifies. Then t* = the u_k of maximum F1, and among equal F1 values
  the **smallest** u_k (the largest k); `rule = max_f1`.
  - F1 values are compared exactly. F1_k > F1_j iff TP_k·(TP_j + FP_j + P) > TP_j·(TP_k + FP_k + P),
    in Python integers or `Fraction`.

**Properties the tests rely on.**
1. Every candidate at or above t* has precision ≥ 0.90. This is what "stable" means: a
   low-threshold spike in a non-monotone curve cannot set the cut-off.
2. The stable t* is ≥ the literal t*, the smallest u_k with precision ≥ 0.90 (the literal
   qualifying set contains the stable one). The ADR prints the literal value when it differs.
   - In spec example 4 the literal value is 0.90 and the stable value is **0.92**. Under lead
     decision 1, example 4's expected t* is 0.92, with `rule = precision_0.90`, precision 1.000
     and recall 8/9.
3. The stable rule depends on the **first** candidate. If the top tie group has precision < 0.90
   (for example 4 positives and 1 negative tied at the maximum), no t qualifies and max-F1 fires,
   even if lower thresholds reach 0.90. This follows lead decision 1 as written. It needs its own
   unit test.
   - On K-12 the top UniOP group is 252 positives and 5 negatives. In 2,000 directon replicates
     the precision rule fired every time.
4. All scores equal: one candidate, AUROC = 0.5. t* is that value under either branch.
5. F1 = 0 at every candidate can only happen if P = 0, which the guard already rejects.

**Recorded value.**
- t* is an observed score, written as Python's shortest round-trip `repr` of the float, and
  **never rounded**. The pipeline then compares the same parsed float with `≥`, so every item tied
  at t* stays called.
  - TerminatorNet's candidates near 1.0 differ only in the 12th to 16th decimal (the design check
    found t* = 0.9999999999961204). Any rounding would move items across the cut-off.
- `threshold_curve` returns (u_k, TP_k, FP_k) as integers. Precision, recall and F1 are derived
  for the TSV only.

### Precision and recall at t*, and the selection effect

**UniOP.**
- Precision(t*) and recall(t*) get 95 % directon-bootstrap CIs with **t* held fixed**. Each
  replicate counts TP, FP and P at the same t*; t* is not re-selected.
- That CI is a CI for "the precision of the threshold t*", **not** for "the precision achieved by
  the rule".
  - Holding t* fixed ignores the selection. The rule chose t* *because* the observed precision
    there is ≥ 0.90, so the point estimate is biased upwards, and the fixed-t CI does not cover
    the variability of the choice.
  - The ADR must say this in those words.

**How large the selection effect is** (design check on the K-12 inputs; preliminary, not the
record):
- In-sample precision at t*: 0.9017.
- 10-fold directon cross-validation of the whole rule, choosing t on 9 folds and counting TP/FP on
  the held-out fold, pooled, 50 repeats: 0.9008 (SD over repeats 0.001).
- The out-of-bag precision of a bootstrap-selected t: 0.9025.
- So the optimism is about 0.001, against a fixed-t CI half-width of about 0.017. The selection
  effect is real but negligible here, because the precision curve is shallow around t*.
- t* itself is **not** stable: 13 % of directon replicates re-select the same candidate, and the
  95 % range of re-selected t* is [0.458, 0.610]. t* is known only to about ±0.08; the
  cut-off's *precision* is what is stable.

Recommended (R1, non-blocking, cheap): add `precision_cv`, the pooled directon 10-fold CV
precision of the rule (folds from the same `rng`, drawn after the bootstrap streams below), and
`threshold_boot_q`, the 2.5 / 50 / 97.5 % quantiles of the re-selected t* over the AUROC
replicates. Together these make the selection statement quantitative.

**What "precision 0.90" means.** The 0.90 is the share of true operon pairs among **all** pairs
called at t*: a cumulative precision. It is not the reliability of a pair near the cut-off. In
the design check, the pairs with p ∈ [0.53, 0.65) were operon pairs 50 % of the time. The ADR
should say this, because users read `pair_threshold` as a per-pair confidence.

**TerminatorNet.**
- Precision(t*) and recall(t*) get Wilson 95 % intervals, as the spec says. Precision is over
  eligible calls; recall is over eligible truth.
- The Wilson interval treats calls as independent. Calls cluster: one window can hold several
  calls, and up to 3 matched calls fall in one window. So Wilson is slightly anti-conservative.
  - In the design check, a bootstrap over windows (each gene's downstream window is a cluster)
    gave [0.194, 0.258], against Wilson's [0.197, 0.257]: about 6 % wider. That is too small to
    block on.
  - Recommended (R2): report the window-cluster interval beside Wilson, or state the 6 %.
- `recall_all` and `recall_by_type` get Wilson intervals over their own denominators.

### The d₀ / s refit

**Model.** p(operon | d) = 1/(1 + exp((d − d₀)/s)), which is y ~ Bernoulli(σ(β₀ + β₁d)) with
**d₀ = −β₀/β₁ and s = −1/β₁**. The fit uses all of L, whatever its UniOP p, with
d = `intergenic_distance(a, b, replicon)` (negative for overlaps, across the origin on a circular
replicon). The fit is unpenalised maximum likelihood.

**Algorithm.** This replaces the spec's plain Newton from (2.5, −0.05); see B1.
1. **Separation check first (exact, one covariate).** The 1-D logistic MLE exists iff the classes
   overlap in d both ways:
   - max{d : y = 1} > min{d : y = 0}, **and**
   - max{d : y = 0} > min{d : y = 1}.

   If the first fails, the data are (quasi-)separated with operons at short distance, and
   β₁ → −∞. If the second fails, the separation runs the other way. Either is a **failure**,
   reported before any iteration.
2. Start at β = (logit(P/|L|), 0). This is the intercept-only MLE, which lies inside the region
   where Newton is well behaved.
3. Damped Newton (IRLS) on the log-likelihood ℓ, computed stably with `logaddexp`.
   - Step δ = H⁻¹g.
   - Halve the step until ℓ(β + hδ) ≥ ℓ(β), with at most 30 halvings. Because ℓ is concave, the
     ascent is monotone.
   - Stop when max|hδ| < 1e-10, after at most 100 iterations.
4. **Failure** is any of the following:
   - the separation check fails;
   - no convergence in 100 iterations;
   - a singular Hessian (`LinAlgError`), which must be **caught** and reported, not raised;
   - a non-finite β;
   - β₁ ≥ 0;
   - |β₁| > 1 (s < 1 bp). This is a backstop: under quasi-separation the damped Newton can
     "converge" on saturated weights, at β₁ ≈ −10.6 in the example-5 separated data.

   A failure gives `converged: false` and keeps d₀ = 50, s = 20 with `calibrated: false`.

**CI.** The same directon bootstrap refits each replicate. A failed replicate is NaN and is
dropped (fewer than B/2 valid gives `null`). The CIs are the 95 % percentile CIs of d₀ and of s,
taken separately. On K-12 no replicate can separate: the largest positive d is 559, and the
smallest negative d is −59.
On synthetic data with known d₀ and s (part A tests, R = 300, B = 200) these percentile CIs cover
at about 0.92–0.93, slightly below 0.95 (0.89 for d₀ at B = 100). Coverage is recorded, not gated.

**Fit quality.**
- The log-loss of the refit and of the placeholder (50, 20) on L are both in-sample. The refit's
  gain is therefore optimistic by about 2/|L| nats (two parameters), which is negligible at
  |L| = 2,739. They are reported without a CI.
- The linear-logit model misfits K-12. It is too flat between 0 and 150 bp, and too confident
  beyond 300 bp, where about 10 % of the labelled pairs are still operon pairs. The model family
  is fixed by §7.1.3. The ADR should show the binned table (observed vs fitted share by distance
  bin) from the run.

### TerminatorNet

**Window eligibility.**
- Each gene of the model has one downstream window, `_window(g, neighbour, replicon, 50)`:
  - it starts at g's 3′ end, moved 50 bp back into g;
  - it ends at the nearest end of g's genome-order neighbour **on either strand**, in the
    direction of transcription;
  - it is unrolled across the origin on a circular replicon, and is 300 bp past the stop on a
    linear contig end;
  - `lo ≥ hi` gives no window.
- A call on strand s is **eligible** iff its whole interval lies in some window of a strand-s
  gene, for a shift in {0, ±length}.
- Containment, rather than overlap, keeps calls that cross into the next gene out of the
  denominator.
- Eligible calls are those in the §15.4 sense, "calls in intergenic regions downstream of genes".
  Only eligible calls enter precision. The design check found 2,447 of 8,294 calls eligible.

**The match.**
- An eligible call c matches when a RegulonDB terminator r on the same contig and strand has an
  interval gap max(0, r.start − c.end, c.start − r.end) ≤ 50 bp, for a shift in {0, ±length}.
  Overlap is a gap of 0.
- An opposite-strand terminator never matches.
- Several calls may match one terminator, and each counts (§15.4 "share of calls").
- Truth eligibility, the recall denominator, uses the same containment rule on the RegulonDB
  intervals.

**The rule.**
- Rule 2 runs on (confidence, matched) over the eligible calls. F1 in the max-F1 branch is
  therefore **call-level**: its recall is TP/(all matched eligible calls), not the truth-level
  recall that is reported.
- The pure function sees only (score, y). Truth-level F1 would need the call-to-terminator
  assignment. That choice is the lead's; the call-level F1 is the specified default.

**Why precision is a lower bound, and when it is not.**
- RegulonDB 14.5.0 lists 375 terminators on K-12, far fewer than the genome's intrinsic
  terminators, and includes Rho-dependent ones that TerminatorNet does not model.
- An unmatched call can therefore be a real terminator RegulonDB lacks. This pushes the observed
  precision **below** the true precision.
- But a false call can also land within 50 bp of a listed terminator by **coincidence**, because
  both concentrate just past stop codons. This pushes the observed precision **above** the true
  precision.
- So the observed precision is a lower bound only when missing terminators outnumber
  coincidences. This is likely, but it is not guaranteed (B2).
- The honest bound: let c be the probability that a false eligible call matches by coincidence.
  Then obs ≤ true + (1 − true)·c, so **true ≥ (obs − c)/(1 − c)**. This bound holds whatever
  RegulonDB is missing.
- c is estimated with a placement null (B2). The ADR states both: the observed precision, which
  is a lower bound to the extent that RegulonDB is incomplete, and the coincidence-corrected
  bound.

**The coincidence null (B2).**
- R = 200 placements. Each placement moves every eligible call to a uniformly chosen
  **same-strand window** that is long enough to hold it.
- The call keeps its length, its strand and its **offset from the window's 3′-side start**
  (`lo` for `+`, `hi` for `−`). Keeping the offset preserves the "just past the stop" position
  distribution that makes coincidences likely.
- The window is drawn uniformly among the same-strand windows that fit, from the same `rng`
  stream, after the bootstrap streams.
- Each placed call is matched with the same 50-bp gap rule.
- c(t) is the mean match rate over placements among calls with confidence ≥ t. The confidence
  travels with the call.
- Reported: `coincidence_rate` at t*, and `precision_lower_bound` = max(0, (obs − c)/(1 − c)).
- In the design check, a cruder uniform-in-window placement gave c ≈ 0.06–0.075 overall. The
  overall observed match rate was 0.120, so roughly half of all matches could be coincidental.
  At t* (obs 0.226) the corrected bound would be about 0.17.

**Recall is low by construction.** Rho-dependent terminators are in the truth, and TerminatorNet
does not call them. `recall_by_type` makes this visible. Recall is not a property of the cut-off
alone.

**What the rule did on K-12** (design check, preliminary):
- No t reaches 0.90 against RegulonDB. The best precision over thresholds with ≥ 20 calls is
  0.229. So the TerminatorNet cut-off is a **max-F1** cut-off, under either reading of the rule.
- t* = 0.9999999999961204, with 727 calls at or above it, precision 0.226 and F1 0.318.
- The F1 curve is flat between 0.9 and 1.0 (0.310–0.318). Re-selection in window replicates
  picks the same candidate 12 % of the time, with a 95 % range of [0.9957, 1.0].
- `terminatornet_yes` will therefore mean "saturated network output". Its 12th decimal is noise.
  - The run records it exactly, which is correct for K-12.
  - On another genome, or another TerminatorNet build, the float pattern near 1.0 can differ.
    Recommended (R3, a question for the lead): record the F1 flatness in the ADR, and consider a
    one-standard-error rule. That rule takes the smallest t whose F1 is within one bootstrap SE
    of the maximum. It would replace a noise-level argmax with a stable one, but it changes the
    §15.4 text, so it is Khoa's call.

**Implementation note for PR-B.** `parse_terminatornet` reads plain text (`Path.read_text`). The
recording is `terms.tsv.gz`, so the driver must decompress it first, or the parser must learn
gzip in a separate change. As the spec section 4 call is written, it fails with
`UnicodeDecodeError`.

### Bootstrap streams and reproducibility

The run has one generator, `rng = giap.rng.for_step("m6_12_cutoffs", 20261014)` (the config
default seed). Each `tu_block_bootstrap` call draws one (B, n_directons) block, B = 2,000. The
draw order is fixed, and the tests pin it.

1. UniOP AUROC.
2. UniOP precision at t*.
3. UniOP recall at t*.
4. Distance AUROC.
5. AUROC difference (R4, if adopted).
6. d₀.
7. s.
8. The iid comparison lines, in the same order as 1–7, with `tu_col="row_id"`.
9. The TerminatorNet coincidence placements (B2).
10. The CV folds (R1, if adopted).

d₀ and s use separate draws (two refits per replicate). That is statistically fine: each CI is a
marginal CI. It costs about 4,000 IRLS fits, which takes seconds. The same inputs and seed give a
byte-identical `cutoffs.json` apart from `date`.

### Design check on the K-12 inputs (statistician, 2026-10-23; preliminary, not the record)

This was a scratch re-implementation, run in the session allocation in under 30 s, so no Slurm
job was needed. It uses the coordinate re-map of the spec (4,372 / 4,372 matched) and the
directons above. PR-B's run is the record and replaces these numbers.
- L: 2,739 pairs (P = 1,730, N = 1,009), no pair missing a p.
  - 773 directons hold labelled pairs, with a median of 2 pairs each and a maximum of 32.
- UniOP:
  - AUROC 0.911: directon CI [0.898, 0.924], iid CI [0.899, 0.922]. The directon CI is about
    12 % wider.
  - t* = 0.529415 under both readings: precision 0.9017 [0.885, 0.918] and recall 0.848
    [0.828, 0.868] at fixed t*, with TP 1,467 and FP 160.
- Distance:
  - AUROC(−d) 0.907 [0.893, 0.920].
  - AUROC(UniOP) − AUROC(−d) = 0.004, paired directon CI [0.001, 0.008].
- Refit:
  - the MLE is d₀ = 102.7 and s = 47.2, reached in 7 steps from (logit ȳ, 0);
  - log-loss 0.407 against 0.616 for the placeholder;
  - **plain Newton from (2.5, −0.05) diverges**: the Hessian is singular at iteration 3 (B1).
- The fallback at p ≥ t* (that is, d ≤ 97 bp under the refit) has precision 0.839 and recall
  0.910. **The distance model does not deliver 0.90 precision at the UniOP cut-off.**
- TerminatorNet: as above. In addition, 112 eligible calls (4.6 %) lie within 50 bp of an
  opposite-strand terminator. That is a crude coincidence-plus-bidirectional rate, consistent
  with c.

### Limitations, and what transfers to other genomes

- **One genome, one truth set.** Every number is about K-12 MG1655 against RegulonDB 14.5.0. The
  CIs are over resampled directons of *this* genome. They describe the uncertainty of K-12's
  numbers, **not** the genome-to-genome variation users will meet. No CI here covers transfer.
- **Precision depends on prevalence.** Precision at t* is π·TPR/(π·TPR + (1 − π)·FPR), where π is
  the share of operon pairs among the scored pairs. K-12's labelled π is 0.632, with TPR 0.848
  and FPR 0.159 at t*. At the same t, a genome with π = 0.5 gets precision ≈ 0.84, and one with
  π = 0.4 gets ≈ 0.78.
  - The labelled set also excludes 10.7 % ambiguous pairs, and only covers genes RegulonDB has
    studied. Its π is not even K-12's all-pairs π.
  - `pair_threshold` is a **K-12-tuned** cut-off. Its 0.90 is not a guarantee elsewhere. AUROC is
    prevalence-free, but it is not genome-free (UniOP's features depend on the genome's
    phylogenetic neighbourhood).
- **d₀ and s transfer worse than AUROC.** β₀ absorbs K-12's prevalence. On a genome with
  different intergenic-length statistics or operon density, the fallback's probabilities are
  miscalibrated, even though its ranking (AUROC of −d) may hold.
- **Truth labels are not perfect.**
  - A "negative" pair (different TUs, none containing both) can be an operon pair under
    conditions RegulonDB has not sampled.
  - A "positive" can come from a weakly supported TU.
  - Label noise of either kind biases both AUROC and precision **towards** 0.5 and the base rate.
    So the UniOP figures are conservative with respect to this noise. They are not conservative
    with respect to study bias, because well-studied operons are probably the easy ones.
- **Possible training overlap.** If UniOP's development or training used K-12 operons, the AUROC
  is partly in-sample for the tool. This was not verified here. The ADR should state UniOP's
  published training data, or say that it is unknown.
- **TerminatorNet.** Precision is against an incomplete truth, with coincidences, as above. The
  cut-off is max-F1 on call-level F1, sits at saturated outputs, and is unstable in its
  decimals. The 0.90 target is unreachable against RegulonDB. So `terminatornet_yes` carries
  **no precision guarantee** on any genome, and the ADR must not describe it as a 90 % cut-off.
- **Selection.** Both cut-offs are chosen and evaluated on the same data. This is negligible for
  UniOP precision (≈ 0.001, see above). For TerminatorNet it is unquantified.

### Statistician review of the design (PR-A, 2026-10-23)

Verdict: **blocked on B1 and B2**, both of which this note already specifies a fix for. Once the
lead accepts the two fixes (as spec amendments or by pointing the spec at this section), the
design is **approved**. Nothing else is anti-conservative at a scale that matters.

**B1: the refit algorithm fails on the real data.** Plain Newton from the spec's start
(β₀, β₁) = (2.5, −0.05) diverges on K-12: the first step goes to β₁ = +0.039, and the Hessian is
singular at iteration 3. Under the spec's failure rule the run would ship `calibrated: false`,
although the MLE exists (d₀ ≈ 102.7, s ≈ 47.2) and is reached in 7 steps from (logit ȳ, 0). The
synthetic example 5 does not show this, because its truth (40, 15) is close to the start.
- **Fix:** the algorithm above: exact separation check, start at (logit ȳ, 0), step-halving, and
  `LinAlgError` caught as a failure. Keep |β₁| > 1 as a backstop.
- **Test (TA), a K-12-like case:** d ∈ {−20, 0, 20, …, 400}, 100 pairs at each d, with
  round(100·σ) positives from d₀ = 100, s = 50.
  - Expected: converged, d₀ = 100.14 ± 0.05 and s = 50.06 ± 0.05 (scratch values 100.136 and
    50.055).
  - Plain Newton from (2.5, −0.05) fails on these data, so the test discriminates.
- **Separation tests:** example 5 (y = 1 iff d < 60) fails at the separation check. Adding 3
  positives at d = 60 (quasi-separation: the boundary value holds both classes, but no
  y = 1 / y = 0 pair is inverted) also fails, without an exception.

**B2: "precision is a lower bound" is not true as stated.** Coincidental matches within 50 bp
push the observed precision up. A crude placement null suggests that about half of all matches
could be coincidental.
- **Fix:** the offset-preserving coincidence null above. Report `coincidence_rate` and
  `precision_lower_bound` = (obs − c)/(1 − c) in `cutoffs.json` and in the ADR. Keep the §15.4
  sentence, but qualify it: "RegulonDB is incomplete, so this precision is a lower bound, except
  for coincidental matches; the coincidence-corrected bound is …".
- **Tests (TA):**
  - On a toy genome whose calls are placed by the null itself, the observed precision equals c
    within Monte Carlo error, and the bound is ≈ 0.
  - With calls placed exactly on truth terminators, c is small and the bound is ≈ obs.

**Settled, not blocking:**
- **Directon resampling** is the right unit (see "Directons"). TUs nest in directons, there are
  773 clusters, and the CI is 12 % wider than iid on K-12. The iid comparison line stays, as the
  lead decided.
- **Selection effect for UniOP** is measured and negligible (optimism ≈ 0.001, against a
  half-width of 0.017). An honest statement with the fixed-t CI is enough. Held-out CV is
  recommended as an extra number (R1), not required.
- **The stable rule** is sound. Its one fragile point, a top tie group below 0.90, triggers max-F1
  by design. It did not occur on K-12 or in any replicate. Spec example 4 now expects 0.92.
- **Degenerate rules:**
  - one-class replicates give NaN, and fewer than B/2 valid gives `null`. This is sound and
    unreachable on K-12, where P and N are in the thousands.
  - "< 500 labelled pairs" and "> 1 % without p" are sound guards.
  - The refit failure set is sound once B1 is applied.
- **The Wilson interval for TerminatorNet** is about 6 % too narrow next to a window-cluster
  bootstrap. Report it (R2), but do not block on it.
- **Multiple testing:** no p or q values, so no §8.1 family. Every reported estimate has a CI,
  except the in-sample log-losses (descriptive), t* itself (R1 adds its bootstrap quantiles) and
  the TerminatorNet AUROC diagnostic. For the diagnostic, recommended: add a window-cluster CI,
  or mark it `ci: null` explicitly.

**Recommendations:**
- **R1:** `precision_cv` and `threshold_boot_q`.
- **R2:** a window-cluster CI for TerminatorNet precision.
- **R3:** the TerminatorNet F1 flatness in the ADR, plus a lead/Khoa question on a one-SE rule.
- **R4:** the paired UniOP − distance AUROC difference.
- **R5:** the ADR states:
  - "cumulative, not per-pair" for precision;
  - that the distance fallback reaches only ≈ 0.84 precision at `pair_threshold`;
  - that `pair_threshold` is K-12-tuned (prevalence);
  - UniOP's training data.

**Attack nulls owed (§8.5; the TA writes them under `tests/calibration/attacks/` with the PR-A
and PR-B tests):**
1. *Directon-clustered null.* Simulate genomes with K-12's directon-size distribution (1–32 pairs)
   with a shared directon random effect on the scores (ICC 0.3) and labels clustered by
   directon (a directon-level prevalence), drawn independently of each other, so the true
   AUROC is 0.5. Over R = 500 replicates, the coverage of the
   95 % directon CI must be ≥ 0.93. The iid CI's coverage is printed; it is expected to
   under-cover, which shows why the unit matters.
2. *Coincidence null* (B2, above).
3. *Top-group spike.* A top tie group at precision 0.85 above a broad region at 0.95: the stable
   rule must fall to max-F1, and the literal value must be reported beside it.
4. *Hard refit start* (B1, above), plus the separated and quasi-separated cases, with no
   exceptions raised.

## Calibration evidence

Fast tier (`tests/calibration/test_operon_fdr_calibration.py`, R = 50; 1,200 targets, 1,200 decoys
per set, 3 bins): all-null, π₀ = 1 → mean FDP 0 at every level (no calls at all); 40 % planted at
σ₁ = 0.02 → mean FDP 0.0000 / 0.0004 / 0.0479 at q = 0.01 / 0.05 / 0.10 (bound 0.1036 at 0.10),
non-vacuous (≈ 480 calls); distance-confounded bin mix, all null → reweighted FDP 0, single-bin FDP
0.02 (R = 50) / 0.03 (R = 300) at 0.10, recorded (the confounding here is too weak to break the
single-bin FDR; attack A with M6-03's generators is M6-09). lfdr under a shared null: mean median
0.981 (R = 50) / 0.985 (R = 300).

Student-t (`tests/calibration/test_operon_tdist_calibration.py`, R = 100 fast): p uniform under a
known fit (KS pass); pooled level under a known fit passes; fit on 400 draws, p of an independent
draw → 0.01 / 0.06 / 0.10 at α = 0.01 / 0.05 / 0.10 (bounds 0.030 / 0.094 / 0.160); plug-in rate as
above; parameter recovery within the asserted tolerances with > 95 % convergence.

Attacks the operon-level gate (M6-09) needs beyond the AR(1) attack A, from this review: (a) Δ noise
that depends on abundance with a different abundance (or detection) mix in NO and CV/SW pairs —
level; at engine level the bound above also caps the damage (a false call needs the null NO pairs
~10× tighter than the decoys); (b) a sparse decoy bin (SW pairs almost absent from the < 0 bp bin)
— power and the saturation of q; (c) a small genome (≈ 300 NO pairs, ≈ 150 decoys per set), where
the +1 floor dominates.

## Sign-off

Statistician (PR #140 review, 2026-10-21): FDR control and the KDE evaluation **approved**; the
zero-MAD / dominant-tie fit (σ → 0 reported as converged) **must be fixed** before merge; the power
limitation above is a property of the specified method and is reported, not a merge blocker.

Statistician (PR #148 review, 2026-10-21): the missed-operon step is at level under every
continuous null tried, including an abundance mix that differs between NO and decoy pairs. Pairs
with both genes imputed share a fill value, so their Δ = 0 is a point mass that can produce
wholesale false calls when it is more common among NO pairs. The imputed state must travel with the
pair rows now, and exclusion on both sides is an ADR for Khoa. Relative-scale rows need
`unestimated` (§7.3.4), and scored rows need `ci_unestimated`.
Statistician (PR #147 review, 2026-10-21): the discordance test is at level with the double use
included, and BH over the both-detected t p-values is the specified family. The flag on tested rows
must say that the CI, not the claim, is unestimated. Power is near zero at realistic effects; that
is reported here, not a merge blocker. Attack D-1 is owed to M6-09.

Statistician (M6-12 PR-A design review, 2026-10-23): **blocked on B1** (the refit diverges on
K-12 from the spec's start; fixed by the separation check, the (logit ȳ, 0) start and step-halving)
**and B2** (the TerminatorNet "lower bound" ignores coincidental matches; fixed by the
offset-preserving coincidence null and the corrected bound). Both fixes are specified in "K-12
cut-off calibration (M6-12)"; once the lead accepts them the design is approved. Directon
resampling, the stable rule, the fixed-t CIs with the selection statement, and the degenerate rules
are sound.
