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
| `missed_operon_power(n, k, s)` | k true NO pairs (UniOP false negatives) with Δ ~ s · SD(NO Δ) · N(0, 1); UniOP keeps its own Beta errors; convergent TU boundaries 10–150 bp | power of the missed-operon step |
| `discordance_power(n, k, e)` | ceil(k/2) OP pairs +e over an internal promoter, floor(k/2) −e over an internal terminator | power of discordance and the §7.4.5 flag rows 1–2 |
| `operon_attack_ar1(n, rho, *, n_planted=0, effect_sd=0)` | = `operon_null(n, rho=rho)`: TU abundance AR(1) along the chromosome; optionally planted pairs as `missed_operon_power` | decoys that ignore position: neighbours look alike without co-transcription |
| `operon_attack_detection(n, type, *, n_planted=0, effect_sd=0)` | the `type` block (a third of the chromosome) has short genes (length × 0.3) and intra-TU gaps between all of them; detection depends on the final (post-walk) abundance, logit = covariates + 2 · (log2 abundance − 20) + 7 (≈ 94 % detected outside the block) and −6 on the block (≈ 49 % inside, ≈ 30 % of block pairs both detected); convergent TU boundaries 10–150 bp; optionally planted pairs as `missed_operon_power` | detection and distance confounded (ADR-0018): the both-detected filter keeps only the block's high-abundance genes, so its NO pairs, the targets of the < 50 bp bins, have smaller \|Δ\| than the CV decoys from outside the block that share those bins |

Decoy geometry: with the merged gaps no decoy pair is closer than 50 bp and SW decoys are ≥ 100 bp
apart. The opt-in `LayoutConfig.short_convergent_gaps` (on in `missed_operon_power` and the
detection attack; `convergent_gap_bp = (10, 150)`) draws every convergent TU-boundary gap uniformly
from 10–150 bp, so CV decoys reach the < 50 bp bin; off, the layouts are byte-identical to before.

Generator checks (tests/calibration/test_sim_operon_calibration.py, R = 20, seed streams
`m6_03_<scenario>`): null NO-vs-decoy small-|Δ| excess 0.003 (z = 0.58); power excess 0.108
(z = 15.3); null OP tail beyond 2 scale 0.104 vs the t tail 0.102 (z = 0.81); discordance power
tail 0.168 (z = 21.9). Known approximation: CV decoys compare two TU-final genes while NO and SW
pairs compare a TU-final with a TU-first gene, so the CV Δ carries one extra within-TU walk; at
`tu_sd = 2` this is ≈ 2 % of the small-|Δ| share, below the gate's resolution.

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
