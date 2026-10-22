# Operons: the within-operon distribution and the missed-operon decoy FDR

Issue M6-01 (GIAP_SPEC.md §7.4.2–§7.4.4, §8.1–§8.3), modules `giap.stats.tdist` and
`giap.stats.decoy_fdr`, with the lead decisions of `docs/specs/M6-01-lead-decisions.md`. The pair
table is M6-04, the discordance tests M6-05, the missed-operon step M6-06 and the operon-level
calibration gate M6-09. This note covers the engine; the sections the step modules add come later.

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
`converged = False` there; **at the reviewed head (326d024) the fit instead returns σ ≈ 10⁻³⁰⁵ with
`converged = True`** (50 equal values; 40 tied of 50), which would give p = 0 for every untied OP
pair downstream. Required fix before merge (statistician review of PR #140).

**Not yet provided (pass 2).** A CI for (μ, σ, ν) (profile likelihood), so the reported fit
carries its uncertainty. The effect size that §8.2 requires next to each discordance p is
z = (Δ − μ)/σ with its CI, owed by M6-05.

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
