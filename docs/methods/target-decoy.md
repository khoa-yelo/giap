# Target–decoy q-values and picked FDR

Implements `giap.stats.tda.{monotone_q, tda_qvalues, picked_fdr}` (GIAP_SPEC.md §7.2.5, §7.3.2,
§8.3, Appendix B.6; issue M4-01). Used by placement FDR (M3-10, per stratum) and the
`decoy_picked` gene-detection tier (M4-11/12).

## Conventions

- Higher score = better. Callers with "lower is better" statistics pass the negation (placement
  passes −E).
- NaN and −inf mean **absent** (no evidence). +inf is a legal, best-possible score (e.g. a rank
  transform of q = 0).
- Thresholds are the distinct observed target scores. At threshold t, T(t) and D(t) count every
  target and decoy with score **≥ t**. A decoy tied with a target is therefore counted against it
  (the conservative reading, and consistent with "ties go to the decoy" in picking).
- Output is one value per input target, in input order. Tied targets get equal q.

## Separate search: `tda_qvalues(target, decoy, pi=1, plus_one=True)`

    FDR(t) = π · (D(t) + 1) / T(t)        (π · D(t)/T(t) if plus_one=False)
    qᵢ     = min_{t ≤ sᵢ} FDR(t), clipped to [0, 1]

`monotone_q` does the second step: given FDR at thresholds sorted strictest → loosest, it returns
the running minimum from the loose end. So a target's q is the smallest FDR among all thresholds
that would still accept it, and accepting q ≤ α is the same as taking the loosest threshold with
FDR ≤ α.

- **The +1.** It is the "knockoff+" correction (Barber & Candès 2015; He et al. 2015): it stops
  D(t) = 0 from producing FDR = 0 when few targets pass. It makes small strata conservative by
  design: with no decoys at all, q = 1/T, so a stratum needs T ≥ 100 targets above the first decoy
  to accept anything at 0.01 (B.6: 150 targets → q = 0.0067 accepted; 50 targets → q = 0.02, none
  accepted).
- **π.** π ∈ (0, 1] scales the decoy count to the expected share of wrong targets (π̂_wrong in
  §7.2.5; 1 for m ≥ 1 strata). π outside (0, 1] (or NaN) raises `ValueError`. π = 1 is always
  conservative when there are as many decoys as targets; π < 1 is valid **only** if π is a
  conservative (upward-biased) estimate of the true wrong share — see the violation table below.
- **Assumption.** Decoys must be a sample from the null-target score distribution with the same
  number of draws as targets (one decoy search per target search). Separate search has no
  finite-sample guarantee; the guarantee is asymptotic plus the +1 margin, and the simulations below
  show it at level.

## Picked competition: `picked_fdr(target_best, decoy_best)`

Savitski et al. (2015), applied to genes (§7.3.2). Per gene g, with target score s_g and paired
decoy score s̃_g:

1. **Pick.** The decoy wins when s̃_g ≥ s_g (**a tie counts as a decoy**) or when the gene has only
   decoy evidence. The target wins when s_g > s̃_g or the gene has only target evidence. A gene with
   neither is in no list.
2. **FDR and q.** q over the picked list = `tda_qvalues(picked targets, picked decoys)` with π = 1
   and +1: FDR(t) = (#picked decoys ≥ t + 1)/#picked targets ≥ t, made monotone.
3. **Output.** `q` is defined only for target winners; it is NaN where the decoy wins or neither
   exists (the gene is then `not_detected`, decided in M4-12). `decoy_wins` is true also for
   decoy-only genes; `n_target_wins`, `n_decoy_wins` are reported counts.

Why this controls FDR: for a null gene whose target and decoy scores are exchangeable, the picked
label is a fair coin independent of the picked score, so picking is target–decoy competition and the
(D + 1)/T rule with the loosest passing threshold is the knockoff+ procedure, which controls FDR in
finite samples. A tie is a case where the coin cannot be read; giving it to the decoy keeps the
coin at most fair for targets. Decoy-only genes count as decoys, which only lowers T/D and is
conservative.

## What is and is not controlled

Controlled (given the assumptions): the expected share of false entries among targets accepted at
q ≤ α — placement FDR within a stratum, gene FDR within the picked list.

Not controlled:

- FDR for anything aggregated afterwards (genes from peptide q, operons from gene q, §8.7).
- FDR when the decoys do not mimic null targets (wrong score distribution, missing decoys,
  different multiplicity). These are input-contract violations; the engine cannot detect them.
- FWER: one false call at q ≤ 0.01 is expected about once per 100 calls.
- Per-gene error: q is a set-level quantity; a gene with q = 0.009 is not "99.1 % certain".

## Simulation evidence

Fast suite (tests/calibration/test_stats_fdr_calibration.py): N = 2000, π₀ = 0.8, null scores
N(0, 1), true scores N(3.5, 1), one N(0, 1) decoy per target; picked adds 5 % of genes without a
target and 5 % without a decoy. Mean FDP ± SE, seed 20261014:

| Method | R | q = 0.01 | q = 0.05 | q = 0.10 |
|---|---|---|---|---|
| separate, π = 1 | 50 | 0.0076 ± 0.0010 | 0.0392 ± 0.0021 | 0.0770 ± 0.0029 |
| separate, π = 1 | 300 | 0.0070 ± 0.0004 | 0.0386 ± 0.0008 | 0.0771 ± 0.0011 |
| separate, π = true π₀ | 50 | 0.0093 ± 0.0011 | 0.0497 ± 0.0025 | 0.0960 ± 0.0029 |
| separate, π = true π₀ | 300 | 0.0088 ± 0.0005 | 0.0485 ± 0.0009 | 0.0965 ± 0.0011 |
| separate, π = 1, no +1 | 300 | 0.0097 ± 0.0005 | 0.0409 ± 0.0008 | 0.0791 ± 0.0011 |
| picked | 50 | 0.0070 ± 0.0010 | 0.0431 ± 0.0022 | 0.0949 ± 0.0030 |
| picked | 300 | 0.0076 ± 0.0005 | 0.0438 ± 0.0009 | 0.0926 ± 0.0012 |

All at level; π = 1 is conservative by ≈ π₀, as expected; without the +1 the q = 0.01 margin
shrinks to the SE.

### Attack nulls (tests/calibration/attacks/test_attack_fdr_engine.py, R = 50 and 300)

| Attack | Built to break | R = 300 mean FDP at q = 0.01 / 0.05 / 0.10 |
|---|---|---|
| picked, gene sizes log-normal (best of k_g), shared target–decoy gene effect (SD 0.7), scores rounded to 0.1, 5 % target-absent, 5 % decoy −inf | heterogeneous multiplicity, correlation, ties | 0.0048 / 0.0408 / 0.0859 |
| picked, integer scores (most genes tie) | tie rule | 0.0017 / 0.0258 / 0.0286 |
| separate, exponential scores on a 0.5 grid + shared batch shift | ties, dependence | 0.0000 / 0.0142 / 0.0536 |
| global null (every call false), picked and separate, n = 2000 and n = 25 | P(any call) = FDR | ≤ 0.0033 everywhere |

All stay at level with the current implementation.

### Assumption violations (documented, not CI gates)

These break the inputs, not the code, and show what upstream modules must guarantee (R = 300,
N = 2000, π₀ = 0.8):

| Violation | q = 0.01 / 0.05 / 0.10 | Consequence for callers |
|---|---|---|
| tie given to the **target** (integer scores) | 0.029 / 0.214 / 0.233 | the ≥ / tie-to-decoy rule is essential; never break ties in the target's favour upstream |
| decoys shifted −0.3 SD below null targets | 0.019 / 0.097 / 0.184 | decoys must mimic targets (§8.2): same length, composition, multiplicity |
| 10 % of decoys lost at random | 0.0085 / 0.052 / 0.106 | borderline |
| 20 % of decoys lost at random | 0.0099 / 0.059 / 0.122 | **the §7.3.2 fallback at > 20 % unpaired decoys is too lax**: at 20 % loss gene FDR is ≈ 1.2 × nominal. Recommendation for M4-11/12: scale D by 1/(paired share) or fall back at a lower share (≈ 5 %) |
| 30 % of null-gene decoys lost | 0.011 / 0.067 / 0.138 | as above |
| π = 0.5 when true π₀ = 0.8 | 0.015 / 0.077 / 0.158 | π̂ must be conservative; M3-10's π̂_wrong (mean PEP, max q, de novo lower bound) needs its own calibration in the placement suite |

## Convergent and sandwich decoys (operons, M6-01)

`giap.stats.decoy_fdr` reuses `monotone_q` for the missed-operon family (§7.4.4); the full method
and its power analysis are in `docs/methods/operons.md`. What differs from the separate search
above:

- **Decoys are not one per target.** CV (convergent) and SW (→ ← →) pairs are whatever the genome
  has, in different numbers and with a different intergenic-distance mix than the NO targets. The
  "one decoy per target, same distribution" assumption of `tda_qvalues` is replaced by
  **within-bin exchangeability**: in each distance bin b, |Δ| of a null NO pair and of a decoy have
  the same distribution. Each bin's decoy tail fraction (D_b + 1)/(N_decoy,b + 1) is scaled by the
  bin's target count N_NO,b, so the decoys are reweighted to the targets' distance mix and their
  number no longer matters.
- **The +1 is per bin**, so the small-list floor is ≈ π₀ Σ_b N_NO,b/(N_decoy,b + 1) phantom false
  targets — about B times the single-bin floor, dominated by the bin with the fewest decoys per
  target. A bin with targets and no decoys contributes all its targets.
- **Two decoy sets, the larger FDR wins**, maximum taken before monotonising. CV partners can be
  *less* alike than same-strand neighbours (transcriptional interference), which is conservative;
  SW pairs are farther apart. Reporting the larger guards against whichever set happens to mimic
  null NO pairs badly in the conservative direction; it does not guard against a set that is
  *more* alike than null targets (anti-conservative), which only an attack null can reveal.
- **π₀ = 1** and N_NO,b counts every target, so the estimate is conservative by ≈ π₀.

## Numerical notes

- Counts use `searchsorted` on sorted arrays (O(n log n)); no floating-point equality beyond the
  exact ≥ comparison of the caller's scores. Callers that compute scores on different paths for
  targets and decoys must round consistently, or a near-tie becomes a target win.
- `monotone_q` assumes no NaN in its input (a NaN propagates to all stricter thresholds); the
  public functions never pass one.

## Sign-off

Statistician (M4-01 review, 2026-10-14): **approved** for `tda_qvalues`, `picked_fdr` and
`monotone_q` as implemented in `src/giap/stats/tda.py`, with the decoy-loss finding forwarded to
M4-11/12.
