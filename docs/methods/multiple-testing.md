# Multiple testing: Benjamini–Hochberg

Implements `giap.stats.multitest.bh` (GIAP_SPEC.md §8.2, §8.3, Appendix B.1; issue M4-01).
Statistician sign-off: see the end of this note.

## Definition

For the m non-missing p-values sorted ascending, p₍₁₎ ≤ … ≤ p₍ₘ₎,

    q₍ᵢ₎ = min( 1, min_{j ≥ i} m · p₍ⱼ₎ / j ).

Rejecting every hypothesis with q ≤ α is the Benjamini–Hochberg (1995) step-up procedure at level α.
q is returned in input order (never sorted).

- **Monotone.** The running minimum from the largest p makes q non-decreasing in p, so the q ≤ α
  set is exactly the BH rejection set.
- **Ties.** Tied p-values get equal q (the running minimum assigns each tie group the value at its
  largest rank). The result does not depend on input order.
- **Missing values.** NaN p-values stay NaN and are **not** counted in m. A NaN means "no test was
  possible" (for example, no permutation was run), not "p = 1". Callers who want such units to
  dilute the family must pass 1.0 explicitly.
- **Range.** Any non-NaN p outside [0, 1] (including ±inf) raises `ValueError`; clipping upstream is
  the caller's job, because silent clipping hides bugs.
- **Empty input** returns an empty array.

## What is and is not controlled

- BH controls FDR ≤ π₀·α ≤ α for independent p-values and under positive regression dependence
  (PRDS; Benjamini & Yekutieli 2001), which covers one-sided tests on positively correlated
  Gaussian statistics. It does **not** guarantee FDR under arbitrary dependence (BY would be
  needed); GIAP's families are chosen so that PRDS is plausible or so that resampling respects the
  dependence (TU blocks, circular shifts).
- No π₀ estimation is done (no Storey adjustment): BH is conservative by the factor π₀.
- Discrete p-values (permutation p = (1 + b)/(1 + B), rounded p) keep control; BH is conservative
  for them.
- BH guarantees nothing across families. Error control is per family (§8.1), so a reader scanning
  several families should expect roughly α × (number of discoveries) false ones per family.

## Families (§8.1, pre-declared)

BH is applied **within** each family and never across them:

| Family | BH unit | Level |
|---|---|---|
| Region effect | region, within region type | 0.05 |
| Autocorrelation | lag, across lags | 0.05 |
| Operon discordance | pair | 0.1 |
| Category detection (ORA) | category, per ontology × direction | 0.05 |
| GSEA | gene set, per collection | 0.05 |

The family keys are fixed in the calling modules (M4-17/18, M5, M6, M8), which must cite this
table. Islands use FWER (maximum statistic) and identification/placement use target–decoy q-values
(see `target-decoy.md`), not BH. All p and q values are reported, whether or not they pass.

## Simulation evidence

Generator (tests/calibration/test_stats_fdr_calibration.py): N = 2000 tests, π₀ = 0.8, null
p ~ U(0, 1), alternative p ~ Beta(0.1, 8). Criterion (§8.4): mean FDP ≤ q + 2 SE. Seed
`calibration_base` = 20261014.

| R | q = 0.01 | q = 0.05 | q = 0.10 | mean discoveries (q = 0.05) |
|---|---|---|---|---|
| 50 (CI) | 0.0094 ± 0.0010 | 0.0428 ± 0.0016 | 0.0790 ± 0.0017 | 331 |
| 300 | 0.0081 ± 0.0003 | 0.0406 ± 0.0006 | 0.0800 ± 0.0008 | 331 |

Observed FDP ≈ π₀·q, as theory predicts.

Attack (tests/calibration/attacks/test_attack_fdr_engine.py, `bh_equicorrelated_ties`): N = 1000,
one-sided z with pairwise correlation 0.5 and p rounded to three decimals (heavy ties), R = 300:
FDP 0.0049 ± 0.0015, 0.0281 ± 0.0046, 0.0592 ± 0.0073 at q = 0.01, 0.05, 0.10 — at level.

Exactness: agrees with `statsmodels.stats.multitest.multipletests(method="fdr_bh")` to rel. 1e-12
(tests/property) and reproduces Appendix B.1.

## Sign-off

Statistician (M4-01 review, 2026-10-14): **approved** for `bh` as implemented in
`src/giap/stats/multitest.py`.
