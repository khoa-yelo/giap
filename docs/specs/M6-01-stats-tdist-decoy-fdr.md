# M6-01: Within-operon Student-t, distance-reweighted decoy FDR, monotone q and local fdr

Milestone M6 · track S (statistics engine) · **pass 1** · depends on: nothing (reuses merged
`giap.stats.bh`) · estimate ≈ 260 code lines (+ ≈ 320 test lines) · **statistical** (label `stats`)

## Spec sections implemented

- §7.4.2: "Fit a Student-t (location, scale and df, with df in [2, 50]) by maximum likelihood to Δ
  over OP pairs where both genes are detected, starting from the median and MAD. Pairs with imputed
  abundance are excluded from the fit."
- §7.4.3 first bullet: "Two-sided p = 2 · min(F_t(Δ), 1 − F_t(Δ)) under the fitted t."
- §7.4.4: the decoy FDR formula with distance bins and the +1 correction, "report the **larger**"
  of the two decoy FDRs, "make the FDR monotone into q-values", and
  "lfdr(s) = min(1, π₀ f₀(s) / f_NO(s))" with "Gaussian KDE, boundary reflection at 0, Silverman
  bandwidth".
- §8.3 (engine functions are pure, typed, with a methods note), §8.2 ("decoys mimic targets").

## Public interface

`src/giap/stats/tdist.py` (pure)
```python
@dataclass(frozen=True, slots=True)
class StudentT:
    loc: float; scale: float; df: float; n: int; converged: bool

def fit_student_t(x: ArrayLike, df_bounds: tuple[float, float] = (2.0, 50.0)) -> StudentT: ...
def student_t_two_sided(x: ArrayLike, fit: StudentT) -> F64: ...
```
`src/giap/stats/decoy_fdr.py` (pure)
```python
@dataclass(frozen=True, slots=True)
class DecoyFdr:
    q: F64                  # per target, monotone, in [0, 1]
    fdr: F64                # the raw step function at each target score
    which: npt.NDArray[np.str_]   # the decoy set that gave the larger FDR, per target
    lfdr: F64               # NaN when the KDE is not defined

def reweighted_fdr(target_s: ArrayLike, target_bin: ArrayLike, decoy_s: ArrayLike,
                   decoy_bin: ArrayLike, n_bins: int, pi0: float = 1.0) -> F64: ...
def max_decoy_fdr(target_s, target_bin, decoys: Mapping[str, tuple[ArrayLike, ArrayLike]],
                  n_bins: int, pi0: float = 1.0) -> DecoyFdr: ...
def local_fdr(target_s, decoy_s, decoy_w, pi0: float = 1.0) -> F64: ...
def bin_index(distance_bp: ArrayLike, edges: Sequence[float]) -> npt.NDArray[np.int64]: ...
```

1. **`fit_student_t`.** Maximum likelihood over (loc, scale, df) with `scale > 0` and df clamped to
   `df_bounds`; start at `loc = median(x)`, `scale = 1.4826 · MAD(x)`, `df = 8`. Optimiser:
   Nelder–Mead on (loc, log scale, df with the bounds given to the optimiser, then both bounds scored), `maxiter = 2000`. NaNs are
   dropped. `converged = False` when the optimiser stops on iterations; the caller decides what to
   do (M6-05 marks the tests `unestimated`).
2. **`student_t_two_sided`.** `2 · min(F, 1 − F)` with `F = t.cdf((x − loc)/scale, df)`; exactly
   1.0 at `x = loc`; NaN in → NaN out.
3. **`bin_index`.** Half-open bins `[edges[i], edges[i+1])` from `operons.distance_bins_bp`
   (`[-inf, 0, 50, 150, 300, inf]` → 5 bins); NaN distance → bin −1, excluded everywhere.
4. **`reweighted_fdr`.** For each distinct target score t, descending,
   `FDR(t) = π₀ · Σ_b N_NO,b · (#{decoy_b : s ≥ t} + 1)/(N_decoy,b + 1) / #{target : s ≥ t}`,
   where `N_NO,b` counts **targets** in bin b (all of them, not only those above t) — the spec's
   `N_{NO,b}`. Bins with no target contribute nothing; bins with no decoy contribute `N_NO,b · 1`,
   which is the maximally conservative value and must be kept (example 4).
5. **`max_decoy_fdr`.** Runs step 4 once per decoy set, takes the **element-wise maximum**, records
   which set won in `which`, then makes it monotone non-increasing in the score by a running minimum
   from the top, and clips to [0, 1] → `q`.
6. **`local_fdr`.** `f₀` = Gaussian KDE of the decoy `|Δ|` values with per-point weights `decoy_w`
   (the bin reweighting of step 4, normalised), reflected at 0; `f_NO` = unweighted Gaussian KDE of
   the target `|Δ|`, also reflected; Silverman bandwidth from the unweighted samples;
   `lfdr = min(1, π₀ f₀ / f_NO)` evaluated at each target, NaN where `f_NO = 0`. Scores are
   `s = −|Δ|`, so the KDE works on `|Δ| = −s`, which is the bounded quantity.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| `n < 5` finite values in `fit_student_t` | `StudentT(nan, nan, nan, n, False)`; no exception |
| zero MAD (all Δ equal) | start `scale = 1e-6`; if the likelihood is unbounded, return `converged = False` |
| df hits a bound | returned at the bound with `converged = True` (a bound is a legitimate optimum) |
| no target scores | `DecoyFdr` with empty arrays |
| no decoy at all in a decoy set | that set is skipped; if no set survives, `q` is all NaN |
| a bin with targets but no decoys | contributes `N_NO,b`; `q` may saturate at 1 (example 4) |
| ties in target scores | all tied targets get the same `fdr`, `q`, `which` |
| < 5 values for either KDE | `lfdr` all NaN |
| `pi0` outside (0, 1] | `ValueError` |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `fit_student_t` on 2,000 draws from t(loc 0.3, scale 0.5, df 5), seed fixed | `loc` 0.3 ± 0.05, `scale` 0.5 ± 0.05, `df` 5 ± 2, `converged` |
| 2 | Δ = 1.2 under `StudentT(0.05, 0.4, 6)` | z = 2.875, p = **0.0282** |
| 3 | Δ = −1.2 under the same fit | p = **0.0205** (asymmetric because `loc ≠ 0`) |
| 4 | targets 200 in bin 0 and 100 in bin 1; CV decoys 100 and 100; at t: 40 targets, 5 and 2 decoys ≥ t | FDR = (200·6/101 + 100·3/101)/40 = **0.3713** |
| 5 | same, at a lower t: 100 targets, 20 and 15 decoys | FDR = **0.5743** (so the raw FDR is not monotone and `q` takes the running minimum from the top) |
| 6 | SW decoys 50 and 20, with 1 and 1 ≥ t of example 4 | FDR = (200·2/51 + 100·2/21)/40 = **0.4342** > 0.3713, so `which = "sandwich"` and `q` uses 0.4342 |
| 7 | one bin holds 10 targets and 0 decoys, 3 targets ≥ t, no other bin | FDR = 10/3 → `q` = **1.0** |
| 8 | `bin_index([-5, 0, 49, 50, 1000, nan], [-inf, 0, 50, 150, 300, inf])` | `[0, 1, 1, 2, 4, -1]` |
| 9 | `local_fdr` with decoy and target drawn from the same distribution (n = 2,000 each) | median lfdr within [0.8, 1.2]; no value > 1 |
| 10 | targets at \|Δ\| ≈ 0 and decoys spread out | lfdr near 0 for the smallest \|Δ\| and rising monotonically in \|Δ\| over the central 90 % |

## Acceptance criteria

- [ ] `tests/unit/stats/test_tdist.py` and `test_decoy_fdr.py` cover examples 1–10 and every
  degenerate row; examples 2–8 are hand-computed constants in the test, not golden files.
- [ ] Property: `q` is non-increasing in the score, lies in [0, 1], and equals the element-wise
  monotone curve of the element-wise maximum of the per-decoy FDRs, which is ≥ the maximum of the
  per-decoy monotone curves (corrected by lead decision 1 of M6-01-lead-decisions.md; both are tested
  so the order of operations is pinned).
- [ ] Property: multiplying every count in one bin by k leaves the FDR unchanged up to the +1 terms
  (checked as a limit at large counts).
- [ ] Property: `student_t_two_sided` is invariant to an affine shift applied to both `x` and `fit`.
- [ ] Calibration (fast, `tests/calibration/test_engine_calibration.py` style): with targets and
  decoys drawn from the same distribution and π₀ = 1, the mean FDP at `q ≤ 0.1` is ≤ 0.1 + 2 SE over
  R = 200 draws. (The full operon-level check is M6-09.)
- [ ] Performance: 50,000 targets × 2 decoy sets × 5 bins < 2 s.

## Hardening owed (pass 2)

- π̂₀ estimation instead of the fixed `operons.pi0 = 1.0`, if attack A in M6-09 shows the FDR is far
  too conservative; needs an ADR (ADR-0019 is the precedent).
- A profile-likelihood CI for the t parameters, so `operons.md` can report the fit's uncertainty.
- Replace Nelder–Mead with an EM for the t if the fit proves unstable on real data (M10 pilot).

## Methods note (statistician sign-off)

`docs/methods/operons.md` §"The within-operon distribution" and §"Decoy FDR for missed operons";
`docs/methods/target-decoy.md` gains §"Convergent and sandwich decoys" (why the larger of the two is
reported, and what exchangeability assumption the reweighting replaces).

## Questions for the lead

1. The spec writes `N_{NO,b}` in the numerator: the count of **all** NO pairs in bin b, not those
   above t. **Proposed default:** read it literally (step 4), which makes `FDR(t)` a ratio of an
   expected count of false targets to the observed count — the standard form.
2. `q` = running minimum from the top (standard) versus the spec's bare "make the FDR monotone".
   **Proposed default:** running minimum from the highest score down, which is the only direction
   that keeps `q` an upper bound on the FDP.
3. Element-wise maximum of two monotone curves, or monotone curve of the element-wise maximum?
   **Proposed default:** maximise first, then monotonise (step 5) — the other order can report a q
   below one of the two decoy FDRs.
