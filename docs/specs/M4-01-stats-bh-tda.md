# M4-01: Stats engine I: Benjamini–Hochberg, target–decoy q-values, picked FDR

Milestone M4 · track S (stats engine) · depends on: none · estimate ≈ 180 code lines (+ ≈ 300 test
lines) · **can start now** · **schedule first:** the same (D + 1)/T monotone q is what M3-10 needs
for placement FDR (§7.2.5); M3-10 currently plans its own `stratum_qvalues` — see Q3

## Spec sections implemented

- §8.3 table rows: "`bh(p) -> q` | Benjamini–Hochberg, monotone, NaN-aware"; "`tda_qvalues(target_scores,
  decoy_scores, pi=1.0, plus_one=True)` | separate-search q-values, (D + 1)/T"; "`picked_fdr(target_best,
  decoy_best)` | picked competition, ties go to the decoy". "Each function is pure, typed and documented,
  with its methods note in `docs/methods/`."
- §7.3.2 picked FDR steps 3–4: "Per gene, keep max(s_g, s̃_g); a tie counts as a decoy." "FDR(t) =
  (#picked decoys ≥ t + 1) / #picked targets ≥ t, made monotone into gene q-values."
- §7.2.5 step 3: "q_i = min over e ≥ E_i of FDR_m(e)" (same monotonisation; M3 passes −E as score).
- §8.2: "BH applies within the families listed in §8.1 … All p and q values are reported."
- Appendix B.1 and B.6.

## Public interface

`src/giap/stats/multitest.py`
```python
def bh(p: ArrayLike) -> npt.NDArray[np.float64]: ...
```
`src/giap/stats/tda.py`
```python
def monotone_q(fdr_at_sorted: ArrayLike) -> npt.NDArray[np.float64]: ...
    # input: FDR at thresholds sorted from strictest to loosest; output: running min from the loose end, clipped to [0, 1]
def tda_qvalues(target_scores: ArrayLike, decoy_scores: ArrayLike, pi: float = 1.0,
                plus_one: bool = True) -> npt.NDArray[np.float64]: ...      # one q per target, input order
@dataclass(frozen=True, slots=True)
class PickedResult:
    decoy_wins: npt.NDArray[np.bool_]      # per gene; true also when only a decoy exists
    q: npt.NDArray[np.float64]             # per gene; NaN where the decoy wins or neither exists
    n_target_wins: int; n_decoy_wins: int
def picked_fdr(target_best: ArrayLike, decoy_best: ArrayLike) -> PickedResult: ...
```
Conventions: higher score = better; NaN = absent. Thresholds are the distinct observed scores;
at threshold t, counts include every score ≥ t (so a decoy tied with a target is counted with it).
`tda_qvalues`: FDR(t) = π·(D(t) + 1)/T(t) (or π·D(t)/T(t) when `plus_one=False`); q_i = min over
t ≤ s_i of FDR(t), clipped to 1. `bh`: NaN entries stay NaN and do not count toward m.

## Examples

| # | Call | Expected |
|---|---|---|
| 1 | `bh(B.1 p)` = `[0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216]` | `[0.01, 0.04, 0.084, 0.084, 0.084, 0.1, 0.105714, 0.216, 0.216, 0.216]` (abs tol 1e-6) |
| 2 | `bh([0.01, nan, 0.04, 0.03])` | `[0.03, nan, 0.04, 0.04]` (m = 3) |
| 3 | `bh([0.5, 0.5, 0.5])`; `bh([])` | `[0.5, 0.5, 0.5]`; empty array |
| 4 | B.6: `tda_qvalues(targets=150 scores 1000…851, decoys=[850])` | all 150 q = 1/150 = 0.0066667 (all accepted at 0.01) |
| 5 | B.6: `tda_qvalues(50 targets, decoys=[])` | all q = 1/50 = 0.02 (none accepted at 0.01) |
| 6 | `tda_qvalues([3, 2, 1], [2], plus_one=False)` | FDR(3)=0/1, FDR(2)=1/2, FDR(1)=1/3 → q = `[0, 1/3, 1/3]` |
| 7 | `tda_qvalues([3, 2, 1], [2], pi=0.5)` | FDR = 0.5·[1/1, 2/2, 2/3] → q = `[1/3, 1/3, 1/3]` |
| 8 | `picked_fdr(target=[5, 4, 3, 3, 2, 1], decoy=[1, nan, 3.6, nan, 2, nan])` (genes G1…G6) | picked list 5T, 4T, 3.6D, 3T, 2D (tie → decoy), 1T; FDR at 5, 4, 3.6, 3, 2, 1 = 1, 1/2, 1, 2/3, 1, 3/4; `q = [0.5, 0.5, nan, 2/3, nan, 0.75]`; `decoy_wins = [F, F, T, F, T, F]` |
| 9 | `picked_fdr([nan, 2], [1, nan])` | gene 0 only a decoy → `decoy_wins[0] = True`, `q[0] = nan`; gene 1 q = (1 + 1)/1 → 1.0 (clipped) |

## Edge cases

- `+inf` scores are allowed (e.g., a rank transform of q = 0); `-inf` is treated as absent.
- All-NaN input → all-NaN output; `n_target_wins = 0`.
- `pi` outside (0, 1] → `ValueError`; mismatched lengths in `picked_fdr` → `ValueError`.
- Ties among targets get equal q. Output order = input order (never sorted).
- `bh` with p > 1 or < 0 → `ValueError` (callers clip upstream; silent clipping hides bugs).

## Acceptance criteria

- [ ] `tests/unit/stats/test_multitest.py`, `test_tda.py`: examples 1–9.
- [ ] Property tests (`tests/property/test_stats_fdr_props.py`): `bh` q is monotone in p, ≥ p, ≤ 1,
  invariant to input permutation; `tda_qvalues` is non-increasing in score; adding a decoy never
  lowers any q; `picked_fdr` with `decoy_best` all NaN equals `tda_qvalues(target, [])` on the targets.
- [ ] `bh` agrees with `statsmodels.stats.multitest.multipletests(method="fdr_bh")` on 100 random
  vectors (rel tol 1e-12) — statsmodels is already a dependency.
- [ ] Calibration (in M4-05): BH FDR on independent nulls.
- [ ] mypy strict; functions pure (no RNG, no I/O).

## Methods note (statistician sign-off)

`docs/methods/multiple-testing.md` (BH, NaN handling, the §8.1 families as fixed in code) and
`docs/methods/target-decoy.md` (separate-search (D + 1)/T, picked competition, tie rule, monotone q).

## Questions for the lead

1. Tie at a threshold between a target and a decoy in `tda_qvalues`: count the decoy at that threshold?
   **Proposed default:** yes (≥ comparison on both), the conservative reading, matching "ties go to
   the decoy" in the picked rule.
2. Gene whose decoy wins: q = NaN (reported `NA`)? **Proposed default:** yes; the gene is
   `not_detected` (M4-12 decides the status).
3. M3-10 (`mapping/placement_fdr.py::stratum_qvalues`) re-implements separate-search q-values with a
   per-peptide π. **Proposed default:** M3-10 calls `tda_qvalues` per stratum (π̂ is constant within
   a stratum for database input; for de novo input with per-bin π̂, `tda_qvalues` gains an optional
   per-target `pi` array — a 10-line extension here), so one implementation is calibrated once.
