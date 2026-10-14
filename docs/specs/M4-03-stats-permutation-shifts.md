# M4-03: Stats engine III: permutation p-values, circular shifts, Besag–Clifford

Milestone M4 · track S · depends on: none · estimate ≈ 220 code lines (+ ≈ 300 test lines) ·
**can start now** (reused by M5 autocorrelation, M6, M8)

## Spec sections implemented

- §8.3: "`perm_p(obs, null, alternative)` | (1 + b)/(1 + B), two-sided centred on the null mean";
  "`besag_clifford(stat_fn, obs, h, max_B, rng)` | sequential Monte Carlo p"; "`circular_shifts(n, idx,
  circular)` | all valid shifts (§7.3.5)".
- §7.3.5(b): "The region's gene-index set I becomes (I + k) mod N, for k = 1 … N − 1 (exhaustive while
  N ≤ `exhaustive_shift_max_genes`, otherwise B random shifts)." "On linear contigs, only shifts that
  keep the set inside [0, N) are used." "Two-sided p = (1 + #{|S_k − S̄| ≥ |S_obs − S̄|}) / (1 + K),
  where S̄ is the null mean and K the number of shifts."
- §8.2: "p = (1 + b)/(1 + B), and seeds recorded"; "sequential Besag–Clifford permutation (stop after
  h = 10 exceedances), which avoid the 1/(B + 1) floor." §2.1 decision 14.
- §6.6: randomness only through `giap.rng.for_step`; functions here take a `np.random.Generator`.
- Appendix B.4: "Exhaustive circular shifts on a 100-gene circular replicon give K = 99 shifts. If the
  observed statistic is more extreme than all 99, p = (1 + 0)/(1 + 99) = 0.01."

## Public interface

`src/giap/stats/permutation.py`
```python
Alternative = Literal["two-sided", "greater", "less"]

def perm_p(obs: float, null: ArrayLike, alternative: Alternative = "two-sided") -> float: ...
def circular_shifts(n: int, idx: ArrayLike, circular: bool) -> npt.NDArray[np.int64]: ...
def apply_shift(idx: ArrayLike, k: int, n: int) -> npt.NDArray[np.int64]: ...      # (idx + k) mod n
def choose_shifts(valid: ArrayLike, n_max_exhaustive: int, n: int, B: int,
                  rng: np.random.Generator) -> tuple[npt.NDArray[np.int64], bool]: ...  # (shifts, exhaustive)
def besag_clifford(stat_fn: Callable[[np.random.Generator], float], obs: float, h: int,
                   max_B: int, rng: np.random.Generator) -> tuple[float, int]: ...   # (p, B_used); "greater"
```
- `perm_p`: NaN entries of `null` are dropped (B = finite count); `obs` NaN or B = 0 → NaN.
  `greater`: (1 + #{null ≥ obs})/(1 + B); `less`: (1 + #{null ≤ obs})/(1 + B); `two-sided`: centred
  on the null mean S̄ as quoted. Comparisons use a tolerance `1e-12 · max(1, |x|)` so float noise in
  exhaustive shifts (identical sums computed in different orders) counts as a tie.
- `circular_shifts`: circular → k = 1 … n − 1; linear → the k in 1 … n − 1 for which (idx + k) mod n
  does not wrap, i.e. `max(idx) + k < n` or `min(idx) + k ≥ n` (the second set are the leftward
  translations). Sorted ascending.
- `choose_shifts`: exhaustive when `n ≤ n_max_exhaustive` (config `context.exhaustive_shift_max_genes`);
  else B shifts drawn **without replacement** from `valid` (all of them when B ≥ len(valid)), sorted.
- `besag_clifford`: draws null statistics one at a time; stops at the h-th null ≥ obs and returns
  p = h / L (L = draws so far); if `max_B` is reached with b < h exceedances, p = (b + 1)/(max_B + 1).
  Two-sided use: the caller passes centred absolute statistics.

## Examples

| # | Call | Expected |
|---|---|---|
| 1 | `perm_p(5.5, [1, 2, 3, 4, 5])` | S̄ = 3, |dev| null = [2, 1, 0, 1, 2], none ≥ 2.5 → 1/6 = 0.166667 |
| 2 | `perm_p(1.0, [1, 2, 3, 4, 5])` | |1 − 3| = 2; two nulls ≥ 2 → 3/6 = 0.5 |
| 3 | `perm_p(4, [1, 2, 3, 4, 5], "greater")`; `"less"` | (1 + 2)/6 = 0.5; (1 + 4)/6 = 0.833333 |
| 4 | B.4: 99 null values all strictly less extreme than obs | 0.01 |
| 5 | `circular_shifts(10, [2, 3, 4], circular=True)` | `[1, …, 9]` (K = 9) |
| 6 | `circular_shifts(10, [2, 3, 4], circular=False)` | `[1, 2, 3, 4, 5, 8, 9]` (K = 7; 8 and 9 shift left by 2 and 1) |
| 7 | `circular_shifts(5, [0, 4], circular=False)` | `[]` (no non-wrapping translation) |
| 8 | `apply_shift([98, 99], 3, 100)` | `[1, 2]`; region size preserved |
| 9 | `besag_clifford(fn, obs=5, h=2, max_B=100, rng)` with `fn` yielding 1, 2, 6, 0, 0, 0, 7, … | stops at L = 7 → p = 2/7, B_used = 7 |
| 10 | same, `fn` never ≥ 5, `max_B = 100` | p = 1/101, B_used = 100 |

## Edge cases

- `circular_shifts(1, [0], …)` → empty. Empty `idx` or any index outside [0, n) → `ValueError`.
- Non-contiguous `idx` (AMR gene sets) is allowed; shifts translate the whole set.
- Duplicate indices in `idx` → `ValueError`.
- `choose_shifts` with `len(valid) == 0` → empty array; callers report `p_value = NA`.
- `besag_clifford` with `h < 1` or `max_B < h` → `ValueError`.
- Determinism: the same `rng` state gives the same shifts and the same B_used.

## Acceptance criteria

- [ ] `tests/unit/stats/test_permutation.py`: examples 1–10.
- [ ] Property tests: shifted sets keep their size and their pairwise gaps (circular distance on
  circular replicons; plain difference on linear); linear shifts never wrap; `perm_p` ∈ [1/(B+1), 1];
  `perm_p` two-sided is invariant to adding a constant to obs and null.
- [ ] Calibration (in M4-05): `perm_p` and `besag_clifford` at level on exchangeable nulls.
- [ ] No global RNG: grep test that `np.random.` appears only as a type in `stats/`.

## Methods note (statistician sign-off)

`docs/methods/permutation-and-shifts.md`: the +1 rule, centring two-sided p on the null mean, the
circular-shift null (what it preserves: region size, spacing, autocorrelation; what it breaks: the
region's position), linear-contig restriction, exhaustive vs random shifts, Besag–Clifford stopping.

## Questions for the lead

1. Random shifts when N > 20,000: with or without replacement? **Proposed default:** without (no
   duplicated null draws; for N > 20,000 genes and B = 1,000 the difference is negligible).
2. Float-tie tolerance `1e-12 · max(1, |x|)` in `perm_p`: acceptable? **Proposed default:** yes; it
   only ever makes p larger (conservative).
