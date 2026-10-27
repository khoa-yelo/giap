# M4-04: Stats engine IV: TU-block bootstrap, circular block bootstrap, Hodges–Lehmann shift

Milestone M4 · track S · depends on: none · estimate ≈ 260 code lines (+ ≈ 300 test lines) ·
**can start now** (circular block bootstrap is used by M5 islands; TU-block bootstrap by M5–M8)

## Spec sections implemented

- §8.3: "`tu_block_bootstrap(df, tu_col, stat_fn, B, rng)` | percentile CI"; "`circular_block_bootstrap(x,
  block_len, B, rng)` | for scan nulls".
- §8.2: "Resample TUs, not genes. Bootstraps and permutations move whole TUs." "Effect sizes with 95 %
  CIs accompany every p or q."
- §7.3.5(a): "O − E, with a 95 % CI from a TU-block bootstrap (resample the group's TUs with
  replacement)". (b): "Abundance statistic. The Hodges–Lehmann shift: the median of pairwise
  differences in log2 abundance between the region's detected genes and the replicon's other
  detected genes." "CI. TU-block bootstrap: resample the region's TUs and, for the HL shift, the
  outside TUs too."
- §7.3.6 unit list: "HL shift".

## Public interface

`src/gippa/stats/bootstrap.py`
```python
@dataclass(frozen=True, slots=True)
class BootCI:
    estimate: float; ci_low: float; ci_high: float; n_valid: int; B: int

def percentile_ci(boot: ArrayLike, level: float = 0.95) -> tuple[float, float]: ...
def tu_block_bootstrap(df: pd.DataFrame, tu_col: str, stat_fn: Callable[[pd.DataFrame], float],
                       B: int, rng: np.random.Generator, level: float = 0.95,
                       strata_col: str | None = None) -> BootCI: ...
def circular_block_bootstrap(x: ArrayLike, block_len: int, B: int,
                             rng: np.random.Generator) -> npt.NDArray[np.float64]: ...   # shape (B, n)
```
- `tu_block_bootstrap`: `estimate = stat_fn(df)`. Each replicate draws, within each stratum
  (`strata_col` value; one stratum when None), as many TU IDs as the stratum has, uniformly with
  replacement, and concatenates their rows (a TU drawn twice contributes its rows twice).
  NaN replicate statistics are dropped; if fewer than B/2 are valid, CI = (NaN, NaN).
  `percentile_ci`: `np.quantile(..., [(1−level)/2, (1+level)/2], method="linear")`.
- `circular_block_bootstrap` (Politis–Romano): ⌈n / block_len⌉ start positions uniform on [0, n);
  each block takes `block_len` consecutive values, wrapping modulo n; concatenated and truncated to n.
- Draw order is fixed (strata in sorted order, then replicates) so results are reproducible.

`src/gippa/stats/robust.py`
```python
def hodges_lehmann(x: ArrayLike, y: ArrayLike) -> float: ...   # median of {x_i − y_j}; NaN-free inputs
```
Exact: for m·n ≤ 10⁶ by direct enumeration; above, by bisection on Δ counting pairs with
x_i − y_j ≤ Δ on sorted y (numba, O(m log n) per step) until the two middle order statistics are
found exactly. Median of an even count = mean of the two middle values.

## Examples

| # | Call | Expected |
|---|---|---|
| 1 | `hodges_lehmann([10, 12], [8, 9, 11])` | differences [2, 1, −1, 4, 3, 1] → sorted [−1, 1, 1, 2, 3, 4] → 1.5 |
| 2 | `hodges_lehmann([5], [5])`; `hodges_lehmann([], [1])` | 0.0; NaN |
| 3 | `hodges_lehmann(x, y)` with x = 2,000 N(1,1) draws, y = 3,000 N(0,1) (seeded) | equals the direct `np.median(np.subtract.outer(x, y))` exactly (path above 10⁶ pairs) |
| 4 | `percentile_ci(np.arange(101))` | (2.5, 97.5) |
| 5 | `tu_block_bootstrap` on 3 TUs `{A: [1, 1], B: [0], C: [1]}` (d − p values), stat = sum, B = 2,000 | estimate 3.0; every replicate sum ∈ {0, 1, 2, …, 6} and is a sum of 3 TU totals from {2, 0, 1} (property check) |
| 6 | same frame, one TU only | CI degenerate at estimate (all replicates equal) |
| 7 | `tu_block_bootstrap(..., strata_col="inside")` with 2 inside and 5 outside TUs | each replicate has exactly 2 inside and 5 outside TU draws |
| 8 | `circular_block_bootstrap(np.arange(10), 3, B=1, rng)` with starts forced to 8, 0, 4, 9 (patched rng) | `[8, 9, 0, 0, 1, 2, 4, 5, 6, 9]` |
| 9 | `circular_block_bootstrap(x, 1, B, rng)` | equals iid resampling with replacement (each row's values ⊆ x) |

## Edge cases

- `block_len` < 1 or > n → `ValueError`; n = 0 → shape (B, 0).
- TU IDs that are NA → `ValueError` (the contract guarantees every gene has a TU).
- `stat_fn` raising on a resample propagates (no silent swallowing).
- `hodges_lehmann` with NaN in inputs → `ValueError` (callers drop missing abundances first).
- `level` must be in (0, 1).

## Acceptance criteria

- [ ] `tests/unit/stats/test_bootstrap.py`, `test_robust.py`: examples 1–9.
- [ ] Property: TU-block replicates contain only whole TUs (row multiset of each drawn TU intact, internal
  order preserved); circular block rows are runs of consecutive (mod n) values of length `block_len`
  except possibly the last; `hodges_lehmann(x + c, y) == hodges_lehmann(x, y) + c` and
  `hodges_lehmann(x, y) == −hodges_lehmann(y, x)`.
- [ ] Recorded (not gated): coverage of the TU-block percentile CI for a sum over 20 TUs of iid
  normal TU effects, R = 200 seeded replicates, written to the test log; flag if outside 0.95 ± 0.04.
- [ ] Benchmark: `hodges_lehmann` with m = 200, n = 4,500 under 2 ms (`benchmarks/test_stats.py`);
  M4-17 evaluates it once per shift.

## Methods note (statistician sign-off)

`docs/methods/bootstrap.md`: TU-block resampling (why TUs), stratified resampling for two-group
statistics, percentile CIs and their known under-coverage with few TUs, circular block bootstrap.
The HL shift and its exact computation go in `docs/methods/region-tests.md` (M4-17).

## Questions for the lead

1. Percentile CIs (as §8.3 says) rather than BCa? **Proposed default:** percentile, with the recorded
   coverage check above; BCa is a later improvement if coverage is poor.
2. Minimum valid replicates before reporting a CI: B/2? **Proposed default:** yes; otherwise CI = NA.

## Lead decisions after statistician review (PR #31, BLOCK)
1. `hodges_lehmann` raises `ValueError` on any non-finite input (NaN, ±inf) on both paths, and falls
   back to direct enumeration (or raises if that is too large) when the bisection bounds overflow
   (|differences| near the float64 limit); it never loops. Regression tests above `_DIRECT_MAX`.
2. The direct-enumeration threshold is m·n ≤ 5×10⁴ (`_DIRECT_MAX`), replacing 10⁶: the bisection path is
   bit-identical to the median of all pairwise differences (statistician verified) and ~6× faster.
3. Forwarded to M4-17: percentile CIs under-cover with few TUs (0.86 / 0.89 at 5 / 10 TUs); report
   `n_tu` and flag CIs from fewer than 10 TUs (docs/methods/bootstrap.md).
