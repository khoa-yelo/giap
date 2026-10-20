# M5-02: Stats engine: Moran's I, join counts, gene-label and TU-block permutations, block length

Milestone M5 · track S (engine) · depends on: none (reuses merged `perm_p`, M4-03) · estimate ≈ 300
code lines (+ ≈ 350 test lines) · **statistical** (label `stats`) · **can start now** · no tools, no
Slurm

## Spec sections implemented

- §8.3 "`morans_i(x, lag, circular)`, `join_count(d, lag, circular)` | §7.3.5".
- §7.3.5(c) quoted: "Moran's I at lag k = 1 … 20. Pairs (i, i + k) along the gene index (circular where
  applicable) where both genes have values: I_k = ((1/W_k) Σ_{(i,j):|i−j|=k} z_i z_j) / ((1/n) Σ_i z_i²),
  z = x − x̄. W_k is the number of such pairs and n the number of genes with values. Join-count for
  detection. BB_k = the number of lag-k pairs of eligible genes where both are detected. Two nulls,
  B permutations each: gene label: shuffle values among the replicon's genes; asks whether there is any
  structure; TU block: shuffle whole TUs along the replicon, keeping internal order and orientation;
  asks whether there is structure beyond operons. Reporting. Two-sided p per lag and null (centred on
  the null mean) … 2.5–97.5 % envelopes for V10."
- §7.3.5(d) step 5: "block length ℓ = the smallest lag at which TU-level Moran's I is not significant at
  0.05 (gene-label null), bounded to [1, 20]" (the rule lives here so M5-06 and M5-08 share it;
  docs/methods/bootstrap.md "Block length").
- §7.3.6 unit "Moran's I on an alternating series (−1)"; property "TU-block permutation preserves TU
  sizes and internal order"; Appendix B.4 first bullet.

## Public interface

`src/giap/stats/autocorr.py` (pure; randomness only from the generator passed in)
```python
def max_lag(n: int, circular: bool) -> int: ...          # circular (n-1)//2, linear n-1, n < 2 → 0
def lag_pairs(n: int, lag: int, circular: bool) -> tuple[I64, I64]: ...  # (i, j = i + lag [mod n]); empty if lag ∉ [1, max_lag]
def morans_i(x: ArrayLike, lag: int, circular: bool) -> float: ...       # NaN marks "no value"
def join_count(d: ArrayLike, lag: int, circular: bool) -> float: ...     # d boolean over eligible genes
def morans_i_matrix(X: F64, lags: Sequence[int], circular: bool) -> F64: ...   # (B, n) → (B, L)
def join_count_matrix(D: B8, lags: Sequence[int], circular: bool) -> F64: ...  # (B, n) → (B, L)
def label_permutations(n: int, B: int, rng, movable: ArrayLike | None = None) -> I64: ...  # (B, n) orders
def tu_blocks(tu_labels: Sequence[str], circular: bool) -> I64: ...      # block id per position
def tu_block_permutations(blocks: I64, B: int, rng) -> I64: ...          # (B, n) orders
@dataclass(frozen=True, slots=True)
class LagTest:
    value: float; null_mean: float; null_lo: float; null_hi: float; p_value: float; n_valid: int
def lag_test(obs: float, null: ArrayLike) -> LagTest: ...
def choose_block_len(z: ArrayLike, circular: bool, B: int, rng, alpha: float = 0.05,
                     max_len: int = 20) -> int: ...
```
- **Lags.** Circular: pairs (i, (i + k) mod n) for i = 0 … n − 1, valid for 1 ≤ k ≤ ⌊(n − 1)/2⌋ (larger
  k would count each pair twice or pair a gene with itself). Linear: pairs (i, i + k), i + k < n,
  valid for 1 ≤ k ≤ n − 1.
- **Moran's I.** x̄ and n over entries with values; z_i = x_i − x̄; W_k = pairs with both values.
  **Does not exist (returns NaN):** lag invalid, n < 2, W_k = 0, or Σ z² ≤ 1e-24·max(1, max x²) (no
  variance).
- **Join count.** BB_k = Σ over valid lag-k pairs of d_i·d_j, returned as a float; invalid lag → NaN. It
  exists for every valid lag (0 is a value).
- **Matrices** equal the scalar functions row by row (permuted series: `X = x[orders]`).
- **Gene-label permutations.** Each row a uniformly random permutation; with `movable` (boolean), only
  the movable positions are permuted among themselves and the others stay fixed. M5-05 passes the genes
  with values, so the missingness pattern stays where it is (see M5-05 Q1).
- **TU blocks.** Block = maximal run of consecutive positions with the same TU label; on a circular
  replicon a run at the end and a run at the start with the same label are one block (an origin-crossing
  TU). A TU interrupted by another TU's gene forms several blocks (documented, not an error).
- **TU-block permutations.** Each row: the blocks in a uniformly random order, each block's positions in
  their original internal order (orientation kept), concatenated; a wrapped block is emitted as its
  end-part followed by its start-part.
- **lag_test.** NaN null entries are dropped (`n_valid` counts the rest); `p_value` =
  `perm_p(obs, null, "two-sided")` (centred on the pooled mean, +1 rule); `null_mean` = mean of the
  valid null; `null_lo`/`null_hi` = linear quantiles 0.025/0.975. **Does not exist:** obs NaN or
  `n_valid` = 0 → everything NaN except `value`; null constant and equal to obs (no variation, e.g. no
  detected gene) → `p_value` NaN (a p of 1 from a degenerate null carries no information and would
  inflate the BH family).
- **choose_block_len.** For lags k = 1 … min(`max_len`, `max_lag(n)`): `lag_test` of Moran's I of z under
  B gene-label permutations; ℓ = the first k whose p > `alpha` or is NaN. If every computed lag is
  significant: ℓ = min(`max_len`, `max_lag(n)` + 1, n). n < 3 (no lag) → ℓ = 1. Lags are evaluated in
  order with one permutation matrix shared by all lags (draws fixed by B and n only).

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | B.4: x = +1, −1, … (n = 10), circular | `morans_i(x, 1, True) = −1.0`; lag 2 → +1.0 |
| 2 | x = [1, 2, 3, 4] | lag 1 linear 0.333333 (W = 3); circular −0.2 (W = 4) |
| 3 | x = [1, NaN, 3, 1, 3] | lag 1 linear −1.0 (W = 2); lag 2 linear 0.0 (W = 2); lag 2 circular 0.333333 (W = 3) |
| 4 | x = [2, 2, 2, 2, 2] | NaN (no variance) |
| 5 | x = [1, NaN, NaN, 2], lag 1 linear | NaN (W = 0) |
| 6 | n = 6, lag 3, circular | `max_lag = 2`; `morans_i` and `join_count` NaN; `lag_pairs` empty |
| 7 | d = [1, 1, 0, 1, 1, 0] | BB lag 1 circular 2, linear 2; lag 2 circular 2, linear 1 |
| 8 | `tu_blocks([A, A, B, C, C, A])` | linear [0, 0, 1, 2, 2, 3]; circular [0, 0, 1, 2, 2, 0] |
| 9 | `lag_test(3, [0, 0, 0, 1])` | p = 0.2, null_mean 0.25, null_lo 0.0, null_hi 0.925, n_valid 4 |
| 10 | `lag_test(0, [0, 0, 0])` | p NaN (degenerate), null_mean 0 |
| 11 | `choose_block_len` on iid N(0, 1), n = 2,000, B = 199 | 1 in ≥ 90 % of 100 seeds |
| 12 | same on AR(1) ρ = 0.9, n = 2,000 | 20 (bounded) in ≥ 95 % of 100 seeds |
| 13 | constant z | 1 (lag 1 Moran NaN → not significant) |

## Edge cases

- Permutations of a series of length 0 or 1: identity rows; all statistics NaN (no lag exists).
- A replicon whose genes are one TU: one block; the TU-block null is the identity → `lag_test` degenerate
  → p NaN (M5-05 reports it).
- `join_count` on all-detected or none-detected genes: BB constant under both nulls → p NaN.
- Matrices for B = 1,000, n = 5,000, 20 lags fit in memory in chunks (implementer chooses the chunking;
  results independent of it).

## Acceptance criteria

- [ ] `tests/unit/stats/test_autocorr.py`: examples 1–10 and 13.
- [ ] Property (`tests/property/test_autocorr_props.py`): every TU-block permutation row is a permutation
  of 0 … n − 1, keeps each block contiguous (cyclically on circular input) in its internal order, and
  preserves the multiset of block sizes; label permutations with `movable` never move fixed positions;
  matrix functions equal scalar loops on random inputs; Moran's I ∈ [−1 − ε, 1 + ε] is **not** asserted
  (I_k can exceed 1 in small samples) — instead I_k is invariant to affine maps x → a·x + b (a ≠ 0).
- [ ] Calibration (`tests/calibration/test_engine_calibration.py`, suite `autocorrelation`, test-author):
  `lag_test` of `morans_i` lag 1 on iid normal series (n = 500) under label permutations: level at
  α ∈ {0.01, 0.05, 0.10}; examples 11–12 as slow tests (`calibration_nightly` if > 60 s).
- [ ] Performance: `morans_i_matrix` with B = 1,000, n = 4,500, 20 lags < 5 s on 1 CPU (benchmark
  recorded).

## Methods note (statistician sign-off)

`docs/methods/autocorrelation.md` §"Statistics and nulls": the two estimands (any structure vs structure
beyond operons), why the label null fixes the missingness pattern, block definition and wrap, degenerate
nulls → NA, two-sided centred p; `docs/methods/bootstrap.md` §"Block length" updated with the bounded
rule implemented here.

## Questions for the lead

1. "Not significant at 0.05" uses the raw per-lag p (not BH-adjusted) and two-sided p? **Proposed
   default:** raw two-sided p > 0.05 — the rule picks a correlation length, it is not a reported test.
2. If every lag up to 20 is significant, ℓ = 20 ("bounded to [1, 20]"). On replicons with fewer than 41
   TUs the lags stop earlier; ℓ = min(20, max_lag + 1, n). **Proposed default:** as stated.
3. A degenerate null gives `p_value` NA rather than 1. **Proposed default:** NA with the flag
   `not_estimable` (vocabulary addition proposed in M5-09).

## Lead decisions (2026-10-20, after the statistician's attack check)

4. Q1: the block-length rule uses the raw two-sided p > 0.05 as written; ℓ̂ enters only through
   ℓ' = min(20, 2ℓ̂ + 1, T) (ADR-0021, accepted).
5. Q2: ℓ̂ = min(20, max_lag + 1, n) as stated.
6. Q3: a degenerate null gives p NA plus `not_estimable`. `perm_p` returns 1/(B + 1) when the null is
   constant and the observation differs, so `lag_test`'s NA rule tests "null constant", not "null
   constant and equal to obs".
7. Residual join count (ADR-0020, accepted): `join_count(d, lag, circular, p=None)` and
   `join_count_matrix(D, lags, circular, P=None)` return Σ(d − p)(d − p) when p is given; p moves with d
   under both nulls; all p equal and all d equal → constant null → NA.
8. Example 12 stays; add a test that giapsim TU-score series at ρ 0.9 give ℓ̂ ≈ 8 (4,500 genes) and ≈ 1.3
   (1,000 genes), pinning the dilution.

## Lead decisions after the test-author's report (2026-10-20)
9. Example 12 amended: AR(1) ρ = 0.9, n = 2,000 gives ℓ̂ ≥ 15 in ≥ 95 % of runs (measured 100 %; ℓ̂ = 20
   in 0.85–0.92 because the sampling SD of I₁₇…I₂₀ ≈ ρ²⁰). The test-author updates the test.
10. `choose_block_len` returns ℓ̂ only. ADR-0021's ℓ' = min(20, 2ℓ̂ + 1, T) is a shared helper
    `giap.stats.autocorr.block_len_prime(l_hat, T)` (M5-06 and M5-08 call it); the test-author adds tests.
11. `lag_test` with NaN obs: `n_valid` is the count of finite null values; every other field NaN.
12. The < 5 s performance target goes to benchmarks/ (implementer adds the benchmark, no assertion in tests/).
13. Decision 11 wording: "finite" reads "non-NaN", matching `perm_p` (±inf null entries count as valid).

## Lead decisions after code review of #116 (2026-10-20)
14. "Constant null" means spread ≤ 1e-12 × max(1, max|null|) over non-NaN entries (the perm_p tie scale),
    not exact min == max; a circular replicon with 2 TU blocks (rotation-only null) gives p = NaN.
15. The property test asserts min ≤ null_mean ≤ max and null_lo ≤ null_hi (not lo ≤ mean ≤ hi).
16. Every documented ValueError has a test (autocorr.py checks: non-contiguous blocks, bad n/B/movable,
    l_hat/T < 1, max_len < 1, …); decision 13 (±inf counts in n_valid) is pinned; `choose_block_len`
    moves NaNs with values (documented and tested); B = 0 raises ValueError.
