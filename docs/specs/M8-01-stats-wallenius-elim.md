# M8-01: Wallenius ORA with p_detect odds, and the generic GO-elim walker

Milestone M8 · track S (statistics engine) · **pass 1** · depends on: nothing (scipy ≥ 1.18 is
locked) · estimate ≈ 230 code lines (+ ≈ 350 test lines) · **statistical** (label `stats`)

## Spec sections implemented

- §7.7.3 "Setup. The universe is the eligible genes. Per category C: M = |universe|;
  n = |C ∩ universe|; N = number of detected genes; k = number detected in C; odds
  ω = mean(p_detect in C) / mean(p_detect outside C), the goseq weighting applied to
  detectability."
- §7.7.3 "p-values: over-representation: `scipy.stats.nchypergeom_wallenius(M, n, N, ω).sf(k − 1)`;
  under-representation: `.cdf(k)`. Wallenius is analytical, so there is no permutation floor."
- §7.7.3 "GO redundancy (elim, Alexa 2006). Test GO terms bottom-up; when a term is significant at
  the raw-p threshold 0.01, remove its genes from all its ancestors before testing those. Report
  `p_classic` and `p_value` (elim)."
- §8.3 rows `wallenius_ora(...)` (§7.7.3); §8.2 "Adjust for detectability"; §7.7.5 "Unit: elim on a
  toy DAG"; §7.7.5 golden "Wallenius p-values against goseq … relative difference ≤ 10⁻⁶" (checked
  in M8-13; this issue must make it achievable).

## Public interface

`src/giap/stats/wallenius.py` (pure)
```python
@dataclass(frozen=True, slots=True)
class OraResult:
    M: int; n: int; N: int; k: int
    odds: float            # ω; NaN when undefined
    p_over: float          # P(X ≥ k); NaN when undefined
    p_under: float         # P(X ≤ k); NaN when undefined

def category_odds(p_detect: ArrayLike, in_cat: ArrayLike) -> float: ...
def wallenius_tails(k: int, M: int, n: int, N: int, odds: float) -> tuple[float, float]: ...
def wallenius_ora(detected: ArrayLike, p_detect: ArrayLike, in_cat: ArrayLike) -> OraResult: ...
```
`src/giap/function/elim.py` (pure; generic over any DAG and any test)
```python
TestFn = Callable[[frozenset[int]], float]          # gene indices of a term -> raw p
def elim(term_genes: Mapping[str, frozenset[int]], parents: Mapping[str, frozenset[str]],
         test: TestFn, threshold: float = 0.01) -> dict[str, tuple[float, float]]: ...
                                                   # term -> (p_classic, p_elim)
```

1. **`category_odds`.** `mean(p[in]) / mean(p[~in])` over finite p. NaN when `in` is empty, when
   `~in` is empty (C = universe), or when `mean(p[~in]) == 0`. `inf` is never returned: ω = 0 is
   allowed (mean inside 0) and makes k > 0 impossible (p_over = 0 for k ≥ 1, exactly).
   Algebraically this is goseq's `avg_weight·(N_genes − n)/(Σpwf − n·avg_weight)` — the test in
   example 3 asserts the identity, so the golden comparison cannot drift on the odds.
2. **`wallenius_tails`** returns `(p_over, p_under)`. Support is `[lo, hi]` with
   `lo = max(0, N − (M − n))`, `hi = min(n, N)`. **Both tails are direct sums of the pmf** over the
   tail (`pmf` evaluated by `scipy.stats.nchypergeom_wallenius(M, n, N, ω)`), summed from the far
   end toward k for accuracy (decision 1). `p_over = Σ_{j=k}^{hi} pmf(j)`,
   `p_under = Σ_{j=lo}^{k} pmf(j)`, each clipped to [0, 1]. k < lo or k > hi is a `ValueError`
   (it cannot happen with consistent counts). ω = 1 must reproduce the central hypergeometric
   tails (scipy `hypergeom`) to 10⁻¹².
3. **`wallenius_ora`** counts M, n, N, k over the finite-p universe it is given (the caller passes
   the eligible genes only; NaN p_detect raises naming the indices, the M4-16 fail-closed rule) and
   returns an `OraResult`. It does **not** apply size filters or BH; M8-04 does.
4. **`elim`.** Terms are visited deepest first (level = longest path to a root; ties by term ID).
   For each term the test sees its genes minus every gene marked removed *for this term*; a term
   whose elim p is < `threshold` marks its current genes removed for **all** its ancestors
   (transitive closure of `parents`). `p_classic` = the test on the term's full gene set. A term
   whose remaining set is empty gets elim p = NaN (the caller decides the label). Cycles raise.
   The walker never calls the test twice with the same set for one term.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| n = 0 or C = universe | `odds` NaN, both p NaN (not tested; M8-04 writes `unestimated`) |
| N = 0 (nothing detected) | support {0}; k = 0; p_over = p_under = 1 |
| N = M (everything detected) | support {n}; k = n; both p = 1 |
| ω = 0 | p_over = 1 for k = 0, exactly 0 for k ≥ 1 (consistent: k ≥ 1 cannot occur) |
| a tail < 10⁻³⁰⁰ | returned as the float underflow (≥ 0), never negative, never NaN |
| k at the support bound | p_over at k = hi equals pmf(hi); p_under at k = lo equals pmf(lo) |
| elim: a term with no genes after removal | elim p NaN; ancestors unaffected by it |
| elim: two roots sharing a child | removal propagates to both roots |
| elim: threshold boundary p = 0.01 | **not** significant (strict <, as topGO's `elim` uses `< cutOff`) |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | M = 4000, n = 60, N = 1500, k = 55, ω = 1.3 | p_over = 3.2049528951 × 10⁻¹⁴ (direct sum; scipy 1.18.1 `.sf(54)` returns 0.0 here — recorded in the session) to rel. 10⁻⁹ |
| 2 | same, k = 40 | p_over = 7.356215112 × 10⁻⁴; `.sf(39)` agrees to rel. 10⁻⁸ (no cancellation at this size) |
| 3 | p_detect = [0.9, 0.9, 0.3, 0.3, 0.3, 0.3], C = first two | ω = 0.9 / 0.3 = **3.0**; goseq's form 0.9·4/(3.0 − 2·0.9) = **3.0** |
| 4 | ω = 1, M = 20, n = 5, N = 8, k = 4 | p_over = `hypergeom(20, 5, 8).sf(3)` to 10⁻¹² |
| 5 | C = universe | `odds` NaN, both p NaN |
| 6 | N = 0 | k = 0, p_over = p_under = 1 |
| 7 | k = n = N = 5, M = 4000, ω = 2 | support [0, 5]; p_over = pmf(5); p_under = **1**; and p_over(k) + p_under(k) − pmf(k) = 1 for every k of the support |
| 8 | toy DAG: root R ← A ← A1, A2; R ← B; genes A1 {1, 2, 3}, A2 {4}, A = A1 ∪ A2 ∪ {5}, B {6}, R = all; test = toy p table where A1 is significant | A1 visited first; its genes {1, 2, 3} removed from A and R; A tested on {4, 5}; R on {4, 5, 6}; `p_classic` of A uses {1…5} |
| 9 | same DAG, A1 p = 0.01 exactly | not significant; A tested on its full set |
| 10 | a cycle A → B → A | `ValueError` naming the cycle |

## Acceptance criteria

- [ ] `tests/unit/stats/test_wallenius.py` covers examples 1–7 and every degenerate row.
- [ ] `tests/unit/function/test_elim.py` covers examples 8–10 (the §7.7.5 "elim on a toy DAG").
- [ ] Property (Hypothesis): for random (M, n, N, ω) with M ≤ 500, `p_over(k) + p_under(k) −
  pmf(k) = 1 ± 10⁻¹²` for every k in the support; `p_over` is non-increasing in k and non-decreasing
  in ω; at ω = 1 the tails equal the central hypergeometric ones.
- [ ] Property: `elim` with `threshold = 0` returns elim p = classic p for every term.
- [ ] Performance (benchmark, `benchmarks/baseline.json`): 3,000 categories × `wallenius_tails` on
  M = 4,500 in < 5 s on 1 CPU.
- [ ] `docs/methods/function.md` §"Wallenius ORA" (statistician sign-off): the odds, the direct
  sums and why, the approximation (mean odds) and the attack it implies.

## Calibration or attack work owed

- The `ora` suite and attack 2 are M8-15's; this issue owes the methods-note paragraph that says
  the mean-odds Wallenius is exact only for constant p_detect within and outside C, and names the
  Poisson-binomial fallback.

## Hardening owed (pass 2)

- Golden: Wallenius vs goseq on the M8-12 inputs at relative 10⁻⁶ (M8-13). If scipy's pmf itself
  differs from BiasedUrn's by more than that (both implement Fog's algorithm, with their own
  precision constants), the fix is a precision argument or a log-space pmf — an ADR, never a
  tolerance change.
- Mutation testing on `wallenius.py` (monthly mutmut covers `stats/`).

## Questions for the lead

1. The spec writes `.sf(k − 1)` / `.cdf(k)`; scipy computes `sf` as 1 − cdf and underflows to 0.0
   for tails below ~10⁻¹⁶ (example 1). **Proposed default:** direct pmf sums with the same
   meaning; no ADR needed (an implementation detail of the same quantity), but the methods note
   says so.
2. elim's significance threshold: strict `< 0.01` (topGO) or `≤`? **Proposed default:** strict.
3. Where does elim live? **Proposed default:** `giap.function.elim` (it walks a DAG, not a stats
   primitive); M8-05's redundancy check reuses its ancestor closure.
