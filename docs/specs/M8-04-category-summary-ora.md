# M8-04: `category_summary` and detection enrichment (conditional Poisson-binomial ORA on TUs, BH families, elim, `parent_group`)

Milestone M8 · track R · **pass 1** · depends on: M8-01 (Wallenius for `p_classic`, the elim
walker; branch `m8/wallenius`, not yet on main), M8-03 (collections, merged) · estimate
≈ 390 code lines (+ ≈ 560 test lines) · **statistical** (label `stats`)

**Revision 2 (2026-10-21) follows ADR-0029** (accepted by Khoa): the ORA p-value comes from the
**conditional Poisson-binomial null counted in TUs**. Wallenius (M8-01) remains only as
`odds_weight` / `p_classic` and for golden parity with goseq (M8-13). Until M8-15's arms pass, ORA
rows carry no calibrated claim (ADR-0029, M8-01 lead decision 5). `category_summary` is unchanged
by the ADR.

## Spec sections implemented

- §7.7.2 "One row per COG category, KEGG pathway and GO slim term, with: n_genes, n_detected;
  expected (Σ p_detect), O − E with a TU-bootstrap CI; median log2 abundance of detected genes;
  coverage. Mass fraction is deferred (improvement 11)."
- §7.7.3 as amended by **ADR-0029**: the setup (universe = eligible genes), "Filters. Category sizes
  are limited to `function.ora_size` ([5, 500]). BH runs within each ontology and direction.",
  "GO redundancy (elim …) Report `p_classic` and `p_value` (elim); BH applies to the elim p-values.
  Group surviving terms under their parents (`parent_group`)." The §7.7.3 Wallenius p-value
  formula is kept only for `p_classic`.
- ADR-0029 "The ORA p-value is computed under the conditional Poisson-binomial null. The unit is
  the TU: a TU is detected when any of its genes is, with probability 1 − ∏(1 − p_detect)."
- §8.1 row "Category detection | category | … | BH per ontology and direction | 0.05" (null
  replaced by ADR-0029); §8.2 "Resample TUs, not genes", "Adjust for detectability", "Effect sizes
  with 95 % CIs accompany every p or q"; §8.3 `poisson_binomial_sf` (extended here, see 1); §9.2
  test columns; §9.3 `category_summary.tsv` and `enrichment_detection.tsv`.

## Public interface

`src/giap/stats/poisson_binomial.py` (pure; extends the merged module)
```python
@dataclass(frozen=True, slots=True)
class CpbTails:
    k: int; N: int
    p_over: float          # P(K >= k | Σ = N)
    p_under: float         # P(K <= k | Σ = N)
    mean: float            # E[K | Σ = N]
    lo: int; hi: int       # support of K given N

def poisson_binomial_pmf_trunc(p: ArrayLike, kmax: int) -> np.ndarray: ...   # P(X = j), j = 0..kmax
def tilt(p: ArrayLike, total: int) -> np.ndarray: ...      # odds × e^θ with Σ p(θ) = total
def conditional_poisson_binomial(p_in: ArrayLike, p_out: ArrayLike, k: int, N: int) -> CpbTails: ...
```
`src/giap/function/ora.py` and `src/giap/function/summary.py` (pure)
```python
@dataclass(frozen=True, slots=True)
class TuUniverse:
    tu_ids: tuple[str, ...]                 # TUs with >= 1 eligible gene, natural-sorted
    q: F64                                  # 1 − ∏(1 − p_detect) over the TU's eligible genes
    detected: B8                            # any eligible gene of the TU detected
    gene_tu: Mapping[str, int]              # eligible locus_tag -> TU index

def tu_universe(genes: pd.DataFrame) -> TuUniverse: ...
def tu_members(u: TuUniverse, genes: frozenset[str]) -> frozenset[int]: ...
def category_summary(genes: pd.DataFrame, colls: Sequence[GeneSets], B: int,
                     rng: np.random.Generator) -> pd.DataFrame: ...
def detection_enrichment(genes: pd.DataFrame, colls: Sequence[GeneSets],
                         size: tuple[int, int], alpha: float) -> pd.DataFrame: ...
```

### 1. The conditional null (`giap.stats.poisson_binomial`)

The TUs are independent units. TU t is detected (D_t = 1) with probability q_t. For a category C
(a set of TUs) with complement Ū in the TU universe U, let K = Σ_{t∈C} D_t and S = Σ_{t∈U} D_t.
The ORA null is the law of K given S = N:

  P(K = j | S = N) = PB_C(j) · PB_Ū(N − j) / Σ_i PB_C(i) · PB_Ū(N − i),

where PB_A is the Poisson-binomial pmf of the q_t in A (the "conditional Bernoulli" law). Then

- `p_over` = Σ_{j ≥ k} P(K = j | N), `p_under` = Σ_{j ≤ k} P(K = j | N). Both include the
  observed k, so `p_over + p_under − P(K = k | N) = 1`. They are sums of non-negative terms, never
  `1 − cdf`.
- `mean` = Σ_j j · P(K = j | N), used for the effect size (4).
- The support is [lo, hi] = [max(0, N − |Ū|), min(|C|, N)], further narrowed to the j with a
  non-zero weight. A k outside it is a `ValueError` (it cannot happen with consistent counts).

**Computation: exact DP, no approximation, no Monte Carlo.**

1. **Truncated DP.** `poisson_binomial_pmf_trunc(p, kmax)` is the merged `_pmf` recurrence,
   stopped at index kmax. PB_C needs kmax = min(|C|, N), and PB_Ū needs kmax = N. The convex-
   combination recurrence keeps relative precision across the whole support: checked in this
   session against exact rationals (`fractions.Fraction`) on 120 Beta(1, 4) probabilities, the
   maximum relative error was **7.1 × 10⁻¹⁶** over the full support, down to pmf = 1.2 × 10⁻¹⁰⁹.
   The inner loop is numba-compiled (`@njit`, a sequential loop, so the §6.6 reduction-order rule
   is met trivially). The merged pure-numpy `_pmf` costs 15 ms at T = 2,500. The numba truncated
   loop (spec-writer check, `$GIAP_ROOT/stat/m8-04-spec/bench.py`) costs ≈ 0.4 ms per category at
   T = 2,500 TUs and N = 1,250, for any |C| in {5, 50, 500}.
2. **Tilting against underflow.** P(K | S) is invariant when every odds q/(1 − q) is multiplied by
   the same e^θ (example 5). `tilt(q, N)` solves Σ_t expit(logit q_t + θ) = N for θ (monotone in θ;
   Newton with a bisection guard to |Σ − N| ≤ 10⁻⁹ N). It is computed **once per run**, because U
   and N do not depend on the category. q_t ∈ {0, 1} stay fixed. All DPs use the tilted q. The mass
   of PB_Ū is then centred near the indices that are read, so no weight underflows even when p̂ is
   badly scaled. That would be possible with raw q, for example q ≈ 10⁻⁶ everywhere and N large.
3. **Normalise by the support sum**, not by a separate PB_U(N). This is the same quantity
   mathematically, it saves a third DP, and it makes the tails sum to 1 by construction. If the sum
   is 0 or not finite, the primitive raises `ValueError`. That happens only when the observation
   has probability 0 under q, for example a detected TU with q = 0.
4. **Not used, and why.**
   - Deconvolution of PB_U by C's factors: unstable for q > 0.5 and over ~500 removals.
   - FFT or characteristic-function inversion: absolute error ~10⁻¹⁶ × the mode, so relative
     precision in the tails is lost.
   - The saddlepoint: an approximation, unnecessary at ≈ 0.4 ms per exact call.
   - Monte Carlo conditional-Bernoulli sampling: the exact DP is cheap. Exhaustive enumeration
     over 2^T configurations (T ≤ 14) is the **test oracle**.
5. **Precision contract.** Every returned tail ≥ 10⁻³⁰⁰ is within relative **10⁻¹²** of the exact
   rational value. Smaller tails are returned as the float underflow (≥ 0, never negative, never
   NaN). A tail of 0.0 is reachable only below ~10⁻³⁰⁰. This is far below the ~10⁻¹⁴ floor of
   Wallenius (M8-01 methods note).

**Reuse.** `poisson_binomial_pmf`, `poisson_binomial_cdf` and `poisson_binomial_sf` stay as they
are (`mapping/chance.py` uses them). The truncated pmf replaces the private `_pmf` loop with the
numba kernel. The existing pmf tests then pin the refactor: same values to 10⁻¹⁵. The conditional
version is **new**. No merged helper conditions on the total, so `conditional_poisson_binomial`
has to be added. The §8.3 engine row `poisson_binomial_sf` is extended, not replaced.

### 2. TU universe and membership (`tu_universe`, `tu_members`)

- **Units.** A TU is in the universe when ≥ 1 of its genes is eligible (§7.3). Only eligible genes
  count, for both q and detection. q_t = 1 − ∏_{i∈t, eligible}(1 − p_detect,i). d_t = any eligible
  gene of t has `detected = true`.
- **Membership rule (consistent with M8-05).** TU t is in set S when **any** of its eligible genes
  is in S. `GeneSets` already hold eligible genes only (M8-03). M8-05 applies the same rule to the
  *detected* genes, because only ranked TUs exist there. ORA must not let membership depend on the
  detection outcome, since that would invalidate the null. So ORA applies the rule to **all**
  eligible genes (Q1).
- **Fail closed** (M4-16 rule, as in `context/summary.py`): an eligible gene with NaN `p_detect` or
  without a `tu_id` raises `ValueError` naming the locus tags.
- TU order is the natural sort of `TU_<n>`, so results do not depend on the input row order.

### 3. `category_summary` (unchanged by ADR-0029)

Collections `cog`, `kegg`, `go_slim` only (§9.3 enum). Per category:
- `n_genes` (eligible members) and `n_detected`;
- `expected` = Σ p_detect;
- `o_minus_e` with a 95 % TU-block bootstrap CI (`giap.stats.tu_block_bootstrap` over the
  category's TUs, B = `stats.bootstrap`, the M4-16 `basic_rows` construction: per-TU sums of
  d − p resampled);
- `median_abundance_log2` over detected members with a non-NA abundance (imputed values
  included, as everywhere under the primary policy);
- the collection's `coverage`.

Every category is reported, whatever its size: this is a descriptive table with no test. The
table stays **gene-level**, because §9.3 names `n_genes`, and its CI is already TU-aware. Rows are
sorted by collection (cog, kegg, go_slim), then category.

### 4. `detection_enrichment`

All six collections. The TU universe U, N = Σ d_t, and the tilted q are computed once. Per term
and direction (`over`, `under`):

| Column | Value |
|---|---|
| `collection`, `term`, `name`, `direction` | as M8-03; `direction` ∈ {over, under} |
| `n_universe` | \|U\| (TUs) |
| `n_in_set` | n = \|C\| in TUs (rule of 2) |
| `n_detected` | N = detected TUs in U |
| `k_detected_in_set` | k = detected TUs in C |
| `n_in_set_genes`, `k_detected_in_set_genes` | **added**: the gene-level n and k that `p_classic` uses (goseq parity, M8-13 hardening; readers) |
| `odds_weight` | ω = M8-01 `category_odds` on genes (mean p_detect in C / outside) |
| `p_classic` | M8-01 Wallenius tail for this direction on genes, full set (never elim'd, never in BH) |
| `p_noelim` | **added**: the conditional-PB tail on C's full TU set |
| `p_value` | the elim conditional-PB p (GO namespaces, `over`); otherwise `p_noelim` |
| `q_value` | BH (5) on `p_value` |
| `stat` | k |
| `effect` | k − E[K \| S = N] for C's **tested** TU set (full set, or the elim-reduced one), in TUs |
| `ci_low`, `ci_high` | NA in pass 1 (accepted Q6 default); flag `ci_unestimated` |
| `parent_group` | see 7 |
| `coverage` | the collection's (gene-level, §7.7.1) |
| `method` | `cpb_tu_ora` |
| `null_model` | `conditional_poisson_binomial_tu` |
| `n_perm` | 0 |
| `flags` | index decision 12 (ADR-A, owed by M8-06): `unestimated`, `null_degenerate`, `odds_out_of_range`, `ci_unestimated`, `elim_empty` |

- **Count columns are in TUs**, because the unit tested is the TU. The units are written in the
  schema description, since §9.2 puts units in names only for physical quantities (Q3).
- **Wallenius is computed on genes.** It is the goseq quantity, and M8-13 compares it against
  goseq. If it fails closed (M8-01 decision 4, ω outside `ODDS_RANGE`), then `p_classic` is NA
  and the flag is `odds_out_of_range`. `p_value` is not affected.

### 5. Size filter and BH

- **Size filter.** Terms with n (TUs) outside `function.ora_size` ([5, 500]) are not tested, and
  no rows are written. This matches GSEA, where sizes are counted in TUs (§7.7.4). Under elim the
  filter applies to the term's **full** TU set (accepted Q3 default). A term that elim reduces
  below 5 TUs is still tested on what remains.
- **BH families** = (collection, direction) (§7.7.3 "within each ontology and direction"; §8.1),
  on `p_value`, with `giap.stats.bh`. A NaN p stays out of BH and carries `flags = unestimated`.
  `p_classic` and `p_noelim` are reported only and are never adjusted.

### 6. elim (`go_bp`, `go_mf`, `go_cc`, direction `over` only, as accepted)

- **Inputs.** M8-01's walker, unchanged, is fed **TU index sets**: `term_genes[term] =
  tu_members(u, set)`. The test is `TestFn(tus) -> conditional_poisson_binomial(q[tus], q[U ∖
  tus], k(tus), N).p_over`, with U, N and the tilt fixed. The threshold is strictly `< 0.01` on
  that p.
- **Removal is of TUs.** A significant child removes its TUs from every ancestor. Removing genes
  instead would leak: a TU of the child that also holds another gene of the ancestor would stay in
  the ancestor and count the same detection event again (example 7; Q2).
- **Outputs.** `p_value` is the elim p, `p_noelim` is the walker's `p_classic` (the full-set CPB
  p), and the Wallenius `p_classic` column is separate.
- **No elim for `under`, `cog`, `kegg` or `go_slim`.** There, `p_value = p_noelim`.
- **R4 label** (M8-01 lead decision 5): "BH on elim p; conditional on removing elim-significant
  descendants; FDR control is empirical". It goes in the methods note and the report legend.

### 7. `parent_group` (GO namespaces only; accepted index decision 11)

For a row with `q_value ≤ stats.alpha.enrichment`, `parent_group` is the most general
**significant** ancestor in the same direction and namespace (by depth, ties by ID), or the term's
own ID if it has none. It is NA otherwise and for non-GO collections.

Rows are sorted by collection, term, then direction (M8-06).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| C covers every TU of U (Ū empty) | not a test: `p_value`, `p_noelim` NA, `flags = unestimated`, out of BH (ω is NaN as well) |
| N = 0 (no detected TU) | support {0}; every p = 1; rows written; nothing significant |
| N = \|U\| (every TU detected) | support {n}; every p = 1 |
| a TU with q = 1 is undetected, or q = 0 is detected | the observation has probability 0: `ValueError` from the primitive. Every row of the run gets `p_value` NA and `flags = unestimated;null_degenerate`, plus a run-level warning. Unreachable from M4-13, whose logistic p lies strictly in (0, 1) |
| all q equal (every TU the same size and constant p) | the null is the central hypergeometric on TUs (example 4) |
| `detectability_constant` (p_detect = the overall rate) | ω = 1 for `p_classic`, but q_t still grows with TU size, so the TU null is **not** hypergeometric. The run flag from M4-13 stays in the manifest |
| ω undefined or out of range (Wallenius only) | `odds_weight` / `p_classic` NA, flag `odds_out_of_range` or `unestimated` on `p_classic` only; the CPB columns stand |
| elim reduces a term to 0 TUs | `p_value` NA, `flags = unestimated;elim_empty`; `p_noelim` and `p_classic` kept |
| a TU partly in C (one member gene, others not) | the whole TU is one unit of C and may be detected through a non-member gene (example 1). This is a power loss, not a level loss; it is measured in M8-15 arm 7 |
| an empty collection | no rows; `category_summary` still reports coverage for cog, kegg and go_slim |
| category of < 10 TUs in `category_summary` | O − E CI NA, flag `ci_few_units` (M4-16) |
| sim100k | most terms are below 5 TUs, so `enrichment_detection.tsv` may be nearly empty. That is a **passing** pass-1 result |

## Examples (numbers computed in the session; `$GIAP_ROOT/stat/m8-04-spec/ex.py`, `ex2.py`)

**Toy universe T6.** Ten eligible genes in six TUs:

| TU | Genes (p_detect) | q |
|---|---|---|
| TU_1 | g1 (0.9), g2 (0.5) | 0.95 |
| TU_2 | g3 (0.5) | 0.5 |
| TU_3 | g4 (0.2), g5 (0.2) | 0.36 |
| TU_4 | g6 (0.5) | 0.5 |
| TU_5 | g7, g8, g9 (0.3 each) | 0.657 |
| TU_6 | g10 (0.1) | 0.1 |

Detected genes are g2, g4, g6 and g9, so the detected TUs are TU_1, TU_3, TU_4 and TU_5.

| # | Input | Expected |
|---|---|---|
| 1 | T6, C = {g1, g4, g7}, so C's TUs = {TU_1, TU_3, TU_5} | n = 3, k = 3, N = 4, \|U\| = 6. `p_over` = P(K = 3 \| N = 4) = 2134593/5159698 = **0.41370503** (enumeration oracle); `p_under` = **1**; E[K \| N] = **2.39086047**; `effect` = **0.60913953**. Gene level: n = 3, k = 1 (g4), N = 4, M = 10; gene CPB p_over 0.97289, p_under 0.56849 (recorded only). ω = 0.46667/0.34286 = **1.36111**; `p_classic` (over) = Wallenius(10, 3, 4, ω) tail k ≥ 1 = **0.89763087**, (under) = **0.54832563**. TU_1 counts as detected through g2, which is not a member (partial-membership row) |
| 2 | **TU-packing.** One TU of 5 genes (p 0.5 each, q = 31/32) = C, plus 20 singleton TUs (p 0.5). All 5 C genes and 10 singletons detected | TU level: n = 1, k = 1, N = 11; `p_over` = **341/351 = 0.97150997**. Gene level: hypergeometric(25, 5, 15), P(k ≥ 5) = **13/230 = 0.05652174**. The gene count treats one detection event as five; this is the teeth of arm 5 |
| 3 | example 1 with `ora_size` [4, 500] | no row (n = 3 TUs) |
| 4 | 20 singleton TUs, q = 0.4 each; C = 5 of them; N = 8, k = 4 | `p_over` = `hypergeom(20, 5, 8).sf(3)` = **0.05779154**, `p_under` = `.cdf(4)` = **0.99638803**, both to 10⁻¹² |
| 5 | q_in = (0.9, 0.5, 0.36), q_out = (0.5, 0.657, 0.1, 0.3, 0.8), k = 2, N = 4; then every odds × 3 | `p_over` = **0.63641697**, `p_under` = **0.91762557**, unchanged after scaling to 10⁻¹² (tilt invariance) |
| 6 | the same term, direction `under` | `p_value = p_noelim` (CPB under tail); separate BH family (collection, under) |
| 7 | **elim on TUs.** M8-01's toy DAG R ← A ← A1, A2; R ← B. TUs: a = {1, 2}, b = {3, 4}, c = {5}, d = {6}. Genes: A1 {1, 2, 3}, A2 {4}, A = {1…5}, B {6}. So A1 → {a, b}, A2 → {b}, A → {a, b, c}, B → {d}. A toy p table makes A1 significant | A1 is visited first and removes TUs {a, b} from A and R. A2 is tested on {b} (not an ancestor of A1, so not reduced). A is tested on {c}, R on {c, d}. Gene removal would have left A = {4, 5} → TUs {b, c}, counting TU b again (rejected design). A's `p_noelim` uses {a, b, c} |
| 8 | A1 and A both q ≤ 0.05 (over, go_bp) | `parent_group(A1)` = **A** (R not significant); `parent_group(A)` = **A** |
| 9 | C = every TU of U | `p_value` NA, `flags unestimated`; the other rows' q equal BH without it |
| 10 | the same input twice with the same seed | identical tables (the bootstrap draws from `function.bootstrap`; ORA and elim are deterministic) |

## Acceptance criteria

- [ ] `tests/unit/stats/test_conditional_poisson_binomial.py`:
  - examples 1, 2, 4 and 5;
  - Hypothesis, T ≤ 14: `p_over`, `p_under` and `mean` equal exhaustive enumeration to 10⁻¹²
    relative;
  - T ≤ 60: equal to a `Fraction` DP to 10⁻¹² relative for tails ≥ 10⁻³⁰⁰;
  - `p_over + p_under − P(K = k) = 1 ± 10⁻¹²` over the support;
  - `p_over` non-increasing in k;
  - odds-scaling invariance;
  - q ≡ 10⁻⁶ with N = T/2 returns finite tails (the tilt);
  - the `ValueError` row;
  - the merged `poisson_binomial_*` tests still pass unchanged.
- [ ] `tests/unit/function/test_ora.py` and `test_summary.py` cover examples 1–10 and every
  degenerate row. The `category_summary` cases from revision 1 are kept: n_genes 5, n_detected 3,
  expected 3.3, o_minus_e −0.3 on the 10-gene universe, CI NA with `ci_few_units`.
- [ ] Property: TU membership never depends on `detected`. Flipping any `detected` value leaves
  every `n_in_set` unchanged.
- [ ] Property: relabelling genes outside C, or permuting TUs outside C that have equal q, leaves
  C's row unchanged. Splitting a C-free TU into singletons with the same product changes only
  `n_universe`.
- [ ] Property: BH families are exactly (collection, direction); no NaN p enters one; `p_classic`
  and `p_noelim` never enter BH.
- [ ] Property: a GO term with no elim-significant descendant has `p_value == p_noelim`. No
  inequality is asserted otherwise, because elim can move p either way.
- [ ] Benchmark (`benchmarks/baseline.json`, implementer): 3,000 categories × both tails plus one
  elim re-test each, at T = 2,500 TUs and N = 1,250, in **< 10 s on 1 CPU**.
- [ ] `docs/methods/function.md` §"Detection enrichment" (statistician sign-off) covers:
  - the TU unit and the q formula with its independence assumption (Q4);
  - the conditional null, the tilt and the precision contract;
  - why Wallenius stays only as `p_classic`;
  - the families, elim on TUs and its R4 label;
  - the effect size;
  - partial-membership dilution.
- [ ] `ora.py` writes no claim: no column or report text says "calibrated" until M8-15 passes
  (ADR-0029).

## Calibration and attack arms (run in M8-15; this issue makes each runnable on a gene frame plus `GeneSets`)

The level bound is §8.4: α + 2√(α(1 − α)/R_eff). It is checked with `check_level_pooled`, per
direction, at α ∈ {0.01, 0.05, 0.10}. p̂ is always the **fitted** M4-13 cross-fitted p_detect
(backlog rule); the truth arm is recorded beside it.

| # | Arm | Generator (M8-14) | Method | Criterion |
|---|---|---|---|---|
| 1 | `ora` (**gate**) | `ora_null(4000)`: gene detection ~ Bernoulli(p_true), categories independent | TU-CPB | level |
| 2 | attack 2, weighted (**gated**) | `ora_attack_confounded(4000, −1.5)` with **TU-aligned** confounded categories (members fill their TUs) | TU-CPB | level, over and under |
| 3 | attack 2, teeth (**gated**) | the same | TU hypergeometric (all q equal) | under rejection at α = 0.05 > 2α |
| 4 | TU-packed attack (**gated**, new) | `ora_attack_tu_packed` (M8-14 to add): logit p = covariates + c · u_TU, u_TU ~ N(0, 1) shared by the TU (abundance), c = 2.5; p̂ from covariates only; categories = 15 random whole TUs | TU-CPB | level |
| 5 | TU-packed teeth (**gated**) | the same | **gene-level** CPB (`conditional_poisson_binomial` on genes with the same p̂) | rate > 2α at α = 0.01 **and** > the level bound at 0.05 (pooled over directions) |
| 6 | size bias (recorded) | `ora_null`, random-**gene** categories | TU hypergeometric | over-rejection recorded: shows q weighting is needed even under the null |
| 7 | partial-membership dilution (recorded) | planted gene-level enrichment in random-gene categories | TU-CPB vs gene-level CPB | power of both recorded |
| 8 | Wallenius (recorded) | arms 1–5 | `p_classic` | rates recorded (why it is not the gate) |
| 9 | truth (recorded) | arms 1, 2 and 4 with p_true | TU-CPB | level recorded |

**Spec-writer pre-check** (illustrative only, not gate evidence; `$GIAP_ROOT/stat/m8-04-spec/
sim.py`, 1,500 genes, TU sizes geometric(0.45) capped at 8, true marginal p as p̂). Rejection
rates are at α = 0.01 / 0.05, pooled over 20 categories per replicate:

| Scenario (R) | TU-CPB over | TU-CPB under | gene-CPB over | gene-CPB under | TU hypergeometric |
|---|---|---|---|---|---|
| null, random 40-gene categories (150) | 0.003 / 0.032 | 0.008 / 0.036 | 0.005 / 0.036 | 0.005 / 0.035 | over **0.168 / 0.400** (arm 6) |
| attack 2, whole-TU categories, shift −1.5 (600, 2 per replicate) | 0.003 / 0.027 | 0.003 / 0.027 | 0.003 / 0.036 | 0.002 / 0.022 | under **0.388 / 0.667** (arm 3) |
| attack 2, random-gene categories (600) | 0.008 / 0.040 | 0.007 / 0.033 | at level | at level | under 0.004 / 0.033: **no teeth**, because of dilution |
| TU-packed, c = 1.5 (150) | 0.004 / 0.025 | 0.006 / 0.039 | **0.021 / 0.070** | **0.020 / 0.065** | — |
| TU-packed, c = 2.5 (150) | 0.007 / 0.038 | 0.014 / 0.043 | **0.037 / 0.093** | **0.035 / 0.093** | — |

The pooled bounds at m = 3,000 are 0.0136 / 0.058.

- **TU-packed teeth.** At c = 2.5 the gene-level CPB rejects above 2α at α = 0.01. At α = 0.05 it
  is above the bound but below 2α (0.093 < 0.10). So arm 5's criterion is written as above, not
  as a flat "> 2α".
- **Attack 2 teeth** need TU-aligned categories. With random-gene categories the confounded genes'
  TUs are detected through their other genes, and the unweighted test does not fail.
- **The TU-CPB under tail at c = 2.5** is 0.014 at α = 0.01, at the bound. This is where the q
  independence assumption (Q4) is weakest, and M8-15 must look at it closely.

## Hardening owed (pass 2)

- The CI on `effect`. Candidates: the TU-block bootstrap of k − E (as §7.7.2), or the conditional
  MLE of a log-odds shift δ for C's TUs, with an exact CI by inverting the CPB test (statistician's
  suggestion, M8-01 note).
- Golden: the frame's gene-level counts (`n_in_set_genes`, `k_detected_in_set_genes`, M, N_genes)
  equal goseq's per-category counts on the M8-12 inputs. The CPB p has no R counterpart. Its
  oracle is the enumeration and rational tests above.
- Downstream text owed (not this issue's files; for the lead to assign):
  - M8-06's `method` enum `wallenius_ora` → `cpb_tu_ora`, plus the schema's new columns and flags;
  - M8-14 gains `ora_attack_tu_packed` and TU-aligned attack-2 categories;
  - M8-15's table gains arms 4–9 and drops its "ω forced to 1" teeth, which arm 3 replaces;
  - the M8-index row 4 and its "Wallenius ω" estimator row.

## Questions for the lead

Revision-1 questions 1–6 were accepted with their proposed defaults (M8-index lead decision 1):
- the universe = all eligible genes;
- no elim for `under`;
- the size filter on the full set;
- six collection values;
- `parent_group` = the most general significant ancestor;
- effect CI NA in pass 1.

They still apply. The new questions:

1. **Membership rule.** The brief says "a TU is in a set when any of its genes is", consistent with
   M8-05. M8-05 applies it to *detected* genes, which would make ORA membership depend on
   detection. **Proposed:** any **eligible** gene. The rule matches M8-05 but runs over the
   eligible genes, as stated in 2.
2. **elim removal unit.** **Proposed:** remove the significant child's **TUs** from its ancestors
   (example 7). The alternative removes genes, which is M8-05's redundancy rule and leaks a
   shared TU's detection into the ancestor. Should M8-05's redundancy re-test switch to TU removal
   too, for symmetry?
3. **Columns.**
   - Count columns `n_universe`, `n_in_set`, `n_detected`, `k_detected_in_set` are in **TUs**.
   - New columns `n_in_set_genes`, `k_detected_in_set_genes` and `p_noelim` are added. Adding
     columns is non-breaking, as with ADR-A.
   - `p_classic` = Wallenius, per the ADR's wording. This drops §7.7.3's sense of `p_classic` as
     "the non-elim p of the reported test"; that sense is now `p_noelim`.

   Accept, or rename (for example `p_wallenius` for Wallenius and `p_classic` for the non-elim CPB
   p)? Either way the §9.3 dictionary note needs an ADR-0029 amendment line.
4. **q_t = 1 − ∏(1 − p) assumes the genes of a TU are detected independently given p̂.** Real
   co-detection is positively correlated, so this overstates q. The pre-check suggests the
   conditioning on N absorbs most of it (arm 4 at level up to c = 2.5, with the under tail at the
   bound). **Proposed:** keep ADR-0029's formula and let arm 4 decide. If arm 4 fails, an ADR on a
   TU-level detectability fit (p̂ for "any gene detected" fitted on TU covariates) goes to Khoa.
5. **BH family.** The brief says "BH per collection", while §7.7.3, §8.1 and revision 1 say per
   (collection, direction). **Proposed:** keep (collection, direction). Was a single family per
   collection, pooling both directions, intended?
6. **Size filter unit.** **Proposed:** [5, 500] in **TUs**, as GSEA counts sizes. Sizes in genes
   would test some 2-TU operon-sized sets of 5 genes, and the TU null has little resolution there.
7. **Effect.** **Proposed:** `effect` = k − E[K | N] in TUs (the mean under the conditional null,
   consistent with the p). Revision 1's gene-level k − Σ_C p_detect stays in `category_summary`.
8. **PR split.** At ≈ 390 code lines this is near the ~400 limit. **Proposed:** one PR. The
   alternative is M8-04a (the `stats` primitive, ≈ 100 lines + 220 test lines) before M8-04b
   (`ora.py`/`summary.py`), which lets the statistician review the null first.
9. **Arm 5 threshold.** Proposed: "> 2α at 0.01 and > level bound at 0.05" (the pre-check shows a
   flat "> 2α" at 0.05 fails at c = 2.5). Is it acceptable as a gated teeth criterion, or should
   M8-14's `c` be raised until "> 2α" holds at both levels?
