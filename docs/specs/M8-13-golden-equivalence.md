# M8-13: Gate, half 1 — GIAP's Wallenius and fgsea against goseq and R fgsea on the golden data

Milestone M8 · track C (calibration and golden) · **pass 2** · depends on: M8-01, M8-11 (and M8-02),
M8-12 (committed golden data) · estimate ≈ 40 code lines (+ ≈ 300 test lines) · **statistical**

## What the gate says

§13 row M8: "GSEA and ORA match fgsea and goseq on golden data". §7.7.5 states the criteria and
`tolerances.yaml` (protected, **unchanged**) holds them:
`golden_r: {wallenius_rel_tol: 1.0e-6, fgsea_es_abs_tol: 1.0e-9, fgsea_nes_abs_tol: 0.02,
fgsea_log10p_k: 3.0}`.

1. "Wallenius p-values against goseq run with a supplied PWF, relative difference ≤ 10⁻⁶";
2. "fgsea ES identical to 10⁻⁹";
3. "NES within 0.02";
4. "|log₁₀ p − log₁₀ p_R| within 3 × the combined `log2err` / log₂ 10".

## Behaviour

`tests/golden/test_golden_wallenius.py`, `tests/golden/test_golden_fgsea.py` (layer "Golden",
§11.2: every PR, < 2 min), reading `tests/data/golden/r/` and `tolerances.yaml`.

1. **Wallenius.** For every golden category: GIAP's counts (M, n, N, k via `wallenius_ora`) equal
   goseq's `numInCat`, `numDEInCat`; ω equals goseq's weight identity (M8-01 example 3) to 10⁻¹²;
   `|p − p_R| / p_R ≤ 1e-6` for **both** tails. p_R = 0 (if BiasedUrn underflows) requires
   p ≤ 10⁻³⁰⁰. Every category is checked; none is skipped.
2. **ES** (against R's `calcGseaStat` column and fgsea's `ES`): `|es − es_R| ≤ 1e-9` for every set,
   on the same integer-scaled, stably sorted stats (M8-02 step 1). The same sets are dropped on
   both sides (names absent, size filter).
3. **NES**: `|nes − nes_R| ≤ 0.02` where both exist; a set where exactly one side is NA fails.
4. **p**: `|log10 p − log10 p_R| ≤ k · sqrt(log2err² + log2err_R²) / log2(10)` with k = 3. Rows
   where either `log2err` is NA (fgsea's `isCpGeHalf` case) are reported, counted, and must be
   ≤ 5 % of sets; they are not silently passed.
5. **Report.** The test writes `local/logs/golden/m8_golden.json`: per criterion the max
   deviation, the count per p decade (10⁰…10⁻³⁰), and the number of sets/categories, so a vacuous
   pass (all p near 1) is visible. The numbers go into PROGRESS.md as gate evidence.
6. **GIAP seeds.** The Python fgsea runs with `giap.rng.for_step("function.gsea", GOLDEN_SEED)`;
   R's seeds are in `PROVENANCE.md`. No attempt is made to match R's RNG.

## Examples (what a passing run looks like)

| # | Check | Expected |
|---|---|---|
| 1 | 150 golden categories × 2 tails | max relative difference ≤ 10⁻⁶, including the ~10⁻⁴⁰ planted tails |
| 2 | 120 sets ES | max \|Δ\| ≤ 10⁻⁹ (in practice ~10⁻¹⁵: same integer arithmetic) |
| 3 | NES | max \|Δ\| ≤ 0.02 |
| 4 | p, the 20 planted sets | every \|Δ log10 p\| within 3 × combined `log2err`/log₂10; decades 10⁻² … 10⁻³⁰ all populated |
| 5 | a deliberately broken port (sign flip in the reversed-ranks ruler, mutation test) | criterion 4 fails on the negative planted sets |

## Acceptance criteria

- [ ] Both golden test files pass in CI (GitHub, no R needed: the data are committed).
- [ ] The JSON report shows every p decade from 10⁻² to 10⁻³⁰ populated by ≥ 1 set.
- [ ] Mutation check (example 5) documented: two injected faults are caught.
- [ ] PROGRESS.md records gate half 1 as met with the four max deviations and the test names.

## Hardening owed / escalation

- If criterion 1 fails because scipy's and BiasedUrn's pmf precision differ (index risk): an ADR on
  the computation (log-space pmf or a precision argument), never on `wallenius_rel_tol` (R3).
- If criterion 4 fails after 3 iterations (R9): the §7.7.4 / §17 fallback ADR goes to Khoa
  (index "Questions that need Khoa" 1); GSEA keeps `fgsea_simple` meanwhile.

## Questions for the lead

1. Rows with NA `log2err` on either side: counted and capped at 5 % (proposed), or excluded
   without a cap?
2. Should the golden test also compare the leading edges (not in §7.7.5)? **Proposed default:**
   yes, as an exact set comparison — it is deterministic given ES and costs nothing.
