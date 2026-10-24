# M6-05: Discordant OP pairs — the t test, the one-detected test and the BH family

Milestone M6 · track R · **pass 1** (both-detected branch) **+ pass 2** (M6-05H: the one-detected
branch) · depends on: M6-01, M6-04; pass 2 also M6-02 · estimate ≈ 220 code lines (≈ 150 of them in
pass 1) (+ ≈ 300 test lines) · **statistical** (label `stats`)

## Spec sections implemented

- §7.4.3 in full: both-detected two-sided t p; the one-detected P(miss) test ("p = P(miss). A small
  value means the missing partner should have been seen"); "With one run and no peptides, this test
  is skipped and marked `unestimated (no_detection_limit_model)`"; "BH across all tested OP pairs
  (both kinds together). Pairs are discordant at q ≤ `stats.alpha.operon_flags` (0.1)";
  "Relative-scale abundance (§7.3.4): abundance tests are marked `unestimated (relative_scale)`".
- §7.4.2 (which pairs enter the fit), §8.1 row "Operon discordance" (family: pair; null: the
  within-operon t distribution and the detection-limit model; correction BH; α 0.1), §9.2 test
  columns, §9.3 `pairs.tsv` `test` enum value `both_detected_t` and `one_detected_pmiss`.

## Public interface

`src/giap/operons/discordance.py` (pure)
```python
@dataclass(frozen=True, slots=True)
class Discordance:
    fit: StudentT
    ratios: F64                  # the oriented within-operon log-ratios r (§7.4.3)
    rows: pd.DataFrame           # pair_id, test, stat, effect, p_value, q_value, method, flags

MIN_OP_PAIRS: int = 20

def fit_rows(pairs: pd.DataFrame) -> F64: ...          # Δ of OP pairs usable for the fit
def discordance(pairs: pd.DataFrame, genes: pd.DataFrame, cfg: GiapConfig,
                p_miss: PMissModel | None) -> Discordance: ...
```

1. **Eligibility.** A pair enters any test only when both genes are eligible (§7.3: not `partial`)
   and the class is OP. Pairs with `rna_between` are tested normally (the flag is informational).
2. **The fit.** Δ over OP pairs with `detected_a and detected_b`, both `abundance_log2` finite and
   neither `abundance_imputed`. `n < MIN_OP_PAIRS` (20) or `not fit.converged` → every OP test is
   `test = "none"`, `flags` gains `too_few_op_pairs`, and **no BH family is formed** (decision 8).
3. **Both detected.** `test = "both_detected_t"`, `stat = delta_log2`, `effect = delta_log2 −
   fit.loc`, `p_value = student_t_two_sided(Δ, fit)`, `method = "student_t"`,
   `null_model = "within_operon_t"`, `n_perm` NA. `ci_low`/`ci_high`: the `fit.scale`-based
   prediction interval is **not** a CI of the effect, so they stay NA with `flags` gaining
   `unestimated` — §8.2 asks for a CI and none exists for a single pair (Q2).
4. **One detected** (pass 2, M6-05H). `test = "one_detected_pmiss"`, `stat` = the detected
   partner's `abundance_log2`, `effect` NA, `p_value = p_miss(u, a_v, ratios)` from M6-02,
   `method = "detection_limit"`, `null_model = "p_detect_model"`. `ratios` are the oriented r of
   §7.4.3: Δ when the undetected gene u is downstream, −Δ when upstream, taken from the same pairs
   as the fit.
5. **Neither detected**, or an abundance missing under `abundance.missing = exclude`:
   `test = "none"`, no p.
6. **Relative scale.** When the run's `abundance_scale` is `relative` (§7.3.4; from the gene table's
   `abundance_type` / the detection summary), every abundance test is skipped: `test = "none"`,
   `flags` gains `relative_scale`. The pair table still carries classes, distances and evidence.
7. **BH.** One family over all OP rows with a non-NA p (both kinds together), `giap.stats.bh`;
   `q_value` written back. `discordant = q_value <= stats.alpha.operon_flags` (0.1) is not a column
   — M6-07 reads the q.
   *Superseded (M6-05H, 2026-10-24): per M6-index lead decision 1 and M6-02 lead decision 4 the family holds only the `both_detected_t` p; P(miss) rows get no q — item 7's "both kinds" and example 8's "one family of 150" stand only once M6-09 attack C and an ADR admit the mixture.*
8. **No RNG.** Everything here is deterministic; the only stochastic part of M6 lives in M6-02.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no OP pair at all | empty family; `operon_flags.tsv` has no discordance row |
| < 20 usable OP pairs | `too_few_op_pairs`; every OP test `none` (so sim100k typically yields none — expected, M6-08) |
| fit did not converge | as above, plus a `WARNING` log line with n and the start values |
| all Δ identical | `fit.scale → 0`; `converged = False` → `too_few_op_pairs` path |
| one run and no peptides | one-detected rows get `test = "none"` + `no_detection_limit_model` |
| relative scale | as step 6 |
| `p_miss` returns NaN for a gene (no observable peptide) | that row `test = "none"`, `flags` gains `no_detection_limit_model` |
| every p is NA | `bh` of an empty family; all `q_value` NA |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | 200 OP pairs, Δ ~ t(0.05, 0.4, 6); one pair with Δ = 1.2 | that pair `p_value` ≈ **0.0282**, `test = both_detected_t`, `method = student_t` |
| 2 | same with Δ = −1.2 | `p_value` ≈ **0.0205** |
| 3 | a pair with Δ = `fit.loc` | `p_value` = **1.0**, `effect` = 0 |
| 4 | 19 usable OP pairs | every OP row `test = none`, `flags` contains `too_few_op_pairs`; no q anywhere |
| 5 | 200 usable pairs + 1 pair whose `abundance_imputed_a` is true | the imputed pair is **not** in the fit but **is** tested (it is detected with a value) |
| 6 | pair where only the upstream gene is detected, a_v = 17.0, r sample {−1.0, −0.2, 0.3, 1.5}, M6-02's toy model | `p_value` = **0.6889**, `test = one_detected_pmiss` (pass 2) |
| 7 | the same pair in pass 1 | `test = none`, `flags` contains `no_detection_limit_model` |
| 8 | 100 both-detected rows with p and 50 one-detected rows with P(miss) | one BH family of 150; `q_value` monotone in the pooled p |
| 9 | `abundance_type = maxlfq` | every OP row `test = none`, `flags` contains `relative_scale`; `pairs.tsv` still has classes and distances |
| 10 | the same gene table twice | identical frame (deterministic) |

## Acceptance criteria

- [ ] `tests/unit/operons/test_discordance.py`: examples 1–10 and every degenerate row on
  hand-built frames; examples 1–3 assert the hand-computed constants from M6-01.
- [ ] Property: the BH family contains exactly the rows with `test != "none"`, and no NA p enters it
  (the M4-17 rule).
- [ ] Property: the fit is invariant to row order and to relabelling `pair_id`.
- [ ] Property: adding pairs whose genes are `partial` changes nothing.
- [ ] Calibration (M6-09 owns the suite): OP pairs simulated without discordance give p-values at
  level for α ∈ {0.01, 0.05, 0.10} (§7.4.7 "discordance p-values uniform for OP pairs simulated
  without discordance"); the one-detected branch is measured separately (risk register).
- [ ] Integration: on sim100k the step runs to completion with `too_few_op_pairs` and writes a
  well-formed empty `operon_flags.tsv`.

## Hardening owed (pass 2)

- **M6-05H**: the one-detected branch (step 4) once M6-02 merges — rows, flags, the ratios sample
  and the shared BH family; example 6 becomes a live test.
- Measure the one-detected branch's level separately in M6-09 and decide, with the statistician,
  whether P(miss) may share a BH family with the t p-values (Q3). If it may not, the ADR options
  are: separate families, or report P(miss) with no q.
- `ci_low`/`ci_high` for the effect: a bootstrap over the OP pairs of the TU, if §8.2's "every p has
  a CI" is to be met literally (Q2).
- KS uniformity of the both-detected p under the full suite (R = 1,000), recorded in
  `docs/methods/operons.md`.

## Questions for the lead

1. Minimum usable OP pairs for the fit: 20 (proposed) or 10? **Proposed default:** 20 — below that
   the df estimate is unusable, and the consequence (no discordance calls) is honest.
2. `ci_low`/`ci_high` on a per-pair test: leave NA with `unestimated` (proposed), or report the
   t-prediction interval (which is not a CI of an effect)? **Proposed default:** NA; §8.2's CI rule
   is written for estimated effects, and pretending otherwise would mislead.
3. One BH family for both kinds, as §7.4.3 says, even though P(miss) is not a p-value?
   **Proposed default:** follow the spec in pass 1 and let M6-09 measure it; any change is an ADR.

## Lead decisions (2026-10-21)
1. Q1: minimum **20** usable OP pairs for the t fit; below it every pair is `unestimated` with a flag.
2. Q2: `ci_low`/`ci_high` NA (no effect is estimated per pair). Amended after the statistician
   review of #147 (F1): the tested row carries `ci_unestimated` (ADR-0020 decision 4 precedent),
   not `unestimated` — the claim has p and q; only the CI is missing (R4).
3. Q3: superseded by M6-index lead decision 1 — `P(miss)` never enters the BH family; in pass 1 the
   one-detected branch is `unestimated (no_detection_limit_model)`, so the family holds only t p-values.
4. The t fit is `giap.stats.tdist.fit_student_t` (M6-01, #144). A degenerate fit (`converged = False`,
   NaN parameters — M6-01 decision 13) means every pair is `unestimated`, never p = 0.
