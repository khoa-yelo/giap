# M9-07: V18 `stat_diagnostics`, panels (a) reliability and (b) p-value histograms

Milestone M9 · **pass 1** · depends on: M9-01 · estimate ≈ 220 code lines (+ ≈ 250 test lines) ·
not statistical (histograms of reported p-values; the reliability points are M4-13's)

## Spec sections implemented

- §10.2 V18: "detectability reliability curve; p-value histograms per test family; null
  distributions against the observed statistic (regions, islands, determinants)". Panel (c), the
  null distributions, needs null draws that no table holds: M9-21 (after M9-17).
- §7.3.3 "Calibration check. 10 equal-count bins, Brier score, and the calibration intercept and
  slope … written to `detection_summary.tsv` (section `detectability`) and plot V18."
- §8.1 (test families), §8.2 (pre-declared BH families), §9.2 (test columns).

## Panel (a): reliability

From `detection_summary` rows `section = detectability`: `group_type = bin` rows give, per bin
`group`, `mean_p` (x), `observed_rate` (y) and `n` (marker area); `group_type = overall` gives
`brier`, `calibration_slope`, `calibration_intercept`. Draw the diagonal y = x (grey dashed), the
points joined in bin order, both axes [0, 1], labels "mean predicted p_detect" and "observed
detection rate", and the three overall numbers in the corner. If the overall row
`detectability_constant` has value 1 (the §9.4 flag of a constant model) or there are no bin rows:
the panel shows "detectability constant: no calibration curve".

## Panel (b): p-value histograms per family

A family = one table and one key, drawn as one small histogram (20 bins of width 0.05 on [0, 1],
counts) with "n = …" and the table name in the panel title. Families, in this order, each drawn
only when it has ≥ 1 finite `p_value`:

| Table | Family key | Rows |
|---|---|---|
| `regions_stats` | `stat` (o_minus_e, hl_shift) | `subset = all`, `p_value` finite |
| `autocorrelation` | (`statistic`, `null_model`) | finite `p_value` |
| `pairs` | `test` (≠ none) | finite `p_value` |
| `enrichment_detection` | (`collection`, `direction`) | finite `p_value` (M8; drawn once the table exists) |
| `enrichment_gsea` | `collection` | finite `p_value` (M8) |

Not drawn as histograms: `islands` (`p_value` is a family-wise maximum-statistic p, not one test
per row), `determinants_model` (`null_p`, one value per model: M9-21 panel c), tables without
`p_value`. A grid of up to 4 columns; more than 16 families → the 16 with the most p-values, the
rest listed in the data TSV only.

## Data TSV `V18_stat_diagnostics.tsv`

`panel` enum{a, b}, `table`, `family`, `bin`, `x`, `y`, `n`, `metric`, `value`. Panel (a) bin rows
(`x` = mean_p, `y` = observed_rate, `n`) and overall rows (`metric`, `value`); panel (b) one row per
family × bin (`bin` = left edge, `n` = count). Sorted by panel, table (order above), family, bin.

## Title

"Detectability calibration slope `<s>` (intercept `<i>`), Brier `<b>`" (`%.2g`), or
"Detectability constant; `<k>` test families shown". Never "calibrated": the slope is a diagnostic,
not a test.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no family has a finite p (validation run: every autocorrelation p NA) | panel (b) "no p-values in this run" |
| one p in a family | a one-count histogram |
| `p_value = 1` | in the last bin [0.95, 1] (closed on the right) |
| `p_value = 0` (a bug: §8.2 permutation p ≥ 1/(B+1)) | drawn in the first bin; not hidden |
| `detection_summary` absent | panel (a) "not available"; panel (b) still drawn |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | bins 1–2: (mean_p 0.1, rate 0.12, n 50), (0.6, 0.55, 50); overall slope 0.97, intercept 0.02, brier 0.11 | 2 bin rows; title "Detectability calibration slope 0.97 (intercept 0.02), Brier 0.11" |
| 2 | `regions_stats` o_minus_e p = 0.01, 0.04, 0.5, 1.0 | family `regions_stats / o_minus_e`, n = 4: bin [0, 0.05) count 2, bin [0.50, 0.55) count 1, bin [0.95, 1.00] count 1, all other bins 0 |
| 3 | `islands` table with p_value 0.03 | no histogram for islands |
| 4 | validation run `OUT/` (overall `detectability_constant` 1) | panel (a) text "detectability constant: no calibration curve" |
| 5 | 20 families | the 16 with most p drawn; all 20 in the TSV |

## Acceptance criteria

- [ ] `tests/unit/viz/test_v18_pass1.py`: examples 1–5; bin edges; the family table above is the
  single source (a test lists tables with `p_value` in their schema and asserts each is either a
  family or in the "not drawn" list, so a new test table cannot be forgotten silently).
- [ ] Renders on `tests/data/golden/sim50k/m4/` and on the validation run.

## Questions for the lead

1. Islands' FWER p in panel (b)? **Proposed default:** excluded (above); they appear in panel (c)
   against their maximum-statistic null (M9-21).
2. The value of the `detectability_constant` overall row: the validation run writes 1 (checked:
   `local/runs/val-mvp-2026-10-22`, where `n_detected = 0`). **Proposed default:** treat value 1 as "constant"; the M4-13 owner confirms the
   encoding in review.
