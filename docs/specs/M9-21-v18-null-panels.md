# M9-21: V18 panel (c): null distributions against the observed statistic

Milestone M9 · **pass 2** · depends on: M9-07, M9-17 (null support files) · estimate ≈ 160 code
lines (+ ≈ 180 test lines) · not statistical (draws the nulls the steps computed)

## Spec sections implemented

- §10.2 V18 "null distributions against the observed statistic (regions, islands,
  determinants)". §7.3.5(b) (circular shifts; two-sided p centred on the null mean), §7.3.5(d)
  (maximum-scan null; FWER), §7.8.2 step 4 (TU-level permutation; p = (1 + b)/(1 + B)), ADR-0034
  (`null_uncalibrated`), ADR-0018 (O − E test statistic is studentised; the null file stores what
  the test compares).

## Panel (c)

Small multiples appended to V18 (same figure, a third row block):

- **Regions:** up to 12 regions with the smallest `q_value` (from `regions_stats`, `stat =
  o_minus_e`), each a histogram (30 bins) of `null_regions` `null_value` with a vertical line at
  `observed` and "p `<p>`, q `<q>`".
- **Islands:** per contig × direction, the histogram of `max_scan_s` with the observed maximum and
  the FWER p of the top island.
- **Determinants:** per model (elastic net, LightGBM), the histogram of permuted `cv_r2` with the
  observed mean CV R²; caption "null p `<p>` (uncalibrated, ADR-0034)" while flagged.

Data TSV: `V18_stat_diagnostics.tsv` gains `panel = c` rows: `table` (regions_stats, islands,
determinants_model), `family` (region_id/stat, contig/direction, model), `bin`, `n`, and one
`metric = observed` row per small multiple.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| support files absent | panel (c) "null draws not available"; (a), (b) unchanged |
| a region with a constant null (all shifts equal) | one bar; the observed line |
| `too_few_rows` determinants | no determinants panel |
| null stored as percentiles (M9-17 Q2 fallback) | drawn as a step CDF instead of a histogram, labelled |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | region nulls for 15 regions | 12 small multiples, those with the smallest q |
| 2 | determinants null 999 values, observed 0.21 | one histogram per model; `observed` row 0.21 |
| 3 | no support | panel (c) text only |

## Acceptance criteria

- [ ] `tests/unit/viz/test_v18_null_panels.py`: examples 1–3; the uncalibrated caption.
- [ ] The observed line equals the table value the p was computed from (regions: the studentised
  statistic stored by M9-17, not the raw O − E; a test asserts the column used).

## Questions for the lead

1. Regions draw the statistic the test used (studentised, ADR-0018), not the reported raw O − E.
   **Proposed default:** yes, with the x label "test statistic (studentised O − E)".

## Lead decisions (2026-10-23)
1. **Q1 yes:** regions draw the studentised test statistic (ADR-0018), x label "test statistic
   (studentised O − E)".
2. **Stacking:** on m9/plot-support-analysis (M9-17 null support files; `null_regions` has a
   `contig` column, `null_islands` is per scanned contig, decision 5).
