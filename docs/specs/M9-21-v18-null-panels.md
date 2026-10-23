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
3. **Test-author (2d7fdb4e): gaps 1–10 accepted as proposed.**
   - **Plots-step wiring** (gap 7): `plots` declares the M9-17 step outputs `context_stats_support`,
     `determinants_support` and `function_support` as optional inputs. Its view finds files in
     those directories directly or under `support/`.
   - **V14:** M9-19 uses the same `function_support` key. Whichever of M9-19 or M9-21 merges second
     reconciles the declaration.
   - **Nextflow:** wiring main.nf → REPORT/PLOTS for these keys follows once M9-17's outputs are
     published (test_report_nf_input_keys.py).
   - **`sa.COLUMNS["null_regions"]`:** the lead adds `contig` to it after the M9-17 implementer
     reports.
4. **Review round (code review of aa50e573).**
   - **B1, split regions:** key regions on (region_id, contig) in the ranking, the null filter and
     `observed`. When a region spans more than one contig, the family is `"<rid> (<contig>) / <stat>"`.
     Tests cover a split region in both the exhaustive and the grid file; the grid step CDF stays
     monotone per contig.
   - **B2, island null:** draw one genome-wide multiple, family `genome / abs`:
     - its values are the replicate-wise maximum of `max_scan_s` over the scanned replicons (NaN skipped as `scan_with_nulls` does);
     - its observed value is the maximum `observed_max`;
     - its caption is "FWER p <min islands p_value>".
     
     Per-replicon multiples are captioned "per-replicon maximum", with no FWER p.
   - **S1, cap:** the genome-wide multiple plus at most 7 per-replicon multiples, ranked by
     `observed_max` descending with ties broken by contig. Write a `metric = n_not_drawn` row.
   - **Nits taken:**
     - N1: "observed NA" in the caption;
     - N2: x labels "max |S_w| over windows" and "mean CV R² (null: permuted response)";
     - N3: read the key columns as str;
     - the ADR-0034 uncalibrated note is unconditional;
     - "p not available" replaces "p NA";
     - tests for `function_support` reaching V14 through the plots step, and for plots/support taking precedence over an analysis directory;
     - simplify `run_view.py:279`;
     - docstring wording "correctly rounded, like the schema reader".
   - **Backlogged:** p ≥ 0.995 printing as "1" (use `.3g` above 0.1).
