# Statistician note: join-count scale in `autocorrelation.tsv` and V10 (2026-10-24)

Trigger: plot survey of the MVP-04 K-12 run (`local/runs/mvp04/20261023_a1`). The `join_count_bb`
rows have `value` = 1,973 at lag 1, while `null_mean` / `null_lo` / `null_hi` = -3.4 / -136 / 122.

## Verdict

The p-values and q-values are correct. The problem is that two scales are mixed in one display.
- In the producer it is by design: ADR-0020 Amendment 1 and lead decision 17 say `value` = raw BB_k,
  while `null_*` and `effect` are on the Pearson-residual statistic S_k.
- In V10 it is a real bug, because the figure compares the two scales.

- **Producer.** `src/giap/context/autocorrelation.py:291-293` sets `value` to the raw count `raw[col]`
  and `observed` to the Pearson S_k `jc_obs[col]`. It then calls `lag_test(observed, null)` with the
  Pearson null. `_row` (lines 201-202) writes `value`, the S-scale `null_*` and
  `effect = S_k - null_mean`.
- **The p is S against S.** S_k itself is in no column; it can be recovered as `effect + null_mean`.
  For K-12 at lag 1 that gives 1117.63 - 3.36 = 494.54 + 619.74 = 1114.3 under both nulls, so the
  rows are consistent. The 855 in the trigger comes from 1,973 - 1,117.6, which mixes the two scales.
- **V10 is the bug.** Two places plot the raw `value` against the S-scale envelope:
  - `src/giap/viz/correlogram.py:143-145`, with the band drawn at 133;
  - the Plotly twin `src/giap/viz/twins_panels.py:77` (band) and `:98` (the hover text "envelope").

  The M9-09 spec (panel (b): "raw BB value, ADR-0020", same layout as Moran's I) asks for this.
  It contradicts docs/methods/autocorrelation.md: "[the envelope is] on the statistic's scale ... not
  the raw count BB_k".

  The effect: under a null, the BB marker sits outside the band at every lag, and the filled or open
  marker (from q) disagrees with its position.

## Replay check (giapsim, real `replicon_autocorrelation`, same rng stream)

1,000 genes, B = 999, 100 replicates per scenario, 20 lags, both nulls. The scenarios are
`autocorr_calib.label_null` and `tu_null`.

| Check | label_null | tu_null |
|---|---|---|
| `value` == raw BB_k | all rows | all rows |
| `effect + null_mean` == S_k | all rows | all rows |
| `p_value` == perm_p(S_k, S null) | all rows | all rows |
| rejection at 0.05, table p, gene_label (lag 1 / all lags) | 0.05 / 0.048 | 0.05 / 0.048 |
| rejection at 0.05, table p, tu_block (lag 1 / all lags) | 0.03 / 0.049 | 0.03 / 0.046 |
| rejection at 0.05, hypothetical p(raw vs S null) | 1.0 / 1.0 | 1.0 / 1.0 |
| raw `value` outside its own `null_lo`–`null_hi` | 100 % of rows | 100 % of rows |
| S_k outside the band | 4.8–5 % | 5 % |

- The raw-count null mean is about 280 here, against an S null mean of about 0. A raw-against-S p
  would be 1/(B + 1) on every row.
- So the calibration suites would have caught a p bug: `tests/calibration/test_autocorrelation_calibration.py`
  gates `p_value` from `autocorrelation_rows` via `autocorr_calib.lag_p`.
- They are blind to the display scale. No calibration, unit or V10 test checks that `value` and
  `null_lo` / `null_hi` share a scale. The V10 unit fixtures, e.g. `tests/unit/viz/test_v10_correlogram.py:329-330`,
  use same-scale synthetic numbers.

## Fix to specify (needs a lead decision; touches the §9.3 table and M9-09)

Preferred option: keep the table and fix the figure.
- Panel (b) of V10 and its twin plot S_k = `effect + null_mean` against the bands.
- The y label becomes "BB join statistic (Pearson residual)".
- The hover shows both S_k and the raw BB_k.
- M9-09 panel (b) is amended to match.
- No §9.3 change, no schema bump.

Alternative: make the table self-consistent.
- Add the raw-count null summaries (`null_mean`, `null_lo` and `null_hi` of BB_k under the same
  permutations) in new columns. Or write S_k to `value` and the raw count to a new
  `value_raw_bb` column.
- This is a schema change (minor version for added columns, major for a redefined `value`), so it
  needs an ADR-0020 amendment.

Either way the `value`, `null_mean` and `null_hi` schema descriptions in
`src/giap/schemas/tables/autocorrelation.schema.json` should state the scale for `join_count_bb` rows.

## M5 gate

No re-run is needed. The gated quantities (`p_value`, the BH `q_value`) use the S statistic against
the S null, which the replay confirms. Only the V10 figure (and its twin) of existing runs is
misleading. Re-render it after the fix; MVP-04 needs no recompute.

Replay script (not committed):
`/tmp/claude-430805/-home-classes-bios-270-khoa-giap/b9f44921-ae38-43a6-9b00-e1601c824b49/scratchpad/jc/check.py`
(`python -I check.py 100`, about 1 min on 2 threads).
