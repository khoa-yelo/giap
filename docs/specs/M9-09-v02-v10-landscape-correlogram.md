# M9-09: V2 `abundance_landscape` and V10 `correlogram`

Milestone M9 · **pass 1** · depends on: M9-01, M9-03 (oriC/ter lines) · estimate ≈ 340 code lines
(+ ≈ 320 test lines) · **stats (light)**: V2 draws a LOESS trend; the statistician confirms it is
the display analogue of the islands detrend

## Spec sections implemented

- §10.2 V2: "per replicon: log2 abundance vs. position (detected genes as points, undetected as a
  rug), LOESS trend, shaded islands, oriC and terminus lines".
- §10.2 V10: "Moran's I vs. lag 1–20 with gene-label and TU-block 95 % envelopes; join-count panel
  for detection".
- §7.3.5(c) (Moran's I, join counts, lags `context.moran_lags`, BH across lags within each null,
  "2.5–97.5 % envelopes for V10"), §7.3.5(d) (islands, `detrend_span_frac`), ADR-0020 (raw `value`,
  `effect` = statistic − null mean, `ci_unestimated`, the envelope `null_lo`/`null_hi` "is reported
  for the V10 figure"; residual join count), ADR-0021 Amendment 1 (no `uncalibrated` flag), M5-01
  (`giap.stats.smooth.loess_detrend`).

## V2 (one panel row per replicon with ≥ 1 gene; double width)

- **Points:** genes with `detected = true`, `abundance_imputed = false`, `partial = false`;
  x = gene midpoint in kb (`(start + end) / 2`, unwrapped for origin-crossing genes then taken
  mod length), y = `abundance_log2`, colour `#0072B2`, size 2 pt.
- **Rug:** genes with `detected = false` (status `not_detected`), `#BBBBBB` ticks at the bottom;
  `single_peptide_unestimated` genes as `#E69F00` ticks.
- **Trend:** `loess_detrend(y, pos, frac = context.islands.detrend_span_frac, circular =
  (topology == circular), period = length)` on the drawn points ordered by midpoint; the `fit`
  is drawn as a black line (circular replicons: the line is drawn in position order and not
  joined across the origin). No trend when the replicon has fewer points than the M5-01 minimum
  (`loess_span_points` < 3).
- **Islands:** spans of `islands.tsv` rows with `p_value ≤ stats.alpha.context` shaded with the
  §10.1 island colour and hatch (alpha 0.25); not-significant rows are not shaded.
- **Landmarks:** vertical solid line "oriC" and dotted line "ter" from
  `plots/support/replicon_landmarks.tsv` when `status = called`.
- Axis labels "position (kb)" (Mb when length > 2 Mb) and "log2 abundance".

Data TSV `V02_abundance_landscape.tsv`: `contig`, `kind` enum{point, rug, trend, island, oric, ter},
`locus_tag`, `position_bp` (1-based midpoint), `start`, `end` (islands), `value`
(abundance or fit), `category` (detection status or island direction). Sorted by contig, kind (enum
order), position, locus_tag.

Title: "`<h>` high and `<l>` low abundance islands at FWER ≤ `<alpha>`" plus, when the lag-1
gene-label Moran's I row has a finite `q_value`, "; lag-1 Moran's I `<I>` (q = `<q>`)". With no
abundance: "No measured abundance: detection only" (points absent; rug drawn for all genes).

## V10 (one column per replicon with finite rows; double width)

- **(a) Moran's I:** `autocorrelation` rows `statistic = morans_i`, x = `lag`, y = `value`; for
  each `null_model` ∈ {gene_label, tu_block} a shaded band `null_lo`–`null_hi` (gene label light
  grey, TU block `#56B4E9` at alpha 0.3) and a dashed `null_mean` line. Points with `q_value ≤
  stats.alpha.context` under the **TU-block** null are filled; others open.
- **(b) Join count:** the same layout for `statistic = join_count_bb` (raw BB value, ADR-0020),
  y label "BB joins (detected–detected)".
- x axis "lag (genes)", ticks 1–20 (`context.moran_lags`).
- Rows with `value` NA (flags `not_estimable`, `few_values`): no point; if every row of a contig is
  NA the contig column shows "not estimable (`<flags>`)".
- `oriter_slope` rows are not drawn here (they are context-section table rows).

Data TSV `V10_correlogram.tsv`: `contig`, `statistic`, `null_model`, `lag`, `value`, `null_mean`,
`null_lo`, `null_hi`, `q_value`, `significant` (bool), `flags`; sorted by contig, statistic,
null_model, lag.

Title: "Lag-1 Moran's I `<I>`: `<above|below|within>` the TU-block envelope (q = `<q>`)" for the
first contig in FASTA order with a finite row; else "Autocorrelation not estimable on this genome".

## Degenerate cases

| Situation | Behaviour |
|---|---|
| validation run: all autocorrelation rows NA, `few_values` | V10 empty-state per contig; title "not estimable" |
| a replicon with < 3 drawn points | V2 panel with points/rug, no trend |
| no islands table | V2 without shading, noted |
| an island crossing the origin (`end > length`) | two shaded spans |
| linear contig | `circular = False` LOESS; no padding |
| no abundance | V2 rug for every gene, no points or trend |
| no `gene_context`-derived oriC (sim without oriC) | no landmark lines |
| > 20,000 genes on a replicon | V2 points binned: 20,000 bins with the bin median of measured values and the bin count in the TSV (`kind = point`, `locus_tag` NA) |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | chr1 L = 100,000; detected measured genes at midpoints 10k, 30k, 50k, 70k, 90k with log2 1, 2, 3, 2, 1; frac 1.0, circular | five `point` rows; five `trend` rows equal to `loess_detrend(...).fit` (test calls the engine as oracle) |
| 2 | gene undetected at 40k | one `rug` row category not_detected |
| 3 | islands ISL_1 high 20,001–40,000 p 0.01; ISL_2 low p 0.30; alpha 0.05 | one `island` row; title starts "1 high and 0 low abundance islands at FWER ≤ 0.05" |
| 4 | morans_i lag 1 value 0.31, tu_block null_lo −0.05, null_hi 0.12, q 0.002 | V10 row `significant = true`; title "Lag-1 Moran's I 0.31: above the TU-block envelope (q = 0.002)" |
| 5 | lag 2 value 0.08 inside both envelopes, q 0.4 | open marker, `significant = false` |
| 6 | validation run `OUT/` | V2 rug only; V10 empty state |

## Acceptance criteria

- [ ] `tests/unit/viz/test_v02_landscape.py`: examples 1–3; origin-crossing island; imputed values
  never drawn; trend equals the engine's fit to 1e-12.
- [ ] `tests/unit/viz/test_v10_correlogram.py`: examples 4–6; envelopes come from the table only
  (no null recomputed: the test patches `giap.stats.autocorr` to raise if called).
- [ ] Statistician note in the PR: the V2 trend is gene-level on measured detected genes; the
  islands detrend is on TU scores; the plot states "trend (LOESS, span `<frac>`)" in the legend and
  makes no claim from it.

## Questions for the lead

1. V2's trend on genes vs the islands' TU-level detrend: **Proposed default:** gene-level LOESS
   with the same span fraction (what the reader sees is genes); legend says LOESS span; no claim.
2. Which null decides the filled markers in V10? **Proposed default:** TU block (the conservative
   null ADR-0020 gates); both envelopes are drawn.

## Lead decisions (2026-10-22)
1. V2's trend is a gene-level LOESS with the islands' span fraction (what the reader sees is genes);
   the legend names the LOESS span; no claim is made from it.
2. V10's filled markers follow the TU-block null (the conservative null ADR-0020 gates); both
   envelopes are drawn.
3. Branch off main: M9-03 (support files, #210) is not merged, so V2's oriC/terminus lines are tested
   through their empty state only; the populated lines are pinned when M9-03 lands.

## Lead decisions on the tests (2026-10-22, after 50bea79)
4. V10's title quotes the TU-block row's q (consistent with decision 2); V2's title clause keeps the
   gene-label q, as the spec says.
5. V10's data TSV keeps every Moran's I and join-count row, NA rows included, with `significant` false.
6. V10 rows sort `statistic` in the schema's enum order (not alphabetical).
7. The `<flags>` in "not estimable (<flags>)" are the row's flags joined by ", "; only the prefix and
   `few_values` are pinned.
8. V2 binning assigns genes to bins by midpoint; the bin count is not a TSV column in pass 1 (backlog).
9. "No measured abundance: detection only" applies whenever no gene qualifies as a point (detected,
   measured, not imputed or partial) — so the validation run, whose undetected genes carry only
   imputed values, gets it.
10. The "no islands table" note wording is free; it mentions islands.
11. An origin-crossing island is one data row with end > length (index decision 11), drawn as two spans.
12. The interface in the test docstrings (V2 axes labelled by contig, V10 axes "a:<contig>"/"b:<contig>",
    `%g` numbers, `position_bp = (start0+end)/2 mod L + 1`) is accepted.

## Lead decisions on the review round (2026-10-22, review of 1f5fba6)

21. **Not-estimable rows are NA.** An `autocorrelation.tsv` row flagged `not_estimable` or
    `few_values` is not drawn, whatever its value (the producer writes the join count as 0), and
    counts as NA for the empty-state check (example 6). The test helper now writes those rows as
    the producer does.
22. **V2 legend** never overlaps the title or the data, with or without the oriC/ter entries
    (measured with FigureCanvasAgg after draw); place it below the axes or reserve the space.
23. **V10 with many replicons** wraps its contig columns into rows of at most 3, the figure height
    growing per row; no text overlaps at 6 replicons.
24. **The detection-only rug** (decision 9 stands): when no gene is a point, the rug holds every
    gene, `all_rug = not measured`. The property test is re-pinned by the lead.
25. Lead decision 6's statistic order is pinned by a new test (the alphabetical-sort mutant C7).
26. NITs (untested lines, the 20k-point memory bound — fit on ≤ 5,000 bin medians later —, rug-only
    y ticks, the V10 legend source, the size trims) go to the backlog; the statistician note (TU
    ordinals in the islands detrend vs bp positions with period L in V2) goes into the PR body.
