# M9-15: Interactive twins of V2, V5, V6, V7, V8, V10, V11, V17, V18

Milestone M9 · **pass 2** · depends on: M9-14 and the static plots (M9-04 … M9-09, M9-11) ·
estimate ≈ 380 code lines (+ ≈ 300 test lines) · not statistical

## Spec sections implemented

§10.2 "Interactive in report" column: Plotly for V2, V6, V7, V8, V10, V11, V17, V18; "Plotly
Sankey" for V5. V14, V15 twins land with their plots (M9-19, M9-20); V16's concordance twin with
M9-23. §10.1 colours apply to the twins (the style module exports a Plotly template).

## Behaviour

1. `giap.viz.style.plotly_template()` returns a `plotly.graph_objects.layout.Template` with the
   §10.1 colours, DejaVu Sans (falls back to sans-serif in the browser), the same font sizes, white
   background; every twin uses it.
2. Each twin reads only its data TSV frame (M9-14 item 1). Hover text shows the IDs and the
   numbers with their error estimates as the table writes them (e.g. "prophage_1 · O − E −4.75
   (95 % CI −7.1 to −2.2) · q 0.003"); `NA` CIs show "CI unestimated".

| Plot | Twin |
|---|---|
| V2 | scatter per replicon (detected points), rug as a marker row, trend line, island shapes; hover locus tag, product not available (not in the data TSV) |
| V5 | `go.Sankey` with the fixed node order and colours of M9-05 |
| V6 | three subplots as the static; the alpha line as a shape |
| V7 | horizontal error-bar scatter in the drawing order; q in hover |
| V8 | box-free strip/violin (`go.Violin` with points when n > 50, `go.Scatter` strips otherwise) |
| V10 | lines with filled envelopes per null |
| V11 | histograms (pre-binned bars from the TSV), lfdr line, partner scatter with flag symbols |
| V17 | five subplots; the warning list as a table trace |
| V18 | reliability scatter + diagonal; p-histograms as pre-binned bars |

3. Twins never recompute: histograms use the TSV's pre-binned rows, not Plotly's binning.
4. Element cap: a twin with more than 20,000 points uses the same binned rows as the static plot.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| a data TSV that is header-only (empty state) | no twin; the static empty state is shown |
| hover fields NA | rendered as "NA" |
| a plot skipped | no twin |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | V05 TSV of M9-05 example 1 | Sankey with 2 + 5 + 4 nodes, 5 links with values 1 |
| 2 | V07 TSV of M9-06 example 1 | 2 region traces; hover of `prophage_1` contains "q 0.003" |
| 3 | V18 TSV with pre-binned counts | bar heights equal the counts; no `histfunc` used |
| 4 | the twin JSON for any plot | colours ⊂ the §10.1 table + viridis/cividis |

## Acceptance criteria

- [ ] `tests/unit/viz/test_twins.py`: per twin, trace values equal the data TSV columns; examples
  1–4; the template's colour set.
- [ ] Every plot marked "Plotly" in §10.2 (except V14, V15, V16) has a registered twin (a test
  reads a list derived from §10.2 kept in `giap.viz.interactive`).

## Questions for the lead

1. V2 hover would be more useful with the gene product. **Proposed default:** add `product` to the
   V2 data TSV (a column, written by M9-09's code path) so the twin stays TSV-only.

## Lead decisions (2026-10-22)

1. V2's data TSV gains a `product` column (written by M9-09's code path), so the twin's hover reads
   the TSV only.
2. M9-15 is stacked on M9-14 (#221, queued); after #221 squash-merges, the branch merges the old
   base tip, then origin/main.
3. Feature loop: test-author → implementer → code review.
