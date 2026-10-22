# M9-14: The interactive report: Plotly.js inlined once, sortable tables, `--no-interactive`

Milestone M9 · **pass 2** · depends on: M9-02 · estimate ≈ 320 code lines (+ ≈ 280 test lines) ·
not statistical

## Spec sections implemented

- §10.3 "Self-contained. One HTML file built with Jinja2. Plotly.js is inlined once, and CSS and
  the small table-sorting script are inline; no network access is needed to view it. The size
  target is < 25 MB: large tables are linked by relative path, not embedded."
- §10.3 "Tables. Every table section shows the top rows, sortable and filterable, with a link to the
  full TSV." §10.3 tests (no external URLs loaded at view time; link checker).
- §10.2 column "Interactive in report" (Plotly for V2, V5 (Sankey), V6, V7, V8, V10, V11, V14,
  V15, V16 (concordance), V17, V18; static SVG for V1, V3, V12, V13). §2.1 decision 13. §5.1
  `--no-interactive`, §5.2 `report.interactive`.

## Public interface

`src/giap/viz/interactive.py`
```python
InteractiveFn = Callable[[pd.DataFrame, str], "plotly.graph_objects.Figure"]   # data TSV frame, title
def interactive(pid: str) -> Callable[[InteractiveFn], InteractiveFn]          # registry by plot ID
def figure_json(pid: str, data: pd.DataFrame, title: str) -> str | None        # None = static only
```
`src/giap/report/templates/` gains `plotly_block.html.j2` and `tables.js` (inline).

## Behaviour

1. **Source of interactive figures.** A Plotly twin is built **from the plot's data TSV** (the
   exact drawn data), never from the tables again, so static and interactive views cannot
   disagree. Twins are registered by M9-15 (and M9-19, -20, -23); this issue builds the mechanism
   and one twin as a working example (V17).
2. **Embedding.** When `report.interactive` is true and at least one twin exists: one `<script>`
   holding `plotly.offline.get_plotlyjs()` (plotly.js 4.1.1 from plotly 7.1.0, MIT), placed once
   in `<head>`; each twin as `<div id="plotly-<pid>[-<suffix>]">` plus an inline JSON spec and a
   `Plotly.newPlot` call with `config = {displaylogo: false, responsive: true, toImageButtonOptions:
   {format: "svg"}}`, `include_mathjax` off, no topojson (no geo traces exist). The static SVG
   stays in a `<details>` element beside the twin ("static figure").
3. **`--no-interactive` / `report.interactive: false`:** no Plotly.js, no twins; static SVGs only.
4. **Tables.** Every embedded table (top 20 rows, M9-02) becomes sortable (click a header: numeric
   sort when every non-`NA` cell parses as a number, else string; `NA` always last) and filterable
   (one text box per table, case-insensitive substring over the row). About 40 lines of vanilla JS,
   written for GIAP (index decision 24). No external JS.
5. **Size budget.** After writing, the step measures `report.html`; > 25 MB raises a warning event
   `report.size_exceeded {bytes}` and, under `--strict`, exit 3. Tables are always linked, never
   embedded beyond 20 rows.
6. **Licence notice.** The page footer names plotly.js and its MIT licence (NOTICE.md gains the
   line if absent).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no twin registered for any rendered plot | Plotly.js not inlined even when `interactive: true` |
| a twin raises | the `report` step fails (exit 4), as for static plots (M9-01 Q2) |
| a table with 0 rows | header only; the filter box still present |
| a 30-MB report (huge run) | written; warning event; exit 3 only with `--strict` |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | report with V17 twin, `interactive: true` | exactly one occurrence of the plotly.js banner `plotly.js v4.1.1`; one `Plotly.newPlot` |
| 2 | same with `--no-interactive` | no `plotly.js` string; V17 as inline SVG |
| 3 | `external_loads(report.html)` (M9-13 checker) | `[]` (plotly.js's internal URL strings are inside the script text, not loads) |
| 4 | table column values `2`, `10`, `NA`, `1.5` sorted ascending | `1.5`, `2`, `10`, `NA` |
| 5 | filter "prophage" on regions table | only rows containing "prophage" visible |
| 6 | sim100k report size | ≈ 5 MB (measured and recorded in the PR) |

## Acceptance criteria

- [ ] `tests/unit/report/test_interactive.py`: examples 1–3 (string checks on built HTML); the
  JSON spec of the V17 twin equals the data TSV values.
- [ ] `tests/unit/report/test_tables_js.py`: examples 4–5 by running the sort/filter logic through
  a small pure-Python mirror of the JS rules **or** by asserting the JS source implements the
  documented rules (no browser in CI; see Q1).
- [ ] Size measured on sim100k and on K-12 (M9-26) and recorded.

## Questions for the lead

1. Testing the JS without a browser: no headless browser is a dependency. **Proposed default:**
   keep the JS minimal, test the rules through a Python mirror plus a structural check of the
   script, and have Khoa exercise sorting during the M9-25 review.
2. Static SVG beside every twin roughly doubles figure bytes. **Proposed default:** keep both (the
   static one is what was approved and snapshotted).
