# M9-19: V14 `enrichment`: which functions are high or low?

Milestone M9 · **pass 2** · depends on: M9-17 (`gsea_running`); PR #181 and the `FUNCTION`
Nextflow wiring merged · estimate ≈ 340 code lines (+ ≈ 320 test lines) · **stats (light)**: the
running sums come from the support file; the plot asserts they match the table's `es`

## Spec sections implemented

- §10.2 V14: "(a) GSEA running-score plots for the top 6 sets; (b) NES vs. −log₁₀ q dot plot;
  (c) detection ORA, O/E by category"; Plotly twin.
- §7.7.3–§7.7.5, ADR-0025 ("The report and V14 must not present flagged rows as calibrated" —
  `genomic_cluster`), ADR-0026 (GSEA "conditional on detection"; flag `detection_confounded`),
  ADR-0027 (NES without a CI), ADR-0029/ADR-0035 (ORA is TU-level: `stat` = k detected TUs,
  `effect` = k − E[K | N], `method = cpb_tu_ora`), §9.3 `enrichment_gsea.tsv`,
  `enrichment_detection.tsv` (branch `origin/m8/function-step` schemas), §5.2
  `stats.alpha.enrichment`.

## Panels (double width)

- **(a) Running scores**, 2 × 3 small multiples, the six sets of `gsea_running.tsv.gz` (M9-17's
  order): x = rank position, y = `running_es`; set-member ticks below; the leading edge shaded up to
  the peak; panel title "`<term>` `<name>` (`<collection>`) NES `<nes>`, q `<q>`". A set flagged
  `genomic_cluster` or `detection_confounded` gets a grey "not calibrated (`<flag>`)" banner and a
  dashed outline (ADR-0025/0026).
- **(b) NES vs −log₁₀ q:** every row with finite `nes` and `q_value`; colour by collection (six
  Okabe–Ito hues, fixed order cog, kegg, go_bp, go_mf, go_cc, go_slim); hollow markers for
  `genomic_cluster` / `detection_confounded` rows; a dashed line at −log₁₀ `stats.alpha.enrichment`;
  `redundant = true` rows drawn small and grey; no error bars (ADR-0027); labels for the 10 most
  significant non-redundant rows.
- **(c) Detection ORA:** rows with `q_value ≤ stats.alpha.enrichment`, else the 20 smallest q;
  O/E = `stat` / (`stat` − `effect`) on a log₂ axis, one row per (collection, term, direction),
  coloured by direction (over `#0072B2`, under `#D55E00`), `coverage` printed at the right; "CI
  unestimated" (the M8-04 pass-1 rows carry `ci_unestimated`); x label "observed / expected
  detected TUs".

Footnote (all panels): "GSEA is conditional on detection, with no p_detect adjustment" (ADR-0026);
"coverage = share of eligible genes with ≥ 1 term".

## Data TSV `V14_enrichment.tsv`

`panel` enum{a, b, c}, `collection`, `term`, `name`, `direction`, `rank_position`, `running_es`,
`in_set`, `nes`, `q_value`, `o_over_e`, `coverage`, `flags`. Panel (a) copies the support rows;
(b), (c) one row per table row. Sorted by panel, collection (enum order), term, direction,
rank_position.

## Title

"`<n>` sets enriched by abundance at q ≤ `<alpha>` (conditional on detection); `<m>` categories
detected more or less than expected" (counts exclude `redundant` rows). None: "No set or category
passes q ≤ `<alpha>`" — never "no enrichment" (ADR-0026 consequence).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| `function` module off | skip `module_off:function` |
| `layer_missing:go;kegg;cog` (no resources) | header-only tables → empty state "No annotation resources installed (gippa setup --go --kegg --cog)" |
| no abundance (`enrichment_gsea` header-only) | panels (a), (b) "no abundance input"; (c) drawn |
| `gsea_running` absent (old run, or `gippa report` before M9-17) | panel (a) "running scores not available"; (b), (c) drawn |
| support `running_es` max deviation ≠ table `es` (> 1e-9) | the plot raises (a data mismatch is a bug, exit 4) |
| q = 0 (cannot happen: fgsea p ≥ floor) | guarded: −log₁₀ capped at the p floor and noted |
| only `under` rows significant | panel (c) shows them; colours as above |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | ORA row stat 12, effect 4 | O/E 12 / 8 = 1.5 |
| 2 | GSEA row flags `genomic_cluster` | (b) hollow; if in the top 6, (a) banner "not calibrated (genomic_cluster)" |
| 3 | 3 non-redundant GSEA rows and 1 redundant at q ≤ 0.05; 2 ORA rows | title "3 sets enriched by abundance at q ≤ 0.05 (conditional on detection); 2 categories detected more or less than expected" |
| 4 | support running sum peak 0.62, table es 0.62 | accepted |
| 5 | support peak 0.60, table es 0.62 | `ValueError` |

## Acceptance criteria

- [ ] `tests/unit/viz/test_v14_enrichment.py`: examples 1–5; ADR-0025/0026 markings; no CI drawn
  on NES.
- [ ] Plotly twin registered (M9-14 mechanism).
- [ ] Coded against the branch schemas; re-checked when #181 merges.

## Questions for the lead

1. Panel (c) when nothing is significant: show the 20 smallest q (with the alpha line) or nothing?
   **Proposed default:** the 20 smallest, so the reader sees how far they are.
2. ORA's `effect` CI is NA in pass 1 (M8-04). **Proposed default:** points without bars, labelled
   "CI unestimated", until M8 hardening adds it.

## Lead decisions (2026-10-23)
1. **Q1 yes:** with nothing significant, panel (c) shows the 20 smallest q with the alpha line.
2. **Q2 yes:** ORA effect points without bars, labelled "CI unestimated", until M8 hardening.
3. **Stacking:** on m9/plot-support-analysis (M9-17; `gsea_running` per its tests a3224acf and
   decision 5). Merge M9-17 into this branch when it lands, then main.
4. **Test-author (be22ed42): pinned choices 1–7 and gaps accepted as proposed.**
   - Axes labels: `a`, `b` and `c*`.
   - Panel (b) is NES against −log10 q; panel (c) uses a log2 O/E axis.
   - O/E = 0 is drawn at the axis floor with a note, and the data keeps 0.
   - A support set missing from the table raises `ValueError`.
   - Both tables header-only gives an empty state with a generic message.
   - The implementer may add "TUs" to `checklist.UNITLESS` (M9-24 decision 7), and may move V14
     out of `LATER` in the twin helpers. That is a lead edit if an existing test must change:
     report it.
   - If a test's own assumption is wrong, report it; the fix goes to the test-author or lead.
5. **Review round (code review of 7410ec5b).** The implementer's choices 1, 2 and 4–7 are accepted. Choice 3 is reworked.
   - **B1, labels overlap:** compute panel (c)'s left edge from the measured width of (b)'s label column
     (renderer extents), so the two can't collide. Add a checklist case with maximum-length names in
     both panels.
   - **Test gaps (surviving mutants):**
     - the negative-set leading edge, which runs from peak to n;
     - the `_REL` boundary: one %.6g step is accepted, 1e-4 relative raises;
     - panel (c) drawn O/E-descending, top to bottom;
     - one `layer_missing` flag gives the generic EMPTY message.
   - **Nits taken:**
     - N1: ORA q = 0 noted in the footnote too.
     - N2: a q = 0 point uses marker `^`.
     - N3: leaders start at the label column with a short elbow, and the legend is pinned to `upper left`, outside the label rows.
     - N4: the "not calibrated" banner goes under the title.
     - N5: panel (c) is capped at 30 rows, plus "+N more in the TSV".
     - N6: ORA ties break on collection, term and direction.
     - N7: "no testable sets" when no GSEA row has a finite NES or q.
     - N8: twin parity. The TSV gains a `redundant` column; the twin keeps q = 0 at the cap and draws the alpha line, (c) in O/E order, and redundant rows small and grey.
     - N11: O/E = NaN rows get a "O/E undefined (E = 0)" note.
   - **Backlogged:** N9 (`x.max()`), N10 (`_draw_b` helper split), and "es" in UNITLESS.
6. **Test-author gaps (0d191489).**
   - **N5:** the data TSV keeps every significant ORA row; the figure draws the top 30. This exception to "panel = drawn rows" is written into the TSV's description.
   - **N2:** the `^` marker applies in both (b) and the (c) q strip.
   - **N3, N4:** read as the tests pin them.
   - **N8:** the twin's alpha line is at the default 0.05.
