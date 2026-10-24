# M9-11: V11 `operon_pairs`: do operon partners look alike?

Milestone M9 · **pass 1** · depends on: M9-01 · estimate ≈ 260 code lines (+ ≈ 260 test lines) ·
not statistical (histograms of table values; lfdr is M6's)

## Spec sections implemented

- §10.2 V11: "(a) Δ densities for OP, NO, CV and SW; (b) local fdr curve; (c) partner scatter
  (upstream vs. downstream log2) coloured by class and flag".
- §7.4 (pair classes OP, NO, CV, DV, SW; DV "reported separately"; decoys CV and SW), §7.4 "Local
  fdr for ranking: lfdr(s) = min(1, π₀ f₀(s) / f_NO(s)) … f₀ the decoy density of |Δ|", §7.4.5
  (flags, `stats.alpha.operon_flags` 0.10), §9.3 `pairs.tsv`, `operon_flags.tsv`.

## Panels (double width, 1 × 3)

- **(a) Δ densities:** pairs with finite `delta_log2` (both partners measured), per class in
  {OP, NO, CV, SW}: density-normalised histograms of |Δ| = |`delta_log2`| on shared bins (index
  decision 16): 30 equal-width bins over [0, max |Δ| across the four classes]; step outlines in
  the Okabe–Ito order OP `#0072B2`, NO `#E69F00`, CV `#009E73`, SW `#CC79A7`; legend with each
  class's n. DV pairs are not drawn (§7.4 reports them separately); their n is in the data TSV.
- **(b) Local fdr:** rows with finite `lfdr` (the `missed_operon_decoy` test), x = |`delta_log2`|,
  y = `lfdr`, points joined in x order; a dashed line at `stats.alpha.operon_flags` labelled "flag
  threshold (q)". Note in the panel: "lfdr ranks pairs; flags use q".
- **(c) Partner scatter:** same-strand pairs (OP, NO) with both abundances measured: x = upstream
  partner's `abundance_log2`, y = downstream's. Upstream = `gene_a` on `+`, `gene_b` on `−`
  (transcription order). Colour by class; pairs in `operon_flags.tsv` with `q_value ≤
  stats.alpha.operon_flags` drawn as larger outlined markers, marker shape by flag (7 shapes for the
  7 §7.4.5 values); identity line y = x.

## Data TSV `V11_operon_pairs.tsv`

`panel` enum{a, b, c}, `pair_id`, `class`, `abs_delta_log2`, `lfdr`, `upstream_log2`,
`downstream_log2`, `flag`, `q_value`, `bin`, `density`. Panel (a) has one row per class × bin
(`pair_id` NA) and per-class count rows (`bin` NA, `density` NA, `abs_delta_log2` NA, n in a
`n` column); (b) and (c) one row per pair. Sorted by panel, class (OP, NO, CV, SW, DV), bin,
pair_id.

## Title

"Median |Δ| OP `<a>` vs NO `<b>` log2; `<k>` pairs flagged at q ≤ `<alpha>`" (medians over panel
(a)'s pairs; descriptive, no comparison claimed). Too few measured pairs (< 1 OP or < 1 NO):
"`<k>` pairs flagged at q ≤ `<alpha>`; too few measured pairs for densities".

## Degenerate cases

| Situation | Behaviour |
|---|---|
| validation run: `operon_flags.tsv` header-only, pairs unmeasured | (a), (c) empty-state text; (b) "no lfdr values"; title "0 pairs flagged …; too few measured pairs" |
| no CV or SW pairs (single-strand genome) | those classes absent from (a), n 0 in TSV |
| `max |Δ| = 0` | bins over [0, 1] |
| `operon_flags.tsv` absent | (c) without flag markers |
| one class with a single pair | a single-bin histogram |
| a pair crossing the origin | no special case (pairs carry no coordinates) |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | OP pairs Δ −0.2, 0.4, 0.6; NO pairs 1.0, −1.6 | medians OP 0.4, NO 1.3; title "Median \|Δ\| OP 0.4 vs NO 1.3 log2; 0 pairs flagged at q ≤ 0.1" |
| 2 | PAIR_1 `+`, a = 20, b = 22 | (c) upstream 20, downstream 22 |
| 3 | PAIR_2 `−`, a = 20, b = 22 | (c) upstream 22, downstream 20 |
| 4 | PAIR_3 flag `unexplained_discordance` q 0.04 | (c) row with flag, outlined marker |
| 5 | PAIR_4 flag q 0.2 | (c) not outlined |
| 6 | DV pairs ×3 | not drawn; TSV count row DV n 3 |

## Acceptance criteria

- [ ] `tests/unit/viz/test_v11_operon_pairs.py`: examples 1–6; shared bins; alpha line from
  `stats.alpha.operon_flags`; upstream rule by strand.
- [ ] Renders on the validation run (empty states) and on an M6 test fixture with flags.

## Questions for the lead

1. Densities (index decision 16): histograms rather than KDEs, so the figure makes no bandwidth
   choice. The lfdr itself uses a Silverman KDE (§7.4). **Proposed default:** histograms.
2. Panel (c) restricted to same-strand pairs (CV and SW have no upstream/downstream). **Proposed
   default:** yes.

## Lead decisions (2026-10-22, after the test author's 16e368b)

1. **Upstream.** `gene_a` is upstream on both strands: pairs.tsv is already in transcription order (M6-04), so V11 never swaps. This overrides the spec's `-` strand rule and its example 3, as pinned by `test_upstream_is_gene_a_on_both_strands`.
2. **Validation run.** It carries finite iBAQ `delta_log2`, so panels (a) and (c) draw. The slow test pins its real counts. The spec's "unmeasured" degenerate row is covered by the synthetic test.
3. **Imputed abundances.** Panels (a) and (c) drop pairs flagged `abundance_imputed` (decision 13: measured values only). The count rows' `n` counts the pairs kept.
4. **k.** k is the number of operon_flags.tsv rows with q ≤ α. This is a documented exception to index decision 5: the title reads operon_flags, not only the data TSV.
5. **Row details.**
   - A pair with an operon_flags row carries both `flag` and `q_value`, whatever q is.
   - There are always five count rows, zeros included.
   - `n` on bin rows is NA.
6. **The rest.**
   - pair_id rows keep pairs.tsv order.
   - Bins are 0-based.
   - Title medians use `%.2g`.
   - A header-only pairs.tsv gives "No gene pairs in this run".
   - Axis labels for (a) and (b) are the implementer's choice, stated in the report.

**Amendment (lead, 2026-10-24, viz batch 2 on the M9-25 fixtures).** When `operon_flags.tsv` has
more rows (`n`) than pass alpha (`k`), the title reads "`<k>` of `<n>` candidate flags at q ≤
`<alpha>`" instead of "`<k>` pairs flagged …". This keeps a 0-flag title true while candidate
rows exist. With `n = k`, the wording above is unchanged.
