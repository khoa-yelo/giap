# M9-06: V7 `obs_vs_exp` and V8 `cog_abundance`

Milestone M9 · **pass 1** · depends on: M9-01 · estimate ≈ 320 code lines (+ ≈ 300 test lines) ·
not statistical (both plots draw estimates and CIs from the tables)

## Spec sections implemented

- §10.2 V7: "forest plot of O − E (or O/E) with 95 % CI for replicons, feature classes, COG
  categories, region types and individual regions; q labels".
- §10.2 V8: "strips (or violins when n > 50) of log2 abundance per COG category, sorted by median,
  annotated with detected fraction".
- §7.3.5(a) (basic summary rows), §7.3.5(b) (region O − E, BH within type, `small_region`,
  plasmids `dosage_plus_regulation` with `p_value = NA`), ADR-0018 (the reported O − E is the raw
  effect; the test is studentised), §9.3 `detection_summary.tsv`, `regions_stats.tsv`,
  `genes.tsv`; §10.1 (titles state the finding, the example "Prophage genes: 12 % detected vs 31 %
  expected").

## V7

**Rows**, top to bottom in blocks:

1. `detection_summary` `section = basic`, `metric = o_minus_e`, blocks by `group_type` in the order
   replicon, feature_class, cog_category; `value`, `ci_low`, `ci_high`; percentages from the same
   group's `n_detected`, `expected` and `n_genes` rows.
2. `regions_stats` rows with `stat = o_minus_e` and `subset = all`, ordered by `region_type` (the
   `context.region_types` order) then `q_value` (NA last) then `region_id`: `effect`
   (= `o_minus_e`), `ci_low`, `ci_high`, `q_value`.

Each row: point and CI bar; a vertical line at 0; the label `group (n = …)`; at the right, `q =
…` for region rows (`%.2g`), "`dosage_plus_regulation` (untested)" for plasmid rows, and
"CI `unestimated`" when `ci_low`/`ci_high` are NA (`ci_few_units`, `ci_unestimated`). Regions
flagged `small_region` are drawn hollow. Region rows are coloured by region type (§10.1); basic
rows are black. The x axis is "observed − expected detected genes".

**Title** (index decision 5): the region with the smallest `q_value ≤ stats.alpha.context`:
"`<region_id>` (`<type>`): `<O/n> %` detected vs `<E/n> %` expected (q = `<q>`)", with O =
`observed`, E = `expected`, n = `n_eligible`; if none passes: "No region departs from detectability
at q ≤ `<alpha>`; `<k>` regions tested".

Data TSV `V07_obs_vs_exp.tsv`: `block` enum{replicon, feature_class, cog_category, region},
`group`, `region_type`, `n`, `observed`, `expected`, `o_minus_e`, `ci_low`, `ci_high`, `q_value`,
`flags`; in drawing order.

## V8

**Rows:** genes with `partial = false`; a gene appears under each letter of its `cog_category`
(M9-01 item 2); `cog_category` NA → the category "none" (drawn last, grey). Per category:
`n_genes`, `n_detected` (`detected = true`), and the values `abundance_log2` of detected genes with
`abundance_imputed = false` (index decision 13).

**Drawing:** categories on the x axis sorted by the median measured log2 abundance (descending; ties
by letter; categories with no measured value last); a strip of jittered points when the category
has ≤ 50 measured values, a violin (matplotlib `violinplot`, default Scott bandwidth) when > 50;
colour by COG super-category; above each category "`n_detected`/`n_genes`"; y "log2 abundance".
Jitter is deterministic: offsets from the gene's rank within the category, not random.

**Title:** "`<letter>` (`<description>`) has the highest median abundance (`<m>` log2, `<n>`
genes)"; with no measured abundance: "No measured abundance; detected fraction per COG category".
*Amended 2026-10-24 (lead, real-data fixes #322):* the claim is made only among categories with at
least 5 measured genes and says so ("… has the highest median abundance among categories with ≥ 5
genes (…)"); with none: "No COG category has 5 or more measured genes (most: `<letter>`, `<n>`
genes)". Gene counts use the singular for 1.

Data TSV `V08_cog_abundance.tsv`: `category`, `super_category`, `locus_tag`, `abundance_log2`,
`detected`, plus per-category summary rows with `locus_tag = NA` carrying `median_log2`,
`n_genes`, `n_detected`, `n_imputed`. Sorted by drawing order then `locus_tag`.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| `regions_stats` absent | V7 draws the basic blocks only; title from the basic rows: "`<group>`: `<O−E>` (95 % CI …)" for the largest \|O − E\| with a CI excluding 0, else "Detection matches detectability in every group with a CI" |
| every `basic` row has `n_detected = 0` (validation run) | rows drawn at 0 with "CI unestimated"; title "No region departs …" |
| `cog_category` NA everywhere (`layer_missing:cog`) | V7 has no COG block; V8 draws one "none" category and the title says "no COG annotation (layer_missing:cog)" |
| no abundance | V8 draws detected fractions only (bars), title as above |
| a category with 1 measured value | a single point, no violin |
| > 20,000 genes in one category | strips are never used above 50 values, so no element issue |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | regions `prophage_1`: n_eligible 25, observed 3, expected 7.75, o_minus_e −4.75, CI (−7.1, −2.2), q 0.003; `prophage_2`: q 0.4; alpha 0.05 | title "prophage_1 (prophage): 12 % detected vs 31 % expected (q = 0.003)"; two region rows in that order |
| 2 | plasmid row `p_value` NA, label `dosage_plus_regulation` | drawn, right label "dosage_plus_regulation (untested)" |
| 3 | basic replicon chr1: o_minus_e 2.1, CI (−1.0, 5.2) | row replicon/chr1 with those numbers |
| 4 | genes: J ×3 measured (20, 22, 24), C ×2 (18, 19), gene `KL` measured 21 | V8 order J (22), K (21), L (21), C (18.5); `KL` appears under K and L |
| 5 | J with 60 measured values | violin for J |
| 6 | V8 run twice | identical TSV (no random jitter) |

## Acceptance criteria

- [ ] `tests/unit/viz/test_v07_obs_vs_exp.py`: examples 1–3; ordering; `small_region` hollow; NA CI
  label; alpha read from `stats.alpha.context`.
- [ ] `tests/unit/viz/test_v08_cog_abundance.py`: examples 4–6; strip/violin switch at > 50; imputed
  values never drawn.
- [ ] Both render on `tests/data/golden/sim50k/m4/`.

## Questions for the lead

1. O − E or O/E (§10.2 "or")? **Proposed default:** O − E, the only one with a CI in the tables
   (index decision 15).
2. "Region types" (index question 3): no table carries a region-type O − E. **Proposed default:**
   omitted; backlog line for M4.
3. Should V7 add `category_summary.tsv` (M8) rows for KEGG and GO slim? **Proposed default:** no:
   V14(c) shows category detection; V7 keeps one COG estimate (the M4 one).
4. Violin bandwidth is a display choice; **Proposed default:** matplotlib's default (Scott), with
   no statistical claim attached.

## Lead decisions (2026-10-22, after the test author's d431aaf)

1. The test author's pinned readings of gaps 1–4, 6–10 are accepted:
   - region `n` = n_eligible;
   - k counts the rows with a q value;
   - a tie goes to the first region in drawing order;
   - an estimate is drawn at 0 when E = 0;
   - the COG "none" row is kept;
   - V08 has one row per non-partial gene and letter, with `abundance_log2` NA unless drawn;
   - "none" goes last.
2. **Gap 5.** With no `regions_stats`, the basic-only title rule applies, even when every row has n_detected = 0.
3. **Gap 9.** The `locus_tag = NA` summary row comes first in its category. `super_category` for "none" is `none`.
4. **Gap 11.** "no COG annotation" appears when every COG letter is NA, whether or not the `layer_missing:cog` flag is present. The flag is named in brackets only when it is present.
5. **Gap 12.** `<m>` is the median log2 abundance, formatted `%.1f`. `<n>` is the number of measured values drawn in that category.
6. **Gaps 13 and 14.** Use the framework's default empty state. V07 is single width up to 20 rows and double above; V08 is double width.
7. **Width (after implementation).** V07 also goes double width whenever an "(untested)" label is drawn, because at single width the label squeezes the forest. Accepted (it amends decision 6).
