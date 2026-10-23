# M9-05: V5 `mapping_sankey` and V6 `mapping_fdr`

Milestone M9 · **pass 1** · depends on: M9-01 · estimate ≈ 340 code lines (+ ≈ 300 test lines) ·
not statistical (V6 draws the M3 estimates; it computes none)

## Spec sections implemented

- §10.2 V5: "input class (peptide or protein) → route (exact, miniprot, tblastn, tolerant,
  unmapped) → outcome (placed, contaminant, unexplained, rejected by FDR)"; custom matplotlib
  alluvial (Plotly Sankey in M9-15).
- §10.2 V6: "(a) targets vs. decoys by length class × differences; (b) stratum FDR with the 0.01
  line; (c) de novo precision by score bin with CI (if de novo)".
- §9.3 `peptide_hits.tsv`, `protein_hits.tsv`, `unmapped.tsv`, `mapping_fdr.tsv` (sections
  `placement` and `denovo_precision`, per-m totals with `length_class = all`), §7.2 (G3 strata),
  §5.2 `stats.alpha.placement`, `mapping.length_classes`.

## V5: inputs and the flow table (index decision 14)

Units are distinct sequences: peptides by sequence, proteins by `protein_id`.

| Source rows | Input class | Route | Outcome |
|---|---|---|---|
| `peptide_hits`, distinct `peptide`, any row with `route = exact` | peptide | exact | placed |
| `peptide_hits`, distinct `peptide`, all rows `route = tolerant` | peptide | tolerant | placed |
| `protein_hits`, distinct `protein_id`, best `rank` row's `route` | protein | miniprot / tblastn | placed |
| `unmapped`, `label = contaminant` | `kind` | unmapped | contaminant |
| `unmapped`, `reason = rejected_placement_fdr` | `kind` | unmapped | rejected by FDR |
| `unmapped`, other reasons (`no_hit`, `below_identity`, `too_short`) | `kind` | unmapped | unexplained |

A sequence present in both `peptide_hits` and `unmapped` is counted once, as placed, and the data
TSV records `n_conflicts`. Node order is fixed: classes (peptide, protein), routes (exact, tolerant,
miniprot, tblastn, unmapped), outcomes (placed, contaminant, rejected by FDR, unexplained). Node
colours: placed `#0072B2`, others from Okabe–Ito in that order; ribbons take the outcome colour.
Routes with zero flow are not drawn but stay in the data TSV with `n = 0`.

Data TSV `V05_mapping_sankey.tsv`: `input_class`, `route`, `outcome`, `n`, sorted by the node
order; a trailing row `input_class = all`, `route = all`, `outcome = conflicts`, `n = n_conflicts`.

Title: "`<x> %` of `<N>` input peptides placed in the genome" (+ "; `<y> %` of `<M>` proteins"
when proteins exist).

## V6: panels (double width)

- **(a)** grouped bars per `length_class` (order of `mapping.length_classes`), one group per
  `n_diff` (0, 1, 2), bars `targets` (`#0072B2`) and `decoys` (`#D55E00`), log scale with 0 drawn
  at the axis floor; `section = placement` rows with `length_class ≠ all`.
- **(b)** `fdr_stratum` of the per-m total rows (`length_class = all`) against `n_diff`, with a
  dashed horizontal line at `stats.alpha.placement` (0.01 by default, read from config, never a
  literal) and `n_accepted` printed above each point; NA strata drawn as an open marker at the
  axis floor labelled "NA (no targets)".
- **(c)** only when `section = denovo_precision` rows exist: `precision` with `ci_low`–`ci_high`
  error bars per `score_bin`; rows flagged `c_assumed` drawn hollow, with a footnote "c assumed
  0.99 (no reference peptides)".

Data TSV `V06_mapping_fdr.tsv`: `panel` enum{a, b, c}, `length_class`, `n_diff`, `metric`, `value`,
`ci_low`, `ci_high`, `flags`. Sorted by panel, length class (config order), n_diff, metric.

Title: from panel (b): "Placement FDR ≤ `<alpha>` in strata `<list>`; `<n>` placements accepted"
or, when a stratum is NA, "… stratum 2 unestimated (no targets)".

## Degenerate cases

| Situation | Behaviour |
|---|---|
| protein-only run (no peptides) | V5 has only the protein class; V6 panels (a), (b) empty-state text ("no peptide placements") |
| no `protein_hits.tsv` rows | protein class absent from the figure; zero rows kept in the TSV |
| every `mapping_fdr` target count 0 (the 2026-10-22 validation run's first rows) | (a) bars at the floor; (b) NA markers; title "unestimated" |
| no de novo rows | panel (c) not drawn, no rows |
| a peptide both placed and unmapped | counted placed; `conflicts` row > 0 |
| > 20,000 sequences | no change: V5 draws ≤ 2 × 5 × 4 ribbons |

## Examples

| # | Input | Expected V05 rows |
|---|---|---|
| 1 | hits: `PEPA` exact (2 loci), `PEPB` tolerant; unmapped peptides `PEPC` no_hit, `PEPD` contaminant, `PEPE` rejected_placement_fdr | peptide/exact/placed 1; peptide/tolerant/placed 1; peptide/unmapped/unexplained 1; peptide/unmapped/contaminant 1; peptide/unmapped/rejected by FDR 1; title "40 % of 5 input peptides placed in the genome" |
| 2 | protein hits `P1` rank 1 miniprot, rank 2 tblastn; unmapped protein `P2` below_identity | protein/miniprot/placed 1; protein/unmapped/unexplained 1 |
| 3 | `PEPA` in hits and in unmapped | placed 1, `conflicts` 1 |

| # | Input `mapping_fdr` | Expected V06 |
|---|---|---|
| 4 | all/0: targets 900, decoys 3, fdr 0.0033, accepted 880; all/1: 120, 2, 0.025, 40; all/2: 0, 0, NA, 0; alpha 0.01 | (b) three rows; title "Placement FDR ≤ 0.01 in strata 0; 920 placements accepted; stratum 2 unestimated (no targets)" |
| 5 | `length_class` 7, n_diff 0: targets 5, decoys 0 | (a) rows `targets` 5, `decoys` 0 |
| 6 | denovo rows bins 1–3, precision 0.8/0.9/0.95, CI given, flag `c_assumed` | (c) three rows; footnote present |

## Acceptance criteria

- [ ] `tests/unit/viz/test_v05_sankey.py`: examples 1–3; node order; zero routes kept in TSV.
- [ ] `tests/unit/viz/test_v06_mapping_fdr.py`: examples 4–6; the alpha line equals
  `cfg.stats.alpha.placement` (test with `--set stats.alpha.placement=0.05`).
- [ ] Both render on the M3 golden tables `tests/data/golden/sim50k/mapping/`.

## Questions for the lead

1. The "stratum FDR" title claims "≤ alpha" only for strata whose `fdr_stratum` is finite and
   ≤ alpha; it says nothing about the others. **Proposed default:** as written above.
2. "rejected by FDR" cannot be split by route (`unmapped.tsv` has no route column). **Proposed
   default:** route `unmapped` for every unmapped row (index decision 14); adding a route column
   to `unmapped.tsv` would be an M3 schema change and is not proposed.

## Lead decisions (2026-10-22, after the test author's a6c601c)

1. The V6 data TSV gets `score_bin` as its 4th column (`panel, length_class, n_diff, score_bin, metric, value, ci_low, ci_high, flags`).
2. Panel (b) has six rows per run: `fdr_stratum` and `n_accepted` for each n_diff. The title reads its count from them (index decision 5).
3. **`decoy_type`.** Panels (a) and (b) use the `reverse` rows (database peptides). The `shuffle` rows (de novo) are used only when the reverse rows have no targets. A run with targets under both decoy types is left for the hardening pass (backlog).
4. **Title wording.** The title says "<n> FDR-tested placements accepted". Exact `same_db` placements are not FDR-tested, so "placements accepted" alone would understate the placed peptides. The test strings were updated by the lead.
5. The test author's readings of gaps 5–8 are accepted as tested:
   - the edge-case titles;
   - the V5 empty state and its 10-row grid;
   - the matching keys: `sequence` for peptides, `sequence_id` for proteins;
   - `contaminant_id` counts as contaminant;
   - placed is `#0072B2`, and the other outcomes take the next Okabe–Ito colours in palette order.
