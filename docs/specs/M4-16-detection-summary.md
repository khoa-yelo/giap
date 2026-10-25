# M4-16: Basic detection summary (O − E) and the `detection_summary.tsv` writer

Milestone M4 · track R (context statistics) · depends on: M4-04 (TU-block bootstrap); consumes the
row producers of M4-12 (`tier`) and M4-13 (`detectability`) through their interfaces · estimate ≈ 330
code lines (+ ≈ 300 test lines) · **can start now** (tests use small hand-built gene frames with the
`genes.tsv` columns)

## Spec sections implemented

- §7.3.5(a) quoted: "Basic detection summary (`detection_summary.tsv`, section `basic`). For each
  replicon, feature class (`cds`, `sorf`, `pseudogene`, `hypothetical`) and COG category, report:
  n_genes, n_detected, % detected; E = Σ p_detect; O − E, with a 95 % CI from a TU-block bootstrap
  (resample the group's TUs with replacement); O/E."
- §9.3 `detection_summary.tsv`: "Long format: `section` enum{tier, basic, detectability, qc},
  `group_type`, `group`, `metric`, `value`, `ci_low`, `ci_high`, `n`, `flags`. Covers: the tier and its
  FDR curve points; the basic O − E rows; detectability calibration bins, Brier score, slope and
  intercept; QC rows: `id_mixups`, `genome_inconsistent_rate`, contaminant counts."
- §8.2 "Adjust for detectability … never against raw counts." "Resample TUs, not genes."
- §5.2 `stats.bootstrap: 1000`.

## Public interface

`src/giap/context/summary.py` (pure)
```python
def basic_rows(genes: pd.DataFrame, B: int, rng: np.random.Generator) -> pd.DataFrame: ...
def qc_rows(peptides, peptide_hits, unmapped) -> pd.DataFrame: ...
def detection_summary(parts: Sequence[pd.DataFrame]) -> pd.DataFrame: ...   # concat, validate, sort
SUMMARY_COLUMNS = ("section", "group_type", "group", "metric", "value", "ci_low", "ci_high", "n", "flags")
```
`src/giap/io/summary_tsv.py`: `write_detection_summary(path, df)`;
`src/giap/schemas/tables/detection_summary.json`.

- `basic_rows` uses eligible genes only (`detection_status != partial_excluded`). Groups:
  `group_type = "replicon"` (each contig), `"feature_class"` (`cds`, `sorf`, `pseudogene`, and
  `hypothetical` = rows with `hypothetical = true`, any feature type), `"cog_category"` (each letter of
  `cog_category`; a gene with `KT` counts in K and in T; genes without a category in `none`).
  Metrics per group: `n_genes`, `n_detected`, `pct_detected`, `expected` (Σ p_detect), `o_minus_e`
  (with `ci_low`/`ci_high` from `tu_block_bootstrap` over the group's TUs; only the group's genes
  inside each drawn TU count), `o_over_e` (NA when E = 0). `n` = n_genes on every row.
- Detected = `detected` column (primary analysis: single-peptide genes are not detected).
- RNG: `giap.rng.for_step("context_stats", seed)`; groups processed in sorted order so draws are
  reproducible.
- `qc_rows` (`section = qc`, `group_type = "run"`, `group = "all"`): `id_mixups` = distinct peptides
  with C11 flag `engine_protein_mismatch`; `genome_inconsistent_rate` = among C9 target peptides with
  `search_db = same`, the share that are not placed exactly (absent from C11, or only tolerant);
  `contaminant_peptides`, `contaminant_proteins` = C13 rows with `label = contaminant` by `kind`;
  `unmapped_peptides`, `unmapped_proteins`.
- `detection_summary` concatenates tier (M4-12), basic, detectability (M4-13) and qc rows; sorts by
  section order (tier, basic, detectability, qc), then group_type, group, metric.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | 4 eligible genes on `chr1`, d = [1, 0, 1, 1], p = [0.9, 0.2, 0.5, 0.4] | `n_genes 4`, `n_detected 3`, `pct_detected 75`, `expected 2.0`, `o_minus_e 1.0`, `o_over_e 1.5` |
| 2 | same 4 genes in TUs {A: g1, g2; B: g3; C: g4}, B = 1,000 | TU sums of d − p: A = 0.1 − 0.2 = −0.1, B = 0.5, C = 0.6 (total 1.0); every replicate is a sum of 3 draws from {−0.1, 0.5, 0.6} → CI ⊆ [−0.3, 1.8] |
| 3 | gene with `cog_category = "KT"` | counted in K and T rows |
| 4 | group with E = 0 (all p_detect 0) | `o_over_e = NA` |
| 5 | partial genes present | excluded from every basic row |
| 6 | C11 with 3 distinct peptides flagged `engine_protein_mismatch` | `id_mixups = 3` |
| 7 | 100 `search_db = same` target peptides, 97 placed exactly, 2 tolerant only, 1 absent | `genome_inconsistent_rate = 0.03` |
| 8 | full fixture written with `write_detection_summary` | validates against `detection_summary.json`; sections in order tier, basic, detectability, qc |

## Edge cases

- A group with one TU: CI degenerate (equal to the estimate) — allowed, reported.
- Empty group lists (no COG annotation at all) → only `none` rows.
- Missing C13 → contaminant/unmapped metrics `NA`.
- `pct_detected` is a percentage (0–100), named `pct_` per §9.2 units-in-names.

## Acceptance criteria

- [ ] `tests/unit/context/test_summary.py`: examples 1–7; `tests/integration/test_detection_summary.py`:
  example 8 on the M4-07/08/09 fixture.
- [ ] Property: relabelling genes outside a group leaves that group's rows unchanged (§7.3.6 property
  "O − E is unchanged by relabelling genes outside the region" applied to groups).
- [ ] Deterministic across runs with the same seed.

## Methods note (statistician sign-off)

`docs/methods/detection-summary.md`: O − E against cross-fitted Σ p_detect, TU-block CI for group sums,
how multi-letter COG categories are counted, the QC metrics' definitions.

## Questions for the lead

1. `genome_inconsistent_rate` is not defined in the spec. **Proposed default:** as above (share of
   same-database target peptides not placed exactly), computed only when `search_db = same` rows exist.
2. Multi-letter COG categories counted in every letter? **Proposed default:** yes (COG convention);
   documented in the methods note.

## Lead decisions (after the test-author's pass)
1. The test-author's resolutions are accepted (`n_tu` row per group; `ci_few_units` flag on `o_minus_e`
   with < 10 TUs; missing `n` → NA, missing `flags` → empty; ValueError on an unknown section or a missing
   core column; `genome_inconsistent_rate` NA without same-database peptides).
2. The schema is `src/giap/schemas/tables/detection_summary.schema.json` (same convention as the M3-16
   table schemas).
3. Stand-ins for unmerged producers (p_detect, COG, tier/detectability columns) are replaced by the real
   M4-12/M4-13 outputs when those merge (backlog, M4-20).

## Lead decisions after the statistician's review (PR #87, changes requested)
4. `ci_few_units` is also set when min(E, n − E) < 10 (E = expected detections in the group), in addition
   to < 10 TUs; a gated attack null for a low-E group (e.g. sORFs) goes to tests/calibration/attacks/.
5. `genome_inconsistent_rate` follows §7.5.1: computed over `search_db = other` peptides (NA for `same`),
   contaminant peptides excluded from the denominator.
6. `id_mixups` (§7.5.1: same-database ALT sites that are not start/stop sites) needs the variant module
   (M7): reported as NA with flag `requires_variants` until then; the engine/protein mismatch count is
   reported as its own metric `engine_protein_mismatch_rate` (no ADR — the spec metric is deferred, not
   changed). Backlog for M7.
7. Fail closed: NaN `p_detect` on an eligible gene, or an eligible gene without `tu_id`, raises
   ValueError naming the genes.
8. COG values `-` and blank map to `none`; `detected` must be boolean dtype (strings raise).
9. Tests kill the surviving mutants (group order; flag only on `o_minus_e`).

## Amendment (M7-16b, 2026-10-24)
`detection_summary.tsv` gains the `qc` placement-kind rows (`placement_kind`, `tolerant_diff_type`,
`placement_flag`, `unmapped_reason`; metric `n_peptides`) from a separate
`placement_kind_rows` part. `qc_rows` is unchanged. See
`docs/specs/M7-16-stop-codon-readthrough.md` rule 17.
