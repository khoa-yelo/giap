# M4-18: `regions_stats.tsv`: plasmid replicon-level rows, prophage subsets, writer, `context_stats` step

Milestone M4 · track R · depends on: M4-15 (gene table), M4-16 (summary writer), M4-17 (region
tests) · estimate ≈ 330 code lines (+ ≈ 300 test lines) · **can start now** on fixtures once its
dependencies merge

## Spec sections implemented

- §7.3.5(b) quoted: "Plasmids are replicon-level effects. Report O − E and the HL shift of the plasmid
  against the chromosome, with TU-block bootstrap CIs. Set `p_value = NA` and `null_model = none`, and
  label the result `dosage_plus_regulation`. With `--replicon-depth`, also report the shift minus
  log2(plasmid depth / chromosome depth)." "Prophage function subsets. Per prophage, report detected /
  total for each subset of §7.1.7 (rows with `subset` ≠ `all`). These rows are descriptive and untested
  in version 1."
- §9.3 `tables/regions_stats.tsv`: "region: `region_id`, `region_type`, `subset` (`all` or a prophage
  function subset), `contig`, `start`, `end`, `n_genes`, `n_eligible`; detection: `observed`,
  `expected`, `o_minus_e`, `o_over_e`; abundance: `hl_shift_log2`; test columns (§9.2): one row per
  statistic (`stat` enum{o_minus_e, hl_shift}); `label` (e.g., `dosage_plus_regulation`), `flags`."
- §9.2 test columns: `stat`, `effect`, `ci_low`, `ci_high`, `p_value`, `q_value`, `method`,
  `null_model`, `n_perm`.
- §8.5 attack 4: "High-copy plasmids → region tests, reported as replicon-level effects only."
- §5.3 "Replicon depth (`--replicon-depth`). TSV `contig`, `mean_depth`."
- §6.2 step `context_stats` ("gene table → `regions_stats`, `autocorrelation`, `islands`,
  `detection_summary`, ori–ter"); M4 implements the regions and detection-summary parts, M5 adds the
  rest.

## Public interface

`src/giap/context/regions_table.py` (pure)
```python
REGIONS_STATS_COLUMNS: tuple[str, ...]      # §9.3 order + the §9.2 test columns
def plasmid_rows(genes, replicons, depth: pd.DataFrame | None, cfg, rng) -> pd.DataFrame: ...
def subset_rows(genes, regions, region_genes) -> pd.DataFrame: ...
def regions_stats(tests: pd.DataFrame, plasmids: pd.DataFrame, subsets: pd.DataFrame,
                  regions: pd.DataFrame) -> pd.DataFrame: ...
```
`src/giap/io/depth.py`: `read_replicon_depth(path) -> tuple[pd.DataFrame | None, list[Issue]]`.
`src/giap/io/regions_tsv.py`: `write_regions_stats(path, df)` (1-based coordinates).
`src/giap/schemas/tables/regions_stats.json`.
`src/giap/steps/context_stats.py`: `run_context_stats(work_dir, cfg, replicon_depth: Path | None)`
→ `tables/regions_stats.tsv`, `tables/detection_summary.tsv` (M4 sections), ready for M5 to extend.

- **Plasmid rows** (one pair per C7 region with `region_type = plasmid`; M2-19 makes one per Bakta
  plasmid replicon or geNOMAD plasmid contig):
  `o_minus_e` over the plasmid's eligible genes; `hl_shift` = HL(plasmid detected log2, chromosome
  detected log2) where "chromosome" = all replicons with `kind = chromosome`; CIs by TU-block bootstrap
  (HL: strata plasmid / chromosome); `p_value = q_value = NA`, `null_model = none`, `n_perm = 0`,
  `method = replicon_contrast`, `label = dosage_plus_regulation`, flag `dosage_plus_regulation`.
  With depth: an extra `hl_shift` row with `label = dosage_adjusted`, `effect` = HL − log2(plasmid
  depth / chromosome depth), CI shifted by the same constant (depth treated as known).
  Chromosome depth = length-weighted mean depth of chromosome replicons.
- **Prophage subset rows**: per prophage and subset present: `subset` = name, `stat = o_minus_e`,
  `observed` = detected count, `n_genes`, `n_eligible`, `expected` = Σ p_detect, `o_minus_e`,
  `o_over_e`; test columns NA except `effect` (= o_minus_e); `null_model = none`, flag `descriptive`.
- **Test rows** (M4-17) get `subset = all`; detection columns are filled on both stat rows of a region
  (`observed`, `expected`, `o_minus_e`, `o_over_e`, `hl_shift_log2`) so each row is self-contained.
- Region coordinates from C7 (AMR: NA). Sort: region_type order of `context.region_types`, then
  contig order, start, region_id, subset (`all` first), stat.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | chromosome + plasmid `p1` (60 genes); plasmid detected log2 median 4 units above chromosome | two `plasmid_1` rows: `stat o_minus_e` and `hl_shift` (≈ 4), `p_value NA`, `null_model none`, `label dosage_plus_regulation` |
| 2 | example 1 + depth `chr1 50`, `p1 400` | extra row `label dosage_adjusted`, `effect` = HL − 3 |
| 3 | prophage_1 with subsets structural 16 genes (8 detected), lysis 4 (0 detected) | rows `subset structural` `observed 8`, `subset lysis` `observed 0`, flag `descriptive`, `p_value NA` |
| 4 | M4-17 output for `bgc_1` | 2 rows, `subset all`, `null_model circular_shift`, `n_perm` = K |
| 5 | genome with no plasmid | no `plasmid` rows |
| 6 | plasmid but no chromosome (`kind` all unknown) | plasmid rows with `hl_shift` NA and flag `no_reference_replicon` (see Q2) |
| 7 | depth file missing a contig | warning issue; adjusted row omitted |
| 8 | fixture end to end through `run_context_stats` | `regions_stats.tsv` validates against its schema; header = `REGIONS_STATS_COLUMNS` |
| 9 | same seed twice | byte-identical TSV |

## Edge cases

- A plasmid with no eligible gene → O − E 0, HL NA.
- Depth ≤ 0 → error issue from the reader (log2 undefined).
- `regions_stats.tsv` with no regions at all: header only.
- Coordinates of origin-crossing regions keep `end > length` (as §6.5) in the 1-based TSV.

## Acceptance criteria

- [ ] `tests/unit/context/test_regions_table.py`: examples 1–7, 9; `tests/integration/test_context_stats_step.py`:
  example 8.
- [ ] Attack null 4 (statistician): with `plasmid_copy_log2 = log2(20)`, no plasmid row ever has a
  p-value; chromosome region tests at level (shared with M4-17's attack test).
- [ ] No row has `p_value` without `null_model ≠ none`; no row has `q_value` without `p_value` (schema
  rule + test).

## Methods note (statistician sign-off)

`docs/methods/region-tests.md` §"Replicon-level effects" (why plasmids get no p-value; the dosage
adjustment as a known offset) and §"Descriptive rows".

## Questions for the lead

1. Plasmid rows: also an O − E row against the chromosome rate, or only the plasmid's own O − E?
   **Proposed default:** the plasmid's own O − E (E already uses p_detect); HL is the contrast.
2. No `chromosome` replicon (draft genome): **proposed default:** plasmid HL = NA with a new flag
   `no_reference_replicon` (schema vocabulary extension, lead to approve), or use all non-plasmid
   replicons as reference — the latter is simpler; choose one.
3. `--replicon-depth` plumbing in `giap run` does not exist yet. **Proposed default:** M4 exposes the
   step parameter and reader; CLI wiring lands with the `run` command.
