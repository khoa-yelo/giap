# M5-09: `autocorrelation.tsv`, `islands.tsv`, gene context columns and the `context_stats` step extension

Milestone M5 · track T (tables and step) · depends on: M5-05, M5-06, M5-08; **M4-18** (the
`context_stats` step and `run_context_stats`, not merged yet) and **M4-15** (`work/genetable/genes.parquet`,
`write_genes_tsv`, not merged yet); the integration test uses the **M4-20** fixture builder · estimate ≈ 320 code lines (+ ≈ 300 test lines) · **not statistical**
(plumbing; the statistician checks that no row loses its error label) · no tools, no Slurm

## Spec sections implemented

- §6.2 step `context_stats` "gene table → `regions_stats`, `autocorrelation`, `islands`,
  `detection_summary`, ori–ter".
- §9.3 quoted: "`tables/autocorrelation.tsv`: `contig`, `statistic` enum{morans_i, join_count_bb,
  oriter_slope}, `lag`, `null_model` enum{gene_label, tu_block, bootstrap}, `value`, `null_mean`,
  `null_lo`, `null_hi`, plus the test columns." "`tables/islands.tsv`: `island_id`, `direction`
  enum{high, low}, `contig`, `start`, `end`, `tu_first`, `tu_last`, `n_tus`, `n_genes`, `genes` (list),
  `mean_abundance_log2`, `detected_frac`, `scan_s`, `block_len`, plus the test columns (`p_value` = FWER,
  `q_value = NA`)." genes.tsv context group: `island_id`, `replichore`, `replichore_position`.
- §9.2 table conventions and test columns; §9.4 "Region and context" flags (`descriptive`, `no_oric`,
  `uncalibrated`) — "extend only through the schema"; §9.5 Frictionless resources; §6.5 coordinates and
  sort; §6.6 / R7 determinism; R4 (every claim carries FWER/q or `unestimated`).

## Public interface

- `src/giap/io/context_tsv.py`: `write_autocorrelation(path, df)`, `write_islands(path, df)` and readers
  (inverse; islands 0-based → 1-based `start + 1`, `end` unchanged, `end > length` kept for origin
  crossers; lists `;`-joined; through M1-04 `write_tsv`, `%.6g`).
- `src/giap/schemas/tables/autocorrelation.schema.json`, `islands.schema.json` (same layout as
  `detection_summary.schema.json`: name, type, unit in the name, description, enums, `coordinates`
  `1-based-closed` for islands, none for autocorrelation).
- `src/giap/context/vocab.py`: `CONTEXT_FLAGS` = §9.4 Region-and-context flags + M4 additions
  (`few_shifts`, `split_contigs`, `hl_variance_nonpositive`, `no_reference_replicon`, `ci_few_units`) +
  **M5 additions `not_estimable`, `few_values`** (Q1); both schemas constrain `flags` items to it.
- `src/giap/steps/context_stats.py` (M4-18's step) gains inputs `tus` (C5, required) and `features`
  (C15, optional) and outputs `autocorrelation` → `tables/autocorrelation.tsv`, `islands` →
  `tables/islands.tsv`, `gene_context` → `work/context/gene_context.parquet` (`locus_tag`, `island_id`,
  `replichore`, `replichore_position`; one row per gene of the gene table).
- genes.tsv: the writer that produces `tables/genes.tsv` (M4-15 `write_genes_tsv`, called after
  `context_stats`) left-joins `gene_context` on `locus_tag`; without it the three columns stay NA (Q2).

**Row assembly.** autocorrelation = M5-05 rows ∪ M5-08 rows, sorted by contig (C1 order), `statistic`
(enum order), `null_model` (enum order), `lag` (NA last), `stat`. islands = M5-06 frame sorted by
contig, start, end, `island_id`. Header order exactly as §9.3 followed by the §9.2 test columns `stat,
effect, ci_low, ci_high, p_value, q_value, method, null_model, n_perm` (with `null_model` appearing once,
in its §9.3 position, for autocorrelation) and `flags` last.

**Schema rules (validated in tests, also as Frictionless constraints where expressible):**
- autocorrelation: `p_value` NA ⇒ `q_value` NA; `statistic = oriter_slope` ⇔ `null_model = bootstrap` ⇔
  `p_value` NA with flag `descriptive`; `lag` NA iff `statistic = oriter_slope`; `p_value`, `q_value` ∈
  [0, 1].
- islands: `q_value` always NA; `p_value` ∈ (0, `stats.alpha.context`]; `method = scan_max_stat`;
  `n_tus` ≥ `context.islands.min_tus`; `island_id` matches `ISL_[1-9][0-9]*`.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | 1,000-gene fixture (M4-20 builder) with an `oric_frac` and a planted island | both TSVs written; validate with `frictionless validate`; headers equal the module constants |
| 2 | island on TUs with 0-based start 1,000, end 9,500 | TSV `start 1001`, `end 9500` |
| 3 | island crossing the origin of a 50,000-bp replicon (0-based 49,000 – 51,200) | TSV `start 49001`, `end 51200` |
| 4 | no island | `islands.tsv` is the header line only; every genes.tsv `island_id` NA |
| 5 | no `features.parquet` | autocorrelation has the `no_oric` row for each circular chromosome; run succeeds; manifest/log records `layer_missing:features` (existing run-flag pattern) |
| 6 | gene in `ISL_2` on the right replichore at position 0.31 | genes.tsv `island_id ISL_2`, `replichore right`, `replichore_position 0.31` |
| 7 | same seed twice | byte-identical `autocorrelation.tsv`, `islands.tsv`, `gene_context.parquet` (TE7) |
| 8 | adding M5 to `context_stats` | `regions_stats.tsv` and `detection_summary.tsv` byte-identical to before (M5 uses its own RNG streams) |
| 9 | an unknown flag reaching the writer | `ValueError` (closed vocabulary) |

## Edge cases

- A replicon with no eligible gene contributes NA autocorrelation rows only (M5-05) — the table is never
  missing a replicon.
- `lag` is an integer column with NA (Frictionless `integer` with `missingValues`).
- Very small floats are written `%.6g` (e.g. `1e-07`), as elsewhere.

## Acceptance criteria

- [ ] `tests/unit/io/test_context_tsv.py`: examples 2, 3, 9 and round trips (dtype-aware) of both tables.
- [ ] `tests/unit/context/test_vocab.py`: `CONTEXT_FLAGS` ⊇ §9.4 Region-and-context flags; schemas list the
  same set.
- [ ] `tests/integration/test_context_stats_m5.py`: examples 1, 4–8 on the fixture builder; rule checks of
  the schema section on the produced tables.
- [ ] R4 check: every autocorrelation row has a `q_value`, NA with `not_estimable`, or `descriptive`;
  every island row has `p_value` (FWER).

## Questions for the lead

1. New flags `not_estimable` (a statistic or p does not exist for the data: no pairs, no variance,
   degenerate null, too few TUs) and `few_values` (< 10 genes with values for Moran). **Proposed default:**
   add both to the schema vocabulary ("extend only through the schema", §9.4), as M4 did for
   `few_shifts`.
2. Who writes the M5 columns of `genes.tsv`? `context_stats` runs after `gene_table`, and §6.2 gives the
   final `genes.tsv` to `finalize`. **Proposed default:** `context_stats` writes `gene_context.parquet`;
   the genes.tsv writer (M4-15's, wherever M4-18/M4-20 call it, later `finalize`) left-joins it. If the
   M4 code writes genes.tsv before `context_stats`, move that call after it in this PR.
3. Schema file names: `<table>.schema.json` (existing files) vs `<table>.json` (M4-15/M4-18 specs).
   **Proposed default:** follow whatever M4-15/M4-18 merged; else the existing `.schema.json`.

## Lead decisions (2026-10-20, after the statistician's attack check)

36. Q1: add `not_estimable`, `few_values` and `ci_unestimated` to the flag vocabulary.
37. Q2, Q3: accept the proposed defaults.
