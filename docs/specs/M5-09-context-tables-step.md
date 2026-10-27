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
  `ci_few_units`; **not** `uncalibrated`, which ADR-0021 Amendment 1 drops) — "extend only through the schema"; §9.5 Frictionless resources; §6.5 coordinates and
  sort; §6.6 / R7 determinism; R4 (every claim carries FWER/q or `unestimated`).

## Public interface

- `src/gippa/io/context_tsv.py`: `write_autocorrelation(path, df)`, `write_islands(path, df)` and readers
  (inverse; islands 0-based → 1-based `start + 1`, `end` unchanged, `end > length` kept for origin
  crossers; lists `;`-joined; through M1-04 `write_tsv`, `%.6g`).
- `src/gippa/schemas/tables/autocorrelation.schema.json`, `islands.schema.json` (same layout as
  `detection_summary.schema.json`: name, type, unit in the name, description, enums, `coordinates`
  `1-based-closed` for islands, none for autocorrelation).
- `src/gippa/context/vocab.py`: `CONTEXT_FLAGS` = §9.4 Region-and-context flags + M4 additions
  (`few_shifts`, `split_contigs`, `hl_variance_nonpositive`, `no_reference_replicon`, `ci_few_units`) +
  **M5 additions `not_estimable`, `few_values`** (Q1); both schemas constrain `flags` items to it.
- `src/gippa/steps/context_stats.py` (M4-18's step) gains inputs `tus` (C5, required) and `features`
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

## Lead decisions after the test-author's report (2026-10-21)
1. **`gippa.context.diagnostics` is accepted as invented.** ADR-0021 Amendment 1 requires the per-replicon
   ℓ̂ / ℓ' / Î₁ to be reported, but `IslandScan` dies with the step, `detection_summary.tsv` must stay
   byte-identical (example 8), `autocorrelation.tsv`'s `statistic` enum is closed, and the §9.4 run report
   does not exist yet (`src/gippa/report/` holds only `__init__.py`). So the step persists
   `work/context/island_diagnostics.parquet` with `DIAGNOSTICS_COLUMNS` and `BLOCK_RULES` as proposed.
   The semantics are the contract and survive any later renaming: ℓ̂ and ℓ' are **NA, never 0**; a fixed
   block gives `block_rule="fixed"`, ℓ̂ NA, ℓ' = the configured value, Î₁ finite; not scanned gives NA/NaN
   with `few_values`; `render_l_hat` never emits the character `0`; `calibration_applies("fixed")` is
   False. When the §9.4 run report is built it **consumes this parquet** rather than recomputing.
2. **`uncalibrated` stays out of `CONTEXT_FLAGS`** as the tests assert, with the legacy
   `regions_stats.schema.json` pattern left alone. (The spec's line 19 correction *is* on main as of
   1189648 — re-read it.)
3. **`tables/genes.tsv` is written by `context_stats`**, not by `run_detectability`. It is the last M4/M5
   step and the only one holding both the gene table and the three context columns, so a single write
   with the left join already applied beats writing it early and patching it. **Known limitation:** a run
   where `context_stats` is skipped (no C7/C8) produces no `genes.tsv`. That is recorded, not hidden, and
   is resolved when §6.2's `finalize` step arrives and takes ownership of the file. M4-20 moves its write
   accordingly.
4. **Example 6's literal numbers are not reproducible on a simulated fixture**; testing the relation (every
   island's genes carry that island's id, replichore ∈ {right, left}, position ∈ [0, 1], the frame equals
   the producers' join) is the right substitute. The exact geometry is M5-08 example 1's job.
5. **The island fixture runs outside the ADR-0021 calibrated regime on purpose** — abundance-driven
   detection, a 0.30 detrend span and `context.islands.block_len = 1` are all needed to make an island
   reachable at this size (decision 50: censoring caps the direction). Acceptable for plumbing, it doubles
   as the fixed-block diagnostics case, and it is stated in the module docstring. Do **not** grow the
   genome to reach the calibrated regime here; that belongs to M5-07.
6. **`stat` is the statistic's *name* in every test table, islands included** (2026-10-21). §9.2 defines
   the test columns with `stat` as an enum of names (`enum{o_minus_e, hl_shift}` for regions), and
   `regions_stats` and `autocorrelation` both follow it. M5-06 shipped `stat` as a float equal to
   `scan_s`, which both breaks the convention and duplicates a column the same row already carries. So
   **islands' `stat` becomes the string `"scan_max"`** and `scan_s` keeps the value. The M5-09
   test-author's fixtures already say `"scan_max"` — they were right and the merged producer is the slip.
   This changes merged code (`src/gippa/context/islands.py`) and one merged oracle
   (tests/unit/context/test_islands.py:134); cheap now, expensive once goldens, the M9 report and any
   generic table reader depend on three tables meaning three different things by `stat`.
7. **The `flag_set` / `NA` contradiction is a test bug**: `test_every_row_keeps_its_own_flags` reads raw
   cells where `flag_set` returns `{"NA"}`, while `test_empty_flags_are_na` requires the cell to be
   exactly `NA`. Drop `"NA"` from `flag_set`; the integration path through `read_tsv` is unaffected.
8. The two changes outside M5-09's nominal surface are **accepted**: `oriter`'s `lag` as nullable `Int64`
   (without it the concat of M5-05's int64 and M5-08's float64 rows writes `1.0`), and
   `registry.optional_outputs`. The latter has no unit test — that is filed and must not be forgotten.
