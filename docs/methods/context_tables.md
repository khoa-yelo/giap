# The context tables: `autocorrelation.tsv`, `islands.tsv` and the genes.tsv context group

Issue M5-09 (`docs/specs/M5-09-context-tables-step.md`, lead decisions 1–5). Spec sections §6.2
(the `context_stats` step), §9.2 (table conventions and test columns), §9.3 (the two column
dictionaries), §9.4 (the flag vocabulary), §9.5 (Frictionless resources), §6.5 (coordinates) and
§6.6 / R7 (determinism).

## Where each row comes from

| Table | Rows | Producer |
|---|---|---|
| `tables/autocorrelation.tsv` | Moran's I and join counts per replicon, lag and null | `giap.context.autocorrelation.autocorrelation_rows` (M5-05) |
| | the pooled / right / left ori–ter slopes per circular chromosome | `giap.context.oriter.oriter_rows` (M5-08) |
| `tables/islands.tsv` | one row per merged significant window | `giap.context.islands.islands_frame` (M5-06) |
| `work/context/gene_context.parquet` | `island_id`, `replichore`, `replichore_position` per gene | `gene_island_map` + `gene_replichores` |
| `work/context/island_diagnostics.parquet` | l-hat, l′ and the lag-1 Moran's I per replicon | `giap.context.diagnostics.island_diagnostics` |

`giap.context.tables` only concatenates and sorts: autocorrelation by contig (C1 order),
`statistic` and `null_model` in §9.3 enum order, `lag` (NA last) and `stat`; islands by contig,
`start`, `end` and the `ISL_<n>` number. Both sorts are stable, so rows that tie keep the order
their producer gave them, and no island is renumbered.

`giap.context.oriter` emits `lag` as the schema's nullable `Int64` rather than `float64`, so the
concatenation of the two row sources keeps integer lags (`1`, not `1.0`) in the written table.

## `stat` is a name, in all three tables

§9.2 defines the test column `stat` as an enum of **statistic names** (`enum{o_minus_e, hl_shift}`
for `regions_stats`), and `autocorrelation` follows it (`morans_i`, `join_count_bb`,
`oriter_slope_pooled`, …). M5-06 originally shipped islands' `stat` as a float equal to `scan_s`,
which both broke the convention and duplicated a column the same row already carries; M5-09 lead
decision 6 fixed it, so islands' `stat` is the string `"scan_max"` (`giap.context.islands.STAT`)
and `scan_s` keeps the value. A generic reader of the three test tables can therefore treat `stat`
as a string everywhere.

## Random streams (§6.6)

`context_stats.autocorrelation`, `context_stats.islands` and `context_stats.oriter`, each from
`giap.rng.for_step(..., cfg.seed)`. The island stream is shared by `scan_replicons` and then
`islands_frame`, in that order. The M4 parts keep their own streams, so changing an M5-only knob
(the detrend span, the block length, the presence of an oriC) leaves `regions_stats.tsv` and
`detection_summary.tsv` byte-identical.

## Writer conventions

`giap.io.context_tsv` hands cells to M1-04 `write_tsv` **by value**: `%.6g` for floats, integers
as integers, `NA` for missing (including a row with no flag), `;` for the `flags` and `genes`
lists. Nothing is coerced to its schema type on the way out, so a frame that violates the schema
is written as it stands and caught by `frictionless validate` rather than silently rewritten.
Only `start` is converted, to 1-based; `end` is unchanged, because a half-open 0-based end equals
a closed 1-based end, and an island that crosses the origin therefore keeps `end > length` (§6.5).
An empty frame writes the header line only.

The readers are the inverses and type each column from the shipped schema: `integer` → pandas
nullable `Int64`, `number` → `float64`, `string` → `object` of `str`, and the list cells → tuples
of `str`.

`flags` is checked against `giap.context.vocab.CONTEXT_FLAGS` (§9.4 is closed) and an unknown flag
raises `ValueError` naming it. Note that `pandas.Series.flags` shadows the `flags` *column*: a
writer reaching for `row.flags` on an `iterrows` row gets a `Flags` object, so every access is by
name or by position.

## The island diagnostics, and why 0 is not a value

ADR-0021 **Amendment 1** (decision 2a) and M5-06 decision 48 dropped the `uncalibrated` flag: no
per-replicon statistic certifies the validity domain, so a flag would be false reassurance.
Instead the run must report, per replicon, the block-length estimate l-hat, the bootstrap block
length l′ and the lag-1 Moran's I of the detrended TU scores. `IslandScan` carries them, but it
dies with the step, so `context_stats` persists them as
`work/context/island_diagnostics.parquet`; the §9.4 run report will consume that file rather than
recompute (lead decision 1).

`IslandScan.l_hat` is `0` in **two** different situations — the replicon was not scanned, and
`context.islands.block_len` fixed the block length so l-hat was never estimated — and
`IslandScan.block_len` is `0` when the replicon was not scanned. Printed as a number, l-hat = 0
reads "no autocorrelation", which is precisely the false reassurance Amendment 1 refuses. The
persisted table therefore carries **NA, never 0**, `render_l_hat` prints a sentence that never
contains the character `0`, and `block_rule` keeps the two reasons apart:

| `block_rule` | l-hat | l′ | lag-1 Moran's I | `flags` | calibration |
|---|---|---|---|---|---|
| `auto` | ≥ 1 | ≥ 1 | finite | none | applies |
| `fixed` | NA | the configured value | finite | none | does **not** apply |
| `not_scanned` | NA | NA | NaN | `few_values` | does **not** apply |

`calibration_applies` is True only for `auto`: Amendment 1's FWER evidence is conditional on
ADR-0021 decision 1's rule l′ = min(20, 2·l-hat + 1, T), so a run with an overridden `block_len`
is outside the regime the calibration covers — as the M5-09 island fixture deliberately is
(lead decision 5).

## genes.tsv

`context_stats` writes `tables/genes.tsv` with `gene_context` left-joined on `locus_tag` (lead
decision 3): it is the last M4/M5 step and the only one holding both the gene table and the three
context columns, so one write with the join already applied beats writing the file early and
patching it. A run that skips `context_stats` (no C7/C8) therefore produces no genes.tsv; that is
recorded in `docs/backlog.md` and resolved when §6.2's `finalize` step takes the file over.
