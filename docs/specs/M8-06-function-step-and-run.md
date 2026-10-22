# M8-06: The `function` step, its three tables, and wiring them into `giap run`

Milestone M8 · track T (tables and integration) · **pass 1, the Module 7 integration milestone** ·
depends on: M8-03 (merged), M8-05 + M8-05b (merged), M8-04a (merged), **M8-04b (PR #179, in
review)** · estimate ≈ 420 lines in `src/` and `workflow/`, of which ≈ 130 are schema/registry JSON
(+ ≈ 520 test lines, + nf-test) · **split into M8-06a and M8-06b** (see "Split") · not statistical

**Revision 2 (2026-10-22).** Revision 1 was written before M8-03, M8-04 and M8-05 were built.
This revision matches the code those issues actually expose: `giap.function.collections` and
`giap.io.function_resources` (main), `giap.function.gsea` (main), and `giap.function.ora` /
`giap.function.summary` (branch `origin/m8/ora`, PR #179). Anything marked **[#179]** depends on
an interface that review could still change.

## The gap this closes

Without this issue, Module 7 repeats the project's recurring failure: tables that exist only in
Python. This issue:

- registers the step;
- writes the three tables;
- persists M4-13's TU folds, because ADR-0035's recalibration needs them;
- stages the five setup-manifest resources under both executors;
- adds the `FUNCTION` subworkflow and calls it from `main.nf`, switched by `params.modules`;
- proves it with a real Slurm run before any of M8's golden or calibration work starts.

## Spec sections implemented

- §6.2 step row "`function` | py | gene table, xrefs, GO/KEGG resources → `category_summary`,
  `enrichment_*` | function". `xrefs` is read through the gene table's function columns: M8-03
  item 1, and revision 1 Q2, accepted under index lead decision 1.
- §6.3 "Optional modules. `mobile`, `bgc` and `function` are switched with `params.modules`.
  Downstream steps accept a missing optional layer and record `layer_missing` in the manifest."
- §6.4: both executors give byte-identical tables.
- §6.6: each random draw has its own stream (index decision 26).
- §7.7.1: the resources; "the manifest notes `bakta_db = light`".
- §7.7.5: the three outputs.
- §9.1:
  - `tables/category_summary.tsv`, `enrichment_detection.tsv` and `enrichment_gsea.tsv`;
  - "Optional modules that are switched off write no files; the manifest lists them under
    `modules_off`".
- §9.2 conventions, §9.3 the three column dictionaries as amended (below), §9.4 run flags
  (`layer_missing:<layer>`, `bakta_db_light`), §9.5 Frictionless resources.
- §7.9.1: the setup items `go_basic`, `goslim_generic`, `kegg_ko_pathway`, `kegg_pathway_names`
  and `cog2020_fun`. Their paths come from `src/giap/data/db_sources.yaml`, and `giap setup --go
  --kegg --cog` installs them. All five are `ok` in FarmShare's `local/db/db_manifest.json`
  (retrieved 2026-10-15).

## What the merged and under-review code exposes (the inputs to this issue)

| Function | Module | Signature (as on main / #179) | Run flags it returns |
|---|---|---|---|
| resources | `giap.io.function_resources` | `read_resources(paths: Mapping[str, Path \| None]) -> Resources` (keys = the five item names; None or absent means not installed) | none |
| collections | `giap.function.collections` | `build_collections(genes, eligible: B8, res) -> (list[GeneSets], list[str])`, in `COLLECTIONS` order (cog, kegg, go_bp, go_mf, go_cc, go_slim) | `layer_missing:cog`, `layer_missing:kegg`, `layer_missing:go`, in that order; `bakta_db_light` when any `annotation_source == "bakta_light"` |
| eligibility | `giap.genetable.eligibility` | `eligible(genes) -> B8` (non-partial, finite `p_detect`) | none |
| summary **[#179]** | `giap.function.summary` | `category_summary(genes, colls, B: int, rng) -> DataFrame` | none (row flag `ci_few_units`) |
| ORA **[#179]** | `giap.function.ora` | `detection_enrichment(genes, colls, size: tuple[int, int], alpha, *, folds: Mapping[str, int] \| None) -> DataFrame` | `df.attrs["run_flags"]`: `tu_recal_skipped`, `p_detect_clipped`, `null_degenerate` |
| GSEA | `giap.function.gsea` | `gsea_table(genes, colls, size, cfg: GseaConfig, alpha, rng, replicons=None) -> DataFrame` | `df.attrs["run_flags"]`: `layer_missing:abundance` or `too_few_ranked_tus` |

Consequences for the step:

- **No `tus` input.** ORA, GSEA and the summary read `tu_id` and `tu_position` from the gene
  table (C5 is already merged into it upstream). Revision 1's required `tus` key is dropped.
- **`replicons` (C1) is required.** GSEA's `genomic_cluster` flag wraps circular replicons only
  (`topology == "circular"`). Without C1, every replicon would be treated as linear.
- **ORA needs M4-13's TU folds [#179, ADR-0035].** `tu_universe(genes, folds)` cross-fits the TU
  recalibration on `DetectabilityResult.fold_of_tu`. With `folds=None` or `{}` it falls back to
  `q_raw` and adds `tu_recal_skipped`. Today the folds are not written anywhere: the
  `detectability` step keeps only `p_detect`. This issue adds the output (item 1).
- **`eligible_genes` fails closed.** An eligible gene with a NaN `p_detect` or no `tu_id` raises
  `ValueError`. The step lets that propagate (exit 4), the same as `detectability`, which already
  raises for a null `tu_id`.
- **No `annotation_meta` input.** `build_collections` detects `bakta_db_light` from the gene
  table's `annotation_source`.

## Public interface

`src/giap/steps/function.py`
```python
@step("function", module="function",
      inputs=["gene_table", "replicons", "tu_folds"],
      optional_inputs=["go_basic", "goslim_generic", "kegg_ko_pathway", "kegg_pathway_names",
                       "cog2020_fun"],
      outputs=["category_summary", "enrichment_detection", "enrichment_gsea"])
def function_step(ctx: StepContext) -> None: ...

def run_function(genes: pd.DataFrame, replicons: pd.DataFrame, folds: Mapping[str, int],
                 res: Resources, cfg: GiapConfig
                 ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    """(category_summary, enrichment_detection, enrichment_gsea, run_flags). Pure."""
```
`src/giap/io/function_tables.py`
```python
def write_category_summary_tsv(path: Path, df: pd.DataFrame) -> None: ...
def write_enrichment_detection_tsv(path: Path, df: pd.DataFrame) -> None: ...
def write_enrichment_gsea_tsv(path: Path, df: pd.DataFrame) -> None: ...
```
`src/giap/setupdb/installed.py`
```python
FUNCTION_RESOURCES: tuple[str, ...] = ("go_basic", "goslim_generic", "kegg_ko_pathway",
                                       "kegg_pathway_names", "cog2020_fun")
def resource(db_dir: Path | None, name: str) -> Path | None: ...   # cog_table generalised
def cog_table(db_dir: Path | None) -> Path | None: ...            # = resource(db_dir, "cog2020_def")
```
Also `src/giap/schemas/tables/{category_summary,enrichment_detection,enrichment_gsea}.schema.json`.

### 1. TU folds: a new `detectability` output

- `detectability` gains the output key `tu_folds`, listed in `optional_outputs`. While a caller
  leaves it unwired (M8-06a, Nextflow), it is not written. The file is
  `genetable/tu_folds.parquet`, with columns `tu_id` (string) and `fold` (int64), sorted by
  natural TU ID and written with `write_parquet_deterministic`.
- Its rows are exactly `DetectabilityResult.fold_of_tu`. When the model is constant
  (`detectability_constant`: fewer than 50 detected or 50 undetected eligible genes, or too few
  TUs), the file is header-only.
- The function step reads it into `dict[str, int]` and passes it as `folds=`. A header-only file
  gives `{}`, so `tu_recal_skipped` is the documented outcome on a constant model. This is what
  sim100k gives (Q1).

### 2. Order inside the step (`run_function`)

1. `el = giap.genetable.eligibility.eligible(genes)`. `res = read_resources({k: ctx.inputs.get(k)
   for k in FUNCTION_RESOURCES})`. A resource file that does not parse (for example a bad OBO)
   raises, and the step exits 4. There is no silent fallback.
2. `colls, flags = build_collections(genes, el, res)`.
3. `summary = category_summary(genes, colls, cfg.stats.bootstrap, for_step("function.bootstrap",
   cfg.seed))` **[#179]**.
4. `detection = detection_enrichment(genes, colls, tuple(cfg.function.ora_size),
   cfg.stats.alpha.enrichment, folds=folds)` **[#179]**. ORA and elim are deterministic and draw
   no randomness.
5. `gsea = gsea_table(genes, colls, tuple(cfg.function.gsea_size), cfg.function.gsea,
   cfg.stats.alpha.enrichment, for_step("function.gsea", cfg.seed), replicons)`.
6. `run_flags` is the union of `flags`, `detection.attrs["run_flags"]` and
   `gsea.attrs["run_flags"]`, without duplicates, in first-seen order.

The streams follow the `context_stats.<part>` precedent: sub-streams from
`giap.rng.for_step(name, cfg.seed)`, not the step's own `ctx.rng`. Adding a draw to one part
therefore never moves the other.

**No eligible gene** (an empty gene table, or every gene partial). The step does not call steps 3
to 5. It writes three header-only files with run flags `flags + ["no_eligible_genes"]`. This guard
exists because the [#179] TU universe would be empty, and `tilt_pair` on an empty array is not
pinned anywhere.

### 3. Run flags and events (§6.7)

- Each run flag is logged once as a WARNING line `function: run flag <flag>`, as `codon_usage` and
  `operons` do.
- The step then emits one event `function.run_flags` with `{"flags": [...]}`.
- It also emits one event `function.collections` with `{"<collection>": {"n_terms", "coverage",
  "n_unknown_ids"}}` for each built collection, for the report's coverage line (§7.7.1).
- `bakta_db_light` is the §9.4 run flag that carries §7.7.1's "the manifest notes `bakta_db =
  light`". The manifest step (§6.7) reads it from the event. No separate event is added.

### 4. Writers (§9.2, ADR-0006)

- Cells use `giap.io.tsv.write_tsv`: `NA`, `%.6g`, `true`/`false`, `;`-joined `leading_edge` and
  `flags`. A NA `redundant` (boolean dtype) or `parent_group` (None) is written as `NA`.
- Column order is each module's `COLUMNS`, verbatim.
- Rows are sorted stably by (collection in `COLLECTIONS` order, term, and direction for ORA).
  ORA and the summary are already in that order. GSEA is sorted by the writer. A header-only
  frame writes the header row only.
- No coordinates are written.

### 5. Columns and schemas (§9.3 as amended by ADR-A / index decision 12, ADR-0029, M8-04 decision 3)

| Table | Columns (in order) | §9.3 departures |
|---|---|---|
| `category_summary.tsv` | `collection` enum{cog, kegg, go_slim}, `category`, `name`, `n_genes`, `n_detected`, `expected`, `o_minus_e`, `ci_low`, `ci_high`, `median_abundance_log2`, `coverage`, `flags` | `flags` added (last); §9.3 has none (Q4) |
| `enrichment_detection.tsv` **[#179]** | `collection` (6), `term`, `name`, `direction` enum{over, under}, `n_universe`, `n_in_set`, `n_detected`, `k_detected_in_set`, `n_in_set_genes`, `k_detected_in_set_genes`, `odds_weight`, `p_classic`, `p_noelim`, `parent_group`, `coverage`, `stat`, `effect`, `ci_low`, `ci_high`, `p_value`, `q_value`, `method`, `null_model`, `n_perm`, `flags` | `n_in_set_genes`, `k_detected_in_set_genes` and `p_noelim` added (M8-04 decision 3); counts in TUs; `flags` last |
| `enrichment_gsea.tsv` | `collection` (6), `term`, `name`, `size_tus`, `es`, `nes`, `log2err`, `leading_edge`, `redundant`, `parent_group`, `coverage`, `stat`, `effect`, `ci_low`, `ci_high`, `p_value`, `q_value`, `method`, `null_model`, `n_perm`, `flags` | `flags` last |

Schema enums and patterns:

- **`collection`:**
  - the enrichment tables take {cog, kegg, go_bp, go_mf, go_cc, go_slim} (index decision 7);
  - the summary takes {cog, kegg, go_slim}.
- **`method`:**
  - detection: {`cpb_tu_ora`};
  - GSEA: {`fgsea_simple`, `fgsea_multilevel`}. `fgsea_multilevel` is reserved for M8-11, so its
    switch needs no schema change.
- **`null_model`:** detection {`conditional_poisson_binomial_tu`}; GSEA {`tu_set_permutation`}.
- **`n_perm`:** `0` in detection.
- **`flags` patterns** (`;`-joined tokens, the pattern style of `pairs.schema.json`):
  - summary: `ci_few_units`;
  - detection: `ci_unestimated`, `unestimated`, `odds_out_of_range`, `elim_empty`,
    `null_degenerate`, `p_detect_clipped`, `ora_uncalibrated`. The last is ADR-0035's switch, in
    `ora.CALIBRATION_FLAGS`, and is empty today;
  - GSEA: `ci_unestimated`, `unestimated`, `size_filtered`, `genomic_cluster`,
    `detection_confounded`.
- **Integer and number constraints:**
  - counts are integers ≥ 0;
  - `p_*` and `q_value` are numbers in [0, 1];
  - `coverage` is in [0, 1].
- **Metadata:**
  - descriptions state the units ("in TUs" for `n_universe`, `n_in_set`, `n_detected`,
    `k_detected_in_set` and `size_tus`);
  - `missingValues ["NA", ""]`;
  - the `coordinates` key follows `pairs.schema.json`, a table without coordinate fields.
- **Not in any table:** the run-only tokens (`layer_missing:*`, `bakta_db_light`,
  `tu_recal_skipped`, `too_few_ranked_tus`, `no_eligible_genes`) are not row flags. They go to
  the event (3).

Revision 1's `schemas/datapackage.json` does not exist on main. M6-08 did not create one either,
so this criterion is dropped (Q7).

### 6. Module switch (§6.3)

- **Local.** `giap.runner.local._enabled` already plans a `module="function"` step only when
  `"function" in cfg.modules`. The default config includes it.
- **Nextflow.** `main.nf` calls `FUNCTION` only when `'function' in (params.modules ?: [])`. This
  is the first module switch in `main.nf`. `mobile` and `bgc` are not switched there yet (Q5).
- **Off.** No file is written, and there is no new event: `modules_off` is derived from the
  resolved `modules` by the §6.7 manifest step, which does not exist yet. Revision 1's
  `function.skipped` event is dropped (Q3).

### 7. Resource staging (both executors)

- `installed.resource(db_dir, name)` returns `db_dir / SOURCES[name].path` when the setup manifest
  lists `name` and the file exists, otherwise None. This is the `cog_table` rule (M4-15
  decision 5). `cog_table` becomes a one-line wrapper.
- **Local.** `cli._local_roots` stages each of `FUNCTION_RESOURCES` that `resource()` resolves,
  under its item name, and only when `"function" in cfg.modules`, so `--modules core` copies no
  32 MB OBO.
- **Nextflow.** `main.nf` adds `functionResources()`. It returns a `[key, path]` channel of the
  installed items, built from `installed(name, relative)` with these relative paths:
  - `go_basic`: `go/go-basic.obo`;
  - `goslim_generic`: `go/goslim_generic.obo`;
  - `kegg_ko_pathway`: `kegg/ko_pathway.tsv`;
  - `kegg_pathway_names`: `kegg/pathway_names.tsv`;
  - `cog2020_fun`: `cog/fun-20.tab`.

  A unit test pins these strings against `db_sources.yaml` (the five paths and `cogTable()`'s).

### 8. Nextflow wiring

`workflow/subworkflows/local/function.nf` follows the OPERONS/DETERMINANTS style:

- the alias is `GIAP_STEP as FUNCTION_STEP`, because an alias cannot share the subworkflow's
  name `FUNCTION`;
- it uses `stepIn`, `kv` and `pick` from `./annotate`;
- it runs at `process_single`, because the step needs < 2 GB and seconds to minutes (index
  feasibility).

```groovy
workflow FUNCTION {
    take:
    gene_table   // genetable/genes.parquet
    annotation   // [genes.parquet, features.parquet, replicons.parquet]
    tu_folds     // genetable/tu_folds.parquet
    resources    // [key, path] for each installed function resource, or channel.empty()
    main:
    def chans = [kv('gene_table', gene_table), kv('replicons', annotation.map { a -> a[2] }),
                 kv('tu_folds', tu_folds), resources]
    def v = stepIn('function', chans, ['gene_table', 'replicons', 'tu_folds'])
    FUNCTION_STEP(v.files, v.keys, ['category_summary=tables/category_summary.tsv',
        'enrichment_detection=tables/enrichment_detection.tsv',
        'enrichment_gsea=tables/enrichment_gsea.tsv'])
    emit: category_summary, enrichment_detection, enrichment_gsea (pick), versions
}
```
Other wiring changes:

- **`genetable.nf`.** DETECTABILITY's out spec gains `'tu_folds=genetable/tu_folds.parquet'`, and
  the subworkflow emits `tu_folds`.
- **`main.nf`.** After `GENETABLE`, beside `CONTEXT` and `OPERONS`, and inside the
  `inputs.peptides` branch:
  `if ('function' in (params.modules ?: [])) FUNCTION(GENETABLE.out.gene_table, gene_model,
  GENETABLE.out.tu_folds, functionResources())`.
- **Publishing (M6-08 decision 7).** Only the three `tables/*.tsv` files and
  `genetable/tu_folds.parquet` are named. `giap.runner.local.PUBLISHED_FILES` gains the same four
  keys, which `tests/unit/test_local_publish_names.py` pins against the Nextflow names.
- **Registry pin.** `workflow/tests/data/annotate/step_registry.json` gains the `function` entry,
  and `detectability`'s outputs gain `tu_folds`.
- **Process-list pin.** `workflow/tests/main_genetable_context.nf.test`'s "every module" test runs
  with `modules` including `function` and `db_dir = data/genetable/db` (COG def only), so
  `FUNCTION_STEP` joins its `want` list. Under `modules = ['core']` (main_operons, the "no
  peptides" case) it must not appear. `genetable.nf.test`'s published-file list gains
  `genetable/tu_folds.parquet`.
- **Conflict note.** #178 (M8-10) edits the same `main.nf` lines, `PUBLISHED_FILES`,
  `step_registry.json` and `want` list (it adds `DETERMINANTS_STEP`). Whichever PR merges second
  merges `origin/main` and keeps both entries.

## Split (PRs ≤ ~400 lines)

| PR | Content | Code | Tests (test-author) |
|---|---|---|---|
| **M8-06a**: step and tables (Python only) | items 1 to 5; `installed.resource`; `step_registry.json`. `tu_folds` stays optional-output and unwired in Nextflow; nothing new is staged by `cli` or published. Under `--executor local` the step is planned (default modules) and writes to `work/` only, with no resources staged | ≈ 210 Python + ≈ 125 JSON | `tests/unit/steps/test_function_step.py` (examples 4 to 9, via `run_step`); writer and schema tests; `tests/unit/steps/test_detectability_folds.py` (example 7); a synthetic function DB helper (Q6) |
| **M8-06b**: `giap run` under both executors | items 6 to 8: `cli` staging, `PUBLISHED_FILES`, `function.nf`, `main.nf`, `genetable.nf`, the pins; the real Slurm run | ≈ 85 Groovy/Python | `tests/integration/test_m8_06_function_run.py` (examples 1 to 3); `workflow/tests/function.nf.test`; `workflow/tests/main_function.nf.test`; pin updates |

M8-06a imports `giap.function.ora` and `giap.function.summary`, so until #179 merges,
`m8/function-step` merges `origin/m8/ora` (stacked, no rebase). M8-06a merges after #179.
M8-06b stacks on M8-06a. One PR at ≈ 420 lines (schema JSON included) is the alternative
(Q8).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| `--modules core` (function off) | step not planned (local) / `FUNCTION` not called (Nextflow); no file; no error; no event |
| no GO/KEGG/COG resource installed (GitHub CI, a fresh machine, `data/genetable/db`) | `colls` empty; three header-only files; run flags `layer_missing:cog`, `layer_missing:kegg`, `layer_missing:go` (+ `tu_recal_skipped` if the model is constant) |
| only some resources (e.g. COG only) | COG rows only; `layer_missing:kegg`, `layer_missing:go` |
| `goslim_generic` missing, `go_basic` present | `go_slim` collection absent (M8-03); `category_summary.tsv` has cog/kegg rows only |
| a resource listed in the manifest but its file deleted | not installed (`resource()` rule) |
| a resource file that does not parse | step fails, exit 4 (no fallback) |
| no abundance | `category_summary.tsv` rows with `median_abundance_log2 = NA`; detection rows as usual; `enrichment_gsea.tsv` header-only; `layer_missing:abundance` |
| < 2 ranked TUs | `enrichment_gsea.tsv` header-only; `too_few_ranked_tus` |
| constant detectability model (sim100k) | `tu_folds.parquet` header-only → `folds = {}` → `tu_recal_skipped`; ORA on `q_raw` |
| no eligible gene / empty gene table | three header-only files; `no_eligible_genes` (step 2 guard) |
| an eligible gene with NaN `p_detect` or no `tu_id` | `ValueError` from `eligible_genes` → exit 4 naming the genes (fail closed) |
| circular vs linear replicon | only `genomic_cluster` uses topology: a set whose TUs sit at both ends of a **circular** replicon is adjacent across the origin, on a **linear** one it is not (M8-05b); every other column is topology-free |
| minus-strand genes | no effect: the tables have no coordinates; TU membership and `tu_position` (transcription order) come from C5 |
| ties in TU scores (e.g. imputed abundances) | ordered by the seeded `function.gsea` permutation, never by TU ID (M8-02 decision 7): reproducible for a seed, may differ across seeds |
| ties for `parent_group` | most general ancestor by depth, then by ID (M8-04/M8-05) |
| `bakta_light` annotation | run flag `bakta_db_light`; tables as usual (low coverage) |
| an upstream channel empty in Nextflow (no peptides) | `stepIn`'s `need` holds the task; `FUNCTION_STEP` does not run |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `giap run --executor local --modules core,function` on sim100k (M3-19 proteomics, precomputed layers, seed 7) with `$GIAP_DB_DIR` = the synthetic function DB (Q6) | `OUT/tables/{category_summary,enrichment_detection,enrichment_gsea}.tsv` exist and validate; `category_summary.tsv` and `enrichment_detection.tsv` have ≥ 1 row; every detection and GSEA row's `flags` contains `ci_unestimated`; the `function.run_flags` event contains `tu_recal_skipped` (constant model on ~100 genes); exit 0 |
| 2 | example 1 through `--executor nextflow --profile farmshare_local` | the three tables and `genetable/tu_folds.parquet` byte-identical to example 1 (§6.4) |
| 3 | example 1 with `--modules core` | no `tables/category_summary.tsv`, `enrichment_detection.tsv` or `enrichment_gsea.tsv`; `tables/operons.tsv` still written; exit 0 |
| 4 | `run_step("function", …)` on a 4-gene toy gene table (g1–g4 in TU_1…TU_4, non-partial, `cog_category` J, J, J, NA; `detected` T, T, F, F; `p_detect` 0.9, 0.5, 0.4, 0.3; no abundance) with only `cog2020_fun` (J = "Translation, ribosomal structure and biogenesis") | `category_summary.tsv` = header + one row `cog  J  Translation, ribosomal structure and biogenesis  3  2  1.8  0.2  NA  NA  NA  0.75  ci_few_units`; `enrichment_detection.tsv` header-only (3 TUs < `ora_size` 5); `enrichment_gsea.tsv` header-only; run flags `layer_missing:kegg`, `layer_missing:go`, `tu_recal_skipped`, `layer_missing:abundance` |
| 5 | example 4 with no resource at all | three header-only files whose header rows equal `summary.COLUMNS`, `ora.COLUMNS`, `gsea.COLUMNS` tab-joined; run flags start `layer_missing:cog`, `layer_missing:kegg`, `layer_missing:go` |
| 6 | example 4 with every gene `partial = true` | three header-only files; `no_eligible_genes`; exit 0 |
| 7 | `detectability` on a 600-gene non-constant fixture (the M4-13 test frame) | `tu_folds.parquet` rows == `fit_detectability(...).fold_of_tu`, sorted by natural TU ID; feeding it to `function` gives no `tu_recal_skipped`; replacing it with a header-only file gives `tu_recal_skipped` and different `p_value`s |
| 8 | the same inputs and seed twice; then seed 7 vs 8 | byte-identical tables; with seed 8, `enrichment_detection.tsv` is unchanged (deterministic), `category_summary.tsv` CIs and GSEA p-values may change |
| 9 | `frictionless validate` on the three outputs of examples 1, 4, 5 | valid; a row with `method = wallenius_ora` or flag `foo` is rejected |
| 10 | `sbatch scripts/slurm/giap_run.sbatch` on sim100k with FarmShare's `local/db` (real GO/KEGG/COG) | completes; the three tables exist and validate (they may be sparse: random simulator IDs); job ID in PROGRESS.md |

Example 4's numbers: n_genes = 3, n_detected = 2, expected = 0.9 + 0.5 + 0.4 = 1.8, O − E = 0.2,
and 3 TUs < `CI_MIN_UNITS` 10, so the CI is NA with `ci_few_units`. Coverage is 3/4 = 0.75,
because g4 is eligible but has no COG.

## Acceptance criteria

M8-06a:
- [ ] `tests/unit/steps/test_function_step.py` covers examples 4, 5, 6, 8 and 9. The example 4 row
  is exact. The order of the run-flag union is pinned.
- [ ] `tests/unit/steps/test_detectability_folds.py` covers example 7. The existing detectability
  tests pass unchanged, and `tu_folds` unwired means the file is not written.
- [ ] Writer tests:
  - list and boolean cells;
  - `NA` for a None `parent_group` and a NA `redundant`;
  - the stable sort of an unsorted GSEA frame;
  - header-only output.
- [ ] Schemas:
  - the enums and flag patterns exactly as item 5;
  - the column order equals each module's `COLUMNS`. A unit test asserts this, so a #179 column
    change fails loudly.
- [ ] `installed.resource`:
  - listed and present → path;
  - listed and absent → None;
  - not listed → None;
  - `cog_table` unchanged in behaviour.
- [ ] `step_registry.json` matches the registry (`test_nf_annotate_registry.py`).
- [ ] Benchmark `function_step_4500_genes` (the 4,500-gene simulated layout with synthetic
  resources, `stats.bootstrap` = 1000) is recorded in `benchmarks/baseline.json`.

M8-06b:
- [ ] `tests/integration/test_m8_06_function_run.py` covers examples 1, 3 and 9. Example 2 is
  marked `needs_nextflow`, as in `test_m6_08_operons_run.py`.
- [ ] `workflow/tests/function.nf.test` (stub):
  - the three keys are emitted with their §9.1 names;
  - the resources are staged under their item names only when given;
  - no `tu_folds` means no run.
- [ ] `workflow/tests/main_function.nf.test`:
  - `modules` with `function` → `FUNCTION_STEP` runs once;
  - `['core']` → it does not run;
  - no peptides → it does not run.
- [ ] Pins:
  - `main_genetable_context.nf.test`'s `want` gains `FUNCTION_STEP`;
  - `genetable.nf.test` gains `tu_folds.parquet`;
  - `test_local_publish_names.py` passes;
  - the `functionResources()` paths are pinned against `db_sources.yaml`.
- [ ] **Evidence that counts (the MVP rule):** example 10, one real `sbatch
  scripts/slurm/giap_run.sbatch` that reaches the end. Its job ID goes in PROGRESS.md.

## Hardening owed (pass 2)

- A meaningful real-data run: the K-12 full-DB Bakta recording
  (`tests/data/recorded/bakta/1.12.2-full/k12`) with M3-22's in-silico peptides and a simulated
  abundance column. Detectability is then non-constant, so the recalibration runs, and ORA and
  GSEA have real GO/KEGG/COG structure (M8-15's real-data sanity check).
- Contract assertions on the three tables inside the pipeline (the MVP-01 pass-2 pattern).
  `-resume` across `FUNCTION`, and the real-tool arm under `-profile farmshare,apptainer`.
- Golden tables for the 50-kb fixture, once M8-14 gives the simulator real-shaped annotations.
- Manifest (§6.7 step):
  - the GO release and the KEGG retrieval date in `run_manifest.json` (ADR-0017);
  - `modules_off`;
  - the run flags read from the `function.run_flags` event.
- `mobile` / `bgc` module switches in `main.nf` (Q5).
- The run-scope flag tokens of item 5 join §9.4's vocabulary when the manifest schema gets a
  run-flag enum.

## Questions for the lead

1. **TU folds source (ADR-0035).** Options:
   - (a) persist `fold_of_tu` as `genetable/tu_folds.parquet`, an optional output of
     `detectability` (proposed);
   - (b) recompute the folds in `function` by replaying `tu_folds(…, for_step("detectability",
     seed))`. This depends on `tu_folds` staying the first draw of that stream and on knowing
     whether the model was constant;
   - (c) pass `folds=None`. Every run would then be `tu_recal_skipped`, which contradicts the
     accepted ADR.

   **Proposed default:** (a). It is small, explicit and executor-neutral, and `tu_recal_skipped`
   on sim100k is expected and pinned.
2. **Stacking on #179.** M8-06a needs `ora.py` and `summary.py`. **Proposed default:** merge
   `origin/m8/ora` into `m8/function-step` now (plain merge, no rebase), and merge M8-06a only
   after #179. Any #179 review change to `COLUMNS`, the run flags or the `folds=` keyword is
   picked up by the schema-equals-`COLUMNS` test.
3. **`modules_off` / the skipped event.** Revision 1 promised a `function.skipped` event. **Proposed
   default:** drop it. The §6.7 manifest derives `modules_off` from the resolved `modules`, which
   is already in `params.json`, and no executor emits per-module events today.
4. **`category_summary.tsv` `flags`.** §9.3 lists no `flags` column, but #179's `COLUMNS` has one
   (`ci_few_units`). **Proposed default:** keep it as an ADR-A extension (non-breaking, last
   column). The lead records it, with ORA's three added columns (M8-04 decision 3), in an
   ADR-0029 / ADR-A amendment line. This spec-writer does not edit ADRs.
5. **First module switch in `main.nf`.** `params.modules` is not consulted anywhere in `workflow/`
   today. **Proposed default:** switch `function` only, here. A backlog line covers
   `mobile`/`bgc`, which need their own steps to exist first.
6. **Test DB for examples 1 and 2.** The M3-17 test DB has no function resources, and sim100k's
   GO, KEGG and COG IDs are random. **Proposed default:** a test helper
   (`tests/helpers/function_db.py`, test-author) writes a synthetic DB with these contents:
   - a toy OBO whose terms are the sim's own GO IDs under the three roots, plus a few parent terms
     sized ≥ 10 TUs;
   - a synthetic `ko_pathway.tsv` mapping the sim's KOs to made-up `path:map9xxxx` IDs;
   - `fun-20.tab` and the sim's `cog-20.def.tab`.

   It contains no KEGG content (index lead decision 4). The real FarmShare `local/db` is used only
   by example 10.
7. **`schemas/datapackage.json`.** It does not exist (M6-08 never created it). **Proposed
   default:** drop the criterion here, and add a backlog line for one package listing every
   table schema.
8. **Split.** **Proposed default:** M8-06a / M8-06b as above. The alternative is one PR at ≈ 420
   lines, if the lead counts the ≈ 125 lines of schema JSON outside the budget (the M6-08
   precedent shipped as one PR).
9. **Run-scope tokens.** `tu_recal_skipped`, `too_few_ranked_tus`, `no_eligible_genes`,
   `p_detect_clipped` and `null_degenerate` are run flags that §9.4's run scope does not list.
   **Proposed default:** emit them in the `function.run_flags` event now. They join the §9.4
   vocabulary through the manifest schema (hardening), as `operon_source_distance_model` did.
10. **sim100k is constant-model.** With about 100 genes, M4-13's `MIN_CLASS` of 50 makes
    `p_detect` constant, so every pass-1 integration run is `tu_recal_skipped`. The recalibration
    is exercised only by example 7 and the benchmark. **Proposed default:** accept this for pass
    1 (the M6-08 "sparse sim100k" precedent). The K-12 hardening run is the first real
    recalibrated run.

## Lead decisions (2026-10-22)
1. **Questions 1–10: the proposed defaults are accepted.**
   - TU folds are persisted as the new `detectability` output `genetable/tu_folds.parquet`.
   - The branch is stacked on #179 (origin/m8/ora merged in), and M8-06a merges after #179.
   - No `function.skipped` event.
   - `category_summary` keeps `flags` as its last column.
   - Only the `function` module switch goes into main.nf.
   - The integration examples use a synthetic test DB built from the simulator's IDs, with no KEGG content.
   - The datapackage criterion is dropped and goes to the backlog.
   - The work splits into 06a and 06b.
   - New run-scope flags are emitted in the event now.
   - sim100k always skipping the recalibration is accepted for pass 1.
2. **Question 4, the columns that depart from §9.3:** the `category_summary` `flags` column and ORA's `n_in_set_genes`, `k_detected_in_set_genes` and `p_noelim`. They go in an ADR-0029 amendment, which the lead drafts as *proposed* for Khoa. 06a goes ahead on the default. If Khoa rejects the amendment, only the schemas and `COLUMNS` change.
3. **Decision 13 of #179 applies:**
   - `folds` is a required keyword, so `function` passes `folds=None` explicitly when `tu_folds.parquet` has only a header.
   - Every ORA row carries `ora_uncalibrated` in pass 1.
   - A skipped recalibration adds `tu_recal_skipped` to the rows as well.
   - The schema-equals-`COLUMNS` test pins the final #179 columns.
4. **The M8-06b test contract is accepted (TA 9a931e5).** Q6's premise was wrong: sim100k writes `db_xrefs: []`, so the simulator has no GO, KEGG or COG IDs. The integration fixture therefore writes deterministic `db_xrefs` into the simulated `bakta_like.json` before either run, and builds the synthetic DB from those IDs. Giving giapsim real xrefs goes to the backlog (M8-14). The lead added `FUNCTION_STEP` to the MVP-02 process-list pin and `genetable/tu_folds.parquet` to the genetable published files.

