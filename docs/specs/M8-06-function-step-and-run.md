# M8-06: The `function` step, its three tables, and wiring them into `giap run`

Milestone M8 · track T (tables and integration) · **pass 1 — the Module 7 integration milestone** ·
depends on: M8-03, M8-04, M8-05 (and merged MVP-01 wiring) · estimate ≈ 320 code lines
(+ ≈ 350 test lines, + nf-test) · not statistical

## The gap this closes

Module 7 would otherwise repeat the project's recurring failure: tables that exist only in Python.
This issue registers the step, writes the tables, adds the `FUNCTION` subworkflow, calls it from
`main.nf`, stages the five setup-manifest resources under both executors, and proves it with a real
Slurm run — before any of M8's golden or calibration work starts.

## Spec sections implemented

- §6.2 step row: "`function` | py | gene table, xrefs, GO/KEGG resources → `category_summary`,
  `enrichment_*` | function".
- §9.1 (`tables/category_summary.tsv`, `enrichment_detection.tsv`, `enrichment_gsea.tsv`; "Optional
  modules that are switched off write no files; the manifest lists them under `modules_off`"),
  §9.2 conventions, §9.3 the three column dictionaries, §9.4 run flags (`bakta_db_light`,
  `layer_missing:<layer>`), §9.5 Frictionless resources.
- §6.3 ("Optional modules. `mobile`, `bgc` and `function` are switched with `params.modules`.
  Downstream steps accept a missing optional layer and record `layer_missing` in the manifest."),
  §6.4 (both executors give identical tables), §7.7.1 ("the manifest notes `bakta_db = light`").

## Public interface

`src/giap/steps/function.py`
```python
@step("function", module="function",
      inputs=["gene_table", "tus"],
      optional_inputs=["go_basic", "goslim_generic", "kegg_ko_pathway", "kegg_pathway_names",
                       "cog2020_fun", "annotation_meta"],
      outputs=["category_summary", "enrichment_detection", "enrichment_gsea"])
def function_step(ctx: StepContext) -> None: ...
```
`src/giap/io/function_tables.py` — three writers; `src/giap/schemas/tables/{category_summary,
enrichment_detection,enrichment_gsea}.json`.

1. **Inputs.** `gene_table` = the detectability output (`genetable/genes.parquet`, as CONTEXT
   takes it); its function columns are the annotation source (M8-03 step 1). `tus` (C5) is
   required for the TU bootstrap and the TU collapse. The five resources are optional root keys
   resolved from the setup manifest exactly like `cog2020_def` (M4-15 decision 5): a new
   `giap.setupdb.installed.resource(name)` generalises `cog_table()`; `cli.py` stages all five for
   `--executor local`, and `main.nf` gains `installed('go_basic', 'go/go-basic.obo')` etc. beside
   `cogTable()` (paths from `db_sources.yaml`).
2. **Order.** collections (M8-03) → `category_summary` + ORA (M8-04) → GSEA (M8-05) → writers.
   RNG streams `function.bootstrap` and `function.gsea` from `giap.rng.for_step`.
3. **Writers.** §9.2: `NA`, `%.6g`, `true`/`false`, `;`-joined `leading_edge` and `flags`, sorted
   by collection then term (then direction for ORA). ADR-0006 encoding. The new `flags` column
   (index decision 12, ADR-A) is the last column of `enrichment_detection.tsv` and
   `enrichment_gsea.tsv`; its values are added to the flag vocabulary resource.
4. **Schemas.** Three Frictionless resources; enums: `collection` per index decision 7, `direction`
   {over, under}, `method` {wallenius_ora, fgsea_simple, fgsea_multilevel}. No coordinates
   (`coordinates: none`). `schemas/datapackage.json` lists them.
5. **Module switch.** `function` not in `cfg.modules` → the step is not planned (local) / not
   called (Nextflow), no files are written, and the module is listed under `modules_off` where the
   manifest is assembled (today: the `function.skipped` event; the manifest step is §6.7 work).
6. **Graceful states** (all exit 0): no resource installed → the three files header-only plus run
   flags `layer_missing:go;layer_missing:kegg;layer_missing:cog`; no abundance →
   `enrichment_gsea.tsv` header-only, `layer_missing:abundance`; `bakta_light` → run flag
   `bakta_db_light` and an event noting `bakta_db = light`.
7. **Nextflow.** `workflow/subworkflows/local/function.nf` in the CONTEXT style (`GIAP_STEP` alias
   `FUNCTION`, `stepIn`, `kv`, `pick`), called from `main.nf` after `GENETABLE`, beside `CONTEXT`,
   only when `'function' in params.modules`. Required keys passed explicitly to `stepIn`
   (`gene_table`, `tus`) — the MVP-01 empty-`need` bug does not apply.
8. **publishDir / local publish.** `tables/*.tsv` published, nothing else (M6-08 decision 7); the
   local executor's published set gains the three keys.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| `--modules core` (function off) | no files; `function.skipped` event; no error |
| no GO/KEGG/COG installed (GitHub CI, a fresh machine) | three header-only files, three `layer_missing` flags |
| no abundance | `category_summary.tsv` written (median abundance NA), `enrichment_gsea.tsv` header-only |
| sim100k (≈ 100 genes, random simulator GO IDs) | sparse tables; most terms below the size floors; **passes** |
| a resource path in the manifest but the file deleted | treated as not installed (the `installed()` / `cog_table()` rule) |
| an upstream channel empty in Nextflow | the step does not run, as CONTEXT |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `giap run --executor local` on sim100k with `$GIAP_DB_DIR` = FarmShare's `local/db` | three files under `tables/`; valid against their schemas; exit 0 |
| 2 | the same through Nextflow (`-profile test`) | byte-identical tables (the §6.4 executor-equality test) |
| 3 | run with `function` removed from `modules` (config or `--set`) | no `tables/category_summary.tsv`; exit 0 |
| 4 | run without `$GIAP_DB_DIR` | header-only files; events `layer_missing:go`, `layer_missing:kegg`, `layer_missing:cog` |
| 5 | run without `--abundance` | `enrichment_gsea.tsv` header-only; `layer_missing:abundance` |
| 6 | `frictionless validate` on the three resources | valid |
| 7 | nf-test stub run of `FUNCTION` | the three output keys emitted and named as in §9.1 |
| 8 | `sbatch scripts/slurm/giap_run.sbatch` on sim100k | completes; the three tables exist (may be sparse) |

## Acceptance criteria

- [ ] `tests/unit/steps/test_function_step.py`: examples 1 (with a toy db dir), 3, 4, 5, 6;
  `tests/integration` adds example 2 beside the existing executor-equality test.
- [ ] `workflow/tests/function.nf.test`: stub coverage in the style of #130, plus a registry test
  that the step's declared keys match the subworkflow's (as `test_nf_annotate_registry.py`).
- [ ] Schemas validate; enums exactly as index decision 7 and §9.3.
- [ ] **Evidence that counts (the MVP rule):** one real `sbatch scripts/slurm/giap_run.sbatch` on
  sim100k with the installed resources that reaches the end and writes the three tables; job ID
  in PROGRESS.md.
- [ ] Runtime on a 4,500-gene simulated genome recorded in `benchmarks/baseline.json`.

## Hardening owed (pass 2)

- A meaningful real-data run: the K-12 full-DB Bakta recording
  (`tests/data/recorded/bakta/1.12.2-full/k12`) with M3-22's in-silico peptides and a simulated
  abundance column, so ORA and GSEA have real GO/KEGG/COG structure (Q3).
- Contract assertions on the three tables inside the pipeline (the MVP-01 pass-2 pattern);
  `-resume` across `FUNCTION`; the real-tool arm under `-profile farmshare,apptainer`.
- Golden tables for the 50-kb fixture once M8-14 gives the simulator real-shaped annotations.
- Manifest: the KEGG retrieval date and the GO release in `run_manifest.json` (ADR-0017 "Every run
  report states the KEGG retrieval date").

## Questions for the lead

1. Resource staging: generalise `cog_table()` into `installed.resource(name)` and stage five more
   roots (proposed), or pass `--db-dir` into the step and let it read the manifest itself?
   **Proposed default:** roots — the step stays a pure function of declared inputs (§6.2).
2. §6.2 lists `xrefs` as an input; the gene table already carries the parsed `go_terms`,
   `kegg_ko`, `cog_category`. **Proposed default:** read the gene table only; `xrefs` is not a
   step input (one source of truth, no second parser).
3. Is a sparse sim100k result an acceptable pass-1 definition of done, with the K-12 run in pass 2?
   **Proposed default:** yes (the M6-08 precedent).
