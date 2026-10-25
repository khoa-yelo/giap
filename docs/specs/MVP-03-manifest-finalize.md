# MVP-03: the §6.2 `manifest` and `finalize` steps (run_manifest.json, provenance, finalising OUT)

Track: **MVP / end-to-end**. Scheduled before M9-16 (M9-index Q2, lead 2026-10-22), because M9's
provenance section (M9-16 §6), `giap compare`'s genome check (M9-22 item 4) and the KEGG retrieval
date (ADR-0017) all read the manifest.
Depends on: nothing beyond main (cf41ca4). Not statistical; no tools; no new dependency
(`jsonschema` 4.26 and `frictionless` 5.20 are already in the venv).
Estimate: **≈ 760 code lines, split into 03a / 03b / 03c** (see "Size and split"; ADR-0007).

Sources: GIAP_SPEC.md §6.2 (rows `finalize`, `manifest`), §6.3 (Launch, Error strategy, Optional
modules), §6.4, §6.6, §6.7, §5.1 (`--publication`), §5.2, §7.1 "Calibration" (terminator cut-off
"flagged `uncalibrated` in the manifest"), §7.3.2 (tier "recorded in the manifest"), §7.7.1
(`bakta_db = light`), §7.9.3, §9.1 (`modules_off`), §9.2 (Versioning), §9.4 (Run row), §9.5;
ADR-0014, ADR-0017, ADR-0021 (Amendment 1), ADR-0036 §2, ADR-0037 rung 2 (on origin/adr/m8-amendments);
backlog lines on M2-15, M4-18, M4-19 d7, MVP-01 pass 2 (sensitivity off), M5-09 d3 (genes.tsv),
M6-06/M6-08 (`bin_counts`, π₀, `±Infinity`), M7-03 review (descriptor location);
M9-index decision 26, M9-02 items 6 and 8, M9-16 §6, M9-22 item 4.

## The gap
- `src/giap/schemas/run_manifest.schema.json` is the M0 skeleton. **Nothing writes
  `run_manifest.json`.** No step is registered under `manifest` or `finalize`.
- Run flags are scattered. Some go to `logs/events.jsonl` (`*.layer_missing`,
  `mapping.input_prefiltered`). Some sit in `annotation_meta.json` (`anti_sd_default`,
  `bakta_db_light`). Some exist only as log lines and never reach a file: `c_assumed` and
  `operon_source_distance_model` (both in the §9.4 Run row), the `layer_missing:*` warnings of
  `region_layer`, `operons` and `context_stats`, and `cai_reference_missing`.
- `context_stats.oriter.multiple_oric` is logged on `giap.context.oriter`. That logger is outside
  `giap.steps`, so `EventsHandler` never sees it and the event is lost under both executors.
- `operons.missed_operons` writes `-Infinity` (not strict JSON). Backlog M6-08 marks this
  "must-do before the manifest writer".
- `logs/events.jsonl` is opened in append mode, and nothing resets it. A re-run or a `-resume` into
  the same OUT appends the preflight events and every cached step's events a second time.
- `tables/genes.tsv` is written by `context_stats` (M5-09 decision 3). A run that skips CONTEXT
  has no master table. §6.2 gives the file to `finalize`.
- `OUT/schemas/datapackage.json` is never written. `src/giap/schemas/datapackage.json` lists 3 of
  the 21 table schemas. Its `tables/...` paths would resolve under `OUT/schemas/` (M7-03 review).

## Goal
Every `giap run` that gets past preflight writes `OUT/run_manifest.json`. That holds for both
executors, for success and for a failed step. The manifest validates against
`src/giap/schemas/run_manifest.schema.json` and holds the §6.7 list. Two manifests of the same inputs
and seed are byte-identical once the declared volatile fields are removed: local against Nextflow,
and a fresh run against a `-resume`. `tables/genes.tsv` always exists when a gene table exists, and
`OUT/schemas/datapackage.json` describes exactly the tables that this run wrote.

## Design

### D-A. Where the two steps run
| Step | Registered | Runs | Why |
|---|---|---|---|
| `finalize` | py, core: `gene_table` (required), `gene_context` (optional) → `genes_tsv` | **in the graph**, after GENETABLE/CONTEXT (Nextflow: `GIAP_STEP` in a new `finalize.nf`; local: topological order) | M9's `plots`/`report` read `tables/genes.tsv` in the graph (M9-02 item 8) |
| `manifest` | py, core, run-tail only: `preflight` (required); optional `provenance`, `events`, `trace`, `annotation_meta`, `sensitivity`, `island_diagnostics` → `run_manifest` | **after the graph**, called by `giap run` under both executors, also after a failed step | the trace is complete only after Nextflow exits; a failed run needs a manifest; it must be the last writer |
| datapackage writer (the descriptor half of §6.2 `finalize`) | a function of the tail, not a separate step | after the graph, before `manifest` | only then is the set of published tables known; it is executor-identical by construction |

`giap _step manifest ...` still works, because the step is registered. The local runner leaves
run-tail steps out of the topological plan (`TAIL_STEPS = ("manifest",)`) and calls them in
`run_tail()`.

### D-B. Run order of `giap run` (both executors)
1. Rotate `OUT/logs/events.jsonl` to `events.previous.jsonl` (one generation kept).
2. Preflight (unchanged). An exit of 1, 2 or 3 here writes **no** manifest: `preflight.json` is the
   record.
3. Write `OUT/work/provenance.json`, the **provenance core** (D-D), before the launch.
4. Run the graph (local loop or `nextflow run`). `append_step_events` runs as today.
5. Run tail: write the datapackage (03b), then the manifest. This runs even when step 4 failed.
6. Exit code: 4 is kept if step 4 failed. If the graph succeeded but the tail failed, the exit is 4
   with the message `step manifest failed; log: OUT/logs/manifest.log`. A tail error after a failed
   graph is logged and never hides the original failure.

### D-C. What the manifest gathers, and from where
| Manifest field | Source | Issue | Volatile |
|---|---|---|---|
| `schema_version` | `cfg.schema_version` (the §9.2 output-standard version) | a | |
| `giap_version`, `git_commit`, `git_dirty` | `giap.__version__`; `git -C <source root> rev-parse HEAD` / `status --porcelain`; `null` when not a git checkout | a | |
| `command_line` | `sys.argv` | a | yes |
| `executor`, `profile`, `resume`, `outdir`, `started_at`, `finished_at` | `giap run` | a | yes |
| `status` ∈ {`completed`, `failed`}, `exit_code`, `failed_step`, `failed_log` | runner result; Nextflow: `failed_step(trace)`; local: `StepError.name` | a | `failed_log` |
| `publication` | the `--publication` flag (it is not in `GiapConfig`) | a | |
| `config` | `cfg.model_dump(mode="json")` | a | the `annotation.precomputed.*` paths |
| `inputs[]` {role, file, path, sha256, size} | `preflight.json` `inputs` (reused, never re-hashed); the manifest hashes the other given files itself: `annotation`, `reference_peptides`, `replicon_depth`, and each set `annotation.precomputed.<key>` as role `precomputed:<key>` | a | `path` |
| `inputs_absent[]` | optional `giap run` inputs that were not given (`proteins`, `abundance`, `annotation`, `reference_peptides`, `replicon_depth`); backlog M4-18 | a | |
| `tools[]` {name, version, container_digest, conda_lock_hash, processes} | Nextflow: `<workdir>/versions.yml` of each COMPLETED/CACHED trace row (no workflow change); the digest comes from the trace column `container` (added to `trace.fields`); the lock hash is the SHA-256 of `workflow/envs/locks/<tool>.explicit.txt` when the task had no container. `giap` itself is left out (it is `giap_version`). Local: `[]` | a | |
| `databases[]` {name, version, checksum, retrieved, size_bytes} | `<db_dir>/db_manifest.json`, items with `status: ok` (KEGG's `retrieved` is the ADR-0017 date); `db_dir` absent → `[]` | a | |
| `annotation` {source, bakta_version, db_version, db_type} | `annotation_meta.json` (C4; §7.7.1 `bakta_db = light`) | a | |
| `detection_tier` | the tier the run used: `genes.tsv` `detection_tier` = `detection_summary.tsv` `tier` (§7.3.2: one tier per run, the same in the manifest and `genes.tsv`). Amended 2026-10-25 (lead, #363): it was preflight's provisional `evidence_tier`, which differs since M4-24 (unscored decoys fall back to `mayu` / `two_peptide`). Preflight's tier stays in `preflight.json`; it is the fallback only when the run wrote neither table (failed or stub runs, a manifest built without tables) | a | |
| `seeds` {`master`, `<step>`: crc32(step)} | `cfg.seed`; one entry per executed py step. The stream of a step is `SeedSequence([master, seeds[step]])` (§6.6). Sub-streams such as `context_stats.islands` are not listed | a | |
| `steps[]` {name, status, n_tasks, n_retries, wall_time_s, cpu_time_s, peak_rss_bytes} | Nextflow: `trace.tsv` grouped by `_step_name` (sum of `realtime`, sum of realtime × %cpu, max `peak_rss`, binary units as in `scripts/measure_runtimes.py`, last attempt only). Local: `perf_counter` / `process_time` around each `run_step`, no RSS | a | yes (whole array) |
| `skipped[]` {name, reason} | registered steps minus executed steps. `reason` ∈ {`module_off:<m>`, `opt_in_off`, `precomputed:<key>`, `not_run`}; the same rule under both executors | a | |
| `modules_off[]` | {mobile, bgc, function} − `cfg.modules` (§9.1) | a | |
| `flags[]` | sorted unique union of: `data.flag` / `data.flags` of every event; `annotation_meta.flags`; `uncalibrated:<key>` (03c) | a (03c adds the calibration source) | |
| `layers_missing[]` | `<layer>` of each `layer_missing:<layer>` flag, minus layers that an off module explains (decision 7) | a | |
| `warnings[]` | `preflight:<rule_id>` for each `warn` check; `<step>:<event>` for each WARNING-level event; sorted unique | a | |
| `volatile[]` | the JSON pointers above (constant `giap.manifest.VOLATILE`) | a | |
| `provenance_digest` | SHA-256 of the canonical JSON of the provenance core (D-D) | a | |
| `sensitivity` {level, n_rows, families: {family: {variant: {n_primary, n_variant, n_shared, jaccard}}}} | `cfg.output.sensitivity` and `tables/sensitivity.tsv`. "Absent" and "header-only" both give `n_rows: 0` (backlog MVP-01 pass 2, M4-19 d7) | c | |
| `context.replicons[]` {contig, scanned, block_rule, l_hat, block_len, morans_i_lag1, flags} | `context/island_diagnostics.parquet` (ADR-0021 Amendment 1) | c | |
| `context.oric[]` | `context_stats.oriter.multiple_oric` events | c | |
| `operons` {pi0, bin_counts} | the `operons.missed_operons` event, with ±∞ → `null` (backlog M6-06/M6-08) | c | |
| `mapping` {search_db, n_queries} | the `mapping.search_db_resolved` event | c | |
| `calibration` {key: {calibrated, provenance}} | `src/giap/data/calibration.yaml` through `giap.data.calibration.load()` (§7.1, backlog M2-15) | c | |

Canonical JSON throughout: `json.dumps(sort_keys=True, indent=2, ensure_ascii=False,
allow_nan=False)` plus a trailing newline. A NaN or ±∞ reaching the writer is a bug (exit 4), not a
`NaN` token.

### D-D. Provenance core, stable view, digest
- **Provenance core:** `schema_version`, `giap_version`, `git_commit`, `git_dirty`, `publication`,
  `config` (each `annotation.precomputed.*` value replaced by its file name), `inputs` (without
  `path`), `databases`, `seeds.master`. All of it is known before launch, so `giap run` writes it to
  `OUT/work/provenance.json` (D-B step 3). M9-16 can then stage the core into `report` and print
  the same digest that the manifest records later.
- **Stable view:** the manifest minus `VOLATILE`, with the same path normalisation. Two
  manifests compare equal on it across executors and across `-resume`.
- `provenance_digest = sha256(canonical(core))`. The digest covers neither `tools` (they follow from
  `git_commit` and the profile), nor `steps`, nor flags.

### D-E. Run flags that reach the manifest (03a converts the log-only ones to events)
| Token | §9.4 / ADR | Emitted by | Reaches a file today | 03a change |
|---|---|---|---|---|
| `layer_missing:<layer>` | §9.4 Run | `context_layers`, `gene_table`, `novel_orfs` | event | none |
| ditto (`genomad`, `antismash`) | §9.4 Run | `region_layer` (`steps/regions.py:111,115`) | log only | add an event |
| ditto (`promoters`, `terminators`, ...) | §9.4 Run | `operons` (`steps/operons.py:196`) | log only | add an event |
| `layer_missing:features` | §9.4 Run | `context_stats` (`steps/context_stats.py:124`) | log only | add an event |
| `anti_sd_default`, `bakta_db_light` | §9.4 Run | `parse_bakta` → `annotation_meta.json` | file | none (read the file) |
| `c_assumed` | §9.4 Run | `placement` (`steps/mapping.py:336`) | log only | add an event |
| `operon_source_distance_model` | §9.4 Run | `operon_pairs` (`steps/operons.py:112`) | log only | add an event |
| `input_prefiltered` | **not in §9.4** | `placement` | event | none; Khoa item 1 |
| `cai_reference_missing` | **not in §9.4** | `codon_usage` (`steps/codon.py:70`) | log only | add an event; Khoa item 1 |
| `no_detection_limit_model` | **not in §9.4** (it is a `pairs` row token) | `operons` (info log) | log only | add an event; Khoa item 1 |
| `tu_recal_skipped`, `too_few_ranked_tus`, `no_eligible_genes`, `p_detect_clipped`, `null_degenerate` | ADR-0036 §2 (accepted-by-khoa) | `function` (`function.run_flags` event, M8-06 branch) | event once merged | none |
| `tu_recal_reduced` | ADR-0037 (accepted-by-khoa) | `function` | event once merged | none |
| `uncalibrated:<key>` | §7.1 says "`uncalibrated` in the manifest"; the parameterised form is new | the manifest, from calibration.yaml (03c) | — | 03c; Khoa item 1 |

The event shape is the existing one (`extra={"event": "<step>.run_flag", "data": {"flag": tok}}`;
the `*.layer_missing` events keep their names). The manifest reads any event whose `data.flag` is a
string or whose `data.flags` is a list of strings. That covers ADR-0036's `function.run_flags`
without a special case. `giap.manifest.RUN_FLAGS` lists the registered vocabulary, including the
`layer_missing:` and `uncalibrated:` prefixes. Tokens pending Khoa are in `RUN_FLAGS_PENDING`.

## Examples (input → expected output)
Examples 1–8 are pure-builder unit tests (`build_manifest(sources) -> dict`, no file I/O).
Examples 9–12 are runs.

| # | Input | Expected |
|---|---|---|
| 1 | events: the same `{"step":"gene_table","event":"gene_table.layer_missing","level":"WARNING","data":{"layer":"cog","flag":"layer_missing:cog"}}` line twice (a resume duplicate) | `flags == ["layer_missing:cog"]`, `layers_missing == ["cog"]`, `warnings == ["gene_table:gene_table.layer_missing"]` |
| 2 | `annotation_meta = {"db_type":"light","flags":["anti_sd_default","bakta_db_light"],"bakta_version":"1.11.4","db_version":"6.0",...}` | `flags ⊇ {"anti_sd_default","bakta_db_light"}`; `annotation.db_type == "light"` |
| 3 | `cfg.modules = [core, mobile]`; events contain `layer_missing:antismash` (from `context_layers`) | `modules_off == ["bgc","function"]`; `layers_missing` has no `antismash`; `flags` has no `layer_missing:antismash` (decision 7) |
| 4 | one trace row `ANNOTATE:PARSE_BAKTA (parse_bakta)`, COMPLETED, realtime `5.5s`, %cpu `65.8%`, peak_rss `269.9 MB` | `steps == [{"name":"parse_bakta","status":"completed","n_tasks":1,"n_retries":0,"wall_time_s":5.5,"cpu_time_s":3.619,"peak_rss_bytes":283010662}]` |
| 5 | trace: `PROMOTECH (chunk_1)` FAILED attempt 1 (exit 137), then `PROMOTECH (chunk_1)` COMPLETED attempt 2, plus `PROMOTECH (chunk_2)` COMPLETED | one entry `promotech`: `n_tasks 2`, `n_retries 1`, `status completed`, metrics from the attempt-2 row and chunk_2 only |
| 6 | `cfg.seed = 7`; executed py steps `detectability`, `gene_table` | `seeds == {"master":7,"detectability":3987542110,"gene_table":393846822}` |
| 7 | the event `operons.missed_operons` with `lower_bp: -Infinity` in bin 0 (03c) | `operons.bin_counts[0].lower_bp is None`; `json.dumps(..., allow_nan=False)` succeeds |
| 8 | `cfg.output.sensitivity = "off"`, sensitivity table header-only **or** absent (03c) | `sensitivity == {"level":"off","n_rows":0,"families":{}}` in both cases |
| 9 | **local**, SIM100K with every precomputed layer, `--seed 7 --modules core,mobile,bgc` (the `test_m2_21_local_vs_nextflow` setup) | exit 0; `run_manifest.json` validates; `status completed`; `executor local`; `modules_off == ["function"]`; `tools == []`; `inputs[role=genome].sha256` equals `preflight.json`'s; `inputs` has a `precomputed:bakta_json` row; `provenance_digest == sha256(canonical(work/provenance.json))` |
| 10 | the same arguments with `--executor nextflow --profile farmshare_local` | `stable_view(nextflow) == stable_view(local)` (this is §6.4 extended to the manifest); `steps[*].peak_rss_bytes` present; `tools == []` (all layers precomputed; `giap` left out) |
| 11 | the M2-21 fault: PROMOTECH exits 137 twice (stub) | exit 4; manifest written; `status failed`, `exit_code 4`, `failed_step "promotech"`, `failed_log` ends in `.command.log`; `steps[promotech].status == "failed"`; no `finalize` entry in `steps` |
| 12 | example 11, then `--resume` without the fault | exit 0; `status completed`, `resume true`; the steps that finished in run 1 have `status cached`; `stable_view` equals that of a fresh run with the same arguments; `logs/events.jsonl` holds each step's events once (the run-1 file is in `events.previous.jsonl`) |

Finalize (03b):

| # | Input | Expected |
|---|---|---|
| 13 | `finalize` with `gene_table` and `gene_context` | `tables/genes.tsv` byte-identical to what `context_stats` writes today (golden: the current `context_stats` output on the M5-09 fixture) |
| 14 | `finalize` with `gene_table` only (CONTEXT skipped: the "precomputed regions layer" nf-test before MVP-02, or no C7/C8) | `tables/genes.tsv` exists, `island_id`/`replichore`/`replichore_position` are `NA` in every row, the columns follow `GENES_COLUMNS` order |
| 15 | the tail on an OUT with `tables/{genes,pairs,sensitivity}.tsv` and `mapping/peptide_hits.tsv` only | `OUT/schemas/datapackage.json` lists exactly those 4 resources (sorted by name), each `schema` copied to `OUT/schemas/<group>/<name>.schema.json`; `frictionless.Package(descriptor, basepath=OUT).validate().valid` |
| 16 | `--modules core` (function off) | no `category_summary`/`enrichment_*` resources and no files; `modules_off == ["bgc","function","mobile"]` (sorted) |

## Edge cases
| Case | Behaviour |
|---|---|
| A step skipped (CONTEXT without C7/C8; OPERONS without inputs) | `skipped[] {name, reason: "not_run"}`; its tables are missing from the datapackage; no flag invented |
| An external step replaced by a precomputed layer | `skipped[] {name, reason: "precomputed:<key>"}`; the layer file is in `inputs` as `precomputed:<key>` with its SHA-256 |
| A missing optional layer, module on (no COG table installed) | `layer_missing:cog` in `flags` and `cog` in `layers_missing`, and a `gene_table:gene_table.layer_missing` warning |
| A module off (`--modules core`) | `modules_off` sorted; the module's steps are `skipped` with `module_off:<m>`; its layer tokens are left out of `layers_missing` (decision 7); no files for it (§9.1) |
| Opt-in off (`bactermfinder`) | `skipped {bactermfinder, opt_in_off}`; not under `modules_off` |
| A failed run (exit 4) | the manifest is written with `status failed`, and the flags and steps known so far; `finalize` may be absent; the datapackage lists only what was published |
| Preflight error, exit 1/2/3 | no manifest and no provenance core; `preflight.json` is the record |
| Exit 5 (internal error) | best effort: the tail is attempted; if it raises too, nothing is written and the exit stays 5 |
| The tail fails after a successful graph | exit 4 `step manifest failed`; `logs/manifest.log` holds the traceback |
| `-resume` | trace rows CACHED → step `status cached` with the original metrics; `resume: true`; flags/warnings from a fresh `events.jsonl` (D-B step 1) are identical to a non-resumed run |
| Re-run of the local executor into the same OUT | same rotation; `steps` from this run only |
| No `--db-dir` and no `$GIAP_DB_DIR` | `databases == []`; the KEGG date is not recorded; M9-16 prints "not recorded" |
| `db_manifest.json` in the provisional `{"<name>": ...}` form (read_manifest's second shape) | `databases[]` gets `{name, version}`, with `checksum`/`retrieved` `null` |
| `trace.tsv` absent (Nextflow died before the first task) | `steps == []`; `status failed`; `failed_step null` |
| A trace metric `-` or empty | that metric is left out of the step entry (the schema makes them optional), never written as 0 |
| `versions.yml` missing in a task dir (failed task, `work/` cleaned) | that tool is left out; no error |
| Two processes report different versions of one tool | one `tools[]` row per (name, version), each with its `processes[]` |
| An event line that is not JSON | the manifest step fails (exit 4): events.jsonl is GIAP's own output, so a bad line is a bug |
| Non-git install (wheel) | `git_commit: null`, `git_dirty: null` |
| Circular/linear replicons, minus strand | do not apply to the manifest. `context.replicons[]` (03c) keeps C1 order. `genes.tsv` (03b) keeps the §6.5 sort |
| Empty inputs (a header-only gene table) | `finalize` writes a header-only `genes.tsv`; the manifest has the `no_eligible_genes` flag only if `function` ran |
| Ties | flags, warnings, layers, `modules_off`, `skipped`, `databases`, `tools` and datapackage resources are all sorted by name; `steps` in trace order (Nextflow) or plan order (local) |

## Acceptance criteria
03a (manifest):
- [ ] `run_manifest.schema.json` extended: every D-C field of 03a, `additionalProperties: false` at
  the top level, `flags` items matching `^[a-z0-9_]+(:[a-z0-9_.-]+)?$`. Existing
  `tests/unit/test_schemas.py` stays green.
- [ ] `tests/unit/manifest/test_build.py`: examples 1–6 against the pure builder, plus every
  edge-case row that does not need a run.
- [ ] `tests/unit/manifest/test_stable_view.py`: `stable_view` removes exactly `VOLATILE`;
  `provenance_digest` is stable under key order and changes when one input SHA-256 changes.
- [ ] `tests/unit/manifest/test_vocabulary.py`: every token emitted anywhere in `src/giap`
  (constants `*_FLAG`, `LAYER_MISSING`, …) is in `RUN_FLAGS ∪ RUN_FLAGS_PENDING`.
- [ ] `tests/unit/steps/test_run_flag_events.py`: each D-E "add an event" site emits an event with
  `data.flag`.
- [ ] Examples 9, 10 (integration, `needs_nextflow`), 11 and 12 (stub variants of
  `test_m2_21_resume_failure.py`): the manifest validates and the stated fields hold. Example 10
  compares `stable_view`s.
- [ ] Evidence: one Slurm `giap_run.sbatch` on sim100k (real path), job ID in PROGRESS.md; the
  manifest validates and `tools[]` lists the containerised tools with digests.

03b (finalize):
- [ ] Examples 13–16. `context_stats` no longer declares `genes_tsv`. The test-author updates the
  tests that pin it: `tests/integration/test_context_stats_step.py`, `context.nf.test`,
  `main_genetable_context.nf.test`. `workflow/tests/finalize.nf.test` is the stub for FINALIZE.
- [ ] `src/giap/schemas/datapackage.json` becomes the full catalogue (every table and annotation
  schema). A test checks that every `schemas/{tables,annotation}/*.schema.json` has exactly one
  resource.
- [ ] Backlog lines closed: M5-09 d3 (genes.tsv when CONTEXT is skipped) and the M7-03 descriptor
  location.

03c (content blocks):
- [ ] Examples 7, 8; `context.replicons` equals `island_diagnostics.parquet` row for row; the
  `oriter` event reaches `events.jsonl` under both executors; `calibration` lists both
  calibration.yaml keys and adds `uncalibrated:operon_distance` and
  `uncalibrated:terminatornet_yes` while both are `calibrated: false`.
- [ ] Backlog lines closed: M2-15 (manifest `uncalibrated`), M4-18 (missing `--replicon-depth`,
  through `inputs_absent` in 03a), M4-19 d7, MVP-01 pass 2 (sensitivity off), and M6-06/M6-08
  (`±Infinity`, the manifest half).

## Test ideas (beyond the examples)
- Property: `build_manifest` is invariant under any permutation of the events lines (Hypothesis).
- Round trip: `json.loads(write(manifest)) == manifest`; the writer output ends in `\n` and has no
  `NaN` or `Infinity`.
- Mutation guard: drop one entry from `VOLATILE` → example 10 fails (it proves that the list is
  needed).
- `giap _step manifest --in preflight=... --out run_manifest=...` works on its own (registry path).

## Size and split (ADR-0007: code under src/, workflow/, scripts/)
| Part | Contents | Code lines |
|---|---|---|
| **03a** manifest | `giap/manifest.py` (pure builder, `VOLATILE`, `stable_view`, digest, `RUN_FLAGS`) ≈ 170; `giap/io/manifest_sources.py` (trace units, versions.yml, events, db manifest, hashing, strict writer, git) ≈ 90; `steps/manifest.py` + `run_tail` ≈ 35; `cli.py` (events rotation, provenance core, tail under both executors and on failure, local step timing) ≈ 45; the 7 event conversions of D-E ≈ 20; `trace.fields += container` 1; schema JSON ≈ 80 | **≈ 440** |
| **03b** finalize | `steps/finalize.py` (moves `_with_context`) ≈ 45; the `context_stats` key removal ≈ 10; `finalize.nf` + `main.nf` ≈ 35; datapackage writer + schema copy in the tail ≈ 70; `datapackage.json` catalogue (JSON) ≈ 30 | **≈ 190** |
| **03c** content blocks | sensitivity / context / oric / operons / mapping / calibration blocks ≈ 80; the `oriter` logger fix and `_plain` ±∞ → null ≈ 10; schema JSON ≈ 40 | **≈ 130** |

03a is above 400 only because of its schema JSON. If the size check fails it, move the four
`tools`/`databases` readers (≈ 40) to 03c. Order: 03a, then 03b and 03c in parallel (they touch
different files apart from the schema; 03c merges main after 03b). M9-16 can start once 03a has
merged.

## Hardening owed (pass 2)
- Validate every published table against its schema inside the tail (Frictionless) and fail with
  exit 4. Pass 1 validates in tests only.
- Executor equality for the step set: pass 1 asserts `stable_view` equality on the
  all-precomputed sim100k case only; real-tool runs compare `skipped`/`steps` names (recorded).
- A per-run `run_id` in every event line, so mixed generations can be told apart without rotation.
- `peak_rss_bytes` under the local executor (process high-water is not per step; needs a
  subprocess per step or `tracemalloc`).
- The Nextflow `-with-report`/timeline HTML are not in the manifest; link them if M9 asks.
- Runs on peptide-free inputs (proteins only) never reach GENETABLE in `main.nf` (`if
  (inputs.peptides)`), so they get no `genes.tsv`; this is an existing gap and is out of scope.

## Proposed defaults (lead decisions)
1. **Proposed default:** split into 03a / 03b / 03c as in the size table; 03a first, because it
   unblocks M9-16 and M9-22.
2. **Proposed default:** the `manifest` step is registered but runs only in the run tail, under
   both executors and after a failed step, never as a Nextflow process (D-A, D-B). This departs from
   §6.3's "one process per Python step"; a lead-level ADR records it (no default, threshold or
   output changes).
3. **Proposed default:** `finalize` is split. The in-graph step owns `tables/genes.tsv` (this
   resolves M5-09 d3; `context_stats` drops the `genes_tsv` key). The descriptor half
   (`schemas/datapackage.json` plus schema copies) runs in the tail, because only then is the set
   of published tables known. The same ADR as decision 2 records this departure from the §6.2 row.
   Alternative: run `finalize`, `plots`, `report` and `manifest` all in the tail. That is simpler and
   executor-identical by construction, but it overturns M9-02's accepted `report.nf` design.
4. **Proposed default:** the datapackage stays at §9.1's `OUT/schemas/datapackage.json`. Resource
   `path`s and `schema`s are relative to `OUT`, and GIAP validates with
   `frictionless.Package(descriptor, basepath=OUT)`. Checked in the session with frictionless
   5.20: a `../tables/x.tsv` path is rejected as "not safe" even with `trusted=True`, and a naive
   `frictionless validate OUT/schemas/datapackage.json` resolves under `OUT/schemas/`. The user
   guide gives the `--basepath OUT` command. The alternative is Khoa item 2.
5. **Proposed default:** `giap run` rotates `logs/events.jsonl` to `events.previous.jsonl` at start.
   This changes today's append behaviour; the test-author checks `test_mapping_nf_run_events.py`.
6. **Proposed default:** `provenance_digest` (D-D) is the "manifest digest" that M9-16 prints, not
   the SHA-256 of the file. The file hash changes on every run (timings) and does not exist yet
   when the in-graph `report` runs. `giap run` writes `OUT/work/provenance.json` before launch so
   that M9-16 can stage it. M9-16 §6 needs a one-line amendment.
7. **Proposed default:** a `layer_missing:<layer>` token that an off module explains (`genomad` ↔
   `mobile`, `antismash` ↔ `bgc`) is left out of `flags` and `layers_missing`. The module is listed
   under `modules_off` instead. `cog`, `go` and `kegg` are always kept: `gene_table` reads COG in
   core.
8. **Proposed default:** per-step metrics are task sums (`wall_time_s` = Σ realtime, task-seconds,
   not elapsed time), `peak_rss_bytes` = max, last attempt only, with `n_tasks`/`n_retries`. The
   trace is parsed in its human form with the `scripts/measure_runtimes.py` rules, moved into
   `giap.io.manifest_sources`; `trace.raw` stays off.
9. **Proposed default:** `tools[]` comes from each task's `versions.yml` in its work dir (through
   the trace), with digests from a new trace column `container`. No subworkflow wiring is needed;
   ANNOTATE and NOVEL emit no `versions` channel today.
10. **Proposed default:** `databases[]` lists every `ok` item of the db manifest, not only those
    the run used. "Used" is not tracked anywhere today.
11. **Proposed default:** `seeds` = `master` plus crc32 per executed py step (example 6); the schema
    keeps integer values. Named sub-streams are not listed.
12. **Proposed default:** a manifest is written for every run that reached D-B step 3; a preflight
    exit (1/2/3) writes none.
13. **Proposed default:** the schema file keeps its name `run_manifest.schema.json`. §9.5 says
    `manifest.schema.json`, and the repo and `test_schemas.py` use the former. The lead's ADR from
    decision 2 notes the name.
14. **Proposed default:** `schema_version` in the manifest is `cfg.schema_version` (the §9.2
    output-standard version). The JSON Schema `$id` is not versioned separately in pass 1.

## Needs Khoa
1. **Run-scope tokens beyond §9.4.** `input_prefiltered`, `cai_reference_missing`,
   `no_detection_limit_model` (as a run flag) and the parameterised `uncalibrated:<key>` (§7.1 says
   "`uncalibrated` in the manifest", and §9.4 lists `uncalibrated` only under Region and context).
   ADR-0036 treated additions to the §9.4 Run row as Khoa-level. One ADR, "run tokens beyond the
   §9.4 Run row", covers all four. Until then they sit in `RUN_FLAGS_PENDING`: the schema pattern
   accepts them and the vocabulary test names them as pending. Nothing blocks.
2. **Only if the lead rejects decision 4:** moving the descriptor to `OUT/datapackage.json`
   (Frictionless's convention, so a naive validate works) changes the §9.1 output tree in
   GIAP_SPEC.md (protected) and needs an ADR accepted by Khoa.
3. **No protected file is touched otherwise.** There are no new config keys (`defaults.yaml`
   unchanged) and no calibration or tolerance changes. The manifest schema is additive, and
   `src/giap/schemas/` is not in `ci/protected.txt`.

## Questions for the lead (gaps in the spec)
1. §6.7 "warnings": should these be the human messages, or the machine IDs proposed here
   (`preflight:<rule_id>`, `<step>:<event>`)? Messages carry paths, which breaks executor
   equality. **Proposed default:** IDs.
2. §6.7 says that per-step metrics come "from the Nextflow trace", and is silent on the local
   executor. **Proposed default:** local records wall and CPU time, no RSS (Hardening owed).
3. §7.9.3 needs the genome SHA-256 "in each manifest". **Proposed default:** `inputs[role =
   genome].sha256` is the raw input file's hash (as `preflight.json` has it), not that of
   `genome.norm.fa`. M9-22 reads that field.
4. Should `finalize` also fold in the §9.3 `genes.tsv` columns of later modules (variants, novel
   ORFs, function) as optional inputs? **Proposed default:** yes, each producing issue adds its own
   optional input to `finalize`; MVP-03b wires only `gene_context`.

## Lead decisions (2026-10-22)

- Decisions 1–14: accepted as proposed. The split is 03a (manifest), then 03b (finalize) and
  03c (content blocks) in parallel.
- Decisions 2 and 3 (the manifest runs after the graph; finalize is split between graph and tail)
  are recorded in a lead-level ADR written in 03a.
- The four spec-gap questions take their proposed defaults.
- Needs Khoa: the run-scope tokens `input_prefiltered`, `cai_reference_missing`,
  `no_detection_limit_model` and `uncalibrated:<key>` (one ADR). They stay in `RUN_FLAGS_PENDING`
  until Khoa answers; nothing is blocked.

## Lead decisions on the 03c tests (2026-10-22, after 0bf8ee2a)

15. The strict `MVP03C` xfail marks in the 03a tests are removed by the lead once 03c's
    implementation lands (they will XPASS); the implementer reports which ones pass.
16. Under the local executor the CLI tests use `fake_giap_root` for a setup-form db; the sim100k
    integration keeps the provisional test db (`retrieved` null there).
17. `operons` and `mapping` are null when no event exists; `context` holds empty lists. A stub
    `mapping.search_db_resolved` event with no values gives `mapping = null`.
18. `context.oric` is sorted by contig; `context.replicons` keeps C1 order.
19. The `calibration` source defaults to None (no `uncalibrated:` flags) so the 03a exact-flag
    tests stay green; the run tail and `giap _step manifest` pass `calibration.load()`.
20. **Lead edits after the 03c implementation (2bb97710):**
    - The 11 strict `MVP03C` xfail marks are removed (decision 15).
    - `test_operons_step.py::test_decision_6_event_records_pi0_and_bin_counts` compares a
      non-finite oracle edge (the −∞ first bin) as null, since `_plain` maps ±∞ to None.
    - A 0-byte output (a `-stub-run` `touch`) reads as absent.
    - Of several distinct `operons.missed_operons` / `mapping.search_db_resolved` payloads, the
      last one in file order wins.
    - `giap _step manifest` still writes `databases: []`; this is a backlog item.
21. **Review on #236 (lead fixes):**
    - B1: an unreadable `versions.yml` (YAML or UTF-8 error) is skipped with a warning, so the run
      tail never fails on it.
    - B2: tests cover the digest winning over a lock, and BaseLoader keeping `1.10` a string.
    - Decision 20 is now honoured for distinct payloads A, B, A (A wins).
    - B3: the backlog lines this issue closes are marked.
    - Tests are in tests/unit/manifest/test_manifest_review_236.py; nits go to the backlog.
