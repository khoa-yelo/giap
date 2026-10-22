# M6-08: The `operons` step, its three tables, and wiring them into `giap run`

Milestone M6 · track T (tables and integration) · **pass 1 — this is the integration milestone** ·
depends on: M6-04, M6-05, M6-06, M6-07 (and merged MVP-01 wiring) · estimate ≈ 320 code lines
(+ ≈ 350 test lines, + nf-test) · not statistical

## The gap this closes

M4 and M5 built tables that `giap run` reaches only because MVP-01 wired them. Module 4 would
repeat the project's recurring failure if its three tables lived in Python only. This issue
registers the step, writes the tables, adds the `OPERONS` subworkflow, calls it from `main.nf`, and
proves it with a real Slurm run — before any of M6's hardening work starts.

## Spec sections implemented

- §6.2 step row: "`operons` | py | pairs, gene table → `pairs`, `operons`, `operon_flags` | core".
- §9.1 output tree (`tables/pairs.tsv`, `tables/operons.tsv`, `tables/operon_flags.tsv`), §9.2
  conventions (TSV, `NA`, 1-based closed coordinates in TSV, `;`-separated flags), §9.3 the three
  column dictionaries, §9.4 the Operons flag scope, §9.5 Frictionless resources.
- §6.3/§6.4 (both executors give identical tables), §6.5 (ordering, coordinate conversion in
  writers only), §5.3 module flags (Module 4 is core; it is skipped when there is no abundance).

## Public interface

`src/giap/steps/operons.py` (extends the merged module)
```python
@step("operons",
      inputs=["gene_table", "genes", "replicons", "tus", "pairs_raw", "features"],
      optional_inputs=["promoter_windows", "terminators", "tx_rates", "peptides",
                       "peptide_hits", "peptide_runs", "detection_summary"],
      outputs=["pairs", "operons", "operon_flags"])
def operons_step(ctx: StepContext) -> None: ...
def run_operons(...) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: ...
```

1. **Inputs.** `gene_table` is the detectability output (`genetable/genes.parquet`), as CONTEXT
   already takes it. The three layer tables are optional (missing → NA evidence + the existing
   `layer_missing:<layer>` run flag). `peptides`, `peptide_hits` and `peptide_runs` are optional and
   used only by M6-02 in pass 2; in pass 1 the step logs `no_detection_limit_model` and proceeds.
2. **Order.** pair frame (M6-04) → discordance (M6-05) → missed operons (M6-06) → flags and the TU
   table (M6-07) → three writers.
3. **Writers** in `src/giap/io/operon_tables.py`: 1-based closed coordinates, `NA` for missing,
   `;`-joined flags and list columns, booleans as `true`/`false`, the §6.5 sort order. Coordinates
   are converted **only here** (the §6.5 rule).
4. **Schemas.** Three Frictionless resources in `src/giap/schemas/`, with the `flag` enum of M6-07,
   the `class` and `test` enums of §9.3, and `coordinates: 1-based-closed`. The new Operons flag
   values of index decision 7 are added to the flag vocabulary resource.
5. **Module skip.** No abundance input → the step still runs and writes `pairs.tsv` (classes,
   distances, evidence) with every test `none`; `operon_flags.tsv` is header-only. Nothing in the
   pipeline fails for lack of quantitation.
6. **Nextflow.** `workflow/subworkflows/local/operons.nf` in the CONTEXT style (`GIAP_STEP`,
   `stepIn`, `kv`, `pick`), called from `main.nf` after `GENETABLE` and beside `CONTEXT`. It needs
   three channels `annotate.nf` does not emit today — `promoter_windows`, `terminators`,
   `tx_rates` — so `annotate.nf` gains `emit: layers_for_pairs` (decision 6). Each is optional, in
   the way CONTEXT already treats `regions`.
7. **publishDir.** `tables/*.tsv` is published; work Parquet is not. This settles the MVP-01 pass-2
   backlog item "which keys are user-facing, which stay in `work/` — with the M6 writers": the rule
   is written into the subworkflow and into `docs/specs/M6-08` so later milestones copy it.
8. **Small refactor named in scope:** `_genes`, `_path` and `_replicons` move out of
   `giap.steps.operons` into a public shared module (the standing M2-16 backlog item), because this
   step and M6-02 both import them.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no abundance | step 5 |
| no peptides and one run | one-detected tests `none` + `no_detection_limit_model` (pass 1 always) |
| layer missing | NA evidence column, run flag, no flags of the affected rows |
| sim100k (≈ 100 genes) | `too_few_op_pairs` and `too_few_no_pairs` are the normal outcome; the run **passes** |
| empty `operon_flags.tsv` | header-only file, valid against its schema |
| an upstream channel empty in Nextflow | the step does not run (the `stepIn`/`need` behaviour, including the known MVP-01 pass-2 bug on empty `need` — this issue passes the real required keys explicitly) |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `giap run --executor local` on the sim100k fixture with all layers | three files under `tables/`; `pairs.tsv` has one row per adjacent pair plus the SW rows; every test column `NA`/`none`; exit 0 |
| 2 | the same through Nextflow (`-profile test`) | byte-identical tables (the §6.4 executor-equality test) |
| 3 | run without `--abundance` | `pairs.tsv` written, `operon_flags.tsv` header-only, no error |
| 4 | run with `annotation.precomputed.terminators` unset and no TerminatorNet | `terminator_between` NA everywhere, run flag `layer_missing:terminators` |
| 5 | the 50-kb golden fixture | `pairs.tsv` matches the committed golden file byte for byte |
| 6 | `frictionless validate` on the three resources | valid |
| 7 | nf-test stub run of `OPERONS` | the three output keys are emitted and named as in §9.1 |
| 8 | `sbatch scripts/slurm/giap_run.sbatch` on sim100k, real tools | completes; the three tables exist and are non-empty (`operon_flags.tsv` may be header-only) |

## Acceptance criteria

- [ ] `tests/unit/steps/test_operons_step.py`: examples 1, 3, 4, 6; `tests/integration` adds
  example 2 beside the existing executor-equality test.
- [ ] `workflow/tests/subworkflows/operons.nf.test`: stub coverage in the style of #130, plus a
  registry test that the step's declared keys match the subworkflow's (as
  `test_nf_annotate_registry.py` does).
- [ ] Golden: `tests/golden` gains the 50-kb `pairs.tsv` (example 5). `operons.tsv` and
  `operon_flags.tsv` join it in pass 2, when the fixture has enough detected pairs to be
  informative.
- [ ] Schema: `frictionless validate` passes; the `flag`, `class` and `test` enums are exactly the
  §9.3 ones; `schemas/datapackage.json` lists three new resources.
- [ ] **Evidence that counts (the MVP rule):** one real `sbatch scripts/slurm/giap_run.sbatch` on
  sim100k that reaches the end and writes the three tables; the job ID in PROGRESS.md.
- [ ] Runtime recorded in `benchmarks/baseline.json` for a 4,500-gene genome.

## Hardening owed (pass 2)

- Contract assertions on the three tables **inside** the pipeline, not only in the Python tests
  (the MVP-01 pass-2 pattern).
- `-resume` correctness across `OPERONS`, and the behaviour when the step fails mid-graph.
- The real-tool arm under `-profile farmshare,apptainer` (pass 1 uses the local executor and
  recorded layers).
- Golden `operons.tsv` and `operon_flags.tsv` once a fixture with ≥ 20 usable OP pairs exists
  (M6-03 can produce one: a 1,000-gene simulated genome with high detection).
- Manifest: record the cut-off provenance (`operon_distance`, `terminatornet_yes`, their
  `calibrated` flags and sources) — **not** an `uncalibrated` flag (ADR-0021 Amendment 1).

## Questions for the lead

1. `annotate.nf` emits the three layer tables (proposed), or `context_layers` is taught to write a
   `pair_evidence.parquet`? **Proposed default:** emit — the alternative changes a merged M2 step
   and forces a fixture re-record.
2. Does Module 4 need a `--modules` flag of its own? §5.3 makes it core. **Proposed default:** core,
   with the graceful no-abundance behaviour of step 5.
3. Is a header-only `operon_flags.tsv` on sim100k an acceptable pass-1 definition of done?
   **Proposed default:** yes — the alternative is to enlarge the integration fixture, which belongs
   to pass 2 (M6-03 hardening item above).

## Lead decisions (2026-10-21)
1. Q1: `annotate.nf` emits `promoter_windows`, `terminators` and `tx_rates` (M6-index decision 3).
2. Q2: Module 4 is core. It has no `--modules` flag and keeps the graceful no-abundance behaviour.
3. Q3: a header-only `operon_flags.tsv` on sim100k is an acceptable pass-1 done.
4. **Inputs:** the M6-04 `pair_frame`, M6-05 `discordance`, M6-06 `missed_operons` (always pass `genes`) and M6-07 `merge_tests`/`flag_pairs`/`operons_table`, as on `m6/flags`.
   - `pairs.tsv` never carries the internal `which` column.
   - `merge_tests(..., missing_layers=...)` is given the annotation layers that are actually absent.
5. **Local/Nextflow parity:** add the three tables to the local executor's `PUBLISHED_FILES` (M4-21). `tests/unit/test_local_publish_names.py` pins them against the Nextflow `'key=path'` names.
6. **Manifest:** record π₀ (`operons.pi0`) and the decoy `bin_counts`, which the #148 review left owed to M6-08.
7. **The pair-scope flag vocabulary is extended through the `pairs.tsv` schema** (§9.4, "extend only through the schema"; #152 review). It adds `abundance_imputed`, `layer_missing:promoters`, `layer_missing:terminators`, `layer_missing:tx_rates` (added after the test step), `ci_unestimated`, `unestimated`, `too_few_op_pairs`, `too_few_no_pairs`, `no_detection_limit_model` and `relative_scale`.
   - `abundance_imputed` stays the pair-scope token. It is distinct from the gene-scope `imputed_abundance`, and the schema documents both.
   - `operon_flags.tsv` and `operons.tsv` take the same values where they carry flags.
8. After the test step (594fbd6): the TA's interface is accepted.
   - `pairs_raw` travels inside `layers_for_pairs` as its required member.
   - The process alias name is free; tasks are found by their tag.
   - Decision 6 is pinned as a §6.7 event `operons.missed_operons`.
   - The `main_genetable_context.nf.test` "precomputed regions" case now asserts that CONTEXT's tables are absent, not that `tables/` is (OPERONS publishes there).
   - Example 5 (a 50-kb golden `pairs.tsv`) moves to pass 2: no fixture has `pairs_raw` yet (M6-03).
   - If example 1's "every test none" fails on real sim100k because ≥ 20 usable OP pairs exist, the implementer reports it and the lead revises the example; never weaken the test silently.
