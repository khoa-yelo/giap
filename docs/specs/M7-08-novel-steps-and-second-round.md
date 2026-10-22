# M7-08: The novel-ORF steps, second-round Promotech and OSTIR, and a real `giap run`

Milestone M7 · track T (tables and integration) · **pass 1 — the Module 6 integration milestone**
· depends on: M7-06, M7-07; merged M2-13/M2-16 steps and processes (`PROMOTECH`, `OSTIR`), the
MVP-01 Nextflow helpers; M6-08's `annotate.nf` layer emission (or this issue adds the
`terminators` emission itself if M6-08 has not merged) · estimate ≈ 380 code lines (+ ≈ 90
simulator lines, + ≈ 350 test lines, + nf-test) · not statistical

## The gap this closes

The §13 deliverable "second-round Promotech and OSTIR" is an integration item: two external tools
run a second time on windows that exist only after mapping. This issue registers the five steps,
runs the second round under both executors, makes the simulator supply the second-round layers for
the local executor, writes the four Module 6 files through `giap run`, and proves it with a real
Slurm run — before the evidence rules (M7-09) and the FDR (M7-11) arrive.

## Spec sections implemented

- §6.2 rows: "`novel_candidates` | py | loci → candidate ORFs, extra Promotech and OSTIR windows";
  "`promotech_novel`, `ostir_novel` | ext | candidate windows → scores"; "`novel_orfs` | py |
  candidates + second-round scores → `novel_orfs.*`, `annotation_patch.gff3`".
- §7.6.2 "Second-round layers"; §7.1.4 (chunks of 100,000 windows run in parallel); §7.1.1 tool
  rows for Promotech RF-HOT and OSTIR 1.1.3 and their gotchas; §6.3 (one process per tool, labels,
  `versions.yml`, stubs, `publishDir`); §6.4 ("Every external layer must come from
  `annotation.precomputed` … a missing layer is an error naming the key"); §11.3 (the simulator
  writes all precomputed layers); §9.1.

## Public interface

`src/giap/steps/novel.py`
```python
@step("novel_candidates",
      inputs=["mapping", "peptides", "genes", "features", "replicons", "genome", "annotation_meta"],
      optional_inputs=["contaminants"],
      outputs=["novel_candidates", "novel_windows", "novel_windows_map", "ostir_input_novel"])
external("promotech_novel", inputs=["novel_windows"], outputs=["promotech_novel"],
         precomputed="promoters")
@step("novel_promoter_layer", inputs=["novel_windows_map", "promotech_novel"],
      outputs=["novel_promoter_windows"], precomputed="promoters")
external("ostir_novel", inputs=["ostir_input_novel"], outputs=["ostir_novel"],
         precomputed="rbs_novel")
@step("novel_rbs_layer", inputs=["ostir_novel", "ostir_input_novel"], outputs=["rbs_novel"],
      precomputed="rbs_novel")
@step("novel_orfs",
      inputs=["novel_candidates", "gene_table", "genes", "replicons", "genome"],
      optional_inputs=["novel_promoter_windows", "rbs_novel", "promoter_windows", "terminators"],
      outputs=["novel_orfs", "novel_orfs_gff3", "novel_orfs_faa", "annotation_patch"])
```
`src/giap/annotate/ostir.py` gains `parse_ostir_starts(path, requested: Mapping[str, set[int]])
-> dict[str, dict[int, RbsFeatures]]` (every requested start of an id, not only one).
`src/giap/config/model.py` + `defaults.yaml`: `annotation.precomputed.rbs_novel: null` (decision 1).
`src/giap/annotate/precomputed.py`: layer `rbs_novel` = `contig, strand, start (pos), tir_ostir,
sd_dg, mrna_dg, spacing_dg, standby_dg, spacer_bp, start_codon`.

## Rules

1. **`novel_candidates`** runs M7-06 and writes `novel_candidates.parquet` (one row per ORF ×
   peptide), the window chunks (`novel_windows/windows_NNN.fa`), `novel_windows_map.parquet` and
   `ostir_input_novel.csv`. It deduplicates its windows against nothing — first-round scores are
   reused by position in `novel_orfs` (rule 4), so a window scored twice costs Promotech time but
   no correctness (decision 2).
2. **`novel_promoter_layer`** = M2-13's `promoter_layer` on the novel map; **`novel_rbs_layer`**
   reads `ostir_novel.csv` with `parse_ostir_starts` and writes `rbs_novel.parquet` keyed by
   `(contig, strand, start)` (0-based start codon position) with the OSTIR features.
3. **Empty second round.** With no candidate start, `novel_candidates` writes an empty window
   directory and a header-only OSTIR CSV; Nextflow skips `PROMOTECH_NOVEL` and `OSTIR_NOVEL` (no
   chunk, no row) and `novel_orfs` runs without the two optional inputs. Under the local executor
   the precomputed layers are used as usual.
4. **`novel_orfs`** joins second-round scores by position: Promotech score at a start = the maximum
   over its windows in `novel_promoter_windows ∪ promoter_windows` (first round, when present);
   OSTIR = the `rbs_novel` row at that start. It fills `start_best_ostir` (highest `tir_ostir`;
   ties → the more upstream start), `ostir_tir`, `promoter_score`; then writes the four M7-07
   files. Evidence, ranking and RBS support are M7-09; class FDR is M7-11. A missing optional
   layer gives NA columns and the run flag `layer_missing:<key>`.
5. **Local executor.** `promotech_novel`/`novel_promoter_layer` are replaced by
   `annotation.precomputed.promoters` (position-keyed already, so candidate windows are looked up in
   the same layer); `ostir_novel`/`novel_rbs_layer` by the new `annotation.precomputed.rbs_novel`.
   A run with peptides and without `rbs_novel` is an error naming the key (§6.4), exactly as a
   missing `rbs` is today. As for every replaced step, `novel_orfs` then receives the layer's TSV
   (M2-02 format, read with `load_precomputed`) instead of the Parquet, and must accept both.
6. **Simulator** (`tests/sim/giapsim/layers.py`): writes `rbs_novel.tsv` — one row per start codon
   (`annotation.start_codons`) on both strands of every replicon, `tir_ostir` drawn from the
   "true" distribution (lognormal 8, as `rbs_row`) at planted novel-ORF and planted N-terminal
   true starts and from the background (lognormal 6) elsewhere — and adds the key to its
   `config.yaml`. Recorded configs under `tests/data/` that run the local executor with peptides
   gain the key in the same PR.
7. **Nextflow** `workflow/subworkflows/local/novel.nf` (CONTEXT style). `take:` genome,
   gene_model, contract (C4 `annotation_meta`), mapping, peptides, contaminants, gene_table,
   layers (first-round `promoter_windows`, `terminators`; each optional). Processes:
   `GIAP_STEP as NOVEL_CANDIDATES`, `PROMOTECH as PROMOTECH_NOVEL` (one task per chunk, the M2-13
   pattern), `GIAP_STEP as NOVEL_PROMOTER_LAYER`, `OSTIR as OSTIR_NOVEL`, `GIAP_STEP as
   NOVEL_RBS_LAYER`, `GIAP_STEP as NOVEL_ORFS`. `main.nf` calls `NOVEL` after `GENETABLE` (it needs
   `gene_table` for M7-09's RBS quantile). Publishes `tables/novel_orfs.tsv`,
   `tables/novel_orfs.gff3`, `tables/novel_orfs.faa`, `tables/annotation_patch.gff3`.
8. **Tools.** No new tool, env or image: the second round reuses `workflow/envs/promotech.yml`
   (Promotech @ `56251ad`, GPL-3.0, conda env with lock; RF-HOT model from `$GIAP_DB_DIR/promotech`)
   and `workflow/envs/ostir.yml` (OSTIR 1.1.3, GPL-3.0, ViennaRNA 2.6.4; lock present). Both stay
   separate processes, never imported (`scripts/ci/check_imports.py`). Neither has an upstream
   image (ADR-0008), so `-profile apptainer` runs them from the conda env.
9. **OSTIR range recording.** One real OSTIR run with a multi-start `start`/`end` window on three
   sim100k candidates (minutes, `sbatch` with the `giap-` prefix, ≤ 4 CPUs) is recorded under
   `tests/data/recorded/ostir/1.1.3/sim100k_novel/` with `PROVENANCE.md`; `parse_ostir_starts`
   gets a contract test on it (the existing recordings only hold single-start rows).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no peptides | the NOVEL block is not entered (as MAPPING) |
| no candidate | rule 3: header-only TSV, pragma-only GFF3s, empty FASTA; exit 0 |
| ORF without a start | no windows, no OSTIR row (M7-06); NA scores |
| OSTIR returns no row at a requested start | that start has NA `tir_ostir`; if every start is NA, `ostir_missing` (M7-07) |
| local run without `annotation.precomputed.rbs_novel` | exit 2 naming `annotation.precomputed.rbs_novel` (rule 5) |
| first-round `promoter_windows` absent (precomputed promoters) | the precomputed layer serves both rounds |
| `terminators` layer missing | ADR-0012 columns NA, run flag `layer_missing:terminators` (M7-09) |
| GENETABLE skipped (no gene table) | NOVEL does not run (gene_table is required); the manifest records it |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `giap run --executor local` on sim100k with planted novel ORFs (M3-18 counts) and the new `rbs_novel` layer | the four files under `tables/`; `novel_orfs.tsv` has one row per reconstructed ORF; `start_best_ostir` set where a start exists; exit 0 |
| 2 | the same through Nextflow (`-profile test`, stubs) | byte-identical tables (executor equality) |
| 3 | stub run of `NOVEL` with two chunks | `PROMOTECH_NOVEL` runs twice, `OSTIR_NOVEL` once; the four outputs are emitted and named as §9.1 |
| 4 | a run whose peptides give no candidate | header-only/pragma-only files; `PROMOTECH_NOVEL` and `OSTIR_NOVEL` do not run |
| 5 | `parse_ostir_starts` on the recorded range run | one feature row per requested start, `start_position` matched 1-based |
| 6 | local run with peptides, sim config without `rbs_novel` | exit 2, message names `annotation.precomputed.rbs_novel` |
| 7 | `sbatch scripts/slurm/giap_run.sbatch` on sim100k (real Promotech and OSTIR in the second round, conda envs) | completes; four files present; trace shows `PROMOTECH_NOVEL` and `OSTIR_NOVEL` tasks |

## Acceptance criteria

- [ ] `tests/unit/steps/test_novel_steps.py`: examples 1, 4, 6; rules 2–4 on hand-built inputs.
- [ ] `tests/contract/test_ostir_range.py`: example 5 on the recording (rule 9).
- [ ] `tests/integration`: example 2 beside the existing executor-equality test; TE7 on the four
  files.
- [ ] `workflow/tests/subworkflows/novel.nf.test` (stub): example 3; the registry test that the five
  steps' declared keys match the subworkflow's.
- [ ] Config sync test passes with the new `annotation.precomputed.rbs_novel` key (defaults.yaml ↔
  model); `precomputed.py` validates the layer.
- [ ] **Evidence that counts:** example 7 as one real `sbatch scripts/slurm/giap_run.sbatch` run;
  the job ID in PROGRESS.md. R14: the second round is ≤ 2 chunks on sim100k; the job stays within
  16 CPUs.
- [ ] Runtime of `novel_candidates` and of the second round on a 4,500-gene genome recorded in
  `benchmarks/baseline.json` (`scripts/measure_runtimes.py`).

## Hardening owed (pass 2, M7-15)

- The real-tool arm under `-profile farmshare,apptainer` on K-12 or Nissle with real peptides
  (pass 1 proves sim100k).
- Contract assertions on the four files inside the pipeline; `-resume` across the second round
  (chunk caching); a failing `OSTIR_NOVEL` task stops the run with exit 4 naming it.
- Second-round window deduplication against first-round windows (decision 2), if the runtime
  measurement shows it matters.

## Questions for the lead

1. The local executor needs a precomputed source for second-round OSTIR (§6.4 makes a missing layer
   an error). **Proposed default:** a new key `annotation.precomputed.rbs_novel` (null by default)
   holding a position-keyed OSTIR layer, written by the simulator for every start codon. This adds
   a config key (a `defaults.yaml` change, ADR-worthy; ADR-0005 item 4 shows Khoa rejected new keys
   once, for thresholds). Alternative without a new key: allow `START_<contig>_<pos1>_<p|m>` rows in
   the existing locus-keyed `rbs` layer — no config change, but it overloads `locus_tag`. The lead
   decides, and whether it needs Khoa.
2. Deduplicate second-round windows against first-round ones? **Proposed default:** not in pass 1
   (rule 1). The first-round scores are merged by position, so duplicates only cost runtime; the
   M7-08 benchmark decides whether pass 2 adds it.
3. `NOVEL` needs `gene_table` (RBS quantile at detected genes, M7-09), so it runs after GENETABLE;
   without a gene table it does not run. **Proposed default:** accept — the RBS threshold is
   undefined without detection calls. The alternative (run without and mark every RBS test
   `unestimated`) adds a branch for a case `giap run` never produces.
4. Is a sim100k run with few or no *called* ORFs an acceptable pass-1 definition of done?
   **Proposed default:** yes, as M6-08 Q3: the acceptance criterion is "the four files exist with
   the right schema"; M7-12 builds the scenario where calls are informative.

## Lead decisions (2026-10-21)

1. **Q1: a new key.** Add `annotation.precomputed.rbs_novel`, null by default: a position-keyed OSTIR layer that the simulator writes for every start codon.
   - Overloading `locus_tag` in the `rbs` layer is rejected because it breaks that column's meaning.
   - This is an input path beside the other `annotation.precomputed.*` keys, not a threshold, so it is not the case ADR-0005 item 4 rejected.
   - It is recorded as ADR-0024 (proposed, Needs Khoa: yes, because it changes the config interface). Pass 1 goes ahead on it.
2. **Q2:** no deduplication against first-round windows in pass 1. Accepted as proposed; the benchmark decides for pass 2.
3. **Q3:** `NOVEL` runs only when a gene table exists. Accepted as proposed.
4. **Q4:** accepted as proposed. A sim100k run with few or no called ORFs is done when the four files exist and validate against their schemas.
5. **Stacking:**
   - The branch m7/novel-steps is stacked on m7/novel-outputs (#165, itself on #161).
   - M6-08 (#155) is not merged yet, so this issue adds the `terminators` emission itself if `annotate.nf` lacks it. The lead resolves the overlap when main is merged.
6. **The test author's choices (commit 5a06f85) are accepted as the contract.**
   - **`rbs_novel` positions:** the first base of the codon in transcription direction, as in `novel_orfs.tsv`. `Orf.starts` and the `scored` `start_best_ostir` stay lowest-coordinate.
   - **`novel_rbs_layer` inputs:** it may take extra inputs from the `novel_candidates` outputs.
   - **Example 6:** exits 1 (`INVALID_INPUT`), as every missing layer does today. The spec's "exit 2" is an error in the spec.
   - **Example 2:** runs real Nextflow on `farmshare_local` with every layer precomputed. This revises ADR-0024: Nextflow honours a given `annotation.precomputed.rbs_novel`, and second-round `promoters`, the same way as the other precomputed keys, and runs OSTIR only when the key is null.
   - **The rest:** the promoter-score join on `(contig, strand, end3p)`, the WARNING `layer_missing:<key>`, and the pinned NOVEL Nextflow interface (`contract`, a `layers` channel of `[key, path]`, `params.stub.novel_windows_chunks`, emits named after the output keys).
   - **MVP-02 nf-test:** the lead updated it to check only CONTEXT's tables, because NOVEL also publishes under `tables/`.
   - **Rule 9 recording:** the implementer records `tests/data/recorded/ostir/1.1.3/sim100k_novel/` from a real OSTIR 1.1.3 run, using the official environment on Slurm.
