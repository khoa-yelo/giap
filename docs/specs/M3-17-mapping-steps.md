# M3-17: Mapping steps, query selection and workflow wiring

Milestone M3 · depends on: M3-01 … M3-16; M1-12 (normalised work files); M2-02 (step registry,
`StepContext`, `giap _step`) and M2-06 (Nextflow skeleton) (if absent, see Q1) · creates `src/giap/steps/mapping.py`,
`workflow/subworkflows/local/mapping.nf`, process modules `workflow/modules/local/giap_mapping/`
· pure Python + Nextflow (external tools only through M3-12/13/15 processes) · estimate ≈ 330
code lines

## Spec sections implemented

- §6.2 steps and their files: `sixframe_index` (genome → six-frame and decoy translations
  `.npz`), `map_exact` (peptides, index → `hits_exact.parquet`), `miniprot` (genome,
  `proteins.long.faa` → GFF), `tblastn` (genome, leftover proteins → tabular hits),
  `map_tolerant` (leftover peptides ≥ 10 aa, index → `hits_tolerant.parquet`, target + decoy),
  `placement` (all hits, decoy hits, input q/PEP → `placement.parquet`, `mapping_fdr.tsv`),
  `loci` (placements, genes → `peptide_hits.tsv`, `protein_hits.tsv`, `mappings.*`,
  `unmapped.tsv`).
- §6.1 data flow ("six-frame index ─► exact peptides ─► tolerant (+ decoys) ─► placement FDR";
  "proteins > 40 aa ─► miniprot ─► tblastn (leftovers)"); §6.3 one process per Python step group,
  `giap _step <name>`, stub blocks; §6.4 local executor (no external tools: protein routes need
  `annotation.precomputed`-style injection or are skipped when there are no long proteins).
- §7.2 "Each sequence takes the cheapest route that can place it; only failures move on."
- §7.3.1 (consumer contract): accepted evidence requires `q_value ≤ input.q_value_max`.

## Behaviour

- **Query peptides** = rows of `peptides.norm.parquet` (aggregated) with `is_decoy = false`,
  `q_value ≤ input.q_value_max` or `q_value` missing, length in `[exact.min_len, exact.max_len]`
  (Q2). Input decoy rows are never mapped (they serve gene-level FDR in M4).
- **Source per peptide** → decoy type: `database` → `mapping.decoys.database_input` (reverse),
  `denovo` → `denovo_input` (shuffle); `auto` is resolved as preflight recorded it (Q3).
- **Routes:** exact for all (target + decoy) → tolerant target search for exact-target leftovers
  ≥ `tolerant.min_len`, tolerant decoy search for every query ≥ `min_len` (M3-10 Q1) →
  frameshift (M3-08) → best costs → de novo precision (if any de novo peptides; reference
  peptides from `--reference-peptides`) → placement → loci → contaminants → outputs.
  Proteins: miniprot → assignment → tblastn on leftovers → assignment → blastp of the final
  leftovers against the contaminant library.
- **Gene model** comes from M3-01 (`gene_model_from_bakta_json`), or from M2-03/M2-04
  (`parse_bakta_json(...).model`, `genes.parquet`) once merged.
- Every step logs counts (queries, hits per route, accepted per stratum) to `logs/<step>.log`.

```python
@step(name="sixframe_index", inputs=[...], outputs=["work/sixframe.npz"], module="core")
def sixframe_index_step(ctx: StepContext) -> None: ...
# map_exact, map_tolerant, placement, loci: same pattern; pure functions from M3-02…16 do the work
def select_queries(peptides: pd.DataFrame, cfg: GiapConfig) -> pd.DataFrame: ...
```

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | peptides: `PEPTIDEK` (q 0.001), `rev_…` decoy row, `AAAK` (4 aa), `LONGPEPTIDE…` (q 0.2) | queries = `PEPTIDEK` only |
| 2 | 50-kb sim run, peptides only, `--executor local` | all seven Module 2 outputs exist; `protein_hits.tsv` header only |
| 3 | same with `proteins.faa` and miniprot unavailable on the local executor | error naming the missing tool (§6.4), exit 4 |
| 4 | nf-test, stub mode, `mapping.nf` | every process emits its declared outputs and `versions.yml` |
| 5 | same inputs, seeds equal, local vs Nextflow executor | byte-identical `mapping/*.tsv` (§6.4) |
| 6 | de novo input without `--reference-peptides` | run flag `c_assumed` recorded; de novo rows in `mapping_fdr.tsv` |

## Acceptance criteria

- [ ] `tests/unit/test_mapping_steps.py`: examples 1, 6 (pure selection and routing logic).
- [ ] `tests/integration/test_mapping_local.py`: examples 2, 3, 5 (5 needs `needs_nextflow`).
- [ ] `workflow/tests/mapping.nf.test`: example 4.
- [ ] Golden (`tests/golden/test_mapping_golden.py`): the four mapping TSVs for the 50-kb fixture
  (`GOLDEN_SEED`) equal the committed files in `tests/data/golden/sim50k/mapping/`.
- [ ] Benchmark: the whole mapping chain on the 100-kb sim genome with 5,000 peptides ≤ 20 s
  (local executor, excluding external tools); recorded in `benchmarks/baseline.json`.

## Questions for the lead

1. If M2-02's step registry (`@step`, `StepContext`, `giap _step`) is not merged when M3 reaches
   this issue: **Proposed default:** M3-17 waits for it (it is the only M3 issue that needs it);
   all other M3 issues are plain functions.
2. Peptides with q above `input.q_value_max`: map them or not? **Proposed default:** not mapped
   (never evidence, §7.3.1); they remain in the work table for M4's gene-level FDR, which does
   not need their genomic placement.
3. `input.peptide_source = auto` with no per-row `source`: **Proposed default:** `database` when
   the table has `q_value` or `pep`, else `denovo` — if M1-12 has not already fixed this.

## Lead decisions (after the test-author's pass)
1. The test-author's resolutions are accepted (root keys `genome`/`peptides`/`bakta_json`; intermediate keys
   named after files; Q3 per table; q at the cut-off is mapped; `c_assumed` flag + WARNING; shuffle decoys
   on the `sixframe_index` stream; MAPPING subworkflow with 5 inputs).
2. The 20 s timing check goes to benchmarks/ (implementer adds `benchmarks/test_mapping_steps_bench.py`).
3. Test load: the full suite and calibration run on Slurm; if the fast suite on GitHub CI exceeds its budget
   because of the chain runs, the implementer reports it and the test-author moves shared chain runs into a
   module-scoped fixture (follow-up).
4. If the work exceeds 800 gated lines, split: (a) steps + local-executor staging, (b) MAPPING subworkflow —
   tell the lead before splitting.

## Lead decisions after code review (PR #77, changes requested)
5. Split (size): #77 keeps the steps and the local executor; the MAPPING Nextflow subworkflow, the
   `giap_mapping` module and workflow/tests/mapping.nf.test move to a follow-up PR (branch
   m3/mapping-nf, after #77), which follows M2-07 (`GIAP_STEP` + `stepIn`, no duplicate processes,
   real stub headers) and wires `contaminants` (review findings 4 (nf part) and 7).
6. `search_db: auto` is resolved in the mapping steps (M3-10 decision 4: `same` when ≥ 95 % of peptides'
   protein_ids are in the database, else `other`), recorded in the run metadata; test.
7. Pruning removes only optional branches; a missing required root input (no evidence, proteins-only
   when peptides are required, a GFF3 annotation where Bakta JSON is needed) exits with the documented
   code and message, never 0 with missing outputs. Tests for the three cases.
8. Optional inputs come only from this run's staging, never from stale `work/<key>` files (clear or
   namespace per run); test with a leftover file.
9. Local `contaminants` staging and wiring (peptide library matching runs when the library exists).
10. `loci.py` `_engine_ids` accepts list/array `protein_ids` (contract C9), like `contaminants._as_ids`;
    the `_ids` join workaround in steps/mapping.py is removed. Test.
11. The local executor publishes mapping outputs to `OUT/mapping` (§9.1), like Nextflow.
12. Tests kill the surviving mutants (tolerant decoys per query, frameshift wiring, best-cost filter,
    de novo π̂, tolerant `min_len` boundary, reference peptides, leftover assignment, tblastn leftover
    filter, protein_ids separator, genome sha); `test_giap_step_cli_exits_4_on_a_corrupt_index` uses
    `fresh_steps`; tests never read the real $GIAP_DB_DIR (tmp DB dir with a test manifest).

## Lead decisions after the re-review (PR #77)
13. When an enabled target step is pruned only because a REQUIRED root is missing (e.g. `bakta_json`:
    `--peptides` without `--annotation` or `annotation.precomputed.bakta_json`), the local run exits 1
    naming the missing key; optional roots (`reference_peptides`, `contaminants`, `proteins_long`) are
    still pruned silently. Test for the review's probe E.
14. The `search_db: auto` resolution is recorded in `logs/events.jsonl` (an event with key and resolved
    value), not only a log line. Test.
15. Tests kill N10 (`loci` uses the resolved queries, so `auto` reaches M3-14's `same_db` mismatch
    check) and M10 (`placed` set).

## Lead decisions after the code review of PR #96 (changes requested)
15. Decision 14 under Nextflow: the per-task `logs/events.jsonl` written by `giap _step` (placement,
    loci, any step) is an optional GIAP_STEP output; the workflow collects them and `giap run` appends
    them to `OUT/logs/events.jsonl` after Nextflow ends (order: task completion as listed in the trace;
    each event keeps its step name). Test on the run directory, not the task directory.
16. Integration example 5 is `xfail(strict=False, reason="needs M2-20 context_layers")` until M2-20
    merges; its backlog line is restored.
17. Tests kill the surviving mutants: `gene_model` channel order (features vs replicons); main.nf wiring
    of `proteins_long` and `reference_peptides`; `params.json` contents written by `giap run`
    (`_NORM_INPUTS` peptides.norm.parquet, proteins_long, contaminants path).
18. The Groovy contaminants fallback reads the current `giap setup` manifest format (db_manifest.json as
    written today), tested with that shape.
