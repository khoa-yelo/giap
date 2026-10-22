# M7-04: The `variants` step, the VARIANTS subworkflow, and a real `giap run`

Milestone M7 · track T (tables and integration) · **pass 1 — the Module 5 integration milestone**
· depends on: M7-01, M7-02, M7-03; merged MVP-01 wiring patterns (`GIAP_STEP`, `stepIn`, `kv`,
`pick`, `inDir` in `workflow/subworkflows/local/*.nf`) · estimate ≈ 300 code lines (+ ≈ 300 test
lines, + nf-test) · not statistical

## The gap this closes

Module 5's three tables must reach `giap run` before any rigour is added (Khoa, 2026-10-20/21).
This issue registers the step, wires it into both executors, resolves `sample_name: auto`, and
proves it with a real Slurm run. Flags (M7-05) and class error rates (M7-10) land later inside the
already-wired step and change no wiring.

## Spec sections implemented

- §6.2 step row: "`variants` | py | hits, genome, genes → `pileup`, `variants.tsv`, `variants.vcf` |
  core".
- §9.1 (`tables/pileup.tsv.gz`, `tables/variants.tsv`, `tables/variants.vcf`), §9.5 (validation),
  §6.3/§6.4 (both executors, identical tables), §6.5, §6.7 (step log and events).
- §5.2 `sample_name: auto  # auto = stem of the peptide (or protein) file` and §7.5.5 "Sample
  column: `config.sample_name` (default: the peptide file stem)" — today nothing resolves `auto`.
- §1.3 table: Module 5 is **core** (no `--modules` switch).

## Public interface

`src/giap/steps/variants.py`
```python
@step("variants",
      inputs=["mapping", "genome", "genes", "replicons", "peptides"],
      optional_inputs=["hits_tolerant", "reference_peptides"],
      outputs=["pileup", "variants", "variants_vcf"])
def variants_step(ctx: StepContext) -> None: ...
def run_variants(peptide_hits: pd.DataFrame, hits_tolerant: pd.DataFrame | None,
                 peptides: pd.DataFrame, gm: GeneModel, genome: Mapping[str, str],
                 cfg: GiapConfig) -> tuple[pd.DataFrame, pd.DataFrame, list[VcfRecord]]: ...
```
`src/giap/cli.py` / `giap.config`: `resolve_sample_name(cfg, peptides_path, proteins_path) -> str`.

## Rules

1. **Inputs.** `mapping` is the LOCI output directory (its `peptide_hits.tsv`, read with the M3-16
   reader; C11). `peptides` is C9 `peptides.norm.parquet` (for `source`, `search_db`). `genes` and
   `replicons` are C2/C1. `hits_tolerant` is optional: without it, frameshift rows are skipped with
   the run flag `layer_missing:hits_tolerant` and a warning (no silent loss). `reference_peptides`
   is declared now and used only by M7-10.
2. **Order.** read → M7-01 `align_hits` → M7-02 `eligible`/`pileup`/`sites` → M7-03
   `vcf_records` → three writers. One pass, no randomness (no RNG stream needed).
3. **Outputs.** Under Nextflow: `pileup=tables/pileup.tsv.gz`, `variants=tables/variants.tsv`,
   `variants_vcf=tables/variants.vcf`, published (`tables/*` is user-facing; work Parquet is not —
   the M6-08 publishDir rule). Under `--executor local`: the same three keys, and the local
   publisher copies them to `OUT/tables/` (add them to `giap.runner.local.PUBLISHED` if the
   table-publishing change the lead mentioned has not landed on main by then; on c4d8841 only
   `annotation/` and `mapping/` are published).
4. **`sample_name`.** `giap run` resolves `auto` before writing `resolved_config.yaml`: the peptide
   file's name without `.gz` and its last suffix (`run_A.tsv.gz` → `run_A`), else the protein
   file's, else `sample`; then mapped to a VCF-safe ID (`[A-Za-z0-9._-]`, other characters → `_`,
   as preflight maps contig IDs, ADR-0005). The manifest records the resolved value.
5. **Nextflow.** `workflow/subworkflows/local/variants.nf`, in the CONTEXT style: `take:` genome,
   `gene_model` (`[genes, features, replicons]`), `mapping` (LOCI directory), `hits_tolerant`,
   `peptides`, `reference_peptides`; one `GIAP_STEP as VARIANTS`; `needed` lists the required keys
   explicitly (avoids the known empty-`need` bug of `stepIn`, MVP-01 pass 2). `mapping.nf` gains
   `emit: hits_tolerant = pick(MT, 'hits_tolerant.parquet')`. `main.nf` calls `VARIANTS` right
   after `MAPPING`, inside the `if (inputs.peptides)` block, independent of GENETABLE.
6. **Stub.** The stub writes header-only `pileup.tsv.gz`, `variants.tsv` and a header-only VCF with
   the §7.5.5 header, so `nf-test` stub runs exercise publishing.
7. **Logging.** `logs/variants.log`: counts of aligned hits, covered codons, sites by call and
   confidence; `events.jsonl` event `variants.summary` with the same counts.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| a peptide table with no accepted placement | the step runs on an empty `peptide_hits.tsv`: header-only tables, header-only VCF, exit 0 |
| no `--peptides` (proteins only, or no evidence) | `main.nf` does not enter the peptide block; no Module 5 files (as MAPPING today); the manifest lists `variants` under steps not run, with reason `no_peptides` |
| `hits_tolerant` missing | rule 1 |
| `peptide_hits.tsv` disagrees with the genome (M7-01 raises) | the step fails; `giap run` exits 4 naming `variants` and its log |
| sample name with spaces | rule 4 maps them to `_` |
| two runs, same seed | byte-identical three files (TE7) |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `giap run --executor local` on sim100k (precomputed layers, de novo peptides with errors) | `tables/pileup.tsv.gz`, `tables/variants.tsv`, `tables/variants.vcf` exist; the VCF reads in pysam; exit 0 |
| 2 | the same through Nextflow (`-profile test`, stubs where tools are absent) | byte-identical three tables (the §6.4 executor-equality test, beside the existing one) |
| 3 | database peptides searched against the strain's own Bakta proteins (`search_db same`) | `variants.tsv` rows only at codons the peptides really contradict (on sim100k: none except annotated non-ATG starts) |
| 4 | `--peptides run_A.tsv.gz`, `sample_name auto` | VCF sample column `run_A`; `resolved_config.yaml` and `run_manifest.json` say `run_A` |
| 5 | `--set sample_name="my sample"` | `my_sample` |
| 6 | nf-test stub run of `VARIANTS` | three outputs emitted and named as §9.1 |
| 7 | `frictionless validate` on the two TSV resources | valid |
| 8 | `sbatch scripts/slurm/giap_run.sbatch` on sim100k | completes; the three tables exist (`variants.tsv` may be header-only) |

## Acceptance criteria

- [ ] `tests/unit/steps/test_variants_step.py`: examples 1, 3, 4, 5, 7 on the simulator, and each
  degenerate row.
- [ ] `tests/integration`: example 2 beside the existing executor-equality test; TE7 on the three
  files.
- [ ] `workflow/tests/subworkflows/variants.nf.test` (stub, the #130 style) and a registry test
  that the step's declared keys match the subworkflow's, as `test_nf_annotate_registry.py` does.
- [ ] `main.nf.test` (stub) shows `VARIANTS` running after `MAPPING` when peptides are given and
  not at all without them.
- [ ] **Evidence that counts (the MVP rule):** one real `sbatch scripts/slurm/giap_run.sbatch` on
  sim100k reaching the end with the three tables; the job ID in PROGRESS.md.
- [ ] Runtime of the step on a 4,500-gene genome with 50,000 peptides recorded in
  `benchmarks/baseline.json`.

## Hardening owed (pass 2, M7-15)

- Contract assertions on the three tables inside the pipeline; `-resume` across VARIANTS; the
  failure path when VARIANTS errors mid-graph.
- The real-tool arm under `-profile farmshare,apptainer` (pass 1 uses the local executor and the
  stub/real mix of `-profile test`).
- Golden `variants.tsv`/`variants.vcf` for the 50-kb fixture.

## Questions for the lead

1. `sample_name: auto` is resolved nowhere today. Resolve it in `giap run` (proposed) or in the
   variants step? **Proposed default:** `giap run`, so the manifest, the VCF and later the report
   all see one value; the step never sees `auto`.
2. Should VARIANTS run when only proteins (no peptides) are given? **Proposed default:** follow
   MAPPING's existing condition (`if (inputs.peptides)`); protein-only runs carry no residue-level
   evidence, so Module 5 has nothing to report and writes nothing.
3. `hits_tolerant` as a step input (M7-01 Q3) — confirm. **Proposed default:** yes, optional.

## Lead decisions

1. **Questions 1–3:** the proposed defaults are accepted.
   - `sample_name: auto` is resolved in `giap run` before any step runs. `write_vcf` rejects `auto`, as M7-03 decision 3 has it.
   - VARIANTS follows MAPPING's `if (inputs.peptides)` condition.
   - `hits_tolerant` is an optional step input.
2. **From M7-03 (#160):**
   - the pileup is passed in memory within the step, as a DataFrame with its `attrs`;
   - the step calls `eligible` before `pileup` (M7-02 review);
   - the run copies `schemas/datapackage.json` so that its `tables/...` paths resolve against the run root.
3. **Parity:** add the published names to `giap.runner.local.PUBLISHED_FILES`, so that they equal the Nextflow `'key=path'` names (`tests/unit/test_local_publish_names.py`).
4. **The test author's choices** (commit 1d51456) are accepted as the contract.
   - **`resolve_sample_name`** lives in `giap.config`. It is pure, and it applies the VCF-safe mapping to explicit names too.
   - **`Support`** is built from the C11 `source`/`placement_q`/`uniqueness` plus C9 `search_db`.
   - **The VCF header** carries `giap.__version__`, contigs in replicon order, and `##reference=` set to the SHA-256 of the genome input.
   - **The `variants.summary` event** has the keys `aligned_hits`, `covered_codons`, `sites_by_call` and `sites_by_confidence`.
   - **Without `hits_tolerant`:** frameshift rows are dropped and a WARNING says `layer_missing:hits_tolerant`.
   - **The stub** must write a real gzip pileup and a full header-only VCF. Extending the GIAP_STEP stub mechanism for this is in scope.
   - **The datapackage** goes to `OUT/schemas/datapackage.json`, with paths resolving under `OUT`. This settles the M7-03 review's N2.
   - The branch has origin/m7/tables (including the datapackage fix) and origin/main merged in (d713b6d).
