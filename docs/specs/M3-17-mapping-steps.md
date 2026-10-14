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
