# M4-15: `genes.tsv` assembly, schema resource and the `gene_table` / `detectability` steps

Milestone M4 · track T (assembly) · depends on: M4-06, M4-10, M4-11, M4-12, M4-13, M4-14 · estimate
≈ 340 code lines (+ ≈ 300 test lines) · **can start now** on fixtures once its dependencies merge

## Spec sections implemented

- §7.3: "All statistics run on one gene table (one row per annotated CDS, sORF or pseudogene)".
- §9.3 `tables/genes.tsv` column groups (identity, context, regulatory and sequence features, evidence,
  abundance, function, flags), in that order.
- §9.2: "TSV, UTF-8, one header row, `NA` for missing, one entity per row, sorted (§6.5)";
  "Coordinates 1-based closed in TSV"; booleans `true`/`false`; flags `;`-separated controlled
  vocabulary. ADR-0006 encoding.
- §9.4 gene flags (all scopes "Gene"); "extend only through the schema".
- §9.5: "`schemas/datapackage.json` is a Frictionless Data Package with one resource per table … The
  descriptors live in `src/giap/schemas/`".
- §6.2 steps `gene_table` ("genes, layers, loci, abundance → `genes_pre.parquet`") and
  `detectability` ("gene table → `p_detect` (cross-fitted)"); §6.5 sort order and partial genes;
  §6.6 seeds (`giap.rng.for_step`); §1.2 rule 6 "Say 'detected', never 'expressed'".

## Public interface

`src/giap/genetable/assemble.py` (pure)
```python
GENES_COLUMNS: tuple[str, ...]        # exact §9.3 order (identity … flags)
GENE_FLAGS: tuple[str, ...]           # §9.4 Gene scope, vocabulary order

def function_columns(xrefs: pd.DataFrame, cog_def: pd.DataFrame | None) -> pd.DataFrame: ...
    # locus_tag, cog_id, cog_category, go_terms, kegg_ko, ec (pivot of C3; cog_category via COG 2020 table)
def assemble_genes(inputs: M4Inputs, cov, evidence, detection, detectability, abundance) -> pd.DataFrame: ...
def gene_flags(row_parts: Mapping[str, Sequence[str]]) -> list[str]: ...   # union, vocabulary order
```
`src/giap/io/genes_tsv.py`: `write_genes_tsv(path, df)` (0-based → 1-based: `start + 1`, `end`
unchanged; lists `;`-joined; via M1-04 `write_tsv`) and `read_genes_tsv(path)` (inverse, for
downstream modules and tests).
`src/giap/schemas/tables/genes.json`: Frictionless resource (name, type, unit, description,
enum constraints, `coordinates: "1-based-closed"`), `flags` items constrained to `GENE_FLAGS`.
`src/giap/steps/gene_table.py`: `run_gene_table(work_dir, cfg) -> Path` writes
`work/genetable/genes_pre.parquet` (everything but `p_detect`); `src/giap/steps/detectability.py`:
`run_detectability(work_dir, cfg) -> Path` adds `p_detect` → `work/genetable/genes.parquet`, writes
diagnostics rows to `work/genetable/detection_summary_parts/detectability.parquet`. Both stop
(exit code 1 path via the step runner) on any contract `error` issue.

Column sources (M4-06 table IDs): identity ← C2 (`aa`/`nt` not written), C14 (`source_id`), M4-10
(`length_aa`), C1 (`replicon_type` = `kind`); context ← C6 `layers.parquet` (`tu_id`, `tu_position`,
`tu_size`, `uniop_p_prev/next`, `region_ids`; M2-20 computes them) with `island_id`, `replichore`,
`replichore_position` = NA (M5) and `gc_local_dev` = NA (no producer yet); regulatory/codon ← C6
(all-NA when a layer is missing; M2 records `layer_missing:<layer>`); protein ← M4-10; evidence ←
M4-11/12 + `p_detect` (M4-13); abundance ← M4-14; function ← C3 via `function_columns`,
`annotation_source` ← C4. Flags: C6 `flags` (`partial`, `crosses_origin`, `pseudogene`,
`cai_reference`) ∪ evidence flags (M4-11) ∪ `two_proteins_one_gene`, `cross_strain_protein`,
`imputed_abundance`, `top3_lt3` (M4-14) ∪ `detectability_constant` (M4-13).

## Examples

| # | Input | Expected in `genes.tsv` |
|---|---|---|
| 1 | M4-07/08/09 default fixture (1,000 genes) | 1,000 rows; header equals `GENES_COLUMNS`; validates against `genes.json` with `frictionless validate` |
| 2 | internal gene `start=0, end=300` | TSV `start = 1`, `end = 300` |
| 3 | origin-crossing gene on a 50,000-bp circular contig, `start=49_900, end=50_200` | TSV `start = 49901`, `end = 50200`, `crosses_origin = true`, flag `crosses_origin` |
| 4 | C6 row for G2: `tu_id TU_7`, `tu_position 2`, `uniop_p_prev 0.9` | copied unchanged into `genes.tsv` |
| 5 | C6 `region_ids = [amr_1, prophage_1]` | `region_ids = amr_1;prophage_1` |
| 6 | C3 rows `(G1, COG, COG0593)`, `(G1, GO, GO:0005524)`, `(G1, EC, 3.6.4.12)`, `(G1, KEGG, K02313)`; COG table maps COG0593 → `L` | `cog_id = COG0593`, `cog_category = L`, `go_terms = GO:0005524`, `ec = 3.6.4.12`, `kegg_ko = K02313`; without the COG table `cog_category = NA` |
| 7 | detected gene with imputed abundance, a tolerant hit and `multi_identical` peptides | `flags = multi_identical;imputed_abundance;tolerant_evidence;strain_difference` (vocabulary order) |
| 8 | partial gene | `detection_status = partial_excluded`, `detected = false`, `p_detect = NA`, flag `partial` |
| 9 | `island_id`, `replichore`, `replichore_position` | `NA` (filled by M5) |
| 10 | run twice with the same seed | byte-identical `genes.tsv` (TE7) |

## Edge cases

- Rows sorted by contig (C1 row order = FASTA order), start, end, strand, locus_tag regardless of input order.
- No column contains the word "expressed" (test greps headers and enum values).
- A gene with no operon neighbour → `uniop_p_prev/next` NA.
- An unknown flag reaching `gene_flags` → `ValueError` (vocabulary is closed).

## Acceptance criteria

- [ ] `tests/unit/genetable/test_assemble.py`: examples 2–9; `tests/integration/test_gene_table_step.py`:
  examples 1 and 10 (fixture written to `tmp_path`, steps run, TSV validated with frictionless).
- [ ] Round trip: `read_genes_tsv(write_genes_tsv(df)) == df` (dtype-aware) on the fixture.
- [ ] Schema resource lists every column of `GENES_COLUMNS` with type and description; a test fails on
  any mismatch.
- [ ] Steps do I/O only through `io/`; `assemble.py` is pure.

## Methods note

None new; the statistician checks that `p_detect` is attached only from cross-fitted predictions
(cross-reference `detectability.md`).

## Questions for the lead

1. Step registration: M2-02 specifies the registry (`@step(name, inputs=…, outputs=…)`).
   **Proposed default:** register `gene_table` and `detectability` with it if M2-02 has merged;
   otherwise ship plain `run_*` functions and register them in M4-20.
2. Schema layout: one resource file per table under `src/giap/schemas/tables/` (as M3-16 does),
   assembled into `datapackage.json` by `finalize` later? **Proposed default:** yes.
3. `cog_category` needs the COG 2020 table from `giap setup --cog` (M2-09). **Proposed default:**
   read it when installed, else `NA` with run flag `layer_missing:cog`.

## Lead decisions (2026-10-20)
1. M2-02's `@step` registry has merged: register `gene_table` and `detectability` with it now (not in
   M4-20), and regenerate `workflow/tests/data/annotate/step_registry.json` if the registry snapshot
   test requires it.
2. Yes: one resource file per table under `src/giap/schemas/tables/` (as M3-16); `finalize` assembles
   `datapackage.json` later.
3. Yes: read the COG 2020 table when `giap setup --cog` installed it (manifest entry `cog2020_def`/
   `cog2020_fun`); otherwise `cog_category = NA` and run flag `layer_missing:cog`. Tests use a tiny
   fixture COG table, never the real `$GIAP_DB_DIR`.

## Lead decisions after the test-author's report (2026-10-20)
4. `gc_local_dev`: copied from C6 / M2-20 layers when the column is present, NA only when absent
   (M2-20 has merged; M4-index decision 2).
5. COG: `assemble_genes(..., cog_def: DataFrame | None = None)` (keyword, default None → `cog_category`
   NA); the step passes the table read from the optional input key `cog2020_def` (raw NCBI
   `cog-20.def.tab`, like `contaminants` reaches `loci`). The workflow resolves that path from the setup
   manifest; this refines decision 3.
6. A contract error raises `StepError` (runner exit code 4), not exit 1.
7. giapsim writes `annotation_source = "bakta_full"` so simulated outputs validate against the §9.3 enum
   (giapsim is tooling: the implementer changes it); the integration test's rewrite can then go.
8. Work-dir layout: `run_gene_table` reads the real `giap run` layout (TSVs under `mapping/`, C14 under
   `annotation/`) through the existing M4 input loader; the flat `default_paths` layout stays accepted
   for unit fixtures.
9. Several COG IDs on one gene: `cog_id` = the lexicographically smallest, `cog_category` from that ID;
   bare category xrefs (`COG:L`) are ignored (category comes from the table). No id map → `source_id` =
   `locus_tag`. Types: `island_id` string, `replichore` string (`right`/`left`, M5-08),
   `replichore_position` float; all NA in M4.
10. "I/O only through io/" is enforced for `genetable/assemble.py` only; the step-level convention is a
    backlog item (merged steps already read parquet directly).

## Lead decisions after the code review of PR #111 (changes requested)
11. (blocking) `read_genes_parquet` keeps missing values missing: float NaN and `pd.NA` in string columns
    stay null (never the text `"nan"`). Tests: a parquet round trip, and
    `pq.read_table(genetable/genes.parquet).column(c).null_count` equals the number of missing values for
    string columns (`island_id`, `replichore`, `cog_id`, `cog_category`).
12. Tests for COG resolution through the setup manifest: a tmp `$GIAP_DB_DIR` with `db_manifest.json` and a
    fixture `cog/cog-20.def.tab` → categories set; file present but no manifest entry → NA + warning.
13. `read_cog_def` falls back to latin-1 (the NCBI file has cp1252 bytes, e.g. 0x92); test with such a
    fixture byte.
14. Tests: a gene missing from a partial id map gets `source_id = locus_tag`; NA `detected_2pep` /
    `abundance_imputed` become False (kill `fillna(True)`).
15. `--executor local` stages `cog2020_def` (like `contaminants`) and passes `abundance` / `preflight`
    when the run produced them, so `cog_category` is set under the local executor too.
16. Rename the schema to `src/giap/schemas/tables/genes.schema.json` (the `<table>.schema.json`
    convention). Per-field `unit` (§9.5) is a project-wide backlog item.
17. The COG path comes from `db_sources.yaml` (`cog2020_def.path`), not a hard-coded constant.
18. Downstream: the stricter loader (parquet/TSV by table name) means M4-18/M4-20 supply C7/C8 as parquet.
19. The input key `preflight` stays; `stage_roots` copies a root file into `workdir/<key>` when that path is
    already a directory, and `gene_table` accepts `preflight` as `preflight.json` or a directory holding
    it (implementer, 9838394). Backlog: give staged root inputs their own namespace so keys never collide
    with executor directories.

## Lead decisions after the re-review of PR #111
20. Tests (blocking): (a) `m4_paths` turns a `preflight` directory into its `preflight.json` and the
    abundance scale is read from it; (b) a local run with an abundance input checks that
    `work/preflight/preflight.json` exists, `abundance` was staged, no `contract.preflight.abundance_scale`
    issue is logged and the abundance columns are filled; (c) a `stage_roots` test for the directory case.
21. `stage_roots`' copy-into-directory branch applies to the `preflight` key only (other keys still fail
    loudly on a directory) until the backlog namespace item lands.
22. `cog_table` moves out of `giap.runner.local` into `giap.setupdb` (or io), using a public accessor
    instead of `setupdb.items._load_sources`; steps do not import the executor.
