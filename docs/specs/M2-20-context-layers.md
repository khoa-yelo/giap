# M2-20: `context_layers` step: `layers.parquet`, `annotation/genes.gff3` and every `annotation/*` file

Milestone M2 · depends on: M2-02, M2-04, M2-11 … M2-19 · estimate ≈ 380 code lines (tests not
counted) · capabilities: none (runs on the simulator's precomputed layers with the local executor)

## Spec sections implemented

- §6.2 `context_layers` | py | "parsed layers → `annotation/*` (§9.1), `layers.parquet`".
- §7.1.8 Outputs of Module 1: "`annotation/genes.gff3`: the merged model with source-tagged layers.
  `annotation/regions.tsv` and `regions.bed`. `annotation/promoters.bed`, `annotation/terminators.bed`.
  `annotation/operon_pairs.tsv`, `annotation/tus.tsv`. `annotation/id_map.tsv`, when `--annotation`
  is a GFF3 or GBFF. `annotation/terminators_crosscheck.tsv`. `annotation/tools/<tool>/`: native outputs, untouched."
- §9.1 tree; §9.2 conventions; §9.3 column dictionaries for `regions.tsv`, `operon_pairs.tsv`,
  `tus.tsv` and the `genes.tsv` context / regulatory / codon columns this layer supplies.
- §6.3 / §9.4 "Downstream steps accept a missing optional layer and record `layer_missing`";
  §6.4 precomputed injection; §6.6 / R7 "byte-identical files under … `annotation/*.tsv`"; no
  timestamps or absolute paths.

## Public interface

`src/giap/annotate/layers.py`; step `context_layers` (inputs: `genes`, `replicons`, `pairs_raw`,
`tus`, optional `promoter_windows`, `tx_rates`, `terminators`, `terminators_crosscheck`, `rbs`,
`codon`, `regions`, `region_genes`, `id_map`, `tool_dirs`; outputs: `annotation_dir`, `layers`).
```python
LAYER_COLUMNS: tuple[str, ...] = (
    "locus_tag", "tu_id", "tu_position", "tu_size", "uniop_p_prev", "uniop_p_next", "region_ids",
    "promoter_score", "promoter_score_tu", "n_promoter_windows_tu",
    "tx_rate_max_tu", "tx_rate_sum_tu", "tss_distance_bp", "terminator_conf_downstream",
    "tir_ostir", "sd_dg", "mrna_dg", "spacing_dg", "standby_dg", "spacer_bp", "start_codon",
    "cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10", "rare_codons_2_10", "flags")

def gene_layers(genes, tus, pairs, promoters, tx, terms, rbs, codon, region_genes, replicons) -> pd.DataFrame
def merged_gff3_features(genes, tus, promoters_bed, terms, regions) -> list[GffFeature]
def copy_native(tool_dirs: Mapping[str, Path], dest: Path) -> None     # annotation/tools/<tool>/, byte copies
def layer_flags(cfg: GiapConfig, present: set[str]) -> tuple[str, ...] # "layer_missing:<layer>"
def write_annotation(outdir: Path, …) -> None                           # writes every §7.1.8 file
```
Column definitions (names match §9.3 `genes.tsv`):
- `tu_position` 1-based in transcription order; `tu_size` genes in the TU; `uniop_p_prev` / `_next`
  = p of the operon pair (M2-12) with the previous / next gene **in transcription order**, whether or
  not that pair is linked into the same TU; `NA` when that neighbour is on the other strand or absent
  (linear contig end).
- `region_ids` list (`;`), sorted by region id.
- `flags`: gene flags from this layer (`partial`, `crosses_origin`, `pseudogene`, `cai_reference`).

`annotation/genes.gff3` (M1-05 writer), one source per layer in column 2:

| Layer | Source | Type | Key attributes |
|---|---|---|---|
| genes | `Bakta` | `CDS` (cds, sorf), `pseudogene` | `ID=<locus_tag>`, `locus_tag`, `gene`, `product`, `feature_type` |
| TUs | `UniOP` or `distance_model` | `operon` | `ID=TU_<n>`, `genes` (list), `upstream=<start1>-<end1>` |
| promoters | `Promotech` | `promoter` | `ID=promoter_<n>`, `score` in column 6 |
| terminators | `TerminatorNet` / `TransTermHP` / `BacTermFinder` | `terminator` | `ID=<source>_<n>`, `confidence` |
| regions | `geNOMAD` / `antiSMASH` / `Bakta` / `AMRFinderPlus` | `prophage`, `plasmid`, `biosynthetic_gene_cluster`, `CRISPR`, — (amr not written) | `ID=<region_id>`, `genes` |

Files in `annotation/`: `genes.gff3`, `regions.tsv`, `regions.bed`, `promoters.bed`,
`terminators.bed`, `terminators_crosscheck.tsv` (when TransTermHP ran), `operon_pairs.tsv`,
`tus.tsv`, `id_map.tsv` (when present), `tools/<tool>/`. Also the Frictionless schemas for these
tables in `src/giap/schemas/annotation/*.json` (incl. `id_map`, `terminators_crosscheck`) (§9.5).

## Examples (simulator SIM100K with precomputed layers, local executor)

| # | Input | Expected |
|---|---|---|
| 1 | `giap run --executor local` on SIM100K's `config.yaml` up to `context_layers` | all §7.1.8 files exist; `annotation/tools/` contains only tools that ran (none with full precomputed) |
| 2 | TU `TU_4` = (A, B, C) on +, p(A,B)=0.9, p(B,C)=0.7 | B: `tu_position=2`, `tu_size=3`, `uniop_p_prev=0.9`, `uniop_p_next=0.7` |
| 3 | − TU (C, B, A) in transcription order (genome order A,B,C); the gene before A in genome order is on + | A: `tu_position=3`, `uniop_p_prev = p(A,B)`, `uniop_p_next = NA` |
| 4 | gene in `prophage_1` and `plasmid_1` | `region_ids = "plasmid_1;prophage_1"` |
| 5 | `modules=[core]` | `region_ids` without prophage/bgc ids; run flags contain `layer_missing:genomad` and `layer_missing:antismash` |
| 6 | `promoter_calculator.enabled=false` | `tx_rate_*_tu = NA`; flag `layer_missing:tx_rates` |
| 7 | `genes.gff3` read back with `read_gff3` | gene features equal the gene table (round trip); `##sequence-region` + `Is_circular=true` for circular replicons |
| 8 | the step run twice | every file under `annotation/` (excluding `tools/`) byte-identical (TE7) |
| 9 | each written TSV validated with its Frictionless schema | valid |

## Edge cases

- Genes without a TU cannot exist (singletons are TUs); a gene missing from `rbs`/`codon` gets `NA`.
- Partial genes keep layer values (statistics exclude them later, §6.5).
- `id_map.tsv` only with a GFF3/GBFF `--annotation`.
- Native copies are byte copies (no gzip, no rename); large native outputs (Promotech chunks) are
  copied as produced (§9.1 says untouched).

## Acceptance criteria

- [ ] `tests/unit/test_layers.py`: examples 2–6.
- [ ] `tests/integration/test_annotation_local.py`: examples 1, 7, 8, 9 on SIM100K (seconds, no tools).
- [ ] Truth check on SIM100K: ≥ 95 % of truth TUs reproduced exactly from the precomputed pairs; every
  planted region present in `regions.tsv`.
- [ ] Schemas for `regions`, `operon_pairs`, `tus`, `id_map`, `terminators_crosscheck` added.

## Questions for the lead

1. SO types for the merged GFF3 are not specified. **Proposed default:** the table above (`operon`
   SO:0000178, `promoter` SO:0000167, `terminator` SO:0000141, `prophage` SO:0001006,
   `biosynthetic_gene_cluster` SO:0001411, `CRISPR` SO:0001459).
2. `uniop_p_prev/next` orientation. **Proposed default:** transcription order (as above); a TU's first
   gene has `uniop_p_prev` = the p to its upstream same-strand neighbour (< threshold) or `NA`.

## Lead decisions after the test-author's 461fa2a (2026-10-16)
1. Step outputs are `annotation` and `layers` (the keys annotate.nf's CONTEXT_LAYERS uses; the spec's
   `annotation_dir` is renamed); optional inputs also take `features` and `promoters`.
2. Precomputed layers (§6.4) arrive as TSVs read through the config; `--modules core` drops prophage/bgc ids
   from precomputed regions; missing optional layers log a `layer_missing` warning event.
3. The run publishes OUT/annotation (`local.PUBLISHED`); `terminators_crosscheck.tsv` only when
   terminator_layer made it; annotation TSV/BED files match the existing writers byte for byte.
4. `gc_local_dev` (§9.3: 10-kb window GC − genome GC) IS produced here (M4-index line 150; GIAP_SPEC §9.3
   "Modules 1 and 3") and belongs in LAYER_COLUMNS — the tests from 461fa2a omit it: add a test (window
   centred on the gene midpoint, circular wrap on circular replicons, edge-truncated on linear ones) before
   or with the implementation.
5. annotate.nf must pass `terminators_crosscheck` to CONTEXT_LAYERS (workflow fix, implementer).
6. Integration example 5 (M3-17) stays xfail until a real run completes past CONTEXT_LAYERS; then re-check.
7. (lead, after c1ac317) `genome` is a required context_layers input (gc_local_dev needs the sequence);
   `gene_layers(..., genome=)` keyword (NA without it). "Genome GC" = whole genome, all replicons, N
   excluded; GC as a fraction, case-insensitive; window [m − 5000, m + 5000) with m = (start + end) // 2,
   modulo length on circular replicons, a circular replicon shorter than the window counted once; no A/C/G/T
   in the window → NA; column right after `region_ids`. annotate.nf stages `genome` and
   `terminators_crosscheck`; step_registry.json and the layers contract (gc_local_dev: optional float) follow.

## Lead decisions after the code review of PR #104 (changes requested)
8. (blocking) GFF3 IDs are unique per file (§9.5): every ID GIAP writes carries its feature kind as a
   prefix — `region:<region_id>`, `operon:<tu_id>`, `promoter:<id>`, `terminator:<id>` — while Bakta's own
   feature IDs are kept as given; a test on sim100k (contig named `plasmid_1`) asserts no duplicate IDs and
   that Parent= references resolve.
9. `layer_missing:<layer>` is a structured event (logger `extra={"event": …}`) that reaches
   OUT/logs/events.jsonl (and the run flags), not only a log line (§9.4, decision 2).
10. `start_codon` comes from genes.parquet; the RBS layer only adds RBS columns — a missing rbs layer never
    blanks start_codon.
11. Tests kill the surviving mutants: `crosses_origin` flag; circular wrap in `_next_same_strand` and in
    `assign_tus`; both `tu_sources` mutants; the precomputed `tx_rates` exemption in `layer_flags`;
    pseudogene phase.
12. GFF3 score for promoters is the tool's own score (Promotech / Promoter Calculator), not the BED score/1000.
13. Not blocking (backlog): annotate.nf wiring of `tool_dirs` → OUT/annotation/tools/ (§7.1.8); no producer
    of `id_map` yet, so `--annotation` GFF3/GBFF runs lack `id_map.tsv`. cai/fop NA on sim100k is expected
    (all genes "hypothetical protein": no ribosomal reference set).
