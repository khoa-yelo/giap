# M4-06: Input contract for M4 (tables consumed from M1, M2 and M3)

Milestone M4 · track C (contract) · depends on: none (M1-04 `Issue`, M1-06 columns already merged) ·
estimate ≈ 330 code lines (+ ≈ 350 test lines) · **can start now**

M2 (annotation) and M3 (mapping) are being specified in parallel. This issue fixes, as table
schemas, exactly what M4 reads from them, so every other M4 issue is developed and tested against
simulator-generated fixtures (M4-07…09) before M2/M3 merge. File and column names below are **aligned
with the M2/M3 drafts** (M2-04 `parse_bakta`, M2-12 TUs/pairs, M2-19 regions, M2-20 `layers.parquet`,
M3-16 mapping outputs; see "Alignment notes"). Any later deviation is resolved by changing this
contract (and its validator) in one PR that the M2/M3 owners review.

## Spec sections implemented

- §6.2: "Intermediate files are Parquet under the work directory." Step rows: `parse_bakta` →
  "`genes.parquet`, `features.parquet` …"; `operon_pairs` → "`pairs_raw.parquet`, `tus.parquet`";
  `context_layers` → "`annotation/*` (§9.1), `layers.parquet`"; `loci` → "`peptide_hits.tsv`,
  `protein_hits.tsv`, `mappings.*`, `unmapped.tsv`"; `gene_table`: "genes, layers, loci, abundance →
  `genes_pre.parquet`".
- §6.1: "The only join between the mapping track and the annotation track is the Bakta locus tag.
  Analysis modules read only the gene table and the layer tables."
- §6.5: "`Interval`: … 0-based half-open. On a circular replicon, a feature crossing the origin is
  stored with `end > length` … and `crosses_origin = True`." "Sort order. All tables sort by contig
  (as in the input FASTA), then start, end, strand and ID." "Partial genes … excluded from detection
  statistics; they still appear in `genes.tsv`."
- §7.1.2 (gene model columns), §7.1.3 ("TUs are maximal runs of same-strand adjacent genes linked by
  p ≥ `uniop.pair_threshold`"), §7.1.7 (regions; "gene membership by ≥ 50 % CDS overlap"; prophage
  function subsets), §7.2.6 (locus classes, uniqueness, peptides → genes), §7.2.9 and §9.3 column
  dictionaries for `peptide_hits`, `protein_hits`, `unmapped`, `regions`, `operon_pairs`, `tus`.
- §5.3 aggregation (`peptides.norm.parquet`, M1-06 `AGG_COLUMNS`).

## The contract

Work files are addressed by the step-output keys of the M2 step registry (M2-02; M2-index "Annotation
step contract"); `load_m4_inputs` takes a key → path mapping (default layout: `OUT/work/<file>`).
Parquet files use internal 0-based half-open coordinates and Arrow `list<string>` list columns;
the M3 files are the **final** `mapping/*.tsv` tables (1-based closed, `;`-joined lists) and are
converted on read. "opt" = file or column may be absent; "NA" = nullable.

| ID | File (producer) | Columns M4 reads (type) | Invariants M4 checks |
|---|---|---|---|
| C1 | `replicons.parquet` (M2-04 `parse_bakta`) | `contig` str, `length` int, `topology` enum{circular, linear}, `kind` enum{chromosome, plasmid, unknown} | rows in FASTA order (row order = contig order); `contig` unique |
| C2 | `genes.parquet` (M2-04) | `locus_tag` str; `gene` str NA; `product` str NA; `feature_type` enum{cds, sorf, pseudogene}; `contig`; `start` int; `end` int; `strand` enum{+, -}; `partial` bool; `crosses_origin` bool; `hypothetical` bool; `aa` str NA (no trailing `*`); `nt` str (5′→3′ on the gene strand); `start_codon` str; opt `ribosomal_protein` bool, `amr` bool | `locus_tag` unique; `0 ≤ start < end`; `end ≤ length` unless `crosses_origin` (then circular and `length < end < 2·length`); `aa` NA only for `pseudogene` |
| C3 | `xrefs.parquet` (M2-04) | `locus_tag`, `db` (e.g. `COG`, `GO`, `EC`, `KEGG`), `id` | long format; M4 builds `cog_id`, `go_terms`, `kegg_ko`, `ec` lists from it |
| C4 | `annotation_meta.json` (M2-04) | `annotation_source`, `flags` | `annotation_source` copied to every `genes.tsv` row |
| C5 | `tus.parquet` (M2-12 `operon_pairs`) | `tu_id` (`TU_<n>`), `contig`, `start`, `end`, `strand`, `genes` list (transcription order) | **every C2 gene is in exactly one TU** (M2-12 partitions CDS, sORF and pseudogene; singletons are TUs) |
| C6 | `layers.parquet` (M2-20 `context_layers`) | M2-20 `LAYER_COLUMNS`: `locus_tag`, `tu_id`, `tu_position`, `tu_size`, `uniop_p_prev`, `uniop_p_next` (transcription order), `region_ids` list, `promoter_score`, `promoter_score_tu`, `n_promoter_windows_tu`, `tx_rate_max_tu`, `tx_rate_sum_tu`, `tss_distance_bp`, `terminator_conf_downstream`, `tir_ostir`, `sd_dg`, `mrna_dg`, `spacing_dg`, `standby_dg`, `spacer_bp`, `start_codon`, `cai`, `tai`, `enc`, `fop`, `gc3`, `cpb`, `at_2_10`, `rare_codons_2_10`, `flags` list | one row per C2 gene; `tu_id` agrees with C5; layer columns may be all-NA (layer missing; run flag from M2) |
| C7 | `regions.parquet` (M2-19 `region_layer`) | `region_id` (`<type>_<n>`), `region_type` enum{prophage, plasmid, bgc, crispr, amr}, `contig`, `start` NA, `end` NA, `source_tool`, `score` NA | AMR: one row per contig, NA coordinates; plasmid: whole replicon |
| C8 | `region_genes.parquet` (M2-19) | `region_id`, `locus_tag`, `overlap_frac`, `subset` enum{integrase, repressor, structural, lysis, cargo} NA, `core_biosynthetic` bool | membership decided by M2 (≥ 0.5); `subset` non-NA iff prophage |
| C9 | `peptides.norm.parquet` (M1-12) | M1-06 `AGG_COLUMNS`: `peptide`, `q_value`, `pep`, `score`, `is_decoy`, `is_contaminant`, `protein_ids` (`;`-joined str NA), `modified_peptides`, `n_rows`, `n_runs`, `intensity` (linear, NA), `source`, `search_db` | one row per stripped peptide; decoy rows included |
| C10 | `abundance.norm.parquet` (M1-12; opt) | `protein_id` (resolved), `value` float NA, `scale` enum{log2, linear} | see Q4 |
| C11 | `mapping/peptide_hits.tsv` (M3-16 `loci`) | §9.3: `peptide`, `peptide_len`, `source`, `best_q_input`, `pep_input`, `score_input`, `contig`, `start`, `end`, `strand`, `frame`, `crosses_origin`, `route` enum{exact, tolerant}, `n_diff`, `diff_types` list, `placement_q`, `length_class`, `locus_class` enum{rna_overlap, pseudogene, sorf, cds_inframe, nterm_extension, cds_outframe, antisense, intergenic}, `uniqueness` enum{unique, multi_identical, shared}, `gene` NA, `novel_candidate_id` NA, `flags` list | one row per target peptide × **accepted** locus; `gene` non-NA exactly for assigning classes; decoys never appear; every `gene` is a C2 locus tag |
| C12 | `mapping/protein_hits.tsv` (M3-16; header-only without proteins) | §9.3: `protein_id`, `protein_len_aa`, `route`, `contig`, `start`, `end`, `strand`, `rank`, `identity`, `query_coverage`, `reciprocal_overlap`, `class` enum{identical, near, divergent}, `gene` NA, `flags` list | `gene` a C2 locus tag when set |
| C13 | `mapping/unmapped.tsv` (M3-16) | §9.3 `unmapped` columns | QC rows only (M4-16) |
| C14 | `annotation/id_map.tsv` (M2-11; opt) | `source_id`, `locus_tag`, … | fills `genes.tsv` `source_id`; else `NA` |

### Alignment notes (differences from the first M4 draft, and gaps)

- M2 names: `aa`/`nt` (not `aa_seq`/`nt_seq`), `replicons.length`, `pairs_raw.parquet` with
  `gene_a` = earlier gene in **genome** order (M4 does not read it: `uniop_p_prev/next` come from
  C6), `region_genes.parquet` with `subset`.
- §9.3 `genes.tsv` columns with **no producer** in the M2 drafts: `gc_local_dev` (context) and
  `cog_category` (needs the COG 2020 `cog-20.def.tab` resource from M2-09 to map `COG####` → letters).
  M4 writes `NA` for `gc_local_dev`; `cog_category` is filled by M4-15 from the COG table when
  `giap setup --cog` has run, else `NA` (see Q6).
- Pseudogene peptide assignment: M3-14 Q3 proposes **no** assignment; this contract's Q3 proposes
  yes. The lead must pick one; M4 code handles both (it trusts C11 `gene`).

Optional input: `--replicon-depth` TSV (`contig`, `mean_depth`, §5.3), read in M4-18.

## Public interface

`src/giap/genetable/contract.py` (pure)
```python
@dataclass(frozen=True, slots=True)
class Column:
    name: str; dtype: Literal["str", "int", "float", "bool", "list"]
    nullable: bool = False; enum: tuple[str, ...] | None = None; optional: bool = False

TABLES: Mapping[str, tuple[Column, ...]]          # keys "replicons", "genes", "xrefs", "tus", "layers", "regions",
                                                  # "region_genes", "peptides", "abundance", "peptide_hits", "protein_hits", "unmapped", "id_map"

def validate_table(name: str, df: pd.DataFrame) -> list[Issue]: ...      # columns, dtypes, enums, NA
def validate_inputs(t: M4Inputs) -> list[Issue]: ...                      # cross-table invariants
```
`src/giap/io/work_tables.py` (the only file I/O)
```python
@dataclass(frozen=True, slots=True)
class M4Inputs:
    replicons: pd.DataFrame; genes: pd.DataFrame; xrefs: pd.DataFrame; annotation_source: str
    tus: pd.DataFrame; layers: pd.DataFrame; regions: pd.DataFrame; region_genes: pd.DataFrame
    peptides: pd.DataFrame | None; abundance: pd.DataFrame | None
    peptide_hits: pd.DataFrame; protein_hits: pd.DataFrame; unmapped: pd.DataFrame
    id_map: pd.DataFrame | None

M4_INPUT_KEYS: tuple[str, ...]     # registry keys: "replicons", "genes", "xrefs", "annotation_meta", "tus", "layers", …
def load_m4_inputs(paths: Mapping[str, Path]) -> tuple[M4Inputs, list[Issue]]: ...
def default_paths(work_dir: Path) -> dict[str, Path]: ...    # fixture / local-executor layout
```
`Issue` is M1's shared type (`src/giap/io/issues.py`); codes are `contract.<table>.<rule>`.
Errors: missing required file or column, wrong dtype, bad enum, invariant broken. Warnings:
extra columns (ignored), missing optional file. `load_m4_inputs` never raises on content; it returns
issues, and the `gene_table` step (M4-15) stops on any error.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | valid fixture written by `giapsim.write_contract_tables` (M4-07) for a 200-gene layout | `load_m4_inputs` → no `error` issues; `len(genes) == 200` |
| 2 | C2 with two rows `locus_tag = GSIM_00005` | error `contract.genes.locus_tag_unique`, value lists `GSIM_00005` |
| 3 | C2 row `feature_type = "CDS"` (upper case) | error `contract.genes.enum`, column `feature_type` |
| 4 | gene `GSIM_00010` absent from every C5 `genes` list | error `contract.tus.gene_coverage` naming `GSIM_00010` |
| 5 | gene in two TUs | error `contract.tus.gene_coverage` |
| 6 | circular 50,000-bp contig, gene `start=49_900, end=50_200, crosses_origin=true` | valid; same with `topology=linear` → error `contract.genes.origin` |
| 7 | C11 row with `uniqueness = unique`, `locus_class = cds_inframe`, `gene = NA` | error `contract.peptide_hits.assignment` (an assigning class needs a gene) |
| 8 | C11 row `gene = GSIM_99999` not in C2 | error `contract.peptide_hits.unknown_gene` |
| 9 | C6 with every promoter column NA | valid (layer missing is not a contract error) |
| 10 | C9 with `protein_ids` as a list column instead of a `;`-joined string | error `contract.peptides.dtype` |

## Edge cases

- Header-only C11 (no peptide placed) is valid.
- A replicon with zero genes (tiny plasmid) is valid; it gets no TUs.
- Pseudogene with `aa` NA must have `nt` (M4-10 translates it).
- `gene` NA is common; `product` NA allowed; C14 absent → `source_id` NA.
- M3 TSVs are read with M1-04 `read_tsv` and converted to 0-based (`start − 1`); `;` lists split.
- Dtype check accepts pandas nullable dtypes (`Int64`, `boolean`, `string`) and numpy equivalents.

## Acceptance criteria

- [ ] `tests/unit/genetable/test_contract.py`: examples 1–10, one test each; every `Column` in
  `TABLES` is exercised by at least one missing-column test (parametrised).
- [ ] `tests/contract/test_m4_input_contract.py` validates the M4-07 simulator output (example 1) and,
  once M2/M3 merge, their real outputs (marked `pending_sample` until then; M4-21 unmarks it).
- [ ] `TABLES` is the single source of the column names; M4-07 (simulator) imports nothing from
  `giap` and a test asserts its column lists equal `TABLES` (catches drift both ways).
- [ ] File I/O only in `io/work_tables.py`; `contract.py` imports no I/O.

## Methods note

None (no statistics). The contract is documented in the module docstring and copied into
`docs/specs/M4-index.md` (this table).

## Questions for the lead

1. Must every gene, including pseudogenes and genes UniOP never saw, belong to a TU? **Proposed
   default:** yes (C5 invariant). M2-12 already partitions CDS, sORF and pseudogene into TUs; M4's
   TU-block resampling and cross-fitting rely on it.
2. M4 reads M3's **final** `mapping/*.tsv` (1-based, converted on read) because M3-16 writes no
   Parquet twin? **Proposed default:** yes; no change asked of M3.
3. Peptides placed in a pseudogene (`locus_class = pseudogene`): assign to the pseudogene?
   **Proposed default (M4):** yes, since pseudogenes are eligible rows with their own covariates
   (§7.3.3, "Pseudogenes. Covariates come from …"); **M3-14 Q3 proposes no** (§7.2.6 names only
   `cds_inframe`, `sorf`, `nterm_extension`). Conflict for the lead. If "no", pseudogenes can only be
   detected through proteins, and M4 should drop them from the eligible set to keep O − E honest.
4. `abundance.norm.parquet` columns are not fixed by M1-12. **Proposed default:** C10 as above
   (`protein_id`, `value`, `scale`); M1-12's writer adopts it.
5. Function columns from `xrefs.parquet` (long) rather than list columns in `genes.parquet`?
   **Proposed default:** yes (M2-04's design); M4-15 pivots them.
6. `cog_category` and `gc_local_dev` have no producer in the M2 drafts. **Proposed default:**
   `cog_category` from the COG 2020 definition table when installed (M4-15), else `NA` (M4-16 then
   skips COG rows with a warning); `gc_local_dev` stays `NA` until an issue (M2 or M5) defines it.
