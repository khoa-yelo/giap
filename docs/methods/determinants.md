# Determinants: what explains abundance

Module 8 track R (GIAP_SPEC.md §7.8). Issues M8-07 onward; this note grows with each.

## Features

Issue M8-07, module `giap.determinants.features` (pure). `feature_frame` builds
`tables/determinants_features.tsv`; `model_matrix` encodes it for M8-08.

**Rows.** One per eligible gene (§7.3: not partial, finite `p_detect`), in §6.5 order (C1 contig
order, start, end, strand, locus tag). Ribosomal-protein genes, the CAI reference set, stay in the
table with `ribosomal_protein = true`; M8-08 excludes them from fitting and from univariate CAI
statistics. `ribosomal_protein` is the C2 column when the gene table carries it, otherwise the
§7.1.2 product / gene-name rule or the `cai_reference` gene flag.

**Columns.** Identity (`locus_tag`, `tu_id`, `detected`, `abundance_log2`, `abundance_imputed`,
`ribosomal_protein`), then the 31 features of §7.8.1 in table order. Name map (feature ← source):

| Feature | Source |
|---|---|
| `promotech_max`, `promotech_tu` | gene table `promoter_score`, `promoter_score_tu` |
| `n_promoters_tu_to_gene` | C6 `n_promoter_windows_tu` (windows ≥ `promotech_yes`, lead decision 1) |
| `mw` | gene table `mw_da` |
| `replichore_position` | M5-09 `context/gene_context.parquet` (NA without it or without an oriC) |
| `in_prophage`, `in_bgc` | C8 `region_genes` with C7 `region_type` `prophage` / `bgc`; NA without C7/C8 |
| `dist_to_tu_start_bp` | computed: bp from the TU's first start codon to this gene's, in the direction of transcription, modulo the length on circular replicons; NA for a gene in no TU |
| `nterm_residue` | computed from C2 `aa` (the Met rule below) |
| every other feature | the gene-table column of the same name |

**The Met rule (decision 19).** `nterm_residue` is residue 2 when residue 1 is M and residue 2 ∈
{A, C, G, P, S, T, V} (methionine aminopeptidase removes the initiator Met before small residues),
else residue 1; NA for proteins shorter than 2 aa or with an empty `aa` (pseudogenes). It is a
covariate only: nothing in mapping or annotation uses it (§16 excludes a recoding catalogue).

**Transforms.** None in the TSV (raw values). In `model_matrix`: `tir_ostir` and `length_aa` →
natural log (lead decision 1; values ≤ 0 → NaN); `start_codon` one-hot over {ATG, GTG, TTG, other},
`replicon_type` over {chromosome, plasmid, unknown}, `nterm_residue` over the 20 residues (every
level kept), named `<feature>_<level>`, an NA value giving NaN in every level; booleans 0/1. A
feature that is NA on every row (a missing layer) stays in the TSV and is dropped from the matrix
with all its encoded columns; the returned column → block map omits it, and a block left without
columns is absent. Run flags (`layer_missing:<layer>`) are the step's (M8-10).
