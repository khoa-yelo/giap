# M8-07: `determinants_features.tsv` — the five feature blocks, one row per eligible gene

Milestone M8 · track R · **pass 1** · depends on: nothing new (merged gene table, C2, C5, C6, C7/C8,
M5-09 gene context) · estimate ≈ 280 code lines (+ ≈ 380 test lines) · not statistical

## Spec sections implemented

- §7.8.1 the feature table verbatim:

  | Block | Features | Source |
  |---|---|---|
  | Transcription | `promotech_max`; `promotech_tu`; `tx_rate_max_tu`, `tx_rate_sum_tu`; `n_promoters_tu_to_gene`; `tu_position`; `dist_to_tu_start_bp` | Promotech, Promoter Calculator, UniOP |
  | Translation initiation | `tir_ostir` (log); `sd_dg`; `spacer_bp`; `start_codon` (one-hot ATG, GTG, TTG, other); `mrna_dg` | OSTIR with the strain's 16S 3′ end |
  | Codon usage | `cai` (ribosomal-protein reference); `tai`; `enc`; `fop`; `gc3`; `cpb`; `at_2_10`; `rare_codons_2_10` | codon-bias |
  | Protein | `length_aa` (log); `mw`; `pi`; `gravy`; `nterm_residue` (Met removed when residue 2 ∈ {A, C, G, P, S, T, V}); `tm_segments` | sequence |
  | Context | `replicon_type`; `replichore_position`; `in_prophage`; `in_bgc`; `gc_local_dev` (10-kb window GC − genome GC) | Modules 1 and 3 |
- §7.8.2 step 1 "Rows … Ribosomal-protein genes, the CAI reference set, are excluded from model
  fitting and from univariate CAI statistics, and are reported separately."
- §7.8.3 "Unit: the Met-processing rule"; §9.1 `tables/determinants_features.tsv`; §7.3 eligibility;
  §7.3.3 ("Expression determinants (CAI, promoter, RBS) are excluded [from detectability] so
  Module 8 is not circular").

## Public interface

`src/giap/determinants/features.py` (pure)
```python
BLOCKS: Mapping[str, tuple[str, ...]]      # block -> model columns (after encoding)
FEATURE_COLUMNS: tuple[str, ...]           # the TSV columns, §7.8.1 order

def nterm_residue(aa: str) -> str | None: ...
def dist_to_tu_start(genes: Sequence[Gene], tus: Sequence[Tu], replicons) -> dict[str, int]: ...
def feature_frame(genes_tbl: pd.DataFrame, genes: Sequence[Gene], tus: Sequence[Tu],
                  layers: pd.DataFrame, regions: pd.DataFrame | None,
                  gene_context: pd.DataFrame | None, replicons) -> pd.DataFrame: ...
def model_matrix(features: pd.DataFrame) -> tuple[pd.DataFrame, Mapping[str, str]]: ...
                                           # encoded columns, column -> block
```

1. **Rows.** One per eligible gene (§7.3; `giap.genetable.eligible()` from M8-03), in §6.5 order.
   Identity columns first: `locus_tag`, `tu_id`, `detected`, `abundance_log2`, `abundance_imputed`,
   `ribosomal_protein` (from C2; the CAI reference set), then the features in §7.8.1 order.
2. **Column sources** (spec name ← existing column):
   `promotech_max` ← `promoter_score`; `promotech_tu` ← `promoter_score_tu`;
   `tx_rate_max_tu`, `tx_rate_sum_tu`, `tu_position` ← gene table; `n_promoters_tu_to_gene` ←
   C6 `n_promoter_windows_tu` (windows ≥ `annotation.promoters.promotech_yes` from the TU's first
   gene to this gene — the same definition); `tir_ostir`, `sd_dg`, `spacer_bp`, `start_codon`,
   `mrna_dg`, the eight codon columns, `pi`, `gravy`, `tm_segments`, `replicon_type`,
   `gc_local_dev` ← gene table; `mw` ← `mw_da`; `replichore_position` ← M5-09
   `context/gene_context.parquet` (NA without it or without an oriC).
3. **Computed here.**
   - `dist_to_tu_start_bp`: bp from the TU's first start codon to this gene's start codon, in the
     direction of transcription (0 for the first gene; circular replicons wrap, as M2-12's
     `intergenic_distance` does). Genes in no TU (should not occur after M2-12) → NA.
   - `nterm_residue` = residue 2 when residue 1 is M and residue 2 ∈ {A, C, G, P, S, T, V}, else
     residue 1 (decision 19: a feature only). Proteins shorter than 2 aa → NA. The input is C2
     `aa` (Bakta's protein; pseudogenes: M4-10's in-frame translation is not re-derived here — NA
     when `aa` is empty).
   - `in_prophage`, `in_bgc` from C8 `region_genes` (types `prophage`, `bgc`); without C7/C8 →
     NA with run flag `layer_missing:regions`.
   - Transforms in the **TSV**: none (the TSV holds raw values). Transforms in `model_matrix`:
     `tir_ostir` → ln (≤ 0 → NA), `length_aa` → ln; one-hot `start_codon` {ATG, GTG, TTG, other}
     and `replicon_type` {chromosome, plasmid, unknown}; `nterm_residue` one-hot over the 20
     residues (rare levels kept: LightGBM and the elastic net both cope).
4. **Missing layers.** A layer column that is NA for every row (e.g. no Promotech) is kept in the
   TSV and **dropped from `model_matrix`**, reported in the returned metadata and as run flag
   `layer_missing:<layer>`; a block with no column left is absent from block statistics (M8-08).
5. **`BLOCKS`**: transcription, translation, codon, protein, context — matching
   `determinants_attribution.tsv`'s `shap_transcription … shap_context`.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| gene `aa` = "M" | `nterm_residue` NA |
| residue 1 not M (non-ATG start translated as V/L by the annotation) | `nterm_residue` = residue 1 (the rule needs an initiator Met) |
| first gene of a TU | `dist_to_tu_start_bp` = 0 |
| TU crossing the origin | distance unwrapped (example 6) |
| linear replicon, no oriC | `replichore_position` NA |
| no C7/C8 | `in_prophage`, `in_bgc` NA + `layer_missing:regions` |
| `start_codon` not in {ATG, GTG, TTG} | one-hot `other` |
| ribosomal protein | row present, `ribosomal_protein = true` (M8-08 excludes it from fitting) |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `aa` = "MAKV…" | `nterm_residue` = **A** (Met removed: A ∈ {A, C, G, P, S, T, V}) |
| 2 | `aa` = "MKLV…" | **M** (K not in the set) |
| 3 | `aa` = "MP…" | **P** |
| 4 | `aa` = "VKL…" | **V** (no initiator M) |
| 5 | TU on `+` with genes at starts 1000 and 1450 (0-based) | second gene `dist_to_tu_start_bp` = **450**; first 0 |
| 6 | TU on `−` with genes [5000, 5600) and [4300, 4950) | first gene in transcription = the right one (start codon at 5599); second at 4949 → **650** |
| 7 | TU on a circular 10,000-bp contig, `+`, genes at 9800 and 200 | **400** (wraps) |
| 8 | gene in `prophage_1` | `in_prophage = true`, `in_bgc = false` |
| 9 | Promotech layer missing | `promotech_max`, `promotech_tu`, `n_promoters_tu_to_gene` NA in the TSV, absent from `model_matrix`; flag `layer_missing:promoters` |
| 10 | `tir_ostir` = 0 | TSV 0; model value NA (ln undefined) |

## Acceptance criteria

- [ ] `tests/unit/determinants/test_features.py` covers examples 1–10 (examples 1–4 are the
  §7.8.3 "Met-processing rule" unit test) and every degenerate row.
- [ ] Property: `model_matrix` columns map to exactly one block each; mirroring the genome
  (strand flip + coordinate mirror) leaves `dist_to_tu_start_bp` unchanged.
- [ ] Conformance: on the merged sim100k fixture, every eligible gene has a row and the feature
  columns equal the gene-table columns they come from.
- [ ] `docs/methods/determinants.md` §"Features": the name map of step 2, the transforms, the Met
  rule and decision 19.

## Hardening owed (pass 2)

- A golden `determinants_features.tsv` for the 50-kb fixture.
- Check on the K-12 recording that `nterm_residue` matches known N-terminal processing for a few
  well-characterised proteins (descriptive, not gated).

## Questions for the lead

1. `nterm_residue` vs §16 "No: recoding and start catalogue (… Met removal)": §7.8.1 asks for this
   rule as a covariate, and §16 excludes a catalogue for mapping. **Proposed default:** feature
   only; nothing in mapping or annotation uses it.
2. `n_promoters_tu_to_gene` = C6 `n_promoter_windows_tu` (windows at or above `promotech_yes`)?
   **Proposed default:** yes; it is the merged M2-13 definition.
3. Logs natural (proposed) or base 10? (Trees are invariant; the elastic net is standardised.)
