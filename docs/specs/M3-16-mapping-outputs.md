# M3-16: Module 2 outputs and their table schemas

Milestone M3 · depends on: M3-10, M3-11, M3-14, M3-15; M1-04 (TSV, BED, proBed), M1-05 (GFF3) ·
creates `src/giap/mapping/outputs.py`, `src/giap/schemas/tables/{peptide_hits,protein_hits,
unmapped,mapping_fdr}.schema.json` · pure Python · estimate ≈ 360 code lines (schemas excluded)

## Spec sections implemented

- §7.2.9 "`mapping/peptide_hits.tsv`, `protein_hits.tsv`, `mappings.probed` (proBed 1.0,
  BED12+13), `mappings.bed`, `mappings.gff3`, `unmapped.tsv` and `mapping_fdr.tsv`. Schemas are in
  §9.3."
- §9.3, exactly these columns and this order:
  - **`peptide_hits.tsv`** ("One row per peptide × accepted locus. Rejected placements stay in
    `unmapped.tsv`."): `peptide`, `peptide_len`, `source`, `best_q_input`, `pep_input`,
    `score_input`, `contig`, `start`, `end`, `strand`, `frame`, `crosses_origin`, `route`
    enum{exact, tolerant}, `n_diff`, `diff_types` (list), `placement_q`, `length_class`,
    `locus_class`, `uniqueness`, `gene`, `novel_candidate_id`, `flags`.
  - **`protein_hits.tsv`**: `protein_id`, `protein_len_aa`, `route` enum{miniprot, tblastn},
    `contig`, `start`, `end`, `strand`, `rank`, `identity`, `query_coverage`,
    `reciprocal_overlap`, `class` enum{identical, near, divergent}, `gene`, `flags`.
  - **`unmapped.tsv`**: `sequence_id`, `sequence`, `kind`, `label`, `reason`,
    `contaminant_match`.
  - **`mapping_fdr.tsv`**: `length_class`, `n_diff`, `targets`, `decoys`, `pi_wrong`,
    `fdr_stratum`, `n_accepted`, `decoy_type`; "Per-m total rows carry `length_class = all`. De
    novo precision rows carry `section = denovo_precision` with `score_bin`, `precision`,
    `ci_low`, `ci_high`, `c_value`, `flags`."
- §9.2 conventions (TSV, `NA`, `%.6g`, `true/false`, `;` lists, sorted, units in names), §6.5
  (TSV/GFF3 1-based closed with `end > length` for origin-crossing; BED/proBed 0-based,
  origin-crossing split into `_part1`/`_part2`), §9.5 Frictionless schemas with a `coordinates`
  property, §6.6 no timestamps or absolute paths in `mapping/`.

## Public interface

```python
PEPTIDE_HITS_COLUMNS: tuple[str, ...]; PROTEIN_HITS_COLUMNS: tuple[str, ...]
UNMAPPED_COLUMNS: tuple[str, ...]; MAPPING_FDR_COLUMNS: tuple[str, ...]
def write_mapping_outputs(outdir: Path, peptide_hits: pd.DataFrame, protein_hits: pd.DataFrame,
                          unmapped: Sequence[UnmappedRow], fdr_rows: pd.DataFrame,
                          replicons: Mapping[str, Replicon], genome_sha256: str) -> None: ...
```

- `mapping_fdr.tsv` column order: the eight §9.3 columns, then `section` (`placement` |
  `denovo_precision`), `score_bin`, `precision`, `ci_low`, `ci_high`, `c_value`, `flags`; columns
  not used by a row are `NA`.
- **Sort:** `peptide_hits` and `protein_hits` by contig (FASTA order), start, end, strand, then
  `peptide`/`protein_id`; `unmapped` by kind, sequence_id; `mapping_fdr` by section, decoy_type,
  n_diff, length class order (`all` last).
- **BED (`mappings.bed`, BED12):** one record per accepted peptide locus and per protein hit;
  `name` = peptide or protein ID; score 1000 for exact peptides and `identical` proteins, 500
  otherwise; one block, or two for frameshift peptides (`blocks` column, M3-08; for `δ = −1` the
  first block is shortened by 1 nt so blocks do not overlap).
- **proBed (`mappings.probed`):** peptides only; `proteinAccession` = `gene` (or `.`),
  `peptideSequence` = peptide, `uniqueness` = `unique` | `not-unique[same-set]`
  (multi_identical) | `not-unique[conflict]` (shared), `psmScore` = `score_input`, `fdr` =
  `placement_q`, `genomeReferenceVersion` = `sha256:<genome hash>` (Q2).
- **GFF3 (`mappings.gff3`):** source `GIAP`, type `nucleotide_to_protein_match` (SO:0000568) for
  peptides and proteins; attributes `ID`, `Name`, `route`, `n_diff`, `diff_types`, `placement_q`,
  `locus_class`, `uniqueness`, `gene`, `flags` (peptides) or `identity`, `query_coverage`,
  `class`, `rank` (proteins); frameshift peptides as two `match_part` children.
- **Schemas:** one Frictionless Table Schema per TSV (types, enums from §9.3/§9.4, `NA` as
  missing value, `coordinates: 1-based-closed`).

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | exact hit `Interval(chr1,130,160,+)`, `PEPTIDEKLLR`, q 0.002 | TSV row `start=131`, `end=160`, `strand=+`, `frame=2`, `crosses_origin=false`, `route=exact`, `n_diff=0`, `diff_types=NA` |
| 2 | origin-crossing hit `Interval(chr1,9990,10020,+)` on L = 10,000 | TSV `start=9991`, `end=10020`, `crosses_origin=true`; BED `chr1 9990 10000 PEP…_part1` and `chr1 0 20 PEP…_part2` |
| 3 | tolerant hit with `diff_types=("substitution","swap")` | TSV `diff_types=substitution;swap`, `n_diff=2` |
| 4 | frameshift hit of M3-08 example 1 | BED12 `blockCount=2`, `blockSizes=18,18`, `blockStarts=0,19` |
| 5 | `placement_q = 0.000123456789` | written `0.000123457` |
| 6 | `mapping_fdr` per-m total | `length_class=all`, `section=placement` |
| 7 | de novo bin row | `section=denovo_precision`, `score_bin=3`, placement columns `NA` |
| 8 | empty inputs | each TSV has its header only; BED/proBed have only the header line |
| 9 | the four TSVs of the 50-kb golden fixture | `frictionless validate` passes against the schemas |

## Edge cases

- Two runs with the same inputs and seed: byte-identical files (TE7).
- `gene` NA → `NA` in TSV, `.` in proBed, attribute omitted in GFF3.
- Peptide with `I` in the sequence: written as in the input (I/L kept as written, §5.3).

## Acceptance criteria

- [ ] `tests/unit/test_mapping_outputs.py`: examples 1–8.
- [ ] Column tuples equal the §9.3 lists exactly (a test pins them).
- [ ] Schema test (`tests/unit/test_schemas.py` extended): example 9; enums match §9.3/§9.4.
- [ ] GFF3 round trip with M1-05's reader; proBed has 25 fields per line.

## Questions for the lead

1. GFF3 feature type. **Proposed default:** `nucleotide_to_protein_match` (SO:0000568), with
   `match_part` children for frameshift peptides.
2. proBed `fdr`: placement q or the engine q? **Proposed default:** `placement_q` (it qualifies
   this genomic placement); the engine q is in `peptide_hits.tsv` as `best_q_input`.
3. Should `mappings.*` include protein hits? **Proposed default:** BED and GFF3 yes; proBed is a
   peptide format, so peptides only.
