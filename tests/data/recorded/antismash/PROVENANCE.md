# antiSMASH 8.0.4 recorded runs (M2-18; GIAP_SPEC.md §7.1.1)

- Tool: antiSMASH 8.0.4 (AGPL-3.0; R5: separate process only, never imported), image
  `sha256:d3f390e6…bda29d` (`versions.txt`).
- Recorded on FarmShare (barley-01, 17 CPUs), 2026-10-15, into
  `$GIAP_ROOT/recordings/antismash/{sim100k,k12_mg1655,nissle1917}/`.
- Command (`command.txt`): `antismash --genefinding-tool none --minimal --cpus 17 --databases
  <db>/antismash --output-dir antismash bakta.gbff` (the recorder passes `--databases`; the spec
  command and `workflow/modules/local/antismash/main.nf` do not).
- Input `bakta.gbff`: the recorded Bakta 1.12.2 full-DB GBFF
  (`$GIAP_ROOT/recordings/bakta_full/<sample>/bakta/<sample>.gbff`; sha256 below). Locus tags are Bakta's (`KOIAIB_…`, `JMJJOG_…`).
- Output JSON is named after the input: `antismash/bakta.json`.

| sample | full `antismash/bakta.json` sha256 | records | areas |
|---|---|---|---|
| sim100k | `317521dd8c7e6f4df7c048d83cc4aceca6592d87b253e5e2fdd7f99424d6e17a` | chr1, plasmid_1 (circular) | 0 |
| k12_mg1655 | `fe8cb6eb9fa01c3d5e273fe69182c4ba7457364be4a655e69249029240d32f3d` | U00096.3 (linear) | 3 |
| nissle1917 | `59b6c76952fc76ea886b347a2f6b383054e12831aacffd3a59c6a5f596c7ff92` | CP171242.1 + 3 plasmids (linear) | 6 (all on CP171242.1) |

GBFF sha256: sim100k `c028b1c43b22260713a77796c6ce79d855773b682a68c2f4a2fe36124481c6fd`,
k12_mg1655 `bc0ad5c0ece48585e913f96ec230f596c9b75238dbf2417451348768b71871d7`, nissle1917
`2d9d799e63173a5b0bb749c39752d6a3f9112d8e8492944b8b407dbc9a27fea3`.

## Files (`8.0.4/<sample>/bakta.json`): trimmed excerpts (M2-index lead decision 4, R11)

Top-level `version`, `input_file`, `schema` kept; per record `id`, `name`,
`annotations.topology`, the complete `areas` list (verbatim, incl. `protoclusters`,
`candidates`, `subregions`), and of `features` only: every `region` / `protocluster` /
`proto_core` / `cand_cluster` feature, every `CDS` overlapping an area, and the nearest `CDS` on
each side outside every area. `seq` and all other record keys are dropped, as are the CDS
qualifiers `translation`, `NRPS_PKS`, `sec_met_domain`. Kept features stay in genome order and
keep their original `location` strings, including the fuzzy `[3082708:>3082864](-)` CDS inside
Nissle area 3. `replicons.tsv`: `contig, length, topology` from the full records.

Core biosynthetic genes (CDS in the area with `gene_kind` containing exactly `biosynthetic`, not
`biosynthetic-additional`), checked by the test author on the full JSON:
K-12 area 1 terpene-precursor `KOIAIB_00460`; area 2 NRP-metallophore/NRPS `KOIAIB_00640`,
`KOIAIB_00647`, `KOIAIB_00650`; area 3 azole-containing-RiPP `KOIAIB_00973`, `KOIAIB_00976`.

No recorded area touches a contig end (`contig_edge` is `False` everywhere) and no record is
circular with an area; origin-crossing areas (`start > end`) are covered by hand-written JSON.
