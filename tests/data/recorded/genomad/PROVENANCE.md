# geNOMAD 1.12.0 recorded runs (M2-18; GIAP_SPEC.md §7.1.1)

- Tool: geNOMAD 1.12.0 (`versions.txt`: "geNomad, version 1.12.0", image
  `sha256:d54652ff…509e5c`), database genomad_db v1.9 (non-commercial licence; R5: run as a
  separate process only, never imported).
- Recorded on FarmShare (barley-01, 17 CPUs), 2026-10-15, into
  `$GIAP_ROOT/recordings/genomad/{sim100k,k12_mg1655,nissle1917}/`.
- Command (`command.txt`): `genomad end-to-end --cleanup --splits 8 --threads 17 genome.fa genomad
  <db>/genomad_db`; output prefix is therefore `genome` (the Nextflow GENOMAD process stages the
  genome as `genome.fa` too).
- Inputs (`genome.fa` in each recording): sim100k sha256 `369d1da8…f2a7f` (chr1 100,000 bp,
  plasmid_1 8,000 bp); k12_mg1655 `918f35c7…f14f2` (U00096.3, 4,641,652 bp); nissle1917
  `f64e01e8…85807` (CP171242.1 5,051,730 bp; plasmids CP171243.1 4,523, CP171244.1 3,173,
  CP171245.1 5,514 bp).

## Files (`1.12.0/<sample>/`)

`genome_summary/genome_{virus_summary,plasmid_summary,virus_genes}.tsv` are the **complete**
recorded files, copied verbatim (they are small). sha256 of the sources
(`$GIAP_ROOT/recordings/genomad/<sample>/genomad/genome_summary/…`):

| sample | virus_summary | plasmid_summary | virus_genes |
|---|---|---|---|
| sim100k | `a14c6fdd1087ae0f1a3c54363f49ab02d091ab0cf3e6d1329b0d26e3ad60588b` (header only) | `1ecc546681777544b5f9c6fd2f2631dcd6a521af3fbf384b35fe9ddb71ad8719` (header only) | `5a1779721b35afa8ebd7a9efe64ed2ff356491187dcffd89f003de3800062a69` (header only) |
| k12_mg1655 | `36f58beca19f70ff0c25c411496603658084f9151d6f44d5c0c62722b544b1be` (5 proviruses) | `1ecc546681777544b5f9c6fd2f2631dcd6a521af3fbf384b35fe9ddb71ad8719` (header only) | `fb17b7825eb3cc325e924aaaa987a7b8abe80e16adda877f8d0ac22775fb2aab` (145 genes) |
| nissle1917 | `41bafd2b2c3a7aa5c8028e445bb0cd19be0bacf8107178e721c6d939177b9199` (3 proviruses) | `50b4b14d99d329889a59869a599b52c9f1e3bff559c344ef745c928e2b8588ea` (CP171244.1, CP171245.1) | `7426e9f756b6049c0cced19a9b172b765674b7ccbeadf5f294e827cfc7ab09c3` (133 genes) |

`replicons.tsv` (written by the test author): `contig, length, topology` of every input contig,
lengths from `genome.fa`, topology from the Bakta GBFF records the same genomes were annotated
with (sim100k circular; K-12 and Nissle linear in these runs).

## What the recordings decide (M2-18 Q2)

geNOMAD 1.12.0 `virus_genes.tsv` `start`/`end` are **absolute, 1-based closed contig
coordinates**, not relative to the provirus: the first gene of
`U00096.3|provirus_563848_584430` (`…_520`) starts at 563848, the provirus start, and every gene
lies inside its provirus only under the absolute reading. Gene numbers (`_520`) are genome-wide
too. Spec example 4 (`start 201` → `12200`) is therefore rewritten in the tests with absolute
coordinates (`start 12201`). No recording has a whole-contig virus call, an FDR value
(all `NA`) or a provirus near a contig end; those cases are covered by hand-written inputs.
