# Bakta JSON excerpts for the M2-19 region tests

- Source: the recorded Bakta 1.12.2 full-DB runs (M2-10) under
  `$GIAP_ROOT/recordings/bakta_full/<sample>/bakta/<sample>.json`, the same annotations whose
  GBFF the recorded antiSMASH 8.0.4 runs used (locus tags `KOIAIB_…`, `JMJJOG_…`) and whose
  genomes the recorded geNOMAD 1.12.0 runs used (`tests/data/recorded/{genomad,antismash}/`).
- sha256 of the full source JSON:
  - k12_mg1655 `13c54a8961651053cda596ffa4cc0c611565588112275faf7f3ddc27e28bacec`
  - nissle1917 `49dbc4dc51e38443e822bf26481a7771c2d95f59315c5cb648cd23808809256f`
- Written by `regions_make_excerpts.py` (in this directory; run by hand with `$GIAP_ROOT` set) on
  2026-10-15, trimmed for R11 (M2-index lead decision 4).

## Selection

Every CDS/sORF overlapping a geNOMAD provirus or an antiSMASH area plus the nearest gene on each
side outside it; every gene within 3 kb of a CRISPR array; every gene on a contig shorter than
10 kb (Nissle CP171243.1, CP171244.1 = pMUT1, CP171245.1 = pMUT2); every gene with an
`amrfinder` expert hit; every CRISPR feature. Kept fields: `type, sequence, start, stop, strand,
locus, gene, product, pseudo` and `expert` reduced to `type`/`source`; CRISPR features keep
`id`. `sequences[]` keeps `id, type, topology, length` (no `nt`; tests pass the replicons).

So every region the recordings define (prophage, BGC, CRISPR with cas genes, plasmid, AMR) has
all its member genes in the excerpt; `tests/contract/test_regions_contract.py` checks that the
excerpt and the full JSON give identical regions when the full recordings are present.

| sample | features kept | sequences (Bakta `type`, topology) |
|---|---|---|
| k12_mg1655 | 286 | U00096.3 chromosome, linear |
| nissle1917 | 450 | CP171242.1 chromosome; CP171243.1, CP171244.1, CP171245.1 `contig`; all linear |
