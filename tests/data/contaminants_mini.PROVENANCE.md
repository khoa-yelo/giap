# contaminants_mini.faa (M3-15 test fixture)

A small "library-like" protein FASTA for the M3-15 unit, property and TE5-proxy tests; it is NOT
an excerpt of the Frankenfield 2022 universal contaminant library (that library is downloaded by
`giap setup --contaminants` and never vendored into tests).

- Entries: porcine trypsin (P00761), human keratin 1 (P04264), human keratin 10 (P13645), bovine
  serum albumin (P02769); UniProt-style headers so `read_protein_fasta` resolves the IDs to the
  accessions.
- Sequences: hand-written test sequences, written by the test author (2026-10-14) to resemble the
  named proteins and to contain the well-known contaminant peptides used by the tests
  (`VATVSLPR` in trypsin; `SLNNQFASFIDK` in keratin 1 only). They are not guaranteed to equal the
  UniProt entries residue for residue and must not be used as reference data.
- Licence: written for this repository (no third-party content).
