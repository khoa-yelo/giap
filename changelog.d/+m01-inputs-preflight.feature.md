<!-- m01 -->
Inputs and preflight: `giap preflight` checks the genome, proteins, peptides, abundance table and
environment before a run and writes `preflight.json` (validated against its JSON Schema) and
`preflight.md`, with exit codes 0, 1, 2 and 3. `giap convert --from generic|diann` turns a search
engine report into GIAP's peptide table. Readers and writers cover FASTA, GFF3, BED, proBed, VCF 4.3
and TSV.
