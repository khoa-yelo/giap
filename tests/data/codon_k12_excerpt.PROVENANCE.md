# codon_k12_excerpt.json.gz (M2-17 codon-usage tests)

- Source: `$GIAP_ROOT/recordings/bakta_full/k12_mg1655/bakta/k12_mg1655.json`, Bakta 1.12.2,
  full DB 6.0 (2025-02-24), E. coli K-12 MG1655 U00096.3, recorded 2026-10-14
  (sha256 `13c54a8961651053cda596ffa4cc0c611565588112275faf7f3ddc27e28bacec`).
- Built by `python -P tests/helpers/codon_excerpt.py <json> tests/data/codon_k12_excerpt.json.gz`
  (deterministic: gzip mtime 0, sorted keys).
- Content: 173 genes (all 59 ribosomal-protein genes, the EF-Tu genes, the 5 pseudogenes, sORFs
  under 46 nt, rpsT/nhaA/carA/lacZ/thrA, every 50th other CDS) with the recorded `nt`,
  `gene`, `product`, type, partial flag and start codon; all 87 tRNA anticodons.
  The CAI/Fop reference set and the tAI weights therefore equal the whole genome's.
- Size about 41 kB (R11). Bakta output data are CC0/public annotation of a public genome.
