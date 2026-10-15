# genomad 1.12.0 on k12_mg1655 (M2-10, GIAP_SPEC.md §11.4)

- Command: `genomad end-to-end --cleanup --splits 8 --threads 17 genome.fa genomad /home/classes/bios/270/khoa/giap/local/db/genomad_db`
- Image: `sha256:d54652ff868af72604ea0e54c541c5c92168c7ac821a422de05dcbad65509e5c` (versions.txt: geNomad, version 1.12.0)
- Database: genomad_db v1.9 (`$GIAP_DB_DIR/genomad_db`)
- Date: 2026-10-15T07:17:44Z (run finished, UTC)
- Host: barley-01 (FarmShare; SLURM_JOB_ID not recorded; 17 CPUs; wall time 105 s)
- Recording: `/home/classes/bios/270/khoa/giap/local/recordings/genomad/k12_mg1655`
- Inputs:
  - `/home/classes/bios/270/khoa/giap/local/recordings/genomad/k12_mg1655/genome.fa`: sha256 `918f35c7bfc8a3c8d824ed281a4185903025a0361242219b6b0e0ea6ecaa14f2`
- `genome_summary/*.tsv` are the complete recorded `genomad/genome_summary/` files, verbatim; `replicons.tsv` was written by the test author (M2-18); details in ../../PROVENANCE.md

## Files (sha256sum)

```
1ecc546681777544b5f9c6fd2f2631dcd6a521af3fbf384b35fe9ddb71ad8719  genome_summary/genome_plasmid_summary.tsv
fb17b7825eb3cc325e924aaaa987a7b8abe80e16adda877f8d0ac22775fb2aab  genome_summary/genome_virus_genes.tsv
36f58beca19f70ff0c25c411496603658084f9151d6f44d5c0c62722b544b1be  genome_summary/genome_virus_summary.tsv
a1141b8a877f2e8644663c93ec95d8835079b115352e1f3e5a98c01f8a2d7660  replicons.tsv
```
