# genomad 1.12.0 on sim100k (M2-10, GIAP_SPEC.md §11.4)

- Command: `genomad end-to-end --cleanup --splits 8 --threads 17 genome.fa genomad /home/classes/bios/270/khoa/giap/local/db/genomad_db`
- Image: `sha256:d54652ff868af72604ea0e54c541c5c92168c7ac821a422de05dcbad65509e5c` (versions.txt: geNomad, version 1.12.0)
- Database: genomad_db v1.9 (`$GIAP_DB_DIR/genomad_db`)
- Date: 2026-10-15T07:15:59Z (run finished, UTC)
- Host: barley-01 (FarmShare; SLURM_JOB_ID not recorded; 17 CPUs; wall time 66 s)
- Recording: `/home/classes/bios/270/khoa/giap/local/recordings/genomad/sim100k`
- Inputs:
  - `/home/classes/bios/270/khoa/giap/local/recordings/genomad/sim100k/genome.fa`: sha256 `369d1da812bd925afe260cad42a51eb1181d4f358f4d94db6ee38e05516f2a7f`
- `genome_summary/*.tsv` are the complete recorded `genomad/genome_summary/` files, verbatim; `replicons.tsv` was written by the test author (M2-18); details in ../../PROVENANCE.md

## Files (sha256sum)

```
1ecc546681777544b5f9c6fd2f2631dcd6a521af3fbf384b35fe9ddb71ad8719  genome_summary/genome_plasmid_summary.tsv
5a1779721b35afa8ebd7a9efe64ed2ff356491187dcffd89f003de3800062a69  genome_summary/genome_virus_genes.tsv
a14c6fdd1087ae0f1a3c54363f49ab02d091ab0cf3e6d1329b0d26e3ad60588b  genome_summary/genome_virus_summary.tsv
ccb4832800acea547b8ac2d045719beaba8c59a40e003215132ac8856bb8817a  replicons.tsv
```
