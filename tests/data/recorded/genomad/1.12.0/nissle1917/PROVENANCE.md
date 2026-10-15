# genomad 1.12.0 on nissle1917 (M2-10, GIAP_SPEC.md §11.4)

- Command: `genomad end-to-end --cleanup --splits 8 --threads 17 genome.fa genomad /home/classes/bios/270/khoa/giap/local/db/genomad_db`
- Image: `sha256:d54652ff868af72604ea0e54c541c5c92168c7ac821a422de05dcbad65509e5c` (versions.txt: geNomad, version 1.12.0)
- Database: genomad_db v1.9 (`$GIAP_DB_DIR/genomad_db`)
- Date: 2026-10-15T07:19:32Z (run finished, UTC)
- Host: barley-01 (FarmShare; SLURM_JOB_ID not recorded; 17 CPUs; wall time 108 s)
- Recording: `/home/classes/bios/270/khoa/giap/local/recordings/genomad/nissle1917`
- Inputs:
  - `/home/classes/bios/270/khoa/giap/local/recordings/genomad/nissle1917/genome.fa`: sha256 `f64e01e8cb859c36221a97f730e6e40067c5abbdf5a861c62982c91d78485807`
- `genome_summary/*.tsv` are the complete recorded `genomad/genome_summary/` files, verbatim; `replicons.tsv` was written by the test author (M2-18); details in ../../PROVENANCE.md

## Files (sha256sum)

```
50b4b14d99d329889a59869a599b52c9f1e3bff559c344ef745c928e2b8588ea  genome_summary/genome_plasmid_summary.tsv
7426e9f756b6049c0cced19a9b172b765674b7ccbeadf5f294e827cfc7ab09c3  genome_summary/genome_virus_genes.tsv
41bafd2b2c3a7aa5c8028e445bb0cd19be0bacf8107178e721c6d939177b9199  genome_summary/genome_virus_summary.tsv
48ae97a1fe49a71926ff2831baf27bf0f5b383efa57676dcdb8c2aee4980d302  replicons.tsv
```
