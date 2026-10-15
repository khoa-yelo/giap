# antismash 8.0.4 on sim100k (M2-10, GIAP_SPEC.md §11.4)

- Command: `antismash --genefinding-tool none --minimal --cpus 17 --databases /home/classes/bios/270/khoa/giap/local/db/antismash --output-dir antismash bakta.gbff`
- Image: `sha256:d3f390e687cfe02f1387be40256dddbf678e91bc56b96c36da07b56ec8bda29d` (versions.txt: antiSMASH 8.0.4)
- Database: antiSMASH 8.0.4 databases (`$GIAP_DB_DIR/antismash`)
- Date: 2026-10-15T07:19:44Z (run finished, UTC)
- Host: barley-01 (FarmShare; SLURM_JOB_ID not recorded; 17 CPUs; wall time 6 s)
- Recording: `/home/classes/bios/270/khoa/giap/local/recordings/antismash/sim100k`
- Inputs:
  - `/home/classes/bios/270/khoa/giap/local/recordings/antismash/sim100k/bakta.gbff`: sha256 `c028b1c43b22260713a77796c6ce79d855773b682a68c2f4a2fe36124481c6fd`
- `bakta.json` is a trimmed excerpt of `antismash/bakta.json` (sha256 `317521dd8c7e6f4df7c048d83cc4aceca6592d87b253e5e2fdd7f99424d6e17a`); `replicons.tsv` was written by the test author (M2-18); trimming rules in ../../PROVENANCE.md

## Files (sha256sum)

```
75ea69289edadf9d77bd7f4d90a12f890bf488a7f2e861c6b4cc8ca99ea32df2  bakta.json
ccb4832800acea547b8ac2d045719beaba8c59a40e003215132ac8856bb8817a  replicons.tsv
```
