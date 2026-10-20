# Runtime-measurement fixtures (M2-21 example 7; `scripts/measure_runtimes.py`)

- `m2_21_trace_sim100k.tsv`: Nextflow `trace.tsv` copied unmodified (test-author, 2026-10-16) from
  `local/scratch/real-run-m2-20/rerun2/out/logs/trace.tsv` (sha256
  `d1fc81cdfb80acc14c2caee992511f021b13490a650d7bf524c5229f653f530f`). Real `giap run` of the
  M2-20 branch (`m2/context-layers` at `9cabfc8`) on the giapsim SIM100K genome
  (`review-pr78/sim/genome.fasta` + `proteins.faa`) with the real `$GIAP_DB_DIR` (Bakta
  `db: auto`, full DB present), `--threads 16 --memory 60.GB --profile farmshare_local,apptainer`,
  Slurm job `1851578` (`giap-real-run`, 16 CPUs, 64 GB). 37 tasks, all `COMPLETED` (ANNOTATE
  only; MAPPING was not reached in that run).
- `m2_21_sacct_1851578.txt`: real output of
  `sacct -j 1851578 --format=JobID,JobName,NodeList,MaxRSS,Elapsed,State -P` (2026-10-16): node
  `barley-01`. The tests feed it to `measure_runtimes.py` through a fake `sacct` on `PATH` that
  answers any `--format` field list from these rows.
