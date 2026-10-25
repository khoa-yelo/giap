# Runtimes

This page will hold the measured wall-clock times of a standard `giap run`, from
`scripts/slurm/benchmark_standard.sbatch` (M10-12). Each row will give the median of three runs,
with the node class, the date and the GIAP commit. The numbers will come from
`benchmarks/baseline.json › standard_run` and the M10-12 decision record.

**Status: pending (M10-12).** No standard-run benchmark has been committed yet. The table below is
empty on purpose: GIAP's documentation states no runtime that has not been measured.

| Run | Median (h:mm) | Range of 3 runs | Peak memory (GiB) | Threads / memory | Node class | Date | Commit |
|---|---|---|---|---|---|---|---|
| E. coli K-12, all modules | pending (M10-12) | pending (M10-12) | pending (M10-12) | 16 / 64 GB | pending (M10-12) | — | — |
| E. coli K-12, `--modules core` | pending (M10-12) | pending (M10-12) | pending (M10-12) | 16 / 64 GB | pending (M10-12) | — | — |

What affects the runtime, qualitatively:

- `--publication` raises permutations and bootstrap draws from 1,000 to 10,000;
- `--promoters genome` predicts promoters genome-wide instead of in upstream windows;
- Bakta runs unless `--annotation bakta.json` is given, and the full database is slower than the
  light one;
- the optional modules: antiSMASH (`bgc`) and geNOMAD (`mobile`).

Nextflow writes per-task times and peak memory to `OUT/logs/trace.tsv` and
`OUT/logs/nf_report.html`. The manifest records per-step time and peak memory.
