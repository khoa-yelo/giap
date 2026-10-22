# fgsea 1.36.2 recorded outputs (M8-17a; GIAP_SPEC.md §11.4, ADR-0038)

Each directory holds one call of `workflow/bin/giap_fgsea.R`: the inputs the wrapper wrote
(`stats.tsv`, `jobs.tsv`, `sets.tsv`), the command (`command.txt`), the outputs (`results.tsv`,
`versions.txt`) and `PROVENANCE.md`.

| Directory | Origin |
|---|---|
| `by_hand/` | **recorded-by-hand** by the test author (2026-10-22), not R output. The numbers are invented; only the file format follows docs/specs/M8-17-r-fgsea-engine.md §3. It pins the parser (`giap.io.fgsea_r.read_results`, `read_versions`) and stays as a format fixture. |
| `sim_main/`, `sim_retest/`, `ties/`, `unbalanced/`, `tiny_p/` | to be recorded by `sbatch scripts/slurm/record_fgsea.sbatch` (implementer, M8-17a §6) in the `giap-fgsea` env; each gets a real `PROVENANCE.md` (command, lock sha256, date, host, Slurm job ID). |

The implementer replaces this table row for every real recording it adds. `by_hand/` must never
be presented as R output.
