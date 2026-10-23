# fgsea 1.36.2 recorded outputs (M8-17a; GIAP_SPEC.md §11.4, ADR-0038)

Each directory holds one call of `workflow/bin/giap_fgsea.R`: the inputs the wrapper wrote
(`stats.tsv`, `jobs.tsv`, `sets.tsv`), the command (`command.txt`), the outputs (`results.tsv`,
`versions.txt`) and `PROVENANCE.md`.

| Directory | Origin |
|---|---|
| `by_hand/` | **recorded-by-hand** by the test author (2026-10-22), not R output. The numbers are invented; only the file format follows docs/specs/M8-17-r-fgsea-engine.md §3. It pins the parser (`giap.io.fgsea_r.read_results`, `read_versions`) and stays as a format fixture. |
| `sim_main/` | real R (Slurm job 1868708): the main call of the M8-06 step fixture (`World`, seed 7, `n_perm_simple` 200), collections cog, kegg, go_bp, go_slim |
| `sim_retest/` | real R (1868708): that ranking's GO re-test call; no fixture GO term reaches q <= alpha, so significance is forced (each tested go_bp term with tested descendants minus its lowest-p tested descendant; see its PROVENANCE.md) |
| `ties/` | real R (1868708): the T5 ranking, 20 of 50 scores tied |
| `unbalanced/` | real R (1868708): set `BOTTOM` has `pval`, `NES` and `log2err` NA (fgsea's `modeFraction < 10`) |
| `tiny_p/` | real R (1868708): set `PACKED` has p = 3.06e-50 with a finite `log2err` |
| `log2err_na/` | real R (1868708), extra: set `BOTTOM` has a finite p = 3.84e-24 with `log2err` NA (`isCpGeHalf`, "likely overestimated"; M8-17 example 4) |

All real recordings come from `sbatch scripts/slurm/record_fgsea.sbatch <worktree>`
(`scripts/dev/record_fgsea.py`); each `PROVENANCE.md` gives the command, the lock sha256, the
date, the host and the job ID.

`by_hand/` must never be presented as R output.
