# R golden references: fgsea and goseq (M8-12; GIAP_SPEC.md §11.4, §11.5)

- Command: `sbatch scripts/slurm/golden_r.sbatch` (runs `scripts/golden/make_golden_inputs.py
  --seed 50000` and `Rscript --vanilla scripts/golden/make_golden_r.R` in env `giap-golden-r`).
- Slurm job: 1868639 on barley-01, 2026-10-22 (19:40 UTC).
- Env lock: `scripts/golden/golden_r.explicit.txt`, sha256 `eb4448ffcc114159f125021649dc3a870bf31a06a796c88117de1a5eeba10b7c`
  (R 4.5.3,
  fgsea 1.36.2, goseq 1.62.0, BiasedUrn 2.0.12; conda-forge and bioconda only).
- Scripts: git commit `88627e4056c88acdf872eced16acca8065fd80b4`.
- Seed: GOLDEN_SEED = 50000 (`tests/calibration/tolerances.yaml`); `set.seed(GOLDEN_SEED)`
  before `fgseaMultilevel` and before `goseq`; `nproc = 1`.
- Reproducibility: the inputs and the R outputs were each generated twice in this job and were
  byte-identical (M8-12 examples 2 and 3).
- Inputs are synthetic (random IDs and memberships): no KEGG, GO or COG content. Only numbers
  are committed; no R package source.
- Regeneration only through this job or `.github/workflows/golden-r.yml`, with a reason (a
  package upgrade needs an ADR); CI fails a change to `tests/data/golden/r/` that does not
  update this file.
