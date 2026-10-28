# ADR-0053: Bit-identical speed-ups of the determinants step; two new benchmark baselines
- Status: proposed
- Date: 2026-10-23
- Spec sections: §7.8.2 (steps 1, 3, 4), §6.6, §8.2, §12
- Needs Khoa: no (no default, threshold, interface, dependency or protected file changes; outputs byte-identical)

## Context
The per-commit fast suite took 1,956 s on 16 CPUs (Slurm job 1879345, barley-01, no coverage). One
test, `tests/calibration/test_gene_fdr_real_mapping.py::test_gene_fdr_one_replicate_through_real_mapping`,
took 1,532 s and was the suite's critical path. Profiling (cProfile; py-spy and scalene are not
installed in the shared venv) showed that `gippa run --executor local` spent most of its time in
the `determinants` step: on the 500-kb replicate, `cross_validate` 682 s, `null_p` (1,000
permutations) about 1,000 s, `univariate` 44 s. On sim50, determinants took 60 s of an 84 s run.

Causes:
1. `ElasticNetCV` calls `enet_path` per inner fold with `check_input=True`, which re-validates the
   Gram matrix once per alpha (`check_array`); that was about 65 % of every elastic-net fit.
2. The elastic net holds the GIL: threading its permutation refits was slower than one thread
   (0.21 s per permutation at 1 thread, 0.35 s at 4).
3. `cross_validate` fitted LightGBM with `n_jobs = threads` (16) on a few hundred rows: one CV
   repeat took 204 s at 16 OpenMP threads and 0.3 s at 1 on a shared node.
4. `tu_block_bootstrap` built a data frame per replicate, and `scipy.stats.spearmanr`/`rankdata`
   have a large per-call overhead (about 90,000 calls per run).

## Decision
- `_ElasticNetCV` (a subclass) routes its inner-fold path through `enet_path(check_input=False)`.
  `_path_residuals` has already validated X, y, Xy and the Gram matrix, so the skipped checks
  returned their inputs unchanged.
- `null_p` runs model by model: LightGBM permutation refits on `threads` workers, then the elastic
  net's serially. They are not overlapped, because that causes a GIL convoy on a loaded node.
- `_repeated_oof` fits every fold single-threaded. With `threads > 1` the LightGBM fits run on a
  pool of `threads` workers and are collected in (repeat, fold) order. The elastic net stays serial.
- `tu_block_bootstrap` gains a keyword-only `columns=` path. It calls `stat_fn(*arrays)` with the
  same draws and the same replicate rows. `univariate` uses it, with exact average ranks in place
  of `rankdata`/`spearmanr` (each rank is (i + j)/2 of integers, and `np.corrcoef` of the stacked
  ranks is what `spearmanr` computes).
- New benchmark `benchmarks/test_determinants_bench.py` and two baseline entries.

## Evidence (before → after: metrics, tests, simulations)
- Byte identity (before = origin/main d29a2ea5):
  - `null_p` (40 permutations) null R² arrays, p and b;
  - `cross_validate` r2/oof/ΔR² at threads 1, 4 and 16, and SHAP blocks;
  - `univariate` tables on the 500-kb replicate and sim50;
  - the three determinants tables of the 500-kb replicate (threads 16).
  `_ranks`/`_rho` equal `rankdata`/`spearmanr` bit for bit on 20,000 random inputs with ties and
  ±0. All determinants unit, property, integration and calibration tests pass. Among them,
  example 4 checks that threads 1 and 8 give byte-identical tables.
- Determinants step on the 500-kb replicate (iron-13, threads 16): `cross_validate` 682 → 9.7 s,
  `null_p` about 1,000 → 305 s, `univariate` 44 → 8 s; total about 1,730 → 323 s.
- Benchmarks (iron-13): `determinants_univariate_b1000_n400` 7.29 → 1.57 s;
  `determinants_null_b50_n400` (4 threads) 33.74 → 11.09 s.
- Fast suite (16 CPUs, barley, `-p no:randomly`, no coverage): 1,956 s (job 1879345) → 1,437 s
  with the first commit (job 1879720). The final run is listed in the commit and PR text.

## Consequences
- `CvResult.boosters` now hold LightGBM models fitted with `n_jobs=1`. Their predictions are
  identical, and SHAP is computed single-threaded per booster.
- The elastic-net null is serial: about 0.2 s per permutation, which is the step's floor without
  processes. A process pool for it is a possible follow-up; it needs an ADR, because of its
  memory and CPU budget under R14.
- Proposed baseline additions: `determinants_univariate_b1000_n400`, `determinants_null_b50_n400`
  in `benchmarks/baseline.json`. No existing entry changes.
