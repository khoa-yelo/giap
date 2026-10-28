# ADR-0063: The determinants null on worker processes and partial dependence on threads; one new benchmark baseline
- Status: proposed
- Date: 2026-10-26
- Spec sections: §7.8.2 (step 4), §6.6, §8.2, §12, §14.4 (speed loop), R14
- Needs Khoa: no (no default, threshold, interface, dependency or protected-file change; outputs
  byte-identical; the step's CPU request is unchanged)

## Context
ADR-0053 left the elastic-net null serial ("a process pool for it is a possible follow-up; it
needs an ADR, because of its memory and CPU budget under R14"). On the real runs the
determinants step sits on the critical path after TerminatorNet: 10–12 min at 14–16 CPUs (pilot
pass 1–3, MVP-04), 829–890 % CPU of 14–16. Replaying the pilot pass-3 A1 step under cProfile
(rye-01, 8 CPUs, job 1897002, 1,440 s) put 712 s in the serial elastic-net null column and
458 s in the LightGBM null on its thread pool; a second profile after the first change (job
1897273) showed `pd_rows` (partial dependence, 25 boosters × 84 grid values, serial LightGBM
predicts) at 143 s.

## Decision
- `null_p`: with `threads` > 1 and B ≥ `_PROCESS_MIN_B` = 200 the refits of both models run on
  one pool of `threads` spawned worker processes (`spawn`, never `fork`; each worker
  `threadpool_limits(1)`), model by model, `pool.map` in b order (chunks of 8). Below B = 200
  (worker spawn ≈ 2 s outweighs the refits) the ADR-0053 path is kept: LightGBM on a thread pool,
  the elastic net serially. `threads` = 1 is serial as before.
- `pd_rows(..., threads=1)`: with `threads` > 1 the boosters predict on a thread pool
  (LightGBM's predict releases the GIL; each predict keeps the booster's single thread),
  collected in booster order. The step passes its `threads`.
- New benchmark `benchmarks/test_determinants_null_processes_bench.py` and one baseline entry,
  `determinants_null_b400_n400_t8`.

## Evidence (before → after: metrics, tests, simulations)
- Pilot pass-3 A1 determinants step replayed from its staged inputs (`gippa _step determinants`,
  `threads` 12, rye-01 AMD EPYC 9754, no profiler); origin/main f1ff5b2f run at the same time:
  | variant | origin/main | this branch | speed-up |
  |---|---|---|---|
  | elastic-net null on processes (285c9b05) | 826.6 s | 603.2 s | 1.37× |
  | both models on one process pool (4510dc94) | 844.1 s | 559.4 s | 1.51× |
  | + PD on threads (1ae0b2c6) | 900.3 s | 414.4 s | 2.17× |
  | final, with `_PROCESS_MIN_B` (aa4abe5a) | 826–900 s | 403.0 s | 2.0–2.2× |
  CPU time 1:52–1:59 h → 1:02 h; MaxRSS 0.51 → 2.15 GB (12 workers; the request is 48 GB).
- Byte identity: all six outputs (`determinants_features/model/attribution.tsv`,
  `determinants_pd.tsv`, `determinants_shap.tsv.gz`, `null_determinants.tsv`) have the same
  sha256 on origin/main and every variant, and the three tables equal the pilot run's own
  outputs (job 1896546, barley-03, 14 threads, pinned commit 4ccce8a5).
- `null_p` on the benchmark world (B 24, 4 threads): pickled null arrays, p and b `cmp`-equal to
  origin/main. New tests: `tests/unit/determinants/test_null_processes.py` (processes vs
  serial at 2 and 3 workers, both models, the small-B path, the pool switch) and
  `tests/unit/determinants/test_pd_rows_threads.py` (2 and 4 threads, log-scale and empty
  features). The determinants unit, property and step tests pass (411).
- Benchmarks (rye-01): `determinants_null_b50_n400` (4 threads, below `_PROCESS_MIN_B`)
  12.15 → 11.87 s (unchanged path); `determinants_null_b400_n400_t8` (8 threads, worker processes) 117.5 s (median of 3, 104.7–120.2; job 1897338) → 45.9 s (45.7–46.4; job 1897339), 2.56×; budget 140 s.

## Consequences
- On a 5-Mb genome the critical path loses about 5 min (10–12 → ~5 min at 14 threads).
- R14: the workers use the step's own CPUs (`threads` = the task's CPUs, unchanged), so no extra
  CPUs; memory ≈ 150 MB per worker.
- Proposed baseline addition: `determinants_null_b400_n400_t8` (45.9 s on rye-01; the nightly
  wheat class has no entry yet, so the nightly gate records it before it can compare). No existing entry changes.
