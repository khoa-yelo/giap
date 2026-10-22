# M6-01 lead decisions (2026-10-21, after the test-author's report)

Recorded here because the M6 design PRs (#134–#136) were still merging; fold into M6-01's own spec when
they land.

1. **Order of operations: maximise first, then monotonise** (spec step 5 and Q3), not the acceptance
   list's "element-wise maximum of the per-decoy monotone curves". The two differ
   (monotone-of-max ≥ max-of-monotones). The acceptance text is **wrong and must be corrected**; the
   test additionally asserts `q ≥ max of the monotone curves`, so the order is pinned either way.
2. **"Monotone" means `giap.stats.tda.monotone_q`** — the running minimum over *looser* thresholds, the
   only direction that keeps q an upper bound and the one the merged engine already uses. The spec's
   "running minimum from the top" reads as a constant equal to the strictest FDR; fix the wording.
3. **`DecoyFdr.lfdr`**: the same "larger wins" rule as the FDR — the element-wise maximum over decoy
   sets, with the per-bin weights `N_NO,b / N_decoy,b` from step 4.
4. **Tie-break in `which`** at the loosest threshold is forced (FDR = π₀ for every decoy set there, the
   +1 terms cancel): **the first key of `decoys` wins**, documented in the docstring.
5. **Skipping an empty decoy set lives in `max_decoy_fdr` only.** `reweighted_fdr` has no skipping rule,
   so example 7's 10/3 stands for it. Example 7 as written contradicts the degenerate row; this is the
   reading.
6. **Calibration R comes from `tolerances.yaml`** — `replicates_fast/full.missed_operons` (50/300) and
   `.discordance` (100/1,000). The acceptance text's "R = 200" invents a number not in that protected
   file; ignore it.
7. **Bandwidth: scipy's `gaussian_kde(bw_method="silverman")`**, named explicitly in the methods note,
   because the lfdr values depend on it and the rule-of-thumb variant differs.
8. **Expose the densities for testing.** Reflection at 0 cancels in the lfdr ratio, so no test through the
   public return value can catch a wrong reflection. The implementation exposes the `f0` / `f_NO`
   densities (a private helper is fine) so a test can pin the reflection directly. A silent error here
   would corrupt every lfdr.
9. **Power warning — for the statistician at M6-09, and recorded now.** With π₀ = 1 and `N_NO,b` counting
   all targets, the called set at q ≤ 0.1 is essentially empty unless the planted |Δ| sits deep inside
   the null (the test-author needed planted sd 0.02 against null sd 1.0, with a 40 % planted share, to
   make the FDP check non-vacuous at 1,200 targets). **At the M6-09 grid's realistic effect sizes
   (0.25–1.0 SD) this family may have almost no power.** That must be measured and stated before the
   M6-09 gate is claimed; if it holds, the missed-operon call is reportable but near-powerless and the
   methods note has to say so plainly rather than implying a working test.
10. **Attack A (distance confounding) is gated on the reweighted FDP only**; the single-bin number is
    *recorded* (`local/logs/calibration/missed_operons_reweighting.jsonl`) with only
    `single_bin >= reweighted` asserted. A convincing version needs M6-03's generators and belongs to
    M6-09 — the inflation is capped by the density ratio, not the bin mix, so the hand-built attack
    produced a false call in only ~14 % of replicates.
11. The performance budget (50,000 targets × 2 decoy sets × 5 bins < 2 s) goes in
    `benchmarks/test_decoy_fdr_bench.py`, added by the implementer — house convention (M5-02 decision 12)
    keeps performance targets out of `tests/`.

## Lead decisions after the statistician's review of #140 (2026-10-21)

12. **Power (Khoa, 2026-10-21: "Label it now, ADR later").** The statistician showed the floor on the
    attainable FDR is set by the score, not by the +1 or the bins: in the large-sample limit it is
    π0 / (π0 + (1 − π0)·L) with L = σ0/σ1, so q ≤ 0.1 needs L ≳ 23–46. Over the M6-09 grid
    (0.25–1.0 SD, L ≤ 4) no variant makes a call. Merge as specified; `docs/methods/operons.md` states
    plainly that the missed-operon family is near-powerless at realistic effects, controls FDR, and its
    lfdr is a ranking output. M6-09 gates FDR control and **claims no power**. A combined score
    (signed-Δ likelihood ratio + UniOP p + distance + promoter/terminator evidence, same decoy FDR, decoys
    matched on every covariate) is a **hardening-pass ADR draft for Khoa**, not M6 pass 1.
13. A degenerate (collapsed-scale) t fit is unestimated (`converged = False`, NaN parameters), fixed in
    7ad17fe; the KDE accuracy stated in code and PR is ~3e-3, as measured.
