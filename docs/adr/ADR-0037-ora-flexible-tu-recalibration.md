# ADR-0037: The ORA TU detection recalibration becomes a flexible cross-fitted ridge logistic (pass 2)
- Status: proposed
- Date: 2026-10-22
- Spec sections: §7.7.3, §8.4, §8.5; amends ADR-0035 decisions 1 and 3; ADR-0034 (flag pattern); ADR-0015 (blocks); docs/specs/M8-04-category-summary-ora.md (lead decisions 12 and 13); docs/methods/function.md "Detection enrichment" (origin/m8/ora)
- Needs Khoa: yes. It changes the model behind every ORA p-value, and it sets the rule for lifting `ora_uncalibrated`.

## Context
**ADR-0035 is in place, but its model is not calibrated.**
- ADR-0035 recalibrates the TU detection probability with a cross-fitted logistic regression on [1, logit q_raw, log n], where n is the TU's eligible-gene count.
- PR #179 implements it with a damped Newton fit (M8-04 decision 13).
- The conditional law is exact for whatever q it is given. Q itself is still misfit, so pass 1 puts `ora_uncalibrated` on every ORA row (decision 13; `ora.CALIBRATION_FLAGS`).

**Pass-1 failures.** The statistician re-checked #179 at c1e3ca0 with R = 2,000 per world, split into ADR-0015 blocks (replicates 0–999 and 1000–1999). Source: `$GIAP_ROOT/stat/m804b/recheck/pr179_recheck_comment.txt` and `summ_blocks.txt`, Slurm 1867941 and 1867942. Rates are at α 0.01 / 0.05 / 0.10, with the flag ignored. Every row below is above the §8.4 bound in both blocks.

| World | Cell | Block 1 | Block 2 | Bound per block |
|---|---|---|---|---|
| 4b c 2.5 (the repo arm) | large, under | 0.0140 / 0.0555 / 0.1082 | 0.0130 / 0.0587 / 0.1098 | 0.0120 / 0.0544 / 0.1060 |
| 4d (4b c 1.5) | large_hiq, under | 0.0299 / 0.1075 / 0.1956 | 0.0301 / 0.1118 / 0.1979 | 0.0129 / 0.0564 / 0.1089 |
| 4e (4b c 1.5) | big_large, under | 0.0235 / 0.1215 / 0.2195 | 0.0260 / 0.1175 / 0.2135 | 0.0144 / 0.0597 / 0.1134 |
| 4e (4b c 1.5) | big_multi, over | 0.0163 / 0.0723 / 0.1340 | 0.0107 / 0.0610 / 0.1293 | 0.0136 / 0.0580 / 0.1110 (fails at 0.05 and 0.10) |
| 4f deep coverage (eta0 +1.5) | large, under | 0.0184 / 0.0700 / 0.1278 | 0.0171 / 0.0684 / 0.1263 | 0.0120 / 0.0544 / 0.1060 |
| 4f deep coverage | loq_any, under | 0.0250 / 0.1000 / 0.1700 | 0.0224 / 0.0974 / 0.1792 | 0.0128 / 0.0562 / 0.1085 |
| 4f deep coverage | big_multi, over | 0.0170 / 0.0827 / 0.1633 | 0.0157 / 0.0877 / 0.1677 | 0.0136 / 0.0580 / 0.1110 |
| combo (4b + 4c: c 2.5, a +0.5, b 0.8) | large_hiq, under | 0.0807 / 0.2202 / 0.3219 | 0.0810 / 0.2165 / 0.3207 | 0.0129 / 0.0564 / 0.1089 |
| csize (c 3 for TUs of ≥ 4 genes, 0.5 otherwise) | large, under | 0.0743 / 0.2116 / 0.3181 | 0.0788 / 0.2144 / 0.3249 | 0.0120 / 0.0544 / 0.1060 |
| small genome (800 genes, 321–410 TUs) | hiq_any, under | 0.0140 / 0.0840 / 0.1475 | 0.0175 / 0.0795 / 0.1530 | 0.0144 / 0.0597 / 0.1134 (fails at 0.05 and 0.10) |

The full table, including more 4b c 2.5 cells and the combo large cell, is in `pr179_recheck_comment.txt`.

What holds at level in both blocks:
- the null world, 4c at a = ±0.5 and slope b = 1.4 (slope b = 0.7 fails one cell at α 0.10; see decision 7);
- random whole-TU and single-gene categories in every world;
- the oracle q (the exact TU marginal, on the same detections) in every cell of all 13 worlds.

So the failures come from the 3-parameter model, not from the conditional law. Any category selected on TU size or on q is at risk (methods note, "Pass-1 limit").

**The flexible comparator.** In the first review the statistician fitted a more flexible cross-fitted model on the same detections and categories. Source: `$GIAP_ROOT/stat/m804b/probe_flex.py`, `flex.sbatch` and `flex_*.jsonl`; Slurm 1867892 (4b c 1.5, 4b c 2.5, hidet, lowdet, small) and 1867900 (combo, csize); worktree head 3f82e0f; R = 3,000 per world. The statistician reported it "at level in every cell of all 8 worlds" (`pr179_comment.txt`). Seven `flex_*.jsonl` world files exist, listed below; the "8" is not reconciled (owed, see the questions at the end). Selected cells, pooled over R = 3,000, at α 0.01 / 0.05 / 0.10, with the pooled bound in parentheses:

| World | Cell | Flexible model |
|---|---|---|
| 4b c 2.5 | large, under | 0.0053 / 0.0303 / 0.0623 (0.0111 / 0.0525 / 0.1035) |
| 4b c 2.5 | large_hiq, under | 0.0029 / 0.0225 / 0.0520 (0.0117 / 0.0537 / 0.1051) |
| 4b c 1.5 (4d) | large_hiq, under | 0.0029 / 0.0207 / 0.0507 (0.0117 / 0.0537 / 0.1051) |
| 4b c 1.5 (4e) | big_large, under | 0.0008 / 0.0075 / 0.0285 (0.0126 / 0.0556 / 0.1077) |
| 4b c 1.5 (4e) | big_multi, over | 0.0046 / 0.0294 / 0.0697 (0.0121 / 0.0546 / 0.1063) |
| hidet c 2.5 (4f) | large, under | 0.0041 / 0.0239 / 0.0531 (0.0111 / 0.0525 / 0.1035) |
| hidet c 2.5 (4f) | loq_any, under | 0.0050 / 0.0253 / 0.0587 (0.0116 / 0.0536 / 0.1049) |
| hidet c 2.5 (4f) | big_multi, over | 0.0038 / 0.0258 / 0.0589 (0.0121 / 0.0546 / 0.1063) |
| combo | large_hiq, under | 0.0043 / 0.0236 / 0.0561 (0.0117 / 0.0537 / 0.1051) |
| csize | large, under | 0.0052 / 0.0298 / 0.0636 (0.0111 / 0.0525 / 0.1035) |
| small c 2.5 | hiq_any, under | 0.0017 / 0.0202 / 0.0523 (0.0126 / 0.0556 / 0.1077) |
| lowdet c 2.5 | big_multi, over | 0.0053 / 0.0343 / 0.0742 (0.0121 / 0.0546 / 0.1063) |

How these numbers were obtained and what they do not show:
- The 4b c 2.5 large-under row is the statistician's own figure. The lead recomputed the other rows by running `$GIAP_ROOT/stat/m804b/summ.py`, unchanged, on the seven `flex_*.jsonl` files. That pass found no flexible-model cell above its bound in any world, in either direction, at any α.
- Splitting the same records at replicate 1000 (`recheck/summ_blocks.py`, unchanged) also gives no exceedance. Those are 1,000 and 2,000 replicates, not two equal ADR-0015 blocks, and no statistician has reviewed that split.
- **The flexible model has no measured rates yet** in the null world, 4b c 3.5, 4c at a = ±0.5, or slope b = 0.7 and 1.4 (owed).
- **Its power has not been measured** against the 3-parameter model on any planted-effect arm (owed).
- Many cells are well below α, for example big_large under at 0.0008 against α 0.01. The model may be conservative there.

## Decision
Adopt the flexible model as the pass-2 replacement for ADR-0035 decision 1's model. The cross-fitting, the conditional law, the tilt and the numerics of ADR-0035 decision 2 do not change. Each detail below is marked **[probe]** when `probe_flex.py` fixes it, or **[proposed]** when it is a default this ADR proposes.

### 1. Design matrix (12 columns)
For TU t, with ℓ_t = logit q_raw,t = log q_raw,t − log_comp,t and n_t = the TU's eligible-gene count (`TuUniverse.n_genes`):

x_t = [ 1, ℓ_t, (ℓ_t − κ₁)₊, (ℓ_t − κ₂)₊, (ℓ_t − κ₃)₊, (ℓ_t − κ₄)₊, 1{n_t = 2}, 1{n_t = 3}, 1{n_t = 4}, 1{5 ≤ n_t ≤ 6}, 1{n_t ≥ 7}, ℓ_t · log n_t ]

- **[probe]** `design()`:
  - the intercept, the linear ℓ term, 4 hinges, 5 size indicators with n = 1 as the reference level, and the interaction ℓ · log n;
  - no separate log n main effect, because the indicators carry the size effect;
  - (a)₊ = max(a, 0).
- **[probe docstring vs code]** The module docstring says "5 quantile knots" and "{1,2,3,4,5-6,7+}". The code uses 4 knots and 5 indicators, with size 1 absorbed into the intercept. This ADR adopts what the code does, because that is what was measured.
- The fitted η_t = x_t·β̂ gives q_t = expit(η_t) and comp_t = expit(−η_t), each computed directly (ADR-0035 decision 2, unchanged).

### 2. Knots (cross-fitted)
- **[probe]** For each fold j, the knots κ₁…κ₄ are the 0.2, 0.4, 0.6 and 0.8 quantiles of ℓ over the training TUs only (fold ≠ j). They use numpy's default linear interpolation (`np.quantile`).
- **[probe]** The held-out fold's design uses the training fold's knots. A held-out TU's ℓ never affects its own knots.
- **[proposed]** Repeated knots, which happen when many training TUs share one q_raw, are de-duplicated. A knot at or above the largest training ℓ is dropped, because its hinge column would be all zero. The probe does neither and relies on the ridge to absorb the collinear or empty columns.

### 3. Size bins
- **[probe]** n_t is the eligible-gene count, not `tu_size`. The bins are {1} (reference), {2}, {3}, {4}, {5, 6}, {7, …}.
- **[proposed] No upper cap.** The probe codes the top bin as `(7, 99)`, so a TU with ≥ 100 eligible genes would land in the reference bin. That cap never mattered in the probe, because the generator caps TU sizes at 8 (`probe.py` `world`). Production uses n ≥ 7 with no cap.
- **[proposed] Thin bins are merged.** A bin whose training TUs number fewer than `MIN_BIN_TUS` = 10 is merged into the next smaller non-reference bin: 7+ into 5–6, then into 4, then into 3, then into 2. A non-reference bin with no smaller bin left is merged into the reference. This happens per fold, and held-out TUs follow the same merge. The probe leaves a thin bin in place, and the ridge shrinks its coefficient to near 0, so a held-out TU in that bin is predicted as a single-gene TU. That is the opposite of the size effect the bin exists to capture.

### 4. Ridge penalty
- **[probe]** The penalised log-likelihood is Σ_t [y_t η_t − log(1 + e^{η_t})] − ½ λ Σ_{k ≥ 1} β_k², with λ = 1e-4. The intercept is not penalised. The columns are not standardised.
- The penalty shrinks the ℓ slope toward 0, not toward the identity slope of 1. At λ = 1e-4 against roughly a thousand training TUs this is negligible, and it is what was measured.
- **[proposed]** λ is a module constant (`RECAL_RIDGE`), not a user parameter. A different value needs a new calibration run.

### 5. Solver
- **[probe]** The probe uses plain Newton on the penalised objective: the gradient is Xᵀ(y − μ) − λβ and the Hessian is XᵀWX + λI′ (I′ has 0 for the intercept). It starts at β = 0, runs at most 100 iterations and stops when max |step| < 1e-9. It has no step halving and no failure return, so a non-converged fit is used silently.
- **[proposed]** Production reuses #179's damped Newton (`ora._logistic_mle`, decision 13) with the penalty added to the objective, the gradient and the Hessian:
  - it starts at the identity: intercept 0, ℓ slope 1, every other coefficient 0;
  - each step is halved until the penalised log-likelihood does not decrease beyond `LL_NOISE`;
  - it has converged when the full step is below `NEWTON_TOL` = 1e-10, within `MAX_NEWTON` = 50 iterations and `MAX_HALVINGS` = 60 halvings;
  - it returns None on a non-finite step or no convergence.

  With λ > 0 the penalised objective is strictly concave once the intercept column is non-zero, so the Hessian is non-singular and the stationary point is unique. The #179 re-check's argument that the slack cannot fake convergence carries over (`recheck/newton_summary.txt`).
- **[proposed]** An independent per-fold oracle test (scipy `minimize` on the same penalised objective, about 1e-8 relative) replaces `test_ora_recal_oracle.py`'s unpenalised oracle.

### 6. Skip rules and the minimum TU count
The model has p = 12 coefficients, against 3 for the pass-1 model. With 5 folds, each fit sees about 4/5 of the TUs. The small-genome world (321–410 TUs, `recheck/res_small_c2.5.jsonl`) therefore fits 12 coefficients on about 260–330 training TUs, and the flexible model was at level there. How many detected and undetected training TUs it had is not logged (owed).

The proposed ladder is evaluated per run. Every fold must meet a rung for the run to use it.

1. **Flexible model.** Used when every fold's training TUs have min(detected, undetected) ≥ `EPV` × (p_j − 1). Here p_j is the column count after the knot and bin reductions of 2 and 3. **[proposed]** `EPV` = 10 events per non-intercept column, the usual logistic rule of thumb, which is 110 per class at full size. Whether the 800-gene world meets 110 per class is unmeasured. If it does not, the statistician should re-run the small world at a lower `EPV` before this rung is fixed (owed).
2. **Reduced model.** Used when rung 1 fails but every fold still has min(detected, undetected) ≥ `MIN_CLASS` = 10, the pass-1 rule. The run uses the ADR-0035 3-parameter model. Every row keeps `ora_uncalibrated` permanently, and the run adds **[proposed]** the new token `tu_recal_reduced` (run flag and row flag; it is registered through the schemas, as in ADR-0036).
3. **Skip.** As in pass 1: `folds` is None or empty, a fold has fewer than `MIN_CLASS` of either class, or a fit fails. q = q_raw, and every row and the run carry `tu_recal_skipped`.
   - **[proposed]** A failed flexible fit does not fall back to the reduced model. It skips. A silent change of model is what B3/B4 warned against.
   - **[proposed]** Complete separation of the training classes by the fitted training η is also a skip. Ridge keeps β finite, but the fitted q would be meaningless.

**Skipped and reduced fits keep `ora_uncalibrated` permanently.** It does not depend on the lift below. As the #179 re-check carry-over requires, the commit that empties `CALIBRATION_FLAGS` must add the flag explicitly on these paths. `test_ora_review_179.py` pins this.

### 7. Lift rule for `ora_uncalibrated` on fitted rows
The flag comes off rows fitted with the flexible model only when all of the following hold in one recorded full-tier run:
- every cell of arms 4b, 4c, 4d, 4e and 4f, of the combination world (4b + 4c), of size-dependent co-detection (csize) and of the small genome passes the §8.4 bound **in both blocks**, at every α in {0.01, 0.05, 0.10} and in both directions. Cells include the low-q (`loq_any`) and ≥ 2-gene (`big_multi`) categories, not only large and high-q ones (re-check "Pass-2 target");
- each block has R ≥ 1,000;
- the rates are logged per cell and per block, with the commit and Slurm job ID, in the calibration records. Full-tier rate logging is owed (backlog).

This is stricter than an ADR-0015 gate, where a failure only in the first block counts as `unconfirmed`, not as a fail. Here an `unconfirmed` cell also blocks the lift.

**[proposed]** Two more conditions:
- The slope b = 0.7 world (gene-level p̂ too steep) joins the cells that must pass. In `recheck/summ_blocks.txt` its big_large under cell fails at α 0.10 in both blocks: 0.1185 and 0.1175 against 0.1134 (Slurm 1867941 and 1867942). The re-check comment's table does not list it.
- The cells at level today must stay at level under the same rule: the null, slope b = 1.4, 4c, and the random whole-TU and single-gene categories.

## Options
- **(a) The flexible model (proposed above).** It is the only candidate with measured rates at level in the cells where pass 1 fails (see Context).
  - Cost: 12 coefficients per fold instead of 3, and still one Newton fit per fold per run.
  - Risks:
    - the small-genome and rule-count limits (decision 6);
    - its power and its null-world rates are not yet measured;
    - its form was chosen by looking at the same attack worlds it was scored on, so the full-tier run has to use fresh seeds.
- **(b) Keep the 3-parameter model and the flag permanently.** Rejected. ORA would never become a discovery claim. Its failures are in routine category kinds within [5, 500], such as high-q operons, categories of 100–200 multi-gene TUs, and deep-coverage runs, and they reach 0.07–0.08 at α 0.01 (combo large_hiq under, 0.0807 / 0.0810; csize large under, 0.0743 / 0.0788). This is ADR-0035's rejected "flag and ship" option, made permanent.
- **(c) A nonparametric recalibration: isotonic regression, or a GAM.** Not chosen.
  - Isotonic regression in ℓ alone cannot express the size effect that csize and combo need. A version per size bin multiplies the small-sample problem of decision 6, and its step function puts ties at the bin edges.
  - A GAM, for example a spline in ℓ by size with a penalty chosen by cross-validation, is close in spirit to (a). It adds nested tuning inside each cross-fit fold and possibly a new dependency (licence check under the import rules: owed).
  - Neither has measured rates (owed). (a) is a fixed-basis, lightly penalised version of the same idea, with evidence.
- **(d) A stratified exact TU marginal: model the oracle q directly.** Not chosen for pass 2.
  - The statistician notes that the oracle q, 1 − E_u ∏_g (1 − expit(η_g + c·u)), is at level in every cell of every world.
  - **It is not computable from GIAP's inputs as they stand.** It needs the per-gene conditional linear predictor η_g and the strength and distribution of the shared TU effect, c and u ~ N(0, 1). M4-13 provides only the marginal p̂, and the probe computes the oracle from the generator's true η and c with 40-node Gauss–Hermite quadrature (`probe.py` `world`).
  - Making it computable needs three things:
    - a TU random-intercept detectability model in M4-13, which is a change to another module and needs its own ADR;
    - an estimate of c, per size stratum, because csize shows c varies with size;
    - a per-gene inversion of p̂ = E_u expit(η_g + ĉ u) for η_g.
  - Its level would then depend on the estimation error in ĉ and on the logit-normal shared-effect form. That form is exactly the generator's, so the attack worlds cannot test it fairly.
  - It is the principled long-term target and a backlog candidate, not a pass-2 fix.

## Evidence (before → after: metrics, tests, simulations)
- **Before (pass 1, 3-parameter model):** the Context table. Source: `recheck/pr179_recheck_comment.txt` and `summ_blocks.txt`, Slurm 1867941 and 1867942, R = 2,000 in two blocks, commit c1e3ca0.
- **Comparator (flexible model, prototype):**
  - the Context table: `probe_flex.py`, `flex_*.jsonl`, Slurm 1867892 and 1867900, R = 3,000 pooled, worktree 3f82e0f, recomputed with `summ.py`;
  - the statistician's figure for 4b c 2.5 large under is 0.0053 / 0.0303 / 0.0623 (`pr179_comment.txt`);
  - the combo and hidet "recal" columns in those files used the old undamped Newton, which skipped 2,898 and 417 replicates. The flex columns were computed on every replicate whatever the PR model did.
- **After:** owed, all of it.
  - The full-tier run of the decision 7 cells with the production implementation (decisions 1–6), in two blocks of R ≥ 1,000, with logged rates.
  - The worlds the flexible model has not been run on (null, 4b c 3.5, 4c ±0.5, slopes).
  - Power against the 3-parameter model on a planted-effect arm.
  - The small-genome class counts under the `EPV` rung.

## Consequences
- **Code:** `ora._recalibrate` gains the design of decisions 1–3, the penalised damped Newton of decision 5 and the ladder of decision 6. `TuUniverse.flags` may hold `tu_recal_reduced`. The `enrichment_detection` schema's flag pattern gains `tu_recal_reduced`. The ORA columns do not change.
- **Tests:**
  - a per-fold penalised oracle;
  - knot and bin edge cases: tied ℓ, an empty top bin, a TU with ≥ 100 genes;
  - the rungs of the ladder;
  - the explicit `ora_uncalibrated` on the skip and reduced paths.
- **Calibration:** M8-15's attack file gains every decision 7 cell at the full tier, with rate logging. The fast tier stays a smoke test, because it cannot resolve 4b or 4f (re-check M5).
- **Methods note:** "Recalibration" and "Pass-1 limit" are rewritten, and the note states that skipped and reduced runs remain uncalibrated.
- **Until the lift:** nothing changes for users. Every row keeps `ora_uncalibrated`.
- **Questions for the lead and Khoa:**
  1. The statistician wrote "8 worlds", but there are 7 `flex_*.jsonl` world files. Which is the eighth? Possibly deep coverage counted twice, with and without the PR-skipped replicates.
  2. Is `EPV` = 10 acceptable, given the risk that real small genomes fall to the reduced model?
  3. Should the reduced rung exist, or should anything below rung 1 skip?
  4. Should the decision 7 lift rule really be stricter than ADR-0015, so that an `unconfirmed` cell also blocks the lift?
