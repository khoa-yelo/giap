# ADR-0030: Detection-limit model (M6-02): the pass-1 calibration choices
- Status: accepted-by-khoa (2026-10-21, in session: "3 c 4 a. The rest keep recommended")
- Date: 2026-10-21
- Spec sections: §7.4 (detection limit), §8.4; docs/specs/M6-02-detection-limit-model.md decisions 5–9
- Needs Khoa: yes (a protected tolerance key and event definitions). Answered: accept all for pass 1.

## Decision
1. **CITL gate.**
   - DL1 gates calibration-in-the-large at |CITL| ≤ 0.2. CITL is the intercept of d on logit(p̂), with the slope fixed at 1.
   - The gate uses the new protected key `detection_limit_citl_max_abs: 0.2` in tests/calibration/tolerances.yaml, added in a follow-up commit with `Gate-Change-Approved: ADR-0030`.
   - CITL is also recorded at k_min = 2.
2. **N1.** The one-detected event in the operon test means "fewer than k_min observable peptides seen" for the partner, not the gene-level FDR call.
3. **Replicate model in pass 1.**
   - With exactly 2 runs the replicate model is `usable = False` (`replicate_two_runs`).
   - The replicate `p_detect` stays the pooled per-run q.
   - For hardening: widening training to every gene seen in at least one run (which makes 2 runs usable), and 1 − (1 − q)^R.

## Evidence (before → after: metrics, tests, simulations)
The M6-02 statistician review and its re-check (68860bd, PR #157).

## Consequences
- `test_detection_limit_citl.py` moves from record-only to gated once the key lands.
- M6-05H uses the N1 definition.
