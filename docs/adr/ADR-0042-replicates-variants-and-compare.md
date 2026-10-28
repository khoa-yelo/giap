# ADR-0042: Calibration replicate counts for variants and compare
- Status: accepted-by-khoa
- Date: 2026-10-22
- Spec sections: §8.4, §12, §7.5.4, §7.9.3; M7-10, M7-13, M9-22
- Needs Khoa: yes (protected `tests/calibration/tolerances.yaml`; decided by Khoa on 2026-10-22)

## Context
Two calibration claims have no `replicates_fast` or `replicates_full` key:
- the variant suites: the M7-10 class error rates and the M7-13 F1 gate (`f1_min.variants_high` = 0.90 already exists);
- the TU-block bootstrap CI for Spearman ρ in `gippa compare` (M9-22).

Without a count neither suite can gate.

## Decision
Add to `tests/calibration/tolerances.yaml`:
- `variants: 20` (fast) and `variants: 100` (full).
  - Each replicate is a full simulation, alignment, pileup and calling, the same cost class as `detectability` and `determinants` (20/100).
  - At F1 ≈ 0.93 with a replicate spread of ≈ 0.03, the fast standard error is ≈ 0.007, enough to catch a drop below 0.90 at `se_multiplier` 2.
- `compare: 100` (fast) and `compare: 1000` (full).
  - Replicates are cheap: table arithmetic on simulated pairs, like `ora` (100/1000).
  - At 1,000 replicates, 95 % coverage is measured to about ±1.4 %.

## Consequences
- The M7-10 and M7-13 suites, and `tests/calibration/test_compare_rho_ci.py` (M9-22), gate instead of recording.
- The values are applied in a follow-up PR carrying `Gate-Change-Approved: ADR-0042`.
