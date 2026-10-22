# ADR-0032: Univariate determinants CIs stay `ci_unestimated` in pass 1; a calibrated SE in hardening
- Status: accepted-by-khoa (2026-10-21, in session: "Take your recommendations")
- Date: 2026-10-21
- Spec sections: §7.8, §8.2, §8.5; docs/specs/M8-08-determinants-models.md decision 7
- Needs Khoa: yes (an uncalibrated interval). Answered: option (a) now, (b) in hardening.

## Context
Under a TU-structured null, the TU-block bootstrap 95 % CIs of Spearman ρ and AUROC miss at 0.065–0.072. The bound is 0.064 (statistician, Slurm 1866223 and 1866311). The bootstrap SE is 5–9 % low, and percentile, normal and basic intervals all miss alike.

## Decision
- **Pass 1:** every univariate row carries `ci_unestimated`. The numbers are descriptive screening statistics with no discovery claim. The attack test records the level and does not gate it while the rows are flagged.
- **Hardening:** a studentised TU-level bootstrap or a cluster-jackknife SE. Once its attack passes, the flag is removed.
- Rejected: TU-collapsed statistics, because they change what the number means.

## Evidence (before → after: metrics, tests, simulations)
- tests/calibration/attacks/test_attack_determinants_tu_null.py, fast tier (Slurm 1866499): ρ 0.09 and AUROC 0.07 at R = 100.
- The gene bootstrap misses at 0.25–0.27.

## Consequences
The report and V15 show these CIs as uncalibrated.
