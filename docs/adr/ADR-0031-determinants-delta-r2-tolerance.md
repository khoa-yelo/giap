# ADR-0031: Determinants ΔR² bias tolerance 0.05; the ΔR² interval stays `ci_unestimated`
- Status: accepted-by-khoa (2026-10-21, in session: "3 c 4 a. The rest keep recommended")
- Date: 2026-10-21
- Spec sections: §7.8, §8.2; docs/specs/M8-08-determinants-models.md decisions 5–6
- Needs Khoa: yes (a protected tolerance key). Answered: as recommended.

## Decision
- Add the protected key `determinants_delta_r2_max_abs_bias: 0.05` to tests/calibration/tolerances.yaml, in a follow-up commit with `Gate-Change-Approved: ADR-0031`.
- On simulated data, a planted block's mean ΔR² must be within 0.05 of truth, and every noise block's within 0.05 of 0.
- The ΔR² t-interval over CV repeats is not calibrated, because the repeats share data. ΔR² rows keep `ci_unestimated`.
- A calibrated ΔR² CI is hardening work.

## Evidence (before → after: metrics, tests, simulations)
tests/calibration/test_determinants_models_calibration.py (`calibration_nightly`) gates on the key.

## Consequences
The M8-08 calibration gate becomes active once the key lands.
