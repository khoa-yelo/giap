# ADR-0015: Two-block confirmation for calibration gates (M4-05)
- Status: accepted-by-khoa (2026-10-14, in session: "Two-block confirmation")
- Date: 2026-10-14
- Spec sections: §8.4 (calibration criteria), §12 (nightly), tests/calibration/tolerances.yaml (unchanged)
- Needs Khoa: yes (changes what a calibration gate means) — answered

## Context
Each §8.4 check run on a method exactly at level raises a false alarm about 2.7 % of the time at R = 1000,
and 8–9 % at α = 0.01 with R = 50–100. With ~46 checks per nightly run, P(at least one red check | every
method correct) ≈ 0.75 (statistician, PR #54, docs/methods/calibration.md).

## Decision
A check that fails on replicates 0..R−1 is re-run once on the fresh block R..2R−1 (same suite seeds,
`SeedSequence([calibration_base, crc32(suite), rep])`, so the second block is independent and reproducible).
The gate fails only if BOTH blocks fail; a first-block failure that the second block does not confirm is
logged as `unconfirmed` (visible in the nightly summary, not red). Applies to the nightly/full runs and to
the per-PR fast subset alike. Thresholds, R and k in tolerances.yaml are unchanged.

## Evidence (statistician)
False alarms per nightly run ≈ 0.75 → ≈ 0.04. Power vs 1.5× FDR inflation: at α = 0.05 0.92 → 0.85
(single check); realistic BH runs at full R = 300 stay at 1.0 for q ≥ 0.05.

## Consequences
The harness implements the confirmation; records carry `block` (1 or 2) and `status` (pass | unconfirmed |
fail). Worst case doubles the replicates of a failing suite.
