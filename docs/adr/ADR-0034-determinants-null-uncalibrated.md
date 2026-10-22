# ADR-0034: The determinants permutation null p is `null_uncalibrated` in pass 1; a structure-aware null in hardening
- Status: accepted-by-khoa (2026-10-22, in session: "Flag it (Recommended)")
- Date: 2026-10-22
- Spec sections: §7.8, §8.4, §8.5; docs/specs/M8-09 (rule 6); statistician review of #173
- Needs Khoa: yes, because a p-value ships uncalibrated. Answered: flag it now, a structure-aware null in hardening.

## Context
The M8-09 null permutes TU response blocks within TU-size groups. It is valid when TUs are exchangeable, and in that world it holds level at R = 200, α 0.01 / 0.05 / 0.10:
- LightGBM: 0.010 / 0.040 / 0.105.
- Elastic net: conservative, 0 / 0.005 / 0.020.

Real proteomes also have inter-TU autocorrelation along the genome and detection selection. Under the statistician's §8.5 attacks (R = 100, B = 99; Slurm 1867304, 1867334, 1867348, 1867349) the null fails. Rejection rates against the bounds 0.030 / 0.094 / 0.160:

| Attack | LightGBM | Elastic net |
|---|---|---|
| inter-TU AR(1) abundance with an AR(1) covariate | 0.02 / 0.12 / 0.23 | 0.05 / 0.12 / 0.13 |
| the same plus `replichore_position` | 0.80 / 0.95 / 0.96 | 0.33 / 0.41 / 0.43 |
| detection depends on abundance and length, y independent of X | 0.21 / 0.53 / 0.69 | 0.64 / 0.83 / 0.86 |

The statistician tried three alternatives, and none fixes every arm:
- Circular shifts within size groups: not enough.
- Spatially contiguous TU folds: fix the AR(1) arm (elastic net 0.02 / 0.06 / 0.06), but not the arm with `replichore_position`.
- Stratifying by detectability: still 0.46 at α 0.05 in the detection arm.

## Decision
- **Pass 1:** every null result carries `null_uncalibrated`, and M8-10 carries it onto the `null_p` row. The p is valid only for exchangeable TUs and is not a discovery claim. A NaN observed R² also gets `unestimated`. The attack test (tests/calibration/attacks/test_attack_determinants_null.py) records the level, and does not gate on it, while every result is flagged (the ADR-0032 pattern).
- **Hardening:**
  - spatially contiguous CV folds for both the observed CV and the null;
  - a decision on positional features (`replichore_position`): drop them, or condition the null on them;
  - a detection-selection model.

  The flag is removed only when every attack arm passes.
- **Rejected:** shipping the p unflagged with a caveat in the methods note only. The p would then read as a discovery claim.

## Evidence (before → after: metrics, tests, simulations)
- Statistician probes: `$GIAP_ROOT/stat/m809-null/` (probe*.py, results*.jsonl, summ.py, the sbatch files). In the original runs the attack arms failed at the fast tier: 5 of 6 cases in Slurm 1867349.
- The attack file runs at the fast tier, R = 20. Unflagged, 5 of its 6 cases fail; flagged, it records the rates.

## Consequences
- The report and V15 show the null p as uncalibrated, alongside the ADR-0032 univariate CIs.
- The ΔR² attribution (ADR-0033) is unaffected: it is an effect size, not a test.
