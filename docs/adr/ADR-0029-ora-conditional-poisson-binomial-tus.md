# ADR-0029: Detection ORA uses the conditional Poisson-binomial null, counted in TUs
- Status: accepted-by-khoa (2026-10-21, in session: "3 c 4 a. The rest keep recommended")
- Date: 2026-10-21
- Spec sections: §7.7.3, §8.4
- Needs Khoa: yes (replaces the Wallenius gate the spec names). Answered: option (a).

## Context
Two problems with the Wallenius approach:
- Wallenius's noncentral hypergeometric only approximates the exact null when each gene has its own detection probability. That null is the conditional Poisson-binomial: the number detected in the set, given the total number detected (M8-01 lead decision 5).
- Genes in one TU are detected together. Counting them as independent genes is anti-conservative.

## Decision
- The ORA p-value is computed under the conditional Poisson-binomial null.
- The unit is the TU: a TU is detected when any of its genes is, with probability 1 − ∏(1 − p_detect).
- Wallenius (M8-01) stays as the classic weighted p (`p_classic`, `odds_weight`) and for golden parity with R.
- M8-04 implements this. Until its calibration passes, ORA rows carry no calibrated claim.

## Evidence (before → after: metrics, tests, simulations)
The M8-01 statistician review. M8-04's calibration arms are the gate.

## Consequences
M8-04's design changes. M8-15 covers TU-level ORA.
