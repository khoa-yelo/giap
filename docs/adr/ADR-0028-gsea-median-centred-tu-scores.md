# ADR-0028 (ADR-B): GSEA ranks median-centred TU scores
- Status: accepted-by-khoa (2026-10-21, in session: "3 c 4 a. The rest keep recommended")
- Date: 2026-10-21
- Spec sections: §7.7.4
- Needs Khoa: yes (a methods choice). Answered: option (a).

## Context
Each TU's score is its mean log2 abundance. All of these are positive, so under `scoreType = std` the test becomes effectively one-sided.

## Decision
- Rank TU scores centred on their median (M8-05 lead decision 1).
- Imputed values count in the TU mean (decision 2); the §8.6 reruns vary the imputation.

## Evidence (before → after: metrics, tests, simulations)
- The statistician's probe: centring changes the weights only, so the level under exchangeability is unchanged.
- Rejection at α = 0.05 is 0.0526 under the null (R = 10,000) and 0.0475–0.049 at 20–60 % median imputation.
- M8-15 runs a centred-versus-raw level and power arm.

## Consequences
None beyond M8-05.
