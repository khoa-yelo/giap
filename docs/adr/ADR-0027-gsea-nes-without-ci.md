# ADR-0027: GSEA reports NES without a confidence interval (§8.2 deviation)
- Status: accepted-by-khoa (2026-10-21, in session: "3 c 4 a. The rest keep recommended")
- Date: 2026-10-21
- Spec sections: §8.2, §7.7.4, §9 `tables/enrichment_gsea.tsv`
- Needs Khoa: yes (a deviation from §8.2). Answered: option (a).

## Context
§8.2 asks every effect for a CI. fgsea's normalised enrichment score has none.

## Decision
- NES is a descriptive effect size. Every GSEA row carries `ci_unestimated`.
- `log2err` is a Monte Carlo precision label only.
- No second effect size with a CI is added.

## Evidence (before → after: metrics, tests, simulations)
The M8-05 tests pin `ci_unestimated` on every row.

## Consequences
The report shows NES without an interval.
