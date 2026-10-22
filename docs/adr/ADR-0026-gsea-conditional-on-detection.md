# ADR-0026: Abundance GSEA is conditional on detection (restricted claim and flag)
- Status: accepted-by-khoa (2026-10-21, in session: "3 c 4 a. The rest keep recommended")
- Date: 2026-10-21
- Spec sections: §7.7.4, §8.2
- Needs Khoa: yes (a claim restriction). Answered: option (c).

## Context
Only detected genes are ranked. When a set's genes are harder to detect, its detected members are selected for high abundance, so the set looks enriched without any true difference.

## Decision
- GSEA results are stated as "conditional on detection, with no p_detect adjustment" in the methods note, the report and V14.
- Rows of sets whose detection rate differs from the background carry a flag. The M8-05 follow-up defines the flag and its threshold, and the statistician signs them off.
- Censored ranks for undetected genes (`blom_scores(x, censored_mask)`) and a p_detect-weighted null are left for a later ADR.
- The detectability-confounding attack is a strict xfail.

## Evidence (before → after: metrics, tests, simulations)
Statistician probe, rejection rate at α = 0.05:
- detection logistic(2(a − 1.5·[in set] + 0.5)): 0.786;
- weaker confounding (shift 1.0, slope 1): 0.137.

## Consequences
- One new flag value.
- "No significant set" is never reported as "no enrichment".
