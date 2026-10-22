# ADR-0025: GSEA on genomically clustered sets: a flag and a restricted claim in pass 1, a rotation null in hardening
- Status: accepted-by-khoa (2026-10-21, in session: "3 c 4 a. The rest keep recommended")
- Date: 2026-10-21
- Spec sections: §7.7.4, §8.2, §8.5
- Needs Khoa: yes (a claim restriction and a new flag). Answered: option (c) now, (a) in hardening.

## Context
The M8-05 TU-collapsed GSEA permutes TU labels, which assumes TUs are exchangeable. Positional autocorrelation between neighbouring TUs breaks that for sets whose TUs are adjacent on the genome. Real COG, KEGG and GO sets often cluster this way (superoperons, flagellar regions).

## Decision
1. **Pass 1.**
   - Rows of sets whose member TUs are genomically clustered carry a `genomic_cluster` flag.
   - The flag is defined in the M8-05 follow-up, from TU order along each replicon. The statistician signs off the definition.
   - The methods note states the claim: p-values are valid under TU exchangeability, and there is no level claim for flagged sets.
2. **Hardening pass.** A circular-shift (rotation) null for set p-values, as the region tests use, replaces the restriction.
3. The positional-autocorrelation attack in `tests/calibration/attacks/test_attack_gsea_structure.py` is a strict xfail until (2) lands. The scattered-set arm is gated.

## Evidence (before → after: metrics, tests, simulations)
Statistician probe (R = 2,000 per arm, Slurm 1865631, 1865650, 1865653), rejection rate at α = 0.05:
- AR(1) φ 0.5, one run of 15 adjacent TUs: 0.186.
- AR(1) φ 0.5, five runs of 3 TUs: 0.130.
- AR(1) φ 0.2, one run of 15 TUs: 0.085.
- AR(1) φ 0.5, 15 scattered TUs: 0.0555.

## Consequences
- A new flag value in the §9.4 vocabulary, through the schema.
- The report and V14 must not present flagged rows as calibrated.
