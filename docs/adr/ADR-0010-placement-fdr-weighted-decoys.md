# ADR-0010: Placement FDR with per-peptide π̂ weights each decoy hit by its own π̂ (M3-10)
- Status: accepted
- Date: 2026-10-14
- Spec sections: §7.2.5, G3, §8.4; docs/specs/M3-10-placement-fdr.md; M3-index lead decision (M3-10 Q2)
- Needs Khoa: no (reduces exactly to the §7.2.5 formula FDR_m(e) = π̂ (D_m(e) + 1) / T_m(e) when π̂ is
  constant; replaces only the lead's provisional rule for the per-peptide-π̂ case the spec leaves open)

## Context
§7.2.5 gives FDR_m(e) = π̂_wrong (D_m(e) + 1) / T_m(e) with one π̂ per stratum. De novo input has a
π̂ per peptide (score bin). The lead's provisional rule (M3-10 Q2) used the mean π̂ over TARGETS with
E ≤ e. The M3-10 test-author found it anti-conservative and the statistician confirmed: targets are
mostly correct, low-π̂ peptides, while decoy hits are chance hits of all searched peptides.

## Decision
FDR_m(e) = (W_m(e) + c) / T_m(e), where W_m(e) = Σ π̂_j over decoy hits with E ≤ e, each weighted by
its own peptide's π̂ (missing → 1; strata m ≥ 1 → 1), and c = max π̂ over the peptides of that decoy
type searched. One list ordered by E per stratum, same q-value monotonisation, ties and clipping as
M4-01. `mapping_fdr` rows report `fdr_stratum` = (W + c)/T and `pi_wrong` = (W + c)/(D + 1).
Signed in docs/methods/placement-fdr.md ("π̂_wrong in the stratum FDR").
Not adopted: separate lists per score bin (valid and more powerful on large inputs, but breaks "q
monotone in E within a stratum"). Database input keeps π̂ = mean PEP (estimator (a)); M3-21 measures
the PEP–length correlation that could make it anti-conservative and revisits this if it is material.

## Evidence (mean FDP at q = .01 / .05 / .10; N = 4000, R = 300; * = above q + 2 SE)
| Scenario | Q2 (old) | mean π̂ over all searched | ADR-0010 |
|---|---|---|---|
| 10 score bins | .0143* / .069* / .130* | .0098 / .050 / .100 | .0095 / .050 / .100 |
| 2 bins (π = .02 / .9) | .035* / .119* / .175* | .0097 / .050 / .100 | .0095 / .050 / .100 |
| 2 bins, chance-hit rate tied to π | .054* / .168* / .253* | .0163* / .085* / .172* | .0093 / .049 / .100 |
| π̂ = 1.5 π | .0099 / .049 / .094 | .0072 / .037 / .073 | .0069 / .036 / .073 |
| constant π̂ | .0100 / .050 / .097 | same | same |
π̂ = 0.5 π fails under every rule (≈ 2q): π̂ must come from the conservative bound (M3-11).
Attack: tests/calibration/attacks/test_attack_placement_pi.py.

## Consequences
- c ≈ 1 when any searched peptide has π̂ ≈ 1: small lists accept nothing at q = 0.01 (N = 150 → 0);
  at N = 4000, 1063 accepted vs 1086 for the (invalid in general) mean-over-all variant.
- `stratum_qvalues(e_target, e_decoy, pi_decoy, pi_plus)` replaces the target-mean interface; M3-10
  tests that pinned Q2 change accordingly; M3-index decision for M3-10 Q2 is superseded.
