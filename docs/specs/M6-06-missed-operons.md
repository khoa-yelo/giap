# M6-06: Missed operons — NO-pair scores, CV and SW decoys, q-values and local fdr

Milestone M6 · track R · **pass 1** · depends on: M6-01 (`max_decoy_fdr`, `local_fdr`, `bin_index`),
M6-04 · estimate ≈ 200 code lines (+ ≈ 280 test lines) · **statistical** (label `stats`)

## Spec sections implemented

- §7.4.4 in full: "s = −|Δ| for NO pairs where both genes are detected"; "Decoys. CV pairs and SW
  pairs, both genes detected, each reweighted to the NO pairs' distribution over intergenic-distance
  bins b (`operons.distance_bins_bp`). With π₀ = 1 …"; "Compute the FDR with each decoy set and
  report the larger one … Make the FDR monotone into q-values, and flag pairs at q ≤ 0.1";
  "Local fdr for ranking".
- §5.2 `operons: {decoys: [convergent, sandwich], distance_bins_bp: [-inf, 0, 50, 150, 300, inf],
  pi0: 1.0}`; §8.1 row "Missed operon" (decoy FDR + local fdr, α 0.1); §9.3 `pairs.tsv` `test` value
  `missed_operon_decoy` and the `lfdr` column.

## Public interface

`src/giap/operons/missed.py` (pure)
```python
MIN_NO_PAIRS: int = 20

@dataclass(frozen=True, slots=True)
class MissedOperons:
    rows: pd.DataFrame          # pair_id, test, stat, p_value, q_value, lfdr, method, flags
    n_targets: int
    n_decoys: dict[str, int]
    bin_counts: pd.DataFrame    # per bin: N_NO, N_cv, N_sw (for the methods note and the report)

def missed_operons(pairs: pd.DataFrame, cfg: GiapConfig) -> MissedOperons: ...
```

1. **Targets.** NO pairs with both genes eligible (§7.3) and detected, both `abundance_log2`
   finite. `stat = s = −|Δ|`. Imputed abundances are kept (unlike the t fit): the spec excludes
   imputation only from the fit. Pairs with `rna_between` are kept.
2. **Decoys.** CV pairs and SW pairs under the same both-detected rule. `operons.decoys` selects
   which sets are used; an empty list is a config error.
3. **Bins.** `bin_index(intergenic_bp, operons.distance_bins_bp)` for targets and decoys alike. The
   reweighting is the one implemented in M6-01 step 4; bins with no decoy contribute conservatively.
4. **FDR and q.** `max_decoy_fdr(..., pi0 = operons.pi0)`: the larger of the two decoy FDRs
   element-wise, then monotone, then clipped → `q_value`. `p_value` is **NA** for this test: a decoy
   FDR produces no p (and §9.2's test columns allow NA). `method = "decoy_fdr"`,
   `null_model = "convergent+sandwich_reweighted"`, `n_perm` NA, `effect = -s = |Δ|`,
   `ci_low`/`ci_high` NA.
5. **lfdr.** `local_fdr(|Δ|_targets, |Δ|_decoys, decoy_weights)` from M6-01, for ranking only; never
   used for a call.
6. **Flagging.** `q_value ≤ stats.alpha.operon_flags` (0.1) marks the pair co-abundant; M6-07 turns
   that into `novel_operon_candidate` or `coabundant_despite_terminator`.
7. **Reporting.** `bin_counts` goes to the §9.4 run report so a reader can see the reweighting that
   was applied, and whether a bin was empty of decoys.
8. **No RNG**: deterministic.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| < 20 scored NO pairs | no `missed_operon_decoy` rows; `flags` gains `too_few_no_pairs`; `bin_counts` still written |
| one decoy set empty (e.g. no SW pair on a 2-gene plasmid) | that set is dropped with a log line; the other alone decides; both empty → no rows + `too_few_no_pairs` |
| a bin with targets and no decoys | conservative contribution (M6-01 example 7); may push every q to 1 |
| all Δ identical | every s equal; one tied q for all |
| < 5 values for a KDE | `lfdr` NA, `q_value` unaffected |
| relative-scale abundance | no rows at all (`relative_scale`, as M6-05 step 6) |
| π₀ from config ≠ 1 | honoured, and recorded in the methods note and the manifest |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | 300 NO targets (200 in bin 0, 100 in bin 1), 200 CV decoys (100, 100), threshold with 40 targets and 5 + 2 decoys above it | FDR = **0.3713** from CV (M6-01 example 4) |
| 2 | plus SW decoys (50, 20) with 1 + 1 above that threshold | reported FDR = **0.4342** (SW is larger), `which = sandwich` |
| 3 | a target with the largest −\|Δ\| (Δ = 0.02) | the smallest q in the frame; `lfdr` ≤ its q is **not** asserted (they are different quantities) |
| 4 | 19 NO pairs | no rows; `too_few_no_pairs` |
| 5 | NO pairs and decoys drawn from the same distribution (n = 2,000 each) | no target at q ≤ 0.1 in ≥ 90 % of seeds; mean FDP ≤ 0.1 (the M6-09 criterion in miniature) |
| 6 | 2,000 decoys and 2,000 NO pairs of which 100 are planted with \|Δ\| ≤ 1e-3 SD (M6-01 lead decision 12: at 5e-3 SD the method makes no calls) | ≥ 50 of the planted pairs at q ≤ 0.1 and FDP ≤ 0.1 (power + control on one frame) |
| 7 | every NO pair in bin 4 (> 300 bp) and every decoy in bin 0 | all q = 1; `bin_counts` shows the mismatch, log line names it |
| 8 | `operons.decoys: [convergent]` | only CV used; `which` is `convergent` everywhere |
| 9 | same frame twice | identical rows (deterministic) |

## Acceptance criteria

- [ ] `tests/unit/operons/test_missed.py`: examples 1–4, 7–9 and every degenerate row on hand-built
  frames; examples 5–6 as `slow` tests with fixed seeds.
- [ ] Property: `q_value` is non-increasing in `stat` and lies in [0, 1].
- [ ] Property: dropping a decoy set can only **lower or keep** the reported q (the maximum over
  sets is monotone in the set collection).
- [ ] Property: the frame is invariant to row order and to `pair_id` relabelling.
- [ ] Hand-check: `bin_counts` sums to the number of scored targets and decoys.
- [ ] Performance: 20,000 NO pairs + 20,000 decoys < 2 s.

## Hardening owed (pass 2)

- The gate itself: `missed_operons` FDR ≤ 0.1 + 2 SE under M6-03's generators, plus attack A
  (autocorrelated abundance, ρ = 0.9) and attack B (detection confounded with distance) — M6-09.
- The power grid (\|Δ\| ∈ {0.25, 0.5, 1.0} SD × 20/50/100 planted pairs), recorded.
- π̂₀ instead of 1.0 if the attacks show the FDR is far too conservative (ADR, M6-01 hardening).
- A golden `pairs.tsv` slice for the 50-kb fixture.

## Questions for the lead

1. `p_value` NA for a decoy-FDR test (proposed), or a pseudo-p from the decoy empirical
   distribution? **Proposed default:** NA — §8.1 lists the correction as "decoy FDR (+ local fdr for
   ranking)", and a pseudo-p would invite a second BH.
2. Keep imputed abundances among the targets? **Proposed default:** yes — §7.4.4 excludes them only
   from the t fit, and the `imputed_abundance` gene flag already travels with the pair.
3. Minimum scored NO pairs: 20 (proposed), matching M6-05's rule.

## Lead decisions (2026-10-21)
1. Q1: `p_value` NA; the family is controlled by the decoy FDR (`q_value`) with lfdr for ranking (§8.1).
2. Q2: keep imputed abundances among the targets; the `imputed_abundance` flag travels with the pair.
3. Q3: minimum **20** scored NO pairs, as M6-05.
4. Power: M6-01 lead decision 12 (Khoa) — the family is near-powerless at realistic effects; example 6
   plants |Δ| ≤ 1e-3 SD. Tests assert FDR control and the honest labelling, never power at 0.25–1 SD.
5. ±inf scores are excluded by the engine (M6-01 B1); a log(0) Δ never reaches it as a target.
