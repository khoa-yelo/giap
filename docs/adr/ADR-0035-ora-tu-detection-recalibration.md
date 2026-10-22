# ADR-0035: The ORA TU detection probability is a cross-fitted TU-level recalibration of q = 1 − ∏(1 − p̂)
- Status: accepted-by-khoa (2026-10-22, in session: "TU recalibration (Recommended)")
- Date: 2026-10-22
- Spec sections: §7.7, §8.4, §8.5; ADR-0029; docs/specs/M8-04-category-summary-ora.md (lead decision 4)
- Needs Khoa: yes (lead decision 4: a failed TU-packed arm sends the TU-detectability formula to Khoa). Answered: TU recalibration.

## Context
ADR-0029 counts ORA in TUs. Each TU has detection probability q = 1 − ∏(1 − p̂) over its eligible genes, and the test is the conditional Poisson-binomial law of K given N (M8-04a). That primitive is exact for any q (statistician, #172). The problem is the q formula itself. When genes of a TU share a detection effect, or p̂ is miscalibrated at gene level, q is biased, and the bias grows with TU size.

The statistician's attack used p̂ equal to the exact marginal over a shared TU effect of strength c, with 3,000 p values per cell. Bounds are 0.0136 at α 0.01 and 0.058 at α 0.05. Slurm 1867453 and 1867454; files in `$GIAP_ROOT/stat/m8-04a-review/`.

| Scenario | Categories | Tail | α 0.01 | α 0.05 |
|---|---|---|---|---|
| c 1.5 | large TUs | under | 0.107 | 0.248 |
| c 1.5 | single-gene TUs | over | 0.019 | 0.079 |
| c 2.5 | large TUs | under | 0.278 | 0.490 |
| c 2.5 | single-gene TUs | over | 0.044 | 0.143 |
| c 0, gene-level p̂ intercept +0.5 | large TUs | under | 0.0137 | 0.067 |

With p̂ = p_true, every category type is at level. Random whole-TU categories (spec arm 4) do not show the bias.

## Decision
1. **q is recalibrated.** M8-04b fits a logistic model of observed TU detection on logit q and log(eligible gene count). The fit is cross-fitted on M4-13's TU folds, so a TU's own detection never enters its q. The fitted value is the TU detection probability passed to the conditional law. The statistician's in-sample probe of this fit restored level.
2. **Numerics are part of the decision.**
   - `TuUniverse` carries the log complement Σ log1p(−p̂), with q = −expm1(Σ).
   - `tilt_pair` accepts a complement input, so it never re-derives 1 − q from a rounded q.
   - A p̂ of exactly 0 or 1 from M4-13 is clipped to [1e-12, 1 − 1e-12] and flagged, never treated as certain.
   - The degenerate-case table drops "unreachable from M4-13".
3. **Gated arms.** M8-14 and M8-15 add two arms, both scored by level in each direction:
   - **4b:** the TU-packed generator, with categories built from TUs of at least 4 genes and from single-gene TUs.
   - **4c:** a gene-level p̂ intercept shift of ±0.5.

   The calibration switches from `tilt` to `tilt_pair`, and a tilt-stress level arm is added. If 4b or 4c fails after recalibration, the ORA rows get `ora_uncalibrated`, following the ADR-0034 pattern, until a hardening fix lands.
4. **Rejected options.**
   - Flag and ship the uncalibrated q: ORA would become descriptive only.
   - Gene-level ORA as primary: it ignores operon co-detection, which ADR-0029 was adopted to handle.

## Evidence (before → after: metrics, tests, simulations)
- Before: the table above.
- After: an in-sample recalibration probe at level (statistician, #172). The cross-fitted version is gated by arms 4b and 4c in M8-15.

## Consequences
- M8-04b gains one logistic fit per run, which is cheap. The ORA columns and the §9.3 schema do not change.
- M4-13's p̂ calibration is no longer assumed at TU level. The methods note must say that conditioning on N cancels only a common TU-level odds factor.
