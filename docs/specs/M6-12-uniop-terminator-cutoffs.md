# M6-12: **BLOCKED** — UniOP AUROC and cut-off, the distance-model refit, the TerminatorNet cut-off, ADR-0022

Milestone M6 · track V (validation) · **blocked on M6-11** · depends on: M6-11, merged M2-12
(UniOP adapter) and M2-15 (TerminatorNet parser) · estimate ≈ 200 code lines in `validation/`
(+ ≈ 200 test lines) · **statistical** (label `stats`) · produces a protected-file change and an ADR

## Status

Blocked for exactly as long as M6-11 is. Nothing in M6-01 … M6-10 depends on it: every cut-off it
would change is read at run time from config or `src/giap/data/calibration.yaml`, which today
carries documented placeholders (`operon_distance {d0: 50, s: 20, calibrated: false}`,
`terminatornet_yes {value: 0.5, calibrated: false}`). **No `uncalibrated` flag is introduced** —
ADR-0021 Amendment 1 dropped it; provenance fields carry the information instead.

## Spec sections implemented (when unblocked)

- §15.4 steps 3–5 in full: UniOP — "Record the AUROC. Choose the smallest threshold with precision
  ≥ 0.90 (if none reaches 0.90, the threshold with maximum F1). Refit the distance-model parameters
  d₀ and s." TerminatorNet — "Precision = the share of calls in intergenic regions downstream of
  genes that lie within ±50 bp of a RegulonDB terminator on the same strand. Choose the smallest
  confidence with precision ≥ 0.90 (else maximum F1). RegulonDB is incomplete, so this precision is
  a lower bound; say so in the ADR." And step 5: "Record the cut-offs in
  `src/giap/data/calibration.yaml`, with RegulonDB version, date, AUROC, precision and recall. The
  ADR goes to Khoa (protected file)."
- §7.4.6 "Truth check" (the AUROC is the M6 gate's third clause; the UniOP-plus-co-abundance
  comparison needs K-12 proteomics and stays deferred — §13 "Adding co-abundance does not lower
  K-12 operon AUROC moves to the full validation").
- §7.1.3 (the distance model p(operon | d) = 1/(1 + exp((d − d₀)/s))).

## Scope

1. `validation/calibrate_cutoffs.py`:
   - join M6-11's `pair_truth.tsv` to UniOP's probabilities on the recorded K-12 run;
   - ROC and AUROC (with a bootstrap CI over pairs), the precision/recall curve, the chosen
     threshold by the §15.4 rule, and the precision and recall at it;
   - refit d₀ and s by maximum likelihood of the logistic distance model on the same labelled
     pairs, reporting the fitted values with CIs and the AUROC of the distance model alone (which
     is the fallback's quality, and worth knowing);
   - TerminatorNet: match calls to `terminator_truth.tsv` within ±50 bp on the same strand,
     restricted to intergenic regions downstream of genes, and choose the confidence by the same
     rule.
2. Write the results into `src/giap/data/calibration.yaml` (protected): values plus `source`,
   `regulondb_version`, `date`, `auroc`, `precision`, `recall`, `n_pairs`, `calibrated: true`.
3. `docs/adr/ADR-0022-uniop-terminatornet-cutoffs.md`: the numbers, the decision, the explicit
   statement that RegulonDB's incompleteness makes the terminator precision a **lower bound**, and
   what changes downstream (the default `annotation.uniop.pair_threshold`, the default
   `annotation.terminators.terminatornet_yes`, the distance-model constants). Goes to Khoa.
4. A regression test that `giap.data.calibration.load()` returns the new values and that the
   defaults in `defaults.yaml` and `calibration.yaml` agree.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no threshold reaches precision 0.90 | the maximum-F1 threshold, and the ADR says so |
| fewer than 500 labelled pairs | the calibration is not claimed; escalate (M6-11 example 5 guards this) |
| UniOP probabilities tie at the chosen threshold | the smallest threshold rule is applied to the sorted distinct values; ties documented |
| the distance refit does not converge | keep d₀ = 50, s = 20 with `calibrated: false`, and say so in the ADR |
| TerminatorNet calls absent (tool not run on K-12) | the terminator half is deferred; the UniOP half can still land — they are separable |
| a new RegulonDB release appears later | re-run; the provenance fields make the change visible |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | 3,000 labelled pairs, UniOP p | an AUROC with a bootstrap CI, recorded to 3 decimals |
| 2 | precision ≥ 0.90 first reached at p = 0.62 | `pair_threshold` = **0.62**, with its precision and recall |
| 3 | precision never reaches 0.90 | the max-F1 threshold, flagged in the ADR |
| 4 | the distance refit | d₀ and s with CIs; the distance model's own AUROC recorded beside UniOP's |
| 5 | TerminatorNet confidence ≥ 0.90 precision first at 0.71 | `terminatornet_yes` = **0.71**, with the "lower bound" caveat |
| 6 | the shipped `calibration.yaml` after the run | `calibrated: true`, with version, date, AUROC, precision, recall, n |
| 7 | a GIAP run after the change | the manifest records the new provenance; no `uncalibrated` flag anywhere |

## Acceptance criteria

- [ ] AUROC recorded in PROGRESS.md and in `docs/methods/operons.md` — the gate's third clause.
- [ ] `calibration.yaml` updated with full provenance; `giap.data.calibration.load()` test updated.
- [ ] ADR-0022 written, with the lower-bound statement, and **accepted by Khoa** (protected file).
- [ ] Unit tests of the threshold-choice rule (examples 2, 3, 5) on synthetic ROC data — these can
  be written **now**, while the issue is blocked, since they need no RegulonDB data.
- [ ] A re-run of M6-09's fast suites after the defaults change (expected: unchanged, because the
  simulations plant their own probabilities).

## Hardening owed

- "Adding co-abundance does not lower K-12 operon AUROC" — deferred to the full validation by §13;
  recorded here so it is not forgotten.
- A second strain's truth (Nissle has no RegulonDB) is out of scope and stays out.

## Questions for the lead / Khoa

1. (Khoa) ADR-0022 is a protected-file change and needs acceptance once numbers exist.
2. (Lead) May the threshold-choice unit tests (synthetic ROC) be merged **before** the data is
   reachable, so the blocked issue shrinks to a data-fetch plus a number? **Proposed default:**
   yes — it is the only part of M6 that can make progress while blocked.
3. (Khoa) If RegulonDB stays unreachable, is M6 declared done on its simulation half with §15.4
   deferred to M10? The milestone table's gate text would then need an amendment.
