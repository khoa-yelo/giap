# ADR-0022: UniOP and TerminatorNet cut-offs and the distance-model refit on RegulonDB 14.5.0 K-12
- Status: accepted-by-khoa (2026-10-23, max-F1; M6-12 lead decision 10; 2026-10-24, UniOP `pair_threshold` 0.5)
- Date: 2026-10-23 (run 2026-10-24 01:05 UTC, Slurm 1879748)
- Spec sections: §15.4 steps 3–5, §7.1.3, §7.1.5, §7.4.5, §7.4.6, §7.4.7; spec
  `docs/specs/M6-12-uniop-terminator-cutoffs.md` with lead decisions 1–10; methods
  `docs/methods/operons.md` › "K-12 cut-off calibration (M6-12)"
- Needs Khoa: yes (a protected file, `src/giap/data/calibration.yaml`, and two defaults change).
  Khoa decided the cut-off rule on 2026-10-23: **max-F1 as §15.4 specifies, no one-SE rule**.

## Context

§15.4 asks for three numbers from RegulonDB on K-12 MG1655, recorded in `calibration.yaml` with the
RegulonDB version, date, AUROC, precision and recall:
1. the UniOP AUROC and pair cut-off: "the smallest threshold with precision ≥ 0.90 (if none
   reaches 0.90, the threshold with maximum F1)";
2. the distance-model refit, p(operon | d) = 1/(1 + exp((d − d₀)/s)) (§7.1.3, placeholder 50 / 20);
3. the TerminatorNet confidence cut-off, by the same rule, with precision against RegulonDB
   terminators within ±50 bp on the same strand.

Release of record: **RegulonDB 14.5.0** (2026-01-28), as Khoa decided on 2026-10-21 (§15.3 says
"v12"). Until now `calibration.yaml` held the placeholders (`d0: 50, s: 20`,
`terminatornet_yes: 0.5`, both `calibrated: false`), every run carried the flags
`uncalibrated:operon_distance` and `uncalibrated:terminatornet_yes`, and
`annotation.uniop.pair_threshold` was the uncalibrated 0.5.

The cut-off rule follows lead decision 1 (the **stable** reading of "smallest t with precision ≥
0.90": every candidate at or above t* has precision ≥ 0.90), with max-F1 as the fallback, ties to
the smallest t. Khoa's decision 10 keeps max-F1 exactly as §15.4 specifies; the one-SE rule the
statistician raised for TerminatorNet's flat F1 curve (methods R3) is **not** adopted, and the
flatness is reported below as context only.

## Decision

Record these values (all from `cutoffs.json`, unrounded; the floats are Python's shortest
round-trip `repr` and must not be rounded, because the pipeline compares with `≥`):

| Key | Value | Rule that fired | Mirrors |
|---|---|---|---|
| `uniop_pair_threshold.value` (new) | **0.5** | Khoa, 2026-10-24 (the `precision_0.90` rule gave 0.529415) | `defaults.yaml` `annotation.uniop.pair_threshold`; `config/model.py` `UniopConfig.pair_threshold` |
| `operon_distance.d0` | **102.68667744326773** | IRLS refit, converged | — (read by `steps/operons.py`) |
| `operon_distance.s` | **47.16037265004493** | IRLS refit, converged | — |
| `terminatornet_yes.value` | **0.9999999999961204** | `max_f1` (no candidate reaches 0.90) | `defaults.yaml` `annotation.terminators.terminatornet_yes`; `config/model.py` |

All three entries become `calibrated: true` with the ADR-0022 provenance. The proposed file is
`$GIAP_ROOT/stat/m612c/calibration.yaml.proposed` (with `calibration.yaml.diff`); the lead
commits it in PR-C with `Gate-Change-Approved: ADR-0022`.

**Downstream changes.**
- TUs and pair classes keep `pair_threshold` 0.5 (`classify`, `build_tus`, the flag code); the
  value is now recorded as calibrated.
- The distance fallback (UniOP missing or failed) uses d₀ = 102.69 bp, s = 47.16 bp instead of
  50 / 20: much more permissive at 40–150 bp (see the binned table).
- Terminator "yes" (§7.1.5, §7.4.5) needs TerminatorNet confidence ≥ 0.9999999999961204, that is a
  saturated network output; 727 of 2,447 eligible K-12 calls (29.7 %) pass, against all 2,447 at
  0.5.
- The run flags `uncalibrated:operon_distance` and `uncalibrated:terminatornet_yes` disappear;
  `manifest.calibration` gains a third key, `uniop_pair_threshold`.
- The §8.6 `full` sensitivity variants stay at UniOP thresholds 0.3 and 0.7 (decision 2, Q7):
  t* = 0.529 < 0.7, so no backlog note is due.

## Evidence

**Run.** `validation/calibrate_cutoffs.py` (part B, #264) at worktree HEAD 3db69d9b, Slurm
**1879748** (wheat-05, 59 s wall), output `$GIAP_DATA_DIR/regulondb/14.5.0/calibration/` (not
committed, R6/R11). Rerun **1879749** with the same seed into `$GIAP_ROOT/stat/m612c/rerun/`:
`cutoffs.json` identical apart from `date`; the three TSVs byte-identical.
- `cutoffs.json` sha256 **edcf9044eafe5f1f3f31aab93c759a035efd25ffa240bd199a637640e843d6eb**
- `uniop_curve.tsv` 87054f2e…11, `terminatornet_curve.tsv` b4de9758…1f,
  `terminatornet_calls.tsv` 705903eb…97
- seed 20261014 (`giap.rng.for_step("m6_12_cutoffs", seed)`), B = 2,000, R = 200, cluster
  `directon`, 773 directons holding labelled pairs.
- Identities checked by the driver: truth SHA-256s `TUSet` c96c8a3b…39fc, `TerminatorSet`
  9d916ffe…f366, `GeneProductSet` 3c2c6be3…9dde; genome (normalised) 918f35c7…14f2; UniOP pin
  3ed14dfc (ADR-0008), `uniop.pred` 7d852a93…5ee0, `uniop_cds.gff` ea587089…2906 (map rebuilt by
  coordinates, 4,372 / 4,372); TerminatorNet commit 72180309 (ADR-0008), `terms.tsv.gz`
  9be61c75…a2ff.

**UniOP** (the M6 gate's third clause, "UniOP-alone AUROC recorded", §13 row M6). 95 % CIs are
directon-bootstrap percentile CIs; the iid (single-pair) line is for comparison only.

| Quantity | Estimate | 95 % CI (directon) | iid CI |
|---|---|---|---|
| labelled pairs L; P / N; without p | 2,739; 1,730 / 1,009; 0 | — | — |
| **AUROC** | **0.911** (0.910877249265283) | [0.898, 0.923] | [0.900, 0.922] |
| t* (`precision_0.90`; literal reading identical) | 0.529415 | bootstrap re-selected t*: 2.5 / 50 / 97.5 % = 0.457, 0.529, 0.610 | — |
| precision at t* (fixed t*) | 0.9017 (TP 1,467, FP 160) | [0.884, 0.918] | [0.887, 0.916] |
| recall at t* (fixed t*) | 0.8480 | [0.827, 0.868] | [0.831, 0.864] |
| F1 at t* | 0.8740 | — | — |
| precision of the whole rule, 10-fold directon CV | 0.9001 | — | — |

- The precision CI holds t* fixed: it is a CI for the precision *of the threshold 0.529415*, not
  for the precision *achieved by the rule*. The rule chose t* because the observed precision there
  is ≥ 0.90, so the point estimate is optimistic; the cross-validated precision of the rule
  (0.9001) puts the optimism at about 0.002, against a CI half-width of 0.017.
- t* itself is known only to about ±0.08 (the re-selection quantiles); the *precision* of the cut-off
  is what is stable.
- Precision 0.90 is **cumulative, not per-pair**: it is the share of operon pairs among all pairs
  with p ≥ t*. Pairs just above the cut-off (p in [0.529, 0.65), 122 pairs) are operon pairs 49 %
  of the time.
- `pair_threshold` is **K-12-tuned**. Precision depends on prevalence (K-12 labelled π = 0.632):
  at the same TPR / FPR, a genome with π = 0.5 would see about 0.84 and one with π = 0.4 about 0.78.
  AUROC is prevalence-free but not genome-free.
- UniOP's training data: **unknown here** (not verified in M6-12). If UniOP's development used
  K-12 operons from RegulonDB, this AUROC is partly in-sample for the tool.
- **Khoa's decision, 2026-10-24: UniOP uses 0.5.** On the K-12 curve 0.5 gives precision 0.896,
  recall 0.853 and F1 0.8739 (TP 1,476, FP 172), against 0.902 / 0.848 / 0.8740 at the
  `precision_0.90` threshold 0.529415: the same operating point, so the round value is kept and
  `rule` is recorded as `khoa_fixed_0.5`. Khoa also asked whether TerminatorNet could use 0.5; it
  cannot without a large change. At 0.5 every eligible candidate is called (recall 1.0), FP rises
  from 563 to 2,154 and precision falls from 0.226 to 0.120, close to the 0.091 coincidence rate,
  so TerminatorNet keeps the max-F1 value.
- **Context for decision 10, not the recorded value.** The max-F1 threshold on the UniOP curve is
  0.25114 (F1 0.8743, precision 0.848, TP 1,562, FP 281), against F1 0.8740 at 0.529415. §15.4
  applies max-F1 only when no threshold reaches precision 0.90; UniOP reaches it, so the recorded
  rule is `precision_0.90`. F1 is flat here too: every candidate from 0.069 to 0.707 is within 1 %
  of the maximum.

**Distance model** (the fallback when UniOP is missing).

| Quantity | Estimate | 95 % CI (directon) | iid CI |
|---|---|---|---|
| d₀ (bp) | 102.69 | [94.59, 111.47] | [96.08, 110.31] |
| s (bp) | 47.16 | [41.19, 53.64] | [41.85, 53.12] |
| converged / reason | true / none | — | — |
| AUROC of −d | 0.907 (0.9067313255841931) | [0.892, 0.920] | [0.894, 0.918] |
| AUROC(UniOP) − AUROC(−d), paired | 0.0041 | [0.0009, 0.0077] | [0.0007, 0.0079] |
| log-loss on L, refit / placeholder (50, 20) | 0.407 / 0.616 | in-sample, no CI | — |

- The refit used the damped Newton from (logit ȳ, 0) after the exact separation check (lead
  decision 6, B1); the plain Newton of the spec's first draft diverges on these data.
- The percentile CIs of d₀ and s under-cover slightly (about 0.92 on synthetic data, part A); they
  are recorded, not gated.
- UniOP beats distance alone by 0.004 AUROC; the paired CI excludes 0, so the gain is real but
  small.
- **The fallback does not deliver 0.90 precision at `pair_threshold`.** p_dist ≥ 0.529415 means
  d ≤ 97.1 bp; on L that calls 1,878 pairs with precision 0.839 and recall 0.910.
- The linear-logit family misfits K-12 (it is fixed by §7.1.3). Observed vs fitted operon share by
  intergenic-distance bin (counts and shares only):

| d (bp) | n | observed | refit | placeholder (50, 20) |
|---|---|---|---|---|
| < −20 | 47 | 0.745 | 0.947 | 0.983 |
| [−20, 0) | 566 | 0.979 | 0.907 | 0.939 |
| [0, 20) | 576 | 0.965 | 0.879 | 0.884 |
| [20, 40) | 222 | 0.914 | 0.830 | 0.750 |
| [40, 60) | 173 | 0.549 | 0.751 | 0.495 |
| [60, 80) | 172 | 0.442 | 0.669 | 0.278 |
| [80, 100) | 136 | 0.449 | 0.566 | 0.122 |
| [100, 150) | 289 | 0.273 | 0.398 | 0.033 |
| [150, 200) | 165 | 0.164 | 0.189 | 0.003 |
| [200, 300) | 218 | 0.119 | 0.057 | 0.000 |
| ≥ 300 | 175 | 0.103 | 0.004 | 0.000 |

  The refit is too flat: it under-predicts at −20 to 40 bp, over-predicts at 40 to 150 bp, and is
  too confident beyond 300 bp, where about 10 % of labelled pairs are still operon pairs. It is still far better
  than the placeholder, which puts almost no operon mass beyond 100 bp.

**TerminatorNet.** Wilson 95 % intervals, plus the window-cluster bootstrap interval (each gene's
downstream window is a cluster) for precision.

| Quantity | Estimate | 95 % CI |
|---|---|---|
| calls; eligible (inside a same-strand downstream window) | 8,294; 2,447 | — |
| RegulonDB terminators; window-eligible | 375; 295 | — |
| t* (`max_f1`; no candidate reaches 0.90 under either reading, `literal_threshold` null) | 0.9999999999961204 | — |
| calls ≥ t*: TP / FP | 164 / 563 (727 calls) | — |
| precision at t* (a lower bound, see below) | 0.226 | Wilson [0.197, 0.257]; window [0.196, 0.258] |
| coincidence rate c at t* (placement null, R = 200) | 0.0906 | MC SE 0.0008 |
| coincidence-corrected precision bound (obs − c)/(1 − c) | 0.148 | — |
| … at the Wilson lower limit / at the window lower limit | 0.117 / **0.115** (`precision_lower_bound_window`, the 97.5 % one-sided bound) | — |
| recall over window-eligible truth | 0.556 (164 / 295) | [0.499, 0.612] |
| recall over all 375 terminators | 0.443 | [0.393, 0.493] |
| recall, rho-independent | 0.580 (156 / 269) | [0.520, 0.637] |
| recall, rho-dependent | 0.308 (8 / 26) | [0.165, 0.500] |
| F1 at t* (call-level) | 0.322 | — |
| AUROC matched vs unmatched eligible calls (diagnostic, not recorded in the yaml) | 0.677 | not estimated |

- **RegulonDB is incomplete, so this precision is a lower bound**, except for coincidental
  matches: a false call can land within 50 bp of a listed terminator by chance, because both
  concentrate just past stop codons. The offset-preserving placement null puts that chance at
  c = 0.091, and the coincidence-corrected bound is 0.148 (0.115 at the window-bootstrap 97.5 %
  lower limit). The true precision of `terminatornet_yes` on K-12 lies above about 0.115 with 97.5 %
  confidence; RegulonDB cannot say how far above.
- The 0.90 target is unreachable against RegulonDB: the best precision at any candidate is 0.229.
  **`terminatornet_yes` is a max-F1 cut-off, not a 90 % cut-off**, and carries no precision
  guarantee on any genome.
- Recall is low by construction: RegulonDB includes Rho-dependent terminators, which TerminatorNet
  does not model.
- **Flat F1 (context only; Khoa decided max-F1, 2026-10-23).** F1 rises from 0.214 at 0.5 to 0.310
  at 0.9 and is flat above: 0.314 at 0.99, 0.313 at 0.999, 0.315 at 0.99999, maximum 0.322 at
  0.9999999999961204. The 45 candidates within 1 % of the maximum all lie in
  [0.9999999999689556, 0.9999999999999982]; the 210 within 5 % lie in [0.8814, 1.0]. The argmax therefore sits on
  saturated network outputs, and its decimals beyond about the 4th are noise; on another genome or
  TerminatorNet build the float pattern near 1.0 can differ. The value is recorded exactly, which is
  correct for K-12. A one-SE rule was considered and not adopted (decision 10).
- The rule's F1 is call-level (recall = TP over all matched eligible calls), as specified; the
  reported recall is truth-level.

**Gate and tests run with the proposed values** (scratch copies under `$GIAP_ROOT/stat/m612c/`;
`calibration.yaml` itself was not edited):
- `tests/calibration -m calibration`, targeted to the cut-off, coincidence, operon-FDR and
  operon-simulation suites with their attacks (7 files, 20 tests, including the directon-clustered
  AUROC null at R = 500): **20 passed** on the unchanged branch (Slurm 1879917), with only the
  proposed `calibration.yaml` (1879918), and with the full PR-C mirrors (1879919). No tolerance
  was changed. The untargeted `tests/calibration -m calibration` hit the 3 h limit on one shared
  16-CPU node in both variants (1879781, 1879783) at about 95 %, with one failure each that the
  log does not name; no calibration suite reads `calibration.yaml` or the two defaults, so that
  run is not evidence for or against these values, and the lead's PR-C fast suite remains the check.
- C0 and the tests that read the defaults: the C0 agreement test
  `test_calibration_and_defaults_agree` passes non-vacuously when `defaults.yaml` and
  `config/model.py` carry the mirrors (Slurm 1879784), and fails, as designed, when only
  `calibration.yaml` changes (Slurm 1879782). `test_model_defaults_equal_defaults_yaml` passes.
  Other failures are listed under Consequences.

**Multiple testing.** M6-12 produces no p or q values; every output is an estimate with a CI or a
cut-off from a fixed rule, so there is no §8.1 family. The exceptions without a CI are the
in-sample log-losses (descriptive), t* (its bootstrap quantiles are given) and the TerminatorNet
AUROC diagnostic.

## Consequences

- One genome, one truth set: every CI here is over resampled directons (or windows) of K-12 and
  describes K-12's uncertainty, not genome-to-genome transfer. d₀ and s transfer worse than AUROC
  (β₀ absorbs K-12's prevalence and intergenic statistics).
- Label noise in RegulonDB biases AUROC and precision towards the base rate (conservative); study
  bias towards well-studied operons may do the opposite.
- A later RegulonDB release reruns the driver into `regulondb/<release>/calibration/` and needs a
  new ADR.
- **PR-C needs test changes beyond C0 (lead work, role guard).** With the proposed values, the
  following tests fail; none is a statistical failure, all pin placeholders or the two-key set:
  1. `tests/unit/test_operons.py::test_calibration_yaml_entries` converts every non-`calibrated`,
     non-`provenance` field to float, so the spec's string fields `regulondb_version`, `date` and
     `rule` raise `ValueError` (spec § "What changes in protected files" vs C0's test). Either the
     test skips string metadata or the yaml folds them into `provenance`.
  2. `tests/unit/manifest/test_content_blocks.py::test_calibration_lists_both_bundled_keys_and_flags_them`
     and `tests/integration/test_mvp03c_manifest_content.py::test_calibration_block[local|nextflow]`
     pin the key set to `{operon_distance, terminatornet_yes}` and (the first) require both
     `uncalibrated:*` flags.
  3. `tests/unit/test_operon_steps.py` (4 tests: `..._falls_back[empty|header_only]`,
     `..._without_uniop_pred_falls_back`, `..._uses_uniop_and_distance_model_for_wrap`) pin the
     distance-model probabilities of d₀ = 50, s = 20.
  4. With the `terminatornet_yes` mirror at 0.9999999999961204, tests that read the default
     implicitly with toy confidences such as 0.9: `tests/unit/operons/test_flags.py` (10 tests),
     `tests/unit/novel/test_novel_evidence.py` (2) and
     `tests/unit/novel/test_novel_evidence_review_180.py` (1). The spec's last acceptance
     criterion asks PR-C to list these; they should pass explicit thresholds.
