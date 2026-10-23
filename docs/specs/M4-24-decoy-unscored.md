# M4-24: Unscored and partly scored decoys in gene detection

Milestone M4 (fix) · track G · depends on: M4-12 (merged) · estimate ≈ 60 code lines (+ ≈ 250 test
lines, + one calibration attack) · **can start now** · defect found by the statistician on main in
`src/giap/genetable/detection.py`. (M4-22 and M4-23 are unused; this issue takes the number the lead
gave it.)

## The defect

The picked tier ranks the best accepted target and the best paired decoy of each gene on the
lexicographic key (q ascending, PEP ascending, score descending), with NA mapped to +∞, so NA ranks
last (`_key_rank`). The ranking uses all three columns, whichever side carries them. Two consequences:

1. **All-unscored decoys.** The generic peptide TSV (§5.3: `q_value`, `pep` and `score` are all
   optional) allows decoy rows with no q, PEP or score. giapsim's M3-19 proteomics writes exactly that
   (`tests/sim/giapsim/proteomics.py`: decoy rows `q_value = pep = score = NA`). Every such decoy ranks
   below every target, so D(t) = 0 at every target threshold, and every target-won gene gets
   q = 1/k. That is a picked FDR computed as if there were no decoys, yet the run reports
   `detection_tier = decoy_picked` and raises no warning. Evidence: Slurm 1869554
   (`local/scratch/qfloor/`). On sim100k (6 seeds × π_wrong ∈ {0.01, 0.05}) the base decoys are all
   NA (`decoy_score_na: true`), `n_decoy_wins` is 0–3, and every decoy win is on a gene with no
   target (`decoy_wins_no_target` = `n_decoy_wins`).
2. **Partly scored decoys.** If the decoys carry a score but no q, and the targets carry q, then
   every decoy sits below every target on the first key column. D is again 0, which is silently
   anti-conservative. The same happens one column down: a decoy with the target's q but NA PEP loses
   the q tie to a target that has a PEP, although on the columns both sides carry it is a tie (and
   ties go to the decoy, §7.3.2 step 3).

DIA-NN input is not affected. The converter (`src/giap/io/converters/diann.py`) requires `q_value`
(`Global.Q.Value`, else `Q.Value`) and maps `pep` from `PEP`, both from the same report columns for
targets and decoys, and it maps no `score`. Decoys therefore carry exactly the targets' key columns.
The Casanovo, InstaNovo and Cascadia converters write no decoy rows, so they are not affected either.

## Spec sections implemented

- §7.3.2 picked steps 1–3: "s_g is the best accepted target peptide of gene g, by −log₁₀ q (ties
  broken by lower PEP, then higher `score`)." "The decoy score s̃_g is the best decoy peptide paired
  with gene g." "Per gene, keep max(s_g, s̃_g); a tie counts as a decoy." Targets and decoys can be
  compared only on a scale both carry. This issue makes the comparison use only the columns both sides
  carry.
- §7.3.2 tier table: `decoy_picked` requires decoy rows that can compete. "Otherwise, q-values or
  PEPs exist" gives `mayu`.
- §8 preamble: "Claims that cannot have one [an error estimate] are labelled `unestimated`, never
  silently trusted"; R4. A `decoy_picked` label on an estimate with D ≡ 0 is a silent claim.
- §8.3 `picked_fdr(target_best, decoy_best)`: "picked competition, ties go to the decoy" (unchanged).
- §5.3 C9 columns `q_value`, `pep`, `score` are all optional; `is_decoy` "requires `protein_ids` for
  gene-level pairing". No column is required on decoys.
- ADR-0009: f = paired share of the decoys that enter pairing; decoy wins count 1/f (unchanged).
- §7.3.1 and §5.2 `input.q_value_max` (the decoy q filter; see Rules › q filter).
- M4-12 lead decision 1, M4-index lead decision 4: "Decoys pass the same engine q threshold as
  targets before picking."

## Public interface

`src/giap/genetable/detection.py` (pure). Additions only, and existing signatures are unchanged:
```python
KeyColumn = Literal["q_value", "pep", "score"]
KEY_COLUMNS: tuple[KeyColumn, ...] = ("q_value", "pep", "score")   # lexicographic order

def shared_key_columns(targets: pd.DataFrame, paired: pd.DataFrame) -> tuple[KeyColumn, ...]: ...
    # the columns of KEY_COLUMNS, in that order, that are "carried" by targets (>= 1 non-NA value)
    # and "complete" on paired (no NA value; paired = pairing rows with gene not NA)
```
- `choose_tier(...)` keeps its signature. It can now also return the warning `decoy_unscored`, and it
  applies the fallback below.
- `detect(...)` keeps its signature. The picked ranking (`_key_rank` / `_picked_q`, private) takes
  the shared key columns.
- New warning string `decoy_unscored` in `Detection.warnings`.

## Rules

Definitions. The **target side** is the `accepted` table when `choose_tier` gets it (`detect` always
passes it), else the C9 non-decoy rows (the existing convention of `choose_tier`). The **paired
decoys** are the `pair_decoys` rows with `gene` not NA: after the q filter, with contaminant decoys
already left out (M4-12 decision 12). A column is **carried** by a side when at least one of its rows
has a non-NA value in it. It is **complete** on a side when no row has NA in it.
C_T = the key columns carried by the target side. S = `shared_key_columns(target side, paired
decoys)` = the columns of C_T that are complete on the paired decoys.

1. **Ranking key.** In the picked tier, targets and paired decoys are ranked by the lexicographic key
   restricted to S, in `KEY_COLUMNS` order and with the same directions (q ascending, PEP ascending,
   score descending). Columns outside S are ignored on both sides. Target rows with NA in a column of
   S keep NA → last, as now: a target with no value cannot beat a decoy that has one, which errs
   conservative. Paired decoys cannot have NA in S, by construction. Dense rank over the union,
   picking, ties → decoy, ADR-0009's 1/f and the monotone q are all unchanged.
2. **Warning.** `choose_tier` adds `decoy_unscored` when the pairing guard passes (share ≥ the M4-12
   cut-off), the target side has ≥ 1 row, and S ≠ C_T or S = ∅. In words, the paired decoys lack,
   wholly or on some rows, a key column that the targets carry, or no key column is shared. A column
   that only the decoys carry (C_T lacks it) does not warn, because dropping it removes no target
   information.
3. **Fallback.** If S = ∅ (under the conditions of rule 2), the decoys are treated as absent for the
   tier choice. `choose_tier` continues with the no-decoy chain unchanged: `mayu` if an accepted target
   has q or PEP, else `two_peptide` if C11 has rows, else `protein_list`. In practice this is `mayu`,
   since the targets then carry q or PEP. The warning `decoy_unscored` is kept. The gene table and
   the `q_threshold` rows then equal those of the same input with its decoy rows removed. The `tier`
   rows still report the actual `n_decoys_paired`, `n_decoys_unpaired`, `decoy_paired_share` and
   `decoy_win_weight`, as they already do under `decoy_pairing_low`. `q_floor` is NA (not the picked
   tier).
4. **Order of checks in `choose_tier`.** (a) Pairing rows exist and share < cut-off → `mayu`,
   warnings exactly `["decoy_pairing_low"]` (no key check, unchanged). (b) Pairing rows exist, share ≥
   cut-off, and the target side has 0 rows → `decoy_picked`, no `decoy_unscored` (unchanged: nothing
   can be detected). (c) Otherwise S is computed: S = C_T ≠ ∅ → `decoy_picked`, no new warning;
   ∅ ≠ S ⊊ C_T → `decoy_picked` + `decoy_unscored`; S = ∅ → rule 3. (d) No pairing rows → unchanged.
5. **Scored on some rows, NA on others (one column).** Completeness is checked per column over all
   paired decoys. A single paired decoy with NA in column c takes c out of S, for every decoy and every
   target. The decoy rows are never dropped or re-ranked on their own. Dropping rows would be decoy
   loss that depends on the input rather than random loss, and ADR-0009's 1/f does not cover that. On
   the target side, a partly NA column stays in S (rule 1, NA last). Unpaired decoys never enter S
   (they do not compete).
6. **q filter (`input.q_value_max`) on decoys.** Intended behaviour: decoys with q > `q_value_max`
   are excluded before pairing and picking; NA q passes, row by row, exactly as for targets in §7.3.1.
   This is M4-12 lead decision 1 and M4-index decision 4, and it keeps ADR-0009's f as "the paired
   share of the decoys passing q ≤ q_value_max" (methods note). When q ∈ S the exclusion changes no
   gene q of a target that has a q. Such decoys rank below every accepted target with q
   (≤ `q_value_max`), so they never enter D(t) at those thresholds. It changes only f, the 80 % guard,
   the pairing counts, and the competition with accepted targets whose q is NA. **New (symmetry):** the
   filter is applied only when at least one C9 non-decoy row carries a q value. Without target q,
   targets are not filtered on q, and filtering the decoys alone would remove low-scoring decoys whose
   targets stay, which is anti-conservative at loose thresholds when the ranking is on PEP or score.
   The comment in `defaults.yaml` / §5.2 ("still used to estimate FDR") is read as "kept in C9 and
   reported by preflight". It does not mean "enter the picked competition". See question 1.
7. **Unchanged.** Pairing (prefix stripping, protein → gene map, contaminant decoys), the 80 % guard,
   1/f, `q_floor`, MAYU, the `two_peptide` and `protein_list` tiers, and every output column. No new
   config key. defaults.yaml is not touched.

## Degenerate cases

- **No paired decoy** but pairing rows exist (all unpaired): share 0 → rule 4(a), `decoy_pairing_low`.
- **No pairing rows** (no decoys, no `protein_ids`, all contaminant decoys, or all filtered): rule
  4(d), no `decoy_unscored`.
- **Zero accepted targets** with paired decoys: rule 4(b). The tier stays `decoy_picked`, every
  eligible gene is `not_detected` with q NA, and there is no `decoy_unscored` (C_T is undefined
  without rows).
- **Targets carry no key column** (every accepted q, PEP and score NA) while decoys are scored:
  C_T = ∅ → S = ∅ → rule 3 with `decoy_unscored`. The chain then gives `two_peptide` (no accepted q or
  PEP), or `protein_list` if C11 is empty.
- **A column carried only by the decoys** (e.g. decoys have `score`, targets do not): ignored, no
  warning. With S = {q}, a q tie goes to the decoy, as today (today the target's NA score ranks last,
  so the decoy also wins).
- **Paired decoys on partial (non-eligible) genes** count for completeness like every paired decoy,
  the same population f is computed on.
- **Ties become more frequent** when S is shorter (e.g. integer scores with S = {score}). They go to
  the decoy, which is conservative.
- **Empty input** (no C9): unchanged (`two_peptide` or `protein_list`).
- Row order: S and the tier depend only on column contents, never on row order.

## Examples

All with f = 1 (every decoy paired), `stats.alpha.identification` = 0.01, `q_value_max` = 0.01,
genes eligible, one accepted peptide per gene. "Main" is the current code. Numbers were computed with
the current `_picked_q`, with the columns outside S set to NA, which is the ranking of rule 1.

| # | Input | Main | Expected (this issue) |
|---|---|---|---|
| 1 | 150 genes G001–G150, targets q = 1e-5·i, PEP and score NA; 150 paired decoys (one per gene), q, PEP, score all NA | `decoy_picked`, no warning; every gene q = 1/150 = 0.006667, 150 detected | C_T = {q}, S = ∅ → tier `mayu`, warnings contain `decoy_unscored`; `table` and `q_threshold` rows equal `detect` on the same input with the 150 decoy rows removed; `n_decoys_paired` = 150, `n_decoys_unpaired` = 0; `q_floor` NA |
| 2 | targets G1–G4 q (0.001, 0.002, 0.003, 0.004), score (50, 40, 30, 20), PEP NA; paired decoys of G1, G3, G4 with score (10, 35, 25), q and PEP NA | `decoy_picked`; all 4 genes q = 0.25 | S = {score} → `decoy_picked` + `decoy_unscored`; G1, G2 q = 0.5; G3, G4 decoy wins → `not_detected`, q NA |
| 3 | 300 genes, targets q = 1e-5·i, score = 301 − i, PEP NA; paired decoys of G201–G300 with score = 301 − i + 0.5 (each beats its target), q and PEP NA | 300 detected, q = 1/300 = 0.003333 | S = {score}, `decoy_unscored`; G001–G200 q = 1/200 = 0.005 (detected); G201–G300 `not_detected`, q NA |
| 4 | 150 genes, targets q = 1e-5·i, PEP = 1e-5·i, score NA; one paired decoy of G001 with q = 1e-5, PEP NA, score NA | 150 detected (G001's target wins on PEP), q = 0.006667 | S = {q}, `decoy_unscored`; G001 q tie → decoy wins → `not_detected`; G002–G150 q = (1 + 1)/149 = 0.013423 → 0 detected |
| 5 | 400 genes, targets q = 1e-5·i, score = 1000 − i, PEP NA; paired decoys: G101 (q 5e-4, score 949.5) and G102 (q NA, score 949.4) | `decoy_picked`; G101 decoy wins, G102's NA-q decoy ranks last; 399 target wins, q = 2/399 = 0.005013 | q not complete on decoys → S = {score}, `decoy_unscored`; G101, G102 decoy wins (`not_detected`); 398 target wins, q = 3/398 = 0.007538 (detected) |
| 6 | targets carry q only (PEP, score NA); 4 paired decoys, q on 3 rows and NA on 1, PEP and score NA | `decoy_picked`; the NA-q decoy ranks last | S = ∅ → tier `mayu`, `decoy_unscored`; table equals the decoy-free run (as in example 1) |
| 7 | 151 genes, targets G001–G150 q = 1e-5·i, G151 accepted with q NA; PEP, score NA; one paired decoy of G151, q 0.005 | G151 decoy wins; 150 detected, q = 1/150 = 0.006667 | identical to main: S = C_T = {q}, no `decoy_unscored` (partial NA on the target side only) |
| 8 | 10 decoys, all q/PEP/score NA, 7 pair (share 0.7) | `mayu`, `decoy_pairing_low` | `mayu`, warnings `decoy_pairing_low` only (rule 4(a)) |
| 9 | 150 genes, no C9 target row has q; accepted targets score = 151 − i, q and PEP NA; one decoy paired to G150 with q = 0.05, score 0.5, PEP NA | decoy dropped by the q filter → no pairing rows, no accepted q or PEP → `two_peptide`, all q NA | q filter not applied (rule 6) → decoy paired, `n_decoys_paired` = 1; S = C_T = {score}, no warning → `decoy_picked`; G150 target wins (1 > 0.5); 150 genes q = 1/150 = 0.006667, detected |
| 10 | DIA-NN-shaped input: targets and decoys carry q and PEP on every row, score absent on both | — | S = C_T = {q, pep}, no `decoy_unscored`; output equals main's. M4-12 examples 1–10 and the `sim50k` M4 golden (M4-09 evidence: decoys carry q and PEP like the targets) are unchanged byte for byte |

## Acceptance criteria

- [ ] `tests/unit/genetable/test_detection_decoy_unscored.py`: examples 1–10. In example 1 and
  example 6, `table` and the `q_threshold` rows equal a second `detect` call on the decoy-free input.
- [ ] Unit tests of `shared_key_columns`: order follows `KEY_COLUMNS`; one NA decoy value drops the
  column; a column with one non-NA target value is carried; empty paired frame → every carried column
  (vacuously complete; `choose_tier` never calls it so, rule 4(a)).
- [ ] Unit tests of the degenerate cases: zero accepted targets (rule 4(b)); targets carry no key
  column (→ `two_peptide` + `decoy_unscored`); a decoy-only column (no warning).
- [ ] Property (`tests/property/`): (a) `decoy_unscored` ∈ warnings ⇔ the guard passes, the target side
  is non-empty and S ≠ C_T or S = ∅; (b) setting every value of a key column outside S to NA on both
  sides leaves `tier`, `table` and the `q_threshold` rows unchanged; (c) when every paired decoy is
  complete in every column of C_T, the output has no `decoy_unscored` and equals the output with
  the decoy-only columns set to NA; (d) row-permutation invariance of tier, warnings and table (as in
  M4-12).
- [ ] Existing M4-12 unit, property and contract tests pass unchanged. Any of them whose fixture has
  partly scored decoys and asserts on the warning tuple is updated by the test-author, not the
  implementer (role guard).
- [ ] **Calibration attack** `tests/calibration/attacks/test_detection_tiers_attack_unscored.py` on the
  `detection_tiers_calsim` `picked` input (R = 50 fast / 300 full, π_wrong ∈ {0.01, 0.05}, seeds from
  the harness, `check_fdr` at q ∈ {0.01, 0.05, 0.10}). Variant `decoys_score_only`: decoy `q_value`
  and `pep` set to NA, score kept on both sides → tier `decoy_picked` with `decoy_unscored` in every
  replicate, FDP passes. Variant `decoys_unscored`: decoy q, PEP and score NA → tier `mayu` in every
  replicate, FDP passes. The test-author records that `decoys_score_only` fails on main (teeth). Runs
  go to Slurm (`test_branch.sbatch`), not the session allocation.
- [ ] Goldens and integration: the `sim100k_real` M4 golden (giapsim M3-19 proteomics, all-NA decoys)
  changes tier `decoy_picked` → `mayu` and gains the warning. It is regenerated in the
  implementation PR, with the diff summarised in the PR body (`n_detected`, `q_floor` → NA, curve
  rows). `tests/calibration/test_gene_fdr_real_mapping.py` (asserts `decoy_picked` on the same
  generator) is resolved per question 3.
- [ ] `docs/methods/detection-tiers.md` gains a "Shared ranking key" paragraph (rules 1–6, the
  examples 3–5 shape, the q filter statement). The statistician signs it off.
- [ ] Lint, types and the fast suite pass (`not slow and not calibration_full and not e2e`).

## Questions for the lead

1. **q filter vs. the "still used to estimate FDR" comment.** §5.2 (GIAP_SPEC.md, read-only) and
   `defaults.yaml` (protected) say rows above `q_value_max` are "still used to estimate FDR". The code,
   M4-12 decision 1, M4-index decision 4 and the methods note filter the decoys before picking.
   **Proposed default:** keep the filter (rule 6). It is the right population for ADR-0009's f, and
   with q in the key it changes no target's gene q that has a q. Read the comment as "kept in C9". A
   short ADR records this reading. If Khoa wants the wording changed, the defaults.yaml comment
   becomes "peptides with q above this are not evidence; kept in C9 and in preflight's q range"
   (protected file: Khoa's approval and the ADR trailer).
2. **Partly NA decoy column: drop the column (rule 5) or keep it with NA decoys ranked best?**
   Ranking NA decoys first (+∞ score, q = 0) is the most conservative choice, but one stray NA then
   forces a decoy win. Dropping the decoy rows would be score-dependent decoy loss that 1/f does not
   cover. **Proposed default:** drop the column (rule 5), with `decoy_unscored`.
3. **Sequencing with the goldens and the real-mapping calibration.** The fix flips every giapsim
   M3-19 proteomics run (sim100k_real golden, `test_gene_fdr_real_mapping.py`, the e2e and report
   runs on sim100k) from `decoy_picked` to `mayu`, because those decoys carry no score.
   **Proposed default:** land this fix first, regenerate the sim100k_real golden to `mayu`, and switch
   the real-mapping calibration's tier assertion to "`mayu` while the simulator's decoys are unscored"
   (test-author/lead edit). The giapsim backlog item (scored null decoys) then restores
   `decoy_picked` on that path with a real D.
4. **Report which key was used?** `detection_summary.tsv` could carry three `tier` rows
   `key_q_value`, `key_pep`, `key_score` (1/0, NA outside the picked tier). **Proposed default:** no new
   rows in this issue (golden churn, and §9.3 lists no such metric). The warning is enough. Revisit
   with M9's detection card.
5. **One warning for both cases?** `decoy_unscored` covers "decoys lack a target column" and also
   "targets carry no key column at all", which is a target-side gap. **Proposed default:** one warning
   (the user's fix is the same: emit scores for both sides). The methods note spells out both cases.
6. **Target-side partial NA (rule 1)** keeps a sparse column such as q in S even when a complete
   shared column such as score exists, so NA-q targets lose to their decoys (conservative, lower power).
   **Proposed default:** keep it (behaviour unchanged from main); revisit only if a pilot input shows
   it.
7. **Preflight.** `candidate_tier` (src/giap/preflight/evidence.py) reports `decoy_picked` whenever
   any row of C9, target or decoy, has a ranking column, so all-unscored decoys still preview
   `decoy_picked`. The value is marked `provisional`. **Proposed default:** out of scope here (M4
   decides the tier); a one-line backlog follow-up if the lead wants preflight to apply
   `shared_key_columns` on C9.

## Lead decisions (2026-10-22)

1. **q filter:** the proposed default. Keep the filter; the reading is recorded in ADR-0047
   (lead-level). The defaults.yaml wording stays until Khoa chooses to change it.
2. **Partly NA decoy column:** drop the column for everyone and warn (rule 5).
3. **Sequencing:** land this fix first.
   - The implementer regenerates the sim100k_real golden with the project's update mechanism.
   - The lead edits the tier assertion in `tests/calibration/test_gene_fdr_real_mapping.py` to
     "`mayu` while giapsim's decoys are unscored".
   - The giapsim backlog item (scored null decoys) restores `decoy_picked` later.
4. **No `key_*` tier rows;** the warning is enough.
5. **One warning,** `decoy_unscored`, for both cases; the methods note explains it.
6. **A sparse target column stays in S** (main's behaviour, conservative).
7. **Preflight `candidate_tier`** is out of scope; this is a backlog item.
8. **Lead edits after the test-author step (9c0b7626):** these are existing tests the fix flips.
   - test_gene_fdr_real_mapping expects `mayu` (decision 3).
   - `test_lexicographic_score_key` case 5 expects `target_wins=True`, since a PEP present only on
     the decoy side is not a key column (rule 5).
   - `test_published_table_threshold_unestimated` becomes `test_published_table_threshold`: sim100k
     now detects genes under `mayu`, so the threshold is estimated.
   - The oracle of `test_tier_depends_on_pairing_share` switches the decoy q filter off when there
     are no targets (rule 6).
   - The implementer regenerates the `sim100k_real` M4 golden.
   - The spec's examples 1, 4 and 6 use target q = 1e-5·i, which keeps all 150 targets under
     `q_value_max`. The test-author's reading is accepted; the example table is amended to match.
