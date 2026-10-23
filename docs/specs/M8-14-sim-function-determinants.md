# M8-14: giapsim — function annotations, ORA and GSEA nulls and attacks, planted determinants

Milestone M8 · track F (simulator; tooling written by the implementer, tested by the test-author) ·
**pass 2** · depends on: nothing (merged giapsim M4-07/08/09, M5-04) · estimate 0 code lines in
`src/` + ≈ 330 lines in `tests/sim/` (+ ≈ 250 test lines) · statistical (light): the statistician
signs off the generators M8-15 gates on

## Spec sections implemented

- §11.3 "Determinants. Planted additive effects of promoter motif score, RBS strength and CAI with
  known coefficients, so the true R² is known."
- §8.4 generator column: "ORA | detection ∼ p_detect only", "GSEA | abundance independent of sets;
  operon-packed sets", "Determinants null | response independent of features, TU structure kept".
- §8.5 attack 2 ("Detectability confounded with region membership … → region tests and ORA, with
  p_detect weighting") and attack 3 ("Gene sets packed into one operon → GSEA after TU collapse");
  §7.7.5 "the attack null 'detectability confounded with category membership' (membrane-like genes
  with low p_detect) stays at level with weights and fails without them".
- §11.4 "Until then, Module 7 tests use the light-DB fixture plus a synthetic cross-reference
  table" (the full-DB fixtures exist now; the synthetic table is still what calibration needs).

## What the simulator must add

`tests/sim/giapsim/function_sim.py` and `tests/sim/giapsim/scenarios.py`
```python
def synthetic_ontology(n_terms: int, depth: int, rng) -> tuple[str, pd.DataFrame]: ...
        # an OBO text (GO-shaped IDs in a reserved range GO:9xxxxxx, is_a + part_of, 3 namespaces)
        # and a gene -> term table
def ora_null(n_genes: int) -> Pair: ...
def ora_attack_confounded(n_genes: int, member_shift: float) -> Pair: ...
def gsea_null(n_genes: int) -> Pair: ...
def gsea_attack_operon_packed(n_genes: int, n_sets: int) -> Pair: ...
def determinants_null(n_genes: int) -> Pair: ...
def determinants_planted(n_genes: int, coefs: Mapping[str, float], noise_sd: float) -> Pair: ...
```

1. **Synthetic ontology + annotations.** IDs `GO:9000000+` (outside real GO, so no simulated
   result can be mistaken for biology), a random DAG with `is_a` and `part_of`, 30–60 % of genes
   annotated; COG letters and two-to-twenty synthetic "pathways" written in the KEGG link/name
   formats (`ko:K9xxxx`, `path:map9xxxx`) — synthetic, never real KEGG. The xrefs land in the
   layout's C3 table, so the gene table's function columns fill through the merged M4-15 path.
2. **`ora_null`.** Detection ~ Bernoulli(p_true) with p_true from the M4-09 logistic on the
   physicochemical covariates only; categories drawn **independently** of detection and of the
   covariates. Truth: no category is enriched.
3. **`ora_attack_confounded`.** One or more categories whose members are "membrane-like": high
   `tm_segments`/GRAVY, hence low p_true (shift `member_shift` on the logit); detection still
   depends only on p_true. Truth: no enrichment *beyond detectability* — the unweighted test must
   reject, the weighted one must not.
4. **`gsea_null`.** Abundance from the merged M4-08 model (TU level + within-TU noise) with sets
   independent of abundance; set sizes spread over [10, 500] TUs.
5. **`gsea_attack_operon_packed`.** Sets each made of the genes of one or a few large TUs (packed
   operons, e.g. 10–30 genes per TU) where the TU has a large random abundance effect: gene-level
   ranking sees 30 correlated genes, TU-level ranking sees 1–3 units. Truth: no set effect.
6. **`determinants_null`.** The response is drawn from the TU structure only — TU random effect +
   within-TU noise, independent of every feature — keeping TU sizes and the detection pattern.
7. **`determinants_planted`.** log2 abundance = Σ_j β_j · z(feature_j) for promoter motif score
   (transcription block), RBS strength (translation) and CAI (codon), plus TU effect and noise;
   the true R² of the features is computed and recorded in `truth/determinants.tsv` (β, block,
   true R², true dominant block per gene). The features come from the simulator's own layers (the
   σ70 motifs, AGGAGG RBS and codon model already planted in §11.3).
8. **No `src/` change.** The truth-file schemas and the reserved-ID constants live in `tests/sim`;
   M8-03's tests import the constants from there.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| `member_shift = 0` | `ora_attack_confounded` reduces to `ora_null` (tested) |
| ontology deeper than the gene count supports | terms with 0 genes are kept in the OBO and dropped by M8-03 |
| `coefs` all 0 | `determinants_planted` equals `determinants_null` in distribution (true R² 0) |
| n_genes < 500 | allowed; M8-15's fast tier uses 400 for determinants (decision 23) |

## Examples

| # | Scenario | Expected property (tested here, not the method) |
|---|---|---|
| 1 | `ora_null(4000)` | detection independent of category membership given p_true (χ² over 50 seeds at level) |
| 2 | `ora_attack_confounded(4000, −1.5)` | members' mean p_true ≤ half the non-members'; detection still depends only on p_true |
| 3 | `gsea_null(4000)` | Spearman(set membership, abundance) ≈ 0 within MC error over seeds |
| 4 | `gsea_attack_operon_packed(4000, 20)` | every set's genes lie in ≤ 3 TUs |
| 5 | `determinants_planted(3000, {promoter: 0.5, rbs: 0.4, cai: 0.3}, 1.0)` | recorded true R² equals the analytic value from the β and feature variances to 10⁻⁶ |
| 6 | `determinants_null(3000)` | Spearman of every feature with the response ≈ 0 over seeds; TU sizes identical to the base layout |
| 7 | same seed twice | identical outputs; < 1 s for 4,000 genes |

## Acceptance criteria

- [ ] `tests/sim/tests/test_function_sim.py` covers examples 1–7 and the degenerate rows.
- [ ] The synthetic OBO parses with goatools and with M8-03's reader; IDs are all ≥ GO:9000000.
- [ ] No real GO/KEGG/COG content in any simulator file (a test greps for the reserved prefixes).
- [ ] Statistician sign-off on the generators in the PR (§14.1).

## Hardening owed

- A 1,000-gene "rich" fixture (high detection, real-shaped annotations) for golden M8 tables
  (M8-06 and M8-10 hardening items), generated from these scenarios.

## Questions for the lead

1. Reserved ID ranges `GO:9xxxxxx`, `K9xxxx`, `map9xxxx` for synthetic annotations (proposed)?
2. `determinants_planted` features from the simulator's planted signals (proposed) or drawn
   directly as Gaussians (simpler, but does not test the feature pipeline)? **Proposed default:**
   planted signals for the recovery run; Gaussian features only in unit tests.

## Lead decisions (2026-10-23)
1. **Question 1 accepted:** synthetic annotations use the reserved ranges `GO:9xxxxxx`, `K9xxxx`
   and `map9xxxx`. A test asserts that no synthetic ID falls outside them and that none collides
   with a real ID in the test GO/KEGG fixtures. No KEGG-derived content is committed.
2. **Question 2 accepted:** `determinants_planted` features come from the simulator's planted
   signals in the recovery run; Gaussian features appear only in unit tests.
3. **Test-author (83d4cb42), interface and gaps 1–10: the test-author's readings are accepted.**
   - Gap 1: `function_sim.draw(pair, rng)` returns the world (layout, expression, OBO, KEGG links,
     sets, features, truth).
   - Gap 2: `synthetic_ontology(..., genes=)`.
   - Gap 3: the synthetic namespace roots stay in the term table. Their sets cover 30–60 % of the
     genes, which is above `max_size`, so M8-03's size filter removes them.
   - Gap 4: five collections; a synthetic `go_slim` goes to the backlog for M8-15.
   - Gap 5: the per-gene dominant block goes in `truth/determinants_genes.tsv`.
   - Gap 6: `determinants_null(noise_sd=1.0)`.
   - Gap 7: confounded members are drawn at random, not by GRAVY (spec wording corrected here).
   - Gap 8: planted features get light tests only.
   - Gap 9: COG keeps its letters; COG IDs are not KEGG-derived.
   - Gap 10: the simulator never emits `GO:9999999`, the absent-term marker of the fixtures.
4. **Statistician review (d2338392, Slurm 1872005): BLOCKED on B1–B3. All three fixes accepted.**
   Spec item 5 and example 4 are amended as below.
   - **B1, attack 2 (`ora_attack_confounded`):** the member shift is carried through a covariate
     the M4-13 detectability model sees, instead of directly on the detection logit.
     - The default covariate is GRAVY. The implementer may choose another M4-13 covariate if it
       gives a cleaner effect, and records the choice.
     - Members are still drawn at random (gap 7).
     - Gate: with the fitted p_detect, the weighted ORA stays at level on the confounded set.
       With the shift on the logit alone, the old behaviour, no valid method could pass.
   - **B2, attack 3 (`gsea_attack_operon_packed`):** the extra N(0, 3²) TU effect is dropped.
     - Packed TUs draw their effect like every other TU. The shared TU level is the attack.
     - Each packed set has at least 10 TUs, inside `gsea_size`: 1–3 packed operons of 10–30
       genes plus single-gene TUs, as in test_gsea_tu_calibration.py.
     - The TU-level GSEA must stay at level, and a gene-level engine must fail (teeth).
   - **B3:** the determinants-null independence check pools the three features into one check
     per replicate (`check_level_pooled`, m = 3). No seed is changed.
   - **Example 2:** `members_mean_p_true_at_most_half` becomes a population-ratio assertion, with
     the shift applied to every gene, as the statistician suggested.
   - **Order:** the test-author updates the PR's own tests (B1–B3, example 2) and turns the
     statistician's strict-xfail arms for B1/B2 into gated tests. Then the implementer. Then the
     statistician re-checks, and the code review runs in parallel.
   - **M8-15 notes:** sets go through M8-03; the GO ancestors of the confounded set stay out of
     the null pool; the signed-argmax comparison; the gene-table column wiring; a decay-0.5 arm;
     `resid` reading `ecfg.tu_sd`. These go to the backlog for M8-15.
5. **Test-author (073c11de, 5432fd47): readings 1, 2, 4 and 5 are accepted.**
   - Reading 1: `function_sim.CONFOUND_COVARIATE` is `"gravy"` or `"n_obs"`, and its column is
     `<name>_true`.
   - Reading 4: one new pre-declared stream for the pooled check, which is not seed shopping.
   - **Reading 3: the covariate must be visible in the sequences.** M8-15 refits p_detect from
     GIAP's gene table, which computes GRAVY and observability from the proteins. So the
     members' recorded `<cov>_true` must equal the value recomputed from their `aa`. The
     implementer edits the members' sequences, e.g. hydrophobic substitutions for GRAVY or K/R
     changes for `n_obs`, rather than shifting the column alone.
     - If the logit shift of `member_shift` = −1.5 cannot be reached through a realistic
       sequence change (GRAVY within [−2, 2]), the implementer reports it. The lead then rules on
       the shift size; the tests are not loosened.
   - Next: the implementer.
6. **Implementer (69ad45f3): B2 done; B1 ruled here (the GRAVY carrier is infeasible).**
   - B2 accepted: packed sets are 1–3 chr1 operons of 10–30 genes plus `PACKED_SINGLES` = 40
     single-gene TUs (≥ 9 on small layouts), whole TUs only, no extra TU effect. Slurm 1872913:
     both calibration files at both tiers pass, apart from attack 2's fitted arm, which waits on B1.
   - **B1: GRAVY cannot carry −1.5.** With coefficient −0.6 the members need GRAVY +2.5, so 97–99 %
     land above 2 (20 seeds × 4,000 genes). A common per-member shift reachable at GRAVY ≤ 2 is
     only about −0.55, varying by seed from −0.87 to +0.02.
   - **Ruling: the carrier is `n_obs` (`CONFOUND_COVARIATE = "n_obs"`), and the size stays −1.5.**
     Observable-peptide count is the strongest real detection determinant, and K/R edits leave
     GRAVY unchanged.
     - Each member's K/R positions are edited so that its recomputed `n_obs` gives a logit shift
       as close to −1.5 as the integer allows. The **mean** member shift is −1.5 ± 0.15.
     - Members with `n_obs` = 0 cannot move. They stay members, they are counted, and the mean is
       taken over the movable members. Members are still drawn at random (gap 7).
     - Every member's recorded `n_obs_true` equals the value recomputed from its `aa`, and its
       logit equals the covariate model's logit from its recomputed covariates (exact, 1e-8).
     - The gate is unchanged: with the fitted p_detect the weighted ORA stays at level, and the
       unweighted teeth must still fail.
   - **Test amendments (test-author, this PR's own tests):**
     - `test_attack_shift_is_carried_by_the_members_covariate` drops the `layout.genes` equality.
       Between shift −1.5 and 0 only the members' `aa` (and their derived columns) may differ.
     - `test_member_shift_moves_the_members_logit_by_the_shift` changes its per-member exact −1.5
       to the logit-from-covariates identity plus the mean over movable members.
     - `test_confounding_covariate_is_recorded` and
       `test_attack_logit_is_computed_from_covariates_only[-1.5]` follow `n_obs`.
     - This is a ruling on how the shift is realised, not a widened gate (R3): the threshold and
       the at-level gate stay.
   - Minor, to the backlog: with fewer than 10 scoring TUs a packed set could, at about 1e-5 per
     set, trip the finite-p assert in `_packed_tu_rep` at a large full-tier R.
   - Next: the test-author, then the implementer, then the statistician re-check and code review.
