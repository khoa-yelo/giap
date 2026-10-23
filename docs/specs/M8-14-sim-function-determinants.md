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
