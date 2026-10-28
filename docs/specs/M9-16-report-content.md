# M9-16: Report content: summary cards, methods and caveats, sensitivity, provenance, wording

Milestone M9 · **pass 2** · depends on: M9-02 · estimate ≈ 360 code lines (+ ≈ 320 test lines) ·
not statistical (cards print table estimates; nothing is computed beyond counts)

## Spec sections implemented

- §10.3 sections 1 (summary cards), 11 (sensitivity), 12 ("methods and caveats (generated from
  `docs/methods/` snippets)"), 13 ("provenance (manifest digest)"); "Summary cards: genes detected
  (by tier), mapping rates, islands, flagged operons, ALT calls, novel ORFs, model R²"; "Wording.
  'detected', not 'expressed'; every number with a statistic carries its error estimate or
  `unestimated`."
- §1 rule 6, R4, §8.6 (sensitivity columns), §7.8.2 "Caveats (printed in the report)", §7.4
  (σ70 promoters and intrinsic terminators only: "The report states this next to the table"),
  §7.7.1 (coverage; `bakta_db = light`), ADR-0017 ("Every run report states the KEGG retrieval
  date"; academic-use terms), ADR-0021 Amendment 1 (ℓ̂, ℓ′, Î₁ per replicon), ADR-0025/0026/0027
  (GSEA claims), ADR-0032/0033/0034 (determinants claims), M6-06 backlog (`bin_counts`, π₀).

## Behaviour

1. **Summary cards** (complete set; each shows the estimate and its error control, or `unestimated`):

   | Card | Number | Error estimate |
   |---|---|---|
   | Genes detected | `detected = true` count, per `detection_tier`; `detected_2pep` count beside it | tier FDR at `stats.alpha.identification` (`detection_summary` tier rows); the `two_peptide` tier says "FDR unestimated" |
   | Mapping | peptides placed / input; proteins placed / input | placement q ≤ `stats.alpha.placement` per family (0 / ≥ 1 difference, per decoy type); `mapping_fdr` stratum rows shown only as descriptive values with their full key (decision 15) |
   | Islands | count with `p_value ≤ stats.alpha.context` | "FWER ≤ `<alpha>`" |
   | Flagged operons | `operon_flags` rows with `q_value ≤ stats.alpha.operon_flags` | "FDR ≤ `<alpha>`" |
   | ALT calls | `variants.tsv` rows by `confidence` | `class_error_rate` range, or `unestimated` |
   | Novel ORFs | `called = true` | `class_fdr` max among called, or `unestimated` |
   | Model R² | elastic-net `cv_r2_mean` ± `cv_r2_sd` (ADR-0033 primary) | "null p `<p>` (uncalibrated)" while rows carry `null_uncalibrated` (ADR-0034) |

2. **Methods and caveats.** Snippets are plain HTML fragments under
   `src/gippa/report/snippets/<topic>.html` (package data, so an installed GIPPA has them). Each is
   generated from a marked block in the matching `docs/methods/<topic>.md`
   (`<!-- report:begin -->` … `<!-- report:end -->`, plain paragraphs and lists only) by
   `scripts/dev/sync_report_snippets.py`; a unit test fails when a snippet and its source differ.
   No Markdown library at run time (`markdown` is not a runtime dependency).
3. **Mandatory caveats** (constant text, each placed in its section): §7.4 σ70/intrinsic-only note
   (operons); §7.8.2 caveats from `gippa.determinants.CAVEATS` (determinants; #178); GSEA "conditional
   on detection, no p_detect adjustment" (ADR-0026), "no level claim for rows flagged
   `genomic_cluster`" (ADR-0025), "NES has no interval" (ADR-0027); determinants "univariate CIs
   uncalibrated" (ADR-0032), "elastic-net ΔR² is the primary attribution; LightGBM and SHAP are
   secondary checks" (ADR-0033), "null p uncalibrated" (ADR-0034); function "annotation coverage
   leans toward well-studied genes" with the coverage numbers (§7.7.1), "Bakta light DB" when the
   run flag `bakta_db_light` is set; KEGG "academic use only (rest.kegg.jp terms)" plus the
   retrieval date (ADR-0017); islands: the validity domain sentence of `docs/methods/islands.md`
   and a table of ℓ̂, ℓ′ (`block_len`) and Î₁ per replicon from `context/island_diagnostics.parquet`
   (ADR-0021).
4. **Operons extras:** the `operons.missed_operons` event's `bin_counts` and `pi0` as a small table
   (backlog M6-06); `-Infinity`/`Infinity` bounds printed as "−∞"/"∞".
5. **Sensitivity section:** `tables/sensitivity.tsv` rendered whole (it is small), grouped by
   `family`, `jaccard` printed `%.2f` or `NA` ("no flagged set in either run" when both n are 0).
6. **Provenance** (index decision 26): with `run_manifest.json`: GIPPA version and commit, command
   line, input SHA-256s, tool versions and digests, database versions and checksums (including the
   KEGG retrieval date and GO release), seed, and the SHA-256 of the manifest itself ("manifest
   digest"). Without it: the M9-02 pass-1 subset and the sentence "run_manifest.json not written by
   this GIPPA version".
7. **Wording check** (index decision 25): `gippa.report.checks.wording(html) -> list[str]` finds
   `express` (any case, word-start) outside quoted citations; the `report` step fails a unit test,
   not the run.
8. **Error-estimate check:** `gippa.report.checks.unlabelled_numbers(cards) -> list[str]`: every card
   value has an attached estimate string or `unestimated`.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| `variants.tsv`, `novel_orfs.tsv`, M8 tables absent | their cards absent, not zero |
| `two_peptide` tier | card shows count and "FDR unestimated (two-peptide rule)" |
| no KEGG installed (`layer_missing:kegg`) | the KEGG caveat is replaced by "KEGG not installed" |
| `island_diagnostics.parquet` absent | ℓ̂ table omitted with a note |
| no `operons.missed_operons` event (local run without OPERONS) | table omitted |
| manifest present but missing the databases block | databases row "not recorded" |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | tier `decoy_picked`, 812 detected, tier FDR at 0.01 = 0.0098 | card "812 genes detected (picked target–decoy FDR 0.0098 at q ≤ 0.01)" |
| 2 | `two_peptide` tier, 640 detected | "640 genes detected (two-peptide rule; FDR unestimated)" |
| 3 | elastic net cv_r2_mean 0.21, sd 0.03, null p 0.001 flagged `null_uncalibrated` | "Model R² 0.21 ± 0.03 (CV); null p 0.001 (uncalibrated)" |
| 4 | a snippet edited without editing its `docs/methods` source | the sync test fails |
| 5 | a template string "expressed genes" | `wording()` returns it |
| 6 | no manifest | provenance shows preflight SHA-256s and "run_manifest.json not written by this GIPPA version" |

## Acceptance criteria

- [ ] `tests/unit/report/test_cards.py`: examples 1–3 and every row of the card table, including
  absent tables.
- [ ] `tests/unit/report/test_snippets_sync.py`: example 4.
- [ ] `tests/unit/report/test_wording.py`: example 5 and a scan of all templates and plot title
  functions on the snapshot fixture.
- [ ] `tests/unit/report/test_caveats.py`: each mandatory caveat present in its section when its
  condition holds.
- [ ] `docs/methods/*.md` gain `report:` blocks (statistician reviews the wording of the methods
  blocks; ADR-mandated sentences are verbatim from the ADRs).

## Questions for the lead

1. The manifest writer (index question 2). **Proposed default:** M9 reads it when present.
2. The KEGG retrieval date lives in `db_manifest.json`, which is not in `OUT/`. Without the run
   manifest, **Proposed default:** print "KEGG retrieval date: not recorded in this run (manifest
   step pending)", never guess.

## Lead decisions (2026-10-22)

1. The run manifest exists since MVP-03a (#196, `gippa.manifest`): M9-16 reads it when present and
   says so when absent; it never recomputes provenance itself.
2. The KEGG retrieval date is printed from the manifest's database entries when recorded there;
   otherwise "KEGG retrieval date: not recorded in this run", never guessed.
3. Cards print table values only (no new statistics); every number on a card is traceable to a
   named table and column, which the tests check.

## Lead decisions on the tests (2026-10-22, after 1c0cea0)

4. **#178 dependency.** The implementer merges origin/main once #178 (M8-10) lands (queued) and
   reads `tables/determinants_model.tsv` through its `PATHS` key; it does not copy the schema. Until
   then it builds the other sections first.
5. **ALT calls card.** Only `call == "ALT"` rows count; `mixed` is not an ALT call. The estimate
   range is over the numeric `class_error_rate` rows; when some rows are `unestimated` the card
   appends "(N unestimated)", and when all are, the estimate is `unestimated`.
6. **Two-peptide wording.** Both places read "(two-peptide rule; FDR unestimated)"; the card table
   follows example 2.
7. **Digests.** Provenance prints the manifest file's SHA-256 (the tested value) and, labelled
   "provenance digest", the manifest's `provenance_digest` field when present.
8. **KEGG date.** It is printed in provenance whenever the manifest records a `kegg*` database,
   whether the function module ran or not, and in the function section as tested.
9. **Ribosomal-gene table.** Not in M9-16; moved to the backlog (M9 pass 2).
10. **Required snippets.** `islands` only, plus sync of every marked doc (as tested).

## Lead decisions after the implementation (2026-10-22, 4c7592d)

11. **Manifest in pipeline runs.** The manifest is written after the run, so the pipeline's report
    step cannot read it. `gippa run` (both executors) re-renders `report.html` once the manifest is
    written; that is a follow-up issue (M9-16b, backlog), not this PR. Until then the absent-manifest
    sentence stays as tested and `gippa report OUT` gives full provenance.
12. **checks.py add/add with M9-13 (#224).** Whichever merges second keeps both modules' functions
    and both `__all__` lists.
13. **Sensitivity tables** carry both `data-table="sensitivity"` and `data-sensitivity-family`;
    accepted. proteins_placed and multi-tier gene estimates read "unestimated"; accepted.
14. **Size.** 610 counted lines (warn > 400, fail > 800) accepted for this integration-pass PR; the
    islands `report:` wording gets a light statistician read in review.

## Lead decisions on the review (2026-10-22, PR #228)

15. **B1 Mapping card.** The estimate reads "placement q ≤ α per family (0 / ≥ 1 difference, per
    decoy type)"; the `mapping_fdr` placement rows follow as "descriptive stratum FDR (not
    q-values)", each keyed by decoy type, differences and length class. Lead fix, tests in
    tests/unit/report/test_cards_review_228.py.
16. **B2 islands wording.** The `report:` block names the supported regime (ρ ≤ 0.5 at every size
    measured, ρ = 0.9 at 4,500 genes), says there is no evidence beyond ρ 0.95 or the 20-TU block
    cap, and how to read a high Moran's I with block_len ≤ 7; snippet regenerated.
