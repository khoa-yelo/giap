# M3-25: Implement ADR-0019 (stratum-0 π̂ from exact-match shares per length class)

Milestone M3 · branch `m3/adr0019-pi-hat` · depends on: M3-21 (#113, merged) · **statistical**
(label `stats`) · the ADR is the spec: docs/adr/ADR-0019-stratum0-pi-hat-exact-share.md, decisions 1–5
and Consequences.

## Lead decisions (2026-10-20, after the test-author's report)

1. **r_ℓ is a share of n_ℓ**, not of k_ℓ: the fraction of the class's peptides with an exact decoy hit.
   This is the reading under which U/(1 − r) estimates π, and it matches the statistician's reference
   rule `_pi_adr0019`.
2. **QC rows sit beside `genome_inconsistent_rate`, they do not replace it.** The two measure different
   populations (that one: C9 `other`, non-contaminant, over accepted exact placements). New metric
   `stratum0_unmapped_share`, value k/n, `n` = n: one row with group_type `run` / group `all`, and one
   per class with group_type `length_class` and group = the class label.
3. **How the class table reaches `context_stats`:** follow M4-18's `detectability_summary` pattern. The
   placement step writes the class table as a new step output `stratum0_classes`
   (`mapping/stratum0_classes.parquet`, the columns of `stratum0_pi_classes`), and `context_stats` takes
   it as a new **optional** input. `mapping_fdr.tsv` columns stay unchanged, as the ADR says. The
   Nextflow wiring joins the existing backlog wiring item; this PR registers the step I/O only.
4. **`input_prefiltered`:** computed on P's collapsed PEP (what `placement` actually sees), pooled over
   all classes (n, k, r), not per class. It is raised through the existing run-warning mechanism, the
   same as `c_assumed`, so a WARNING log plus the run warning; no extra table.
5. **`pi_wrong_database` is no longer read for the stratum-0 π̂** (**superseded by decision 11**: the
   decision-4 warning means its own mean PEP over P, so the function is dead in `src/` and goes).
6. **π̂[3] (`null_share_stratum0`): adopt the statistician's per-class form** as the primary assertion —
   compare each class's π̂_ℓ with that class's realised null share — and keep the pooled
   (W + c)/(D + 1) comparison alongside it. (Statistician note, #113 round 1.)
7. **Contaminant rows are excluded from P**, like `same_db` exact peptides: they are not genome-derived,
   so they carry no information about the chance-hit rate. Add a test.
8. **Peptides outside every length class take π̂ = 1 and do not enter c** (ADR decision 2). Reachable only
   if `mapping.exact` admits lengths outside the classes, so a pure-function test is enough.

## Definition of done
ADR-0019 Consequences satisfied; both strict xfails flipped with their ΣV0 floors intact (≈ 79 and
≈ 275 measured under the test-author's stand-in, floors 30 and 100); methods notes
docs/methods/placement-fdr.md ("Cross-strain input and run collapse") and docs/methods/denovo-precision.md
(estimand) updated; statistician sign-off.
9. **The Nextflow out spec is in this PR, not the backlog** (amending decision 3). Registering
   `stratum0_classes` makes three static tests bind the registry to `workflow/`, so the implementer also
   edits `workflow/subworkflows/local/mapping.nf` (the PLACEMENT out spec gains
   `'stratum0_classes=stratum0_classes.parquet'`) and regenerates both
   `workflow/tests/data/mapping/step_registry.json` and `workflow/tests/data/annotate/step_registry.json`.
   What stays on the backlog is the *consumer* wiring: joining the new file into `context_stats` in the
   Nextflow graph, with the other wiring items.
10. **`is_contaminant` on the costs frame** (filled by `best_costs`) is the accepted interface for
    decision 7. Whether a contaminant peptide's own hits stay in the stratum lists is out of scope here;
    the ADR only speaks about P.
11. (code review of #124) Resolve `pi_wrong_database`: `_stratum0_pooled` computes the mean PEP itself,
    correctly, because `pi_wrong_database`'s max-q / 1.0 fallbacks would misfire. It is therefore dead in
    `src/`. **Delete it** and drop the two sentences in lead decision 5 and in
    docs/methods/placement-fdr.md that keep it "because the warning needs it".
12. (code review of #124) Add the two missing edge-case tests; both behave correctly today, both mutants
    survive, and the wrong direction is silently anti-conservative: a contaminant falling in a class
    empty of P weighs 1 (not 0), and an empty P gives c = 1 (not 0).
13. (code review of #124) Guard the degenerate configuration: `mapping.length_classes: []` now raises
    IndexError where main returned a table. Add a config validator rejecting an empty list, and make
    `stratum0_pi_classes` raise a named ValueError, not a bare KeyError, on a frame without `m_target`.
14. (code review of #124) Non-blocking, recorded only: `_stratum0_pooled` takes n/k/r from the class table
    but the mean PEP from all of P, unreachable with the default config; the test helper's `prefilter=`
    silently discards `inject_short=`; and the per-class gate has zero tolerance for a class seen in one
    replicate. Backlog, not this PR.
