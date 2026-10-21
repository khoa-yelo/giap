# M6 "Operons": issue index

Spec: GIAP_SPEC.md §13 row M6. Deliverables: "pair table, discordance, missed operons, flags;
RegulonDB calibration on the K-12 genome (§15.4)". **Gate:** "decoy FDR calibrated on simulations;
UniOP and TerminatorNet cut-offs recalibrated (ADR); UniOP AUROC recorded", i.e. the `missed_operons`
suite (R = 50 per PR / 300 nightly) meets FDR ≤ 0.1 + 2 SE, the `discordance` suite (R = 100 / 1,000)
stays at level, and §15.4 runs on K-12 with its ADR. **The §15.4 half of the gate is blocked** (see
"RegulonDB" below); everything else is reachable now.

Sources: §7.4 (whole module: §7.4.1 classes and columns, §7.4.2 the within-operon t, §7.4.3
discordance and the detection-limit model, §7.4.4 missed operons and the decoy FDR, §7.4.5 the flag
table, §7.4.6 the truth check, §7.4.7 outputs and tests), §7.1.3 (pairs and TUs, the distance model),
§7.1.4–§7.1.5 (`promoter_between`, `tx_rate_between`, `terminator_between` — all three merged in M2),
§7.3.4 (molar vs relative scale), §8.1 rows "Operon discordance" and "Missed operon", §8.2, §8.3,
§8.4 suites `discordance` and `missed_operons`, §8.6 sensitivity family "operon flags", §9.2 test
columns, §9.3 `pairs.tsv`, `operons.tsv`, `operon_flags.tsv`, §9.4 flag scope "Operons", §6.2 step
`operons`, §6.5 coordinates and IDs, §6.6 RNG, §11.3 simulator, §15.3–§15.4 (K-12 calibration),
§10.2 V11 and V3 (M9 consumers). ADRs read: ADR-0007 (size counts code), ADR-0008 (tool pins; UniOP
and TerminatorNet are never imported), ADR-0015 (two-block confirmation), ADR-0018 (detectability
confounding), ADR-0021 **Amendment 1** (the `uncalibrated` flag is dropped — M6 must not bring it
back), ADR-0019 (π̂₀ discipline on decoy FDR).

One issue = one PR; size counts **code** only (ADR-0007; `tests/sim/` lines are tooling and listed
separately). "Stats" = statistical change: label `stats`, the statistician reviews the methods note
and the simulations before the code (§14.1).

## How M6 is sequenced (Khoa, 2026-10-20/21)

"Write the spec as usual, but prioritize end to end integration first, then work on the rigour
after." So this index carries the **whole** of M6 — every table, test, calibration and ADR the spec
owes — and splits it into two passes, as MVP-01 does:

- **Pass 1 (integration):** 01, 04, 05 (both-detected branch), 06, 07, 08. End state: `giap run` on
  sim100k writes `tables/pairs.tsv`, `tables/operons.tsv` and `tables/operon_flags.tsv`, with the
  one-detected test emitted in its spec-sanctioned `unestimated (no_detection_limit_model)` state.
- **Pass 2 (hardening and the gate):** 02 (+ the M6-05H follow-up that turns the one-detected test
  on), 03, 09 (**the gate**), 10, and the per-issue "Hardening owed" sections.
- **Blocked:** 11 and 12 (RegulonDB). They are the §15.4 half of the gate and nothing else depends
  on them: every default they would change already exists as a documented placeholder.

Nothing is dropped; it is scheduled. Every issue that produces a table or a step states how it
reaches `giap run`: that is M6-08 for the three tables, and M6-10 for the sensitivity rows.

## How this design landed (three PRs, docs only)

At ~1,600 lines the design exceeds the ~800-line review budget, so it is split by pass. The files
are disjoint, so the three PRs are independent and all target `main`; the links above resolve once
all three have merged.

| PR | Branch | Files |
|---|---|---|
| 1/3 | `m6/design` | this index, M6-01, M6-04 |
| 2/3 | `m6/design-core` | M6-05, M6-06, M6-07, M6-08 (the rest of pass 1) |
| 3/3 | `m6/design-hardening` | M6-02, M6-03, M6-09, M6-10, M6-11, M6-12 (pass 2 and the blocked work) |

## Issues

| Order | Issue | Track | Pass | Scope | Depends on | Code lines | Stats |
|---|---|---|---|---|---|---|---|
| 1 | [M6-01](M6-01-stats-tdist-decoy-fdr.md) | S | 1 | `stats/`: Student-t MLE + two-sided p; distance-reweighted decoy FDR, monotone q, local fdr (KDE) | none | 260 | yes |
| 2 | [M6-04](M6-04-pair-table.md) | R | 1 | pair table: adjacency, OP/NO/CV/DV/SW, orientation, distances, `rna_between`, evidence columns | none | 340 | no |
| 3 | [M6-05](M6-05-discordance.md) | R | 1 (+H) | discordance: both-detected t test, BH family, `unestimated` states; **pass 2**: the one-detected P(miss) branch | 01, 04 (H: 02) | 220 | yes |
| 4 | [M6-06](M6-06-missed-operons.md) | R | 1 | missed operons: s = −\|Δ\|, CV and SW decoys, larger FDR wins, q and lfdr | 01, 04 | 200 | yes |
| 5 | [M6-07](M6-07-flags-and-operons-table.md) | R | 1 | evidence "yes" rules, the seven flags, `operons.tsv` per-TU aggregation | 05, 06 | 220 | no |
| 6 | [M6-08](M6-08-operons-step-and-run.md) | T | 1 | the `operons` step, three writers + schemas, `operons.nf`, `main.nf`, nf-test, a real `giap run` | 04–07 | 320 | no |
| 7 | [M6-02](M6-02-detection-limit-model.md) | S | 2 | detection-limit model: peptide logistic + Poisson-binomial `P_det`, replicate model, `p_miss` | none | 300 | yes |
| 8 | [M6-03](M6-03-sim-operon-scenarios.md) | F | 2 | giapsim: operon null, planted co-abundant NO pairs, planted discordance, AR(1) attack, detection truth | none | 40 (+260 sim) | yes (light) |
| 9 | [M6-09](M6-09-calibration-gate.md) | C | 2 | **the gate**: `missed_operons` FDR, `discordance` level, attack null, power recorded | 03, 05, 06 (+02) | 60 (+tests) | yes |
| 10 | [M6-10](M6-10-sensitivity-operon-flags.md) | R | 2 | sensitivity family `operon_flags`; `full` variants UniOP 0.3 / 0.7 | 07, 08 | 110 | yes (light) |
| 11 | [M6-11](M6-11-regulondb-fetch.md) | V | **blocked** | `validation/fetch_regulondb.py`, RegulonDB v12 TU and terminator truth on K-12 | network (see below) | 180 | no |
| 12 | [M6-12](M6-12-uniop-terminator-cutoffs.md) | V | **blocked** | UniOP AUROC + cut-off, d₀/s refit, TerminatorNet cut-off → `calibration.yaml` + **ADR-0022** | 11 | 200 | yes |

Total ≈ 2,410 code lines in `src/` and `validation/` (+ ≈ 260 simulator lines, ≈ 3,000 test lines);
about 1,600 of those lines are pass 1. The spec's 1.5 weeks fits if pass 1 runs while M6-02 and
M6-03 are written in parallel.

## Dependency DAG

```text
Pass 1   01 ─┬─► 05 ──┐
             └─► 06 ──┴─► 07 ─► 08 ─► (giap run writes the three tables)
         04 ─────────────┘  ┘
Pass 2   02 ─► 05H (one-detected branch) ─┐
         03 ───────────────────────────────┴─► 09  (the gate)
         07, 08 ─► 10
Blocked  11 ─► 12 (ADR-0022) ─► calibration.yaml defaults ─► re-run 09 is **not** needed
```

04 has no engine dependency and can start first; 01 is two independent engine functions and can be
split across two authors if needed. 05, 06 and 07 all need 04's frame, so 04 is the critical path.

## Waves

| Wave | Issues | Notes |
|---|---|---|
| A (now) | 01, 04 | two independent PRs; 04 is the critical path |
| B | 05, 06 | after A; both read one frame |
| C | 07, 08 | 08 is the integration milestone: a real `sbatch scripts/slurm/giap_run.sbatch` on sim100k |
| D (hardening) | 02, 03 in parallel; then 05H, 09, 10 | 09 submits `calibration_full.sbatch` tasks `discordance`, `missed_operons` (≤ 2 jobs, ≤ 32 CPUs) |
| E (blocked) | 11, 12 | only when RegulonDB is fetchable; ADR-0022 goes to Khoa |

## What M6 consumes (checked against merged main, 1189648)

| ID | Artefact | Producer (status) | What M6 relies on |
|---|---|---|---|
| G | `genetable/genes.parquet` | M4-15 + M4-13 `detectability` (merged; wired by MVP-01) | `locus_tag`, `contig`, `start`, `end`, `strand`, `partial`, `feature_type`, `tu_id`, `tu_position`, `tu_size`, `uniop_p_prev/next`, `detected`, `detection_tier`, `p_detect`, `abundance_log2`, `abundance_imputed`, `abundance_type`, `n_observable_peptides`, `gravy`, `length_aa`, `flags` |
| C1 | `replicons.parquet` | M2-04 (merged) | `contig`, `length`, `topology`, `kind` |
| C2 | `genes.parquet` | M2-04 (merged) | `aa`, intervals, `crosses_origin` — M6 rebuilds `Gene` objects for the M2 `*_between` helpers |
| C5 | `tus.parquet` | M2-12 (merged) | `tu_id`, `genes` (transcription order), coordinates — `operons.tsv` is one row per TU |
| C15 | `features.parquet` | M2-04 (merged) | rRNA / tRNA / ncRNA / tmRNA rows for `rna_between` |
| P | `pairs_raw.parquet` | M2-12 (merged) | `gene_a`, `gene_b`, `uniop_p`, `source`; **same-strand pairs only** (see decision 1) |
| L1 | `promoter_windows.parquet` | M2-13 (merged) | per-window Promotech scores for `giap.annotate.promoters.promoter_between` |
| L2 | `terminators.parquet` | M2-15 (merged) | `giap.annotate.terminators.terminator_between` (called with `yes = 0.0` so the column is the raw maximum) |
| L3 | `tx_rates.parquet` | M2-14 (merged) | `giap.annotate.promoter_calculator.tx_rate_between` |
| C9/C11/C16 | `peptides.norm.parquet`, `mapping/peptide_hits.tsv`, `work/peptide_runs.parquet` | M1-12, M3-16 (merged) | M6-02 only; `peptide_runs` is **not** in the M4-06 contract (decision 9) |
| E | `giap.stats` | merged | `bh`, `tda_qvalues` (shape reference only), `poisson_binomial_sf`, `perm_p` |
| T | calibration harness | M4-05 (merged) | `calibrate`, `check_level_pooled`, `check_fdr`; `replicates_fast/full.{discordance, missed_operons}` already exist in `tolerances.yaml` (protected, unchanged) |
| K | config | merged `OperonsConfig` (`decoys`, `distance_bins_bp`, `pi0`), `annotation.uniop.pair_threshold`, `annotation.promoters.promotech_yes`, `annotation.terminators.terminatornet_yes`, `stats.alpha.operon_flags`, `abundance.observable_len`, `detection.min_peptides_sequences_only` | no `defaults.yaml` change needed |
| Sim | giapsim | M4-07/08/09 (merged): TU layout, convergent/divergent frequencies, AR(1) TU abundance, `gene_sd`, `stoich_sd`, peptide evidence | M6-03 adds the operon scenarios |
| W | `workflow/subworkflows/local/{annotate,context}.nf`, `GIAP_STEP`, `stepIn/kv/pick` | MVP-01 (merged, #129/#130) | M6-08 copies the CONTEXT pattern; `annotate.nf` must also **emit** L1–L3 (decision 6) |

## RegulonDB: the blocked work, and what is needed from Khoa

PROGRESS.md logged `regulondb.ccg.unam.mx` unreachable on 2026-10-14 and the lead re-probed on
2026-10-21 ("connection refused"). A TLS-level probe from the session allocation on 2026-10-21 gives
a different and more actionable diagnosis:

- the host answers and `GET /` returns **HTTP 200** when the official GlobalSign intermediate
  `GlobalSign GCC R46 OV TLS CA 2025` is added to the CA bundle;
- the server certificate was issued **2026-10-14** (the day the blocker was logged) and the
  deployment **omits the intermediate** from the handshake, so every client that does not follow the
  AIA URL fails with `unable to get local issuer certificate`;
- the matching root, `GlobalSign_Root_R46.pem`, is already in FarmShare's trust store.

So this looks like a server-side incomplete chain, not a block. **This is a question for Khoa, not a
decision M6 may take:** is fetching the intermediate from GlobalSign's own CA endpoint
(`http://secure.globalsign.com/cacert/gsgccr46ovtlsca2025.crt`) and pinning it in
`scripts/` acceptable, or must RegulonDB fix its chain first? R6 is not at issue — the data would
still come from RegulonDB itself and no mirror is involved — but trusting a hand-fetched
intermediate is a configuration change, and only Khoa approves those. Until Khoa answers:

- M6-11 and M6-12 stay **blocked** and carry the PROGRESS.md "Blocked" entry;
- the defaults they would set stay as they are and are **documented, not flagged**:
  `calibration.yaml` keeps `operon_distance {d0: 50, s: 20, calibrated: false}` and
  `terminatornet_yes {value: 0.5, calibrated: false}`, and the manifest records those provenance
  fields. **No `uncalibrated` flag is introduced anywhere** (ADR-0021 Amendment 1).
- the M6 gate is claimed in two parts: the simulation half (M6-09) now, the §15.4 half when 11–12
  land. PROGRESS.md must say so rather than record a green gate.

Nothing else in M6 touches RegulonDB. In particular the UniOP threshold (0.5), the Promotech "yes"
threshold (0.5) and the TerminatorNet "yes" threshold (0.5) are read from config and
`calibration.yaml` at run time, so M6-12 changes data files only — no code in 01–10 changes when the
cut-offs are recalibrated, and **M6-09 does not have to be re-run** (its simulations plant their own
UniOP-like probabilities).

## What happens when an estimator does not exist (summary; details per issue)

| Estimator | Does not exist when | Output |
|---|---|---|
| intergenic distance | genes on different contigs (cannot happen: adjacency is per contig) | — |
| UniOP p | pair absent from `pairs_raw` (CV, DV, SW rows) | `uniop_p` NA; class from geometry only |
| evidence columns | the layer is missing, or the class is CV / DV | NA + run flag `layer_missing:<layer>`; no flag row |
| within-operon t | < `MIN_OP_PAIRS` (20) OP pairs with both genes detected and neither imputed, or zero scale | every OP test `unestimated`; `flags` gains `too_few_op_pairs`; no BH family |
| both-detected p | either gene undetected or without abundance | that pair takes the one-detected branch, or `test = none` |
| one-detected p | no peptides and one run (§7.4.3), or the logistic is separated / constant | `test = none`, `flags` gains `no_detection_limit_model` |
| any abundance test | `abundance_scale = relative` (§7.3.4) | `test = none`, `flags` gains `relative_scale`, `operon_flags.tsv` empty of abundance-driven rows |
| missed-operon FDR | < `MIN_NO_PAIRS` (20) scored NO pairs, or no decoy pair at all | no `missed_operon_decoy` rows; `flags` gains `too_few_no_pairs` |
| lfdr | KDE undefined (< 5 decoy or < 5 NO values) | `lfdr` NA; q-values still reported |
| `operons.tsv` row values | a TU with no detected gene | `detected_frac` 0, `median_abundance_log2` NA, `max_abs_delta_log2` NA, `n_flags` 0 |

No silent fallbacks (the M4-17 rule): every NA carries a flag, and no NA p enters a BH family.

## Calibration suites owned by M6 (replicates from `tolerances.yaml`, unchanged)

| Suite (key) | Issue | Generator (M6-03 scenario) | Criterion | R fast / full |
|---|---|---|---|---|
| `discordance` | M6-09 | `operon_null` — Δ of OP pairs drawn from the fitted t, no internal promoters or terminators | level, pooled over α ∈ {0.01, 0.05, 0.10} | 100 / 1,000 |
| `missed_operons` | M6-09 (**gate**) | `missed_operon_power` — independent NO pairs plus planted co-abundant NO pairs | FDR ≤ 0.1 + 2 SE | 50 / 300 |
| attack (new, A) | M6-09 | `operon_attack_ar1(ρ = 0.9)` — autocorrelated abundance along the chromosome, no operon effect | FDR | as `missed_operons` |
| attack (new, B) | M6-09 | detection confounded with distance (short intergenic gaps in a prophage) | FDR / level | as above |
| power (recorded) | M6-09 | planted co-abundance at \|Δ\| ∈ {0.25, 0.5, 1.0} SD × 20/50/100 planted pairs | none (recorded) | 50 per cell, Slurm |

Attack A is §8.5 attack 1 applied to Module 4 and is the statistically dangerous one: adjacent
same-strand genes share a local abundance level, while convergent decoys do too — the distance
reweighting is the only thing separating them. M6-09 must report whether it holds at ρ = 0.9.

## Methods notes (statistician sign-off, §14.1, §14.9)

New: `docs/methods/operons.md` (pair classes and orientation, the within-operon t, the two
discordance tests and why mixing a t p-value with a P(miss) in one BH family is conservative at best
— M6-05 Q3, the decoy FDR and its reweighting, the lfdr, the flag table and what a missing promoter
or terminator call does *not* prove), `docs/methods/detection-limit.md` (M6-02: the peptide model,
the Poisson-binomial, the replicate model, the oriented r sample). Updates: `multiple-testing.md`
(two new families), `sensitivity.md` §"Operon flags" (M6-10), `calibration.md` (two new suites),
`target-decoy.md` §"Convergent and sandwich decoys" (M6-01).

## ADRs M6 owes

- **ADR-0022 "UniOP and TerminatorNet cut-offs from RegulonDB v12"** (M6-12, **blocked**): the AUROC,
  the chosen thresholds, the refitted d₀ and s, the statement that RegulonDB's incompleteness makes
  the terminator precision a lower bound (§15.4 step 4), and the provenance fields written into
  `calibration.yaml`. Protected file; goes to Khoa.
- A second ADR is needed **only if** a calibration in M6-09 fails: a failure is fixed by an ADR on
  the method, never by loosening `tolerances.yaml` (R3).

## Backlog items folded in (docs/backlog.md, not duplicated)

| Backlog line | Where it lands |
|---|---|
| "M6 operons; RegulonDB calibration (source currently unreachable)" | this index + M6-11/M6-12; the checkbox closes when both merge |
| "M6: the simulator plants stop-codon readthrough (C-terminal extensions, ADR-0013) for `cterm_extension` calibration" | **not M6** in substance — it is an M3/M7 calibration need parked under the M6 heading. M6-03 does not implement it; the index records that it should move to M7's backlog line (lead to confirm, Q6) |
| "MVP-01 pass 2 (#129 review): settle publishDir discipline … with the M6 writers" | M6-08 decision 4: `tables/*.tsv` is published, work Parquet is not |
| "M2: move the step helpers (`_genes`, `_path`, `_replicons`) from `giap.steps.operons` to a public shared module before more step modules import them" | M6-08 does import them, so M6-08 performs the move (one small refactor PR-internal, named in its scope) |
| "M2: one shared IUPAC-aware reverse complement (… `operons.py` …)" | not touched by M6; stays in the backlog |
| "M2-15: test the manifest `uncalibrated` flag …" | **superseded**: ADR-0021 Amendment 1 drops the flag. M6-08 records cut-off provenance in the manifest instead and the backlog line is rewritten, not deleted (Q5) |
| "M4-15 decision 10: step-level I/O convention (steps read Parquet directly)" | M6-08 follows the merged convention and does not change it |

## Risks

- **Decoys are not exchangeable with targets.** Convergent pairs share no promoter but do share a
  local abundance level, and transcriptional interference can make them *less* alike than NO pairs;
  the spec's answer is to take the larger of the CV and SW FDRs. If attack A inflates the FDR anyway,
  the fix is an ADR (a π̂₀ estimate, or a third decoy), not a tolerance change.
- **P(miss) is not a p-value.** It is a posterior-ish probability put into a BH family with genuine
  t p-values (§7.4.3). It cannot be uniform under any null. The methods note must say so and M6-09
  must measure the one-detected branch's level separately; if it is anti-conservative, the ADR
  options are to BH the two kinds separately or to report P(miss) without a q.
- **Small genomes.** sim100k has ≈ 100 genes, so ≈ 40 OP pairs and far fewer with both genes
  detected: the t fit will often be `unestimated` on the integration fixture. Pass 1 must therefore
  accept an empty `operon_flags.tsv` as a *passing* end-to-end result, and M6-08's acceptance
  criterion is "the three tables exist with the right schema", not "flags are non-empty".
- **Evidence columns depend on three optional layers.** A run without Promotech or TerminatorNet
  produces no flags at all. That is correct behaviour, and the report text must not read as "no
  internal promoters exist".
- **`uncalibrated` must not come back.** Four places invite it (the distance model, the UniOP
  threshold, the TerminatorNet cut-off, the t fit). All four use provenance fields and explicit
  `unestimated` states instead.

## Feasibility on FarmShare

No external tool, database or network access is needed by 01–10; memory < 2 GB; the pair table is
O(n_genes) rows and the heaviest loop is the per-pair evidence lookup. The full calibration
(1,000 + 300 replicates plus the attacks and the 9-cell power grid) fits one `giap-calibration`
array within R14. Only 11–12 need the network.

## Decisions taken in these specs (for the lead to confirm)

1. `giap.annotate.operons.adjacent_pairs` gains a `same_strand: bool = True` parameter so adjacency
   and origin wrapping have **one** definition; M6-04 calls it with `same_strand=False` (M6-04 Q1).
2. The pair chain contains every C2 gene (cds, sorf, pseudogene), including `partial` ones; `partial`
   genes keep their geometry but are never eligible for a test (M6-04 Q2).
3. `pair_id = PAIR_<n>`, ordered by contig (C1 order), then `start_a`, then `start_b`, then class —
   an ordering that does not depend on the UniOP threshold, so IDs survive the M6-10 reruns
   (M6-04 Q3).
4. Evidence columns are computed for OP, NO and SW pairs only; CV and DV get NA (M6-04 Q4).
5. `terminator_between` is called with `yes = 0.0` so `terminator_between` is the raw maximum
   confidence, and the "yes" decision applies `annotation.terminators.terminatornet_yes` in M6-07
   (M6-04 Q5).
6. `annotate.nf` emits the three layer channels so the `operons` subworkflow can read them; the
   alternative (teaching the merged `context_layers` step to write a `pair_evidence.parquet`) is
   rejected because it re-records M2 fixtures (M6-08 Q1).
7. New `flags` values in the Operons scope, added through the schema: `too_few_op_pairs`,
   `too_few_no_pairs`, `no_detection_limit_model`, `unestimated`. `relative_scale` and `rna_between`
   already exist in §9.4 (M6-07 Q2).
8. The within-operon t needs ≥ 20 usable OP pairs; the missed-operon FDR needs ≥ 20 scored NO pairs
   (M6-05 Q1, M6-06 Q1).
9. `work/peptide_runs.parquet` becomes contract table **C16** (optional), so M6-02's replicate model
   has a validated input (M6-02 Q1).
10. RNG streams: `operons.detection_limit` (the only stochastic part — fold assignment and any
    bootstrap). The t fit, the decoy FDR and the lfdr are deterministic (M6-02, M6-06).
11. The M6 gate is claimed in two parts, and PROGRESS.md records the §15.4 half as outstanding until
    M6-12 merges.

## Questions that need Khoa (not the lead)

1. **RegulonDB TLS** (above): may the GlobalSign intermediate be fetched from GlobalSign and pinned,
   or should GIAP wait for RegulonDB to publish a complete chain? Nothing else unblocks M6-11.
2. ADR-0022 is a protected-file change and needs Khoa's acceptance once M6-12 produces numbers.

## Lead decisions on the designer's questions (2026-10-21)
1. **M6-05 Q3 — a t p-value and a `P(miss)` probability in one BH family.** The designer is right that
   `P(miss)` cannot be uniform under any null, so BH over the mixture is not justified. Do **not** ship
   that family: in pass 1 the one-detected branch is already `unestimated (no_detection_limit_model)`, so
   nothing mixes. When M6-02 turns it on, `P(miss)` is reported **with no q-value** and stays out of the
   BH family, unless M6-09 attack C shows the mixture is safe — in which case an ADR records why. This
   reverses the spec's reading deliberately; the spec predates the observation.
2. **M6-04 Q1 — extend `giap.annotate.operons.adjacent_pairs` with `same_strand=False`.** Accepted. One
   definition of adjacency and origin wrapping beats two that can drift apart.
3. **M6-08 Q1 — `annotate.nf` emits `promoter_windows`, `terminators`, `tx_rates`.** Accepted; re-recording
   the M2 fixtures to teach `context_layers` a new output costs more and buys nothing.
4. **M6-07 Q1 — keep flag rows 3 and 5 when there is no terminator layer**, with `layer_missing:terminators`
   in the row's flags. Suppressing them would silently lose every `novel_operon_candidate` on runs without
   TerminatorNet; a visible flag beats a missing row.
5. **M6-02 Q1 — add `work/peptide_runs.parquet` as optional C16.** Accepted.
6. **New flag values accepted** (`too_few_op_pairs`, `too_few_no_pairs`, `no_detection_limit_model`,
   `unestimated`), added through the schema. No `uncalibrated` anywhere.
7. **Backlog hygiene:** the stop-codon readthrough / `cterm_extension` line is M7 work; re-file it there.

## RegulonDB is unblocked (2026-10-21)
Not an outage: the server omits its TLS intermediate. The GlobalSign intermediate verifies against
`GlobalSign Root R46`, already in FarmShare's trust store, and with the chain completed the site and
`/datasets` return HTTP 200 (verified in-session). **Khoa approved completing the chain.** So **M6-11 and
M6-12 are no longer blocked** — `validation/fetch_regulondb.py` completes the chain itself and pins the
intermediate (sha256 d160e2de…63f9, valid to 2029-06-23); never disable verification. Note for ADR-0022:
RegulonDB is incomplete, so the terminator precision is a lower bound (spec §15.4).
