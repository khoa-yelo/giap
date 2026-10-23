# M7-13: The M7 gate — variant and novel-ORF F1, the class-FDR suite, attack 6

Milestone M7 · track C (calibration) · **pass 2 — the gate** · depends on: M7-05, M7-09, M7-10
(#239), M7-11 (#237, the pure `giap.novel.fdr`), **M7-11b** (the `novel_orfs` step wiring of
`giap.novel.fdr`; backlog, not yet specified), M7-12 (#225), the wired steps M7-04 and M7-08 ·
estimate ≈ 0 `src/` lines, ≈ 60 harness lines (`tests/calibration/_criteria.py`), ≈ 750 test
lines, 1 sbatch line · **stats** · split in three PRs (§ Split)

Refreshed 2026-10-23 against main 9da78562 (M7-10, M7-12 merged), origin/m7/novel-class-fdr
(M7-11) and origin/m7/te2-proxy-id-mixups (M7-14). The 2026-10-21 version predates those specs.

## Spec sections implemented

- §13 row M7 gate: "F1 ≥ 0.9 on simulated variants and novel ORFs; TE2 proxy passes" (TE2 is
  M7-14; its simulated twin is G2 here).
- §7.5.6 "**Calibration:** classification F1 on simulated de novo errors (substitutions, swaps,
  blocks, modifications), with variants planted between 'sample' and 'reference' genomes"; "**Gate
  M7:** simulated F1 ≥ 0.9 for `high` ALT calls; no `high` ALT calls on in-silico same-strain
  peptides except at `start_codon` or `stop_codon` sites."
- §7.6.6 "**Gate M7:** F1 ≥ 0.9 for planted novel ORFs with ≥ 2 peptides"; "**Calibration:** class
  FDR on simulations with planted novel ORFs and de novo errors".
- §7.6.4 (three regimes; a call needs class q ≤ 0.01), §7.5.4 (class error rates).
- §8.4: FDR criterion "Mean FDP over R replicates ≤ q + 2 SE(FDP)", FWER criterion "share of
  replicates with ≥ 1 false discovery is ≤ 0.05 + 2 SE", fixed seeds; suite row "Novel-ORF class
  FDR | planted ORFs + false peptides + de novo errors | FDR | 50 | 300"; the 10-minute per-PR
  budget ("set R and the genome sizes to fit"); "The full suite runs nightly on FarmShare as a
  Slurm job array, one suite per task."
- §8.5 attack 6: "Isobaric swaps everywhere → `mass_equivalent` flags present; no `high` ALT calls
  from swaps alone at q ≤ 0.01." "A failing attack blocks the merge."
- §11.6 TE2 ("in-silico proxy on same-strain peptides").
- `tolerances.yaml`: `f1_min: {variants_high: 0.90, novel_orfs: 0.90}`, `se_multiplier: 2.0`,
  `levels`, `denovo_precision_max_abs_bias`, `replicates_fast/full.variants` (20 / 100, ADR-0042),
  `replicates_fast/full.novel_orfs` (50 / 300). ADR-0015 (two-block confirmation), ADR-0042.

## Interfaces used (as merged)

| Piece | Interface |
|---|---|
| scenarios | `giapsim.scenarios.variant_f1(seed, size)`, `novel_orf_power(seed, size, n_per_class=None)`, `same_strain_insilico(seed)` (100 kb only), `isobaric_attack(seed)` (100 kb only); `SIZES = {"100k", "300k", "1M"}`; `NOVEL_PER_CLASS = {100k: 8, 300k: 16, 1M: 40}`; **`NOVEL_PER_CLASS_GATE = 300`** (M7-12 decision 26); `giapsim.simulate(cfg, seed, out)` |
| variant truth | `truth/alt_sites.tsv`: `contig, codon_start (1-based lowest base), strand, gene, ref_aa, alt_aa, n_true_peptides, frameshift_upstream`; `truth/indel_sites.tsv`: `… frameshift, junction` (`pos_reference + 1`); `truth/denovo.tsv` (error type per spectrum, M7-12 decision 6) |
| novel truth | `truth/novel_orfs.tsv`: `orf_id, class, contig, start, end, strand, aa, stop (1-based lowest base of the stop codon), true_start, n_true_peptides, rbs_planted`; `truth/novel_orfs.json` `{requested, realised, shortfall}` per class (intergenic saturates at ≈ 100–145 at 1M) |
| peptide files | `peptides.tsv` (database, six-frame wrong peptides with scored `rev_GIAPORF_*` decoys in `novel_orf_power`), `peptides_denovo.tsv`, `reference_peptides.tsv` (paired `spectrum_id`) |
| harness | `_criteria.calibrate(suite, criterion, level, fn, R, stream=…)` (criteria `level`, `uniform`, `fdr`, `fwer`; ADR-0015), `replicates(suite, full)`, `suite_seed`, `replicate_rng(stream, i)`, `log_dir()`, `replicate_range`, shards/`--merge`; `conftest.SUITE_REGISTRY`, `NIGHTLY` (fast tests > 60 s) |
| array | `calibration_full.sbatch --array=0-12` today; task i = i-th key of `replicates_full`: … 11 `novel_orfs`, 12 `denovo_precision`, **13 `variants`**, **14 `compare`** |
| GIAP outputs | `variants.tsv` (`confidence`, `alt_type`, `flags` incl. `mass_equivalent`, `start_codon`, `stop_codon`; `class_error_rate`, `class_error_n`); `novel_orfs.tsv` (`class`, `called`, `class_fdr`, `flags`, starts); after M7-11b: per-class T, D, regime in `events.jsonl` `novel_orfs.class_fdr` |

## Checks

Every gated check runs through the harness with ADR-0015 two-block confirmation and fixed seeds
(`replicate_rng(stream, i)`); the stream names below are the `stream=` argument.

| Check | Scenario, size, input | Prediction set | Truth set | Criterion | R fast / full |
|---|---|---|---|---|---|
| **G1 variant F1** (stream `variants.g1`) | `variant_f1`, 100k fast / 1M full; `--peptides peptides_denovo.tsv --reference-peptides reference_peptides.tsv --genome-relation related` | `variants.tsv` rows `confidence = high` (the ALT part of `mixed` included), `alt_type = substitution`, one prediction per (site, ALT residue), minus indel-masked rows (rule 4) | precision: all `alt_sites` rows; recall: rows with `n_true_peptides ≥ 2` and MS-observable (rule 3) | pooled F1 (rule 2) ≥ `f1_min.variants_high`; non-vacuity guard | `variants`: 20 / 100 |
| **G2 same-strain** (`variants.g2`) | `same_strain_insilico` (100k), database input, `search_db = same` | `high` ALT rows | none | zero `high` rows lacking both `start_codon` and `stop_codon`, **in every replicate** (deterministic: no error source exists) | `variants`: 20 / 100 |
| **G3 novel F1** (`novel_orfs.g3`) | `novel_orf_power(seed, "1M", n_per_class=NOVEL_PER_CLASS_GATE)`; database input (regime `picked_decoy`) | `called = true` rows of `novel_orfs.tsv`, classes intergenic, cds_outframe, antisense, nterm_extension | precision: all truth rows; recall: rows with `n_true_peptides ≥ 2` | pooled F1 over the four classes ≥ `f1_min.novel_orfs`; per-class F1 recorded | Q1: 10 / 300 |
| **G4 class FDR** (`novel_orfs.g4.<regime>.<class>`) | as G3, two regimes: database input (`picked_decoy`) and `peptides_denovo.tsv` (`denovo_subgroup`) | ORFs with class q ≤ q (`class_fdr`, float) and the class rule met | truth rows (any `n_true_peptides`) | per class and per regime, mean FDP ≤ q + 2 SE at q ∈ `levels` (0.01, 0.05, 0.10): 2 × 4 × 3 = 24 FDR checks; vacuity status (rule 6) | Q1: 10 / 300 |
| **G4-k3 clustered false ORFs** (`novel_orfs.g4k3.<class>`) | G4 de novo regime plus injected false ORFs with 3–5 wrong de novo peptides each (rule 8) | as G4 | as G4 | as G4 (M7-11 decision 3 risk) | as G4 |
| **A6a isobaric flags** (`variants.a6`) | `isobaric_attack` (100k), de novo + reference peptides, `--genome-relation same` | ALT rows (any confidence), `alt_type ∈ {substitution, block}`, every supporting read an isobaric error in `truth/denovo.tsv` | — | every such row carries `mass_equivalent`, in every replicate | `variants`: 20 / 100 |
| **A6b swaps, random** (`variants.a6b`) | as A6a | `high` rows with `alt_type = swap`, I/L-folded token (rule 9) | none (no strain) | FWER: share of replicates with ≥ 1 such row ≤ 0.05 + 2 SE | `variants`: 20 / 100 |
| **A6c swaps, systematic** (`variants.a6c`) | as A6a with the systematic transform (rule 9) | as A6b | none | as A6b; expected to fail (risk 1) | `variants`: 20 / 100 |
| attacks a, b, d, GC-poor | § Attacks | | | | |
| recorded | see rule 10 | | | | recorded only |

## Rules

1. **One run per replicate feeds every check of that replicate.** A replicate is one
   `giapsim.simulate` call seeded `int(replicate_rng(stream_base, i).integers(0, 2**31 - 1))` (the
   M7-10 pattern), then the GIAP chain through `variants` or `novel_orfs`. G3 and both G4 regimes
   share the simulation; a per-replicate count record (JSON under `log_dir()`) holds every count a
   check needs, so a check never reruns GIAP and shards can be merged.
2. **Pooled F1.** Over the R replicates of a block: P = ΣTP_P / Σ|pred| (TP_P = predictions
   matching any truth row); Rc = ΣTP_R / Σ|truth_R| (TP_R = recall-truth rows matched by a
   prediction); F1 = 2·P·Rc / (P + Rc), and F1 = 0 when Σ|pred| = 0 or Σ|truth_R| = 0. A replicate
   with no truth still contributes its predictions; a replicate with no prediction still
   contributes its truth. The check passes when F1 ≥ `f1_min.<key>` (point estimate; Q4); block 2
   runs only when block 1 fails (ADR-0015). The jackknife-over-replicates SE is recorded.
3. **Variant matching.** A prediction matches a truth row when `(contig, codon_start, strand)` are
   equal (both 1-based lowest base) and the ALT residues are equal after folding I, L and J to one
   letter. Truth rows with ref_aa and alt_aa both in {I, L} are not MS-observable and leave the
   recall denominator (they stay in the precision set). Rows with `frameshift_upstream = true` stay
   in both sets.
4. **Indel mask (G1).** A prediction whose codon contains, or lies within 3 nt of, an
   `indel_sites.tsv` `junction` is neither TP nor FP; the masked count is recorded (Q3).
5. **Novel matching.** A called ORF matches a truth row when `(contig, strand, stop)` are equal and
   the class is equal; the start is ignored (start accuracy is recorded). `nterm_extension` calls
   match by `(contig, strand, stop)` of the extended gene. A call at a truth row's key but of
   another class is an FP in its class, and the truth row stays unmatched (Q6). Two calls at one
   key: one TP, the others FP. Origin-crossing keys are compared modulo the replicon length.
6. **FDP and vacuity (G4).** FDP = false calls / max(calls, 1) per replicate (0 when nothing is
   called). A (regime, class, q) check whose mean number of target ORFs entering the class list
   is below 1/q over the block is logged with status `vacuous` and not counted as a pass: the class
   q floor `(D + 1)/T ≥ 1/T` (M7-12 decision 26) means nothing can be called. At the gate count
   this is expected only for `intergenic` at q = 0.01 (gap supply), and for de novo classes with
   few candidates. The gate record (PROGRESS.md and the calibration log in
   `docs/methods/novel-orfs.md`) states: "a class with fewer than about 1/α = 100 candidates cannot
   be called at α = 0.01" (Q9).
7. **Non-vacuity guards.** G1: Σ|truth_R| ≥ 50 and Σ|pred| ≥ 50 over the block; G3: Σ|truth_R| ≥
   100 per class except intergenic (≥ 50). A guard miss fails the check (it is a broken scenario,
   not a pass).
8. **G4-k3 injection** (test-side, no simulator change). Per replicate and per class, 20 false
   ORFs: a non-coding six-frame segment of the class's territory (the M7-12 pool), 3–5 wrong de
   novo peptides drawn from it (half of them as reversed stretches of the segment, §8.5 attack 5),
   appended to `peptides_denovo.tsv` with fresh `spectrum_id`s. The methods note's 1–5 generator
   failed at q = 0.01 (0.0185 vs bound 0.0117), so this check may fail (risk 2).
9. **A6 rows.** Swap and block tokens are compared after folding I/L/J, and a folded token reaching
   ≥ 2 distinct peptides counts as `high` even if GIAP split it (backlog line 360), so the I/L
   split cannot hide a failure; the unfolded count is recorded too. The systematic transform
   (A6c): per replicate, 30 protein sites are drawn; every de novo read covering one is rewritten
   with the same adjacent swap at that site (statistician NIT "systematic error keyed by protein
   site"). It is test-side, keyed by `replicate_rng(stream, i)`.
10. **Recorded only** (diagnostics jsonl, PROGRESS.md):
    - G1 also on the non-`mass_equivalent` subset, the high non-substitution rows (block, swap) and
      the masked count;
    - G3 per class; the `denovo_subgroup` and `unestimated` (decoys stripped, ≥ 2-peptide rule)
      regimes' F1; start accuracy (truth `true_start` ∈ {`start_longest`, `start_best_ostir`});
      the share of planted RBS starts that pass the threshold; the realised planted and
      false-candidate counts per class;
    - **class error rate, reference-proteome arm** (moved from M7-10 B2; Q8): `variant_f1` 100k
      with every reference identification whose peptide is absent from the reference genome's
      tryptic digest removed; per class the rate and n. Asserted: rate ≥ 0.95 on every class with
      `class_error_n ≥ 10` (the degeneracy of methods note "Estimand" 1, end to end);
    - G2's `start_codon` high count (at 100k the 11 non-ATG N-terminal peptides stay unmapped, so
      the exception is not exercised; M7-14's Nissle proxy exercises it: 80 sites).
11. **Tiers.** Fast tests over 60 s go into `conftest.NIGHTLY` (GitHub CI skips them; the nightly
    Slurm job runs them). The per-PR subset keeps one smoke test per test file: one replicate,
    asserting the count record is complete, no verdict. Full tier: `calibration_full` markers, array
    tasks of § Array.
12. **Failures.** A failure is fixed by an ADR on the method, never by loosening `tolerances.yaml`
    (R3). The pre-registered candidates are in § Risks.

## Attacks (tests/calibration/attacks/test_attack_m7_gate.py)

Each attack has a mutant that must fail, so the check has teeth (statistician, PR #225 rounds 1–3).

| # | Attack | Passes when | Mutant that must fail |
|---|---|---|---|
| a | G4 picked regime with the class-FDR filter removed (every candidate meeting the class rule is called) | — (mutant only) | FDP > q + 2 SE at q = 0.10 in ≥ 1 class (statistician: calling every candidate gives FDP ≈ 0.7) |
| b | matched-TIR arm (`planted_tir = MATCHED_TIR`), RBS threshold quantile 0.10 vs a mutant at 0.50 | G3 intergenic + nterm recall at 0.10 recorded | recall drops by > 2 SE at 0.50 |
| d | `nterm_extension` FDP (inside G4, locus-tag pairing, M7-12 N3); pooled novel list | per-class FDP at level (G4) | one list over all classes: FDP > bound in ≥ 1 class |
| GC | false peptides drawn preferentially from GC-poor intergenic segments (M7-11 owed attack 2), de novo regime | intergenic FDP ≤ q + 2 SE | — (recorded mean FDP) |
| 5i | reversed stretches concentrated in intergenic space (M7-11 owed attack 1) | covered by G4-k3 | — |
| 6 | A6a–c above | as above | A6c doubles as the teeth of A6b |

## Array extension (ADR-0042 backlog line 614)

- `scripts/slurm/calibration_full.sbatch`: `--array=0-12` → `--array=0-14` and the header comment
  "15 suites". Task 13 (`variants`) runs M7-10's full tests plus G1, G2, A6a–c; task 14 (`compare`)
  runs M9-22's (merged). If Q1 is accepted, `novel_orfs_gate` becomes task 15 (`--array=0-15`).
- **Cost (to be measured by the first test PR, recorded in `benchmarks/baseline.json`):** assumed
  ≤ 45 s per `variant_f1` 1M replicate and ≤ 90 s per `novel_orf_power` 1M gate replicate and
  regime. Task 13 ≈ 100 × 45 s + 3 × 100 × 5 s ≈ 1.7 h on 4 CPUs (fits `--time=02:00:00`, else
  `--shard K/2`). G3/G4 full ≈ 300 × 2 × 90 s ≈ 15 h CPU: submitted as `--array=<task> --shard K/10`
  then `--merge` (40 CPUs; sum the running `giap-*` CPUs first, R14), the job IDs in PROGRESS.md.
- The harness gains a pooled criterion (≈ 60 lines with self-tests): `check_pooled_min` /
  `calibrate_pooled(suite, stream, key, bound, records, R)` over the per-replicate count records,
  shard-mergeable, ADR-0015 blocks, statuses `pass` / `unconfirmed` / `fail` / `vacuous` (Q4).

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | G1, one replicate: 10 predictions, 8 match a truth row (any n); 9 recall-truth rows, 7 predicted | P 0.8, Rc 0.778, F1 **0.789** |
| 2 | G1 pooled: replicate A (TP_P 8 of 10, TP_R 7 of 9), replicate B (TP_P 0 of 0, TP_R 0 of 3) | P 8/10 = 0.8, Rc 7/12 = 0.583, F1 **0.675**; B's truth counts though it predicts nothing |
| 3 | truth `alt_aa = I`, prediction `alt_aa = L` at the same minus-strand codon (`codon_start` 4201) | TP (I/L/J folded) |
| 4 | truth row `ref_aa = L`, `alt_aa = I` | not in the recall denominator; a prediction there is a TP for precision |
| 5 | G1 prediction N>K, codon 2 nt from an indel `junction` | masked: neither TP nor FP; masked count 1 |
| 6 | G2 replicate: `high` M at a GTG start flagged `start_codon`; nothing else | pass; another replicate with one `high` N>D at an internal codon → G2 **fails** |
| 7 | G3: truth antisense ORF, strand −, stop 5021; call antisense, strand −, stop 5021, other start | TP; start accuracy false |
| 8 | G3: call `intergenic` at the key of a truth `cds_outframe` ORF | FP in intergenic, FN in cds_outframe |
| 9 | G4, cds_outframe, q 0.05: 120 called, 4 not in truth | FDP **0.0333** for that replicate |
| 10 | A6b, R = 20: one replicate has a `high` swap `IK>KL`; others none | share 0.05 ≤ 0.05 + 2·√(0.05·0.95/20) = 0.147 → pass; a `high` block `GG>N` flagged `mass_equivalent` is allowed |

## Degenerate cases

| Situation | Behaviour |
|---|---|
| replicate with no prediction | contributes its truth to the recall denominator |
| replicate with no truth | contributes its predictions to the precision denominator |
| Σ predictions or Σ recall-truth 0 over a block | F1 = 0 → guard (rule 7) fails the check |
| G4 class with no call | FDP 0 for that replicate; vacuity (rule 6) if T < 1/q |
| intergenic below 100 targets at q = 0.01 | `vacuous`, reported, not a pass and not a fail |
| `class_fdr = unestimated` on a class (paired share < 80 %, `decoys_unscored`, `no_makedb_pair`) | the class's G4 checks fail: the scenario provides scored makedb-tagged decoys, so `unestimated` is a bug |
| ORF across the origin (circular replicon) | key compared modulo the replicon length (none planted; M7-12 pool excludes them) |
| open-3′ ORF on a linear contig | key `stop ± 3` as M7-11 decision 13 (not expected in the circular giapsim genomes) |
| minus strand | every coordinate is the lowest 1-based base, codon and stop alike |
| multi-allelic site (two ALT residues) | one prediction per residue |
| two calls with one key | one TP, the others FP |
| `mixed` call | its ALT part counts when `high` |
| no `--reference-peptides` (any arm) | not a gate input; the class error arm requires it |

## Risks (decided before the numbers)

1. **A6c is expected to fail.** A systematic swap read by every spectrum of two overlapping
   peptides is two distinct peptides with one ALT → `high`, and swaps carry no `mass_equivalent`
   (M7-05 decision 1), and §7.5.3 says flags never change the call. ADR options (Khoa: a call-rule
   change): (i) a swap ALT is never `high`; (ii) `high` also needs `class_error_rate` ≤ a bound when
   it is estimated; (iii) restate the attack as random swaps only (A6b). Until decided, A6c is
   committed as `xfail(strict=True)` naming the ADR, as M7-11 decision 3 was.
2. **G1 may fail on `mass_equivalent` substitutions.** Statistician preview (seed 7, 100k): P 0.82,
   Rc 0.84 (F1 ≈ 0.83); 5 of 6 FPs are A>S / N>D modification errors, 1 an indel junction (rule 4
   now masks it). G1 gates the spec's literal set (all `high` substitutions, Q2) and records the
   non-`mass_equivalent` subset. If it fails: an ADR on the method (Khoa), e.g. `high` needs ≥ 2
   peptides not explained by a modification, or the gate restated on that subset.
3. **G4-k3** may fail (methods note, decision 3): ADR between a picked ORF-level competition on the
   genome decoy and a per-ORF Bonferroni over peptides.
4. **G4 de novo regime is weak**: 2–9 false candidates per class per replicate at 1M (statistician
   round 2). The FDR check stays valid; its teeth come from G4-k3 and attack a.
5. **G3 depends on the RBS threshold**; `test_sim_m7_tir_calibration.py` keeps the simulator honest.

## Acceptance criteria

- [ ] `tests/calibration/test_variants_gate.py` (G1, G2, A6a–c, the reference-proteome arm),
  `tests/calibration/test_novel_orfs_gate.py` (G3, G4, G4-k3), `attacks/test_attack_m7_gate.py`;
  fast and full markers; registered in `SUITE_REGISTRY` (`variants`, `novel_orfs` or Q1's key);
  slow fast tests in `NIGHTLY`; examples 1–10 and every degenerate row as unit tests of the
  matching and pooling helpers (`tests/calibration/_m7_gate.py`).
- [ ] `_criteria` pooled criterion with self-tests (pass, unconfirmed, fail, vacuous, shard merge).
- [ ] `calibration_full.sbatch` array extended (§ Array); `test_criteria_review.py` still finds no
  orphan.
- [ ] Full-tier results with Slurm job IDs and the recorded metrics in PROGRESS.md and in the
  calibration logs of `docs/methods/variants.md` and `docs/methods/novel-orfs.md` (§8.4), including
  the class q floor sentence (rule 6).
- [ ] The M7 gate in PROGRESS.md lists G1, G2, G3, G4 (per regime and class), G4-k3, A6a–c and TE2
  (M7-14) separately; any `unconfirmed`, `vacuous` or strict-xfail check is visible.
- [ ] Backlog lines 360 (A6 part), 614, 665–673, 681–688 (M7-13 parts) and 778 (first half) are
  closed or re-pointed.

## Split (PRs ≤ ~400 counted lines)

- **M7-13a — variants gate** (starts now): pooled criterion + helpers + G1, G2, A6a–c, the
  reference-proteome arm, the array line. ≈ 390 lines.
- **M7-13b — novel gate** (after M7-11b merges): G3, G4, G4-k3, attacks a, b, d, GC. ≈ 380 lines.
- **M7-13c — evidence**: full-tier runs, methods-note calibration logs, PROGRESS gate record.
  ≈ 60 lines of docs.

## Questions for the lead

Old questions (2026-10-21): Q1 (`replicates_*.variants`) is resolved by ADR-0042 (20 / 100); Q2
(substitutions only), Q3 (match by stop and class) and Q4 (A6 as FWER) were accepted with the
M7-index defaults on 2026-10-21 and are rules 2–5 and 9 above.

1. **Replicates for the 1M novel-ORF gate — pending Khoa (protected file; ADR-0048 draft).**
   `replicates_fast.novel_orfs` = 50 is shared with the cheap simulator checks; 50 × 2 regimes ×
   ≈ 90 s at 1M does not fit a nightly fast test, and the 300 full replicates need their own
   sharded array task. **Proposed:** new key `novel_orfs_gate: 10` (fast) / `300` (full).
   **Reversible default until Khoa answers:** the module constant
   `GATE_R = MappingProxyType({"fast": 10, "full": 300})` in `tests/calibration/_m7_gate.py`, the
   checks logged under suite `novel_orfs` (task 11); when the key lands, the constant is deleted,
   `C.replicates("novel_orfs_gate", full)` replaces it and the suite moves to task 15.
2. **G1 and `mass_equivalent`** (backlog 683). **Proposed default:** gate the spec's literal set
   (all `high` substitutions), record the non-`mass_equivalent` subset; a failure goes to a method
   ADR (risk 2), not to a narrower prediction set.
3. **G1 indel mask** (rule 4: within 3 nt of a junction). **Proposed default:** yes; indels are
   outside the truth by M7-12 decision 2, so calls they cause are neither TP nor FP.
4. **F1 criterion: point estimate.** `f1_min` has no SE term in §8.4. **Proposed default:** pooled
   F1 ≥ `f1_min` on block 1, else on block 2 (ADR-0015); the SE is recorded only. The harness
   extension is test code, written by the test-author in M7-13a.
5. **G3 pooled over classes** (§7.6.6 says "planted novel ORFs"). **Proposed default:** gate the
   pooled F1, record per class; intergenic saturation would otherwise dominate.
6. **Class mismatch counts as FP + FN** (rule 5). **Proposed default:** yes, strict; the mismatch
   count is recorded.
7. **A6c expected to fail** (risk 1). **Proposed default:** strict xfail naming the ADR; the M7 gate
   stays open on A6c until Khoa picks (i)–(iii). A6b (random swaps) gates now.
8. **Reference-proteome arm** (backlog 778). **Proposed default:** test-side transform (drop
   reference IDs absent from the reference digest), recorded plus the ≥ 0.95 degeneracy assertion;
   the "every class rate is 1" run flag stays an M7-10 backlog item.
9. **Class q floor wording** (backlog 688). **Proposed default:** the gate record and the
   methods-note calibration log carry the sentence of rule 6; a GIAP run warning when a class has
   T < 1/α belongs to M7-11b (backlog).
10. **M7-11b is not specified.** G3/G4 need `class_fdr` and the class-q call rule in
    `novel_orfs.tsv`. **Proposed default:** M7-13a proceeds; M7-13b's test-author starts once the
    M7-11b spec fixes the step outputs, so the lead should spec M7-11b next.
11. **Cost assumptions** (45 s / 90 s per replicate) are unmeasured. **Proposed default:** M7-13a's
    first commit measures one 1M replicate of each scenario on Slurm and the lead revisits R and
    sharding if either is > 2× the assumption.

## Lead decisions (2026-10-23)
1. **Q1: pending Khoa (ADR-0048).** Until Khoa decides, the module constant `GATE_R = {"fast": 10,
   "full": 300}` in `tests/calibration/_m7_gate.py` is used. The question is listed in PROGRESS.md ›
   Questions for Khoa.
2. **Q2–Q11:** the proposed defaults are accepted as written.
   - G1 gates the literal set and records the non-`mass_equivalent` subset.
   - The indel mask is 3 nt.
   - F1 is a point estimate with two-block confirmation.
   - G3 is pooled over the classes, and a class mismatch counts as a false positive plus a false
     negative.
   - A6c is a strict xfail.
   - The reference-proteome arm is a test-side transform.
   - The floor wording goes into the gate record.
   - M7-13a starts now. M7-13b waits for M7-11b, which the spec-writer is speccing.
   - M7-13a measures costs first.
3. **Order:** M7-13a test-author, then the implementer for the harness extension (pooled, shard-
   mergeable criterion), then statistician review, because the PR is labelled stats.
