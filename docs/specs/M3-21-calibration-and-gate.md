# M3-21: Placement calibration suites, attack 5, recall gate and golden tables

Milestone M3 · depends on: M3-17 (whole chain), M3-19, M3-20 · creates tests only
(`tests/calibration/test_placement_fdr_sim.py`, `tests/calibration/test_denovo_precision_sim.py`,
`tests/calibration/attacks/test_attack5_reversed_stretches.py`, `tests/integration/test_m3_gate.py`,
`tests/data/golden/sim50k/mapping/`), plus `scripts/slurm/calibration_nightly.sbatch` array entries
· pure Python · estimate ≈ 80 code lines (scripts) + tests · **statistical**: evidence for
`docs/methods/placement-fdr.md` and `docs/methods/denovo-precision.md`; the statistician owns
`tests/calibration/attacks/` (§14.1) and signs the calibration results

## Spec sections implemented

- §13 Gate M3: "≥ 99.9 % recall of planted peptides; mapping FDR calibrated; TE1 and TE8 proxies
  pass" (TE proxies: M3-22). §7.2.9 "**Calibration:** placement FDR ≤ q + Monte Carlo error on
  simulated wrong peptides (database and de novo error models)"; "**Golden:** the hit tables for
  the 50-kb fixture."
- §8.4 "**FDR.** Mean FDP over R replicates ≤ q + 2 SE(FDP)." Rows: "Placement FDR | wrong
  peptides at known rate; database and de novo error models | FDR | R = 50, 100 kb | R = 300,
  1 Mb"; "De novo precision | known precision per bin | |bias| ≤ 0.05 | 50 | 300". "The per-PR
  subset MUST finish in under 10 minutes on a 4-vCPU GitHub runner …; record the timings in
  `benchmarks/baseline.json`. The full suite runs nightly on FarmShare as a Slurm job array, one
  suite per task."
- §8.5 attack 5 "De novo errors that reverse a stretch of sequence → placement FDR with shuffled
  decoys (a reversed decoy would inflate D)". "A failing attack blocks the merge."
- `tests/calibration/tolerances.yaml` (protected): `levels`, `se_multiplier: 2.0`,
  `mapping_recall_min: 0.999`, `denovo_precision_max_abs_bias: 0.05`, `replicates_fast.placement:
  50`, `replicates_full.placement: 300`, `seeds`.

## Definitions

- **False placement** (FDP numerator): an accepted peptide (placement q ≤ q) none of whose
  accepted loci is its true origin (wrong peptides have no origin; a correct peptide placed only
  elsewhere is also false). FDP = false / accepted (0 when none accepted), per stratum m and pooled.
- **Recall:** among correct, planted peptides from annotated CDS/sORFs with exactly one genomic
  locus (`n_loci_genome = 1`) and q ≤ `input.q_value_max`, the share with an accepted placement at
  the true locus with `locus_class` ∈ {`cds_inframe`, `sorf`}.

## Suites

| Suite | Generator | Replicates (fast / full) | Criterion |
|---|---|---|---|
| database errors, shuffled wrong peptides | M3-19, `pi_wrong` ∈ {0.01, 0.05} | 50 @ 100 kb / 300 @ 1 Mb | mean FDP ≤ q + 2 SE at q ∈ {0.01, 0.05, 0.10}, pooled and per stratum with ≥ 1 accepted per replicate on average |
| database errors, other-strain wrong peptides, `related` genome | M3-18 strain + M3-19 | same | same (exercises strata 1–2, π̂ = 1) |
| de novo error model | M3-20 | same | same, shuffled decoys |
| de novo precision | M3-20 | 50 / 300 | mean \|bias\| ≤ 0.05 per score bin vs `truth/denovo_precision.tsv` |
| attack 5 (statistician) | M3-20 with `reversed_stretch` weight 1.0 | 50 / 300 | placement FDR stays ≤ q + 2 SE with shuffled decoys; the test also shows that reversed decoys inflate D (documentation, not a gate) |
| recall gate | M3-19, `pi_wrong = 0.01`, 100 kb, 5 seeds | — | recall ≥ 0.999 (`mapping_recall_min`) |
| golden | 50-kb fixture, `GOLDEN_SEED` | — | the four TSVs equal `tests/data/golden/sim50k/mapping/` |

All seeds derive from `seeds.calibration_base`; each replicate uses `base + r`.

## Acceptance criteria

- [ ] The fast placement, de novo precision and attack-5 suites pass under `-m calibration` and
  together take < 6 minutes on 4 CPUs (leaving room for other milestones' suites in the 10-minute
  budget); timings in `benchmarks/baseline.json`.
- [ ] `-m calibration_full` versions run as tasks of the nightly Slurm array
  (`scripts/slurm/calibration_nightly.sbatch`, job name `giap-cal-placement` etc., R14 caps).
- [ ] Recall gate and golden tests pass; failures print the missed peptides with their truth
  loci (first 20).
- [ ] The statistician's sign-off (PR review) cites the suite outputs (FDP table per q, per
  stratum) and is linked from both methods notes.
- [ ] `PROGRESS.md` records the gate evidence (test names, FDP, recall, PR links).

## Questions for the lead

1. FDP when a correct peptide is accepted at its true locus **and** at a wrong one (multi-locus):
   **Proposed default:** not a false placement (uniqueness handles it; the peptide becomes
   `shared`/`multi_identical`).
2. Strata with very few accepted peptides per replicate (m = 2 at 100 kb): **Proposed default:**
   assert per-stratum FDR only where the mean accepted count per replicate is ≥ 1; always assert
   the pooled FDR.

## Lead decisions (2026-10-20, after the statistician's attack check)
1. Q1 accepted: a peptide is a false placement only when none of its accepted loci overlaps a true
   origin; a correct peptide also accepted at an equal-cost wrong locus is not false (uniqueness makes
   it `shared`/`multi_identical`).
2. True origin per row type: correct database peptide → all its truth `locus_tag` copies; de novo row →
   the loci of its spectrum's `true_peptide`, whatever the error type; other-strain peptide → the lifted
   ortholog locus in the sample; shuffled/random peptide → none. "Overlap" = same contig, strand and
   frame, covering ≥ 50 % of the true interval.
3. Anti-gaming guard: assert no accepted peptide with a non-origin accepted locus is labelled `unique`
   unless all its loci lie in its true gene; report the locus-level false share as a diagnostic.
   Include the mutant "placement adds all loci of cost ≤ m_T + 1" and show this check kills it.
4. Q2 rejected: assert the FDR for every list (decoy type × m) the suite exercises, regardless of the
   accepted count, plus the pooled FDR, under the existing two-block confirmation. "Exercised" is decided
   from the truth (mean ≥ 1 wrong or non-origin peptide with a target hit in that stratum per
   replicate); a claimed stratum that is not exercised is a suite-configuration error, not a skip.
5. Pooled FDP: over the union of accepted peptides per replicate, each peptide once (its best-cost list).
6. ADR-0010's PEP–length correlation: giapsim has Cov = 0 by construction; record in the methods note
   that the measurement moves to A2/M3-24 (no new simulator variant here).
7. PENDING (needs an ADR from Khoa): the related-genome suite with `search_db = other` and correct
   peptides drawn from the database strain underestimates stratum-0 false placements (statistician B2).
   Until the ADR is decided, write that case as a strict `xfail` with reason "ADR pending: stratum-0 π̂
   for search_db=other"; the rest of M3-21 proceeds.

## Lead decisions after ADR-0019 (accepted by Khoa, 2026-10-20)
8. Decision 7 is replaced: the genome-level related suite keeps a normal FDP gate (it passes under the
   current rule, ≤ 0.0005 at 1 Mb, so a strict xfail would XPASS). Add the estimate check the statistician
   specified (related suite, `search_db = other`, correct peptides from the M3-18 reference strain's lifted
   proteins via test-local construction as in the statistician's b2sim.py `--mode strain`; per replicate and
   decoy type at m = 0, assert mean(Ŵ0 = fdr_stratum × targets) ≥ mean(V0) − 2 SE at 1 Mb with SNV ≥ 2/kb),
   `calibration_full`, strict xfail "ADR-0019 implementation pending".
9. `test_database_stratum0_pi_hat_covers_the_wrong_share[3]` stays as the min-collapse regression under the
   same strict xfail; the [1] control stays a normal test.
10. The statistician's attack file tests/calibration/attacks/test_attack_placement_strain_span.py is
    committed by the lead (strict xfail through `placement` until ADR-0019 is implemented); its edge-case
    unit test gets the calibration marker or moves under tests/unit.
11. De novo precision gates against the realised share of correct distinct peptides per bin (ADR-0019
    decision 5).
12. ADR-0019 is implemented in a separate PR after M3-21 merges; it flips the xfails and updates the tests
    that encode π̂ = mean PEP (listed in the ADR's Consequences).
