# M7-15: Integration hardening for Modules 5 and 6 (the MVP-01 pass-2 items)

Milestone M7 · track T (integration) · **pass 2** · depends on: M7-04, M7-08 (and M7-09 for the
golden novel table) · estimate ≈ 200 code lines (+ ≈ 300 test lines, + recorded fixtures) · not
statistical

## Spec sections implemented

- §11.2 Golden layer: "gene, pair and variant tables for the 50-kb fixture"; §11.4 (the 50-kb
  golden fixture under `tests/data/golden/sim50k/`, `GOLDEN_SEED`); §11.2 Integration layer
  ("resume and failure paths … local and Nextflow executors give identical tables"); §6.3 error
  strategy ("`giap run` exits 4 naming the failed step and its log"); §9.5 ("Every integration and
  end-to-end test validates outputs with `frictionless validate`, and checks GFF3 and VCF with
  format validators"); §11.7 TE7; the MVP-01 "Hardening owed" list, applied to VARIANTS and NOVEL.

## Scope (each item a checkbox; the issue may split into two PRs at the lead's call)

1. **Golden tables.** `tests/data/golden/sim50k/` gains `pileup.tsv.gz`, `variants.tsv`,
   `variants.vcf`, `novel_orfs.tsv`, `annotation_patch.gff3` from the 50-kb fixture run with a
   strain difference and de novo peptides (M7-12 settings at 50 kb). A golden test diffs them byte
   for byte; regenerating them is a reviewed PR (§11.4).
2. **In-pipeline contract checks.** The `variants` and `novel_orfs` steps validate their TSVs against
   the Frictionless resources and their VCF/GFF3 with the M1-05 readers before returning; a failure
   is a step error (exit 4 naming the step), not a test-only check.
3. **Resume and failure.** `-resume` after a completed run re-executes nothing in VARIANTS/NOVEL;
   after changing `variants.high_min_peptides` only VARIANTS re-runs; a forced failure of
   `OSTIR_NOVEL` (stub that exits 1) makes `giap run` exit 4 naming `OSTIR_NOVEL` and its log.
4. **Real-tool arm.** One `sbatch scripts/slurm/giap_run.sbatch` under `-profile
   farmshare,apptainer` (Promotech and OSTIR from their conda envs, ADR-0008) on K-12 with the
   in-silico TE8 peptide set (cross-strain, `related`), so the second round runs on real candidates;
   runtimes per process into `benchmarks/baseline.json` via `scripts/measure_runtimes.py`.
5. **VCF FILTER.** Close the backlog item "declare GIAP's VCF FILTER values with descriptions" per
   M7-03 decision 6: PASS only; remove or keep the writer's generic `##FILTER` branch (unused) as
   the lead prefers.
6. **`novel_candidate_id`** in `mapping/peptide_hits.tsv` (M7-06 Q6), if the lead takes the
   back-fill option.
7. **Determinism.** TE7 across the two executors and two seeds-equal runs for all seven new files.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | golden run at `GOLDEN_SEED` | the five files match the committed goldens |
| 2 | `-resume` with no change | trace shows VARIANTS and NOVEL tasks `CACHED` |
| 3 | forced `OSTIR_NOVEL` failure | exit 4, message names `OSTIR_NOVEL` and `logs/` path |
| 4 | a corrupted `variants.tsv` row (test hook) | the step fails its own schema check |
| 5 | K-12 real-tool run | completes; `tables/novel_orfs.tsv` non-empty; Promotech/OSTIR second-round runtimes recorded |

## Acceptance criteria

- [ ] Items 1–7 each have a test or a recorded run (job IDs in PROGRESS.md).
- [ ] `nf-test test workflow/tests` passes with the new resume/failure tests.

## Questions for the lead

1. Golden fixture with planted strain differences: the current 50-kb golden has none. **Proposed
   default:** a second golden directory `tests/data/golden/sim50k_m7/` with the M7 settings, so the
   existing goldens are not regenerated.
2. Real-tool arm genome: K-12 with in-silico TE8 peptides (proposed) or Nissle with A2 once the pilot
   data exist (M10)? **Proposed default:** K-12 now; Nissle in M10.

## Lead decisions (2026-10-22)

1. **Golden with planted strain differences:** the proposed default. Add a second golden directory,
   `tests/data/golden/sim50k_m7/`, with the M7 settings; the existing sim50k goldens are not
   regenerated. Numeric columns that depend on BLAS (if any) follow the M5-11 decision 44/46/47
   tolerance pattern rather than byte equality.
2. **Real-tool arm genome:** K-12 with in-silico TE8 peptides now; Nissle with A2 in M10. The
   real-tool arm runs on Slurm only (never in the session), from `scripts/slurm/`.
3. **Test-author resolutions (2d5f534c) accepted:**
   - the fixture uses SNVs only;
   - `BLAS_NUMERIC` is empty;
   - example 3 uses the M2-21 `.command.log` path;
   - item 6 is the back-fill reading;
   - the writer names stay module-level, as the damage hook needs.
4. **Sequencing:**
   - Items 2, 3, 4 and 6 are implemented now.
   - The `sim50k_m7` goldens are generated only after M7-12 (#225) merges; the implementer then
     merges origin/main and runs `GIAP_UPDATE_GOLDEN=1` on Slurm.
   - The PR does not merge before that.
5. **Real-tool arm (item 4):**
   - The implementer submits the K-12 `farmshare,apptainer` run with
     `scripts/slurm/giap_run.sbatch`.
   - It records `novel_second_round_k12` in `benchmarks/baseline.json` and the job ID in
     PROGRESS.md (the lead adds that line), and sets `$GIAP_M7_15_K12_OUT` for the e2e check.
   - No genome FASTA is committed.
6. **Resume (example 2):** finding why an unchanged `--resume` re-runs 8 of 39 tasks is part of
   this issue; the cause goes in the PR body.
7. **Item 3:** each step's config slice becomes part of its task hash, so changing one step's key
   re-runs only that step and its dependants.

## Lead decisions after implementer step 3 (876ee82f)
8. **Config fixes accepted in this PR:** the farmshare 12 h time request for process_medium/high,
   VARIANTS_STEP 8 GB, and `maxRetries = 2` for 137/143 (Khoa 2026-10-23). The lead adds
   tests/unit/test_base_config_retry.py and records Slurm 1879271 in PROGRESS.md.
9. **The variants crash on X/U residues** (`minimal_alt_codon`, "no table-11 codon encodes 'X'") becomes a
   separate issue, M7-18, kept out of this PR. Reversible default for its spec:
   - an ambiguous residue (B, J, X, Z) is never an ALT residue;
   - U over a TGA codon and O over a TAG codon count as matches (selenocysteine and pyrrolysine readthrough), not ALTs;
   - any other U or O is never an ALT;
   - each skipped site is counted in one event.
10. **Goldens:** #225 has merged. The implementer merges origin/main, regenerates them with
    `GIAP_UPDATE_GOLDEN=1` on Slurm and reruns the fast suite. Then the PR opens.
11. **Review round (code review of 5c9695ce, #278).** Blocking:
    - **B1, retry growth:** a retry must change the request.
      - `time` grows with `task.attempt`: `12.h * task.attempt` on medium/high, and an explicit `2.h * task.attempt` on single/low under farmshare.
      - A 137 retry happens only while the memory request still grows under the profile's cap (`errorStrategy` closure). Otherwise the task terminates.
    - **B2, unit tests for the new pure functions:** one hand-made input per branch for `check_tsv`, `check_vcf`, `config_slice`, `config_hashes`, `candidate_ids` and `backfill_candidate_ids`, killing the seven surviving mutants listed in the review.
    - **B3, the retry-policy test:** assert on `nextflow config -flat` per profile (base, farmshare, farmshare_local) rather than a regex over base.config.
    - **B4:** the PR body explains the 8-of-39 re-run on an unchanged `--resume` (decision 6), plus items 2 and 6 and the preflight change.

    Nits taken:
    - `backfill_candidate_ids` writes every row, NA when there is no id (review 8);
    - `candidate_ids` also matches on position within the ORF span (review 9);
    - `check_tsv` honours the schema `pattern`, and `check_vcf` checks POS ≤ contig length (review 7);
    - the damage test matches the exact exception type.

    Backlogged:
    - a static guard for the step-to-section map (review 5);
    - the GIAP version in the task hash (review 6, predates this PR);
    - the preflight keep-unchanged hunk and an all-cached assertion (review 10);
    - recording the `cache 'deep'` hashing cost (review 11).

    Paused (Khoa, 2026-10-23): the round starts when agents resume. The test-author goes first (B2, B3 and the nit tests), then the implementer (B1 and the nits), then the lead (B4).
