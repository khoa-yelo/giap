# M2-21: Wave-B recordings, real-tool integration, local ≡ Nextflow, `nightly.sbatch` / `weekly.sbatch`, measured runtimes

Milestone M2 · depends on: M2-06, M2-07, M2-10, M2-20 (and every wrapper issue) · estimate ≈ 300 code
lines (sbatch, Python helpers, nf-test; tests not counted) · capabilities: `needs_slurm`,
`needs_nextflow`, `needs_apptainer`, `needs_conda`, `needs_bakta_db`, `needs_bakta_full`,
`needs_genomad_db`, `needs_promotech_models` — the M2 gate issue

## Spec sections implemented

- §13 M2 gate: "a contract test per tool; the full Bakta DB downloaded and the full-DB fixtures
  recorded (G17); `nightly.sbatch` and `weekly.sbatch` proposed to Khoa for `scrontab`; measured
  FarmShare runtimes per tool (with node class) replace the estimates".
- §12 `scripts/slurm/nightly.sbatch`: "`scrontab` daily at 02:00 (installed by Khoa), or submitted by
  the agent on demand | full calibration as a job array; tool integration with `-profile
  farmshare,apptainer` on the 100-kb genome (Bakta light and full DB); `record-fixtures` mode
  (`--export=MODE=record`); benchmarks pinned to one node class with `--constraint` (median of 5; fail
  on > 15 % slowdown against `baseline.json` on the same node class). Tests `main` and the newest open
  PR branches, each in its own `git worktree` under `$GIAP_ROOT/ci/`". `weekly.sbatch`: "end-to-end
  runs on K-12 and Nissle with live full-DB Bakta and with the recorded JSON; antiSMASH live run; …".
- §11.2 Integration: "the 100-kb synthetic genome through every wrapper (Bakta light DB plus recorded
  full-DB fixtures for Module 7); resume and failure paths; stub mode for fast workflow tests; local
  and Nextflow executors give identical tables"; §6.4 "an integration test checks that they produce
  byte-identical tables on the synthetic genome".
- §3.6 "every marker a job declares as required MUST be satisfied" (`--require-capabilities`,
  `ci/gates.yaml` `farmshare_nightly`/`farmshare_weekly`); §3.10 results posted with `gh issue comment`
  to "Nightly (FarmShare)"; §3.8 measured numbers into `PROGRESS.md`; §12 `containers.yml`.

## Public interface

- Wave-B recordings via `record_fixtures.sbatch TOOLS=uniop,promotech,promoter_calculator,ostir,transtermhp`
  on sim100k (and `uniop` on k12 for M6 calibration), inputs built by the GIAP builders from the
  recorded full-DB Bakta genes; flips the remaining `pending_sample` contract tests.
- `scripts/slurm/nightly.sbatch` (`MODE=test|record`; job array: `0` = integration, `1` = benchmarks,
  `2..` reserved for calibration from M3): worktrees under `$GIAP_ROOT/ci/<branch>`; runs
  `uv run pytest tests/integration tests/contract --require-capabilities <farmshare_nightly list>`,
  `nf-test test workflow/tests --profile farmshare,apptainer`, and
  `giap run … --executor nextflow -profile farmshare,apptainer` on sim100k with Bakta light and with
  the recorded full-DB JSON; posts a summary with `scripts/ci/post_summary.py` (`gh issue comment`).
- `scripts/slurm/weekly.sbatch`: live full-DB Bakta + annotation track on K-12 and Nissle, antiSMASH
  live run, TE6 on the live K-12 JSON; later milestones append their e2e steps.
- `scripts/measure_runtimes.py TRACE.tsv [--sacct JOBID]` → markdown rows `step | genome | wall | peak
  RSS | node class (job ID)` for `PROGRESS.md` › Measured numbers.
- `tests/integration/test_local_vs_nextflow.py`: same SIM100K config with every layer precomputed,
  `--executor local` vs `--executor nextflow -profile farmshare_local` → identical `annotation/*.tsv`.
- `.github/workflows/containers.yml` (§12): builds `workflow/containers/<tool>/` images, pushes
  `ghcr.io/khoa-yelo/giap-<tool>`, opens a PR updating `containers.config` digests.
- `benchmarks/test_annotation_steps.py` + `benchmarks/baseline.json` (first baseline, per node class).

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `sbatch --export=MODE=test scripts/slurm/nightly.sbatch` with all DBs present | all required capabilities satisfied; integration + contract tests pass; summary comment posted to the "Nightly (FarmShare)" issue |
| 2 | same with `local/db/genomad_db` renamed away | the job fails (not skips): `required capability needs_genomad_db is missing` |
| 3 | `MODE=record` | runs `record_fixtures.sbatch` logic for every tool; `git status` shows no fixture diff when versions are unchanged (byte-identical gz) |
| 4 | `test_local_vs_nextflow` | `cmp` equal for `operon_pairs.tsv`, `tus.tsv`, `regions.tsv` and every BED |
| 5 | `giap run` on sim100k, then kill one Promotech task and `-resume` | completes; finished tasks cached (trace shows `CACHED`) |
| 6 | a tool exits 1 (simulated with a broken env var) | `giap run` exits 4 naming the step and its `.command.log` |
| 7 | `measure_runtimes.py trace.tsv --sacct 1850000` | one row per tool, e.g. `bakta (full) | sim100k | 0:06:12 | 9.8 GB | barley (1850000)` |

## Edge cases

- Jobs named `giap-nightly`, `giap-weekly`; Nextflow `queueSize = 20`; GIAP never exceeds 128 CPUs (R14).
- Nightly checks out the newest open PR branches in their own worktrees and never pushes.
- If a §3.2 source is down, the affected recording sub-task is recorded under "Blocked" (R6) and the
  rest of the nightly continues.
- `scrontab` entries are only **proposed** (text in the PR and `PROGRESS.md`); Khoa installs them.

## Acceptance criteria

- [ ] Every tool in §7.1.1 except miniprot/BLAST+ (M3) has an active contract test on a committed
  recorded fixture (no `pending_sample` left in `tests/contract/` for annotation tools).
- [ ] `nightly.sbatch` ran green once on FarmShare (`MODE=test`) and once with `MODE=record`; job IDs in `PROGRESS.md`.
- [ ] `weekly.sbatch` ran once on K-12 and Nissle (annotation track only at M2).
- [ ] Measured runtimes per tool × genome with node class in `PROGRESS.md` replace the estimates.
- [ ] Examples 4–6 pass; proposed `scrontab` lines written for Khoa (§14.8 checkpoint list unaffected).

## Questions for the lead

1. Which node class pins the benchmarks? **Proposed default:** `wheat` (8 nodes, most available) via
   `--constraint`, if FarmShare exposes node features; otherwise `--nodelist` of one wheat node.
2. `containers.yml` needs GHCR write access for Actions. **Proposed default:** use `GITHUB_TOKEN` with
   `packages: write`; until the first images exist, the conda fallback of M2-01 Q3 stays.

## Lead decisions (2026-10-16, Khoa: "run those jobs; don't worry about uploading images; have the tools work first")
1. Q1: benchmarks pin the wheat node class (`--constraint` if FarmShare exposes the feature, else
   `--nodelist` of one wheat node); baseline per node class.
2. Q2: `.github/workflows/containers.yml` and GHCR pushes are DEFERRED (backlog); tools keep the conda /
   prebuilt-env path (M2-01 Q3 fallback). Not part of the M2-21 acceptance.
3. The existing `scripts/slurm/nightly.sbatch` (PR #105: calibration_nightly on origin/main, scrontab
   installed by Khoa 2026-10-16 01:30) is extended, not replaced: it becomes the §12 job array
   (0 = integration + contract with `--require-capabilities`, 1 = benchmarks on wheat, 2 = the
   calibration_nightly tests, 3 = real `giap run` on sim100k with Bakta light and with recorded full-DB
   JSON); the scrontab line keeps pointing at it. MODE=record runs scripts/slurm/record_tools.sbatch logic.
4. Results: one line per sub-task in local/logs/nightly.tsv plus a comment on a "Nightly (FarmShare)"
   GitHub issue via scripts/ci/post_summary.py (create the issue once if missing).
5. Testing open PR branches nightly is deferred (backlog): main only for now.
