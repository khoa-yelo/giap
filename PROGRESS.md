# GIAP progress

> **Resuming after a timeout:** start an allocation, run `bash scripts/resume_status.sh`, then use the
> resume prompt (GIAP_SPEC.md §0.3). Every worktree's next feature-loop step is printed; uncommitted
> files mean an agent stopped mid-step — inspect them, then re-dispatch that step (test-author →
> implementer → code-reviewer). Slurm jobs survive the session; their logs are in `local/logs/`.
_Last updated: 2026-10-15 13:48 PDT · session https://claude.ai/code/session_01LEyHKvwXcvHWvkxzhtD54J_

## Milestone progress
<!-- progress:start -->
🟩 done (merged) · 🟨 running (open branch) · ⬜ pending; designed issues: 66/79 done (84 %), 5 running · generated 2026-10-15 13:48 by `scripts/dev/progress_bars.py`

| Milestone | Progress | Done | Running | Pending |
|---|---|---|---|---|
| M0 Scaffold | 🟩 | 1 | 0 | 0 |
| M1 Inputs and preflight | 🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩 100 % | 12 | 0 | 0 |
| M2 Annotation wrappers | 🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩⬜⬜ 90 % | 19 | 0 | 2 |
| M3 Mapping and decoy FDR | 🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟨⬜⬜ 88 % | 21 | 1 | 2 |
| M4 Gene table and region tests | 🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟨🟨🟨🟨⬜⬜⬜⬜ 62 % | 13 | 4 | 4 |
| M5 Autocorrelation and islands | _not designed yet_ | - | - | - |
| M6 Operons | _not designed yet_ | - | - | - |
| M7 Variants and novel ORFs | _not designed yet_ | - | - | - |
| M8 Function and determinants | _not designed yet_ | - | - | - |
| M9 Visualisation and report | _not designed yet_ | - | - | - |
| M10 Pilot validation | _not designed yet_ | - | - | - |
| M11 Release 1.0 | _not designed yet_ | - | - | - |
<!-- progress:end -->

## Current milestone
M0: Scaffold. Status: done (PR #1, merged 2026-10-14 under ADR-0004)
M1: Inputs and preflight. Status: done (12 issues; PRs #2–#11, #13, #14)
M2–M4: designed (docs/specs/M{2,3,4}-index.md); implementation in progress (≈ 60 % of M2–M4 issues merged)

## Done (with evidence)
- [M0] `gh` authenticated as khoa-yelo (SSH protocol); repository is public (Actions minutes unlimited).
- [M0] Repository connected to git@github.com:khoa-yelo/giap.git; branch `m0/scaffold` pushed.
- [M0] Bootstrap ran (`scripts/bootstrap_farmshare.sh`): uv 0.13.1, Nextflow 26.10.0, nf-test 0.9.5,
  gh 2.102.0, OpenJDK 21, all under `local/` (ADR-0002).
- [M0] `uv.lock` resolved (131 packages); venv in `local/venv`.
- [M0] Probe in the session allocation (iron-01) and as Slurm job 1845260 (`giap-probe`, wheat-08,
  COMPLETED in 31 s): batch path proven.
- [M0] Appendix A files, ADR-0001 (imports §2.1/§2.4), CI guards (`scripts/ci/check_*.py`), `ci.yml`,
  capability-marker plugin with `--require-capabilities` (tests/unit/test_capabilities.py, 2 passed).

## In progress
Merged since 10:45: M3-17 (#77), M4-13 (#83), M2-17 (#86), M2-10 (#93), M2-19 (#91), M2-14 (#94), waitmerge
conflict stop (#90), antiSMASH compare ignores set order (#97). Run `bash scripts/resume_status.sh`.
Real `giap run` on main fa19c16 (Slurm 1849451, sim100k, real DBs; local/scratch/real-run-main): every
ANNOTATE process completes (Bakta … Promoter Calculator, UniOP, codon usage, regions, tx_layer); stops at
CONTEXT_LAYERS (`unknown step 'context_layers'` = M2-20, not built yet) → next: M2-20.
Long runs on main (all green): calibration_full 89/89 (1849366), integration 90 passed (1849367); real
per-process arrays K-12 (1849282) and Nissle (1849365): all ok after #97/#98.
Open PRs (state at 13:48):
- #98 real Promotech/UniOP inputs + recordings for K-12/Nissle: merging on green CI (its backlog entry is
  removed here).
- #92 M4-12 detection tiers: fixes in (suite 1849404, calibration 1849405 green; decision 11) → needs the
  statistician's sign-off.
- #89 M4-17 region tests: ADR-0018 (accepted-by-khoa) implemented (1186184) but A2-HL calibration still fails
  (0.144 at α 0.05, R=1000, job 1849427) → statistician to reconcile decision 10 (see spec "Open").
- #87 M4-16 detection summary: fixes in (lead pushed once with Khoa's approval; test 3/6 → 3/5) → statistician
  re-check.
- #96 MAPPING subworkflow (M3-17 follow-up): needs a code review.
Session mode: Khoa set light mode (≤ 2 agents; token-light, Slurm-heavy work first).

Follow-ups queued (docs/backlog.md has the full list):
- Per-process real-tool Slurm job array for workflow modules (inputs from recordings) — after #65 (Khoa).
- calibration_full variants for migrated suites (statistician); GO data-version pin (M2-09).
- M3-21: c_ref − c_true, chance rate of wrong peptides, PEP–length correlation (ADR-0010).
- M4-12/M4-17 forwarded items; M6: simulate readthrough (cterm_extension).
Lessons (keep):
- A PR that conflicts with main gets NO CI run; waitmerge now stops on it (#90). Merge main, then retry.
- Fixture genomes need `[topology=circular]` (spec §~770); plain NCBI FASTA → linear replicons. Nissle
  GCA_043228065.1 has 3 plasmids (p15A, pMUT1, pMUT2).
- Generated golden outputs (tests/data/golden/) are exempt from the role check, like recorded fixtures.
- Shared test files across tool-layer PRs: scratchpad merge_main_nf.sh resolves BUILDERS + step_registry.
- Real `giap run` (SIM100K, farmshare_local,apptainer, real $GIAP_DB_DIR): preflight ok, BAKTA (13.5 min),
  GENOMAD, TERMINATORNET ok, `--resume` caches; stops at the first unregistered py step (codon_usage, M2-17).
- Tests built on invented fixtures hid contract bugs (M4-13, M4-14): fixtures must come from the real
  producers' outputs. Test module basenames must be unique (scripts/ci/check_test_names.py in CI).
- Slurm "user env retrieval failed / requeued held" happens repeatedly on batch nodes: release held giap-*
  jobs (`scontrol release <id>`); a 2-min watchdog does this during sessions.
- Shared files across tool-layer PRs (tests/unit/test_nf_cmd_sync.py BUILDERS,
  workflow/tests/data/annotate/step_registry.json): merge main, keep both builders, regenerate the registry.
- Never rebase + force-push: merge origin/main into the branch and push plainly (force-push is blocked;
  Khoa, 2026-10-15). Recover a locally rebased branch by cherry-picking onto origin's branch.
- Slurm work dirs live under $GIAP_ROOT (scratch/…), never node-local /tmp or the session scratchpad.
- Rebuilding a shared micromamba env from one branch breaks `test_envs_built` on all other branches until
  that branch merges (promotech/PyQt5, #65).
- Every tool wrapper gets a real-tool run before code review; stub fixtures hid ~10 bugs in M2-07.
- Keep ≤ 128 CPUs: Slurm rounds each job up by 2 CPUs; chain extra jobs with --dependency.
- Slurm jobs run from the main checkout only; never remove a worktree a job may be using.
- Subagents never install packages; set PYTHONPATH per worktree; unset PYTHONHOME after sourcing the env.

## Next actions (ordered)
1. Merge #66, #67; finish #68, #69, #65 (+ run-wiring PR), M4-08.
2. Then: M4-09 (after M4-08), M3-17 mapping steps (after M3-16), M2-13…M2-20 tool layers (after #65,
   using the recordings), M2-10, per-process real-tool job array, first full ANNOTATE `giap run`.

## Blocked
| Since | What | Host or tool | Needed from Khoa |
|---|---|---|---|
| 2026-10-14 | RegulonDB v12 download (needed at M6) | `regulondb.ccg.unam.mx` unreachable from FarmShare (probe) | nothing yet; retry at M6, then ask |

## Answered by Khoa
- 2026-10-15: force-push blocked → merge main into feature branches instead of rebasing (option 3).
- 2026-10-15: light mode (no new agents) 01:30–03:30 for the usage reset; per-process Slurm testing for M2-07 wanted.
- 2026-10-14: ADR-0013 cterm_extension, ADR-0015 two-block confirmation, ADR-0016 size gate excludes workflow/tests (all accepted-by-khoa).
- 2026-10-14: novel ORFs are checked for a nearby promoter and terminator — evidence columns, flags and a
  ranking key, not a calling requirement (ADR-0012, accepted-by-khoa; implemented in M6).
- 2026-10-14: ADR-0011 (`low_complexity` placement flag) accepted as proposed.
- 2026-10-14: M3-13 tblastn HSP chaining limits (max gap 300 nt, query overlap ≤ 10 aa) are fixed
  constants documented in the methods note — no config keys.
- 2026-10-14: replace the slow Bakta full-DB download with a resumable one (done, #33).
- 2026-10-14: §7.3.2 picked FDR rescales decoy wins by the paired share (ADR-0009, accepted-by-khoa).
- 2026-10-14: subagent commit blocked by the permission system → Khoa approved the lead committing the M2-02 test fix.
- 2026-10-14: ADR-0005 item 4 (new preflight warning thresholds) rejected: decoy share and free RAM/disk are Info only.
- 2026-10-14: ADR-0004 (M0 size exemption; hook realpath fix) and ADR-0007 (size guard counts code only) accepted; "merge the stack".

## Questions for Khoa
| # | Question | Options | Default if no answer | Asked |
|---|---|---|---|---|

## Environment notes
- capabilities: allocation cpus=8, mem=31.2 GB, time left at 16:32 UTC=4:18; local/ used=1.75 GB;
  apptainer=1.4.4, micromamba=2.3.3, nextflow=26.10.0, nf-test=0.9.5, gh=2.102.0, sbatch=ok,
  gpu partition=visible
- databases present: Bakta light+full v6.0, geNOMAD, Promotech RF-HOT, antiSMASH (4.1 GB); genomes K-12,
  Nissle; UP000000625 proteome
- sources: all §3.2 sources reachable except regulondb.ccg.unam.mx

## Measured numbers
| Step | Genome | Wall | Peak RSS | Node class (Slurm job ID) |
|---|---|---|---|---|
| bootstrap (micromamba env giap-wf) | — | ~3 min | — | iron-01 (session) |
| uv sync --all-extras | — | 1 min 48 s | — | iron-01 (session) |
| Bakta full DB | Nissle 1917 | 13.5 min | — | wheat (1846846) |
| TerminatorNet (-n 17) | K-12 MG1655 | 60 min | — | wheat-04 (1847778) |
| giap setup GO/KEGG/COG/contaminants | — | — | go-basic 32.2 MB, KEGG links 2.3 MB | session (review #67) |
| Bakta 1.12.2 full DB 6.0, circular-tagged, `--replicons` (M2-10) | K-12 MG1655 | 9.9 min (592 s) | 20.6 GB (job, with SIM100K) | barley-01, 16 CPU (1849096) |
| Bakta 1.12.2 full DB 6.0, circular-tagged (M2-10) | SIM100K | 5.3 min (317 s) | (same job) | barley-01, 16 CPU (1849096) |
| Bakta 1.12.2 full DB 6.0, circular-tagged (M2-10) | Nissle 1917 | 9.9 min (595 s) | 20.5 GB | barley-01, 16 CPU (1849097) |
| Bakta 1.12.2 light DB 6.0 (M2-10) | SIM100K | 34 s | 5.9 GB | barley-01, 16 CPU (1849098) |
