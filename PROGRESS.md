# GIAP progress

> **Resuming after a timeout:** start an allocation, run `bash scripts/resume_status.sh`, then use the
> resume prompt (GIAP_SPEC.md §0.3). Every worktree's next feature-loop step is printed; uncommitted
> files mean an agent stopped mid-step — inspect them, then re-dispatch that step (test-author →
> implementer → code-reviewer). Slurm jobs survive the session; their logs are in `local/logs/`.
_Last updated: 2026-10-15 06:30 PDT · session https://claude.ai/code/session_01LEyHKvwXcvHWvkxzhtD54J_

## Milestone progress
<!-- progress:start -->
🟩 done (merged) · 🟨 running (open branch) · ⬜ pending; designed issues: 57/79 done (72 %), 8 running · generated 2026-10-15 06:29 by `scripts/dev/progress_bars.py`

| Milestone | Progress | Done | Running | Pending |
|---|---|---|---|---|
| M0 Scaffold | 🟩 | 1 | 0 | 0 |
| M1 Inputs and preflight | 🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩 100 % | 12 | 0 | 0 |
| M2 Annotation wrappers | 🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟨🟨⬜⬜⬜⬜⬜ 67 % | 14 | 2 | 5 |
| M3 Mapping and decoy FDR | 🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟨🟨⬜⬜ 83 % | 20 | 2 | 2 |
| M4 Gene table and region tests | 🟩🟩🟩🟩🟩🟩🟩🟩🟩🟩🟨🟨🟨🟨⬜⬜⬜⬜⬜⬜⬜ 48 % | 10 | 4 | 7 |
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
Merged (all with code review; stats PRs with statistician sign-off): M2-02…M2-06, M2-08 (#57), M2-11 (#64),
M2-12 (#49), M3-04…M3-15 (incl. #41 placement FDR, #48 tolerant seeds, #50 loci, #53 tblastn, #60
frameshifts, #63 contaminants), M3-19 (#55), M4-02, M4-04…M4-07 (#51, #54 calibration harness, #62),
M4-10 (#61), ADR-0016 size gate (#56), Slurm helpers (#46, #58, #59). Run `bash scripts/resume_status.sh`.
Open PRs:
- #65 `m2/nf-tools` (M2-07): real-tool run failed review (registry keys, conda off, no `ps`, PyQt5,
  stub headers); fixes committed locally, real ANNOTATE run = Slurm 1847900; then re-review.
  Stacked follow-up `m2/run-wiring` (giap run → preflight, replicons, db_dir; `_step` loads steps) pushed.
- #66 M3-20 sim de novo errors (approved after test pins) — merging. #67 M2-09 setup resources (ADR-0017) —
  merging. #68 M3-16 mapping outputs — review fixes in progress. #69 calibration migration to
  `calibrate()` — statistician review.
- `m4/sim-m4-08` (M4-08 sim expression): test-author.
Slurm: tool recordings split per tool (giap-rec-*: 1847974 terminatornet, 1847975 genomad, 1847977
transtermhp, 1847979 promoter_calculator; chained 1847981 antismash, 1847982 ostir, 1847983 promotech);
the sequential job 1847778 was cancelled (TerminatorNet on K-12 alone took 60 min).
Follow-ups queued (docs/backlog.md has the full list):
- Per-process real-tool Slurm job array for workflow modules (inputs from recordings) — after #65 (Khoa).
- calibration_full variants for migrated suites (statistician); GO data-version pin (M2-09).
- M3-21: c_ref − c_true, chance rate of wrong peptides, PEP–length correlation (ADR-0010).
- M4-12/M4-17 forwarded items; M6: simulate readthrough (cterm_extension).
Lessons (keep):
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
