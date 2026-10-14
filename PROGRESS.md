# GIAP progress
_Last updated: 2026-10-14 · session https://claude.ai/code/session_01LEyHKvwXcvHWvkxzhtD54J_

## Current milestone
M0: Scaffold. Status: done (PR #1, merged 2026-10-14 under ADR-0004)
M1: Inputs and preflight. Status: in progress — 9 of 12 issues merged (M1-01…09 incl. this PR); M1-10…12 left

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
- Merged to main (squash, one commit per PR): #1 M0, #2 M1-01, #3 M1-02, #4 M1-03, #5 M1-04,
  #6 M1-09, #7 M1-06, #8 M1-05, #9 M1-08; this PR: M1-07 (`giap convert --from generic|diann`).
- Untracked specs M1-10, M1-11, M1-12 in the main checkout's docs/specs/ (commit with each PR).
- Slurm jobs: `giap-fetch-bakta-db` (scripts/slurm/fetch_bakta_db.sbatch; log local/logs/giap-fetch-bakta-db.<id>.out).
- Lessons:
  - Review fixes can contradict older tests → the test-author lists and aligns them first; the
    implementer never edits tests and stops on contradictions.
  - Subagents must never run `uv sync` / `pip install`: the venv is shared, and an install in a
    worktree repoints the editable `giap` at that worktree. Always set
    `PYTHONPATH=$PWD/src:$PWD/tests/sim` when running tests.
  - Merging a stack of squash-merged PRs: retarget the PR to main *before* pushing the rebased
    branch, and merge only from one snapshot showing every check present and passed.

## Next actions (ordered)
Parallel orchestration (Khoa, 2026-10-14: "design Mx in parallel to speed up development"):
- Wave 1 (now): finish M1-10/11/12 (sequential: 11 builds on 10, 12 on 11). In parallel:
  spec-writers for M2, M3, M4 (design only, docs/specs/M{2,3,4}-*), and the long-latency
  Slurm prefetch of the Bakta light+full DB (job `giap-fetch-bakta-db`).
- Wave 2: implement M2 and M3 side by side, each split into parallel tracks from its index
  (one worktree per track, test-author → implementer → code-reviewer per issue); M4 issues that
  only need the stats engine or fixtures start alongside.
- Wave 3: M4 integration, then M5–M8 in parallel (§13), each with its own spec-writer first.
- Concurrency limits: ≤ 6 subagents at once; inside the 8-CPU session allocation, subagents run
  targeted tests with `-n 2` and the full suite once before committing; anything long → sbatch.
  R14 caps apply to GIAP's Slurm jobs (≤ 128 CPUs, ≤ 20 jobs); the account already runs other work.

## Blocked
| Since | What | Host or tool | Needed from Khoa |
|---|---|---|---|
| 2026-10-14 | RegulonDB v12 download (needed at M6) | `regulondb.ccg.unam.mx` unreachable from FarmShare (probe) | nothing yet; retry at M6, then ask |

## Answered by Khoa
- 2026-10-14: ADR-0005 item 4 (new preflight warning thresholds) rejected: decoy share and free RAM/disk are Info only.
- 2026-10-14: ADR-0004 (M0 size exemption; hook realpath fix) and ADR-0007 (size guard counts code only) accepted; "merge the stack".

## Questions for Khoa
| # | Question | Options | Default if no answer | Asked |
|---|---|---|---|---|

## Environment notes
- capabilities: allocation cpus=8, mem=31.2 GB, time left at 16:32 UTC=4:18; local/ used=1.75 GB;
  apptainer=1.4.4, micromamba=2.3.3, nextflow=26.10.0, nf-test=0.9.5, gh=2.102.0, sbatch=ok,
  gpu partition=visible
- databases present: none yet (Bakta light/full, geNOMAD, Promotech models come in M2)
- sources: all §3.2 sources reachable except regulondb.ccg.unam.mx

## Measured numbers
| Step | Genome | Wall | Peak RSS | Node class (Slurm job ID) |
|---|---|---|---|---|
| bootstrap (micromamba env giap-wf) | — | ~3 min | — | iron-01 (session) |
| uv sync --all-extras | — | 1 min 48 s | — | iron-01 (session) |
