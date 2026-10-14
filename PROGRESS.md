# GIAP progress
_Last updated: 2026-10-14 · session https://claude.ai/code/session_01LEyHKvwXcvHWvkxzhtD54J_

## Current milestone
M0: Scaffold. Status: gate pending (PR #1 green except the PR-size guard; ADR-0004 awaits Khoa)
M1: Inputs and preflight. Status: in progress (12 issues, docs/specs/M1-index.md)

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
- PR #1 `m0/scaffold` (M0): all CI green except PR size (8,585 > 800) → ADR-0004 Q1.
- PR #2 `m1/coordinates` (M1-01, base m0/scaffold): reviewed, findings fixed; 206 tests.
- PR #3 `m1/genome-fasta` (M1-02, base m1/coordinates): in code review.
- `m1/tsv-io` (M1-04, worktree local/worktrees/tsv-io): implementer working; rebase onto m1/coordinates before PR.
- `m1/protein-fasta` (M1-03, worktree local/worktrees/protein-fasta, base m1/genome-fasta): test-author working.
- Untracked specs M1-05…M1-12 in the main checkout's docs/specs/ (committed with each issue's PR).
- Worktrees share the venv; set `PYTHONPATH=$PWD/src:$PWD/tests/sim` inside a worktree.
- Slurm jobs: none running.

## Next actions (ordered)
1. Finish M1-02/03/04 loops (review → fixes → PR); then M1-05 (after 04), M1-06 → 07/08, M1-09 (after 02).
2. When Khoa answers ADR-0004: merge #1, retarget and merge the stacked PRs in order.
3. M1-10/11 need ADR-0005 item 4 (new preflight thresholds) accepted by Khoa.

## Blocked
| Since | What | Host or tool | Needed from Khoa |
|---|---|---|---|
| 2026-10-14 | RegulonDB v12 download (needed at M6) | `regulondb.ccg.unam.mx` unreachable from FarmShare (probe) | nothing yet; retry at M6, then ask |

## Answered by Khoa
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
