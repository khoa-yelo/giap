# GIAP progress

> **Resuming after a timeout:** start an allocation, run `bash scripts/resume_status.sh`, then use the
> resume prompt (GIAP_SPEC.md §0.3). Every worktree's next feature-loop step is printed; uncommitted
> files mean an agent stopped mid-step — inspect them, then re-dispatch that step (test-author →
> implementer → code-reviewer). Slurm jobs survive the session; their logs are in `local/logs/`.
_Last updated: 2026-10-14 · session https://claude.ai/code/session_01LEyHKvwXcvHWvkxzhtD54J_

## Current milestone
M0: Scaffold. Status: done (PR #1, merged 2026-10-14 under ADR-0004)
M1: Inputs and preflight. Status: 11 of 12 merged; M1-12 (`giap preflight`, PR #14) in review fixes
M2–M4: designed (docs/specs/M{2,3,4}-index.md with lead decisions); wave 2 started

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
- PR #14 `m1/preflight-cli` (M1-12): review findings pinned by tests; implementer next.
- Wave 2 worktrees (local/worktrees/<name>, branch from main, each issue: test-author → implementer → reviewer):
  `m2/tool-envs` (M2-01), `m2/step-registry` (M2-02), `m3/gene-model` (M3-01, shared with M2-03),
  `m3/sixframe` (M3-02), `m4/stats-core` (M4-01, reused by M3-10).
- Slurm: `giap-fetch-bakta-db` (job 1845787) downloading Bakta DB v6.0 light then full into local/db/bakta.
- Lessons (keep): test-author aligns contradicted older tests before the implementer starts; subagents
  never install packages; retarget before pushing rebased stacks; merge only on one all-green snapshot.

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
| 5 | M3-13: tblastn HSP chaining needs two parameters (max gap 300 nt, query overlap ≤ 10 aa). New config keys, or fixed constants documented in the methods note? | config keys (ADR) / constants / other | M3-13 waits | to ask when M3-13 starts |

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
