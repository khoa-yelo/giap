# GIAP progress

> **Resuming after a timeout:** start an allocation, run `bash scripts/resume_status.sh`, then use the
> resume prompt (GIAP_SPEC.md §0.3). Every worktree's next feature-loop step is printed; uncommitted
> files mean an agent stopped mid-step — inspect them, then re-dispatch that step (test-author →
> implementer → code-reviewer). Slurm jobs survive the session; their logs are in `local/logs/`.
_Last updated: 2026-10-14 · session https://claude.ai/code/session_01LEyHKvwXcvHWvkxzhtD54J_

## Current milestone
M0: Scaffold. Status: done (PR #1, merged 2026-10-14 under ADR-0004)
M1: Inputs and preflight. Status: done (12 issues; PRs #2–#11, #13, #14)
M2–M4: designed (docs/specs/M{2,3,4}-index.md); wave 2 implementation in progress

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
Merged since wave 2: M3-04 (#27), M3-09 (#28), M2-03 (#30), M4-04 (#31), env locks (#32), resumable Bakta
download (#33), M4-02 (#34), M3-05 (#35), M2-04 (#36), M2-06 (#37), M2-05 (#38), M3-06 (#39), role guard
fix (#40), M3-11 (#42), timing budgets → benchmarks (#43). Run `bash scripts/resume_status.sh` for each
branch's next step.
Open:
- #41 `m3/placement-fdr` (M3-10, ADR-0010 + ADR-0011): statistician BLOCK fixed (same_db shortcut for
  database rows only; π̂ over remaining rows); statistician re-check running.
- `m2/nf-tools` (M2-07): tests + lead CLI test update pushed; implementer running.
- `m3/tolerant-seeds` (M3-07), `m3/miniprot` (M3-12), `m4/input-contract` (M4-06), `m2/uniop-tus` (M2-12):
  test-authors running.
Slurm: `giap-fetch-bakta-db` job 1846085 (aria2c, resumable, md5-checked, then `bakta_db install` in the
pinned container); resubmit `scripts/slurm/fetch_bakta_db.sbatch` if it times out.
Follow-ups queued:
- M3-21: measure c_ref − c_true (reference vs searched de novo precision) and the target-vs-decoy chance
  rate of wrong peptides; minimum reference peptides per class (or resample c_ℓ); count unreportable
  bootstrap resamples as precision 0 (M3-11 statistician). Measure the PEP–length correlation for the
  database π̂ (ADR-0010).
- M4-12: M4-02 forwarded items (input validation, report fp/n_genes, unspecific-search warning, `fdr_2pep`
  label); M4-17: report n_tu and flag bootstrap CIs from < 10 TUs (M4-04).
- M3-08: never call `verify` with sides ≤ 4 aa (M3-06 decision 1); anchored ends defined in M3-08.
- M3-07: specifies a public njit batch entry point over `tolerant_dp` (M3-06 decision 6).
- M3-06 benchmark at 3.7–4.2 s of the 5 s budget after input normalisation — watch in the speed loop.
- Skip pip freeze for pip-less envs when built by name; `_step` failure test should chdir to tmp_path;
  waitmerge.sh should read the PR head commit's check runs directly; background merges must use the tool's
  background mode (a plain `&` job dies with its shell).
Lessons (keep):
- Slurm jobs run from the main checkout only; never remove a worktree a job may be using.
- Clean up a worktree only after a confirmed merge; resolve rebase conflicts file by file.
- Subagents never install packages; set PYTHONPATH per worktree; unset PYTHONHOME after sourcing the env.
- Wall-clock budgets live in benchmarks/, unit tests keep only loose 10× guards.

## Next actions (ordered)
1. Finish #41 (statistician re-check → merge), M2-07, M3-07, M3-12, M4-06, M2-12.
2. Then: M3-08 frameshifts (after M3-07), M3-13 tblastn (after M3-12; HSP limits are constants — Khoa),
   M3-14 loci (after M3-10), M3-19 sim peptides, M4-05 calibration harness, M4-07…09 simulator (after
   M4-06), M2-08 setup core, M2-11 user annotation, M2-13…M2-20 tool layers (after M2-07/M2-12).
3. After the Bakta full DB installs: `giap setup` verification (M2-08) and recorded runs (M2-10).

## Blocked
| Since | What | Host or tool | Needed from Khoa |
|---|---|---|---|
| 2026-10-14 | RegulonDB v12 download (needed at M6) | `regulondb.ccg.unam.mx` unreachable from FarmShare (probe) | nothing yet; retry at M6, then ask |

## Answered by Khoa
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
- databases present: none yet (Bakta light/full, geNOMAD, Promotech models come in M2)
- sources: all §3.2 sources reachable except regulondb.ccg.unam.mx

## Measured numbers
| Step | Genome | Wall | Peak RSS | Node class (Slurm job ID) |
|---|---|---|---|---|
| bootstrap (micromamba env giap-wf) | — | ~3 min | — | iron-01 (session) |
| uv sync --all-extras | — | 1 min 48 s | — | iron-01 (session) |
