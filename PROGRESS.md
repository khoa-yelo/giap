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
Merged in wave 2: M2-01 tool envs (#21, pins in ADR-0008), M2-02 step registry (#20), M3-01 gene model
(#18), M3-02 six-frame index (#17), M4-01 BH/TDA/picked FDR (#19), ADR-0009 (#22).
Open (run `bash scripts/resume_status.sh` for each branch's next step):
- #23 `m4/permutation` (M4-03): statistician BLOCKED (two-sided centring); tests pinned; implementer
  fixing; then statistician re-check and merge.
- #24 `m3/decoys` (M3-03): statistician reviewing; writing docs/methods/placement-fdr.md (decoy part).
- #25 `m3/sim-features` (M3-18): review findings pinned by tests (start codons, reference aa, VCF-style
  indels); implementer next.
Slurm: `giap-build-envs` (job 1845936; all 12 envs created incl. Promotech TF 2.1.0, now pulling images;
rerun with the merged script afterwards to get pip locks and markers); `giap-fetch-bakta-db` (job 1845787;
light DB nearly done; full DB ~14 h at Zenodo's ~0.7 MB/s; resubmit if it times out — it skips present DBs).
Follow-ups queued: skip pip freeze for pip-less envs (transtermhp, blast, miniprot) when built by name;
`_step` failure test should chdir to tmp_path; index save/load lives in sixframe.py (revisit at M3-17);
waitmerge.sh should read check runs of the PR head commit directly.

## Next actions (ordered)
1. Finish #23, #24, #25 (loops above) and merge on green.
2. Next parallel issues (specs on main): M3-04 exact matcher (needs 02, 03), M3-05…08 tolerant matcher
   track, M3-09 chance-hit model, M2-03 Bakta JSON parser (extends M3-01), M2-05 simulator layers (needs
   M3-18), M2-06/07 Nextflow skeleton, M4-02 MAYU, M4-04 bootstrap, M4-07…09 simulator.
3. Ask Khoa when M3-13 starts: tblastn HSP-chaining parameters (question 5).

## Blocked
| Since | What | Host or tool | Needed from Khoa |
|---|---|---|---|
| 2026-10-14 | RegulonDB v12 download (needed at M6) | `regulondb.ccg.unam.mx` unreachable from FarmShare (probe) | nothing yet; retry at M6, then ask |

## Answered by Khoa
- 2026-10-14: §7.3.2 picked FDR rescales decoy wins by the paired share (ADR-0009, accepted-by-khoa).
- 2026-10-14: subagent commit blocked by the permission system → Khoa approved the lead committing the M2-02 test fix.
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
