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
Stack (each PR's base is the previous branch): #1 m0/scaffold → #2 m1/coordinates (M1-01) →
#3 m1/genome-fasta (M1-02) → #4 m1/protein-fasta (M1-03) → #5 m1/tsv-io (M1-04) →
{#6 m1/preflight-genome (M1-09), #7 m1/peptide-table (M1-06), #8 m1/gff3-vcf (M1-05)} (siblings on m1/tsv-io).
- #2–#6: reviewed, findings fixed, done.
- #7 M1-06 and #8 M1-05: review findings pinned by tests; implementers fixing.
- Worktrees: local/worktrees/<topic> per branch; shared venv; inside a worktree set
  `GIAP_REPO=/home/classes/bios/270/khoa/giap PYTHONPATH=$PWD/src:$PWD/tests/sim`.
- Untracked specs M1-07, 08, 10, 11, 12 in the main checkout's docs/specs/.
- Slurm jobs: none running.
- Lessons: review fixes can contradict older tests → the test-author lists and aligns them before
  the implementer starts; the implementer never edits tests and stops on contradictions.

## Next actions (ordered)
1. Finish #7/#8 fixes; then M1-07 (convert generic/DIA-NN) and M1-08 (abundance) on m1/peptide-table;
   then M1-10 (needs 03, 06, 08, 09), M1-11, M1-12.
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
