# GIAP progress
_Last updated: 2026-10-14 · session https://claude.ai/code/session_01LEyHKvwXcvHWvkxzhtD54J_

## Current milestone
M0: Scaffold. Status: in progress

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
- m0/scaffold: failing tests for config, CLI and simulator (test-author), then implementation.
- Slurm jobs: none running.

## Next actions (ordered)
1. Implement config model, CLI skeleton, schema skeleton and simulator skeleton to the M0 tests.
2. Code review; open the "M0 scaffold" PR once `gh` is authenticated; get CI green.
3. M1: io/ readers and writers, generic peptide table, `giap convert` (generic, DIA-NN), preflight rules.

## Blocked
| Since | What | Host or tool | Needed from Khoa |
|---|---|---|---|
| 2026-10-14 | RegulonDB v12 download (needed at M6) | `regulondb.ccg.unam.mx` unreachable from FarmShare (probe) | nothing yet; retry at M6, then ask |

## Questions for Khoa
| # | Question | Options | Default if no answer | Asked |
|---|---|---|---|---|
| 1 | ADR-0004: may PR #1 (M0 scaffold, ~4,800 lines + the spec) merge with the PR-size guard failing, as a one-off? | yes / split into ~6 PRs | wait; no merge of PR #1 | 2026-10-14 |
| 2 | ADR-0004: may I apply the two-line `realpath` fix to the protected `.claude/hooks/protect_paths.py` (closes a symlink bypass)? | yes / no | leave as is; CI `check_protected.py` still catches protected-file changes | 2026-10-14 |

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
