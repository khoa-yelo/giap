# ADR-0002: FarmShare toolchain as installed in M0
- Status: accepted
- Date: 2026-10-14
- Spec sections: §3.3, §3.4
- Needs Khoa: no (versions match the spec's pins; uv is unpinned in the spec)

## Context
The bootstrap (§3.4) installs the developer toolchain in user space under `local/`.

## Decision
Record the versions the bootstrap installed: Nextflow 26.10.0 (bioconda, the spec's pin),
nf-test 0.9.5, gh 2.102.0, OpenJDK 21 (conda-forge, in env `gippa-wf`), uv 0.13.1
(standalone installer, `local/bin`). The Python lock (`uv.lock`) resolves numpy 2.5.4,
scipy 1.18.1, numba 0.68.0, pandas 3.0.6 for Python 3.12.

## Evidence (before → after: metrics, tests, simulations)
`scripts/probe_env.py` on iron-01 (session allocation) and on wheat-08 (Slurm job 1845260):
all tools present; every §3.2 source reachable except `regulondb.ccg.unam.mx` (needed at M6).

## Consequences
Upgrades go through a new ADR with contract-test diffs (§14.6).
