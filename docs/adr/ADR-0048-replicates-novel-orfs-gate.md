# ADR-0048: Replicate counts for the 1-Mb novel-ORF gate (M7-13 G3/G4)
- Status: proposed
- Date: 2026-10-23
- Spec sections: §8.4 (suite "Novel-ORF class FDR", 50 / 300; "set R and the genome sizes to fit"
  the 10-minute budget; one suite per array task), §7.6.4, §7.6.6, §13 gate M7; M7-12 decision 26,
  M7-13 Q1
- Needs Khoa: yes (protected `tests/calibration/tolerances.yaml`)

## Context
- G3 (novel-ORF F1) and G4 (class FDR) can only be measured at 1 Mb with
  `NOVEL_PER_CLASS_GATE = 300` ORFs per class: the picked class q `(D + 1)/T` cannot fall below
  1/T, so the default counts (16 at 300 kb, 40 at 1 Mb) call nothing at α = 0.01 (M7-12
  decision 26, statistician round 3 on #225).
- Each gate replicate is one 1-Mb simulation and a GIPPA run per regime (database and de novo),
  assumed ≈ 90 s each (to be measured in M7-13a).
- `replicates_fast.novel_orfs` = 50 is shared with the cheap simulator checks of the same suite
  (`test_sim_m7_*`, `test_attack_m7_sixframe_picked.py`). At 50 × 2 × 90 s ≈ 2.5 h single-core the
  gate cannot run as a fast (nightly) test, and its 300 full replicates (≈ 15 h CPU) need a
  sharded array task of their own, not task 11's 2-hour slot shared with the other tests.

## Decision (proposed)
Add to `tests/calibration/tolerances.yaml`:
- `replicates_fast.novel_orfs_gate: 10`: pooled over ≈ 600 recall-truth ORFs per replicate, F1 has
  an SE well below 0.01; at q = 0.10 with ≈ 150 calls per class the FDP SE over 10 replicates is
  ≈ 0.008.
- `replicates_full.novel_orfs_gate: 300`: the §8.4 full count, unchanged.
- The key is appended after `compare`, so it becomes array task 15 (`--array=0-15`); existing task
  numbers do not move.

Until Khoa decides, M7-13 uses the module constant
`GATE_R = MappingProxyType({"fast": 10, "full": 300})` in `tests/calibration/_m7_gate.py` and logs
under suite `novel_orfs` (task 11). Accepting this ADR replaces the constant by
`C.replicates("novel_orfs_gate", full)` in a PR carrying `Gate-Change-Approved: ADR-0048`.

## Alternatives
- Keep 50 / 300 under `novel_orfs`: the fast gate runs ≈ 2.5 h every night; the full gate needs
  hand-submitted shards of task 11, mixed with the other `novel_orfs` tests.
- Run G3/G4 fast at 300 kb: vacuous at q = 0.01 (class q floor), so not a gate.

## Evidence
To be added by M7-13a: measured per-replicate cost at 1 Mb, and F1 / FDP SE at R = 10.

## Consequences
- The fast tier for G3/G4 is a 10-replicate nightly check; the per-PR subset keeps only a
  one-replicate smoke test.
- The full tier gets its own array task, submitted sharded (`--shard K/10`, then `--merge`).
