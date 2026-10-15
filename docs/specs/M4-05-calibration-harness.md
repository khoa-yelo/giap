# M4-05: Calibration harness and stats-engine calibration tests

Milestone M4 · track S · depends on: M4-01, M4-03 (for the engine tests; the harness itself has no
dependency) · estimate ≈ 0 `src/` lines; ≈ 150 lines of test infrastructure + ≈ 200 test lines ·
**can start now** · authored by the **test-author** (tests/ only), reviewed by the statistician

This issue has no implementer step for the harness: it is test infrastructure (§14.1 test author owns
`tests/`). The engine calibration tests it adds go red until M4-01/03 merge, then green.

## Spec sections implemented

- §8.4 quoted in full for the criteria: "Level. For α ∈ {0.01, 0.05, 0.10}, the null rejection rate is
  ≤ α + 2√(α(1 − α)/R), where R counts independent null replicates." "Uniformity (continuous p only).
  KS test against U(0, 1), passing at p > 0.001. Discrete permutation p-values use the level check
  only." "FDR. Mean FDP over R replicates ≤ q + 2 SE(FDP)." "FWER. The share of replicates with ≥ 1
  false discovery is ≤ 0.05 + 2 SE." "Fixed seeds. Seeds are fixed in the tolerances file, and tests
  are deterministic."
- §8.4: "The per-PR subset MUST finish in under 10 minutes on a 4-vCPU GitHub runner"; "The full
  suite runs nightly on FarmShare as a Slurm job array, one suite per task."
- §11.2 markers `calibration` (fast, per PR) and `calibration_full` (nightly).
- `tests/calibration/tolerances.yaml` (protected): `levels`, `se_multiplier: 2.0`, `ks_min_p: 0.001`,
  `seeds.calibration_base: 20261014`, `replicates_fast`, `replicates_full`.

## Interface (test infrastructure)

`tests/calibration/_criteria.py`
```python
@dataclass(frozen=True)
class Tolerances: ...                                   # parsed tolerances.yaml (read-only)
def load_tolerances() -> Tolerances: ...
def replicates(suite: str, full: bool) -> int: ...      # replicates_fast / replicates_full[suite]
def suite_seed(suite: str, rep: int) -> int: ...        # SeedSequence([calibration_base, crc32(suite), rep]) → int
def level_bound(alpha: float, R: int) -> float: ...     # α + se_multiplier·√(α(1−α)/R)
def check_level(p: ArrayLike, alpha: float) -> CheckResult: ...       # rate of p ≤ α vs level_bound
def check_uniform(p: ArrayLike) -> CheckResult: ...                   # KS p > ks_min_p
def check_fdr(fdp: ArrayLike, q: float) -> CheckResult: ...           # mean ≤ q + k·sd/√R
def check_fwer(any_false: ArrayLike, alpha: float = 0.05) -> CheckResult: ...
@dataclass(frozen=True)
class CheckResult: ok: bool; observed: float; bound: float; R: int; detail: str
```
`tests/calibration/conftest.py`: fixture `full` = true under `-m calibration_full`; every result is
appended to `local/logs/calibration/<suite>.jsonl` (suite, criterion, α or q, observed, bound, R,
seed, duration) for the PR summary and `benchmarks/baseline.json` timings.

## Engine calibration tests added here

| Test | Null generator | Criterion | R (fast / full) |
|---|---|---|---|
| `test_bh_fdr` | 1,000 p per replicate: 900 U(0,1) nulls + 100 Beta(0.1, 1) alternatives | FDR at q ∈ levels | `gene_fdr`: 50 / 300 |
| `test_perm_p_level` | two groups of 30 iid N(0,1); stat = difference of means; B = 199 label permutations | level (discrete) | `regions`: 100 / 1,000 |
| `test_circular_shift_level_iid` | 1,000 iid N(0,1) values on a circular index; region = 20 consecutive indices; stat = sum; exhaustive shifts | level | `regions`: 100 / 1,000 |
| `test_besag_clifford_level` | same as `test_perm_p_level`, h = 10, max_B = 1,000 | level | `regions`: 100 / 1,000 |

## Examples

| # | Call | Expected |
|---|---|---|
| 1 | `level_bound(0.05, 100)` | 0.05 + 2·√(0.0475/100) = 0.093589 |
| 2 | `level_bound(0.01, 1000)` | 0.01 + 2·√(0.0099/1000) = 0.016293 |
| 3 | `check_level([0.001]*3 + [0.5]*97, 0.01)` | observed 0.03 > bound 0.0299 → `ok=False` |
| 4 | `check_level([0.001]*2 + [0.5]*98, 0.01)` | observed 0.02 ≤ 0.0299 → `ok=True` |
| 5 | `check_fdr([0.0]*49 + [0.5], 0.01)` | mean 0.01, sd 0.0707, SE 0.0100 → bound 0.0300 → ok |
| 6 | `check_fwer([True]*6 + [False]*94)` | 0.06 ≤ 0.05 + 2·√(0.0475/100) = 0.0936 → ok |
| 7 | `replicates("regions", full=False)`; `replicates("regions", full=True)` | 100; 1,000 |
| 8 | `suite_seed("regions", 0)` called twice | identical ints; differs from `suite_seed("regions", 1)` |

## Edge cases

- `check_fdr` with all FDP = 0 → SE 0, bound = q, ok.
- An unknown suite name → `KeyError` (no silent default R).
- `check_uniform` refuses discrete p (more than 5 % ties) with a clear error, per §8.4.

## Acceptance criteria

- [ ] `tests/calibration/test_criteria.py` (unit tests of the harness, unmarked so they run every
  commit): examples 1–8.
- [ ] The four engine tests above, marked `calibration` (fast R) and parametrised to `calibration_full`.
- [ ] Fast subset of this file finishes in < 60 s on 4 CPUs (recorded in the jsonl log).
- [ ] The harness never writes to `tolerances.yaml`; no constants duplicated outside it.

## Methods note (statistician sign-off)

`docs/methods/calibration.md`: the four criteria, the meaning of R (independent replicates), how SE is
computed, the discrete-p rule, seed derivation, and how suites map to `tolerances.yaml` keys.

## Questions for the lead

1. `tolerances.yaml` has no suite for the engine itself. **Proposed default:** engine tests reuse the
   R of the closest consumer suite (`gene_fdr` for BH, `regions` for permutation/shift/BC) rather
   than adding keys to a protected file (which needs an ADR).
2. Results log path `local/logs/calibration/*.jsonl` (gitignored) for PR summaries? **Proposed
   default:** yes.

## Lead decisions after statistician review (PR #54, BLOCK)
1. Fail closed: `check_level`, `check_fwer` (and `check_fdr` for FDP) raise ValueError on non-finite or
   out-of-[0,1] input; `check_fdr` rejects R < 2. The statistician's two unmarked tests must pass.
2. `--shard K/N` is honoured: suites read `GIAP_CALIBRATION_REPLICATES=start:stop` through the harness and
   run only that range; a merge step (`merge_shards` / `--merge`) combines shard records and evaluates the
   criterion on the full R. Without the variable, behaviour is unchanged.
3. ADR-0015 (accepted by Khoa): two-block confirmation — a failing check is re-run on replicates R..2R−1;
   the gate fails only if both blocks fail; records carry `block` and `status` (pass | unconfirmed | fail).
4. Suite assignment is explicit (a marker or a registry), not matched by name in source text, so a test runs
   in exactly one array task.
5. Sharded runs honour ADR-0015: `--merge` evaluates block 1 from the shard records and, for every check
   that fails, runs block 2 itself (replicates R..2R−1, sequentially in the merge task) before the verdict;
   it never reports `fail` without block 2. `--block 2` stays available for manual reruns.
6. Accepted implementer choices: `check_fwer` accepts only 0/1 outcomes; `check_uniform` fails closed;
   suite assignment by `SUITE_REGISTRY` or a `calsuite` marker; sim_v2 tests in suite `missed_operons`.
7. Follow-up (backlog, test-author): migrate the existing calibration tests from direct
   `run_replicates`/`check_*` calls to `calibrate()`, so ADR-0015 and sharding apply to every suite.
