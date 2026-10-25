# ADR-0057: the calibration merge honours the xfail mark of the test that owns a check
- Status: accepted-by-khoa (2026-10-24, in session: "Yes and yes", to "should the merge step honour
  expected-failure marks?")
- Date: 2026-10-24
- Spec sections: §8.4 (calibration criteria, nightly suites), §12 (calibration job array), §14.6;
  ADR-0015 (two-block confirmation)
- Needs Khoa: yes (it changes when the calibration gate's exit code is a failure). Decided.

## Context
A calibration test may carry a strict `pytest.mark.xfail(raises=AssertionError)` naming a pending
decision. ADR-0055 Amendment 1 (accepted on m6/gate-09) is the motivating case:
`test_attack_d1b_step[full]` is measured to fail its pooled level at full R, and the claim is narrowed
by the amendment, never by widening a bound (R3).

Unsharded, pytest applies the mark. A confirmed `fail` raises `AssertionError` in the test, which is
reported XFAIL, and the task exits 0. A pass is a strict XPASS, and the task exits non-zero.

Sharded (`--shard K/N`), no shard has a verdict. Each shard logs its block and skips. The verdict
comes from `_criteria.py --merge`, which until now knew nothing about the test that wrote a record.
So a confirmed `fail` of a strict-xfail check exited `MERGE_FAIL`, and a pass exited 0. That is the
opposite of the unsharded run, so the same suite passed or failed depending on how it was run.

## Decision
1. **Owner in the record.** `tests/calibration/conftest.py` names the running test around every
   test call (`pytest_runtest_call` → `_criteria.set_current_test`). It does this with pytest's own
   evaluation of the xfail marks (`_pytest.skipping.evaluate_xfail_marks`, which handles conditions
   and `xfail_strict`). Every record that `calibrate`, `calibrate_pooled` and `calibrate_fdr_counts`
   write then carries `test` (the node ID) and `xfail`. `xfail` is `{"strict", "reason"}` when the
   mark would turn an `AssertionError` into XFAIL, and None otherwise: no mark, a false condition,
   `run=False`, or a `raises` that excludes `AssertionError`.

   This mapping comes from the tests themselves, not from a hand-kept list. Records of one check
   that disagree on `xfail`, or an xfail check written by two tests, are a merge error
   (`MERGE_ERROR`, fail closed). Logs written before this ADR, which have no `test` / `xfail`
   fields, merge as before.
2. **Per-test outcome, as pytest's.** A test with an xfail mark *failed* if any of its checks is a
   confirmed `fail` (both blocks, ADR-0015), or if one of its direct assertions failed in shard 0.
   In the second case conftest writes an `xfailed` line to `<suite>.outcomes.jsonl`
   (`_criteria.log_test_outcome`).
   - If the test failed, every one of its checks gets `outcome = "xfail"`.
   - If every check is `pass` / `unconfirmed` and the test did not fail, every check gets
     `outcome = "xpass"`.
   - A `vacuous` check leaves the test's outcome unset.
3. **Exit codes of `--merge`.**
   - A confirmed `fail` with `outcome = "xfail"` is an expected failure. It does not make the exit
     `MERGE_FAIL`.
   - A strict XPASS exits `MERGE_FAIL` (1), and the merge prints `XPASS(strict) <test>`.
   - A non-strict xfail is reported (`xfail` / `xpass` in the summary's status column and in
     `<suite>.merged.jsonl`) but never changes the exit code.
   - Any other confirmed `fail` is `MERGE_FAIL`, as before. Merge errors are unchanged.
4. **Unsharded task exit** is pytest's own exit code and already honours the mark. The teeth test
   (`tests/calibration/test_criteria_shard_equivalence.py`) checks that the unsharded and the
   sharded exits agree in every case.

No bound, alpha, R, level, seed or gate condition changes. Every check is evaluated exactly as
before (ADR-0015). The only change is how `--merge` turns a check's verdict into an exit code,
which now matches what pytest does with the same test unsharded.

## Evidence (before → after: metrics, tests, simulations)
- `tests/calibration/test_criteria_merge_xfail.py` covers the four cases on sharded logs:
  - strict-xfail fail → exit 0, `xfail`;
  - strict-xfail pass → exit 1, `xpass`, `XPASS(strict)`;
  - non-strict xfail → exit 0, reported only;
  - unmarked fail → exit 1.

  It also checks:
  - the per-test verdict (one failing and one passing check → `xfail`);
  - a shard-0 direct failure → `xfail`;
  - disagreeing marks → `MERGE_ERROR`;
  - pre-ADR logs merge unchanged.
- `tests/calibration/test_criteria_shard_equivalence.py` runs a fixture suite unsharded and as
  3 shards + `--merge`, through `_criteria.main` as `calibration_full.sbatch` does. It finds the same
  verdicts, the same pooled rates and the same exit in the `ok`, `xpass` and `fail` scenarios.
- Before: the same strict-xfail check exited 1 sharded and 0 unsharded (or 0 and 1 for a pass).

## Consequences
- A strict xfail now works the same in the nightly array (unsharded tasks), in task 13 (sharded in
  the task) and with `scripts/slurm/calibration_sharded.sh` (N shard jobs + a merge job).
- A test's xfail mark is the single source of truth. Removing the mark, for example once an ADR
  amendment narrows the claim, takes effect in sharded runs without touching the harness.
- A check that is shared by two tests, one with a mark and one without, can no longer be merged.
  Give each test its own stream.
- `<suite>.outcomes.jsonl` is a new file in the calibration log directory. It is read only by
  `--merge`.
