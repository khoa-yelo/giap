# ADR-0007: The PR-size guard counts code, not tests or docs
- Status: proposed
- Date: 2026-10-14
- Spec sections: §14.2, §14.6 (PR size), `ci/gates.yaml` (`pr_lines_warn`, `pr_lines_fail`)
- Needs Khoa: yes (changes how a gate is measured; the 400/800 values stay)

## Context
§14.2 asks for PRs of about 400 changed lines "excluding lock files and fixtures", and CI warns
above 400 and fails above 800. In practice the test-first loop (R2) plus independent review makes
tests the bulk of every PR. The first four M1 issues, each one small module:

| PR | src | tests | docs/other | total |
|---|---|---|---|---|
| M1-01 coordinates (#2) | 157 | 462 | 211 | 830 |
| M1-02 genome FASTA (#3) | 169 | 527 | 92 | 788 |
| M1-03 protein FASTA (#4) | 123 | 469 | 81 | 673 |
| M1-04 TSV/BED/proBed | 329 | 976 | 126 | 1,431 |

The guard fails #2 and M1-04 for their tests (review findings added 100–300 test lines each),
which pushes toward thinner tests, the opposite of R2 and the calibration culture of §8.

## Decision (proposed)
Count only code under `src/`, `workflow/` and `scripts/` toward `pr_lines_warn` (400) and
`pr_lines_fail` (800). Tests, docs, ADRs, specs, lock files and fixtures are reported in the PR
summary but not gated. The values in `ci/gates.yaml` do not change; only `ci.yml`'s counting does.

## Evidence (before → after: metrics, tests, simulations)
With code-only counting, all four PRs above pass (157, 169, 123, 329 lines), and M1-04 sits under
the 400-line warning.

## Consequences
Reviewers still see the full diff. A PR that adds 1,000 lines of tests for 50 lines of code is
fine; a PR that adds 900 lines of code still fails.
