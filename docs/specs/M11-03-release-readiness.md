# M11-03: release readiness (`scripts/release/readiness.py`)

Milestone M11 · **pass 1 (MVP)** · depends on: M11-01 (version source), M11-02 (soft: changelog
rules) · estimate ≈ 340 code lines + ≈ 320 test lines · not `stats` (it compares committed numbers
with fixed thresholds; no inference) · roles: implementer (`scripts/release/readiness.py`),
test-author (`tests/unit/test_m11_03_*.py`), validation (runs it on FarmShare and commits the output
in M11-13) · Khoa: Q-K5, Q-K6, Q-K12 (criteria readings).

## Spec sections implemented

- §13 row M11 gate: "the acceptance criteria below"; §13.2, all seven bullets (quoted in the table
  below).
- §12 `release.yml`: "verify gates, changelog and the version".
- §11 "a real E. coli run gates every release"; §11.2 e2e "weekly and before release".
- §13 "A milestone is done when its gate passes and `PROGRESS.md` records the evidence".
- R3: every threshold below is §13.2's or §11.6's own; the checker never reads a threshold from a
  config file it could be tuned through. R11: the outputs hold aggregates only.

## Why a checker

The evidence lives in five places (committed pilot metrics, `benchmarks/baseline.json`, the FarmShare
`nightly.tsv`/`weekly.tsv` logs, PROGRESS gate blocks, the repository itself). Khoa approves one
table, not five sources. `release.yml` runs on GitHub and cannot read `$GIAP_ROOT`, so the table is
computed on FarmShare, committed, and re-verified by the workflow (Q-L4).

## Interface

```text
python scripts/release/readiness.py compute --version 1.0.0 [--repo .] [--logs $GIAP_ROOT/logs]
        [--freeze-run DIR_A2] --out validation/release/1.0.0/
python scripts/release/readiness.py verify validation/release/1.0.0/readiness.json [--repo .]
```

`compute` writes `readiness.json` and `READINESS.md` (rendered from the JSON only). Exit 0 when
`overall = pass`, 1 when any row is `fail`, 3 when no row fails but some are `pending`.
`verify` re-checks, without FarmShare, the rows marked `repo` below plus the commit rule; exit 0/1.

### `readiness.json`

```json
{"version": "1.0.0", "commit": "<40-hex R>", "computed": "<ISO-8601 UTC>",
 "overall": "pass|fail|pending",
 "rows": [{"id": "R1", "criterion": "...", "status": "pass|fail|pending|deferred|waived",
           "value": "...", "threshold": "...", "evidence": ["path or log line", "..."],
           "waiver": "ADR-NNNN or null", "source": "repo|farmshare"}]}
```

`overall` = `fail` if any row fails; else `pending` if any is pending; else `pass` (`deferred` and
`waived` count as pass).

## Rows

| id | Criterion (spec words) | Rule | Evidence read | source |
|---|---|---|---|---|
| R1 | "≥ 99.5 % of A2 peptides map exactly" | `metrics.json › te1.share_exact_inframe ≥ 0.995` and `te1.passed` | `validation/pilot/metrics/metrics.json` | repo |
| R2 | "the decoy-estimated mapping FDR is ≤ 1 % at defaults" | the A2 freeze run's overall placement FDR at the default threshold ≤ 0.01 | `metrics.json` key if present (gap G1), else `--freeze-run/mapping/mapping_fdr.tsv` | farmshare |
| R3 | "Every statistical method passes its calibration suite (§8.4)" | the newest `calibration_full` and `calibration` nightly rows for commit R (or an ancestor with no change under `src/`, `tests/calibration/`, `ci/`) have exit 0 and 0 failed | `logs/nightly.tsv` | farmshare |
| R4 | "every attack null (§8.5) stays at level" | as R3, restricted to `tests/calibration/attacks` (its own row in the run summary) | `logs/nightly.tsv` | farmshare |
| R5 | "Every true-example assertion that the available data allow (§11.6) passes" | TE1 (R1), TE2–TE5 real (`te_real`), TE7 (`te7`), TE8 (newest weekly `te8` row for R) each `pass`; a TE without data is `deferred` with the reason; `reported (not passed)` is **fail** (Q-K5) | `metrics.json`, `logs/weekly.tsv` | both |
| R6 | "The same seed gives byte-identical tables" | `metrics.json › te7.passed` and the TE7 integration tests in R3's run | as named | both |
| R7 | "A standard run on 16 threads and 64 GB takes about 35 min (range 16–55 …)" | `baseline.json › standard_run.wheat.k12_standard.median_s` in [960, 3300] (Q-K6); value printed as `h:mm`, node class shown | `benchmarks/baseline.json` | repo |
| R8 | "Core modules alone take about 25 min" | reported, never fails (spec gives no range): `k12_core.median_s` shown, status `pass` if present else `pending` | `benchmarks/baseline.json` | repo |
| R9 | "The pilot passes (§15.6)" | `PROGRESS.md` has `### Gate M10` with a Khoa review line **and** `validation/pilot/REPORT.md` status table has no `fail`/`pending` | PROGRESS, REPORT.md | repo |
| R10 | "co-abundance … moves to the full validation" | always `deferred` (spec text quoted) | none | repo |
| H1 | milestone gates | `PROGRESS.md` has a `### Gate M<k>` block for k = 1…10 (Q-K12) | PROGRESS | repo |
| H2 | version agreement | `giap.__version__` = `--version` = `CITATION.cff › version`; tag `v<version>` does **not** exist (`git tag -l`) | repo | repo |
| H3 | changelog | `CHANGELOG.md` has `## GIAP <version> (` and `changelog.d/` holds no fragment (built) | repo | repo |
| H4 | containers by digest | every `container =` in `workflow/conf/containers.config` is `docker://…@sha256:<64 hex>`; every process of `src/giap/annotate/tools.py`'s registry is either pinned or in the documented conda-only list (UniOP, fgsea, and those M11-06 has not yet published) | repo | repo |
| H5 | licences | every tool in the registry appears in `NOTICE.md`'s table; `LICENSE` is MIT | repo | repo |
| H6 | CI on R | GitHub check runs for R all `success` (`gh api repos/{owner}/{repo}/commits/R/check-runs`) | GitHub | farmshare (also re-run by `release.yml`) |

A row with an accepted waiver ADR (`Status: accepted-by-khoa` and the row id in its body, e.g.
"Waives readiness R7") becomes `waived` and cites the ADR. Nothing else turns a fail into a pass.

### Commit rule (`verify`)

`verify` passes the commit rule when HEAD = `commit`, or when every path changed in `commit..HEAD`
is in the paperwork allow-list: `validation/release/**`, `PROGRESS.md`, `docs/**` (not
`docs/methods/**`), `README.md`, `CITATION.cff`, `CHANGELOG.md`, `changelog.d/**`. Any other path →
`fail (stale: <first path>)`.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | fixture repo: te1 share 0.9971, passed; baseline k12_standard 2,280 s; all gates; nightly rows exit 0 | R1 `pass` value `99.71 %`; R7 `pass` `0:38`; overall `pass`, exit 0 |
| 2 | te1 share 0.9949 | R1 `fail` (`99.49 % < 99.5 %`); overall `fail`, exit 1 |
| 3 | te1 share exactly 0.995 | R1 `pass` (inclusive bound) |
| 4 | k12_standard 3,400 s, no ADR | R7 `fail` (`0:57 > 0:55`) |
| 5 | as 4, plus `docs/adr/ADR-0060-…md` `Status: accepted-by-khoa` "Waives readiness R7" | R7 `waived`, evidence cites ADR-0060; overall not failed by R7 |
| 6 | as 4, ADR status `proposed` | R7 stays `fail`; the ADR is named in `evidence` as "not accepted" |
| 7 | `metrics.json` te_real: TE3 `reported (not passed, explained)` | R5 `fail` naming TE3 (Q-K5 default) |
| 8 | no `### Gate M9` in PROGRESS | H1 `fail (missing: M9)` |
| 9 | `nightly.tsv` newest calibration row is for commit X, and X..R changes `src/giap/stats/bh.py` | R3 `pending (no calibration run for R)` |
| 10 | X..R changes only `docs/user-guide/install.md` | R3 uses X's row: `pass`, evidence `nightly.tsv 2026-11-02 <job> X` |
| 11 | `verify` with HEAD = R + one commit touching `validation/release/1.0.0/READINESS.md` | commit rule `pass` |
| 12 | `verify` with HEAD = R + a commit touching `src/giap/cli.py` | exit 1, `fail (stale: src/giap/cli.py)` |
| 13 | `containers.config` has `container = 'docker://oschwengers/bakta:v1.12.2'` | H4 `fail (tag, not digest: BAKTA)` |

## Edge cases

- **Missing evidence files** (`metrics.json` absent before M10-15b): the row is `pending
  (missing: <path>)`, never a crash; exit 3 if nothing fails.
- **Malformed log lines** in `nightly.tsv` (fewer than 5 tab fields, a killed task): skipped and
  counted in `READINESS.md` › Notes; never treated as a pass.
- **Several log rows for one commit**: the newest by date then job ID; a later failure beats an
  earlier pass.
- **Ties at a bound**: all bounds inclusive as the spec words them (≥ 99.5 %, ≤ 1 %, 16–55 min).
- **Tag already exists** (H2): `fail (tag v1.0.0 exists)`; the checker never creates or deletes tags.
- **Row-level data**: `READINESS.md` prints counts, shares, durations and paths only; a test greps
  the rendered file for 7+ letter amino-acid strings (the M10-05 guard rule) and fails on any.
- **Detached HEAD / shallow clone** in Actions: `verify` fetches nothing; `release.yml` checks out
  with full history (M11-04).
- Genomic edge cases (strand, circular): not applicable; R2 uses whatever the run's FDR table says.

## Acceptance criteria

1. Every row of the table is implemented with a unit test on a fixture repo in `tmp_path`
   (fake PROGRESS, metrics.json, baseline.json, nightly/weekly TSVs, a git repo with two commits);
   examples 1–13 are tests.
2. `compute` and `verify` never touch the network except H6, which is skipped (`pending (offline)`)
   when `--offline` is passed; tests always pass `--offline`.
3. `READINESS.md` is a pure function of `readiness.json` (a test renders twice: byte-identical).
4. No threshold appears in a config file: thresholds are module constants with the spec quote in a
   comment (R3); a test asserts the four constants (0.995, 0.01, 960, 3300).
5. `ruff`, `mypy` (scripts are outside the strict `src/giap` target; type hints anyway) pass.

## Gaps (for the lead)

- **G1.** No committed key holds A2's "decoy-estimated mapping FDR at defaults". Default: R2 reads
  `metrics.json › mapping_fdr` if the validation role adds it (an aggregate, allowed by Khoa's data
  rule), else the freeze run's `mapping/mapping_fdr.tsv` passed with `--freeze-run`; which row is
  "at defaults" must be confirmed against M3-10's table layout by the implementer and stated in a
  comment.
- **G2.** The nightly summary rows carry a free-text summary; R3/R4 parse `N passed`, `N failed`,
  `N error` from it. If attack nulls are not a separate row today, R4 is `pending (no attack row)`
  until nightly reports them separately (a nightly change, not this issue).

## Estimate

Code ≈ 340 (rows ≈ 220, commit rule ≈ 40, render ≈ 50, CLI ≈ 30). Tests ≈ 320.

## Questions for the lead (reversible defaults in use)

Q-K5, Q-K6, Q-K12 (index); Q-L4 (committed readiness + commit rule).
