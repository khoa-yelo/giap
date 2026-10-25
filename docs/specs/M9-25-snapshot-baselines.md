# M9-25: Snapshot baselines and the approval packet (the M9 gate)

Milestone M9 · **pass 2 — the gate** · depends on: M9-13, M9-24, and every plot issue (it can open
with the pass-1 plots and add the rest as they merge) · estimate ≈ 150 code lines (+ ≈ 350 test
lines, + baseline images) · not statistical · **Khoa checkpoint (§14.8 item 3)**

## Spec sections implemented

- §10.1 "Snapshot tests. Every plot is snapshot-tested with pytest-mpl on the simulated genome
  (tolerance RMS ≤ 2). Khoa approves the baselines at M9, and after that a style change is a
  reviewed change."
- §10.3 tests "A snapshot of the HTML structure (headings and IDs)."
- §11.2 snapshot layer (`tests/snapshot`, every PR, < 2 min; `ci/gates.yaml` `budgets_s.snapshot`
  120), §12 `ci.yml` (snapshot tests on every push), §13 gate, §14.4, §14.8, R11 (file sizes).

## Work

1. **Frozen fixture** `tests/data/snapshot/sim_plots/`: an `OUT/` subset (only the files `RunView`
   reads; gzipped TSVs, the preflight JSON, the support files) produced by
   `scripts/dev/make_snapshot_fixture.py` from a `giap run --executor local` on a giapsim
   configuration rich enough that every plot is non-empty (planted island, prophage and BGC with
   low detection, flagged operon pairs, ALT sites, novel ORFs, function and determinants tables at
   ≥ 100 model rows), plus a second, header-only fixture for empty states. Target < 2 MB; a test
   checks the R11 limits. The fixture is regenerated only by a PR that says why; plot snapshots are
   thereby decoupled from upstream table changes.
2. **Snapshot tests** `tests/snapshot/test_plots.py`: one `@pytest.mark.mpl_image_compare(
   tolerance=2, savefig_kwargs={"dpi": 100})` test per plot instance (V03 and V13 locus plots: the
   first two instances only), on both fixtures; the figure comes from the registered plot function
   with `style.apply()`. Baselines under `tests/snapshot/baseline/` (PNG at 100 dpi to keep each
   < 300 kB).
3. **Report structure snapshot** `tests/snapshot/test_report_structure.py`: the ordered list of
   `(tag, id)` for every `<section>`, `<h2>`, `<h3>` and plot `<div>`/`<svg>` id of `report.html`
   built on the fixture, compared with `tests/snapshot/report_structure.json`; plus `compare.html`.
4. **CI.** One ubuntu-latest / Python 3.12 job runs `pytest tests/snapshot --mpl`; the existing
   matrix legs keep running `tests/snapshot` without `--mpl` (rendering only) (index decision 22;
   which job, Q1). `ci.yml` is not a protected file; `ci/gates.yaml` is unchanged.
5. **Approval packet** `scripts/dev/approval_packet.py` → `local/review/m9-approval/`: every
   baseline, the same plots on the K-12 run (M9-26), their deuteranopia simulations, the checklist
   table and score (M9-24), `report.html` and `compare.html`, and an index page. The lead posts the
   path and a summary to Khoa (chat and the PR); the PR carries `needs-khoa`.
6. **After approval:** the lead records Khoa's decision in the PR and PROGRESS.md, removes
   `needs-khoa`, merges. Any later baseline change is a reviewed PR whose description shows before
   and after images (§10.1).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| a plot issue merges after the baselines are approved | it adds its baseline in its own PR, marked `needs-khoa` for that plot only |
| Khoa rejects a plot | §14.4 loop: change, re-run checklist, new packet; after 3 rejected rounds on one plot, R9 (question for Khoa in PROGRESS.md) |
| macOS/3.11 render differs | not compared (decision 22) |
| fixture > 2 MB | shrink the simulated genome; never commit a file > 5 MB (R11) |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `pytest tests/snapshot --mpl` on ubuntu 3.12 with approved baselines | pass, < 120 s |
| 2 | change `style.DETECTION["detected"]` to another hex | V01, V02, V03, … snapshot failures (RMS > 2) |
| 3 | rename a report section id | structure snapshot fails |
| 4 | `pytest tests/snapshot` (no `--mpl`) on macOS | pass (figures render) |

Amendment (lead, 2026-10-24): example 2 checks V01, V02 and V16; V03 colours genes by abundance
and does not use `DETECTION["detected"]`.

## Acceptance criteria

- [ ] Baselines for all 16 plots (V03 and V13 locus instances included) on the rich fixture and the
  empty-state fixture; structure snapshots for both HTML files.
- [ ] CI runs the comparison on one leg; the snapshot layer stays under its 120 s budget.
- [ ] **Gate evidence:** Khoa's approval of the baselines and the report layout recorded in the PR
  and PROGRESS.md (date and words), the checklist score at approval, the packet path.

## Questions for the lead

1. Which job runs `--mpl`: add `tests/snapshot --mpl` to the ubuntu 3.12 `calibration-and-coverage`
   job, or add a 3.12 leg to the `test` matrix? **Proposed default:** a separate small `snapshot`
   job on ubuntu 3.12 (fastest feedback; no matrix growth on macOS).
2. Baseline DPI 100 rather than `report.dpi` 300 (size; RMS is computed on the baseline
   resolution). **Proposed default:** 100.
