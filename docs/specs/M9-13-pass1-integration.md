# M9-13: Pass-1 closer: every available plot in `report.html` from a real sim100k run

Milestone M9 · **pass 1 — the M9 integration milestone** · depends on: M9-02 … M9-12 · estimate
≈ 100 code lines (+ ≈ 250 test lines) · not statistical

## The gap this closes

Each plot issue is tested on in-memory frames. This issue proves the whole chain on a real run: a
`giap run` on sim100k through Slurm and Nextflow writes `plots/` and `report/report.html` with
every pass-1 plot rendered or skipped for a recorded reason; the local executor gives the same
data TSVs; `giap report` rebuilds them identically; the report opens offline. It is the MVP
evidence before any hardening.

## Spec sections implemented

- §7.9.4 (twice → identical `plots/data/*.tsv`), §6.4 (both executors identical), R7/TE7 (the
  determinism rule, extended to `plots/data/`), §10.3 tests "A link checker finds no broken
  relative links. A test asserts that no external URLs are loaded at view time."

## Work

1. `tests/integration/test_plots_report_run.py`: local-executor `giap run` on sim100k (the
   giapsim default seed), then: every registered pass-1 plot is either rendered or has a
   `plots.skipped` event; `giap report OUT` twice → identical `plots/data/*.tsv` and identical
   `report.html`; executor equality against the Nextflow `-profile test` run (beside the existing
   equality test).
2. `src/giap/report/checks.py` (pure; used by tests and later by M9-14/M9-25):
   `relative_links(html) -> list[str]`, `external_loads(html) -> list[str]` (any `src`, `href` on
   `<link>`/`<script>`/`<img>`/`<iframe>`, CSS `url(…)` or `@import` that is absolute, protocol-
   relative or `file:`). Anchor `href`s to documentation are allowed only if M9-16 needs them and
   are not "loaded at view time".
3. **Evidence that counts:** one real `sbatch scripts/slurm/giap_run.sbatch` on sim100k
   (`-profile farmshare,apptainer`) reaching the end with `report/report.html`; job ID in
   PROGRESS.md; a short note listing which plots rendered non-empty, which drew empty states, and
   which skipped with which reason.
4. A preview copy for Khoa under `local/review/m9-pass1/` (gitignored; not the M9-25 approval).
5. Any wiring gap found (a key missing from `REPORT`'s channel, a published file missing from
   `RunView`'s map) is fixed here.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| the sim100k run with zero detections (as on 2026-10-22) | report written; V7/V18/V10 empty states; recorded in the note |
| a pass-2 plot's table already merged (e.g. #164) | its plot still absent (not registered yet); no error |
| Nextflow run after a module is switched off | the corresponding plots skip with `module_off` |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | sim100k, local executor | `plots/V01_circos.svg`, `V02`, `V03_locus.tsv` (index), `V05`, `V06`, `V07`, `V08`, `V10`, `V11`, `V13`, `V17`, `V18` + data; `report/report.html` |
| 2 | `external_loads(report.html)` | `[]` |
| 3 | every `relative_links` target | exists on disk |
| 4 | `diff -r plots/data` between local and Nextflow runs | no difference |
| 5 | `giap report OUT` twice | identical `plots/data/*.tsv` |

## Acceptance criteria

- [ ] The integration test (examples 1–5) passes on FarmShare (`tests/integration`, run in a Slurm
  job via `test_branch.sbatch`, not in the session allocation).
- [ ] `tests/unit/report/test_checks.py` covers the link and external-load detectors on crafted
  HTML (absolute `src`, `//cdn`, CSS `url(http…)`, `file:` and data URIs allowed).
- [ ] The Slurm job ID and the rendered/empty/skipped list in PROGRESS.md.

## Questions for the lead

1. Should M9-13 also run on the K-12 full-DB recording to get non-empty plots? **Proposed
   default:** no; that is M9-26 (pass 2), as M8 deferred its meaningful real-data run.
