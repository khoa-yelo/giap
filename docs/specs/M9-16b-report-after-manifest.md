# M9-16b: `gippa run` re-renders `report.html` after writing `run_manifest.json`

Milestone M9 · **pass 2** · depends on: M9-16 (#228, open), MVP-03a (#196, merged); the KEGG half
of the acceptance also needs MVP-03c (the `databases[]` reader, not started) · estimate ≈ 70 code
lines (+ ≈ 220 test lines) · not statistical (no number is computed; the report is re-rendered
from files already in `OUT/`)

## The gap this closes

`gippa run` writes `run_manifest.json` in the run tail, after the graph has finished (ADR-0040,
MVP-03 D-B step 5). The in-graph `report` step (M9-02) runs before that, so it never sees the
manifest. Every pipeline report therefore prints "run_manifest.json not written by this GIPPA
version" and "KEGG retrieval date: not recorded in this run". Only a later `gippa report OUT` shows
the full provenance (M9-16 lead decision 11; backlog "M9-16b"). With this issue, `gippa run` under
both executors renders `OUT/report/report.html` again as the last step of the tail, after the
manifest. That render reads the plots already on disk and draws none, and it is deterministic. The
issue also states which file records which, so that the manifest and the report never need each
other's checksums.

## Spec sections implemented

- §10.3 section 13 "provenance (manifest digest)"; §1 rule 5 ("`run_manifest.json` records tool
  versions, database checksums, parameters, input SHA-256 hashes, seeds and per-step time and peak
  memory"); §6.7 (the manifest's contents); §6.2 rows `report` and `manifest`.
- §6.1 data flow "`plots ─► report.html ─► run_manifest.json`": this issue keeps that order and
  adds one more step, "`run_manifest.json ─► report.html` (re-render)" (see Q1).
- §7.9.4 "Rebuilds `plots/` and `report/` from existing tables without recomputing statistics":
  the tail re-render is the `report/` half of `gippa report`, and it leaves `plots/` alone.
- §6.4 (executor equality), §6.3 Error strategy (exit 4), ADR-0040 (run tail, exit rules),
  ADR-0017 ("Every run report states the KEGG retrieval date"), M9-16 behaviour 6 and lead
  decisions 1, 2, 7, 8, 11.

## Behaviour

1. **Run order of `gippa run`** (both executors; MVP-03 D-B with one step added):
   rotate events → preflight → `work/provenance.json` → graph (including the in-graph `plots` and
   `report` steps) → tail: datapackage (03b) → `run_manifest.json` → **report re-render** → exit.
   The re-render is the last write to `OUT/` outside `logs/`.
2. **When it runs.** Only when the manifest was written in this tail **and** its `status` is
   `completed` (exit 0 so far). It does not run after a failed graph (exit 4/5), after a failed
   manifest write, or after a preflight exit (1/2/3, when there is no tail). The condition is read
   from the manifest just written, so both executors use the same rule.
3. **What it does.** `rerender_report(OUT, cfg)` builds `RunView.from_outdir(OUT)`, reads the
   existing `OUT/plots/` with `read_plot_files`, and calls the M9-02/M9-16 `build_report` with the
   run's `cfg`. It writes `OUT/report/report.html` **atomically**: it writes
   `report/.report.html.tmp` in the same directory, then `os.replace`. It draws no plot. It never
   calls `render_plots` or `write_plot_support` and never writes under `plots/`. Under the
   Nextflow executor, `nextflow run` has already returned, so every `publishDir` copy has
   finished.
4. **One code path with `gippa report`.** `gippa report OUT` keeps its plot rendering, then calls the
   same `rerender_report`. The page written by the tail is therefore byte-identical to the page
   that `gippa report OUT` writes right after the run (M9-13 example 3 already asserts that
   `gippa report` is idempotent).
5. **Which file records which (the chicken-and-egg).** Dependencies go one way only: the report
   depends on the manifest, never the reverse.
   - The report prints the SHA-256 of the `run_manifest.json` file (M9-16 lead decision 7) and
     the manifest's `provenance_digest`. Nothing rewrites the manifest after the re-render, so the
     printed file digest equals `sha256(OUT/run_manifest.json)` once `gippa run` has exited.
   - The manifest records no checksum of `report/report.html` or of anything under `plots/`. The
     03a schema has no output-checksum field. The datapackage lists tables only. `tools[]` leaves
     out `gippa` itself (D-C). **Invariant for later work:** a future output-checksum field, in the
     manifest or as a datapackage `hash`, must leave out `report/` (its bytes depend on the
     manifest), or it must be added after the re-render with an ADR.
   - `steps[]` still records the in-graph `report` step (its status and timings), which really
     ran. The re-render is not a registered step. It is not in `steps[]`, has no seed entry, and
     emits no event to `logs/events.jsonl`, because that file is a manifest input and the
     manifest has already been written. It logs one INFO line to the `gippa` logger, and on
     failure it writes `OUT/logs/report_tail.log`.
   - The GIPPA version that renders the final page is the host `gippa run` process, which is the
     `gippa_version` that the manifest records. Under `-profile …,apptainer` the in-graph page may
     come from the container's GIPPA. The final page does not.
6. **Failure of the re-render** (after a completed graph and a written manifest): exit 4,
   `gippa run: report re-render failed; log: OUT/logs/report_tail.log`. Because the write is
   atomic, `OUT/report/report.html` keeps the in-graph page, and `run_manifest.json` is unchanged.
   This is the same rule as ADR-0040's "a tail error after a successful graph is exit 4" (see Q3).
7. **Work-directory copies stay as they were.** `OUT/work/report/report.html` (local) and the
   Nextflow task's `report/report.html` keep the in-graph page. `-resume` therefore still caches
   `REPORT_HTML`, and the tail re-renders again on every `gippa run`.
8. **Drop the #178 guards (conditional).** If #178 (M8-10) and #228 (M9-16) have both merged on
   `main` when this issue starts, and the guards are still there, remove
   `if "determinants_model" in PATHS` in `src/gippa/report/build.py` and
   `getattr(gippa.determinants, "CAVEATS", ())` in `src/gippa/report/content.py`, and use
   `gippa.determinants.CAVEATS` and the `determinants_model` key directly. Otherwise leave them
   alone (see Q4: the M9-16 review backlog gives this to the #178 merge instead).

## Public interface

`src/gippa/steps/plots.py`
```python
def rerender_report(outdir: Path, cfg: GippaConfig) -> Path:
    """OUT/report/report.html from OUT/ and the existing OUT/plots (no plot is drawn); atomic
    replace; returns the path. Shared by the gippa run tail and gippa report."""
```
`src/gippa/cli.py`: `_run_tail(...)` calls `rerender_report(outdir, cfg)` after `run_tail` when the
written manifest has `status == "completed"`. This needs the run's `cfg` in `_run_tail`: either a
new keyword argument or `base["cfg"]`, which `launch_sources` already carries. `gippa report` swaps
its `write_report(view, …)` call for `rerender_report(out, cfg)`. `write_report` stays for the
in-graph `report` step. No new CLI option, config key, schema field or Nextflow process.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| graph failed (exit 4) or internal error (exit 5) | manifest written with `status failed`; no re-render; exit code unchanged |
| manifest write failed | exit 4 "step manifest failed" (ADR-0040); no re-render; `report.html` is the in-graph page, which says "run_manifest.json not written by this GIPPA version" (see Q2) |
| preflight exit 1/2/3 | no tail, no manifest, no re-render (MVP-03 decision 12) |
| re-render raises | exit 4, `logs/report_tail.log`; `report.html` keeps the in-graph bytes (atomic replace); manifest unchanged |
| `OUT/plots/` absent or empty (e.g. every plot skipped) | rendered with `plots_dir = None` / empty-state slots, as `gippa report`; exit 0 |
| `OUT/report/report.html` absent after a completed graph | written fresh (the tail does not need the in-graph page) |
| `OUT` reused from an earlier run, new run fails | the earlier run's `report.html` stays next to the new failed manifest (OUT is never cleaned; out of scope, as today) |
| Nextflow `-resume`, everything cached | `PLOTS`/`REPORT_HTML` cached; the tail writes a new manifest (new timings) and re-renders, so the printed manifest digest changes and stays correct |
| a published `report.html` that is a symlink (non-`copy` publish mode) | `os.replace` replaces the link and never writes through it into `work/` |
| `report.interactive: false` / `--no-interactive` on `gippa run` | the re-render uses the run's `cfg`, so the page stays static (M9-14 once merged) |
| manifest `databases: []` (03a today: no db reader) | provenance shows the manifest but "KEGG retrieval date: not recorded in this run"; databases row "not recorded" (M9-16 degenerate row) |
| header-only tables, zero detected genes (empty inputs) | same page that `gippa report` gives (empty states); the re-render adds no condition of its own |
| circular/linear replicons, minus strand, ties | do not apply: no coordinate or ordering logic here; the page is whatever `build_report` gives on `OUT/` |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `gippa run --executor local` on the tiny sim fixture (the M9-02 example 1 setup), exit 0 | `OUT/report/report.html` provenance has GIPPA version, commit, a command line containing `--executor local`, seed; "manifest digest (SHA-256 of run_manifest.json)" equals `sha256(OUT/run_manifest.json)`; "provenance digest" equals the manifest's `provenance_digest`; the sentence "run_manifest.json not written by this GIPPA version" is **absent**; `OUT/work/report/report.html` still contains it |
| 2 | the same run with the Nextflow executor (`-profile test`) | the same assertions as example 1 on `OUT/report/report.html` |
| 3 | example 1's `OUT`, then `gippa report OUT` | `report.html` bytes identical to the tail's page |
| 4 | `rerender_report(OUT, cfg)` called twice on example 1's `OUT`, with every registered plot function monkeypatched to raise | both calls succeed and give identical bytes; the SHA-256 of every file under `OUT/plots/` is unchanged |
| 5 | the tail with `run_tail` stubbed to write a valid manifest whose `databases` has `{"name": "kegg_ko", "version": "2026-10-01", "retrieved": "2026-10-01", …}` and function tables present | the page prints "KEGG retrieval date: 2026-10-01" in provenance and in the function section |
| 6 | local run whose `operons` step raises (stub) | exit 4; `run_manifest.json` `status failed`; `rerender_report` not called (spy); no `logs/report_tail.log` |
| 7 | `run_tail` monkeypatched to raise after a successful graph | exit 4 with "step manifest failed"; `rerender_report` not called; `OUT/report/report.html` bytes equal `OUT/work/report/report.html` |
| 8 | `build_report` monkeypatched to raise only in the tail call | exit 4, stderr "report re-render failed; log: …/logs/report_tail.log"; `OUT/report/report.html` bytes equal the in-graph page; `run_manifest.json` bytes unchanged; no `.report.html.tmp` left |
| 9 | the manifest of example 1 | no key or value of `run_manifest.json` names `report/report.html` or a `plots/` file; the MVP-03 `stable_view` local-vs-Nextflow equality test still passes unchanged |
| 10 | `OUT/report/report.html` replaced by a symlink to a scratch file before `rerender_report` | afterwards `report.html` is a regular file and the scratch file's bytes are unchanged |

## Acceptance criteria

- [ ] `tests/unit/steps/test_rerender_report.py`: examples 4 and 10, atomicity (example 8's
  temp-file check), `plots_dir` absent.
- [ ] `tests/unit/manifest/test_manifest_cli_local.py` (extended) or a new
  `tests/unit/report/test_report_tail.py`: examples 1, 5, 6, 7, 8 under the local executor, and
  the run-order spy (manifest written before `rerender_report`, nothing written to `OUT/` after it
  apart from `logs/`).
- [ ] `tests/unit/manifest/test_manifest_cli_nextflow.py` (extended): the tail calls
  `rerender_report` after `run_tail` on exit 0, and not on exit 4 (stubbed `nextflow`).
- [ ] `tests/integration/test_report_run.py` (extended): examples 1, 2 (sim100k, both executors,
  on Slurm) and 3; example 9's no-report-hash check.
- [ ] `tests/unit/test_cli_report.py`: `gippa report` goes through `rerender_report`; its existing
  examples stay green.
- [ ] `grep -n '"determinants_model" in PATHS\|getattr(gippa.determinants' src/gippa/report` is
  empty **if** behaviour 8's condition held; otherwise the PR body says why the guards remain.
- [x] ADR-0040 gains an amendment (lead-level) for behaviours 1, 5 and 6 (Q1): Amendment 1.
- [ ] The backlog line "M9-16b" is closed; the KEGG end-to-end check (a real run with a `kegg*`
  database in `databases[]`) goes to the MVP-03c acceptance as a backlog line, if 03c has not
  merged.

## Questions for the lead

1. **ADR for the extra tail step.** §6.1's data flow ends at `report.html ─► run_manifest.json`.
   The re-render adds a step after the manifest, and the tail no longer contains only the
   manifest (ADR-0040 decision 1). **Proposed default:** a short "Amendment 1" to ADR-0040
   (lead-level, no Khoa: no output path, default or config key changes), written in this PR. It
   records the one-way rule (behaviour 5) and the exit-4 rule (behaviour 6).
2. **The sentence shown when the manifest is absent.** After this issue, "run_manifest.json not
   written by this GIPPA version" is wrong in the two cases where a user still sees it (the
   manifest write failed, or the page is the in-graph copy under `work/`). It is right only for
   `gippa report` on an `OUT` from before MVP-03a. **Proposed default:** keep the M9-16 constant
   (it is tested, and the role guard means changing it needs a test-author step). Add a backlog
   line to reword it to "run_manifest.json absent (see logs/manifest.log, or an older GIPPA
   version)" in the hardening pass. Alternative: change it here, with a test-author commit first.
3. **Exit code when the re-render fails.** **Proposed default:** exit 4, matching ADR-0040's tail
   rule. The in-graph page stays in place. Alternative: exit 0 with a warning, because every
   table and the manifest are complete and only the provenance section is stale. That would make
   a provenance gap silent, which goes against the spirit of R4.
4. **#178 guards: which change removes them?** The task brief and the backlog line "M9-16b" give
   the removal to M9-16b. The M9-16 review backlog says "Removing the two #178 guards belongs to
   the #178 merge into this branch, not M9-16b". **Proposed default:** behaviour 8. Whichever of
   #178 or #228 merges second removes them. M9-16b removes them only if they are still on `main`
   when it starts.
5. **The KEGG date cannot reach a real run yet.** `build_manifest` writes `databases: []` and
   `tools: []` until MVP-03c reads `db_manifest.json` and `versions.yml`. So after M9-16b, a real
   pipeline report shows the manifest digest, commit, command line, seed and inputs, but still
   prints "KEGG retrieval date: not recorded in this run" and an empty tools table. **Proposed
   default:** M9-16b tests the KEGG line with a stubbed manifest (example 5). The end-to-end KEGG
   date becomes an MVP-03c acceptance item. Should MVP-03c be scheduled before M9-25's approval
   packet, so that Khoa reviews a report with a full provenance section?
6. **Considered alternative: stage `work/provenance.json` instead.** ADR-0040 notes that the
   in-graph `report` could stage the provenance core. The core holds the version, commit, config,
   inputs, databases (so the KEGG date, once 03c lands) and seed. It does not hold the command
   line, `tools[]`, step timings, or the manifest file's SHA-256 (M9-16 lead decision 7).
   **Proposed default:** re-render as decided in M9-16 lead decision 11. Staging the core would
   still leave those fields missing and would add a `provenance` input to `REPORT`.

## Line estimate

| Piece | Code lines |
|---|---|
| `rerender_report` (view, plot files, atomic write) | ≈ 25 |
| `_run_tail` hook: condition, cfg threading, failure log and exit 4 | ≈ 30 |
| `gippa report` switched to `rerender_report` | ≈ 5 |
| #178 guard removal (conditional) | ≈ −4 / +2 |
| docstrings and ADR-0040 amendment (docs, not counted) | — |
| **Total** | **≈ 60–70** (tests ≈ 220) |

## Lead decisions (2026-10-22)

1. **ADR-0040 Amendment 1** (lead-level, written in the implementation PR): the run tail ends
   datapackage → run_manifest.json → report re-render; dependencies run one way (the report reads
   the manifest, the manifest records nothing under `report/` or `plots/`); a failed re-render
   exits 4.
2. **Absent-manifest sentence.** M9-16's tested text stays; rewording goes to the backlog (it needs
   a test-author step first).
3. **Exit 4** on a failed re-render, keeping the in-graph page (the proposed default).
4. **#178 guards.** Whichever of #178 / #228 merges second removes them (the lead does it in that
   merge); M9-16b removes them only if they are still on main.
5. **KEGG date.** Tested with a stubbed manifest here (example 5); the end-to-end KEGG check moves
   into MVP-03c's acceptance. MVP-03c is scheduled before the M9-25 approval packet, so Khoa
   reviews a report with full provenance.
6. Staging `work/provenance.json` into the in-graph step is rejected (the proposed default).
7. M9-16b starts after #228 merges; feature loop test-author → implementer → code review.
8. **Test-author (3d6222e9): gaps accepted.** The implementer starts only after #178 and #228 merge
   and the lead merges origin/main into this branch; the merge brings #236's `databases[]` reader
   and M9-14's twins. `rerender_report` builds through the same twins path as `write_report`, and
   its size warning does not reach events.jsonl. Example 2 (Nextflow) runs in the implementer's
   Slurm check. The 16 failing test_cards/test_caveats tests wait on #178.
