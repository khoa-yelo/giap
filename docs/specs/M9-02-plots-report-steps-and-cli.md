# M9-02: The `plots` and `report` steps, a minimal static report, and `giap report`

Milestone M9 · **pass 1 — the report integration milestone** · depends on: M9-01 · estimate
≈ 380 code lines (+ ≈ 380 test lines, + nf-test) · not statistical

## The gap this closes

Without this issue, plots exist only as functions. It registers the two §6.2 steps, adds the
`REPORT` subworkflow at the end of `main.nf`, publishes `plots/` and `report/` under both executors,
implements `giap report OUT`, and renders a static `report/report.html` with all 13 §10.3
sections. A section shows every registered plot that rendered, and says why the others did not.
Plots added by later issues appear in the report with no change here.

## Spec sections implemented

- §6.2 rows "`plots` | py | tables → `plots/*` | core" and "`report` | py | tables, plots →
  `report/report.html` | core"; §6.1 data flow (plots after the tables).
- §7.9.4 "Rebuilds `plots/` and `report/` from existing tables without recomputing statistics.
  Running it twice MUST produce identical `plots/data/*.tsv`."
- §5.1 `giap report OUT [--no-interactive] [--plots V01,V02,...]`, shared options, exit codes 0, 1,
  4, 5. §5.2 `report.*`.
- §10.3 bullet 1 (one HTML file, Jinja2, inline CSS, no network), bullet 2 (the 13 sections),
  bullet 4 (top rows with a link to the full TSV); the interactive half is M9-14, the full content
  M9-16.
- §6.3 (process conventions, optional modules), §6.4 (executor equality), §9.1 (`plots/`,
  `report/`).

## Public interface

`src/giap/steps/plots.py`
```python
@step("plots", inputs=["preflight"], optional_inputs=[*ALL_PUBLISHED_KEYS, "events", "genome",
      "features", "replicons", "genes"], outputs=["plots"])
def plots_step(ctx: StepContext) -> None: ...

@step("report", inputs=["plots", "preflight"], optional_inputs=[*ALL_PUBLISHED_KEYS, "events"],
      outputs=["report"])
def report_step(ctx: StepContext) -> None: ...

def render_plots(view: RunView, cfg: GiapConfig, plots_dir: Path,
                 only: Sequence[str] | None) -> list[PlotEvent]
```
(`genome`, `features`, `replicons` and `genes` are used only by the support writers of M9-03 and
M9-18.)

`src/giap/report/build.py` (pure: data in, HTML string out)
```python
SECTIONS: tuple[Section, ...]   # 13 entries: id, heading, plot IDs, tables
def build_report(view: RunView, plots: Mapping[str, PlotFiles], skipped: Mapping[str, str],
                 cfg: GiapConfig) -> str
```
`src/giap/report/templates/report.html.j2`, `src/giap/report/templates/report.css` (inlined).
`workflow/subworkflows/local/report.nf` (`REPORT`, `GIAP_STEP as PLOTS; GIAP_STEP as REPORT_HTML`).

## Behaviour

1. **`plots` step.** `style.apply()`; build `RunView.from_inputs`; for each registered plot in ID
   order (or the `--plots` subset): skip check (M9-01 item 6) → call → `write_instance` for each
   instance. Events: `plots.rendered {plot, files, instances}` and `plots.skipped {plot, reason}`.
   Output `plots` is a directory: `Vnn_*.{svg,png,pdf}`, `data/`, `support/`.
2. **`report` step.** Reads `plots/` and the tables; writes `report/report.html` only. The
   template inlines the CSS; each rendered plot is embedded as its SVG text (the XML prolog
   stripped, IDs prefixed with the plot ID so two SVGs never share an `id`); each plot is followed
   by relative links to its `.png`, `.pdf` and data TSV (`../plots/…`). A skipped plot shows a
   one-line note with the reason in words ("not run: the `function` module is off").
3. **Sections** (fixed `id`s; §10.3 order): `summary`, `inputs-qc` (V17), `mapping` (V5, V6),
   `detection` (V7, V18), `genome-context` (V1, V2, V10, islands and regions tables, V3),
   `operons` (V11), `variants` (V12), `novel-orfs` (V13), `function` (V8, V14), `determinants`
   (V15), `sensitivity` (`tables/sensitivity.tsv`), `methods` (placeholder paragraph in pass 1),
   `provenance`. Every section always exists, so the structure snapshot is stable (M9-25).
4. **Tables in the report.** For each section's tables: the first 20 rows in the table's own sort
   order (§6.5), all columns, values as written in the TSV, and a link `../tables/<name>.tsv`. Absent
   tables → "not produced in this run" with the reason when known.
5. **Summary cards (pass-1 subset).** Each card shows a number and its error estimate or the word
   `unestimated` (§10.3, R4): genes detected (from `genes.tsv` `detected`, by `detection_tier`, with
   the tier FDR from `detection_summary` row `tier/q_threshold/<alpha.identification>/fdr`), peptides
   placed / input (with the mapping FDR at `stats.alpha.placement`), islands (count with
   `p_value ≤ stats.alpha.context`, "FWER-controlled"), flagged operon pairs (count with `q_value ≤
   stats.alpha.operon_flags`), novel ORFs called (with `class_fdr` or `unestimated`). Cards whose
   table is absent are not shown. M9-16 completes the set.
6. **Provenance (pass 1).** GIAP version (`giap.__version__`), the preflight input file names and
   SHA-256s, and the SHA-256 of `resolved_config.yaml` when readable. No timestamps, no absolute
   paths (§6.6).
7. **`giap report OUT`.** Resolves the config as other commands do; `RunView.from_outdir(OUT)`
   (exit 1 on `NotAGiapRun`); for the selected plot IDs, `remove_plot_files` first, so a plot that
   is now skipped leaves no stale files; runs `render_plots` then `build_report` in process;
   `--no-interactive` sets `report.interactive = false` for this call. Unknown `--plots` ID → exit
   1. A plot exception → exit 4 naming the plot. Prints the report path.
8. **`giap run`.** Local executor: both steps join the topological order (all their inputs are
   optional except `preflight`, so they run last); `PUBLISHED` gains `plots` → `plots/` and
   `report` → `report/`. Nextflow: `REPORT` takes one channel of `[key, file]` pairs mixed from the
   subworkflows' published outputs plus their `events` emissions, `collect()`s it, and runs
   `PLOTS` then `REPORT_HTML` (label `process_low`); `main.nf` calls it last. A subworkflow that did
   not run contributes nothing; the plot skips itself.
9. **Config.** No new keys. `report.formats`, `report.dpi`, `report.max_locus_plots`,
   `report.interactive` are read.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| `giap report` on a directory that is not a GIAP run | exit 1, "OUT has no tables/ or preflight/" |
| `giap report` on a run with `function` off | V14 note "not run: module `function` off"; exit 0 |
| every analysis table header-only (sim100k validation run) | each plot draws its empty state; exit 0 |
| `--plots V07` only | V07 files rebuilt; other plot files left as they are; report rebuilt from all present plot files |
| an `OUT/plots/` written by an older GIAP | rebuilt plots replace their own files only |
| `report.formats = [svg]` | the report still links PNG/PDF only when they exist |
| Nextflow run where OPERONS did not run | no operon keys in the channel; V11 skips with `input_absent` |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `giap run --executor local` on sim100k (giapsim seed 20261014) | `plots/` with every pass-1 plot that has input, `plots/data/*.tsv`, `report/report.html`; exit 0 |
| 2 | the same through Nextflow (`-profile test`) | identical `plots/data/*.tsv` bytes (executor equality) |
| 3 | `giap report OUT` twice on example 1's `OUT` | identical `plots/data/*.tsv` (§7.9.4) and identical `report.html` |
| 4 | `giap report OUT --plots V1,V99` | exit 1, "unknown plot ID V99" |
| 5 | `giap report /tmp/empty_dir` | exit 1 |
| 6 | report of example 1 parsed with `html.parser` | 13 `<section>` elements with the ids of item 3, in order |
| 7 | report of example 1 | no `src=`/`href=` starting with `http`, `//` or `file:`; every relative link resolves to an existing file |
| 8 | `giap report OUT --plots V17` after deleting `OUT/plots/V07_obs_vs_exp.svg` | V07 not rebuilt; its section notes "plot file absent" |
| 9 | nf-test stub of `REPORT` | `plots/` and `report/report.html` emitted |

## Acceptance criteria

- [ ] `tests/unit/steps/test_plots_step.py`: registry run with two fake plots (one renders, one
  skips): files, events, skip reasons; `--plots` subset; stale-file removal.
- [ ] `tests/unit/report/test_build_report.py`: examples 6, 7 on a tiny hand-built `OUT/`; cards
  omitted when their table is absent; every card carries an error estimate or `unestimated`.
- [ ] `tests/unit/test_cli_report.py`: examples 3, 4, 5, 8 (exit codes 0, 1, 4).
- [ ] `tests/integration/test_report_run.py`: example 1, and example 2 beside the existing
  executor-equality test.
- [ ] `workflow/tests/report.nf.test` (stub, example 9) and the registry key test for `plots` /
  `report` (as `test_nf_annotate_registry.py`).
- [ ] `giap report` no longer calls `_not_implemented`.

## Test ideas

- A property test: for any subset of registered plots, `render_plots` writes files only for that
  subset.
- Feed the report a plot SVG whose internal ids collide with another; assert prefixed ids are
  unique in the page.

## Hardening owed (pass 2)

- Plotly figures and sortable tables (M9-14, M9-15); full cards, methods, caveats, provenance
  (M9-16); the HTML-structure snapshot (M9-25); `-resume` across `REPORT` (M9-26).

## Questions for the lead

1. `--plots` accepts `V1` or `V01` (index decision 1). **Proposed default:** both, normalised.
2. The `plots` step takes **every** published key as an optional input, which makes it the
   largest fan-in in the workflow. Alternative: plots read `OUT/` after publishing. **Proposed
   default:** staged inputs (Nextflow's `publishDir` is asynchronous, so reading `OUT/` from a task
   races; and inputs keep the step a function of declared files, §6.2).
3. `giap report` writes into `OUT/plots` and `OUT/report` in place (the spec's wording "rebuilds").
   **Proposed default:** in place; no `--outdir` option (§5.1 has none).
4. Is a run whose report shows empty states for most plots (sim100k) an acceptable pass-1 done?
   **Proposed default:** yes; M9-13 records which plots rendered non-empty.
