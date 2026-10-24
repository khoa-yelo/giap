# ADR-0040: The manifest runs after the graph; finalize is split between the graph and the run tail
- Status: accepted (lead-level, 2026-10-22)
- Date: 2026-10-22
- Spec sections: §6.2 (rows `finalize`, `manifest`), §6.3 ("one process per Python step", Error strategy), §6.7, §9.1, §9.5; docs/specs/MVP-03-manifest-finalize.md (D-A, D-B, lead decisions 2, 3 and 13)
- Needs Khoa: no. No default, threshold, output file, config key or protected file changes. The run order of `giap run` gains a tail; every §9.1 output keeps its path.

## Context
§6.2 lists `finalize` and `manifest` as ordinary steps, and §6.3 runs every Python step as a Nextflow process. Three facts make that impossible for the manifest and awkward for finalize:
- The Nextflow `trace.tsv` is complete only after Nextflow exits, and the manifest's `steps[]` and `tools[]` come from it.
- A failed run must still get a manifest (§6.3 exit 4; MVP-03 example 11). A process downstream of the failed task never runs.
- The set of tables that a run published, which the Frictionless descriptor must list, is known only once the graph has ended.

M9's in-graph `plots` and `report` processes (M9-02 item 8) read `tables/genes.tsv`, so that table cannot wait for the tail.

## Decision
1. **`manifest` runs in the run tail.** The step is registered (`giap _step manifest` works on files alone), but neither executor plans it: `giap.runner.local.TAIL_STEPS = ("manifest",)`, and `workflow/` has no MANIFEST process. `giap run` calls the tail under both executors after the graph, also after a failed step. The order is: rotate `logs/events.jsonl` → preflight → write `OUT/work/provenance.json` (the provenance core) → graph → tail. A preflight exit (1, 2, 3) writes no manifest. A tail error after a successful graph is exit 4 (`step manifest failed; log: OUT/logs/manifest.log`). After a failed graph a tail error is only logged and the original failure is kept.
2. **`finalize` is split.** The in-graph `finalize` step owns `tables/genes.tsv`, and `context_stats` drops its `genes_tsv` output (MVP-03b). The descriptor half of the §6.2 row (`OUT/schemas/datapackage.json` plus the schema copies) is a function of the run tail, run before the manifest.
3. The manifest schema keeps its repository name `src/giap/schemas/run_manifest.schema.json`. §9.5 calls it `manifest.schema.json`, but the repo and `tests/unit/test_schemas.py` use the former name (lead decision 13).

Alternative rejected: run `finalize`, `plots`, `report` and `manifest` all in the tail. That is simpler and executor-identical by construction, but it overturns M9-02's accepted `report.nf` design.

## Evidence
- MVP-03a tests: `tests/unit/manifest/test_manifest_cli_local.py` and `test_manifest_cli_nextflow.py` pin the run order, the rotation, the failure exits and the "no manifest after preflight" rule under both executors. `tests/integration/test_mvp03_manifest_executors.py` (example 10) checks that `stable_view` is equal across executors. `test_mvp03_manifest_resume.py` (examples 11 and 12) covers the failed run and `-resume`.

## Consequences
- §6.3's "one process per Python step" has one exception, `manifest`. The tail is identical under both executors because it runs in the `giap run` process.
- The manifest's `provenance_digest` is the SHA-256 of `OUT/work/provenance.json`, written before launch, so the in-graph `report` (M9-16) can stage it and print the same digest (lead decision 6).
- MVP-03b adds `finalize.nf` and the datapackage writer of the tail. Later tail work (for example the pass-2 Frictionless validation of published tables) goes in the same tail.

## Amendment 1 (M9-16b, 2026-10-23; lead-level, no output path, default or config key changes)
- **The tail re-renders the report after the manifest.** After `run_manifest.json` is written
  with `status == "completed"`, `giap run` (under either executor) re-renders
  `OUT/report/report.html` from `OUT/` and the existing `OUT/plots/`. It draws no plot, and it
  writes atomically (`report/.report.html.tmp`, then `os.replace`). `giap report` uses the same
  `rerender_report`. The tail is therefore no longer the manifest alone (decision 1).
- **One-way dependency.** The report depends on the manifest, never the reverse.
  - The report prints the manifest's file SHA-256 and its `provenance_digest`.
  - The manifest records no checksum of `report/` or `plots/`.
  - A future output-checksum field must leave out `report/`, or be added after the re-render
    with an ADR.
  - The re-render is not a registered step: it does not appear in `steps[]` and emits no event
    to `logs/events.jsonl`.
- **Exit 4 on failure.** A re-render failure after a completed graph and a written manifest is
  exit 4 (`logs/report_tail.log`): the same rule as a tail error after a successful graph.
  `report.html` keeps the in-graph page and the manifest is unchanged.
- **Skipped when the run did not complete.** There is no re-render after a failed graph (exit
  4/5), a failed manifest write, or a preflight exit (1/2/3). Work-directory copies keep the
  in-graph page, so `-resume` still caches `REPORT_HTML`.
