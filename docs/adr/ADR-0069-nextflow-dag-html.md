# ADR-0069: Nextflow DAG HTML in OUT/logs
- Status: accepted-by-khoa (Khoa, 2026-10-29, in session: "Yes add with-dag as well")
- Date: 2026-10-29
- Spec sections: §6.3 (the `nextflow run` command line); docs/specs/M2-06-nextflow-skeleton.md
- Needs Khoa: answered 2026-10-29 (§6.3 pins the command line)

## Context
§6.3 fixes the launch command with `-with-trace`, `-with-timeline` and `-with-report`. Every run
therefore writes `OUT/logs/trace.tsv`, `timeline.html` and `nf_report.html`, but no diagram of the
run itself. The curated process diagram (docs/user-guide/run.md, "Pipeline processes") is drawn
for every module on; a user's run may switch modules, inputs and tools off.

## Decision
1. `launch_command` adds `-with-dag OUT/logs/dag.html` after `-with-report`. Nextflow 26.10
   writes the DAG of that run as an HTML page with a Mermaid diagram (processes, operators and
   the channels between them).
2. `workflow/nextflow.config` sets `dag.overwrite = true`, as for the trace, timeline and report,
   so a `-resume` into the same OUT replaces the file instead of failing.
3. The docs list `dag.html` among the logs and point to the curated diagram for reading it.

## Evidence (before → after: metrics, tests, simulations)
- tests/unit/test_runner_nextflow.py pins the new command line (spec M2-06 example 1 + the flag).
- No output table, statistic or manifest field changes.

## Consequences
- One more HTML file in OUT/logs (tens of kB). The raw DAG draws Nextflow operators between
  processes, so it is busier than the curated diagram.
- `src/gippa/runner/` and `workflow/` are on the never-reusable list (ADR-0060 §4): merged before a
  release, it needs that release's evidence to be re-run.
