# ADR-0043: Plot support files under `OUT/plots/support/`
- Status: accepted (lead)
- Date: 2026-10-22
- Spec sections: §7.9.4, §9.1, §9.2, §10.2 (V1, V2); docs/specs/M9-03-plot-support-run-inputs.md
- Needs Khoa: no (a non-breaking addition to §9.1; no gate, threshold or protected file changes)

## Context
V1 draws GC skew and oriC/terminus marks, and V2 draws oriC and terminus lines. None of these is in a
§9.1 file: GC skew needs the genome sequence, and the oriC call is in the intermediate
`features.parquet`. `giap report OUT` must rebuild every plot from `OUT/` alone (§7.9.4) and has no
genome argument (§5.1).

## Decision
1. A **support file** holds data that a plot needs, that no §9.1 file holds, and that only a step with
   the run's inputs (or a fitted model, M9-17) can compute. It lives under `OUT/plots/support/`.
2. Format: TSV per §9.2. A table larger than 2^20 bytes is written as `.tsv.gz` with gzip `mtime=0`,
   so the bytes are deterministic. Exactly one of `<name>.tsv` and `<name>.tsv.gz` exists; a rewrite
   removes the other.
3. Each support file is one Frictionless resource in `datapackage.json`, with the custom property
   `role: plot_support` and `schema: plot_support/<name>.json` (relative to `src/giap/schemas`).
4. Positions are 1-based as written, like every published table (spec §6.5).
5. The producing step writes it during `giap run`; `giap report` reads it and never writes or
   recomputes it. Missing support data gives the plot's empty state, not an error.
6. A support file is not a result table: it has no §9.2 test columns and makes no claims.

The first two support files are `replicon_landmarks` (oriC, terminus) and `gc_skew` (10 kb windows),
both written by the `plots` step (M9-03).

## Evidence (before → after)
- Before: V1/V2 could draw GC skew and oriC/ter only during `giap run`, so `giap report OUT` would
  produce different plots from the same run.
- After: both entry points read the same support files; M9-03's step tests pin that `giap report`
  leaves them byte- and mtime-identical.

## Consequences
- §9.1 gains an optional `plots/support/` directory; existing consumers are unaffected.
- M9-17 (plot support from fitted models) follows the same rules.
