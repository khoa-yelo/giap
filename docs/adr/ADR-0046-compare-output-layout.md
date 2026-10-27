# ADR-0046: `gippa compare` output layout (ADR-M9-B)
- Status: accepted (lead-level, 2026-10-22)
- Date: 2026-10-22
- Spec sections: §5.1, §7.9.3, §9.2, §9.5
- Needs Khoa: no (layout of new outputs; no default, threshold or protected file changes)

## Context
§7.9.3 names three comparison tables but places only `compare_detection.tsv` under `tables/`,
does not say how a gene x run matrix is laid out for any run count, how genes are joined, or
where the genome digest comes from before `run_manifest.json` existed (M9-22 spec, ADR-M9-B).

## Decision
1. All three tables live under `CMP/tables/`: `compare_detection.tsv`, `compare_abundance.tsv`,
   `compare_stability.tsv`, each with a Frictionless schema in `gippa/schemas/tables/`.
2. `compare_detection.tsv` is long, with `section` in {matrix, intersection}; per-run columns
   would make the schema depend on the run count. Intersections are exact (UpSet) sets.
3. Genes are joined by `locus_tag`; every run must have the same locus-tag set (exit 1 with the
   counts and the symmetric-difference size otherwise).
4. Genome identity: the `sha256` of the `role = genome` input of `run_manifest.json` when that
   file exists, else of `preflight/preflight.json`; no digest is exit 1 ("cannot verify the
   genome"); differing digests are exit 1 naming both labels and hashes.
5. The ADR number replaces the placeholder "ADR-M9-B" of the M9-22 spec (0044 and 0045 are taken
   on an open branch).

## Evidence (before → after: metrics, tests, simulations)
`tests/unit/test_cli_compare.py` (every exit-1 row, manifest over preflight, schema validation,
byte-identical reruns); `tests/unit/compare/`.

## Consequences
V16 and `report/compare.html` (M9-23) read these three files, plus each run's `tables/genes.tsv`
for the per-gene log2 values of the V16 panel (b) scatters (amended 2026-10-23, M9-23 lead decision 3;
the three tables hold no per-gene abundance). Runs annotated separately
(different locus tags) cannot be compared until a locus-tag mapping exists (backlog).
