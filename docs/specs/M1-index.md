# M1 "Inputs and preflight": issue index

Spec: GIPPA_SPEC.md §13 row M1. Deliverables: `io/` readers and writers (FASTA, GFF3, BED, proBed,
VCF, TSV); generic peptide table; `gippa convert` (generic, DIA-NN from its documentation); every
preflight rule. **Gate:** every preflight rule unit-tested (`T-PRE-01 … T-PRE-25`, one pass and one
fail test each); `preflight.json` validated against `src/gippa/schemas/preflight.schema.json`.

One issue = one PR, ≤ ~400 changed lines excluding fixtures (AGENTS.md; `ci/gates.yaml` warns at 400,
fails at 800). Estimates include tests.

| Order | Issue | Scope | Depends on | Est. lines |
|---|---|---|---|---|
| 1 | [M1-01](M1-01-coordinates-model.md) | `model/coords.py`: `Interval`, `Replicon`, 0-based ↔ 1-based, origin-crossing split, sort key | none | 280 |
| 2 | [M1-02](M1-02-genome-fasta.md) | genome FASTA reader/writer: gzip by magic, first-token IDs, circular tags, safe-ID mapping, IUPAC facts | 01 | 320 |
| 3 | [M1-03](M1-03-protein-fasta.md) | protein FASTA: UniProt IDs, `A;B` groups, `*` handling, alphabet counts, length routing | 02 | 260 |
| 4 | [M1-04](M1-04-tsv-bed-probed.md) | TSV conventions (`%.6g`, `NA`, `true/false`, reproducible `.gz`), `read_tsv`, BED and proBed writers | 01 | 360 |
| 5 | [M1-05](M1-05-gff3-vcf.md) | GFF3 writer + minimal reader (round trip, `Is_circular`), VCF 4.3 writer (Appendix B.7) | 01, 04 | 380 |
| 6 | [M1-06](M1-06-peptide-table.md) | generic peptide table: typed read, validation `Issue`s, aggregation, `peptide_runs` | 04 | 380 |
| 7 | [M1-07](M1-07-convert-generic-diann.md) | modification stripping; `gippa convert --from generic\|diann`; `pending_sample` contract test | 04, 06 | 390 |
| 8 | [M1-08](M1-08-abundance-table.md) | abundance table: shape, numeric parsing, scale detection, missing counts | 04, 06 | 220 |
| 9 | [M1-09](M1-09-preflight-genome.md) | `Check` model; genome rules T-PRE-01…10; topology; replicon classification | 01, 02 | 400 |
| 10 | [M1-10](M1-10-preflight-evidence-rules.md) | inputs, proteins, peptides, abundance, evidence-tier rules T-PRE-11…22 | 03, 06, 08, 09 | 400 |
| 11 | [M1-11](M1-11-preflight-env-report.md) | environment rules T-PRE-23/24; `preflight.schema.json` (T-PRE-25); `preflight.json`/`.md` writers | 09, 10 | 360 |
| 12 | [M1-12](M1-12-preflight-cli.md) | `run_preflight`, normalised work files, `replicons.tsv`, `gippa preflight` CLI, exit codes 0/1/2/3 | 02–11 | 360 |

Total ≈ 3,900 lines. Parallelisable after M1-01: {02 → 03}, {04 → 05}, {04 → 06 → 07, 08};
09 needs only 01–02, so it can start alongside 04–08.

## Rule ID map (T-PRE)

| Range | Input | Issue |
|---|---|---|
| T-PRE-01 … 10 | genome: parse, IDs unique, IDs renamed, count, total length, contiguity, gaps, GC, topology, small replicons | M1-09 |
| T-PRE-11 | inputs: proteins or peptides present (§1.1) | M1-10 |
| T-PRE-12 … 15 | proteins: IDs, alphabet, length classes, duplicate sequences | M1-10 |
| T-PRE-16 … 18 | peptides: schema, content, evidence tier | M1-10 |
| T-PRE-19 … 22 | abundance: shape, join, values, scale | M1-10 |
| T-PRE-23 … 24 | environment: tools and databases (exit 2), resources | M1-11 |
| T-PRE-25 | `preflight.json` validates against its schema | M1-11, M1-12 |

## Decisions taken in these specs (for the lead to confirm)

- GFF3 **reader** is in M1 (minimal; needed for the §9.5 round trip and reused by M2).
- A shared `Issue` type (`src/gippa/io/issues.py`) is produced by io validators; preflight rules turn
  issues into `Check`s, so `gippa convert --from generic` and preflight share one validator.
- `preflight.json` lists every evaluated rule, including passing ones (`level: "ok"`).
- Coordinates in `preflight.json` are 1-based closed (human-facing, like TSV).
- `pending_sample` is a new pytest marker registered in `pyproject.toml` (not a protected file).
- Each issue file ends with "Questions for the lead" holding the proposed defaults; the most
  consequential are M1-06 Q1 (median across runs), M1-07 Q1–Q3 (DIA-NN abundance), M1-10 Q1–Q3
  (thresholds for decoy share and q range; abundance join before annotation) and M1-11 Q1–Q2
  (environment checks before M2 exists).
