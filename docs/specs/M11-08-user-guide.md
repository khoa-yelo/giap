# M11-08: the user guide and tutorials

Milestone M11 · **pass 2** · depends on: M11-07 (nav, generated references) · estimate 0 code lines
+ ≈ 60 test lines + ≈ 1,200 prose lines, as up to three PRs (get-started; user guide; tutorials) ·
not `stats`, but the statistician reviews `interpreting-results.md` · roles: docs (all pages),
test-author (`tests/unit/test_m11_08_docs.py`), statistician (review) · Khoa: none.

## Spec sections implemented

- §13 row M11 "MkDocs site"; §14.1 docs role: "user guide, tutorials, changelog".
- Appendix A.5: "Document every CLI option and output table from the code and schemas, not from
  memory. Include the licence notes for third-party tools (geNOMAD non-commercial, KEGG academic
  use, GPL tools run as separate programs)."
- R5 "geNOMAD is non-commercial: document that in the user guide"; §7 KEGG "free for academic use;
  the user guide says so"; R4 (every claim carries a q-value, FDR, FWER or `unestimated`).
- §5.1–5.4 (CLI, configuration, inputs, converters), §9 (outputs, conventions, flags), §10.3
  report, §7.9.3 `gippa compare`, `gippa makedb`, `gippa setup` (§3.2 sources, DB sizes §3.9), §6.4
  (local executor → Developer guide only), §13.2 runtimes (measured, from M10-12).

## Definition: pages and what each must contain

Paths: `docs/get-started/` (install, setup-databases, quickstart), `docs/user-guide/` (the rest,
including M10-12's `runtimes.md` and M11-09's `platforms.md`), `docs/tutorials/`,
`docs/developer.md`.

| Page | Must contain (checked by tests where marked ✓) |
|---|---|
| `install.md` | pip/uv install of `gippa` from PyPI; Bioconda (after M11-10); Nextflow ≥ the pinned version + Java 17+; profiles `apptainer` (tested), `conda`, `docker` (not tested, see platforms); Python ≥ 3.11 ✓ |
| `setup-databases.md` | `gippa setup` flags (`--bakta light,full`, `--genomad`, `--promotech-models`, `--go`, `--kegg`, `--cog`, `--contaminants`), sizes from §3.9, where `db_manifest.json` goes; RegulonDB is not downloaded by users ✓ (every flag exists in `cli.md`) |
| `quickstart.md` | one end-to-end command block on a tiny public genome with a peptide table, the expected output tree (§9.1), and how to open `report/report.html` |
| `inputs.md` | genome FASTA (circular tags), peptide table, proteins, abundance (scale detection rule of §5.3, one value column), annotation (Bakta JSON, GFF3, GBFF), reference peptides; plain NCBI FASTA preflight info (M2-10) |
| `convert.md` | `gippa convert --from diann|generic` (and which `--from` values exit 5) |
| `run.md` | the common options, `--modules`, `--publication`, `--seed`, `--resume`, `--strict`, exit codes (§5.1 table) ✓ |
| `outputs.md` | the §9.1 tree; conventions (§9.2: 1-based closed in TSV/GFF3/VCF, 0-based half-open in BED/proBed); link to `reference/tables.md` ✓ |
| `interpreting-results.md` | what q-values, FDR, FWER and `unestimated` mean in GIPPA; detection tiers; why absence of detection is not absence of expression; placement FDR vs. search-engine FDR |
| `plots-and-report.md` | the 16 plots (V1–V3, V5–V8, V10–V18) by ID, one line each, with the data TSV next to each ✓ |
| `compare.md`, `makedb.md` | as their commands; `makedb` writes a FASTA only, GIPPA never searches spectra (§16) |
| `licences.md` | NOTICE.md's table, **geNOMAD non-commercial**, **KEGG academic use**, GPL/AGPL tools run as separate programs, UniOP fetched at a pinned commit and never redistributed, DIA-NN is the user's own licence ✓ |
| `runtimes.md` | M10-12's table (measured, node class, date, commit), not §13.2's estimates |
| `tutorials/ecoli-k12.md` | K-12 from public inputs (genome accession, UniProt-derived or user peptides), each command, expected key numbers stated as ranges from committed evidence only |
| `tutorials/nissle-pilot.md` | the pilot (PXD078459) as a worked example from `validation/pilot/README.md`; aggregates only (Khoa's data rule) |
| `developer.md` | `--executor local` and its `annotation.precomputed` requirement, test layers (AGENTS.md), FarmShare notes, link to `release.md` |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `licences.md` | contains the words "non-commercial" next to geNOMAD and "academic" next to KEGG |
| 2 | `run.md` mentions `--threads` | `--threads` exists in `docs/reference/cli.md` under `gippa run` |
| 3 | `run.md` mentions `--executor local` | test fails: the local executor belongs in `developer.md` only (§6.4) |
| 4 | a tutorial states "about 35 min" | allowed only if `runtimes.md` (measured) supports it; the test fails on the literal "35 min" outside a quote of §13.2 |

## Edge cases

- **Option drift**: every `` `--option` `` written in a user-guide page must exist in `cli.md`
  (test). Options of other programs (`nextflow -resume`, DIA-NN flags) are written with their
  program name in the same code span and skipped.
- **Coordinates**: every page that shows a coordinate says which convention it uses; examples
  include one minus-strand gene and one feature crossing the origin of a circular replicon (§6.5:
  TSV and GFF3 keep `end > length`, with `Is_circular=true` in GFF3; BED and proBed write two records
  `_part1`/`_part2`; linear contigs never wrap). The page shows one row of each form.
- **Empty results**: `interpreting-results.md` explains an empty `novel_orfs.tsv` (header only) and
  a module switched off (`modules_off` in the manifest, no files).
- **Row-level pilot data**: never shown (no peptide sequences ≥ 7 aa; the M10-05 guard regex runs
  on `docs/tutorials/`).

## Acceptance criteria

1. Every page in the table exists and is in the nav; `mkdocs build --strict` passes.
2. `tests/unit/test_m11_08_docs.py`: the ✓ checks, examples 1–4, the option-drift rule, the
   local-executor rule, the row-level-data guard.
3. The statistician's review of `interpreting-results.md` is recorded as a PR comment.

## Estimate

0 code; ≈ 1,200 prose; tests ≈ 60.

## Questions for the lead (reversible defaults in use)

- **Gap: no public demo dataset.** The 100-kb test genome and its fixtures come from `gippasim`, which
  is test tooling and is not installed with `gippa`, so a user without the repository cannot run a
  quickstart on them. Default in use: the quickstart downloads the K-12 MG1655 genome (NCBI RefSeq,
  a §3.2 source) and runs `gippa preflight` and `gippa run` with the genome alone (`--peptides` is
  optional), showing the annotation and context outputs; the K-12 tutorial adds a user's own
  peptide table. Alternative for the lead: ship a small demo peptide table (in-silico tryptic
  peptides of a few K-12 proteins, < 100 KB) under `docs/tutorials/data/`.
