# M11-11: the preprint outline (`docs/preprint/outline.md`)

Milestone M11 · **pass 2** · depends on: M11-03 (readiness rows, for evidence pointers), M11-08
(soft: user-guide page names) · estimate 0 code lines + ≈ 50 test lines + ≈ 250 prose lines ·
not `stats` (the statistician reviews the claims column) · roles: docs (writes it), validation
(evidence pointers), statistician (review), test-author (`tests/unit/test_m11_11_outline.py`) ·
Khoa: Q-K11 (venue, authors, AI statement).

## Spec sections implemented

- §13 row M11 "preprint outline". An **outline**, not a draft: headings, one-sentence intents,
  the figure and table plan, and a claim → evidence map. No results prose beyond what committed
  evidence states.
- R4: every quantitative claim names its error estimate (q-value, FDR, FWER, CI) or says
  `unestimated`. R11 and Khoa's data rule: aggregates only. §16 (limitations are what v1 defers).

## Definition

`docs/preprint/outline.md` (excluded from the site nav, kept in the repo; Q-L6 style):

1. **Front matter**: working title; venue (Q-K11 default: bioRxiv, Bioinformatics); authors
   `TBD (Khoa)`; a placeholder "AI-assistance statement (Khoa to write)".
2. **Sections**, each with 1–3 intent sentences: Abstract (structure only), Introduction (the
   question: what genome context explains about detection and abundance), Design (inputs, one gene
   table keyed by Bakta locus tags, the modules), Error control (placement FDR, gene-level FDR,
   circular-shift nulls, attack nulls — pointing at `docs/methods/`), Validation (simulation
   calibration; the Nissle 1917 pilot A1/A2/A3; TE1–TE8), Performance (measured runtimes),
   Availability (PyPI, Bioconda, containers by digest, Zenodo DOI, docs URL, licences), Limitations
   (§16 "Later" items; single sample; holdout not used; co-abundance AUROC deferred).
3. **Figure plan**: Fig. 1 workflow schematic (new, drawn by Khoa or from the docs); Fig. 2 V01
   circos of the pilot; Fig. 3 calibration (QQ / FDR-vs-q panels from the calibration evidence);
   Fig. 4 V07/V08 detection; Fig. 5 V11 operon pairs; Fig. 6 V13 novel ORFs; Supp: V16 compare.
   Each row: plot ID, run directory it comes from, the claim it supports.
4. **Claim → evidence map** (a table): claim (one line), number or `TBD`, error estimate, evidence
   (path in the repo: `validation/release/1.0.0/readiness.json › R1`, a `docs/methods/*.md`
   section, an ADR, `benchmarks/baseline.json` key), status (`supported` / `pending`).
   At least the §13.2 bullets and the §15.6 pilot bullets appear.
5. **Data and code availability** checklist: PXD078459 accession; genome accessions; what is not
   redistributed (RegulonDB, DIA-NN output rows) and why.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | the claim "≥ 99.5 % of A2 peptides map exactly" | row with evidence `validation/release/1.0.0/readiness.json › R1` and `validation/pilot/metrics/metrics.json › te1`, status `pending` until readiness exists |
| 2 | a claim "operon calls are calibrated" | error-estimate column `decoy FDR ≤ q (M6 gate)`, evidence `docs/methods/operons.md` + the M6 gate block in PROGRESS |
| 3 | a claim with no error estimate | test fails unless the estimate column says `unestimated` (R4) |

## Edge cases

- **Evidence path that does not exist yet** (readiness before M11-13): allowed only with status
  `pending`; a `supported` row with a missing path fails the test.
- **Numbers**: any number with 4+ digits in the outline must appear in a cited evidence file (the
  M10-15 notes rule), so the outline cannot invent results.
- **Row-level data**: no 7+ letter amino-acid strings (M10-05 guard regex).
- Genomic edge cases: not applicable.

## Acceptance criteria

1. The file has items 1–5; every §13.2 and §15.6 bullet is a row of the claim map.
2. `test_m11_11_outline.py`: the error-estimate rule, the `supported`-path rule, the number rule,
   the row-level guard.
3. The statistician's review comment on the claim map is on the PR.

## Estimate

0 code; ≈ 250 prose; tests ≈ 50.

## Questions for the lead (reversible defaults in use)

- Q-K11 (Khoa) venue, author list, AI statement, figure permissions. Defaults as in item 1.
- Is a preprint outline in the public repository acceptable before submission? Default: yes (the
  repository is public and the outline holds no unpublished numbers beyond committed evidence);
  alternative: keep it under `local/` and give Khoa the path.
