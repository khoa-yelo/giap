# M9-18: V12 `variants`: where do peptides support ALT residues?

Milestone M9 · **pass 2** · depends on: M9-03 (support writer), M9-10 (shared locus helpers); PRs
#164 (variants step) and #175 (variant flags) merged · estimate ≈ 320 code lines (+ ≈ 300 test
lines) · not statistical

## Spec sections implemented

- §10.2 V12: "(a) counts by call × confidence × flag; (b) peptide-vs-genome text alignments at up
  to 30 top ALT sites"; static SVG + table.
- §7.5 (REF/ALT/mixed calls, `high` / `candidate`, the §7.5.3 flags `mass_equivalent`, `snv`,
  `multi_nt`, `strain_difference`, `database_derived`, `start_codon`, `stop_codon`, class error
  rates), §9.3 `variants.tsv` (branch `origin/m7/variant-flags`: `call` enum{ALT, mixed},
  `alt_type` enum{substitution, swap, block, frameshift}), `peptide_hits.tsv`; §6.5.

## Support: `plots/support/variant_windows.tsv` (written by the `plots` step under `giap run`)

The reference protein around each of the top-30 sites is not in `OUT/` (`genes.parquet` with its
`aa` column is an intermediate). The `plots` step, given `genes` (`genes.parquet`) as an input,
writes `variant_id`, `gene`, `aa_pos`, `window_start_aa` (1-based), `ref_window` (the gene's
protein from `aa_pos − 10` to `aa_pos + 10`, clipped to the protein). Under `giap report` the file
is read back; when absent, panel (b) shows peptides only.

## Panel (a)

Grouped bars: x = `call` (ALT, mixed), groups = `confidence` (high, candidate), one dodged bar per
flag (*amended 2026-10-24, lead, #322:* stacking counted a multi-flag site several times). A
site with several flags is counted once in each flag and once in "all"; sites without flags count
under "no flag". The "all" total is printed above each bar. Colours: Okabe–Ito in the §9.4 flag
order.

## Panel (b): text alignments

Top sites: `call = ALT`, ordered `confidence` high first, then `n_alt` descending, then
`class_error_rate` ascending (NA last), then `variant_id`; at most 30. For each, monospace text:

```
GSIM_00105 aa 214  Q→K  substitution  high  n_alt 3 / n_ref 0  class error 0.004
ref   …LVAEQTRGLK…
pep   …LVAEKTRGLK…   (3 peptides)
```

Peptides are the accepted `peptide_hits` rows covering the codon (same contig and strand, codon
inside `[start, end]`); each distinct peptide is placed by its offset from the gene start
(in-frame), the ALT residue in bold `#D55E00`. `alt_type` ∈ {block, frameshift}: peptides listed
unaligned with the note "length-changing: shown unaligned". Lines are truncated to 21 residues
around the site.

## Data TSV

- `V12_variants.tsv`: `panel` = a: `call`, `confidence`, `flag`, `n`; `panel` = b: `rank`,
  `variant_id`, `gene`, `aa_pos`, `ref_aa`, `alt_aa`, `alt_type`, `confidence`, `n_alt`, `n_ref`,
  `class_error_rate`, `ref_window`, `peptides` (list). Sorted by panel then the orders above.

## Title

"`<h>` high-confidence and `<c>` candidate ALT sites (class error `<min>`–`<max>`)", or
"… (class error unestimated)" when every `class_error_rate` is NA; "No ALT site supported by
peptides" when the table is header-only. "ALT site", never "mutation" or "modification" (§7.5
never names modifications).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| `variants.tsv` absent (variants step not run) | skip `input_absent` |
| header-only | empty state |
| a site on the `−` strand | the window is the protein sequence (already in transcription order); peptide offsets from the gene's protein start |
| a site in a gene crossing the origin | no special case: positions are protein positions |
| a site near the protein start (`aa_pos` 3) | window clipped at 1; `window_start_aa` 1 |
| `--genome-relation same` (no ALT calls; QC only) | empty state; title "No ALT calls (genome relation `same`)" |
| > 30 ALT sites | 30 drawn; all in panel (a) |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | v1 ALT high flags `snv`; v2 ALT high flags `snv;start_codon`; v3 mixed candidate no flags | (a): ALT/high/snv 2, ALT/high/start_codon 1, ALT/high/all 2, mixed/candidate/no flag 1, mixed/candidate/all 1 |
| 2 | v1 n_alt 3, v2 n_alt 5, v4 ALT candidate n_alt 9 | (b) order v2, v1, v4 |
| 3 | gene protein of 300 aa, `aa_pos` 214, Q→K | `ref_window` = residues 204–224; `window_start_aa` 204 |
| 4 | v5 `alt_type` frameshift | (b) row present; peptides unaligned |
| 5 | header-only `variants.tsv` | empty state |

## Acceptance criteria

- [ ] `tests/unit/viz/test_v12_variants.py`: examples 1–5; strand handling; truncation.
- [ ] `tests/unit/viz/test_support.py` extended for `variant_windows` (clipping, `−` strand).
- [ ] Coded against `origin/m7/variant-flags` `variants.schema.json`; re-checked at merge.

## Questions for the lead

1. Top-site order (above): **Proposed default:** confidence, n_alt, class error, ID.
2. Window ±10 aa: a display choice. **Proposed default:** ±10.

## Lead decisions (2026-10-22)

1. Top-site order: confidence, then n_alt, then class error, then ID (the proposed default).
2. The display window is ±10 aa (the proposed default).
3. Feature loop: test-author → implementer → code review (no stats label).

## Lead decisions on the tests (2026-10-22, after e01a6a9)

4. Title counts and the class-error range cover `call == "ALT"` sites only (as M9-16 decision 5).
5. The `same` title reads "No ALT calls (genome relation same)" — no literal backticks.
6. A top site whose gene is missing from `genes.parquet` or has an empty `aa` keeps its row with
   an NA `ref_window` (and NA window start); panel (b) shows it unaligned.
7. A covering peptide must be in frame with the codon; out-of-frame placements do not cover.
8. Panel (b) prints one peptide line per site plus "(N peptides)", the distinct-peptide count.
9. A header-only `variants.tsv` gives a header-only `variant_windows.tsv` (schema-valid).
10. Decisions 6–9 are implemented but untested here; the hardening pass adds tests (backlog).

## Lead decision on the review (2026-10-22, PR #227 B1)

11. A site with no gene (`gene` / `aa_pos` NA: sixframe or novel-ORF peptides) has no
    `variant_windows` row and is shown unaligned in panel (b) with "(no gene: not aligned)".
    Lead fix with tests in tests/unit/viz/test_v12_review_227.py.
