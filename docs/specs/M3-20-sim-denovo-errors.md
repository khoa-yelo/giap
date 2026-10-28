# M3-20: gippasim, de novo errors and reference peptides

Milestone M3 · depends on: M3-19 · changes `tests/sim/gippasim/` (new `denovo.py`) · simulator
tooling · pure Python · estimate ≈ 300 lines (not gated by ADR-0007)

## Spec sections implemented

- §11.3 "**De novo errors:** substitutions, mass-equivalent swaps and blocks (N ↔ GG, Q ↔ AG, …),
  I/L, adjacent swaps, reversed stretches, modifications, and random peptides. Scores are
  correlated with correctness. A `spectrum_id` pairs each de novo row with a reference
  identification."
- §5.3 / G2 reference peptides "from a same-strain database search on the same spectra" (used for
  c_ℓ, M3-11).
- §8.4 "De novo precision | known precision per bin" and "Placement FDR | … de novo error
  models"; §8.5 attack 5 "De novo errors that reverse a stretch of sequence".

## Behaviour

- For each spectrum of M3-19's correct and wrong identifications (per run), emit one de novo row
  (`source = denovo`, same `run`, `spectrum_id`) into `peptides_denovo.tsv` (generic table). With probability `p_correct(score)` the row is
  the true peptide (I and L re-drawn at random: "I/L"), else an **error** of one type, drawn with
  configurable weights: `substitution`, `swap` (adjacent), `block` (the five §7.2.4 mass-equivalent
  compositions, either direction), `near_block` (W ↔ AD/EG, R ↔ GV), `reversed_stretch` (3–6
  residues reversed), `modification` (A+Oxidation read as S, N+Deamidated read as D, G+Methyl read
  as A), `random` (a random 7–20-mer).
- **Scores** correlated with correctness: score = logistic-linked draw so that the true precision
  per score decile is known and written to `truth/denovo_precision.tsv` (`score_bin`,
  `precision_true`, `n`).
- The block list is gippasim's own copy of §7.2.4 (independent from `gippa.mapping.diffops`); a test
  checks both equal the spec list.
- `reference_peptides.tsv`: the database identifications of the same spectra (generic table,
  `source = database`, `search_db = same`, q-values from M3-19).
- **Truth:** `truth/denovo.tsv`: `run`, `spectrum_id`, `peptide`, `true_peptide`, `error_type`
  (`none` or one of the above), `correct` (bool: equal to the true peptide up to I/L).

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | true `PEPTNDEKAR`, error `block` (N → GG) | `PEPTGGDEKAR`, `error_type=block` |
| 2 | true `AKLMNPQRSTVW`, error `swap` at 4 | `AKLMPNQRSTVW` |
| 3 | true `PEPTIDEKAR`, I/L redraw | `PEPTLDEKAR` counted `correct = true` |
| 4 | true `MKAESTVLLR`, `reversed_stretch` of 4 at 2 | `MKTSEAVLLR` |
| 5 | true `PEPTADEKAR` with an oxidised A (A + 15.9949 = S) | observed `PEPTSDEKAR`, `error_type=modification` |
| 6 | 10⁴ spectra, configured precision 0.9 in the top decile | empirical precision of the top decile in `truth/denovo.tsv` within 0.9 ± 0.03 |
| 7 | every de novo row | its `(run, spectrum_id)` exists in `reference_peptides.tsv` when the spectrum had a correct database ID |
| 8 | `peptides_denovo.tsv` / `reference_peptides.tsv` | pass M1-06 validation with no error-level issues |

## Acceptance criteria (test author)

- [ ] `tests/unit/test_sim_denovo.py`: examples 1–8.
- [ ] Property: every error type, applied to a peptide, changes the sequence (except I/L) and keeps
  it within 7–40 aa; `block`/`near_block` errors change the monoisotopic mass by ≤ 0.05 Da.
- [ ] `truth/denovo_precision.tsv` matches the empirical precision per bin within 3 SE.

## Questions for the lead

1. Default error-type weights. **Proposed default:** substitution 0.3, block 0.2, swap 0.15,
   near_block 0.05, reversed_stretch 0.1, modification 0.1, random 0.1 (configurable; calibration
   tests sweep them).

## Lead decisions
1. Default error-type weights accepted as proposed (substitution 0.3, block 0.2, swap 0.15, near_block 0.05,
   reversed_stretch 0.1, modification 0.1, random 0.1); configurable, calibration sweeps them.
2. The test-author's API choices are accepted: `gippasim.denovo` (`DenovoConfig`, `ERROR_TYPES`,
   `DEFAULT_ERROR_WEIGHTS`, `apply_error`), switched on by `SimConfig.denovo` (requires `proteomics`, own RNG
   stream), deciles pooled over runs, one de novo row per non-decoy database row, `reference_peptides.tsv` =
   every non-decoy `peptides.tsv` row with `search_db = same`.
3. The calibration file is registered in the M4-05 SUITE_REGISTRY by the lead once both are on main.

## Lead decisions after code review (PR #66, changes requested — tests only)
1. Redraw rule: when the drawn type has no site, redraw among the remaining positive-weight types; only if
   none of them has a site, fall back to `substitution` (always has a site) even when its weight is 0 —
   documented here, tested.
2. Tests kill the surviving mutants: wrong-ID spectra's `true_peptide` is a correct peptide ≠ the DB peptide
   and varies across spectra; with weights {block: 0.5, swap: 0.5} `substitution` never appears; every
   stretch length 2–6 and both block orders occur; de novo `protein_ids` are NA; the de novo RNG stream is
   distinct from proteomics (one pinned hash of a de novo output).
