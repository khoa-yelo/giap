# M6-07: Promoter and terminator evidence, the seven flags, and `operons.tsv`

Milestone M6 · track R · **pass 1** · depends on: M6-05, M6-06 (and M6-04 for the evidence columns)
· estimate ≈ 220 code lines (+ ≈ 300 test lines) · not statistical, but it is where the module's
claims become words

## Spec sections implemented

- §7.4.5 in full: the "yes" definition ("a Promotech score ≥ `annotation.promoters.promotech_yes`
  (0.5), or a TerminatorNet terminator at or above the calibrated confidence cut-off, between the
  two genes"), the six-row flag table, the four abundance-pattern definitions, "Any other
  combination for a discordant or co-abundant pair is labelled `other`", and the caveat that
  "Promotech and Promoter Calculator see σ70 promoters only, and TerminatorNet and TransTermHP see
  intrinsic terminators only. A missing call therefore never proves a missing element."
- §7.4.7: `operons.tsv` ("per TU: genes, size, detected fraction, median log2 abundance, maximum
  |Δ| between adjacent members, flag counts") and `operon_flags.tsv` ("flagged pairs only, with
  evidence columns").
- §9.3 the column lists of both tables; §9.4 flag scope "Operons".

## Public interface

`src/giap/operons/flags.py`, `src/giap/operons/tu_table.py` (pure)
```python
FLAGS = ("internal_promoter_supported", "internal_terminator_supported", "unexplained_discordance",
         "conflict_terminator_in_operon", "novel_operon_candidate", "coabundant_despite_terminator",
         "other")

def evidence_yes(pairs: pd.DataFrame, cfg: GiapConfig) -> pd.DataFrame: ...   # two bool columns
def abundance_pattern(pairs: pd.DataFrame, alpha: float) -> pd.Series: ...    # see step 2
def flag_pairs(pairs: pd.DataFrame, cfg: GiapConfig) -> pd.DataFrame: ...     # operon_flags.tsv
def operons_table(tus: pd.DataFrame, pairs: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame: ...
```

1. **Evidence "yes".** `promoter_between_yes = promoter_between >= annotation.promoters.
   promotech_yes`; `terminator_between_yes = terminator_between >= annotation.terminators.
   terminatornet_yes` (or `/100` when `--terminators transtermhp`, as M2-15 already normalises).
   NA evidence → `false` **and** the pair is excluded from flagging, with `flags` gaining
   `layer_missing:<layer>`: a missing layer is not evidence of absence.
2. **Abundance pattern** (α = `stats.alpha.operon_flags`, 0.1):
   `downstream_higher` — OP, discordant (`q ≤ α`), Δ > 0, or one-detected with the **downstream**
   gene detected; `upstream_higher` — the mirror; `discordant` — OP with `q ≤ α` in either
   direction; `concordant` — OP, both detected, `q > α`; `coabundant` — NO with `q ≤ α` from
   M6-06; otherwise `none` (no flag row).
3. **Flagging**, first match wins, in this order:

   | # | class | promoter yes | terminator yes | pattern | flag |
   |---|---|---|---|---|---|
   | 1 | OP | true | false | downstream_higher | `internal_promoter_supported` |
   | 2 | OP | false | true | upstream_higher | `internal_terminator_supported` |
   | 3 | OP | false | false | discordant | `unexplained_discordance` |
   | 4 | OP | any | true | concordant | `conflict_terminator_in_operon` |
   | 5 | NO | false | false | coabundant | `novel_operon_candidate` |
   | 6 | NO | any | true | coabundant | `coabundant_despite_terminator` |
   | 7 | OP or NO | any | any | discordant or coabundant | `other` |

   A pair whose pattern is `none` gets **no** row. CV, DV and SW pairs are never flagged.
4. **`operon_flags.tsv`** is `pairs.tsv` filtered to flagged rows, with the §9.3 columns:
   `pair_id`, `gene_a`, `gene_b`, `uniop_class` (= the pair `class`), `promoter_between_yes`,
   `terminator_between_yes`, `abundance_pattern`, `flag`, `q_value`, `lfdr`, `evidence_note`.
   `evidence_note` is a short generated sentence naming the evidence actually used, e.g.
   `"promoter 0.81 between; no terminator; q=0.02; downstream higher by 1.4 log2"`; it never
   asserts absence of an element (the §7.4.5 caveat).
5. **`operons.tsv`.** One row per TU of C5 (every TU, including single-gene TUs):
   `tu_id`, `contig`, `start`, `end`, `strand`, `genes` (list, transcription order), `n_genes`,
   `detected_frac` = detected / eligible genes, `median_abundance_log2` = median over the detected
   genes with a value, `max_abs_delta_log2` = max |Δ| over the TU's internal adjacent pairs with
   both genes detected, `n_flags` = flagged pairs whose **both** genes lie in the TU.
6. **Sorting.** `operon_flags.tsv` by `pair_id`; `operons.tsv` by contig (C1 order) then `start`
   (§6.5).
7. **No RNG, no new statistics.**

## Degenerate cases

| Situation | Behaviour |
|---|---|
| promoter layer missing | no `internal_promoter_supported` ever; the run flag already says why; the report text must repeat it |
| terminator layer missing | no row 2, 4 or 6; rows 3 and 5 still possible (they need "no") — **this is the risk** (Q1) |
| one-detected discordant pair | direction from which gene is detected; `conflict_terminator_in_operon` is impossible (it needs both detected) |
| single-gene TU | `operons.tsv` row with `max_abs_delta_log2` NA, `n_flags` 0 |
| TU with no detected gene | `detected_frac` 0, `median_abundance_log2` NA |
| TU spanning the origin | coordinates unwrapped (`end > length`), §6.5 |
| pair whose genes lie in two different TUs (an NO pair, by construction) | counted in neither TU's `n_flags`; it appears in `operon_flags.tsv` only |
| every test `unestimated` (sim100k) | `operon_flags.tsv` is a header-only file; `operons.tsv` still has one row per TU |

## Examples

OP pair with q = 0.02: Δ = +1.4, promoter 0.81, terminator 0.10 (cut-off 0.5).

| # | Input | Expected flag |
|---|---|---|
| 1 | the pair above | `internal_promoter_supported` (row 1) |
| 2 | same but Δ = −1.4, promoter 0.10, terminator 0.80 | `internal_terminator_supported` (row 2) |
| 3 | same, promoter 0.10, terminator 0.10, Δ = +1.4, q = 0.02 | `unexplained_discordance` (row 3) |
| 4 | OP, both detected, q = 0.60, terminator 0.80 | `conflict_terminator_in_operon` (row 4) |
| 5 | NO pair, q = 0.04 from M6-06, promoter 0.10, terminator 0.10 | `novel_operon_candidate` (row 5) |
| 6 | NO pair, q = 0.04, terminator 0.80 | `coabundant_despite_terminator` (row 6) |
| 7 | OP, q = 0.02, Δ = +1.4, promoter 0.81 **and** terminator 0.80 | `other` (rows 1 and 4 both fail) |
| 8 | OP, q = 0.02, Δ = −1.4, promoter 0.81, terminator 0.10 | `other` (row 1 needs downstream higher) |
| 9 | NO, q = 0.04, promoter 0.81, terminator 0.10 | `other` (row 5 needs promoter no) |
| 10 | OP, both detected, q = 0.60, no terminator | no row (pattern `concordant`, no rule matches) |
| 11 | TU of 3 genes, 2 detected at 17.0 and 18.4, internal Δ = +1.4 and NA | `detected_frac` 0.667, `median_abundance_log2` 17.7, `max_abs_delta_log2` 1.4 |
| 12 | one-detected OP pair, downstream gene detected, q = 0.03, promoter 0.81, terminator 0.10 | `internal_promoter_supported` (pattern `downstream_higher`) |

## Acceptance criteria

- [ ] `tests/unit/operons/test_flags.py` asserts **every** row of the §7.4.5 table and at least
  three `other` combinations (examples 7–9) — the spec's unit-test item "every row of the flag
  table".
- [ ] Property: the flag is a pure function of (class, promoter_yes, terminator_yes, pattern); the
  full 2 × 2 × 5 × 2 truth table is enumerated in the test and compared against the table above.
- [ ] Property: every row of `operon_flags.tsv` exists in `pairs.tsv` with the same `q_value`, and
  no flagged pair has pattern `none`.
- [ ] Property: `Σ n_flags` over TUs ≤ the number of flagged OP pairs (NO pairs never count).
- [ ] `operons.tsv` has exactly one row per C5 TU; `genes` is in transcription order.
- [ ] Schema: the enum of `flag` is exactly `FLAGS` and matches §9.3.

## Hardening owed (pass 2)

- `evidence_note` wording reviewed with Khoa when M9's V11 and the report text land (it is
  user-facing prose).
- The §7.4.5 caveat rendered once in the report next to the table (M9), tracked here so it is not
  lost.
- A sensitivity rerun of the flag set under the UniOP thresholds 0.3 / 0.7 — M6-10.
- A golden `operon_flags.tsv` for the 50-kb fixture once the fixture has enough detected pairs.

## Questions for the lead

1. When the terminator layer is missing, rows 3 and 5 (`… no terminator …`) still fire, which reads
   as "no terminator here" when the truth is "we never looked". **Proposed default:** keep firing
   (the flag names are about the abundance pattern), but add `layer_missing:terminators` to the
   row's `flags` and say it in `evidence_note`. The alternative — suppress rows 3 and 5 entirely
   without a terminator layer — loses every `novel_operon_candidate` on runs without TerminatorNet.
2. `n_flags` counts pairs with both genes in the TU (proposed), or any flagged pair touching the
   TU? **Proposed default:** both genes, so the count is a property of the TU.
3. `evidence_note` free text in a controlled-vocabulary table: acceptable? **Proposed default:**
   yes, §9.3 defines it as a note column; the schema marks it as non-enumerated.

## Lead decisions (2026-10-21)
1. Q1: keep rows 3 and 5 without a terminator layer, adding `layer_missing:terminators` to the row's
   flags and saying it in `evidence_note` (M6-index decision 4).
2. Q2: `n_flags` counts flagged pairs with both genes in the TU.
3. Q3: `evidence_note` free text is accepted; the schema marks it non-enumerated.
4. Inputs are the merged M6-05 (`Discordance.rows`) and M6-06 (`MissedOperons.rows`) frames as they
   are on `m6/discordance` and `m6/missed-operons`. The internal `which` column of M6-06 is not a
   §9.2 column: drop it before anything is written (#148 review N5).
5. A tested discordance row carries `ci_unestimated`, not `unestimated` (#147 review F1). A flag
   whose row is `unestimated` never makes evidence "yes".
6. Imputed abundances (#148 review N2): a pair row whose gene has `abundance_imputed` gains the flag
   `abundance_imputed` here, where the gene table is joined, so a co-abundance or discordance call
   driven by imputation is visible. The statistician's answer on #148 may tighten this (exclusion
   from the targets); if so, that lands in M6-06 and this rule stays as the label.
