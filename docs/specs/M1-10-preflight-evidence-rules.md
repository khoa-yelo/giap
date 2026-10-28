# M1-10: Preflight rules part B: inputs, proteins, peptides, abundance, evidence tier

Milestone M1 · depends on: M1-03, M1-06, M1-08, M1-09 · estimate ≈ 400 changed lines

## Spec sections implemented

- §1.1: `proteins.faa` is required "yes, unless peptides are given".
- §7.0 rows for Proteins (IDs, alphabet, length classes, duplicate sequences), Peptides (schema,
  content, evidence tier) and Abundance (shape, join, values, scale), quoted per rule below.
- §7.3.2 tiers: "`decoy_picked` — decoy rows with `protein_ids` exist and ≥ 80 % of them pair with a
  gene; `mayu` — otherwise, q-values or PEPs exist; `two_peptide` — otherwise, peptide sequences only;
  `protein_list` — proteins only, no peptides".
- §2.1 decision 3: protein-level abundance values "are used only when no peptide intensities are given".

## Public interface (`src/gippa/preflight/evidence.py`, pure)

```python
Tier = Literal["decoy_picked", "mayu", "two_peptide", "protein_list"]

def inputs_checks(has_peptides: bool, has_proteins: bool) -> list[Check]: ...
def protein_checks(records: Sequence[ProteinRecord], parse_error: str | None,
                   peptide_len: tuple[int, int]) -> list[Check]: ...
def peptide_checks(issues: Sequence[Issue], df: pd.DataFrame | None, cfg: InputConfig) -> list[Check]: ...
def candidate_tier(df: pd.DataFrame | None, has_proteins: bool) -> Tier: ...
def tier_check(tier: Tier) -> Check: ...
def abundance_checks(table: AbundanceTable | None, issues: Sequence[Issue],
                     proteins: Sequence[ProteinRecord] | None, scale: Scale,
                     requested_scale: str, peptides_have_intensity: bool,
                     cfg: GippaConfig) -> list[Check]: ...
```
`candidate_tier`: `decoy_picked` if any row has `is_decoy` true and non-empty `protein_ids`; else
`mayu` if a `q_value` or `pep` column has any non-missing value; else `two_peptide` if peptides
exist; else `protein_list`. The ≥ 80 % pairing test needs the gene model and is applied in M4; the
check value says `"provisional": true` for `decoy_picked`.

## Rules, with one passing and one failing example each

| Test | `rule_id` | Rule (spec) | Level | Passing example | Failing example |
|---|---|---|---|---|---|
| T-PRE-11 | `inputs.evidence_present` | proteins required unless peptides given (§1.1) | error | `--peptides` only | neither `--peptides` nor `--proteins` |
| T-PRE-12 | `proteins.ids` | "non-empty, unique; UniProt `sp\|ACC\|NAME` and `A;B` groups parsed" | error | `>sp\|P1\|X`, `>P2;P3` → IDs `P1`, `P2;P3` | `>P1`, `>sp\|P1\|Y` (duplicate `P1`); header `>` |
| T-PRE-13 | `proteins.alphabet` | "20 aa + U; X, B, Z, J counted; internal `*`" | warn | `MKUVL*` | `MKX*VL` → value `{"nonstandard":{"X":1},"internal_stops":1}` |
| T-PRE-14 | `proteins.length_classes` | "≤ 40 aa → peptide path, > 40 aa → proteins; < 7 aa dropped" | info / warn | 41-aa and 20-aa records → `info`, value `{"protein":1,"peptide":1,"dropped":0}` | a 6-aa record → `warn`, `dropped: 1` with its ID |
| T-PRE-15 | `proteins.duplicate_sequences` | "same sequence under several IDs" | info | distinct sequences → `ok` | `P1` and `P2` both `MKV…` → `info`, value `[["P1","P2"]]` |
| T-PRE-16 | `peptides.schema` | "required `peptide`; types and ranges of optional columns (§5.3)" | error | `peptide\|q_value` / `PEPTIDEK\|0.01` | no `peptide` column; `q_value=1.2`; `source=dna` |
| T-PRE-17 | `peptides.content` | "non-letters after stripping; length < 7 dropped; share of decoys; q-value range" | warn | 3 valid peptides, no decoys, q ≤ 0.01 | `PEPK` (short, dropped); `PEPT(UniMod:21)IDEK` (non-letters, dropped); value includes `decoy_share`, `q_min`, `q_max`, `n_dropped` |
| T-PRE-18 | `peptides.evidence_tier` | "which detection tier the columns allow (§7.3.2)" | info | always `info`; four cases: decoy rows with `protein_ids` → `decoy_picked`; `q_value` only → `mayu`; `peptide` only → `two_peptide`; proteins only → `protein_list` | (no failing case; the four cases are the tests) |
| T-PRE-19 | `abundance.shape` | "first column = ID, one numeric column; extra columns warn" | error / warn | `protein_id\|ibaq` | one column → error; `protein_id\|v\|w` → warn |
| T-PRE-20 | `abundance.join` | "IDs missing from proteins, and proteins missing from abundance" | warn | proteins `P1,P2`, abundance `P1,P2` | abundance `P1,P9` → value `{"abundance_not_in_proteins":["P9"],"proteins_not_in_abundance":["P2"]}` |
| T-PRE-21 | `abundance.values` | "non-numeric → error; negatives → error only on a linear scale …; zeros and NA counted; > 20 % missing → warn" | error / warn | `1e6, 2e5` linear | `abc` → error; `-1, 100` forced linear → error; 3 NA of 10 → warn |
| T-PRE-22 | `abundance.scale` | "log or linear auto-detected; `--abundance-scale` overrides" | info | `21.3, 22.0` auto → `{"scale":"log2","source":"auto"}` | `--abundance-scale linear` → `{"scale":"linear","source":"cli"}` (info, both) |

Join rule details: when `proteins` is `None`, the join check is `info` with value
`"deferred: no proteins.faa; IDs are matched to locus tags after annotation"`. Group IDs match
when any member matches. If peptides carry intensities, `abundance.scale` adds
`"used": false` (decision 3: the table is ignored), and the join check is still run.

## Edge cases

- `peptide_checks` with `df=None` (no peptides given) returns no checks; tier comes from proteins.
- Decoy share is computed over rows after validation; `decoy_share = 1.0` (all decoys) → warn.
- `q_value` column all NA → treated as absent for the tier.
- Abundance with zero rows → `abundance.shape` error.

## Acceptance criteria

- [ ] `tests/unit/preflight/test_evidence_rules.py` has `test_T_PRE_11_pass` / `_fail` … through
  `T_PRE_22`, with the four tier cases for T-PRE-18 and the three failing variants listed for
  T-PRE-16, T-PRE-17 and T-PRE-21.
- [ ] All inputs are in-memory literals (records, DataFrames); no fixture files.

## Questions for the lead

1. "share of decoys" has no threshold. **Proposed default:** always report it; warn when
   `decoy_share ≥ 0.5` (more decoys than targets suggests mislabelled rows).
2. "q-value range" has no threshold. **Proposed default:** warn when no row has
   `q_value ≤ input.q_value_max` (no evidence would pass), and report min and max.
3. Abundance IDs may be locus tags or gene IDs (§5.3), unknown before annotation. **Proposed default:**
   join warnings only when ≥ 1 abundance ID matches a protein ID (otherwise the table is assumed to be
   keyed by locus tags and the join is deferred as `info`).

## Lead decisions after code review
- Peptides given but every row dropped, and no proteins file → `inputs.present` error "no usable evidence"
  (the `protein_list` tier means a proteins file and no peptides, §7.3.2).
- `decoy_picked` needs decoy rows **and** a ranking column (q_value, pep or score); decoy rows without
  one fall back to the next tier, with a warning naming the reason (§7.3.2 picking ranks by q, PEP, score).
- "No row passes `input.q_value_max`" counts target rows only.
- No separate `peptides.decoy_share` rule; the share stays an Info value inside `peptides.content`.
- `peptides.schema` is only for errors; `unknown_column` and other warn-level findings go to
  `peptides.content`.
- Protein records whose groups share a member are listed as Info under `proteins.ids` (capped).
- `abundance_checks` resolves the scale once and reports that one value.
