# M1-07: Modification stripping and `gippa convert --from generic|diann`

Milestone M1 · depends on: M1-04, M1-06 · estimate ≈ 390 changed lines (fixtures excluded)

## Spec sections implemented

- §5.4: "Each converter is a pure function from input files to the generic table, with a contract test
  against a recorded sample file (§11.4). Until a real sample exists, the converter is written from the
  tool's documentation and its contract test is marked `pending_sample`".
- §5.4 DIA-NN 2.x (`report.parquet` or `.tsv`): "`peptide` ← `Stripped.Sequence`; `modified_peptide` ←
  `Modified.Sequence`; `protein_ids` ← `Protein.Ids`; `q_value` ← `Global.Q.Value`, else `Q.Value`;
  `pep` ← `PEP`; `is_decoy` ← `Decoy` (present with `--report-decoys`); `intensity` ←
  `Precursor.Normalised`, else `Precursor.Quantity`; `run` ← `Run`; `charge` ← `Precursor.Charge`;
  `source` = database; `--abundance-out` writes `PG.MaxLFQ` per protein group, with type `maxlfq`".
- §5.4 generic: "validates, normalises and reports problems".
- §5.4 stripping: "removes bracketed or parenthesised modifications (`(UniMod:35)`, `[+15.995]`,
  `[UNIMOD:35]`, `+15.995`, `n[42]`) and lower-case residues, and keeps only uppercase letters."
- §5.1: `gippa convert --from … INPUT [INPUT ...] --out peptides.tsv [--abundance-out abundance.tsv]
  [--source database|denovo] [--search-db same|other]`. Casanovo/InstaNovo/Cascadia are not M1.

## Public interface

`src/gippa/io/converters/strip.py`
```python
def strip_modifications(seq: str) -> str: ...
```
Order: remove `(…)` and `[…]` groups (non-nested), then `[+-]\d+(\.\d+)?` mass deltas, then every
character that is not `A–Z` (this drops lower-case residues and `n`/`c` terminus markers, `.`, `_`).

`src/gippa/io/converters/base.py`
```python
@dataclass(frozen=True, slots=True)
class ConvertResult:
    peptides: pd.DataFrame                 # generic columns only, §5.3 order, only columns with data
    abundance: pd.DataFrame | None         # columns protein_id, maxlfq
    issues: list[Issue]
```
`src/gippa/io/converters/generic.py`: `def convert_generic(paths: Sequence[Path], *, source: str |
None, search_db: str | None, cfg: InputConfig) -> ConvertResult`.
`src/gippa/io/converters/diann.py`: `def convert_diann(paths: Sequence[Path], *, search_db: str |
None, cfg: InputConfig) -> ConvertResult` (reads `.parquet` with pyarrow, else TSV).
`src/gippa/io/converters/__init__.py`: `CONVERTERS: Mapping[str, Callable[..., ConvertResult]]`.

Rules: several INPUT files are concatenated in argument order. Output rows sorted by `run`,
`peptide`, `modified_peptide`, `charge`. DIA-NN `Decoy` 0/1 → `false/true`; intensity `0` → `NA`
(DIA-NN writes 0 for unquantified precursors). `--abundance-out` (DIA-NN): one row per
`Protein.Group`, value = median over runs of `PG.MaxLFQ` (NA and 0 excluded), column name `maxlfq`.
`--source denovo` with `--from diann` → error. Generic: strips `peptide` when it contains
modification syntax (issue `stripped`, info, with count), canonicalises booleans to `true/false`,
drops unknown columns (issue `unknown_column`, warn), and fails on any error-severity issue.
CLI: writes with `write_tsv`; prints issues to stderr; exit 0, or 1 on error issues (no file written).

## Examples

| # | Input | Expected output |
|---|---|---|
| 1 | `strip_modifications("_PEPM(UniMod:35)IDEK_")` | `PEPMIDEK` |
| 2 | `strip_modifications("n[42]AC[+57.021]DEK")` | `ACDEK` |
| 3 | `strip_modifications("PEPS[UNIMOD:21]TIDE+15.995K")` | `PEPSTIDEK` |
| 4 | `strip_modifications("PEPmIDEK")` | `PEPIDEK` (lower-case residue removed) |
| 5 | DIA-NN TSV row `Run=r1, Stripped.Sequence=PEPTIDEK, Modified.Sequence=PEPTIDEK, Protein.Ids=P1;P2, Q.Value=0.002, Global.Q.Value=0.004, PEP=0.01, Precursor.Charge=2, Precursor.Normalised=1500, Precursor.Quantity=1400` | generic row `PEPTIDEK, PEPTIDEK, P1;P2, q_value=0.004, pep=0.01, intensity=1500, run=r1, charge=2, source=database` |
| 6 | same without `Global.Q.Value` and `Precursor.Normalised` | `q_value=0.002`, `intensity=1400` |
| 7 | DIA-NN rows `PG.MaxLFQ` 10 (r1), 30 (r2), 0 (r3) for `Protein.Group=P1;P2` | abundance row `P1;P2  20` |
| 8 | DIA-NN report lacking `Stripped.Sequence` | error issue `missing_column:Stripped.Sequence`; CLI exit 1, no output file |
| 9 | generic input with `peptide=PEPT[+79.966]IDEK`, `is_decoy=1` | `peptide=PEPTIDEK`, `is_decoy=true`, issue `stripped` (1 row) |

## Edge cases

- Empty DIA-NN report (header only) → `empty_table` error.
- `Decoy` column absent → no `is_decoy` column in the output.
- Same file passed twice → rows duplicated (aggregation handles it); no dedup in the converter.
- Output path ending `.gz` → gzipped (M1-04 rules).

## Acceptance criteria

- [ ] `tests/unit/test_strip.py` (examples 1–4 + Hypothesis: output matches `^[A-Z]*$`; stripping a
  plain uppercase sequence is the identity).
- [ ] `tests/unit/test_convert_diann.py` uses a hand-written TSV fixture built from the DIA-NN 2.x
  documentation (examples 5–8) and a parquet written at test time from the same frame.
- [ ] `tests/contract/test_diann_contract.py` is marked `pending_sample` and skips while
  `tests/data/recorded/diann/` is absent; the marker is registered in `pyproject.toml`.
- [ ] `tests/unit/test_cli_convert.py`: example 9 end to end through `CliRunner`; exit codes 0 and 1.

## Questions for the lead

1. `abundance.tsv` cannot carry the type. **Proposed default:** the converter prints
   "use `--abundance-type maxlfq`" and names the value column `maxlfq`; M1-08 does not infer the type
   from the column name (no hidden behaviour).
2. Multi-run `PG.MaxLFQ` collapse to one value (v1 allows one value column): median across runs.
   Accept?
3. DIA-NN `Precursor.Normalised == 0` → NA. Accept?

## Lead decisions after code review
- Stripping removes `(...)`, `[...]`, `{...}` and `<...>` groups repeatedly until none remain (nested
  modification names such as MaxQuant's `(Acetyl (Protein N-term))` leave no letters), then
  flanking residues: a string matching `^[A-Z-]\.(.+)\.[A-Z-]$` keeps only its middle.
- The DIA-NN reader loads only the mapped columns and converts numbers directly (500k rows < 15 s).
- MaxLFQ abundance excludes decoy rows and drops protein groups with no value.
- Converters report `validate_peptides` issues but keep warn-level rows in the output: any error
  already makes `gippa convert` exit 1 without writing, and preflight re-validates the table.
