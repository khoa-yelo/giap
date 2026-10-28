# M1-08: Abundance table reader and scale detection

Milestone M1 · depends on: M1-04, M1-06 (`Issue`) · estimate ≈ 220 changed lines

## Spec sections implemented

- §5.3 Abundance table: "The first column is a protein ID that matches `proteins.faa`, a locus tag, or
  a gene ID. Exactly one numeric value column is allowed in version 1; extra columns warn and are
  ignored. The scale is detected automatically (`log2` if any value is negative or all values are
  ≤ 50; otherwise `linear`) unless `--abundance-scale` is given. Negative values are an error only on a
  linear scale. `--abundance-type` declares iBAQ, Top3, MaxLFQ or raw intensity (default `unknown` =
  relative)."
- §7.0 Abundance rows: shape "first column = ID, one numeric column; extra columns warn" (Error/Warn);
  values "non-numeric → error; negatives → error only on a linear scale (checked after scale
  detection); zeros and NA counted; > 20 % missing → warn"; scale "log or linear auto-detected;
  `--abundance-scale` overrides" (Info). Config: `abundance.scale`, `preflight.abundance_missing_warn_frac`.

## Public interface (`src/gippa/io/abundance.py`)

```python
Scale = Literal["log2", "linear"]

@dataclass(frozen=True, slots=True)
class AbundanceTable:
    id_column: str; value_column: str
    ids: tuple[str, ...]; values: tuple[float | None, ...]     # None = NA
    ignored_columns: tuple[str, ...]

@dataclass(frozen=True, slots=True)
class AbundanceStats:
    n_rows: int; n_na: int; n_zero: int; n_negative: int; missing_frac: float

def read_abundance(path: Path) -> tuple[AbundanceTable | None, list[Issue]]: ...
def detect_scale(values: Sequence[float | None]) -> Scale: ...
def resolve_scale(values: Sequence[float | None], requested: Literal["auto","log2","linear"]) -> Scale: ...
def abundance_stats(table: AbundanceTable) -> AbundanceStats: ...
def validate_abundance(table: AbundanceTable, scale: Scale, missing_warn_frac: float) -> list[Issue]: ...
```
Value column = the second column. A second column that is not numeric (after excluding NA/empty)
→ `non_numeric_values` error with example rows. Fewer than two columns → `no_value_column` error.
Columns after the second → `extra_columns` warn (names listed). Duplicate IDs → `duplicate_ids`
error. Empty ID → `empty_id` error. `detect_scale` ignores NA; all-NA → `linear` with issue
`scale_undetermined` (warn). `missing_frac = n_na / n_rows`.

## Examples

| # | Input (`\|` = tab) | Expected output |
|---|---|---|
| 1 | `protein_id\|ibaq` / `P1\|1e6` / `P2\|2.5e5` | table ok; `detect_scale → linear`; no issues |
| 2 | `protein_id\|log2_abund` / `P1\|21.3` / `P2\|-0.5` | `log2` (a negative value) |
| 3 | values `12.1, 30.0, 50` | `log2` (all ≤ 50) |
| 4 | values `12.1, 50.01` | `linear` |
| 5 | `protein_id\|v\|w` / `P1\|3\|x` | issue `extra_columns` warn listing `w`; value column `v` |
| 6 | `protein_id\|v` / `P1\|abc` | `non_numeric_values` error, example row `1` |
| 7 | values `-1, 100` with `requested="linear"` | `negative_on_linear` error |
| 8 | values `-1, 100` with `requested="auto"` | scale `log2`, no error |
| 9 | 10 rows, 3 `NA`, 1 `0` | `n_na=3, n_zero=1, missing_frac=0.3`; `too_many_missing` warn (> 0.20) |
| 10 | 10 rows, 2 `NA` | `missing_frac=0.2`, no warn (strictly greater than) |

## Edge cases

- Header only → `empty_table` error. A single column → `no_value_column` error.
- `0` on a log2 scale is a valid value (log2 abundance 0); on linear it is counted in `n_zero`.
- Scientific notation parses as a float. Non-finite text (`inf`, `INF`, `-inf`, `nan`, `NaN`, and
  overflow such as `1e400`) is an error (`non_finite`); only `NA` or an empty cell is missing
  (ADR-0006; lead decision after review).
- IDs are resolved as in `proteins.faa` (ADR-0005); the IDs as written are kept in `raw_ids`.
  Protein groups keep their written order (`P2;P3` and `P3;P2` are different IDs).
- A table with valid shape is always returned; cells with numeric errors become `None` and are
  flagged in `bad_rows`, so preflight still reports counts and scale. Only shape errors return no
  table.

## Acceptance criteria

- [ ] `tests/unit/test_abundance.py`: examples 1–10, one test each.
- [ ] Boundary tests for `detect_scale` at exactly 50 and 50.000001.
- [ ] Gzipped input reads identically.

## Questions for the lead

1. Does a UniProt-style ID (`sp|P1|X`) in the abundance file resolve to `P1` like in `proteins.faa`?
   **Decision:** yes, using `resolve_protein_id` from M1-03, so joins succeed; the IDs as written are kept in `raw_ids`.
2. Linear zeros: count as missing for the 20 % rule? **Proposed default:** no; only NA/empty are
   missing (§7.0 says "zeros and NA counted" separately). Zeros are reported in `preflight.json`.
3. Header required? **Proposed default:** yes, the first line is always the header.
