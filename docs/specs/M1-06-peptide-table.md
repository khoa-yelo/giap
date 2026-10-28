# M1-06: Generic peptide table: read, validate, aggregate

Milestone M1 · depends on: M1-04 (`read_tsv`) · creates `src/gippa/io/issues.py` · estimate ≈ 380 changed lines

## Spec sections implemented

- §5.3 generic peptide table (all 13 columns; quoted rules below), aggregation:
  "**Group by stripped `peptide`.** I/L are kept as written; equivalence is applied only in matching."
  · "`q_value` and `pep` take the minimum over rows; `score` takes the maximum; `is_decoy` is true only
  if all rows are decoys; `protein_ids` is the union." · "**Intensity:** summed over charges and
  modified forms within a run, then the median across runs (`abundance.run_aggregation`)." ·
  "**Per-run detail:** per-run presence is kept in `work/peptide_runs.parquet`".
- Column rules: `peptide` str **required**, "stripped sequence, uppercase A–Z only"; `protein_ids`
  "`;`-separated … Decoy and contaminant prefixes are recognised from config"; `q_value`, `pep` float
  [0, 1]; `is_decoy` bool "requires `protein_ids` for gene-level pairing"; `intensity` float ≥ 0;
  `charge` int; `source` enum `database|denovo` "default from `input.peptide_source`"; `search_db`
  enum `same|other` "default from `input.search_db`"; `spectrum_id` "unique within `run` (G2)".
- §7.0 Peptides/schema (Error) and content (Warn) are evaluated in M1-10 from the `Issue`s below.

## Public interface (`src/gippa/io/peptides.py`)

```python
PEPTIDE_COLUMNS: tuple[str, ...]   # the 13 columns in §5.3 order
Severity = Literal["error", "warn", "info"]

@dataclass(frozen=True, slots=True)
class Issue:                        # shared by io validators and preflight (src/gippa/io/issues.py)
    code: str; severity: Severity; message: str
    n_rows: int = 0; examples: tuple[str, ...] = ()     # ≤ 5 examples, 1-based data-row numbers

def read_peptide_table(path: Path) -> tuple[pd.DataFrame, list[Issue]]: ...
def validate_peptides(df_raw: pd.DataFrame, cfg: InputConfig) -> tuple[pd.DataFrame, list[Issue]]: ...
def aggregate_peptides(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]: ...  # (peptides, peptide_runs)
```
`read_peptide_table` reads every column as string (`NA` and empty → missing), keeps unknown columns
(Issue `unknown_column`, info) and returns typed columns. Typing rules: bools accept
`true/false/TRUE/FALSE/True/False/1/0`; ints must be integral; floats parse with `float()`.

`validate_peptides` (pure) returns rows that are evidence-eligible plus issues:
| code | severity | condition |
|---|---|---|
| `missing_peptide_column` | error | no `peptide` column |
| `bad_type:<col>` | error | a value cannot be parsed as the column type |
| `out_of_range:<col>` | error | `q_value`/`pep` ∉ [0,1], `intensity` < 0, `charge` < 1 |
| `bad_enum:<col>` | error | `source`/`search_db` outside the enum |
| `duplicate_spectrum_id` | error | `(run, spectrum_id)` repeated with different `peptide` |
| `non_letters` | warn | `peptide` not matching `^[A-Z]+$` (row dropped) |
| `short_peptide` | warn | `len(peptide) < input.peptide_len[0]` (row dropped) |
| `long_peptide` | warn | `len(peptide) > input.peptide_len[1]` (row dropped) |
| `decoy_without_protein_ids` | warn | `is_decoy` true and `protein_ids` missing |
| `empty_table` | error | zero data rows after the header |

Derived columns: `is_decoy` missing → true when every `protein_ids` entry starts with a
`decoy_prefixes` entry, else false; `is_contaminant` = every entry starts with a
`contaminant_prefixes` entry; `source`/`search_db` filled from config when missing (`auto` stays
`auto`, resolved later). Rows with q > `q_value_max` are kept (needed for FDR estimation).

`aggregate_peptides` returns one row per `peptide` with `q_value`=min, `pep`=min, `score`=max,
`is_decoy`=all, `is_contaminant`=all, `protein_ids`=sorted union, `modified_peptides`=sorted
union, `n_rows`, `n_runs`, `intensity`=median over runs of the per-run sums (runs with no
intensity excluded; all missing → NA). Missing `run` → one run named `""`. `peptide_runs` has
`peptide, run, intensity, n_rows` (sum per run). Both sorted by `peptide` (then `run`).

## Examples (tables as TSV, `\t` shown as `|`)

| # | Input | Expected output |
|---|---|---|
| 1 | `peptide` / `PEPTIDEK` / `LLLLIIIK` | 2 peptides, no issues; `L` and `I` kept as written |
| 2 | no `peptide` column (`sequence\|q_value`) | `missing_peptide_column` error |
| 3 | `peptide\|q_value` / `PEPTIDEK\|1.2` | `out_of_range:q_value` error, example row `1` |
| 4 | rows `PEPTIDEK\|0.004\|r1\|2\|100`, `PEPTIDEK\|0.001\|r1\|3\|50`, `PEPTIDEK\|0.02\|r2\|2\|90` (`peptide\|q_value\|run\|charge\|intensity`) | one row: `q_value=0.001`, per-run sums `r1=150`, `r2=90`, `intensity=120`, `n_runs=2` |
| 5 | `PEPTIDEK` with `is_decoy` `true` and `false` on two rows | `is_decoy=false` |
| 6 | `protein_ids` `P1;P2` and `P2;P3` for one peptide | `protein_ids=["P1","P2","P3"]` |
| 7 | `peptide\|protein_ids` / `AAAAKPEP\|rev_P1` (no `is_decoy` column) | `is_decoy=true` |
| 8 | `PEPK` with `peptide_len=[7,40]` | row dropped, `short_peptide` warn, `n_rows=1` |
| 9 | `PEPT(UniMod:21)IDEK` in `peptide` | `non_letters` warn, row dropped |
| 10 | header only | `empty_table` error |

## Edge cases

- Ties in `score`/`q_value` do not matter (min/max are values). Median of an even number of runs is
  the mean of the middle two (example 4).
- `intensity = 0` is a value (summed); NA is missing. A peptide with all intensities NA has NA.
- `pep` present but `q_value` absent (and vice versa) is fine. Gzipped TSV is accepted.
- `spectrum_id` without `run`: uniqueness is checked over the whole table.

## Acceptance criteria

- [ ] `tests/unit/test_peptide_table.py`: examples 1–10, one test each.
- [ ] Property test: aggregation is invariant to row order (shuffle rows → identical frames).
- [ ] `peptide_runs` sums equal the per-run intensity sums of the input (property test).

## Questions for the lead

1. "Median across runs": over runs where the peptide was quantified, or all runs (absent = 0)?
   **Proposed default:** runs with a non-missing intensity only; absence is detection-limit data kept in
   `peptide_runs`.
2. Should `peptide` values with modification syntax be stripped by the reader (as converters do)?
   **Proposed default:** no; the generic table promises stripped sequences, so such rows warn and are
   dropped; `gippa convert --from generic` (M1-07) is the place that strips.
3. Peptides longer than `peptide_len[1]` (40 aa): drop, or route to the protein path? **Proposed
   default:** warn and drop (no de novo or DIA-NN peptide is that long in practice).
4. Rows of one peptide with different `source`: **Proposed default:** `database` wins; issue
   `mixed_source` (warn) with the count.
