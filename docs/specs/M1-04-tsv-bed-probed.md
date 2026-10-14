# M1-04: TSV reader/writer conventions, BED and proBed writers

Milestone M1 · depends on: M1-01 · estimate ≈ 360 changed lines (≈ 160 code, ≈ 200 tests)

## Spec sections implemented

- §9.2: "TSV, UTF-8, one header row, `NA` for missing, one entity per row, sorted (§6.5)" ·
  "snake_case, with units in the name" · "Booleans `true` / `false`" · "`flags` column: `;`-separated".
- §9.3: "Types: `str`, `int`, `float`, `bool`, `enum{…}`, `list` (`;`-separated)."
- §6.6: "Floats are written with 6 significant digits (`%.6g`); integers are written as integers;
  booleans as `true` and `false`." · "No timestamps or absolute paths" in tables.
- §6.5: "BED and proBed are 0-based half-open; origin-crossing features are written as two records
  with `name` suffixes `_part1` and `_part2`."
- §7.2.9: "`mappings.probed` (proBed 1.0, BED12+13)". R7: byte-identical outputs.

## Public interface

`src/giap/io/tsv.py`
```python
NA = "NA"
COLUMN_NAME = re.compile(r"^[a-z][a-z0-9_]*$")
Cell = str | int | float | bool | None | Sequence[str]

def format_value(v: Cell) -> str: ...
def write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Sequence[Cell]]) -> None: ...
def read_tsv(path: Path) -> tuple[list[str], list[list[str | None]]]: ...   # header, cells (NA/"" → None)
```
`format_value`: `None`, `float('nan')`, `pd.NA` → `NA`; `bool` → `true`/`false` (checked before
`int`); `int`/`numpy.integer` → `str(int)`; float → `%.6g`, with `-0` → `0`, `inf` → `inf`,
`-inf` → `-inf`; list/tuple → `;`-join of `format_value` of items (empty list → `NA`); `str` as is,
but a tab, `\n` or `\r` inside raises `ValueError`. `write_tsv` checks every column name against
`COLUMN_NAME` and row width; writes `\n` endings, no trailing tab. A path ending in `.gz` is written
with `gzip.GzipFile(mtime=0, filename="")` so bytes are reproducible. `read_tsv` uses `open_text`
(gzip by magic), strips a UTF-8 BOM, errors on ragged rows with the 1-based line number.

`src/giap/io/bed.py`
```python
@dataclass(frozen=True, slots=True)
class BedRecord:
    interval: Interval; name: str; score: int = 0
    thick: Interval | None = None; rgb: str = "0"
    blocks: tuple[tuple[int, int], ...] = ()    # (size, start relative to interval.start); () = BED6

def bed_lines(rec: BedRecord, replicon: Replicon) -> list[str]: ...  # 1 line, or 2 if origin-crossing
def write_bed(path: Path, records: Iterable[BedRecord], replicons: Mapping[str, Replicon]) -> None: ...
```
`src/giap/io/probed.py`
```python
PROBED_EXTRA = ("proteinAccession", "peptideSequence", "uniqueness", "genomeReferenceVersion",
    "psmScore", "fdr", "modifications", "charge", "expMassToCharge", "calcMassToCharge",
    "psmRank", "datasetID", "uri")
@dataclass(frozen=True, slots=True)
class ProBedRecord:
    bed: BedRecord; extra: Mapping[str, Cell]      # keys ⊆ PROBED_EXTRA; missing → "."
def write_probed(path: Path, records: Iterable[ProBedRecord],
                 replicons: Mapping[str, Replicon], genome_version: str) -> None: ...
```
BED rows are BED6 when `blocks == ()`, BED12 otherwise; strand `.` written as `.`. proBed is
always BED12 + 13 columns (25 tab-separated fields), one block if none given, with a
`# proBed-version 1.0` header comment. Records are sorted by contig order, start, end, strand,
name after splitting; for an origin-crossing record, `_part2` (at 0) therefore precedes `_part1`.

## Examples (C = circular `c1` length 1000)

| # | Input | Expected output |
|---|---|---|
| 1 | `format_value(0.000123456789)`, `(1234567.0)`, `(3)`, `(True)` | `0.000123457`, `1.23457e+06`, `3`, `true` |
| 2 | `format_value(None)`, `(float("nan"))`, `(-0.0)`, `(["a","b"])` | `NA`, `NA`, `0`, `a;b` |
| 3 | `write_tsv(p, ["gene_id","length_aa"], [["g1", 300]])` | file bytes `gene_id\tlength_aa\ng1\t300\n` |
| 4 | `write_tsv(p, ["GeneID"], …)` | `ValueError` (not snake_case) |
| 5 | BED6 of `Interval("c1", 99, 200, "-")`, name `g1` | `c1\t99\t200\tg1\t0\t-` |
| 6 | BED6 of `Interval("c1", 950, 1050, "+")`, name `g2`, on C | `c1\t0\t50\tg2_part2\t0\t+` then `c1\t950\t1000\tg2_part1\t0\t+` |
| 7 | same interval on linear `c1` length 1000 | `ValueError` from `check_interval` |
| 8 | proBed for `Interval("c1", 9, 30, "+")`, `peptideSequence="PEPTIDEK"`, `fdr=0.001` | 25 fields; fields 10–12 = `1`, `21`, `0`; field 14 = `PEPTIDEK`; field 18 = `0.001`; unspecified extras `.` |
| 9 | `write_tsv("x.tsv.gz", …)` twice | identical bytes both times |

## Edge cases

- Zero rows: header line only. Empty-string cell `""` is written as an empty field, not `NA`
  (callers use `None` for missing).
- `numpy.float32` values format like `float`. Large ints (> 2⁵³) stay exact.
- BED12 blocks on an origin-crossing record: blocks are split at the origin with the interval and
  re-based on each part; `blockCount` per part counts only the blocks in that part.
- `thick` defaults to the full interval (`thickStart = chromStart`, `thickEnd = chromEnd`).

## Acceptance criteria

- [ ] `tests/unit/test_tsv.py`, `test_bed.py`, `test_probed.py` cover examples 1–9.
- [ ] Property test: `read_tsv(write_tsv(rows))` returns `format_value` strings, `NA` → `None`.
- [ ] Property test: for random intervals on circular replicons, BED parts cover exactly
  `iv.length_bp` bases and all coordinates lie in `[0, length]`.
- [ ] Every proBed line has exactly 25 fields (checked on 50 random records).

## Questions for the lead

1. proBed `uniqueness`, `psmRank` etc. come from Module 2 (M3). **Proposed default:** M1 ships the
   generic writer; M3 decides values. Missing proBed values are `.` (proBed convention), unlike TSV `NA`.
2. Are `%.6g` exponents (`1.23457e+06`) acceptable in TSV? **Proposed default:** yes, it is the literal
   §6.6 rule; Frictionless parses them as numbers.
3. BED `score` must be 0–1000; **Proposed default:** writer raises `ValueError` outside that range.
