# M1-01: Coordinates, `Interval` and `Replicon`

Milestone M1 · depends on: none · estimate ≈ 280 changed lines (≈ 120 code, ≈ 160 tests)

## Spec sections implemented

- §2.1 decision 16: "1-based closed in TSV, GFF3 and VCF; 0-based half-open in BED and proBed. Code
  works 0-based internally and converts only in the file writers."
- §6.5: "**`Replicon`:** `id`, `length`, `topology` (`circular` | `linear`), `kind` (`chromosome` |
  `plasmid` | `unknown`; from Bakta's replicon table or size)."
- §6.5: "**`Interval`:** `contig`, `start`, `end`, `strand`, 0-based half-open. On a circular replicon, a
  feature crossing the origin is stored with `end > length` (an unwrapped coordinate) and
  `crosses_origin = True`; helpers split it into `[start, length)` and `[0, end − length)`. Linear
  contigs never wrap."
- §6.5: "TSV and GFF3 use 1-based closed coordinates (`start + 1`, `end`). Origin-crossing features
  keep `end > length` … VCF positions are 1-based within `[1, length]`."
- §6.5 sort order: "All tables sort by contig (as in the input FASTA), then start, end, strand and ID."
- §11.2 unit: "origin-crossing coordinates; BED ↔ GFF"; property: "coordinate round trips".

## Public interface (`src/giap/model/coords.py`, re-exported from `giap.model`)

```python
Strand = Literal["+", "-", "."]
Topology = Literal["circular", "linear"]
RepliconKind = Literal["chromosome", "plasmid", "unknown"]

@dataclass(frozen=True, slots=True)
class Replicon:
    id: str; length: int; topology: Topology; kind: RepliconKind = "unknown"
    # __post_init__: length >= 1 else ValueError

@dataclass(frozen=True, slots=True)
class Interval:
    contig: str; start: int; end: int; strand: Strand = "."
    # 0-based half-open; __post_init__: 0 <= start < end else ValueError
    @property
    def length_bp(self) -> int: ...                       # end - start
    def crosses_origin(self, replicon: Replicon) -> bool: ...   # end > replicon.length

def check_interval(iv: Interval, replicon: Replicon) -> None: ...      # raises ValueError
def split_at_origin(iv: Interval, replicon: Replicon) -> tuple[Interval, ...]: ...
def to_one_based_closed(iv: Interval) -> tuple[int, int]: ...           # (start + 1, end)
def from_one_based_closed(contig: str, start1: int, end1: int, strand: Strand,
                          replicon: Replicon) -> Interval: ...
def wrap_position(pos0: int, replicon: Replicon) -> int: ...            # 0-based, in [0, length)
def interval_sort_key(iv: Interval, contig_order: Mapping[str, int], id_: str = "") -> tuple: ...
```

Validity (`check_interval`): `iv.contig == replicon.id`; `start < replicon.length`; on linear
`end <= length`; on circular `end - start <= length` (so `end < 2·length`).
`split_at_origin` returns `(iv,)` when the interval does not cross, else
`(Interval(c, start, length, s), Interval(c, 0, end − length, s))`, always in that order (the part
at the contig end first, for both strands). `from_one_based_closed` also accepts a wrapped input
`start1 > end1` on a circular replicon and unwraps it (`end = end1 + length`); on a linear one this
is a `ValueError`. Strand ranks for sorting: `+` < `-` < `.`.

## Examples (L = linear length 1000, C = circular length 1000)

| # | Input | Expected output |
|---|---|---|
| 1 | `to_one_based_closed(Interval("c1", 99, 200, "+"))` | `(100, 200)` |
| 2 | `from_one_based_closed("c1", 100, 200, "+", L)` | `Interval("c1", 99, 200, "+")` |
| 3 | `Interval("c1", 10, 19, "-")` → `to_one_based_closed` | `(11, 19)`, strand stays `-`, `length_bp == 9` |
| 4 | `from_one_based_closed("c1", 951, 1050, "-", C)` | `Interval("c1", 950, 1050, "-")`, `crosses_origin(C) is True` |
| 5 | `from_one_based_closed("c1", 951, 50, "+", C)` (wrapped) | `Interval("c1", 950, 1050, "+")` |
| 6 | `split_at_origin(Interval("c1", 950, 1050, "-"), C)` | `(Interval("c1", 950, 1000, "-"), Interval("c1", 0, 50, "-"))` |
| 7 | `check_interval(Interval("c1", 950, 1050, "+"), L)` | `ValueError` ("linear contigs never wrap") |
| 8 | `wrap_position(1003, C)` / `wrap_position(1003, L)` | `3` / `ValueError` |
| 9 | `Interval("c1", 5, 5)` | `ValueError` (empty interval) |
| 10 | sort `[("c2",0,9,"+"), ("c1",5,9,"-"), ("c1",5,9,"+")]` with FASTA order `c2, c1` | `c2:0-9:+`, `c1:5-9:+`, `c1:5-9:-` |

## Edge cases

- An interval covering the whole circular replicon (`start=0, end=length`) does not cross the origin.
- `end == length` on circular is not origin-crossing; `end == length + 1` is.
- `end − start == length` with `start > 0` on circular is valid (whole replicon, unwrapped);
  `end − start > length` is a `ValueError`.
- Negative `start` is a `ValueError`. Contig-ID mismatch with the replicon is a `ValueError`.
- `split_at_origin` on a linear replicon always returns `(iv,)` after `check_interval`.

## Acceptance criteria

- [ ] `tests/unit/test_coords.py` covers examples 1–10, each a separate test.
- [ ] `tests/property/test_coords_roundtrip.py` (Hypothesis): for any valid interval on a random
  replicon (both topologies, both strands), `from_one_based_closed(*to_one_based_closed(iv), …) == iv`,
  and the parts of `split_at_origin` have total length `iv.length_bp` and each lies in `[0, length]`.
- [ ] No file I/O in `giap.model`; `mypy --strict` and `ruff` clean.

## Questions for the lead

1. §6.5 lists `crosses_origin = True` as a stored field, but the flag depends on the replicon length.
   **Proposed default:** compute it (`iv.crosses_origin(replicon)`) instead of storing it, so it can
   never disagree with `end > length`; `Gene`/`Feature` in later milestones may cache it.
2. Is the `.` strand needed internally? **Proposed default:** yes, for regions and BED-only layers.
3. `Gene` and `Feature` (§4 layout lists them in `model/`) are not needed by M1. **Proposed default:**
   defer them to M2 (`parse_bakta`).
