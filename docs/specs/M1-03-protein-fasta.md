# M1-03: Protein FASTA reader

Milestone M1 · depends on: M1-02 (`open_text`, `FastaFormatError`) · estimate ≈ 260 changed lines

## Spec sections implemented

- §5.3 Proteins FASTA: "The ID is the first header token. UniProt `sp|ACC|NAME` and `tr|ACC|NAME`
  resolve to `ACC`. A protein group written `A;B` is split, keeping the group ID." · "Records of ≤ 40 aa
  are routed to the peptide path, and records < 7 aa are dropped with a warning." · "Allowed letters:
  the 20 amino acids plus U. X, B, Z and J are counted. An internal `*` warns."
- §7.0 rows Proteins/IDs, alphabet, length classes, duplicate sequences (rules in M1-10).
- §5.2 `input.peptide_len: [7, 40]` ("shorter are dropped; longer protein records go to miniprot").

## Public interface (`src/gippa/io/protein_fasta.py`)

```python
STANDARD_AA: frozenset[str]          # ACDEFGHIKLMNPQRSTVWY + U
COUNTED_AA: frozenset[str]           # X B Z J
Route = Literal["peptide", "protein", "dropped"]

@dataclass(frozen=True, slots=True)
class ProteinRecord:
    id: str                      # group ID (resolved accessions joined by ";"), or the single ID
    members: tuple[str, ...]     # resolved accessions; len 1 for non-groups
    header: str
    sequence: str                # upper case; one trailing "*" removed
    internal_stops: int          # "*" remaining inside the sequence
    nonstandard: Mapping[str, int]   # counts of X/B/Z/J
    other: Mapping[str, int]     # any other character (e.g. O, digits, "-") → count

def resolve_protein_id(token: str) -> tuple[str, tuple[str, ...]]: ...   # (group_id, members)
def read_protein_fasta(path: Path) -> list[ProteinRecord]: ...
def route_by_length(length_aa: int, peptide_len: tuple[int, int] = (7, 40)) -> Route: ...
```

`resolve_protein_id`: split on `;`, then for each part, if it matches `^(sp|tr)\|([^|]+)\|` take
group 2; strip whitespace; empty parts are dropped. `group_id = ";".join(members)`.
`route_by_length`: `< min` → `dropped`; `min ≤ L ≤ max` → `peptide`; `> max` → `protein`.

## Examples

| # | Input | Expected output |
|---|---|---|
| 1 | header `>sp\|P0A7V8\|RS4_ECOLI 30S ribosomal protein S4` | `id="P0A7V8"`, `members=("P0A7V8",)` |
| 2 | header `>tr\|A0A0K\|A0A0K_ECOLI` | `id="A0A0K"` |
| 3 | header `>NISS_00012;NISS_00013 group` | `id="NISS_00012;NISS_00013"`, `members=("NISS_00012","NISS_00013")` |
| 4 | header `>sp\|P1\|X_ECOLI;sp\|P2\|Y_ECOLI` | `id="P1;P2"`, `members=("P1","P2")` |
| 5 | sequence `MKV*` | `sequence="MKV"`, `internal_stops=0` |
| 6 | sequence `MK*VLX` | `internal_stops=1`, `nonstandard={"X": 1}` |
| 7 | sequence `mkub` | `sequence="MKUB"`, `nonstandard={"B": 1}`, `other={}` (U is standard) |
| 8 | `route_by_length(6)`, `(7)`, `(40)`, `(41)` | `dropped`, `peptide`, `peptide`, `protein` |
| 9 | gzipped file with two records | same records as the plain file |

## Edge cases

- Header `>` with no token, or `>;` (no members after splitting) → `FastaFormatError`
  (reported by M1-10 as an ID error).
- Duplicate IDs are returned as-is (M1-10 reports them); duplicate sequences under different IDs are
  returned as-is (M1-10 lists them as Info).
- A sequence that is only `*` → empty sequence, length 0 → route `dropped`.
- Length used for routing is `len(sequence)` after removing the single trailing `*`.
- `sp|ACC` without a third field: does not match the UniProt pattern → ID stays `sp|ACC`.
- Empty file → `[]` (M1-10 decides whether proteins were required).

## Acceptance criteria

- [ ] `tests/unit/test_fasta_protein.py`: examples 1–9, one test each.
- [ ] Hypothesis property: for any sequence over `STANDARD_AA`, `sum(nonstandard) == 0`,
  `other == {}`, `internal_stops == 0`.
- [ ] Boundary tests for `route_by_length` at 6/7/40/41 using a non-default `peptide_len` too
  (e.g. `(8, 30)`).

## Questions for the lead

1. "A protein group written `A;B` is split, keeping the group ID": is the group ID the raw token or
   the resolved accessions? **Proposed default:** resolved accessions joined by `;` (example 4), so
   the group ID matches DIA-NN `Protein.Ids` written with bare accessions.
2. `O` (pyrrolysine) is not in "20 amino acids plus U" nor in the counted set. **Proposed default:**
   count it under `other` and warn like X/B/Z/J (no error).
3. Characters that are not letters (digits, `-`, `.`): **Proposed default:** counted under `other`,
   Warn in M1-10; the record is kept (the mapper ignores it until M3 decides otherwise).
