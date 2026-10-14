# M1-02: Genome FASTA reader and writer

Milestone M1 · depends on: M1-01 · estimate ≈ 320 changed lines (≈ 140 code, ≈ 180 tests)

## Spec sections implemented

- §5.3 Genome FASTA: "The first token of each header is the contig ID." · "`circular=true` or
  `[topology=circular]` in a header marks a circular replicon. Contigs without such a tag are linear
  unless `annotation.bakta.complete: true` declares every contig a complete, circular replicon." ·
  "Unsafe ID characters are mapped to `[A-Za-z0-9._-]`; the mapping is saved." · "gzip input is
  accepted."
- §7.0 rows "parses as FASTA: ≥ 1 record, IUPAC nucleotides only" and "contig IDs: unique first
  tokens; unsafe characters mapped …, map saved" (the rules themselves are checked in M1-09; this
  issue supplies the parsed facts they need).
- §6.2 `preflight` outputs `genome.norm.fa` and `contig_map.tsv` (written by M1-12 using this writer).

## Public interface

`src/giap/io/_open.py`
```python
def open_text(path: Path) -> TextIO: ...      # gzip detected by magic bytes 1f 8b, not by suffix; UTF-8
```

`src/giap/io/fasta.py`
```python
IUPAC_DNA: frozenset[str]   # ACGT N RYSWKMBDHV (upper case; input is upper-cased first)
SAFE_ID = re.compile(r"[A-Za-z0-9._-]")

class FastaFormatError(ValueError): ...      # carries .line_no (1-based)

@dataclass(frozen=True, slots=True)
class GenomeRecord:
    id: str               # safe ID (after mapping)
    original_id: str      # first header token as written
    header: str           # full header line without '>'
    sequence: str         # upper case, no whitespace
    circular_tag: bool    # header carries circular=true or [topology=circular]

def sanitize_id(raw: str) -> str: ...        # every char not in SAFE_ID → "_"
def read_genome_fasta(path: Path) -> list[GenomeRecord]: ...   # file order kept
def invalid_letters(seq: str) -> dict[str, int]: ...           # letters outside IUPAC_DNA → count
def write_fasta(path: Path, records: Iterable[tuple[str, str]], width: int = 80) -> None: ...
```

Reader rules: lines before the first `>` that are non-blank → `FastaFormatError`; blank lines
skipped; `\r\n` accepted; sequence whitespace removed; lower case upper-cased (soft-masking is not
kept). A header with an empty first token → `FastaFormatError`. The topology tag match is
case-insensitive on the key and value (`circular=true`, `Circular=True`, `[topology=circular]`).
The reader does not reject duplicates or invalid letters; it reports facts and M1-09 decides.
`write_fasta` writes `>id` only (no description), 80 columns, `\n` endings, uncompressed.

## Examples

| # | Input | Expected output |
|---|---|---|
| 1 | `>chr1 Escherichia coli circular=true\nACGT\nacgn\n` | `GenomeRecord(id="chr1", original_id="chr1", sequence="ACGTACGN", circular_tag=True)` |
| 2 | `>p1 [topology=circular] [gcode=11]\nGG\n` | `circular_tag=True`, `id="p1"` |
| 3 | `>ctg_7 len=5000\nAT\n` | `circular_tag=False` |
| 4 | `>NZ_CP0\|chr:1 x\nA\n` | `id="NZ_CP0_chr_1"`, `original_id="NZ_CP0\|chr:1"` |
| 5 | same content as #1 gzipped, file named `g.fa` (no `.gz`) | identical records to #1 |
| 6 | `ACGT\n>c\nA\n` | `FastaFormatError` with `line_no == 1` |
| 7 | `>c\nACGTXU-\n` | record parsed; `invalid_letters("ACGTXU-") == {"X": 1, "U": 1, "-": 1}` |
| 8 | empty file | `[]` (the "≥ 1 record" error is raised by M1-09, not here) |
| 9 | `write_fasta(p, [("c1", "A"*85)])` | `>c1\n` + 80 `A` + `\n` + 5 `A` + `\n` |

## Edge cases

- Record with a header but no sequence: parsed with `sequence=""` (M1-09 reports it as an error).
- Two headers mapping to the same safe ID (`a|b`, `a:b`): both returned; the duplicate is detected on
  `id` by M1-09.
- `circular=false` and `[topology=linear]` mean linear. A tag in the sequence part of the file is
  irrelevant (only headers are scanned).
- Trailing `*` or digits in sequence lines are invalid letters, not stripped.
- A gzip file that is truncated → `FastaFormatError` (wraps `EOFError`/`gzip.BadGzipFile`).

## Acceptance criteria

- [ ] `tests/unit/test_fasta_genome.py`: examples 1–9, one test each; example 5 builds the gzip at
  test time (no committed binary).
- [ ] Reading `write_fasta` output gives back the same `(id, sequence)` pairs (round trip test).
- [ ] The sim genome (`giapsim`, header `[topology=circular]`) reads as one circular record.
- [ ] Only `io/` touches files; `mypy --strict` clean.

## Questions for the lead

1. Is `U` acceptable in a genome? It is an IUPAC code but marks RNA. **Proposed default:** not in
   `IUPAC_DNA`, so it fails the "IUPAC nucleotides only" rule; `-` and `.` (gap characters) also fail.
2. Unsafe-ID mapping replaces each unsafe character with `_`. **Proposed default:** `_`; collisions
   after mapping are duplicate-ID errors (M1-09), never silently suffixed.
3. `annotation.bakta.complete: true` makes every contig circular (§5.3). **Proposed default:** the
   reader only reports `circular_tag`; topology is resolved in M1-09 (`resolve_topology`), so the
   reader stays config-free.
