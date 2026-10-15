# M3-13: Route 2b, tblastn for protein leftovers (HSP merging)

Milestone M3 · depends on: M3-12 (`ProteinAlignment`, `assign_alignments`) · creates
`src/giap/mapping/tblastn.py`; reuses `workflow/envs/blast.yml` from M2-01 (creates it if M2-01 has not merged),
`workflow/modules/local/tblastn/main.nf` · **external tool: BLAST+ `tblastn`, `makeblastdb`**
(conda `blast` pinned, or Apptainer by digest) · estimate ≈ 260 code lines

## Spec sections implemented

- §7.2.3 "**tblastn leftovers.** Proteins with no assignment go to
  `tblastn -db_gencode 11 -seg no -evalue 1e-5`, with the same thresholds after merging collinear
  HSPs per locus." Config `mapping.proteins.tblastn_args`.
- §7.2.8 reasons `below_identity` and `no_hit` for proteins that fail both routes (written by
  M3-15).
- §6.2 step `tblastn`: "genome, leftover proteins → tabular hits".

## Public interface

```python
OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen nident positive"
def makeblastdb_cmd(genome: Path, db: Path) -> list[str]: ...          # -dbtype nucl -parse_seqids
def tblastn_cmd(db: Path, proteins: Path, threads: int, args: str) -> list[str]: ...  # -outfmt OUTFMT
@dataclass(frozen=True, slots=True)
class Hsp:
    qseqid: str; sseqid: str; qstart: int; qend: int; sstart: int; send: int   # as BLAST reports (1-based)
    length: int; nident: int; positive: int; bitscore: float; qlen: int
def parse_tblastn(path: Path) -> list[Hsp]: ...
def merge_hsps(hsps: Sequence[Hsp], replicons: Mapping[str, Replicon],
               max_gap_nt: int = 300, max_q_overlap_aa: int = 10) -> list[ProteinAlignment]: ...
```

- **Strand** from `sstart > send` (minus). **Collinear chain** per (protein, contig, strand):
  HSPs sorted by genome position in transcription order, each next HSP starting after the previous
  one on the protein (overlap ≤ `max_q_overlap_aa`) and within `max_gap_nt` on the genome.
  Greedy chaining from the best-bitscore HSP. On a circular contig, an HSP ending within
  `max_gap_nt` of the contig end may chain with one starting near 0 (origin-crossing locus,
  unwrapped interval).
- **Merged statistics:** identity = Σ nident / Σ length; coverage = |union of query ranges| /
  qlen; `rank` = order of chains by summed bitscore (1 = best). `frameshift`/`stop_codon` = 0.
- Merged alignments go through `assign_alignments` (same thresholds, `route = tblastn`).

## Examples

| # | Input (tabular rows) | Expected |
|---|---|---|
| 1 | `tblastn_cmd(db, p, 8, "-db_gencode 11 -seg no -evalue 1e-5")` | `["tblastn", "-db", "db", "-query", "p", "-num_threads", "8", "-db_gencode", "11", "-seg", "no", "-evalue", "1e-5", "-outfmt", OUTFMT]` |
| 2 | one HSP `P1 chr1 … qstart 1 qend 120 sstart 1001 send 1360 nident 110 length 120 qlen 120` | one alignment `Interval(chr1,1000,1360,+)`, identity 0.917, coverage 1.0 |
| 3 | two HSPs of P1: q 1–60 at s 1001–1180 and q 61–120 at s 1184–1363 (a frameshift between) | merged: `Interval(chr1,1000,1363,+)`, coverage 1.0 |
| 4 | same two HSPs but the second at s 9001–9180 | two chains, coverage 0.5 each → not assigned |
| 5 | minus-strand HSP `sstart 2360 send 2001` | `Interval(chr1,2000,2360,-)` |
| 6 | circular chr1 (L = 50,000): q 1–50 at s 49851–50000, q 51–100 at s 1–150 | one chain `Interval(chr1,49850,50150,+)`, crosses origin |
| 7 | same on a linear contig | two chains |
| 8 | contract: recorded tblastn output for leftover proteins of the 100-kb sim genome (`tests/data/recorded/blast/<version>/sim100k/`) | parses; diverged sim proteins (80–95 % identity) are assigned `divergent` |

## Edge cases

- `qseqid` containing `|` (UniProt headers): proteins are written with GIAP's resolved IDs before
  the search (M1-03), so tabular IDs match.
- No hits file (empty): every leftover → unmapped.
- Two chains on different CDSs of one protein: both rows kept (multi-locus), ranks 1 and 2.

## Acceptance criteria

- [ ] `tests/unit/test_tblastn.py`: examples 1–7.
- [ ] `tests/contract/test_tblastn_contract.py`: example 8 with `PROVENANCE.md`.
- [ ] Property: merging is invariant to the input order of HSPs; merged coverage ≤ 1.
- [ ] Integration (`needs_conda`/`needs_apptainer`): live `makeblastdb` + `tblastn` on the 100-kb
  sim genome match the recorded fixture.

## Questions for the lead

1. HSP chaining limits are not in the spec. **Proposed default:** `max_gap_nt = 300`,
   `max_q_overlap_aa = 10`, as new keys `mapping.proteins.hsp_max_gap_nt` and
   `hsp_max_q_overlap_aa` (ADR; additive config).
2. BLAST+ version pin. **Proposed default:** the bioconda `blast` version current when the env is
   built, recorded in `workflow/envs/blast.yml` and an ADR (spec says "pinned" without a number).

## Answered by Khoa (2026-10-14)
Question 1: `max_gap_nt = 300` and `max_q_overlap_aa = 10` are fixed module constants (keyword defaults of
the chaining function), documented with their rationale in the methods note; no config keys, no ADR.

## Lead decisions (test-author questions)
1. Coverage of a tblastn chain = |union of its HSP query ranges| / qlen (spec), not the span. `ProteinAlignment`
   (M3-12) gains an optional explicit coverage that M3-13 sets; when absent it stays the span (miniprot's
   single alignments, unchanged). Test: a chain with a query gap reports the union, not the span.
2. Accepted: both chaining limits inclusive; genome gap = `next.start − prev.end`, unwrapped across the
   origin; genome overlap between HSPs allowed (real frameshifts give a few-nt overlap); rank per protein;
   `alt_start` False; unknown contig → ValueError; process `TBLASTN` with outputs `hits` (tblastn.tsv) and
   `versions`, its own container selector, `tblastn_args` only via `ext.args`; makeblastdb checked by flags.

## Lead decisions after code review (PR #53, changes requested)
1. `makeblastdb` runs WITHOUT `-parse_seqids` (interface amended): with it, BLAST reports accession-style IDs
   as `gb|U00096.3|` / `ref|NC_000913.3|`, which never match GIAP's contig IDs. Without it BLAST reports the
   first header token unchanged. Builder, main.nf and tests drop the flag. A recorded fixture with an
   accession-style contig ID (e.g. a K-12 excerpt named `NC_000913.3`, pinned BLAST image) guards it.
2. Chaining picks the BEST chain, not the nearest extension: per protein × contig × strand, a dynamic program
   over compatible HSPs (same strand, query order, gap ≤ 300 nt, query overlap ≤ 10 aa, small genome
   overlaps allowed, unwrapped across the origin) maximises summed bitscore; ties → larger union coverage,
   then lower start. Test: the reviewer's A/B/C case → chain A+B, coverage 1.0, gene assigned.
3. Tests kill the surviving mutants: seed/best-chain choice, a chain exceeding a small circular contig's
   length (both checks), the `b0 > a0` order check, ValueError for out-of-contig hits and bad query ranges.
4. The parser rejects rows that do not have exactly 15 columns.
