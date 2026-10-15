# M2-11: `--annotation` GFF3/GBFF: Bakta `--regions` input and `annotation/id_map.tsv`

Milestone M2 · depends on: M2-03 · estimate ≈ 240 code lines (tests not counted) · capabilities: none
(one integration test `needs_conda` + `needs_bakta_db` runs Bakta `--regions` on sim100k)

## Spec sections implemented

- §5.3: "*GFF3 or GBFF from another annotator:* Bakta runs with `--regions` on that file's CDS, so
  the gene coordinates match. `annotation/id_map.tsv` maps the source IDs to Bakta locus tags by
  identical CDS coordinates." "*Bakta JSON:* used as is, and Bakta is not run."
- §7.1.1 Bakta gotcha: "`--regions` accepts CDS only, with at most 30 bp overlap with de novo CDS."
- §7.1.8: "`annotation/id_map.tsv`, when `--annotation` is a GFF3 or GBFF." §6.2 `bakta` inputs
  "[, `regions.gff3`]".
- §5.3 / M1-02: contig IDs are renamed to safe IDs; "the mapping is saved" (`contig_map.tsv`).

## Public interface

`src/giap/annotate/user_annotation.py`
```python
AnnotationFormat = Literal["bakta_json", "gff3", "gbff"]

@dataclass(frozen=True, slots=True)
class UserCds:
    source_id: str; interval: Interval; pseudo: bool

def detect_format(path: Path) -> AnnotationFormat     # by content: '{' → JSON, '##gff-version' → GFF3, 'LOCUS' → GBFF (gzip aware)
def read_user_cds(path: Path, fmt: AnnotationFormat, replicons: Mapping[str, Replicon],
                  contig_map: Mapping[str, str]) -> list[UserCds]
def write_bakta_regions(path: Path, cds: Sequence[UserCds], replicons: Sequence[Replicon]) -> None
@dataclass(frozen=True, slots=True)
class IdMapRow:
    source_id: str; locus_tag: str | None; contig: str; start: int; end: int; strand: str  # 0-based internally
    match: Literal["identical", "unmatched"]
def build_id_map(cds: Sequence[UserCds], genes: Sequence[Gene]) -> list[IdMapRow]
def write_id_map(path: Path, rows: Sequence[IdMapRow]) -> None
```
Rules:
- GFF3: only `type == CDS` lines (via M1-05 `read_gff3`); multi-line CDS sharing one `ID` (split
  across exons is impossible in bacteria, but origin-split rows happen) are merged into one interval.
  `source_id` = `locus_tag` attribute, else `ID`, else `Name`, else `<contig>_<start1>_<strand>`.
- GBFF: Biopython `SeqIO.parse(..., "genbank")` (BSD-3, an allowed in-process dependency); CDS
  features; `source_id` = `/locus_tag`, else `/protein_id`, else `/gene`; `join(…,…)` across the
  origin on a circular record → one interval with `end > length`.
- Contig names are translated through `contig_map` (original → safe); an unknown contig →
  `ValueError` (preflight turns it into an error, exit 1).
- `/pseudo` or `pseudo=true` CDS are written to the regions file (Bakta decides) but flagged.
- `write_bakta_regions`: GFF3 with `##sequence-region` lines and CDS lines only, `ID=<source_id>`.
- `build_id_map`: identical `(contig, start, end, strand)` → `identical`; otherwise `unmatched` with
  `locus_tag = NA`. One row per user CDS, sorted per §6.5.

`annotation/id_map.tsv` columns (Q1): `source_id`, `locus_tag`, `contig`, `start`, `end`, `strand`,
`match` (1-based closed).

## Examples (chr1 circular, 100,000 bp)

| # | Input | Expected |
|---|---|---|
| 1 | GFF3 with `gene`, `mRNA`, `CDS` (ID `cds-1`, `locus_tag=ABC_0001`, 101–400 +) | one `UserCds("ABC_0001", Interval("chr1",100,400,"+"), False)`; regions file has exactly one CDS line `chr1\t…\tCDS\t101\t400\t.\t+\t0\tID=ABC_0001` |
| 2 | GBFF CDS `complement(1001..1300)` `/locus_tag="ABC_0002"` | `Interval("chr1",1000,1300,"-")` |
| 3 | GBFF CDS `join(99901..100000,1..150)` | `Interval("chr1",99900,100150,"+")` |
| 4 | Bakta gene `GSIM_00005` at `(chr1,100,400,+)`; user CDS same coords | `IdMapRow("ABC_0001","GSIM_00005","chr1",100,400,"+","identical")`; TSV line `ABC_0001\tGSIM_00005\tchr1\t101\t400\t+\tidentical` |
| 5 | user CDS `(chr1,97,400,+)` (start differs) | `locus_tag=None`, `match=unmatched`; TSV `NA` |
| 6 | GFF3 seqid `chr 1|x` mapped by `contig_map` to `chr_1_x` | interval contig `chr_1_x` |
| 7 | GFF3 seqid `plasmidZ` not in the genome | `ValueError("annotation contig plasmidZ not in genome")` |
| 8 | `detect_format` on a gzipped GBFF | `"gbff"` |

## Edge cases

- Empty GFF3 (no CDS) → regions file with header only; Bakta runs de novo; id_map has a header only.
- Two user CDS with identical coordinates (duplicates) → both map to the same locus tag.
- A user CDS on a linear contig with `end > length` → `ValueError`.
- Bakta JSON given as `--annotation` → this module is not used (`detect_format` only).

## Acceptance criteria

- [ ] `tests/unit/test_user_annotation.py`: examples 1–8.
- [ ] Round trip: simulator GFF3 → `read_user_cds` → `build_id_map` against the simulator genes →
  every row `identical`.
- [ ] `tests/integration/test_bakta_regions.py` (`needs_conda`, `needs_bakta_db`, `slow`): Bakta light
  with `--regions` on sim100k keeps ≥ 99 % of the given CDS coordinates (id_map `identical` share).

## Questions for the lead

1. `id_map.tsv` columns are not in §9.3. **Proposed default:** the 7 columns above; added to the
   schema in M2-20.
2. Should `unmatched` rows try a looser match (same stop, different start)? **Proposed default:** no
   in v1 (spec says "identical CDS coordinates"); the count of unmatched rows goes to the manifest.

## Lead decisions (test-author questions)
1. Accepted: Bakta's `pseudogene=<type>` (GFF3) and `/pseudogene` (GBFF) count as pseudo, alongside
   `pseudo=true` / `/pseudo`; fallback IDs for CDS without one start `<contig>_<start>_` and differ by strand;
   `detect_format` raises ValueError on unrecognised content; `id_map.tsv` has exactly the 7 columns proposed
   in Q1 and keeps `end > length` for origin-crossing CDS; an origin-crossing CDS is written once in the
   regions file.
2. `annotation_source` stays in `bakta_layers` (M2-04); M2-11 defines no separate API for it.
