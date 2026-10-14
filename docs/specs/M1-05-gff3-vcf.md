# M1-05: GFF3 writer and reader, VCF writer

Milestone M1 · depends on: M1-01, M1-04 · estimate ≈ 380 changed lines (≈ 170 code, ≈ 210 tests)

## Spec sections implemented

- §6.5: "TSV and GFF3 use 1-based closed coordinates (`start + 1`, `end`). Origin-crossing features
  keep `end > length`, as the GFF3 convention for circular sequences allows, and the
  `##sequence-region` line and region feature carry `Is_circular=true`." · "VCF positions are
  1-based within `[1, length]`."
- §9.5: "checks GFF3 and VCF with format validators (`pysam` for VCF; a GFF3 parser round trip)."
- §7.5.5 VCF 4.3 header: "`##fileformat=VCFv4.3`; `##source=giap-<version>`; `##reference=` with the
  genome SHA-256; `##contig` lines with lengths; `##INFO` and `##FORMAT` definitions." ·
  "FORMAT: `GT:AD:DP`" · "Sample column: `config.sample_name`". Appendix B.7 record example.

Decision: M1 includes a **minimal GFF3 reader** (needed for the §9.5 round trip, and reused in M2 for
`--annotation` GFF3); it reads what the writer writes plus standard third-party GFF3 feature lines.
Parsing other annotators' CDS for `--regions` and `id_map.tsv` is M2.

## Public interface

`src/giap/io/gff3.py`
```python
@dataclass(frozen=True, slots=True)
class GffFeature:
    interval: Interval; source: str; type: str
    score: float | None = None; phase: int | None = None
    attributes: Mapping[str, str | Sequence[str]] = field(default_factory=dict)  # order kept

def write_gff3(path: Path, features: Iterable[GffFeature],
               replicons: Sequence[Replicon], fasta: Mapping[str, str] | None = None) -> None: ...
def read_gff3(path: Path, replicons: Mapping[str, Replicon]) -> list[GffFeature]: ...
def escape_attr(value: str) -> str: ...        # percent-encode ; = & , TAB LF CR % and control chars
```
Writer: `##gff-version 3`; for each replicon in the given order a `##sequence-region <id> 1 <len>`
line, followed for circular replicons by a region feature line
`<id>\tgiap\tregion\t1\t<len>\t.\t.\t.\tID=<id>;Is_circular=true` (also before features).
Features sorted per §6.5; columns: start+1, end (unwrapped), score `.`/`%.6g`, strand, phase
`.`/0–2. Attributes in insertion order, list values joined by `,` after escaping. Optional
`##FASTA` section at the end. Reader: unescapes attributes; `start > end` is a parse error (GFF3
requires `start ≤ end`); `end > length` is accepted only on circular replicons. Unknown `seqid` → `ValueError`.

`src/giap/io/vcf.py`
```python
@dataclass(frozen=True, slots=True)
class VcfRecord:
    contig: str; pos0: int; ref: str; alts: tuple[str, ...]
    info: Mapping[str, str | int | float | Sequence[str]]    # order kept; flag lists joined "|"
    gt: str; ad: tuple[int, ...]; dp: int; id: str = "."; qual: str = "."; filter: str = "PASS"

@dataclass(frozen=True, slots=True)
class VcfHeaderField: id: str; number: str; type: str; description: str

def write_vcf(path: Path, records: Iterable[VcfRecord], replicons: Sequence[Replicon], *,
              giap_version: str, genome_sha256: str, sample: str,
              info_fields: Sequence[VcfHeaderField], format_fields: Sequence[VcfHeaderField]) -> None: ...
```
POS = `wrap_position(pos0, replicon) + 1`; records sorted by contig order then POS. INFO values are
escaped for `: ; = % ,` and whitespace (VCF 4.3 percent-encoding); list values in `FLAGS` use `|`.

## Examples (C = circular `chr` length 2000)

| # | Input | Expected output |
|---|---|---|
| 1 | feature `Interval("chr", 99, 300, "+")`, CDS, `ID=g1;Name=dnaA` | `chr\tgiap\tCDS\t100\t300\t.\t+\t0\tID=g1;Name=dnaA` |
| 2 | feature `Interval("chr", 1950, 2050, "-")` on C | `…\t1951\t2050\t.\t-\t…`; file has `##sequence-region chr 1 2000` and the `Is_circular=true` region line |
| 3 | attribute `product="tRNA-Leu; anticodon=CAG"` | `product=tRNA-Leu%3B anticodon%3DCAG`; reader returns the original string |
| 4 | round trip: write 3 features (one minus, one origin-crossing) then `read_gff3` | equal `GffFeature` list |
| 5 | reader on `chr\t.\tCDS\t1951\t2050…` with linear `chr` | `ValueError` |
| 6 | Appendix B.7 record: `pos0=1000, ref="GCT", alts=("TCT",)`, INFO as in B.7, `gt="1", ad=(0,3), dp=3`, sample `pilot_A3` | line equals the B.7 record (tab-separated, `FLAGS=mass_equivalent\|snv`) |
| 7 | `VcfRecord(pos0=2003, …)` on C | `POS` = `4` |
| 8 | written VCF opened with `pysam.VariantFile` | parses; `rec.pos == 1001`; header has `##contig=<ID=chr,length=2000>` and `##reference=sha256:<hex>` |

## Edge cases

- No features: header lines, sequence-region lines and region lines only. No records: valid VCF with
  header and `#CHROM` line only (pysam accepts it).
- `mixed` call: `gt="."`. Multi-allelic: `alts=("TCT","ACT")`, `ad` has 3 values.
- A feature whose phase is `None` on a CDS is written `.` (writer does not invent phases).
- GFF3 lines starting with `#` other than `##sequence-region`/`##gff-version`, and blank lines,
  are skipped by the reader; reading stops at `##FASTA`.

## Acceptance criteria

- [ ] `tests/unit/test_gff3.py` (examples 1–5) and `tests/unit/test_vcf.py` (examples 6–8).
- [ ] Example 8 is marked to run only when `pysam` is importable (dev extra; present in CI).
- [ ] Hypothesis round trip: random features on random replicons (both topologies, both strands,
  attribute values with `;=,&%`) survive `write_gff3` → `read_gff3` unchanged.
- [ ] Two writes of the same input give identical bytes.

## Questions for the lead

1. Exact `##reference=` form is not specified. **Proposed default:** `##reference=sha256:<64 hex>`.
2. A codon that crosses the origin on a circular replicon cannot be one VCF REF allele within
   `[1, length]`. **Proposed default:** the writer raises `ValueError` when `pos0 + len(ref)` crosses
   the contig end after wrapping; M7 decides how to report such codons (e.g. `variants.tsv` only).
3. Should the region line be written for linear replicons too (without `Is_circular`)? **Proposed
   default:** no; only `##sequence-region` for linear ones.
