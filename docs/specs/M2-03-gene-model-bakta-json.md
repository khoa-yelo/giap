# M2-03: Gene model and defensive Bakta 1.12.2 JSON parser; Bakta command builder

Milestone M2 · depends on: M1-01, M1-02; shares `giap.model.genes` with M3-01 · estimate ≈ 330 code lines (tests not counted) ·
capabilities: none (contract test on recorded JSON is `pending_sample` until M2-10)

## Spec sections implemented

- §7.1.1 Bakta: "`bakta --db DB --output bakta --prefix P --threads N --keep-contig-headers
  (--complete | --replicons replicons.tsv) [--regions regions.gff3] genome.norm.fa`"; keeps "JSON
  (primary), GFF3, GBFF (for antiSMASH), FAA, TSV". Gotcha: "Parse the JSON defensively: the
  contig field may be `sequence` or `contig`, and pseudogenes may be flagged CDS. A contract test on
  recorded 1.12.2 JSON fixes the field names."
- §7.1.2: "`parse_bakta` produces `genes.parquet`, one row per CDS, sORF or pseudogene, with: locus
  tag, gene, product, type, interval, strand, `partial`, `crosses_origin`, `hypothetical`, amino-acid
  sequence, nucleotide sequence, start codon and `db_xrefs`. Non-coding features go to `features.parquet`."
- §5.2 `annotation.bakta: {db: auto, complete: auto, extra_args: []}` "complete=auto: --complete if
  all contigs circular"; §5.3 "*Bakta JSON:* used as is … The Bakta version and DB version are taken
  from the JSON."
- §6.5 `Replicon.kind` "from Bakta's replicon table or size"; origin-crossing features stored with
  `end > length`; IDs "genes: Bakta locus tags"; partial genes.

## Public interface

**Shared with M3-01.** M3-01 (`docs/specs/M3-01-gene-model.md`) defines `src/giap/model/genes.py`
(`Gene`, `NcFeature`, `GeneModel`, `gene_model_from_bakta_json`) as the interface Module 2 consumes,
and says "M2's `parse_bakta` MUST return the same `GeneModel`". This issue therefore **extends** that
module and never forks it: whichever of M2-03 / M3-01 merges second rebases onto the other. Fields
below marked *(M2)* are additions with defaults, so M3 code is unaffected.

`src/giap/model/genes.py` (additions to M3-01)
```python
@dataclass(frozen=True, slots=True)
class Gene:                          # M3-01 fields: locus_tag, feature_type, interval, aa, nt,
    ...                              # start_codon, gene, product, partial, hypothetical, db_xrefs
    expert: tuple[str, ...] = ()     # (M2) e.g. ("amrfinder",)
    # crosses_origin is derived: interval.crosses_origin(replicon) (M1-01), not stored

@dataclass(frozen=True, slots=True)
class NcFeature:                     # M3-01 fields: id, type, interval
    ...
    gene: str | None = None; product: str = ""                    # (M2)
    amino_acid: str | None = None; anticodon: str | None = None   # (M2) tRNA
    nt: str | None = None                                         # (M2) rRNA (16S → anti-SD)
```
M3-01's `NcType` gains `"oriv"`. In this file the names `Gene`, `NcFeature`, `GeneModel` refer to
`giap.model.genes`.

`src/giap/annotate/bakta.py`
```python
@dataclass(frozen=True, slots=True)
class BaktaAnnotation:
    model: GeneModel                              # replicons, genes, features (M3-01 type)
    sequences: Mapping[str, str]
    bakta_version: str; db_version: str; db_type: Literal["full", "light", "unknown"]
    skipped_types: Mapping[str, int]              # unknown feature types, counted
class BaktaParseError(ValueError): ...
def parse_bakta_json(path: Path, replicons: Mapping[str, Replicon] | None = None,
                     genome: Mapping[str, str] | None = None) -> BaktaAnnotation
# M3-01's gene_model_from_bakta_json(path, replicons) becomes `parse_bakta_json(path, replicons).model`
```
`src/giap/annotate/bakta_cmd.py`
```python
def bakta_db_path(db_setting: str, db_dir: Path) -> tuple[Path, Literal["full", "light"]]
def bakta_command(*, genome: Path, outdir: Path, prefix: str, db: Path, threads: int,
                  all_circular: bool, complete: bool | Literal["auto"], replicons_tsv: Path,
                  regions_gff3: Path | None, extra_args: Sequence[str]) -> list[str]
```

Parsing rules (each a test):
- Contig of a feature: `feature["sequence"]`, else `feature["contig"]`, else `BaktaParseError`.
- Sequence records: `data["sequences"]` (fall back to `data["contigs"]`); nucleotides from `nt`
  or `sequence`; `type` chromosome/plasmid → kind, anything else → `unknown`; `topology`
  circular/linear (missing → linear). When `replicons` is given (preflight's `replicons.tsv`), it
  wins and each JSON sequence must match its length (else `BaktaParseError`); a feature on an
  unknown contig is an error. When nucleotides are absent, `genome` must be given.
- Coordinates: Bakta 1-based closed `start`, `stop` → `Interval(contig, start-1, stop, strand)`.
  If `start > stop` (Bakta `edge` feature) and the replicon is circular → `end = stop + length`,
  `crosses_origin = True`; on a linear replicon → `BaktaParseError`.
- `type` `cds` → `cds`, or `pseudogene` when `pseudo` is true or a `pseudogene` key is present;
  `sorf` → `sorf`. RNA types are lower-cased (`tRNA`→`trna`, `ncRNA-region`→`ncrna_region`, `oriC`→`oric`).
- `partial`: `truncated` key present and non-empty, or the interval touches a linear contig end.
- `hypothetical`: `feature["hypothetical"]` if present, else product starts with `hypothetical protein`.
- `nt`: from the feature, else sliced from the sequence (strand-aware, origin-aware); `start_codon`
  = `start_type` if present, else `nt[:3]`. `aa` from the feature; missing on a pseudogene → `""`.
- tRNA `anti_codon` → upper case, `U`→`T`; `amino_acid` kept.
- AMR: `expert` entries whose `type`/`source` is `amrfinder` (recorded field name fixed by contract).
- Versions: `data["version"]["bakta"]`, `data["version"]["db"]["version"]`, `["db"]["type"]`.
- Locus tags unique among genes; output sorted by `interval_sort_key` (§6.5).

## Examples (chr1 circular, length 100,000)

| # | Feature (JSON) | Expected |
|---|---|---|
| 1 | `{"type":"cds","sequence":"chr1","start":101,"stop":400,"strand":"+","locus":"GSIM_00005",…}` | `Gene(… interval=Interval("chr1",100,400,"+"), feature_type="cds")`, not origin-crossing |
| 2 | same with `"contig":"chr1"` instead of `"sequence"` | identical `Gene` |
| 3 | `{"type":"cds",…,"start":1001,"stop":1300,"strand":"-"}`, no `start_type` | interval `(1000,1300,"-")`; `start_codon` = first 3 nt of `revcomp(seq[1000:1300])` |
| 4 | `{"type":"cds",…,"start":99901,"stop":150,"strand":"+","edge":true}` | interval `(99900, 100150)`; `interval.crosses_origin(chr1)` is true |
| 5 | example 4 on a linear `chr1` | `BaktaParseError` naming the locus |
| 6 | `{"type":"cds","pseudo":true,…}` | `feature_type="pseudogene"` |
| 7 | `{"type":"tRNA","anti_codon":"cag","amino_acid":"Leu",…}` | `NcFeature(type="trna", anticodon="CAG", amino_acid="Leu")` |
| 8 | `"version":{"bakta":"1.12.2","db":{"version":"6.0","type":"light"}}` | `bakta_version="1.12.2", db_version="6.0", db_type="light"` |
| 9 | two CDS with locus `X_1` | `BaktaParseError("duplicate locus tag X_1")` |
| 10 | `bakta_command(all_circular=True, complete="auto", regions_gff3=None, threads=16, …)` | `["bakta","--db",DB,"--output",OUT,"--prefix","P","--threads","16","--keep-contig-headers","--complete",GENOME]`; with `all_circular=False` → `--replicons replicons.tsv` instead of `--complete` |

## Edge cases

- Features of an unknown type (`oriT` future types, `cds` with no `locus`) are counted in
  `skipped_types`, not fatal; a CDS without `locus` is fatal (locus tags are the join key).
- `strand` `?` or missing on CRISPR/gap features → `Interval` strand `.`; on a gene → error.
- Empty `features` list → valid annotation with zero genes.
- The M0 simulator JSON (fields `contig`, `start`, `stop`, `locus`, `start_type`, `nt`, `aa`, no
  `version`) parses with `db_type="unknown"`.

## Acceptance criteria

- [ ] `tests/unit/test_bakta_json.py`: examples 1–9 on hand-written JSON snippets in the test.
- [ ] `tests/unit/test_bakta_cmd.py`: example 10, `--regions` appended when given, `extra_args` last
  before the genome, `db: auto` → `db` if `$GIAP_DB_DIR/bakta/db` exists else `db-light`.
- [ ] `tests/contract/test_bakta_contract.py` on `tests/data/recorded/bakta/1.12.2-light/sim100k/`
  and `1.12.2-full/{k12,nissle,sim100k}/`: parses; gene count equals the GFF3 CDS count; every
  locus tag of the FAA appears; `pending_sample` until M2-10 records them.
- [ ] `parse_bakta_json` on the M0 simulator output reproduces `truth/genes.tsv` coordinates.

## Questions for the lead

1. How does Bakta 1.12.2 encode origin-spanning features? **Proposed default:** `start > stop`
   with `edge: true` (as above); the contract test on a recorded circular genome confirms or the
   rule is changed in this PR's follow-up.
2. Should pseudogenes keep a translated `aa`? **Proposed default:** keep whatever Bakta gives; empty
   string when absent (mapping uses six-frame, not `aa`).
