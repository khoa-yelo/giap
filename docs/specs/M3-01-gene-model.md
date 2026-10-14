# M3-01: `GeneModel`, the gene-model interface Module 2 consumes

Milestone M3 · depends on: M1-01 (`Interval`, `Replicon`), M1-02 (genome FASTA) · creates
`src/giap/model/genes.py` · pure Python · estimate ≈ 230 code lines

## Why this issue exists

Module 2 (mapping) joins to the annotation only through Bakta locus tags (§6.1: "The only join
between the mapping track and the annotation track is the Bakta locus tag"). The full Bakta parser
(`parse_bakta` → `genes.parquet`, `features.parquet`, §6.2, §7.1.2) belongs to **M2** (the M2 issue
that implements `parse_bakta`; M2 issue files do not exist yet, so the lead should name it when M2
is split). M3 must not wait for it. This issue fixes the **interface** that every M3 issue
consumes and gives it one minimal constructor that reads the simulator's Bakta-like JSON (§11.3:
"a Bakta-like JSON containing only the fields GIAP parses (validated by GIAP's own parser)").
M2's `parse_bakta` MUST return the same `GeneModel` (contract: `GeneModel` is the API; parquet is
its serialisation), so nothing in M3 changes when M2 lands.

## Spec sections implemented

- §7.1.2 "`genes.parquet`, one row per CDS, sORF or pseudogene, with: locus tag, gene, product,
  type, interval, strand, `partial`, `crosses_origin`, `hypothetical`, amino-acid sequence,
  nucleotide sequence, start codon and `db_xrefs`. Non-coding features go to `features.parquet`."
- §7.1.1 Bakta gotcha: "Parse the JSON defensively: the contig field may be `sequence` or
  `contig`, and pseudogenes may be flagged CDS."
- §6.5 `Interval` 0-based half-open, origin-crossing stored unwrapped (`end > length`); partial
  genes "cut by a linear contig end".
- §7.2.6 needs, per feature: type (rRNA/tRNA/ncRNA for `rna_overlap`; `pseudogene`; `sorf`; CDS),
  interval, strand, frame and the protein sequence (for `multi_identical`).

## Public interface (`src/giap/model/genes.py`)

```python
GeneType = Literal["cds", "sorf", "pseudogene"]
NcType = Literal["rrna", "trna", "tmrna", "ncrna", "ncrna_region", "crispr", "oric", "orit", "gap"]

@dataclass(frozen=True, slots=True)
class Gene:
    locus_tag: str
    feature_type: GeneType
    interval: Interval            # strand "+" or "-"; may cross the origin (end > length)
    aa: str                       # translated protein, no terminal "*"; "" for pseudogenes without one
    nt: str | None
    start_codon: str | None
    gene: str | None = None
    product: str = ""
    partial: bool = False
    hypothetical: bool = False
    db_xrefs: tuple[str, ...] = ()

@dataclass(frozen=True, slots=True)
class NcFeature:
    id: str; type: NcType; interval: Interval

@dataclass(frozen=True)
class GeneModel:
    replicons: Mapping[str, Replicon]          # input FASTA order
    genes: tuple[Gene, ...]                    # sorted by interval_sort_key, then locus_tag
    features: tuple[NcFeature, ...]
    def by_locus(self, locus_tag: str) -> Gene: ...
    def overlapping(self, iv: Interval) -> list[Gene | NcFeature]: ...   # both strands; origin-aware

def gene_model_from_bakta_json(path: Path, replicons: Mapping[str, Replicon]) -> GeneModel: ...
```

`overlapping` uses a sorted-start array per contig plus a max-end sweep (or `numpy.searchsorted`);
an origin-crossing query or feature is split with `split_at_origin` before testing.

**JSON subset read** (per feature): `type` (`cds`, `sorf`, `tRNA`, `tmRNA`, `rRNA`, `ncRNA`,
`ncRNA-region`, `crispr`, `oriC`, `oriT`, `gap`; case-insensitive), contig from `sequence` or
`contig`, `start`/`stop` 1-based closed (`start > stop` on a circular contig = origin-crossing,
via `from_one_based_closed`), `strand` (`+`/`-`), `locus` (genes) or `id` (others), `aa`, `nt`,
`start_type`, `gene`, `product`, `db_xrefs`. A `cds` with `pseudo: true` or a `pseudogene` key is
`pseudogene`. `truncated` present (Bakta) → `partial = True`. `hypothetical` = `hypothetical`
key, else product matches `^hypothetical protein`. Unknown types are ignored (counted in a
returned `Issue`, info). Replicons come from the caller (preflight's `replicons.tsv`); a feature
on an unknown contig is an error.

## Examples

| # | Input feature (JSON) | Expected `Gene`/`NcFeature` |
|---|---|---|
| 1 | `{"type":"cds","contig":"chr1","start":101,"stop":400,"strand":"+","locus":"GSIM_00005","aa":"M…","nt":"ATG…","start_type":"ATG"}` | `Gene("GSIM_00005","cds",Interval("chr1",100,400,"+"),…)` |
| 2 | same with `"sequence":"chr1"` instead of `"contig"` | identical result |
| 3 | `{"type":"cds",…,"pseudo":true,"aa":""}` | `feature_type="pseudogene"`, `aa=""` |
| 4 | `{"type":"sorf","start":2001,"stop":2075,"strand":"-",…}` | `feature_type="sorf"`, `Interval("chr1",2000,2075,"-")` |
| 5 | circular chr1 of length 50,000, cds `start=49901, stop=200` | `Interval("chr1",49900,50200,"+")`, `crosses_origin` true |
| 6 | linear contig, cds `start=49901, stop=200` | `ValueError` naming the feature (linear contigs never wrap) |
| 7 | `{"type":"rRNA","start":5001,"stop":6542,"strand":"+","id":"r1"}` | `NcFeature("r1","rrna",Interval("chr1",5000,6542,"+"))` |
| 8 | the M0 simulator's `bakta_like.json` (seed 1, 50 kb) | `len(genes) == SimResult.n_genes`; every `aa` equals `truth/genes.tsv` length |
| 9 | `overlapping(Interval("chr1",49990,50010,"+"))` with the gene of example 5 | returns that gene |

## Edge cases

- Duplicate locus tags → `ValueError`. Genes are returned sorted (contig order, start, end, strand,
  locus tag) whatever the JSON order.
- Feature exactly at the contig end (`stop == length`) does not cross the origin.
- Minus-strand origin-crossing feature (Bakta writes `start > stop` for both strands) → unwrapped.
- `aa` ending in `*` is stripped of the terminal `*` (internal `*` kept for pseudogenes).
- Empty `features` list → a valid empty model.

## Acceptance criteria

- [ ] `tests/unit/test_gene_model.py`: examples 1–9, one test each.
- [ ] Property: for simulated genomes (giapsim, random seeds), `overlapping(g.interval)` contains
  `g` for every gene, and `nt` translated with table 11 equals `aa` (+ `*`).
- [ ] mypy strict; `GeneModel` is immutable (frozen).
- [ ] Docstring states that M2's `parse_bakta` must construct the same type (contract test lives
  with M2 against the recorded Bakta 1.12.2 JSON).

## Questions for the lead

1. Where does `GeneModel` live? **Proposed default:** `giap.model.genes` (spec §4: `model/` holds
   "Interval, Replicon, Gene, Feature"); M2's `parse_bakta` imports it and adds
   `to_parquet`/`from_parquet`.
2. Should the minimal JSON reader stay after M2? **Proposed default:** no; M2's `parse_bakta`
   replaces it, and this function becomes a thin wrapper around it (same tests keep passing).
3. Bakta `start_type` may be absent for partial genes. **Proposed default:** `None`.
