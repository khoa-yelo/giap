# M2-19: Regions: prophage (with function subsets), plasmid, BGC, CRISPR, AMR; `regions.tsv`/`.bed`

Milestone M2 · depends on: M2-04, M2-18 · estimate ≈ 350 code lines (tests not counted) · capabilities: none

## Spec sections implemented

- §7.1.7 table: "`prophage` | geNOMAD proviruses | gene membership by ≥ 50 % CDS overlap. Prophage
  genes get a function subset by keyword rules on the Bakta product and geNOMAD markers: `integrase`
  (integrase, recombinase, excisionase), `repressor` (repressor, cI, Cro, antirepressor),
  `structural` (capsid, portal, tail, baseplate, head, terminase, scaffold, tape measure, fiber),
  `lysis` (holin, lysin, lysozyme, endolysin, spanin, Rz), `cargo` (everything else)"; "`plasmid` |
  Bakta replicon table, geNOMAD plasmid score | the whole replicon"; "`bgc` | antiSMASH regions |
  core biosynthetic genes flagged"; "`crispr` | Bakta CRISPR arrays | with flanking `cas` genes within
  1 kb, if any"; "`amr` | genes Bakta annotates via AMRFinderPlus | a gene set rather than a contiguous region".
- §5.2 `region_min_overlap_frac: 0.5` "gene in region if ≥50 % of its CDS lies inside";
  `context.region_types: [prophage, plasmid, bgc, crispr, amr]`.
- §6.5 "regions: `<type>_<n>`, numbered by contig order then start"; §9.3 `annotation/regions.tsv`:
  "`region_id`, `region_type`, `contig`, `start`, `end`, `source_tool`, `score`, `genes` (list), `notes`".
- Module 1 test: "region membership at the 50 % boundary".

## Public interface

`src/giap/annotate/regions.py`; step `region_layer`.
```python
RegionType = Literal["prophage", "plasmid", "bgc", "crispr", "amr"]
@dataclass(frozen=True, slots=True)
class Region:
    region_id: str; region_type: RegionType; contig: str; interval: Interval | None   # None = gene set (amr)
    source_tool: str; score: float | None; genes: tuple[str, ...]; notes: str

@dataclass(frozen=True, slots=True)
class RegionGene:
    region_id: str; locus_tag: str; overlap_frac: float
    subset: str | None            # prophage only: integrase|repressor|structural|lysis|cargo
    core_biosynthetic: bool       # bgc only

def overlap_frac(gene: Interval, region: Interval, replicon: Replicon) -> float   # origin-aware
def members(genes: Sequence[Gene], region: Interval, replicon: Replicon, min_frac=0.5) -> list[tuple[str, float]]
def prophage_subset(product: str, marker_description: str | None) -> str
def build_regions(genes, features, replicons, genomad: GenomadResult | None,
                  antismash: Sequence[BgcRegion] | None, cfg) -> tuple[list[Region], list[RegionGene]]
def write_regions_tsv(path, regions) / regions_bed(regions, replicons) -> list[BedRecord]
```
Rules:
- Prophage: one region per geNOMAD provirus; `source_tool = geNOMAD`, `score = virus_score`,
  `notes = taxonomy`; whole-contig virus → interval = whole contig, note `whole_contig`.
- Subset keywords, case-insensitive, whole-word, checked in the order integrase → repressor →
  structural → lysis → cargo on `product + " " + marker_description`; `cI` and `Cro` matched
  case-sensitively as whole words (to avoid `micro…`, `circ…`).
- Plasmid: each replicon with `kind == plasmid` (Bakta) or geNOMAD `plasmid_score ≥ 0.7` (Q2) →
  whole replicon `[0, length)`; `source_tool` `Bakta`, `geNOMAD` or `Bakta+geNOMAD`
  (no `;` inside a str cell, ADR-0006); `score` = plasmid score or `NA`; genes = all genes on it.
- BGC: one region per antiSMASH area; `source_tool = antiSMASH`, `score = NA`, `notes =
  products=<p1>,<p2>`; `core_biosynthetic` from `core_genes`.
- CRISPR: each Bakta `crispr` feature; `cas` genes = gene name `^cas` (case-insensitive) or product
  containing `CRISPR-associated`, whose interval lies within 1,000 bp of the array on either side;
  region interval = union of array and those genes (Q3); `source_tool = Bakta`.
- AMR: one region per contig holding ≥ 1 gene with `expert` containing `amrfinder`; `interval =
  None` (TSV `start`/`end` = `NA`); `source_tool = AMRFinderPlus`; not written to BED (Q1).
- Membership: `overlap_frac ≥ min_frac` (≥, so exactly 0.5 is in). Region numbering per type by
  contig order, then start (AMR by contig order).

## Examples (chr1 circular, 100,000 bp)

| # | Input | Expected |
|---|---|---|
| 1 | gene `(1000, 1300)` (300 bp), region `(1150, 5000)` | overlap 150 / 300 = 0.5 → member |
| 2 | gene `(1000, 1300)`, region `(1151, 5000)` | 149 / 300 = 0.496667 → not a member |
| 3 | gene `(99900, 100200)` (crosses origin), region `(0, 5000)` | overlap 200 / 300 = 0.666667 → member |
| 4 | products `phage integrase` / `CI repressor` / `cI` / `tail fiber protein` / `holin` / `DNA methylase` / `microcin immunity protein` | `integrase` / `repressor` / `repressor` / `structural` / `lysis` / `cargo` / `cargo` |
| 5 | product `tail fiber assembly integrase-like protein` | `integrase` (order of rules) |
| 6 | proviruses on chr1 at 50 kb and on p1 at 1 kb; antiSMASH area on chr1 at 20 kb | `prophage_1` (chr1), `prophage_2` (p1), `bgc_1` |
| 7 | CRISPR array `(60000, 60800)`, `cas1` at `(61500, 62400)`, `cas2` at `(63000, 63300)` | genes `(cas1,)` (gap 700 ≤ 1,000; cas2 gap 2,200 excluded); interval `(60000, 62400)` |
| 8 | two AMR genes on chr1, one on p1 | `amr_1` (chr1, NA, NA, genes = 2), `amr_2` (p1) |
| 9 | `regions.tsv` row for example 6 `prophage_1` at `(50000, 62000)`, score 0.97 | `prophage_1\tprophage\tchr1\t50001\t62000\tgeNOMAD\t0.97\t<genes ;-joined>\t<taxonomy>` |
| 10 | origin-crossing prophage | TSV `end > length`; BED two records `prophage_1_part1`, `prophage_1_part2` |

## Edge cases

- `modules` without `mobile` → no prophage regions, plasmid regions from Bakta only, run flag
  `layer_missing:genomad`; without `bgc` → no BGC regions, `layer_missing:antismash`.
- `context.region_types` restricts which types are built.
- A gene may belong to several regions (prophage inside a plasmid); `region_ids` is a list (M2-20).
- Empty everything → `regions.tsv` with a header only.

## Acceptance criteria

- [ ] `tests/unit/test_regions.py`: examples 1–10.
- [ ] `tests/property/test_region_membership.py`: `overlap_frac` ∈ [0, 1]; equals brute-force base
  counting on random intervals of both topologies.
- [ ] Contract-level test on recorded Nissle Bakta + geNOMAD + antiSMASH (`pending_sample` until M2-10):
  ≥ 1 prophage, ≥ 1 BGC, 2 plasmid regions (pMUT1, pMUT2).

## Questions for the lead

1. AMR "gene set" in a coordinate table. **Proposed default:** one row per contig, `start`/`end` `NA`,
   not in `regions.bed`.
2. geNOMAD plasmid-score cut-off for calling a contig a plasmid region is not specified. **Proposed
   default:** 0.7 (geNOMAD's default classification threshold), as a new config key
   `annotation.regions.genomad_plasmid_min_score` (ADR, default change).
3. CRISPR region bounds with `cas` genes. **Proposed default:** the interval is extended to cover the
   listed `cas` genes, so membership (≥ 50 %) includes them.

## Lead decisions (after the test-author's pass)
1. The test-author's resolutions are accepted: every geNOMAD plasmid row counts (M2-index decision 6
   overrides Q2's 0.7); CRISPR cas genes within 1 kb of the array, region genes = ≥ 50 % members;
   prophage gene descriptions from geNOMAD only on identical intervals; the step reads the single
   top-level antiSMASH JSON and warns `layer_missing:genomad`/`layer_missing:antismash`; parquet columns
   include the `Region`/`RegionGene` fields.
2. Test files carry the `regions_` prefix (M4-17 owns `tests/unit/context/test_regions*.py`).

## Lead decisions after code review (PR #91, changes requested)
3. CRISPR cas-gene search is origin-aware on circular replicons: the gap is the circular distance
   between the array and the gene (also for origin-crossing arrays); region extents stay in one frame
   (unwrapped `end > length` per §6.5). Tests for both review examples.
4. Nissle has THREE plasmids (CP171243.1 4.5 kb, CP171244.1 pMUT1, CP171245.1 pMUT2; chromosome
   CP171242.1) — the circular-replicon recording (bakta_full_circ, M2-10 lead decision 1) is correct and
   gives 3 plasmid regions; the linear recording typed CP171243.1 as `contig` only for lack of topology
   tags. The contract test uses the circular recording (falls back to skip if missing).
5. Contract C7 (src/giap/genetable/contract.py) gains optional `genes` and `notes` columns (no
   `extra_column` warning on real regions tables).
6. Tests kill the surviving mutants: the step with an AMR gene (`expert`); AMR rows in FASTA contig
   order (contig names that sort differently); `cas` as a prefix (`re.match`).
