# M9-08: V1 `circos`: where in the genome is the proteome?

Milestone M9 · **pass 1** · depends on: M9-01, M9-03 (optional tracks) · estimate ≈ 380 code lines
(+ ≈ 320 test lines) · not statistical

## Spec sections implemented

- §10.2 V1: "contigs as sectors with ticks; genes by strand coloured by detection; abundance heat
  track; peptide density (10-kb bins); TU boundaries; promoter and terminator ticks; prophage,
  plasmid, BGC and CRISPR arcs; islands; GC skew; markers for ALT sites and novel ORFs"; pyCirclize;
  static SVG in the report.
- §10.1 (colours, 20,000 elements per track, double width), §6.5 (circular/linear, origin
  crossing, 1-based TSV coordinates), §9.3 `genes.tsv`, `peptide_hits.tsv`, `tus.tsv`,
  `regions.tsv`, `islands.tsv`, `variants.tsv`, `novel_orfs.tsv`; BED files (0-based).

## Tracks (outer → inner radius; each optional except sectors and genes)

| # | Track | Source | Drawing |
|---|---|---|---|
| 1 | sectors + ticks | `preflight.replicons` (id, length, topology) | one sector per replicon in FASTA order, angle ∝ length; ticks every round 10ⁿ giving 4–10 ticks on the largest sector, labels in kb or Mb; linear contigs show a gap mark at both ends |
| 2 | region arcs | `annotation/regions.tsv` | prophage, plasmid, BGC, CRISPR, AMR arcs in the §10.1 colours (all five §9.3 types; §10.2 names four, AMR is drawn too) |
| 3 | genes, + strand | `genes.tsv` `strand = +` | rectangles coloured by `detection_status` (index decision 9) |
| 4 | genes, − strand | `genes.tsv` `strand = -` | same |
| 5 | abundance heat | `genes.tsv` measured `abundance_log2` | one cell per gene, viridis over the 2nd–98th percentile of measured values; NA and imputed drawn white |
| 6 | peptide density | `peptide_hits.tsv` | accepted placement rows per 10-kb bin (by `start`), bars, per-sector max scale; bins tile from position 1 and do not wrap |
| 7 | TU boundaries + promoter / terminator ticks | `annotation/tus.tsv`, `promoters.bed`, `terminators.bed` | TU starts as thin black ticks; promoters up-ticks `#009E73`, terminators down-ticks `#D55E00` |
| 8 | islands | `islands.tsv` rows with `p_value ≤ stats.alpha.context` | arcs, high `#0072B2` `///`, low `#D55E00` `\\\` |
| 9 | GC skew | `plots/support/gc_skew.tsv` (M9-03) | line of cumulative GC skew, per-sector symmetric scale |
| 10 | markers | `variants.tsv` (`call ∈ {ALT, mixed}`), `novel_orfs.tsv` (`called = true`) | ALT sites `▲` black (filled `confidence = high`, open `candidate`); novel ORFs `●` `#CC79A7` |
| — | oriC / ter | `plots/support/replicon_landmarks.tsv` | radial line labelled "oriC" and dotted "ter" on circular sectors with `status = called` |

**Coordinates.** TSV inputs are 1-based closed (converted once to 0-based half-open on read); BED
is 0-based. On circular replicons an interval with `end > length` (origin crossing, §6.5) is drawn
as two pieces, [start, length) and [0, end − length). On linear contigs nothing wraps.

**Binning** (index decision 12). A track with more than 20,000 elements switches to 20,000
equal-width bins per sector (M9-01 `binned`): gene tracks draw the share of detected genes
per bin on viridis (the legend says "detected share"), abundance draws the bin median, ticks draw counts. Peptide density is already binned at 10 kb.

## Data TSV `V01_circos.tsv`

Long format: `track` enum{sector, region, gene_plus, gene_minus, abundance, peptide_density,
tu_start, promoter, terminator, island, gc_skew, alt_site, novel_orf, oric, ter}, `contig`,
`start`, `end` (1-based closed; origin-crossing rows keep `end > length`), `value`, `category`,
`binned` (bool). Sorted by track (enum order), contig (FASTA order), start, end, category.

## Title

"`<d>` of `<n>` genes detected (`<x> %`); `<k>` islands at FWER ≤ `<alpha>`" (+ "; `<a>` ALT
sites, `<m>` novel ORFs" when those tables exist). Counts only.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| one linear contig | one sector with a gap mark; no wrap |
| a contig with no genes | sector drawn, empty tracks |
| a gene crossing the origin | two rectangles; one data TSV row with `end > length` |
| no `peptide_hits` (protein-only) | track 6 omitted, noted |
| `variants.tsv` absent (variants step not on main yet) or `novel_orfs.tsv` absent | the marker kind omitted |
| no support files (`giap report` on an old `OUT/`) | tracks 9 and oriC/ter omitted, noted |
| no regions (`layer_missing:genomad`, `layer_missing:antismash`) | prophage/BGC arcs absent; CRISPR/plasmid arcs from Bakta still drawn; the figure notes the missing layers |
| a 20-kb replicon among 5-Mb ones | sector drawn at its true angle (≈ 0.4 % of the circle) and labelled outside |
| > 20,000 genes (a large genome) | gene tracks binned; `binned = true` rows |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | replicons chr1 100,000 circular; plasmid_1 8,000 circular | two `sector` rows (`start` 1, `end` = length) |
| 2 | gene g1 chr1 99,801–100,300 `+` (crosses origin), `detected` | one `gene_plus` row `start` 99,801, `end` 100,300, `category` detected; two drawn pieces |
| 3 | peptide placements at chr1 starts 5, 9,999, 10,001, 25,000 | `peptide_density` rows: 1–10,000 → 2; 10,001–20,000 → 1; 20,001–30,000 → 1 |
| 4 | islands ISL_1 high p 0.01, ISL_2 low p 0.2, alpha 0.05 | one `island` row (ISL_1, category high) |
| 5 | variants v1 ALT high, v2 mixed candidate, v3 call REF (not in variants.tsv by schema) | two `alt_site` rows with categories high, candidate |
| 6 | 25,000 genes on one 5-Mb contig | gene rows `binned = true`, ≤ 20,000 per gene track |
| 7 | landmarks chr1 oriC 1, ter 50,001 | `oric` and `ter` rows |
| 8 | the validation run's `OUT/` | renders; no islands, no markers; title "0 of 87 genes detected (0 %); 0 islands at FWER ≤ 0.05" |

## Acceptance criteria

- [ ] `tests/unit/viz/test_v01_circos.py`: examples 1–7 on the data TSV; origin crossing; linear
  contig; the element count per track ≤ 20,000 (count of drawn patches per track via the figure).
- [ ] Example 8 as a slow-marked test or on a copied mini fixture.
- [ ] Render time on a 5-Mb simulated genome recorded (M9-26 benchmarks).

## Questions for the lead

1. §10.2 lists "prophage, plasmid, BGC and CRISPR arcs"; §10.1 also colours AMR. **Proposed
   default:** draw all five types.
2. Islands drawn: only those at FWER ≤ `stats.alpha.context` (the claimed ones). **Proposed
   default:** yes; non-significant scan windows are not islands.
3. Peptide density counts placement rows (a multi-locus peptide counts at each accepted locus).
   **Proposed default:** rows, because the question is "where in the genome"; the data TSV says so
   in the `category` column (`placements`).

## Lead decisions (2026-10-22)
1. All five feature-arc types are drawn: prophage, plasmid, BGC, CRISPR and AMR.
2. Only claimed islands (FWER ≤ `stats.alpha.context`) are drawn.
3. Peptide density counts placement rows; the data TSV's `category` column says `placements`.
4. Branch off main: M9-03 (support files) is not merged, so the oriC/ter and GC-skew tracks are tested
   through their empty states only; the populated tracks are pinned when M9-03 lands.

## Lead decisions on the tests (2026-10-22, after d82bebf)
5. **AMR without coordinates.** `regions.tsv` AMR rows with NA start/end are not drawn and get no data
   row in pass 1; the figure says "AMR sets without coordinates not drawn" when any exist. AMR rows with
   coordinates are drawn. (Backlog: draw per member gene from genes.tsv in pass 2.)
6. **TU start tick** is at the transcription start: the low coordinate on `+`, the high one on `−`.
7. **Abundance rows.** One row per gene; NA and imputed genes have `value = NA` (drawn white).
8. **Title plurals.** Singular at 1 ("1 island", "1 ALT site", "1 novel ORF"); the suffix lists only the
   marker tables present, in the order ALT sites, novel ORFs.
9. Omission-note wording is free; each names what is missing ("peptide", "GC skew", "oriC", tool names).
10. **Density bins.** Every bin gets a row, zero counts included, so the bins tile each sector.
11. **Novel ORFs with NA `start_longest`** use the called ORF's own `start`/`end`.
12. **islands.tsv absent:** the island clause is left out of the title.
13. **genes.tsv with no rows:** the framework empty state, "No genes in this run" (index decision 4(b)).
13a. (amends 13, 2026-10-22, after 934902e) A genes.tsv with no rows still draws the replicon circle
   (sectors, ticks, region arcs), titled "0 of 0 genes detected (NA)": the sector layout comes from
   preflight, and several tests render from a zero-gene fixture. Sector rows carry `value` = replicon
   length. (Backlog: revisit with M9-16's empty-state rules.)
