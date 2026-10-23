# M9-10: V3 `locus_<id>`: the locus renderer and locus selection

Milestone M9 · **pass 1** · depends on: M9-01 · estimate ≈ 340 code lines (+ ≈ 320 test lines) ·
not statistical

## Spec sections implemented

- §10.2 V3: "What happens at this prophage, BGC, island, flagged operon or novel ORF? gene arrows
  coloured by abundance (cividis; undetected grey); peptide pileup track; promoter (bent arrow) and
  terminator (lollipop) marks; region shading. At most `report.max_locus_plots` (50), ordered by q";
  pyGenomeViz; static SVG.
- §10.1 multi-instance naming `V03_locus_<region_or_flag_id>`; ADR-0012 (promoters and terminators
  at novel ORFs; the renderer is shared with V13, M9-12); §6.5 (origin crossing; strand); §7.4.5
  (operon flags, q at `stats.alpha.operon_flags`); §7.3.5(b),(d).

## Public interface

`src/giap/viz/locus.py`
```python
@dataclass(frozen=True)
class Locus:
    locus_id: str          # region_id | island_id | pair_id | orf_id
    kind: Literal["prophage", "bgc", "island", "operon_flag", "novel_orf"]
    contig: str
    start0: int; end0: int # 0-based half-open, end0 may exceed length on circular contigs
    error: float | None    # the ordering estimate (q, FWER p, class FDR); None = NA
    label: str             # e.g. "prophage_1, O − E −4.8, q 0.003"

def select_loci(view: RunView, cfg: GiapConfig) -> list[Locus]
def render_locus(locus: Locus, view: RunView, cfg: GiapConfig, scale: tuple[float, float]) -> Instance
```

## Selection and order (index decision 18)

| Kind | Rows | Span | Error |
|---|---|---|---|
| prophage, bgc | `regions_stats`, `stat = o_minus_e`, `subset = all`, `region_type` ∈ {prophage, bgc} | region | `q_value` |
| island | `islands` | island | `p_value` (FWER) |
| operon_flag | `operon_flags` rows with `q_value ≤ stats.alpha.operon_flags` | the TU containing `gene_a` (from `genes.tsv` `tu_id` → `annotation/tus.tsv`), or the two genes' span when they are in different TUs | `q_value` |

Novel ORFs are V13's (index decision 19). The pool is sorted by `error` ascending (NA last), then
kind in the table order, then `locus_id`; the first `report.max_locus_plots` are rendered. An
operon flag whose TU is already rendered for another flag is merged into that plot's label.

## Rendering

- **Window** = span ± 2 kb (index decision 20); clipped to [0, length) on linear contigs; on
  circular contigs the window is kept unwrapped (it may extend past `length` or below 0, drawn
  in unwrapped coordinates with the origin marked by a vertical dashed line).
- **Gene arrows** (pyGenomeViz feature track): every `genes.tsv` gene overlapping the window,
  arrow by strand, fill = cividis of measured `abundance_log2` on one scale shared by all locus
  plots of the run (`scale` = 2nd–98th percentile of measured detected abundance in the run);
  undetected `#BBBBBB`; `single_peptide_unestimated` `#E69F00` outline; imputed values drawn grey
  with an "imputed" hatch (index decision 13); gene name or locus tag as the label (≤ 12 labels,
  the rest unlabelled to avoid overlaps, M9-24).
- **Peptide pileup:** depth per base of accepted placements from `peptide_hits.tsv` in the window,
  as a filled step plot above the genes, y "peptides".
- **Promoters** (bent arrow, `#009E73`) from `annotation/promoters.bed`; **terminators**
  (lollipop, `#D55E00`) from `terminators.bed`; strand-aware.
- **Region shading** behind everything: `regions.tsv` regions overlapping the window in the
  §10.1 region colour (alpha 0.2); islands with their colour and hatch.
- Single width (3.5 in) for windows ≤ 20 kb, double width otherwise; height fixed.

Data TSV per instance `V03_locus_<locus_id>.tsv`: `track` enum{gene, pileup, promoter, terminator,
region, island}, `id`, `contig`, `start`, `end` (1-based closed, unwrapped), `strand`, `value`
(abundance, depth or score), `category`. Plus a single index file `V03_locus.tsv` listing the
selected loci: `rank`, `locus_id`, `kind`, `contig`, `start`, `end`, `error`, `rendered` (bool).

Title per instance: the `label`, e.g. "prophage_1: 12 % detected vs 31 % expected (q = 0.003)";
for islands "ISL_3 (low): mean log2 `<m>`, FWER p `<p>`"; operon flags "`<flag>` between
`<gene_a>` and `<gene_b>` (q = `<q>`)".

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no prophage/BGC rows, no islands, no flags (validation run) | no instance files; `V03_locus.tsv` header-only; report text "no loci to show" |
| pool larger than 50 | 50 rendered, the rest listed with `rendered = false` |
| `max_locus_plots = 0` | none rendered; index lists all |
| window crossing the origin | unwrapped drawing; origin line |
| a linear contig end inside the window | clipped |
| a region with `q_value` NA (`small_region`) | ordered last within NA |
| `peptide_hits` absent | no pileup track, noted |
| a locus with 200 genes (large BGC) | arrows drawn (≤ 20,000), labels capped at 12 |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | prophage_1 q 0.003; ISL_1 FWER 0.01; flag PAIR_7 q 0.02; prophage_2 q NA | order prophage_1, ISL_1, PAIR_7, prophage_2 |
| 2 | 60 candidate loci, `max_locus_plots` 50 | 50 instance file sets; 60 index rows |
| 3 | prophage_1 0-based [10,000, 35,000) on 100-kb circular chr1; genes A [7,500, 8,200), B [37,100, 38,000) | window [8,000, 37,000); gene A drawn (overlaps), gene B not drawn |
| 4 | island 0-based [99,000, 101,500) (crosses origin) | window [97,000, 103,500) unwrapped; origin line at 100,000 |
| 5 | gene measured log2 25 with scale (18, 24) | colour clipped to the top of cividis |
| 6 | two placements of one peptide covering 1-based 101–120 and 111–130 | pileup depth 1 on 101–110, 2 on 111–120, 1 on 121–130 |
| 7 | flags PAIR_7 and PAIR_8 in the same TU | one instance, label naming both |

## Acceptance criteria

- [ ] `tests/unit/viz/test_locus_select.py`: examples 1, 2, 7; NA ordering; the cap reads
  `report.max_locus_plots`.
- [ ] `tests/unit/viz/test_locus_render.py`: examples 3–6; origin crossing; linear clipping;
  shared colour scale across instances.
- [ ] File names `V03_locus_<id>.{svg,png,pdf}` and `data/V03_locus_<id>.tsv`; IDs sanitised.

## Questions for the lead

1. Locus pool (index decision 18): novel ORFs excluded from V3 because V13(b) renders them.
   **Proposed default:** excluded.
2. Window flank 2 kb (index decision 20). **Proposed default:** 2 kb.
3. Ordering across kinds mixes q-values, FWER p-values and (none here) class FDRs. **Proposed
   default:** a single list by the estimate itself, as §10.2's "ordered by q" implies; the index
   TSV records which estimate each row used (`kind`).

## Lead decisions (2026-10-22)
1. Novel ORFs are excluded from the V3 pool; V13(b) renders them (index decision 18).
2. Window flank is 2 kb (index decision 20), a module constant.
3. One list ordered by each row's own estimate (q or FWER p), ascending, ties broken by the
   index TSV's stable order (kind, then id); the index TSV's `kind` column records which estimate each
   row used.

## Lead decisions on the tests (2026-10-22, after 39f9ead)
4. **Index channel.** The index is an `Instance` with `suffix=None`, written by the framework to
   `plots/data/V03_locus.tsv` (not `plots/V03_locus.tsv`), columns exactly `rank, locus_id, kind,
   contig, start, end, error, rendered`. No loci gives the framework empty state "no loci to show".
   The index figure it also writes is a plain ranked list; its content is not pinned (backlog: decide in
   M9-16 whether the report links it).
5. `select_loci` returns the whole ordered pool; the plot renders the first `report.max_locus_plots`.
6. **Labels.** The Title-section form wins over the `Locus` docstring: "prophage_1: 12 % detected vs
   31 % expected (q = 0.003)", `%g` numbers, "(q = NA)" when missing; islands "ISL_1 (high): mean log2
   21.5, FWER p 0.01". A merged operon locus keeps the best flag's `pair_id` and q and lists each clause,
   joined by "; ".
7. Data-TSV rows carry each feature's own coordinates (unwrapped on circular contigs, so starts can be
   negative), not clipped to the window. Pileup rows are maximal runs of constant non-zero depth.
8. NA loci sort last, then by kind and id; `small_region` gets no extra rule.
9. pyGenomeViz may be used only for its patches/tracks drawn onto our own `Figure`; no pyplot figure
   may survive a render.

## Lead decisions on the review round (2026-10-22, #209)
10. An operon flag across the origin spans the shorter arc (b unwrapped relative to a).
11. The title wraps and stays inside the figure; every clause with its q stays in the drawn text.
12. A colorbar labelled "log2 abundance" shows the run-wide cividis scale, plus a legend outside the
    gene and pileup axes with undetected, imputed and single_peptide_unestimated swatches.
13. Strands: arrow heads point downstream. The promoter bent arrow stands at its TSS, the downstream end
    of the upstream Promotech window (high coordinate on +, low on −); + marks above the gene-track
    centre, − below. (Corrects the review-round brief, which said low on +.)
14. Imputed genes keep their data-TSV row with `value` NA and `category = imputed`.
