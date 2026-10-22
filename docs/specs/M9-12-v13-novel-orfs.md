# M9-12: V13 `novel_orfs`: which new genes or upstream extensions are supported?

Milestone M9 · **pass 1** · depends on: M9-10 (the locus renderer) · estimate ≈ 160 code lines
(+ ≈ 200 test lines) · not statistical

## Spec sections implemented

- §10.2 V13: "(a) class counts by call status; (b) locus views of the top 20 (V3 renderer)";
  static SVG + table in the report.
- §7.6 (classes, `class_fdr` "float or `unestimated`", `called`, `rank` §7.6.5), ADR-0012
  (`promoter_supported`, `terminator_downstream`, drawn by the V3 renderer), ADR-0013
  (`cterm_extension`), §9.3 `novel_orfs.tsv` (main; M7-09 evidence columns already in the schema,
  values filled by PR #180).

## Panel (a): class counts by call status

Bars per `class` in the schema enum order (nterm_extension, cterm_extension, intergenic,
cds_outframe, antisense, pseudogene, sorf), stacked by status: `called` (`called = true`,
`#0072B2`), `not called` (`called = false` and `class_fdr` finite, `#BBBBBB`), `unestimated`
(`class_fdr` NA or flag `class_fdr_unestimated`, `#E69F00`). Each bar labelled with its class FDR
range among called rows ("FDR ≤ 0.008") or "unestimated".

## Panel (b): locus views (`V13_novel_orfs_locus_<orf_id>`)

The 20 `called = true` rows with the smallest `rank` (§7.6.5), each rendered with
`giap.viz.locus.render_locus` as a `Locus(kind="novel_orf")`: span from `start_longest` and `end`, which
are in transcription direction (on `−`, `start_longest > end`; the schema says so), converted to
0-based half-open [min − 1, max) on read; `start_longest` NA → the span starts at the first peptide, window ± 2 kb; the ORF drawn as a highlighted
arrow (`#CC79A7` outline) above the annotated genes; a second tick at `start_best_ostir` when it
differs from `start_longest`; peptides from `peptides` (list) highlighted in the pileup track;
promoter/terminator marks as V3. Title: "`<orf_id>` (`<class>`): `<n>` peptides, best q `<q>`,
class FDR `<f>`|unestimated".

## Data TSVs

- `V13_novel_orfs.tsv`: `class`, `status` enum{called, not_called, unestimated}, `n`,
  `class_fdr_max_called`; sorted by class (enum order), status.
- `V13_novel_orfs_locus_<orf_id>.tsv`: the M9-10 locus TSV layout, plus `track = novel_orf`.
- `V13_novel_orfs_locus.tsv`: index of the top-20 selection (`rank`, `orf_id`, `rendered`).

## Title (main figure)

"`<c>` novel ORFs called (class FDR ≤ `<f>`), `<u>` unestimated" or "No novel ORF called; `<u>`
candidates unestimated".

## Degenerate cases

| Situation | Behaviour |
|---|---|
| `novel_orfs.tsv` header-only | empty state "No novel ORF candidates in this run"; no locus files |
| fewer than 20 called | that many locus views |
| no called ORF but unestimated candidates | panel (a) only; no locus views |
| ORF crossing the origin | unwrapped window (M9-10) |
| a `−`-strand ORF read as if `start < end` | the bug this table guards against: example 3 |
| `class_fdr` written as the literal `unestimated` (§9.3 "float or `unestimated`") | read as NA + status unestimated |
| PR #180 not merged (evidence columns NA) | panel (a) unchanged; locus views without promoter/terminator emphasis |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | 3 intergenic rows: called (FDR 0.004), called (0.008), not called (0.2); 2 antisense with `class_fdr` `unestimated` | (a) intergenic called 2, not_called 1; antisense unestimated 2; title "2 novel ORFs called (class FDR ≤ 0.008), 2 unestimated" |
| 2 | 25 called rows ranks 1–25 | 20 locus files, ranks 1–20; index has 20 rows |
| 3 | a novel ORF on chr1 (100 kb, circular), strand `−`, `start_longest` 1,500, `end` 1,201 (transcription direction) | span 0-based [1,200, 1,500); window [−800, 3,500) on a circular contig, drawn unwrapped |
| 4 | header-only table | empty state; no `V13_novel_orfs_locus_*` files |

## Acceptance criteria

- [ ] `tests/unit/viz/test_v13_novel_orfs.py`: examples 1–4; the literal `unestimated` parsed;
  locus files named from `orf_id` (sanitised).
- [ ] Uses M9-10's renderer (no second locus implementation).

## Questions for the lead

1. Top 20 by `rank` among `called` rows only (index decision 19). **Proposed default:** yes; the
   uncalled are not claims.
2. Example 3: a minus-strand ORF near position 1 on a circular contig has a window below 0.
   **Proposed default:** draw unwrapped (negative coordinates shown as `length + x` tick labels);
   linear contigs clip at 0.
