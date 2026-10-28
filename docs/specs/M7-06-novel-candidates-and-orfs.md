# M7-06: Novel-ORF candidates, ORF reconstruction, grouping and the second-round windows

Milestone M7 · track N (novel ORFs) · **pass 1** · depends on: merged M3-14 (`LocusClassifier`,
`in_frame`, the `cterm_extension` class of ADR-0013), M3-15 (`match_contaminant_peptides`),
M2-13 (`gene_windows`, `dedupe_and_chunk`), M2-16 (`ostir_row`, `write_ostir_input`) · estimate
≈ 360 code lines (+ ≈ 550 test lines) · not statistical

## Spec sections implemented

- §7.6.1: "A candidate is an accepted placement whose locus class is `nterm_extension`,
  `intergenic`, `cds_outframe`, `antisense`, `pseudogene` or `sorf`. Peptides that also map in
  frame to an annotated gene (`cds_inframe` at any locus), or to a contaminant, never support a
  novel ORF." ADR-0013 decision 2 extends the list by `cterm_extension`.
- §7.6.2: "**Downstream.** Each candidate ORF extends in frame to the next stop codon.
  **Upstream.** Candidate start codons (`annotation.start_codons`) lie between the previous in-frame
  stop and the first peptide. **Two starts are reported:** `start_longest`: the furthest upstream
  start; `start_best_ostir`: the candidate with the highest OSTIR initiation score.
  **Second-round layers.** Every candidate start gets targeted Promotech windows (−200 to +20) and
  an OSTIR window (−50 to +50; for several candidates, one window spans them all with `start`/`end`
  ranges). **Grouping.** Peptides are grouped per ORF, and every ORF contains all its peptides."
- §7.1.4 (targeted windows: 40 nt, 3′ ends at S + k, k ∈ [−200, +20], dedupe, 100,000-window
  chunks), §7.1.6 and the §7.1.1 OSTIR gotchas ("`start_position` is 1-based in the input
  sequence. Restrict … with CSV `start`/`end` columns"), §6.2 step `novel_candidates` ("loci →
  candidate ORFs, extra Promotech and OSTIR windows"), §6.5 (IDs `NOVEL_<contig>_<start1>_<strand>`),
  §7.6.3 (the class column only; the evidence rules are M7-09).

## Public interface

`src/gippa/novel/candidates.py` (pure)
```python
NovelClass = Literal["nterm_extension", "cterm_extension", "intergenic", "cds_outframe",
                     "antisense", "pseudogene", "sorf"]

@dataclass(frozen=True, slots=True)
class Orf:
    orf_id: str                      # NOVEL_<contig>_<start1>_<p|m>
    cls: NovelClass
    contig: str
    strand: Literal["+", "-"]
    frame: int
    stop: int                        # lowest 0-based coordinate of the stop codon (wrapped)
    starts: tuple[int, ...]          # candidate start codons, 5'→3' (lowest coordinate of each)
    start_longest: int | None
    end_5p_open: bool                # upstream walk hit a linear contig end
    end_3p_open: bool                # downstream walk hit a linear contig end
    peptides: tuple[str, ...]        # distinct, I/L-folded, sorted
    best_q: float | None             # min input q of its peptides (NA for de novo)
    overlapping_gene: str | None

def candidate_rows(peptide_hits: pd.DataFrame, peptides: pd.DataFrame,
                   contaminants: Sequence[ProteinRecord] | None) -> pd.DataFrame: ...
def reconstruct(rows: pd.DataFrame, gm: GeneModel, genome: Mapping[str, str],
                cfg: GippaConfig) -> list[Orf]: ...
def promoter_windows(orfs: Sequence[Orf], gm: GeneModel, genome, cfg) -> list[Window]: ...
def ostir_rows(orfs: Sequence[Orf], gm: GeneModel, genome, anti_sd: str, cfg) -> list[OstirRow]: ...
```

## Rules

1. **Candidate rows.** Rows of C11 `peptide_hits` with `locus_class` in the seven classes above.
   A peptide is excluded entirely when any of its rows is `cds_inframe` or `rna_overlap`, when C9
   marks it `is_contaminant`, or when it is an exact substring (I = L) of a contaminant protein
   (decision 1: placed peptides were never checked against the library, M3-15 only checks unplaced
   ones). Peptides with more than one candidate locus are excluded too (decision 2).
2. **Downstream.** From the most downstream peptide, extend in the peptide's frame to and including
   the first stop codon (table 11). Circular replicons wrap (unwrapped coordinates, §6.5); a linear
   contig end without a stop sets `end_3p_open` and the ORF ends at the last whole codon.
3. **ORF key and grouping.** ORF identity is `(contig, strand, stop)`. All candidate peptides with
   the same key form one ORF; that makes "every ORF contains all its peptides" true by
   construction, because the start is searched upstream of the most upstream one (rule 4).
4. **Upstream.** Walk in frame from the most upstream peptide's first codon back to the previous
   in-frame stop (exclusive). Candidate starts: every codon in `annotation.start_codons` in
   `(previous stop, first peptide codon]` — the first peptide's own first codon is eligible
   (an N-terminal peptide that begins with the initiator M). `start_longest` = the furthest
   upstream. A linear contig start without a stop sets `end_5p_open`. No candidate start →
   `starts = ()`, `start_longest = None`, and the flag `no_start_codon` (decision 5).
5. **Class.** The highest-priority M3-14 class among the ORF's peptides (order `pseudogene` >
   `sorf` > `nterm_extension` > `cterm_extension` > `cds_outframe` > `antisense` > `intergenic`;
   decision 3). For `nterm_extension`, the ORF is the annotated gene's own reading frame and stop;
   `overlapping_gene` = that gene. For `cds_outframe`/`antisense`/`pseudogene`/`sorf`/
   `cterm_extension`, `overlapping_gene` = the gene of the peptides' locus; for `intergenic`, NA.
6. **ID.** `NOVEL_<contig>_<start1>_<p|m>` with `start1` = the 1-based lowest coordinate of the
   **longest** ORF interval (§9.2's `<start>`; on `−` this is the stop codon's first base, so the ID
   is stable on both strands). Without a start codon the interval starts at the first peptide
   (decision 4).
7. **Second-round Promotech windows.** For every candidate start S of every ORF, the 40-nt windows
   on the ORF's strand whose 3′ ends lie at S + k, k ∈ `annotation.promoters.window_3p_range`
   (`gene_windows` generalised to a start position), deduplicated by `(contig, strand, end3p)` and
   chunked as §7.1.4. Header `w<index>`, the `windows_map` columns of M2-13 plus `orf_id`.
8. **Second-round OSTIR row**, one per ORF with ≥ 1 start: `id = orf_id`; `seq` from 50 nt upstream
   of the most upstream start to 50 nt downstream of the most downstream start
   (`annotation.ostir.upstream_nt/downstream_nt`), strand-aware, wrapping or clipping as
   `ostir_row`; `start` and `end` = the 1-based positions of the first and last candidate start in
   `seq`; anti-SD from C4 `annotation_meta.json`.
9. **Peptide fields.** `peptides` distinct after I/L folding; `best_q` = min C9 `q_value` (the gene
   table's `best_peptide_q` rule, M4-11); NA when no peptide has one (de novo).
10. **Minimal call (pass 1 only, replaced by M7-09 and M7-11):** `called` = distinct peptides ≥
    `novel_orfs.min_peptides_without_decoys` (2); `class_fdr = unestimated`; flag
    `class_fdr_unestimated`. This is §7.6.4's "Neither" branch applied to every input until the
    FDR branches exist, so pass 1 never claims error control it does not have.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no candidate rows | `[]`; no windows, no OSTIR rows |
| peptide also `cds_inframe` elsewhere | excluded (rule 1) |
| contaminant substring | excluded (rule 1) |
| ORF crossing the origin | unwrapped coordinates; ID from the 1-based lowest coordinate (may be near `length`) |
| upstream walk reaches a linear contig start | `end_5p_open`, starts searched from the contig start (flag `partial`) |
| downstream walk reaches a linear contig end | `end_3p_open` (flag `partial`) |
| no start codon before the first peptide | rule 4; no Promotech windows, no OSTIR row for that ORF |
| two peptides, same frame and strand, an in-frame stop between them | two ORFs (different stops) |
| `nterm_extension` peptide but an in-frame stop between it and the annotated start | cannot occur (M3-14 requires "after the previous in-frame stop") |
| `pseudogene` peptides in two frames of one pseudogene | two ORFs (two stops) |
| candidate start window running off a linear contig | Promotech window dropped (as M2-13); OSTIR window clipped (`ostir_row`) |
| circular replicon with no stop in the frame for a whole turn | `ValueError` naming contig and frame (impossible in a real genome) |

## Examples

Circular `c1`, 1,000 nt, `+` strand, frame 0 unless stated.

| # | Input | Expected |
|---|---|---|
| 1 | in-frame stop TAA at [30, 33); ATG at [42, 45); GTG at [60, 63); peptides P1 at [72, 102) and P2 at [90, 111), intergenic; next stop TGA at [120, 123) | one ORF, key (c1, +, 120); `starts` (42, 60); `start_longest` 42; ID **`NOVEL_c1_43_p`**; length from 42: (120 − 42)/3 = **26 aa**; 2 peptides → `called` (pass 1) |
| 2 | as 1 plus P3 in frame 1 at [151, 181), its stop at [190, 193) | a second ORF, key (c1, +, 190) |
| 3 | as 1, but P2 also has a `cds_inframe` locus elsewhere | P2 excluded; ORF holds P1 only → not called |
| 4 | the mirror of 1 on `−` (reverse complement of the contig, coordinates x → 1000 − x) | key (c1, −, 877); the ORF interval [877, 958); `start_longest` = the ATG at [955, 958); ID **`NOVEL_c1_878_m`** |
| 5 | linear 300-nt contig, peptide at [6, 36) in frame 0, no stop upstream | `end_5p_open`; starts searched in [0, 6]; flag `partial` |
| 6 | gene G `+` [300, 900); in-frame stop at [240, 243); TTG at [249, 252); peptide at [255, 285) classed `nterm_extension` | class `nterm_extension`; stop = G's stop at [897, 900); `start_longest` 249; `overlapping_gene` G |
| 7 | no start codon in (33, 72] for example 1's frame | `starts` (), `start_longest` NA, flag `no_start_codon`; ID from the first peptide: `NOVEL_c1_73_p` |
| 8 | example 1's Promotech windows | 221 windows for S = 42 and 221 for S = 60 (3′ ends 42 − 200 … 42 + 20, wrapped), deduplicated: the union has **239** distinct 3′ ends (the ranges [−158, 62] and [−140, 80] overlap) |
| 9 | example 1's OSTIR row | `seq` = genome[−8, 110) wrapped (118 nt); `start` **51**, `end` **69** (60 − 42 + 51) |
| 10 | a contaminant peptide (exact substring of a cRAP keratin) in example 1 | excluded; the ORF has 1 peptide |

## Acceptance criteria

- [ ] `tests/unit/novel/test_candidates.py`: examples 1–10 and every degenerate row; §7.6.6 tests
  "ORF reconstruction on both strands, including circular wrap", "exclusion of annotated-gene and
  contaminant peptides".
- [ ] Property: mirroring the genome maps every ORF to its mirror (strand, stop, starts) and its ID
  to the mirrored ID.
- [ ] Property: every peptide of an ORF lies between `start_longest` (or the first peptide) and the
  stop, in the ORF's frame, with no in-frame stop between start and stop.
- [ ] Property: a planted ORF on a random genome with ≥ 2 peptides cut from it is reconstructed
  with the planted stop, and its planted start is in `starts`.
- [ ] The second-round windows equal M2-13's `gene_windows` for a start S when S is an annotated
  gene's start (one generalised function, not two).

## Calibration and attacks owed

None here; class FDR (M7-11) and the novel-ORF F1 (M7-13) measure it.

## Hardening owed (pass 2)

- Back-fill `novel_candidate_id` in `mapping/peptide_hits.tsv` (§9.3 "Assignment: `gene`,
  `novel_candidate_id`"; M3-14 left it NA "novel-locus clustering is Module 6") — M7-15, Q6.
- Golden candidate table for the 50-kb fixture.

## Questions for the lead

1. Placed peptides are never checked against the contaminant library (M3-15 checks only unplaced
   ones), so §7.6.1's contaminant exclusion needs its own check. **Proposed default:** exact
   substring (I = L) against the same library, plus C9 `is_contaminant` (engine IDs with a
   contaminant prefix), inside `novel_candidates`.
2. Multi-locus candidate peptides (two intergenic loci are `shared`, M3-14). **Proposed default:**
   excluded from every ORF — a peptide that fits two unannotated loci equally well cannot raise
   either above the ≥ 2-peptide or FDR bar. The alternative (count at both) double-counts.
3. ORF class when its peptides carry different locus classes (e.g. an ORF running out of a CDS into
   intergenic sequence). **Proposed default:** highest priority among the peptides' classes in the
   M3-14 order (rule 5). Classifying the ORF interval itself instead would make the class depend on
   the start choice.
4. `NOVEL_<contig>_<start1>_<strand>`: which "start"? **Proposed default:** the 1-based lowest
   coordinate of the longest ORF interval (rule 6) — independent of OSTIR (so IDs do not move when
   the second round changes), and on `−` it is the stop side, which is fixed.
5. New flags in the Novel-ORF scope: `no_start_codon`, `partial` (reused from the gene scope), and
   M7-09's `rbs_threshold_unestimated`, `promoter_supported`, `terminator_supported` (ADR-0012),
   `readthrough_candidate` (ADR-0013). **Proposed default:** add through the schema as M6 did
   (M6 lead decision 6); note ADR-0011/0013 took vocabulary additions to Khoa, so the lead decides
   whether this needs Khoa.
6. `novel_candidate_id` in `peptide_hits.tsv`: back-fill it (that file is published by LOCI, before
   Module 6 runs) or leave it NA? **Proposed default:** leave NA in pass 1; in M7-15 the
   `novel_orfs` step rewrites the published `mapping/peptide_hits.tsv` with the column filled (one
   writer, `gippa.mapping.outputs`), or `finalize` does it when it exists.

## Lead decisions

1. **Questions 1–4 and 6:** the proposed defaults are accepted.
   - Contaminants are checked by exact I = L substring match against the library, plus C9 `is_contaminant`, inside `novel_candidates`.
   - Multi-locus peptides are excluded from every ORF.
   - An ORF's class is the highest-priority class among its peptides' classes.
   - The ID uses the 1-based lowest coordinate of the longest interval.
   - `novel_candidate_id` stays NA in pass 1.
2. **Question 5:** the Novel-ORF flags are added through the schema, as M6 did (M6 lead decision 6). This means `no_start_codon`, `partial` (reused), `class_fdr_unestimated`, and M7-09's flags. Khoa is told in the integration-pass summary. They are recorded in the run report, not taken as an ADR, because they describe existing conditions and add no new claims.
3. **The test author's choices** (commit 50b2d45) are accepted as the contract.
   - **`minimal_call(orfs, cfg)`** is the rule-10 entry point. It returns the columns `orf_id, n_peptides, called, class_fdr = "unestimated", flags`, with flags sorted and joined by `;`.
   - **`candidate_rows`** returns the retained C11 rows with all of their columns.
   - **Coordinates:** the 0-based lowest codon coordinate, wrapped. `starts` run 5'→3' on the strand.
   - **`frame`** is the peptides' C11 frame label.
   - **`best_q`** is `None` when no peptide has a q-value, otherwise the minimum over the raw strings.
   - **Multi-locus** is decided by distinct (contig, start, end, strand).
   - **`overlapping_gene`** for `cds_outframe` comes from the gene model, on the same strand.
   - **Promotech windows** are unique by (contig, strand, end3p) across ORFs, with `owner = orf_id`.
   - **Whole turn without a stop:** `ValueError` naming the contig and the frame.
   - **Exclusion reasons:** no output column in pass 1.
