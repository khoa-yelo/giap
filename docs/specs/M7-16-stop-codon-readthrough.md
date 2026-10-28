# M7-16: Peptides that span a stop codon (Sec, read-through, `cterm_extension`)

Milestone M7 · track R (mapping) · **pass 2 · option B accepted by Khoa (2026-10-24, "B"); ADR-0058**
· depends on: merged M3-06/M3-07/M3-08 tolerant routes, M3-10 placement FDR, M3-14 loci, M3-21
placement suites, M3-23 makedb; M4-16 detection summary; M7-01, M7-05, M7-07, M7-09, M7-11b,
M7-12, M7-18 · estimate ≈ 544 code lines in three PRs (+ ≈ 185 simulator lines, + ≈ 1,700 test
lines; 16a merged as #331) · **stats**: the
statistician reviews this spec and the `placement-fdr.md` addition **before any code** (§14.1);
reviewed 2026-10-24 at e7f7f744 (APPROVE WITH CHANGES), and the required changes 1–9 are applied
in the rules, test plan and acceptance criteria below · **amended 2026-10-24 evening** with
Khoa's decisions for 16b and the 16c labels (section "Amendment: 16b reporting and the 16c
labels") · **16c amendment 2026-10-25** after the statistician's pre-code check (section "16c
amendment (2026-10-25)")

The decision-gated text (the gap, the options, the original questions) is kept unchanged below as
history. The buildable spec starts at "Accepted: option B".

## History: the decision-gated version (2026-10-21)

### The gap

Three accepted parts of the design assume that a peptide can be placed across a genome stop codon:

- §7.5.3 G8 `stop_codon`: "The REF residue is a stop (e.g., Sec read through UGA)", and §7.5.1 "Sec
  sites … appear as ALT calls" under `search_db = same`;
- ADR-0013: `cterm_extension` evidence "a peptide spanning the annotated stop"; the
  `readthrough_candidate` label; the patch row with a new end;
- §11.6 TE2's exception list (`stop_codon` sites).

The merged matcher cannot do it. The six-frame index writes `*` at stops (§7.2.1), the exact
automaton never matches across one, and the M3-06 DP states "A target `*` (`STOP`) … is never
consumed". M3-22 found exactly this: the 4 TE1 misses are selenocysteine (fdnG) peptides. So
`stop_codon` and `cterm_extension` calls are unreachable in v1 as built (M7-05 rule 6, M7-09
rule 4).

### Options (for the decision; the spec-writer does not choose)

| Option | What changes | Cost | Risk |
|---|---|---|---|
| A. **Leave it** | nothing; document `stop_codon`/`cterm_extension` as reachable only through frameshift splits; TE2's `stop_codon` clause is vacuous | 0 | the G8 flag and ADR-0013 stay partly dormant |
| B. **Annotated-stop substitution** | the tolerant DP may consume a target `*` **only at the annotated stop codon of a CDS/sORF**, at cost 1, recorded as `substitution` with REF `*`; decoys get the same rule at the mirrored positions | ≈ 220 lines + re-calibration of M3-21 placement suites | small: one position per gene; decoys keep symmetry |
| C. **Any-stop substitution** | `*` consumable anywhere at cost 1 | ≈ 150 lines | large: every stop becomes a chance-hit site; the M3 chance-hit table (B.2) and placement FDR need re-derivation |

**Proposed default: B**, as an ADR (it changes Module 2's search space and therefore placement
FDR; Khoa accepts it). With B: `makedb --patch` writes the read-through residue of a
`cterm_extension` patch as `X` (decision for the ADR); the M7-12 simulator plants Sec/read-through
(the backlog line); TE1's fdnG misses become `stop_codon` ALT sites and TE2 exercises its second
exception.

### Acceptance criteria (if B is accepted), superseded by the section below

- [ ] The M3 tolerant property tests (zero false negatives vs brute force) extended to the
  annotated-stop rule; the M3-21 placement calibration re-run green (fast + full).
- [ ] `stop_codon` ALT sites on the Nissle TE1/TE2 set at fdnG (and any other Sec genes).
- [ ] `cterm_extension` reachable and patchable; ADR-0013 consequences tested end to end.

### Questions for Khoa (through the lead), answered

1. Option A, B or C? **Proposed default:** B, by ADR.
2. With B, how does `makedb` translate the read-through codon of a patch? **Proposed default:** `X`
   (search engines treat it as any residue), since Sec (`U`) and Trp (`W`) read-through both occur.

   **Answers (Khoa, 2026-10-24):** 1. B. 2. `X`. Recorded in ADR-0058.

---

## Accepted: option B (Khoa, 2026-10-24)

ADR-0058 is the decision record. Its items 1–5 are in force. Item 6 (annotated internal Sec/Pyl
codons) was **approved by Khoa on 2026-10-24 (evening)** (Q1 below). PR 16c holds it, and it is
next after 16b. ADR-0058's status line and item 6 were updated in 16b (lead, db7d92ab, 933f0486).

**One finding changes the motivation.** The 4 TE1 misses are **not** at an annotated stop. Bakta
annotates both fdnG copies as full-length selenoproteins with `U` over an *internal* TGA (residue
196 of 1,015 and of 1,016 aa). Option B as accepted therefore does not place them (example 2). It
places C-terminal read-through (`cterm_extension`) and any residue read at a terminal stop. Placing
fdnG needs ADR-0058 item 6.

## Spec sections implemented

- §7.2.4 cost table ("Each of these costs one difference: a substitution …") extended by ADR-0058
  item 1. G4 completeness ("zero false negatives against a brute-force matcher … for every
  difference type", §11.2) now holds **with** the consumable-stop rule.
- §7.2.1 decoys ("Decoy hits are only ever counted, never placed") with the same rule at mirrored
  positions (ADR-0058 item 2).
- §7.2.5 G3: estimator unchanged. R3: no tolerance, bound, α or R changes.
- §7.5.1 "Sec sites … appear as ALT calls", §7.5.3 G8 `stop_codon` (M7-05 rule 6 becomes reachable).
- ADR-0013 items 2–3: `cterm_extension` evidence, `readthrough_candidate`, the patch row; §7.9.2
  `--patch` with the residue `X` (ADR-0058 item 5).
- §11.6 TE1/TE2 (re-run, expectations below); backlog lines "the simulator plants stop-codon
  readthrough … for cterm_extension calibration" (M7-12 Q1) and "M7-16: decide `stop_codon` and
  `mass_equivalent` for blocks and swaps whose REF contains `*`" (closed: never reachable).

## Definitions

- **Stop site.** `StopSite(contig, strand, codon_start, locus_tags, codon, kind)`. `codon_start` is
  the lowest 0-based coordinate of the codon, wrapped modulo the replicon length on circular
  replicons. `codon` is one of `TAA`/`TAG`/`TGA` on the gene strand. `kind` is `terminal`, or
  `recoded` (16c only). For a `+` gene with interval `[s, e)` (unwrapped, `e` may exceed the
  length `n`) the terminal codon is `[e − 3, e)`, so `codon_start = (e − 3) mod n`. For a `−` gene
  it is `[s, s + 3)`. Features: `feature_type` in {`cds`, `sorf`}; pseudogenes never.
- **Marked cell.** A frame position `p` (aa index) whose residue is `*` and whose
  `index.to_interval(fid, p, p + 1)` is the stop site's codon. In code the tolerant search text
  carries `READTHROUGH = 22` there (a new target-only code next to `STOP = 19`, `AMBIG = 20`,
  `NOMATCH = 21`). The six-frame frames themselves and the exact route keep `19`.
- **Mirrored cell.** `reverse` decoy: `n_aa − 1 − p`. `shuffle` decoy: `p`. Both always hold `*`.

## Rules

1. **Building the stop sites** (`annotated_stops`, pure). Take every CDS/sORF of the gene model.
   Read the terminal codon on the gene strand, across the origin when it crosses it. If it is not
   a table-11 stop, skip the feature and count it under `not_a_stop` (for example a 3′-partial gene
   at a linear contig edge). Several features with the same `(contig, strand, codon_start)` form one
   site, with `locus_tags` sorted. Output sorted by (contig order, strand, codon_start).
   **16c:** with `recoded=True` the internal `U`-over-TGA and `O`-over-TAG codons are added as
   `recoded` sites (rule 21: residue → codon map, length check, skip-and-count, merge).
2. **Locating the cells** (`frame_cells`, pure). Let `x` range over the codon's positions in the
   frame's extended sequence: `codon_start`, plus `codon_start + n` when that fits inside `ext_nt`
   on a circular replicon.
   - On `+`, the cell is `p = (x − offset) / 3` in the frame with `offset = x mod 3`.
   - On `−`, it is `p = (ext_nt − offset − 3 − x) / 3` in the frame where that is an integer ≥ 0.
   - Every cell is checked: `to_interval(fid, p, p + 1)` must equal the site's codon (start wrapped)
     and `frame.aa[p] == STOP`. A failure raises `ValueError` naming the site; it is a bug, not
     data. Linear replicons have one cell per site, and circular ones one or two.
3. **Decoy cells** (`mirror_cells`). For every decoy type built, the cells of each frame are mirrored
   (Definitions). Cell counts per frame equal the target's. Every mirrored cell must hold `STOP`, or
   the function raises. `mirror_cells` is deterministic and draws nothing from any RNG (§6.6): the
   shuffle decoy keeps its positions, and the decoy frames themselves are built exactly as in
   format 1, so for the same seed they are byte-identical across index formats 1 and 2.
4. **The DP** (`tolerant_dp._verify_core`, `verify`, `verify_batch`). Only the (1, 1)
   move gains one case. When `t1 == READTHROUGH` it is a substitution, cost `w_sub`, rank 0, for
   any query code, `NOMATCH` included. Every other move keeps `t0 < 19` and `t1 < 19`. So a
   swap, block or near-block next to or over a marked cell is impossible. A marked cell never
   matches. Without `READTHROUGH` codes the kernel is unchanged, and so are its results.
5. **Seeds.** No change. The q-gram index skips codes ≥ 19, and a marked cell is one substitution
   that destroys at most q of the (q + 1)-per-difference q-grams the lemma allows. The threshold
   `t = (L − q + 1) − k(q + 1)` stays lossless.
6. **The tolerant route** (`match_tolerant`). It gets a new keyword `cells: StopCells | None =
   None`. `Search(frames, marked=…)` writes `READTHROUGH` into its concatenated text at the target
   cells (target search) or at the mirrored cells of the search's decoy type (decoy search). Rows,
   sorting and de-duplication are as today. `diff_types` stays `substitution`, with no new column.
   `cells=None` or empty cells give byte-identical output to main.
7. **The frameshift route** (`match_frameshift`). It gets the same keyword. `_FSearch` builds its
   `Search` with the marks, so the side DP `_verify_anchor` sees code 22 and never consumes it
   (22 ≥ 19). The precedence check `nonfs_overlap` gains rule 4's (1, 1) case. Effect: a split is
   suppressed by an overlapping non-frameshift read-through alignment at a cost ≤ the split's cost.
   A decoy frameshift search builds its `_FSearch` with **that decoy type's mirrored cells**, never
   the target cells, so `nonfs_overlap` suppresses a decoy split exactly as it suppresses the
   corresponding target split.
8. **Mapping steps.** `sixframe_index` gains input `genes` (the M2-04 `genes.parquet`; mapping.nf
   passes `gm[0]`). It builds the stop sites and cells, and writes them into the index (format 2):
   arrays `target_rt`/`target_rt_bounds` and `{kind}_rt`/`{kind}_rt_bounds` per decoy type, plus a
   `stop_sites` JSON list in `meta`. `load_index` keeps its signature. A new `load_stop_cells(path)
   -> StopCells` reads the cells, and format-1 files raise `ValueError` ("rebuild the six-frame
   index"). `map_tolerant` loads the cells and passes them to both routes. `sixframe_index` emits
   one INFO event `mapping.stop_sites` with `data = {"terminal": n, "recoded": 0, "not_a_stop": n,
   "cells_target": n, "cells_decoy": {kind: n}}` (16c: `recoded` counts the recoded sites, the
   cell totals include them, and `recoded_skipped` is added; rule 22). An empty gene model gives zero sites, an
   unchanged search and the event with zeros.
9. **Placement FDR, loci and outputs.** No code change. Read-through hits are tolerant rows at
   `n_diff ≥ 1` and enter strata 1–2 with π̂ = 1. M3-14 classifies them: crossing the annotated 3′
   end in frame → `cterm_extension` (`_cterm`, unchanged), and a window that ends on the stop
   codon itself → `cds_inframe` (Q3).
10. **M7-01 re-alignment** (`variants.align`). `align_hit`/`align_hits` gain `stop_sites:
    frozenset[tuple[str, str, int]] = frozenset()`, keyed `(contig, strand, codon_start)`.
    `_residues` accepts a `*` only at a codon in that set and still raises on any other `*` or `X`.
    The traceback DP `_dp` encodes such a `*` as code 22 and allows only move 0 at cost 1 against
    it. The recomputed cost must still equal `n_diff` (`_check`, unchanged). `variants_step` builds
    the set with the same `annotated_stops` from its gene model.
11. **Variants calls** (`variants.pileup`, `vcf_records`, `flags`).
    - A substitution op with REF `*` and a **standard** ALT residue is an ordinary ALT (M7-02 rules
      4–7). `flags_for` already adds `stop_codon` when `ref_aa == "*"`, and `nt_class` comes from
      the stop codon (W over TGA → `snv`; Q over TAA → `snv`; Y over TAG → `snv`; R over TAA →
      `multi_nt`).
    - **`U` over TGA and `O` over TAG** (M7-18 rule 4c, narrowed) are ALT `U`/`O`. Flags: `stop_codon`
      and `vcf_unrepresentable`. They get neither `snv` nor `multi_nt` (`nt_class` is not called;
      zero changes). They have no VCF record, because the minimal codon equals REF.
      `minimal_alt_codon` gains `U → TGA`, `O → TAG`, used only for that zero-change check. Every
      other U/O stays `rare`-skipped. B/J/X/Z stay `ambiguous`-skipped, including an `X` read from
      a patched database (example 9).
    - Blocks and swaps never carry REF `*` (rule 4).
    - **Error-rate class (§7.5.4).** A REF `*` site, and a U/O zero-change ALT, is in class
      `(substitution, mass_equivalent = false)`. `mass_equivalent` is false for REF `*` (`_mass("*")` is
      None, and `left_block_read` sees the `*` in its pair). Its `class_error_rate` is therefore
      the rate of that class, pooled over all its sites (the non-mass-equivalent
      substitutions). That rate is dominated by ordinary SNV sites and does not
      measure read-through reads, and with a reference search on Bakta proteins (which end at the
      stop) every read-through read is in the degenerate case of `variants.md` bias 1: the
      reference cannot contain the peptide. A second-pass search on a `makedb --patch` database
      does not lift this: its identification carries `X` at the stop, which differs from every
      de novo residue after I→L folding, so the pair still counts as an error. The rate feeds
      only `class_error_rate`/`CLASS_ERR`; `confidence` (`high` = `n_alt` ≥ k, M7-02)
      and therefore TE2 do not depend on it. 16b states this in `variants.md` (Q11: pooled,
      answered).
12. **Novel ORFs.**
    - M7-09 `class_rules`: `cterm_extension` is called when `spans_annotated_stop` is true and the
      in-frame continuation exists (M7-06: a next stop or `end_3p_open`). The label is
      `readthrough_candidate`. The M7-11 decision-2 rule then applies (`unestimated`, ≥ 2 distinct
      peptides). M7-09 decision 4's `ok = False` is removed.
    - `spans_annotated_stop` means: some peptide interval strictly contains the gene's terminal stop
      codon and extends past it on the 3′ side (along the origin on circular replicons).
    - M7-07 `patch_features` writes called `cterm_extension` rows: `CDS`, `ID` = `locus_tag` = the
      gene, `Note=cterm_extension`, gene start to the ORF's stop-codon end (`partial=true` when
      `end_3p_open`, as for new ORFs).
13. **`makedb --patch`** (`makedb.build.build_strain_db`). For a patch CDS with
    `Note=cterm_extension` and a `locus_tag` resolving to a CDS/sORF on the same strand, the residue
    at that gene's terminal stop codon translates as `X`. It must lie in frame inside the patch
    interval; otherwise raise `ValueError` naming the CDS. After that substitution, any remaining
    internal `*` raises exactly as today. The terminal `*` is stripped as today. Non-cterm patch
    rows are unchanged.
14. **Simulator** (gippasim, tooling). `SimConfig.n_readthrough: int = 0`. With the default 0, every
    existing scenario and seed is byte-identical. With `n > 0` the simulator picks `n` annotated
    genes (seeded through the generator's stream). Each must have its next in-frame stop 5–60 codons
    after the annotated stop, inside the replicon, wrapping on circular ones. The expressed protein
    reads the stop as `W` (TGA) or `Q` (TAA/TAG) and continues to the next stop. The annotation
    keeps the original stop. Truth goes to `truth/readthrough.tsv`: locus_tag, contig, strand,
    codon_start, codon, residue, ext_stop (0-based start of the new stop codon),
    `n_true_peptides_spanning`. In 16c, `n_selenoproteins: int = 0` plants an internal TGA read as
    `U` with `U` in the annotated protein; the full rule is rule 23.

## Examples

Toy contig `c1`, linear, frame `+0`. Residue `p` covers nt `[3p, 3p + 3)`:

```text
p      0         1         2         3         4
       012345678901234567890123456789012345678901
aa     MSDKPLTEAGRWQLNEFK*VPAGHTLDMRK*AAGLKEPSTR*
```

`G1` = CDS `+` `[0, 57)`: residues 0–18, and its terminal stop is the TGA at p = 18, which is
marked. The stop at p = 30 (TAA) and the one at p = 41 (TAG) are **not** annotated stops.
`max_diff`: 1 at 10–11 aa, 2 at 12–40 aa (§7.2.4). Unless stated otherwise the input is de novo
with `genome_relation = same`.

| # | Input | Expected |
|---|---|---|
| 1 | Query `WQLNEFKWVPAGHTL` (15 aa, k = 2) | one target row: `route` tolerant, `aa [11, 26)`, interval `c1:[33, 78) +`, `n_diff` 1, `diff_types` `substitution`; `locus_class` `cterm_extension`, gene G1 (crosses G1's 3′ end in frame). Variants: codon 54 REF `*` (TGA) ALT `W`, flags `snv;stop_codon` (`snv;strain_difference;stop_codon` under `related`); VCF POS 55, REF TGA, ALT TGG |
| 2 | **Sec peptide over fdnG's UGA**, Nissle, `search_db = same`: `SLGMLAVDNQARVUHGPTVASLAPTFGR` (28 aa); JMJJOG_02511, `−` `[2576167, 2579215)`; frame window `SLGMLAVDNQARV*HGPTVASLAPTFGR` over `[2578585, 2578669) −`, the `*` = internal TGA `[2578627, 2578630)` | **16a/16b (ADR-0058 items 1–5): no target row**: the TGA is not the gene's last codon (`[2576167, 2576170)`), so it is not marked; the peptide stays `no_hit` exactly as on main. **16c (item 6, approved; rules 21–26):** one row `[2578585, 2578669) −`, `n_diff` 1, `cds_inframe` JMJJOG_02511; variants: codon 2578627 ALT `U`, `candidate` (one distinct unique peptide), flags `stop_codon;vcf_unrepresentable`, no VCF record (would be POS 2578628) |
| 3 | **Non-annotated in-frame stop:** query `HTLDMRKQAAGLKE` (14 aa, k = 2) over `aa [23, 37)` = `HTLDMRK*AAGLKE` | no target row: p = 30 is not marked, so the `*` is never consumed; no other window within 2; `unmapped` `no_hit` |
| 4 | **Annotated stop + one more mismatch, over budget:** `NEFKWVPAGHA` (11 aa, k = 1) vs `aa [14, 25)` = `NEFK*VPAGHT` | no row (cost 2 > 1) |
| 4b | Same window, `NEFKWVPAGHT` (11 aa, k = 1) | one row `[42, 75) +`, `n_diff` 1 (the stop is the only allowed difference) |
| 4c | `NEFKWVPAGHAL` (12 aa, k = 2) vs `NEFK*VPAGHTL` | one row `[42, 78) +`, `n_diff` 2, `diff_types` `substitution;substitution`, stratum m = 2 |
| 5 | **Decoy symmetry:** `reverse` decoy of the 42-aa frame | mirrored cell 42 − 1 − 18 = 23 holds `*` and is marked; p = 30 → 11 is not. Query `LTHGAPVWKFENLQW` (the reverse of example 1): no target row; one decoy row at decoy `aa [16, 31)`, cost 1, mirrored target interval `[33, 78) +`. With the decoy cells removed (a deliberately broken build) that decoy row disappears: the test that pins symmetry. `shuffle` decoy: cell 18 marked, `*` at 18, 30, 41 in place |
| 6 | **`cterm_extension` past the annotated stop into the 3′ flank:** example 1's peptide plus `KWVPAGHTLDMR` (12 aa, `aa [17, 29)` = `K*VPAGHTLDMR`, cost 1) | both `cterm_extension` on G1. M7-06 ORF: G1's start to the next in-frame stop, `[0, 93) +` (stop TAA at p = 30). M7-09: `spans_annotated_stop` true, label `readthrough_candidate`, class FDR `unestimated`, 2 distinct peptides → `called` true. M7-07 `annotation_patch.gff3`: `c1 GIPPA CDS 1 93 . + 0 ID=G1;locus_tag=G1;Note=cterm_extension` |
| 7 | **makedb patch with `X`:** `gippa makedb --patch` on example 6's patch | record `>G1_patch …` with sequence `MSDKPLTEAGRWQLNEFKXVPAGHTLDMRK` (p = 18 → `X`, terminal `*` stripped); `db_id_map.tsv` row `G1_patch patch G1 … c1 1 93 +` |
| 7b | Patch row `Note=cterm_extension`, G1, `[0, 126) +` (runs over p = 30) | `ValueError`: "patch CDS G1_patch has an internal stop codon" (only the annotated stop becomes `X`) |
| 8 | **Minus strand:** contig `c2` linear, 300 nt; G2 CDS `−` `[120, 210)`, stop codon `[120, 123)` (gene-strand TAA) | site `(c2, −, 120)`; cell in frame `−0` (`ext_nt` 300, offset 0) at p = (300 − 0 − 3 − 120)/3 = 59; `to_interval(−0, 59, 60)` = `c2:[120, 123) −`. A read-through peptide's flank lies at coordinates < 120 |
| 9 | Second-pass peptide over the patched `X`: the engine reports `…NEFKXVPAG…` | placed as example 1 (any query letter over a marked cell costs 1); variants skip the X op (`ambiguous`, M7-18), so no ALT site; the ORF evidence (example 6) is unchanged |
| 10 | Nissle, ADR-0058 items 1–5, `sixframe_index` | event `mapping.stop_sites`: `terminal` 4,689, `not_a_stop` 3, `recoded` 0, `cells_target` 4,689 (all four replicons are circular in the recording; the 3 are partial genes at plasmid edges, JMJJOG_04967/04968/04975), equal decoy counts. **16c column:** example 29 (`recoded` 3, `cells_target` 4,692) |

## Edge cases

- **Minus strand.** Example 8. The terminal codon is the gene's lowest codon. The 3′ flank lies
  at lower coordinates. Variants read the codon reverse-complemented (M7-01 `translate_locus`).
- **Replicon end, linear.** A stop in the last codon of a contig has no 3′ flank. Only windows
  ending on it can consume it, and they classify as `cds_inframe` (Q3). A gene whose last codon
  runs off the contig is `not_a_stop` (rule 1).
- **Wrap, circular.** Length 300, `ext_nt` 426. A `+` gene `[240, 302)` has its stop at `[299, 302)`,
  which crosses the origin: `codon_start` 299, a single cell in frame `+2` at p = 99 (599 > 426).
  A `+` stop at `[0, 3)` has two cells, frame `+0` at p = 0 and p = 100, and the hits deduplicate
  modulo the length. A `−` stop at `[0, 3)` has two cells in frame `−0`, at p = 141 (x = 0) and
  p = 41 (x = 300). Every cell passes rule 2's round trip.
- **Overlapping genes sharing a stop.** Same strand, same frame (nested or alternative-start
  annotations, e.g. CDS `[0, 90)` and `[30, 90)` `+`): one site, `locus_tags` `G4;G5`, one cell.
  Genes in other frames or strands never share a codon. Their stops are separate cells, and a
  marked cell of gene A lying inside gene B in another frame is still only A's frame's cell.
- **sORFs.** Included like CDS.
- **Pseudogenes and frame-shifted CDS.** Pseudogenes are excluded entirely. A non-pseudo CDS with
  internal stops gets only its terminal codon marked. Its internal stops stay unconsumable
  (16c's recoded rule needs `U`/`O` in the annotated protein, never a `*`).
- **Alternative start annotation.** Irrelevant: only the stop is used. Two annotations differing
  only in their start share a site.
- **TAA / TAG / TGA.** The mapping rule ignores the codon. The variants rule distinguishes
  (rule 11): `U` is ALT only over TGA and `O` only over TAG. Standard residues get `snv`/`multi_nt`
  from the actual codon.
- **Peptides shorter than 10 aa** over a stop are never placed (exact route and `min_len`). This is
  documented in the methods note.
- **Two marked cells in one peptide** (adjacent short genes): each costs 1, so it is possible only
  at k = 2 with no other difference.
- **Peptide over a marked cell and an X codon:** the X is never consumed, so there is no row
  (unchanged).
- **Empty annotation** (no CDS/sORF): zero sites, byte-identical hits.
- **Old index file** (format 1) read by the new `map_tolerant`: `ValueError` asking for a rebuild
  (no silent fallback).

## Amendment: 16b reporting and the 16c labels (Khoa, 2026-10-24 evening)

Khoa's decisions (lead brief, 2026-10-24 evening; backlog "2026-10-24 evening (lead)" on #333).
They are written in here, not re-opened:

- **D1 (16b).** A `readthrough` flag on every `peptide_hits.tsv` row whose accepted alignment
  consumes an annotated terminal stop cell, a matching `mappings.gff3` attribute, the flag in the
  §9.4 vocabulary and in the schema. This answers Q6 (the earlier default was "no marker").
- **D2 (16b).** A run-level count of peptides by placement kind: exact; tolerant by diff type;
  frameshift; read-through; unmapped by reason. Only what no existing output gives is added.
- **D3 (16c, approved).** ADR-0058 item 6 is in force: annotated internal Sec/Pyl codons (`U`
  over TGA, `O` over TAG in the Bakta protein) are consumable, mirrored in the decoys. Rows that
  consume them get the flags `sec_site` and `pyl_site`. Only the label part is written here (rule
  19); the 16c matcher part is rule 1's `recoded` kind, written out in "16c amendment (2026-10-25)" (rules 21–26). 16c is next after 16b.
- **D4.** Q11 stays pooled: REF `*` sites (read-through) and U/O zero-change ALTs are in class
  `(substitution, mass_equivalent = false)`. Checked: rule 11, the 16b `variants.md` touch-point row, the
  16b acceptance criterion and Q11 already say so. Q11 is now marked answered.
- **D5 (statistician, #331 final review).** In `related_other_strain` at 1 Mb, 0–5 true
  read-through rows per replicate are accepted. 16b's e2e test must show that they become
  `stop_codon` ALT sites or `cterm_extension` evidence (rule 20).

### What already exists (checked on main c523c1c0)

| Output | Gives | Missing for D2 |
|---|---|---|
| `mapping/mapping_fdr.tsv` | per decoy type, stratum `n_diff` (and length class): `targets`, `decoys`, `n_accepted` | diff type, frameshift, read-through, unmapped. `n_diff` is not a diff type, and its counts are per stratum list, not distinct placed peptides |
| `run_manifest.json` › `mapping` | `search_db`, `n_queries` only | everything |
| V05 data TSV (`plots/data/V05_mapping_sankey.tsv`, M9-05) | distinct peptides placed exact / tolerant; unmapped as `contaminant` / `rejected by FDR` / `unexplained` | diff type, frameshift, read-through. `no_hit` and `too_short` are merged into `unexplained`. It is a plot twin, written only when plots run, with no table schema |
| `tables/detection_summary.tsv` QC rows (M4-16, `context.summary.qc_rows`) | `contaminant_peptides`, `unmapped_peptides` (by label) | everything by kind and by reason |
| event `mapping.stop_sites` | cell counts, not peptides | read-through peptides |

So nothing gives the D2 breakdown. It goes into **existing** `detection_summary.tsv` rows (rule
17): section `qc`, new `group_type` values, no new column and no new file. The exact count and the
unmapped totals overlap V05 and the existing QC rows on purpose, so that the breakdown sums within
one table. Acceptance checks that they agree.

### Rules (16b, 16c labels)

15. **Which rows consumed a cell** (`readthrough.consumed_sites`, new, pure). A `peptide_hits` row
    consumes the site `S` when all of these hold:
    - the row is a target row with no `blocks` (frameshift rows never consume a cell, rule 7;
      exact rows never contain a `*`);
    - `S.contig == contig` and `S.strand == strand`;
    - with `n` the replicon length and the row's unwrapped 0-based interval `[start, end)`, some
      `x ∈ {S.codon_start, S.codon_start + n}` (the second only on circular replicons) has
      `start ≤ x` and `x + 3 ≤ end` and `(x − start) mod 3 == 0`.

    For a non-frameshift row this is exactly "the window holds a cell of `S`". The tolerant DP
    can cross a marked cell only by consuming it (rule 4), so holding a cell and consuming it are
    the same thing. This equals the test helper `tests/helpers/readthrough_counts.holds_cell` on
    the row's `frame_id`/`aa_start`/`aa_end`, which the acceptance criteria pin as a property. The
    function returns, per row, the consumed sites sorted by `codon_start`, and an empty tuple for
    the rest. It reads the sites the index was searched with (`load_stop_cells(six).sites`), never
    a fresh `annotated_stops` call, so the flag can never disagree with the search (16c changes
    the site set). Signature:

    ```python
    def consumed_sites(rows: pd.DataFrame, sites: Sequence[StopSite],
                       lengths: Mapping[str, int], circular: Mapping[str, bool]
                       ) -> list[tuple[StopSite, ...]]
    ```

16. **The `readthrough` flag and the GFF3 attribute** (`steps.mapping.loci`,
    `mapping.loci.annotate_peptide_hits`, `mapping.outputs`).
    - `loci` gains the **optional** input `sixframe`; `mapping.nf` always wires it
      (`stepIn('loci', [g, pep, six, hx, ht, …])`). With it, `loci` calls
      `load_stop_cells(six)`, so a format-1 file raises as in `map_tolerant`. It then runs
      `consumed_sites` on the `_peptide_rows` frame, which still carries `blocks`, and passes the
      result to `annotate_peptide_hits(…, consumed=…)`.
    - A row with at least one consumed site of kind `terminal` gets the flag `readthrough`. Flags
      stay `;`-joined in alphabetical order: `readthrough;strain_difference` under `related`,
      where every tolerant row is `strain_difference`.
    - Without the `sixframe` input (old work directories, unit callers) no row is flagged, and
      `loci` logs one WARNING, event `mapping.readthrough_unchecked`. The flag is never guessed
      from geometry against a gene model.
    - `peptide_hits.tsv` gains **no column**. The consumed codon positions travel in a work column
      `readthrough_codon` that `write_mapping_outputs` does not write to the TSV.
    - `mappings.gff3`: a flagged peptide feature carries the attribute
      `readthrough_codon=<p1>[,<p2>]`, the **1-based** lowest coordinate of each consumed terminal
      codon (`codon_start + 1`, the VCF POS convention of example 1), ascending, `,`-separated (the
      GFF3 multi-value separator). Its `flags` attribute already carries `readthrough`. Frameshift
      `match_part` children never carry it. The attribute comes after `flags`; the order of the
      existing attributes is unchanged.
    - BED, proBed and `unmapped.tsv` are unchanged.
17. **Placement-kind counts** (`context.summary.placement_kind_rows`, new, pure; added as one more
    part in `steps.context_stats._context_stats`, next to `qc_rows`, which is **not** changed).
    The rows are section `qc`, metric `n_peptides` (the count is in `value`), `ci_low` and
    `ci_high` NA, and empty `flags`. Units are distinct peptide sequences: `peptide_hits.peptide`
    and `unmapped.sequence` for `kind = peptide`. The groups, in this order:

    | `group_type` | `group` | Counted | `n` |
    |---|---|---|---|
    | `placement_kind` | `exact`, `tolerant`, `frameshift` | placed peptides, **exclusive**: `exact` if any row has `route = exact`; else `tolerant` if any row is a non-frameshift tolerant row; else `frameshift` (every row has `frameshift` in `diff_types`) | placed peptides |
    | `tolerant_diff_type` | `substitution`, `swap`, `block`, `near_block` | `tolerant`-kind peptides with at least one non-frameshift row whose `diff_types` contains the type. **Not exclusive**: `substitution;swap` counts in both | `tolerant`-kind peptides |
    | `placement_flag` | `readthrough` (16c adds `sec_site`, `pyl_site`) | placed peptides with at least one row carrying the flag | placed peptides |
    | `unmapped_reason` | `no_hit`, `rejected_placement_fdr`, `too_short`, `contaminant_id` | `unmapped.tsv` rows with `kind = peptide` by `reason`, both labels | peptide rows of `unmapped.tsv` |

    - Every group is written, zeros included, so the row set does not depend on the data.
      `below_identity` is protein-only (`contaminants.unmapped_table`) and is not a peptide group.
    - With an empty `peptide_hits` the placement rows are 0 with `n` 0. Without `unmapped.tsv`
      (`unmapped` None) the `unmapped_reason` rows are NA, as `qc_rows`' C13 counts are.
    - These are counts, not rates or FDRs. No row divides read-through by anything, and none
      carries D_rt/T_rt (methods note item 5, R3). `n` is the population the count comes from,
      not a denominator to report: no metric name ends in `_rate`, `_share` or `_fdr`, and no
      writer or report presents `value / n` of a `placement_flag` or `tolerant_diff_type` row
      as an error rate.
    - The flagged rows are a subset of a list accepted at q. The list FDR does not control that
      subset (methods note items 5–6): the 16b paragraph of `placement-fdr.md` says so, and no
      output describes the `readthrough` rows (or, in 16c, the `sec_site`/`pyl_site` rows) as
      "accepted at q" on their own.
18. **Vocabulary and schema.** §9.4 *Placement* gains `readthrough` (16b), and `sec_site` and
    `pyl_site` (16c). §9.4's own heading ("extend only through the schema") allows this. The
    `flags` field of `src/gippa/schemas/tables/peptide_hits.schema.json` gains a `pattern`
    constraint in the form `variants.schema.json` uses (`(?:a|b|…)(?:;(?:a|b|…))*`). It lists every
    value the code writes today plus the new ones: `engine_protein_mismatch`, `low_complexity`,
    `multi`, `multi_identical`, `nterm_extension_evidence`, `origin_spanning`, `readthrough`,
    `same_db`, `strain_difference` (16c adds `pyl_site`, `sec_site`). Its description names
    `readthrough`. `multi` and `low_complexity` are written today but are not in §9.4 (Q16).
19. **16c labels** (approved; PR 16c, after 16b). A row consuming a site of kind `recoded` gets
    `sec_site` when the site's codon is TGA (Bakta `U`) and `pyl_site` when it is TAG (Bakta `O`),
    **instead of** `readthrough`. A row consuming a terminal and a recoded site gets both flags.
    The GFF3 attribute for recoded sites is `recoded_codon=<p1>[,<p2>]`, with the same coordinate
    rule as `readthrough_codon`. The `placement_flag` groups `sec_site` and `pyl_site` are added
    to rule 17. `StopSite` already carries `kind` and `codon`, so 16c changes only the mapping
    from kind to flag in `annotate_peptide_hits`, the GFF3 writer and the summary groups. Decoys
    mirror the recoded cells exactly as the terminal ones (rule 3). The TE1 exception label is
    `sec_site` (rule 24). The flag depends on the consumed site's `kind` and `codon`, never on the
    query letter: a `W` read over fdnG's TGA still gets `sec_site` in `peptide_hits`, and
    `variants` calls it an ordinary `W` ALT with flags `snv;stop_codon` (rule 26, example 30).
20. **Strain stop-loss in gippasim (D5).**
    - **Finding (the strain model, checked).** `gippasim.strain.derive_reference` puts SNVs at
      uniform positions, stop codons included. The mapping genome is the sample. The reference
      carries the SNVs, and `_lift` translates each reference feature over the lifted interval
      only (`protein(g.nt)`, which stops at the first stop or at the interval end). When an SNV
      turns the sample's terminal stop into a sense codon in the reference, the reference protein
      is the sample protein plus **one** residue at the old stop. The feature is then flagged
      pseudo there, but the `other_strain` wrong-peptide pool (`proteomics`, `f.kind in
      ("cds", "sorf") and f.aa`) still digests it. Its C-terminal peptide ends on the sample's
      marked cell. That explains the 0–5 "true read-through" rows per replicate at 1 Mb: they are
      placed at their origin locus, with one substitution at the stop and the window **ending on
      the stop codon**.
    - **Consequence.** Under Q3 they classify as `cds_inframe`, not `cterm_extension`, and
      rule 12 gives them no `spans_annotated_stop` (nothing extends past the stop). So they must
      become **`stop_codon` ALT sites**: REF `*`, ALT the reference residue, flags
      `snv;strain_difference;stop_codon` (every single-base sense change of a stop is one nt away;
      `nt_class` comes from the codon), confidence `candidate` with one distinct peptide. The
      gippasim model never makes the other-strain peptides extend past the stop, so the
      `cterm_extension` path is covered only by `n_readthrough` (rule 14).
    - **Simulator knob** (tooling, implementer). `StrainConfig.n_stop_loss: int = 0`. With 0,
      every existing strain build is byte-identical and no new file is written. With `k > 0`,
      `derive_reference` additionally places `k` SNVs, drawn from a separate generator (for
      example `default_rng([seed, 2, seed_offset, 1])`) so the existing draws are unchanged. Each
      SNV lies in the terminal stop codon of a distinct non-pseudo sample CDS, turns it into a
      sense codon in the reference, and keeps `_SPACING` from the other events. It is also
      written as an ordinary `snv` row of `strain_differences.tsv`. A gene is eligible only if
      the reference protein then has a C-terminal digest product (`digest_len`,
      `missed_cleavages`) that contains the new residue, has **at least 10 aa** (shorter
      peptides over a stop are never placed) and holds no other strain event in its codons
      (so it places at `n_diff` 1). Fewer than `k` eligible genes, or fewer wrong-peptide slots
      than `k`, raise `ValueError` (no silent shortfall). Truth: `truth/stop_loss.tsv` with
      `locus_tag`, `contig`, `strand`, `codon_start` (sample, 0-based), `sample_codon`,
      `reference_codon`, `reference_residue`, `cterm_peptide` (the reference protein's C-terminal
      digest product that contains that residue) and `planted` (true for the `k` knob SNVs). It
      lists **every** stop-loss in the build: the `k` planted ones and any that the uniform SNV
      draws produced (`planted` false, `cterm_peptide` empty when none qualifies). With
      `wrong_source = "other_strain"`, every planted `cterm_peptide` that is absent from the
      sample (I = L) joins the wrong-peptide set first; the remaining slots are filled as
      today.

### Examples (16b, 16c labels)

Toy contig `c1` and gene G1 as in the examples above (G1's stop: TGA, codon `[54, 57)`, cell
p = 18).

| # | Input | Expected |
|---|---|---|
| 11 | Example 1's accepted row (`WQLNEFKWVPAGHTL`, `[33, 78) +`, `cterm_extension`) | `peptide_hits.tsv` `flags` `readthrough` (`readthrough;strain_difference` under `related`); no new column. `mappings.gff3` peptide feature `… ;flags=readthrough;readthrough_codon=55` |
| 12 | **Window ending on the stop (Q3):** `GRWQLNEFKW` (10 aa, k = 1) vs `aa [9, 19)` = `GRWQLNEFK*` | one row `[27, 57) +` (TSV `start` 28, `end` 57), `n_diff` 1, `cds_inframe` G1, flag `readthrough`, `readthrough_codon=55`; variants: codon 54 ALT `W`, `snv;stop_codon`; no `cterm_extension` evidence (`spans_annotated_stop` false) |
| 13 | Example 3's query (over the unannotated stop at p = 30) and any exact row | not placed / no flag; `consumed_sites` returns `()` |
| 14 | **Minus strand:** example 8's site `(c2, −, 120)`; an accepted tolerant row `c2:[111, 141) −` | `(120 − 111) mod 3 = 0` and `[120, 123) ⊂ [111, 141)`: flag `readthrough`, `readthrough_codon=121`. The same row at `[112, 142) −` (another frame): no flag |
| 15 | **Origin copy:** circular length 300, `+` site at `codon_start` 0; a row with unwrapped `[282, 312) +` | `x = 0 + 300` lies in it with `(300 − 282) mod 3 = 0`: flag, `readthrough_codon=1`, plus `origin_spanning` (as today) |
| 16 | **Frameshift row** whose span covers a site codon (blocks set, `diff_types` `frameshift`) | never flagged (rule 15), even if the codon lies inside its unwrapped span in frame |
| 17 | **Placement-kind rows.** Placed: 3 exact peptides; tolerant `A` (`substitution`, consumes G1's cell); `B` (`substitution;swap`); `C` (`block`); `D` with two loci, one row `substitution`, one row `frameshift;substitution`; `E` with only a frameshift row. Unmapped peptides: `no_hit` 2, `too_short` 1, `rejected_placement_fdr` 1, `contaminant_id` 1 (label `contaminant`) | `placement_kind`: exact 3, tolerant 4 (A, B, C, D), frameshift 1 (E), each with `n` 8. `tolerant_diff_type`: substitution 3 (A, B, D), swap 1, block 1, near_block 0, each with `n` 4. `placement_flag` readthrough 1 (`n` 8). `unmapped_reason`: no_hit 2, rejected_placement_fdr 1, too_short 1, contaminant_id 1 (`n` 5). The existing `unmapped_peptides` = 4 and `contaminant_peptides` = 1 are unchanged |
| 18 | Empty run (no peptides): `peptide_hits.tsv` header only, `unmapped.tsv` header only | every rule-17 row is written with `value` 0 and `n` 0 |
| 19 | **Strain stop-loss (rule 20):** sample gene with stop TAA at `[c, c + 3)`; the reference SNV makes it TAC (`Y`); the wrong peptide is the reference's C-terminal peptide `…LVKAEY` | one tolerant row ending at `c + 3` (+ strand), `n_diff` 1, `cds_inframe`, flags `readthrough;strain_difference`; variants: REF `*` (TAA) ALT `Y`, flags `snv;strain_difference;stop_codon`, `candidate`; no novel ORF row from it |
| 20 | **16c:** example 2's fdnG row (`SLGMLAVDNQARVUHGPTVASLAPTFGR`, `[2578585, 2578669) −`, recoded TGA at `[2578627, 2578630)`) | flag `sec_site` (no `readthrough`); GFF3 `recoded_codon=2578628`; `placement_flag` `sec_site` counts its peptide. An `O`-over-TAG site gives `pyl_site` |

### Edge cases (16b, 16c labels)

- **Multi-locus peptide** (`multi_identical` or `shared`): the flag is per row. Only rows that
  consume a cell carry it, and the peptide counts once in `placement_flag` if any row does.
- **Two cells in one row** (adjacent short genes, k = 2): one flag; `readthrough_codon` lists both
  positions, ascending.
- **Ties between kinds:** one peptide's best-cost rows can mix a frameshift and a non-frameshift
  locus at equal cost. `placement_kind` takes the first match in the order exact, tolerant,
  frameshift (example 17, `D`). Exact and tolerant never mix, since tolerant rows exist only
  for exact leftovers.
- **Linear replicon:** one candidate `x` per site. A site in the last codon of a linear contig
  can only be consumed by a window ending on it (example 12's shape).
- **Circular, origin-crossing codon** (`codon_start` = n − 1 or n − 2): `x = codon_start` with
  `x + 3 > n`; the row's unwrapped interval must reach `x + 3`. Positions in the GFF3 attribute
  use the wrapped `codon_start + 1`.
- **Minus strand:** containment and `(x − start) mod 3` work unchanged because a
  non-frameshift row's length is a multiple of 3 (example 14).
- **Decoy rows** never reach `peptide_hits.tsv`, so they are never flagged. The D_rt counts stay
  in the calibration tests only.
- **Empty annotation:** zero sites, no flag, `readthrough` count 0.
- **No `sixframe` input to `loci`:** no flag, one WARNING event (rule 16). The summary then
  counts 0 read-through peptides; the WARNING is the only signal (Q14).
- **`same_db` exact skipping** never touches read-through rows (they are tolerant, m ≥ 1).
- **16c:** **per feature**, a site is either `terminal` or `recoded`, never both (an internal codon
  is never the last codon). Across features, a gene's terminal codon can be another gene's
  recoded codon; they merge into one `recoded` site with one cell (rule 21 (f), example 27), so
  the flags still never conflict.

### ADR and spec deviations

- **ADR-0058** needs an amendment by the lead, recording Khoa's evening answers:
  - status: item 6 accepted (D3);
  - Consequences › Interfaces says "no new column, no new diff type, **no new flag value**". The
    flags `readthrough`, `sec_site` and `pyl_site`, the GFF3 attributes `readthrough_codon` and
    `recoded_codon`, the optional `loci` input and the rule-17 summary rows replace that sentence;
  - Consequences › TE1: under item 6 the 4 misses become `sec_site` exceptions.
- **GIPPA_SPEC.md, no ADR needed in my reading** (Q15):
  - §9.4 Placement gains three flags, which §9.4's heading allows "through the schema";
  - §9.3 `peptide_hits.tsv` columns are unchanged;
  - §9.3 `detection_summary.tsv` gains rows of the existing `qc` section with no new column. The
    QC list in §9.3 ("id_mixups, genome_inconsistent_rate, contaminant counts") is a summary, not a
    closed list; M4-16 already added `engine_protein_mismatch_rate` and the ADR-0019 rows;
  - §7.2.9 / M3-16's GFF3 attribute list gains `readthrough_codon` (and `recoded_codon`), as M3-16
    lead decision 2 added `frameshift_reread`.
- No tolerance, bound, α, R or gate changes (R3). No protected file is touched.

## 16c amendment (2026-10-25)

This section makes 16c buildable. It writes out the items the statistician's pre-code check
(`local/stat/m7-16c-quick-check.md`, verdict GO WITH CHANGES, item Q5) found missing. It also
replaces the 16c rows, the 16c acceptance items and attack (d) in place. Inputs: ADR-0058 item 6,
`docs/methods/placement-fdr.md` "16c: annotated recoded codons". The method does not change: no
tolerance, bound, α or R changes (R3), and no protected file is touched. Where this section and an
older 16c sentence disagree, this section wins.

### Rules (16c)

21. **Rule 1, recoded part** (`annotated_stops(gm, genome, recoded=True)`, pure).
    - **(a) Which features.** Every gene with `feature_type` in {`cds`, `sorf`} whose terminal
      codon is a stop (the features rule 1 turns into terminal sites). Pseudogenes
      (`feature_type = pseudogene`) are never scanned. Features counted under `not_a_stop` (for
      example 3′-partial genes) are not scanned either (Q24).
    - **(b) Length check.** For the unwrapped interval `[s, e)`, the protein must satisfy
      `len(gene.aa) == (e − s) / 3 − 1`, and `(e − s)` must be divisible by 3. Bakta's `aa` has no
      trailing `*`. A feature that fails is not scanned. It still gets its terminal site, and the
      failure is counted once per feature under `length_mismatch`.
      A feature with an empty `aa` (no protein in the annotation) is not scanned and not
      counted under either skip (lead decision, #346 review).
    - **(c) Residue → codon.** For a 0-based residue index `i` with `1 ≤ i ≤ len(aa) − 1`:
      - on `+` the codon is `[s + 3i, s + 3i + 3)`;
      - on `−` the codon is `[e − 3i − 3, e − 3i)`;
      - `codon_start` is the lowest coordinate, wrapped mod the replicon length `n` on circular
        replicons (the same convention as terminal sites).

      The codon is read on the gene strand with rule 1's reader (reverse complement on `−`).
    - **(d) Which residues.** Only an upper-case `U` or `O` at `i ≥ 1` is a candidate.
      - Residue 0 is the start codon and is never scanned: a `U`/`O` there is ignored and not
        counted.
      - The terminal codon is never part of `aa`, so it is never recoded.
      - `*`, `X` and every other letter inside `aa` are ignored (no site, no count).
    - **(e) Codon check.** A `U` must sit over `TGA` and an `O` over `TAG`. **Default: skip and
      count** (Q5 1b). Any other codon gives no site, and the event counts it under
      `codon_mismatch`, once per residue. It is annotation data, not a bug, so it never raises.
    - **(f) Merging.** Sites are keyed by `(contig, strand, codon_start)` across all features,
      terminal and recoded together. Features that share a key form **one site with one cell**,
      with `locus_tags` the sorted union. Its `codon` is the codon read there (identical for every
      feature, since the key fixes the strand and the bases). Its `kind` is:
      - `recoded` when at least one feature recodes that codon. This covers the nested-gene case,
        where one gene's terminal stop is another gene's internal `U`-over-TGA (Q22). The flag
        follows the kind (rule 19), so such a row gets `sec_site` and never `readthrough`;
      - `terminal` otherwise.

      The 16a edge case "a site is either `terminal` or `recoded`, never both" holds **per
      feature**. Across features the merge applies.
    - **(g) Order.** Recoded and terminal sites are sorted together by (contig order, strand,
      `codon_start`), as in rule 1. A protein with several `U`s gives one site per `U`.
    - **(h) `recoded=False`** (the default) gives exactly the 16a/16b output. No `aa` is read.

22. **Signature and callers.**
    ```python
    def annotated_stops(gm: GeneModel, genome: Mapping[str, str], recoded: bool = False
                        ) -> tuple[tuple[StopSite, ...], int]       # unchanged return type
    def recoded_skips(gm: GeneModel, genome: Mapping[str, str]) -> dict[str, int]
        # {"length_mismatch": n, "codon_mismatch": n}; pure; the rule-21 (b)/(e) counts
    ```
    - The return type is unchanged, so every 16a/16b caller and test stays valid.
    - `sixframe_index` and `variants_step` both call `annotated_stops(gm, genome, recoded=True)`.
      They must build the same site set; attack (d4) pins that.
    - The `mapping.stop_sites` event (rule 8) has these values:
      - `recoded` is the number of `kind == "recoded"` sites;
      - `terminal` is the number of `kind == "terminal"` sites, so a merged site is counted once,
        as recoded;
      - `cells_target` and each `cells_decoy[kind]` include the recoded sites' cells;
      - a new key `recoded_skipped` holds `recoded_skips(...)`, for example
        `{"length_mismatch": 0, "codon_mismatch": 0}` (Q21).

23. **Simulator `n_selenoproteins`** (gippasim tooling; replaces rule 14's last sentence).
    - **Config.** `SimConfig.n_selenoproteins: int = 0`.
    - **Default stays byte-identical.** With 0, nothing is drawn and nothing is edited.
      `truth/selenoproteins.tsv` is not written, so every existing scenario and seed is
      byte-identical, including the committed gippasim goldens.
    - **Generator.** All draws come from a separate generator,
      `np.random.default_rng([seed, _SEC_STREAM])`, with a new SeedSequence key that no other
      gippasim module uses. Keys in use on 267d0130 are 1–9 and 40, so 12 is suggested. No other
      stream moves.
    - **Gene choice.** `n` distinct genes are drawn uniformly without replacement from the
      eligible genes. A gene is eligible if it:
      - is an annotated `cds` (not sORF, not pseudogene);
      - does not cross the origin;
      - has `len(aa) ≥ 21`;
      - has at least one eligible codon (next item).
    - **Codon choice.** One codon per gene, drawn uniformly from the gene's eligible residues `i`
      with `10 ≤ i ≤ len(aa) − 11`. That is at least 10 codons from the start codon and at least
      10 sense codons before the terminal stop, so a tryptic peptide of ≥ 10 aa can span it. The
      codon must not overlap any other annotated feature (any type, either strand) or any other
      simulator edit (stop-loss, N-terminal extension, novel ORF).
    - **Edit.** The codon is replaced by `TGA` on the gene strand (`TCA` on `+` coordinates for a
      `−` gene). Only that codon changes, so no other in-frame stop appears in the gene.
    - **Annotation and expression.** Residue `i` of the gene's `aa` becomes `U` wherever gippasim
      writes the annotated protein: the gene table, the Bakta JSON and `genes.parquet`. The
      expressed protein carries `U` at `i` as well, so true peptides spanning it contain `U` and
      digestion treats it as an ordinary residue. The implementer checks that gippasim's peptide
      mass and digestion code accept `U`; gippa's `residue_masses.tsv` already has
      `U 150.95364`.
    - **Truth.** `truth/selenoproteins.tsv`, one row per planted site, sorted by (contig,
      `codon_start`). Columns: `locus_tag`, `contig`, `strand`, `codon_start` (0-based, lowest
      coordinate), `residue` (0-based `i`).
    - **Shortfall.** If fewer than `n` genes are eligible, raise
      `ValueError("n_selenoproteins=<n> but only <k> eligible genes")` before any edit.
    - **Strain.** `n_selenoproteins > 0` together with `strain` set raises `ValueError` (not
      supported: the strain SNVs could hit the planted codon). Attack (d1) and the 16c
      integration test are same-strain.
    - **Pyl. Default: Sec only. No `O`/TAG is planted** (Q23). `pyl_site` is therefore reached
      only by unit tests. That is a stated limitation of 16c's evidence: no genome in the
      recordings has a Pyl site either.

24. **TE1 `sec_site` label** (`validation/proxies/te_proxies.py`, `te1_result`). **Default:** the
    exception row keeps `locus_tag`, `peptide` and `hits` as today, and its `kind` is:
    - `"sec_site"` for a miss whose peptide has at least one `peptide_hits` row and every one of
      its rows carries `sec_site` in `flags` (`peptide_hits.tsv` holds accepted placements only);
    - `"miss"` for every other miss. A peptide with no row keeps `hits = "unplaced"`.

    The other parts of the result:
    - `counts` gains the diagnostic `sec_site`: the number of misses labelled `sec_site`, counted
      before the exception cap;
    - `misses`, `exact_inframe`, the share and the pass rule (`share ≥ 0.999`) are unchanged.
      Tolerant rows are never exact, so the 4 fdnG peptides stay misses;
    - `te1_result` reads the `flags` column when it is present. A frame without it is treated as
      having no flags (Q25), so callers that omit it keep today's result;
    - **q-values.** `te1_result` gains `per_peptide: pd.DataFrame | None = None`, the run's
      per-peptide placement table (the one TE8 reads). When it is given, every `miss` or
      `sec_site` exception row whose peptide appears there also gets `placement_q`, `m_target`
      and `accepted`. `te1_proxy` passes it. With `None` the rows are as today (Q26).

25. **TE2 note and count** (`te2_result`).
    - `TE2_NO_STOP_NOTE` becomes exactly:
      `"no high stop_codon site: annotated stops and recoded Sec/Pyl codons are consumable
      (M7-16, ADR-0058 items 1-6), but no stop_codon ALT site reached high confidence (on Nissle
      the two fdnG Sec sites are candidate, one unique peptide each; fdhF has no tryptic peptide
      over its Sec codon)"`
      It is one string, wrapped here only for width. The module docstring's TE2 sentence ("…
      stay unmapped; ADR-0058") is changed to match.
    - `counts` gains `alt_stop_codon`: the ALT rows flagged `stop_codon`, at any confidence (Q27).
      The pass rule is unchanged.

26. **Rule 19, added sentence.** The flag depends on the consumed site's `kind` and `codon`, never
    on the query letter. A `W` read over fdnG's TGA still gets `sec_site` in `peptide_hits`, and
    `variants` calls it an ordinary `W` ALT with flags `snv;stop_codon`.

### Examples (16c)

Toy genomes: one contig `c`. Only the named codons matter, and every other base is filler that
makes no stop on the gene's strand and frame. Sites are written
`(contig, strand, codon_start, locus_tags, codon, kind)`.

| # | Input | Expected output |
|---|---|---|
| 21 | **`+`, linear.** `c` linear, length 60. G1 `cds` `+` `[3, 33)`, bases `[3, 33)` = `ATG GCT GCT GCT TGA GCT GCT GCT GCT TAA`, `aa = "MAAAUAAAA"` (9 = 30/3 − 1). | `recoded=True`: `((c,+,15,(G1,),TGA,recoded), (c,+,30,(G1,),TAA,terminal))`, `n_not` 0; `recoded_skips` `{length_mismatch 0, codon_mismatch 0}`. `recoded=False`: `((c,+,30,(G1,),TAA,terminal),)`. Residue 4 → `[3 + 12, 3 + 15)` |
| 22 | **`−`, linear.** `c` length 60. G2 `cds` `−` `[3, 33)`, `aa = "MAAAUAAAA"`. On `+` coordinates: `[30, 33)` = `CAT` (ATG), `[18, 21)` = `TCA` (TGA), `[3, 6)` = `TTA` (TAA). | `((c,−,3,(G2,),TAA,terminal), (c,−,18,(G2,),TGA,recoded))`. Residue 4 → `[33 − 15, 33 − 12)` |
| 23 | **`+`, circular wrap.** `c` circular, length 300 (`ext_nt` 426, the 16a edge case). G3 `cds` `+` `[270, 330)` (crosses the origin), 19 aa, `U` at residue 12, `TGA` at `[306, 309)` ≡ `[6, 9)`, terminal at `[327, 330)` ≡ `[27, 30)`. | `((c,+,6,(G3,),TGA,recoded), (c,+,27,(G3,),TAA,terminal))`. `frame_cells`: the recoded site has two cells in frame `+0`, p = 2 (x = 6) and p = 102 (x = 306), and both pass the round trip |
| 24 | **`−`, circular wrap.** Same contig. G4 `cds` `−` `[270, 330)`, 19 aa, `U` at residue 2. | Codon `[330 − 9, 330 − 6)` = `[321, 324)` ≡ 21. Sites `((c,−,21,(G4,),TGA,recoded), (c,−,270,(G4,),TAA,terminal))` (terminal `[270, 273)` reads TAA on `−`). Two cells in frame `−0`: p = 134 (x = 21) and p = 34 (x = 321) |
| 25 | **Codon mismatch.** Example 21, but `[15, 18)` = `TAA` with `aa` still `U` at 4; and a second gene with `O` over `TGA`. | No recoded site for either. Terminal sites unchanged. `recoded_skips` `codon_mismatch` 2. No exception |
| 26 | **Length mismatch.** Example 21 with `aa = "MAAAUAAAA*"` (10 residues). | Only `(c,+,30,(G1,),TAA,terminal)`; `length_mismatch` 1, `codon_mismatch` 0 |
| 27 | **Nested merge.** Example 21's G1, plus G5 `cds` `+` `[3, 18)`, `aa = "MAAA"`, whose terminal codon `[15, 18)` is G1's recoded TGA. | `recoded=True`: `((c,+,15,(G1,G5),TGA,recoded), (c,+,30,(G1,),TAA,terminal))`: one site, one cell. Event: `terminal` 1, `recoded` 1, `cells_target` 2. A row over `[15, 18)` gets `sec_site`, not `readthrough`. `recoded=False`: `((c,+,15,(G5,),TGA,terminal), (c,+,30,(G1,),TAA,terminal))` |
| 28 | **Pyl, excluded and ignored letters.** G6 `+` with `O` at residue 3 over `TAG`; G7 a `pseudogene` with `U` over TGA; G8 with `U` at residue 0; G9 with `X` and `*` inside `aa`. | One recoded site for G6 with codon `TAG` (rule 19: `pyl_site`). Nothing for G7, G8 or G9, and no skip counted for them |
| 29 | **Nissle recording** (`1.12.2-full/nissle`), `sixframe_index` (example 10's 16c column). | `mapping.stop_sites`: `terminal` 4,689, `recoded` 3, `not_a_stop` 3, `cells_target` 4,692, `cells_decoy` 4,692 for each decoy type, `recoded_skipped` `{0, 0}`. The recoded sites are `(contig of JMJJOG_02511, −, 2578627, (JMJJOG_02511,), TGA, recoded)`, JMJJOG_04733 `+` 4825395 and JMJJOG_04490 `+` 4576497 (residues 195, 195, 139) |
| 30 | **Rule 26.** Nissle `same` search with a strain-like peptide `SLGMLAVDNQARVWHGPTVASLAPTFGR` over JMJJOG_02511's recoded TGA. | `peptide_hits` row tolerant, `n_diff` 1, flags `sec_site` (no `readthrough`), GFF3 `recoded_codon=2578628`. `variants`: ALT `W` at codon 2578627, flags `snv;stop_codon`, a VCF record (TGA → TGG on the gene strand) |
| 31 | **TE1 label** (`te1_result`, unit). Sources `{pS: {G1}, pT: {G1}, pU: {G1}}`, no exclusion. `peptide_hits`: pS two tolerant rows, flags `sec_site` and `sec_site;multi`; pT one tolerant row, flags `""`; pU no row. | `misses` 3, `sec_site` 1. Exceptions: `{kind: sec_site, peptide: pS, hits: "tolerant:…;tolerant:…"}`, `{kind: miss, peptide: pT}`, `{kind: miss, peptide: pU, hits: "unplaced"}`. The share and `passed` are as without the flags column |
| 32 | **gippasim** sim100k, same strain, `n_selenoproteins = 5`, seed 1. | `truth/selenoproteins.tsv` has 5 rows. Each genome codon there reads TGA on the gene strand, and the gene's `aa` has `U` at `residue`. The event reports `recoded` 5. `n_selenoproteins = 0` gives an output tree byte-identical to the config without the field. `n_selenoproteins` above the eligible count raises `ValueError` |

### Edge cases (16c)

- **Minus strand.** Rule 21 (c), examples 22 and 24. The residue index counts from the gene's
  5′ end, which is the highest coordinate.
- **Circular and linear.** On a circular replicon the codon wraps mod `n`, and a codon near the
  origin can have two cells (example 23). On a linear replicon the codon always lies inside the
  gene's interval, so it never runs off the contig.
- **A codon that crosses the origin** (`[n − 1, n + 2)`): `codon_start = n − 1`, read with rule
  1's wrapped reader. `frame_cells` treats it like a terminal codon at the origin (16a edge case).
- **Empty input.** An empty gene model gives `((), 0)` for any `recoded` value, and
  `recoded_skips` gives zeros. A gene model whose proteins have no `U`/`O` gives exactly the
  `recoded=False` output, which is the 16a/16b toys and gippasim by default.
- **Ties.** Two features recoding the same codon (alternative starts): one recoded site, tags
  merged. A terminal and a recoded site at the same key: one `recoded` site (example 27). Sites
  never tie in the sort after merging, because the key is unique.
- **Several `U`s in one protein:** one site each, all sorted with the rest.
- **Adjacent recoded cells in one peptide** cost 1 each, so they are possible only at k = 2 with
  no other difference (as for terminal cells).
- **`O` over TAG on `−`:** `+` coordinates read `CTA`.
- **Peptides spanning fdhF's Sec codon:** none of 10–30 aa, since there is no K/R within 30 aa
  (statistician Q3). So fdhF has no TE1/TE2 row. This is expected, not a miss.
- **Stale index.** An index built by 16b code is format 2 but has no recoded cells. Every run
  rebuilds `sixframe_index`, so the default is no format bump (Q28).

### TE2 and TE1 expectations on Nissle (16c)

- **TE2** (same strain, `te2_proxy`):
  - `alt_sites` 233 → **235**, `high` **80**, `high_start_codon` **80**, `high_stop_codon`
    **0**, `violations` **0**, `alt_stop_codon` **2**;
  - the two new sites are ALT `U` at JMJJOG_02511 `[2578627, 2578630) −` and JMJJOG_04733
    `[4825395, 4825398) +`. Each is `candidate` with `n_alt` 1 (its 28-mer). The 16-mer and
    27-mer are `shared`, because the copies differ in length (1,015 vs 1,016 aa), and
    `pileup.eligible` drops them;
  - fdhF (JMJJOG_04490) has no row: no K/R within 30 aa of its `U`, so no placeable tryptic
    peptide spans it;
  - the note is rule 25's text.
- **TE1:**
  - `misses` stays 4 and the share stays 0.99998;
  - the 4 exceptions become `kind = sec_site`, with `hits` `tolerant:cds_inframe:<copy>` (the
    16-mer and 27-mer list both copies), and `counts["sec_site"]` = 4;
  - the 4 rows record `placement_q`, `m_target` (1) and `accepted`.
- **Rejection is a finding.** The expectations above hold only if all 4 Sec peptides are accepted
  (`placement_q ≤ α`). If one is rejected, it stays `kind = miss` with `hits = "unplaced"` and
  `accepted = False`, and TE2 has one site fewer. That goes to the lead as a finding. It is
  never a reason to change α or any criterion (R3), and the test then records the value rather
  than asserting the count.

### Test plan (16c; test-author; new files unless noted)

Following the statistician's Q5 item 6. Slow and Slurm cases are marked.

- **C1** `tests/unit/test_readthrough_recoded.py`: examples 21–24 (both strands, linear and
  circular, two cells at the origin, the round trip of every recoded cell).
- **C2** same file: example 25 (U over TAA, O over TGA: skipped, counted, no raise) and
  example 26 (length mismatch, terminal kept).
- **C3** same file: example 27, the nested merge (one site, one cell, `kind` recoded, tags
  merged, event counts), and two features recoding one codon.
- **C4** same file: example 28 (O over TAG gives a TAG site; pseudogene, residue 0, `X` and `*`
  ignored), the origin-crossing codon, several `U`s in one protein, and the sort order.
- **C5** same file: `recoded=False` equals the 16a output on every C1–C4 input and on the 16a
  toys, and empty input gives `((), 0)` and zero skips.
- **C6** property (`tests/property/test_recoded_props.py`): the 16a decoy-symmetry property (a)–(c),
  re-run with random genomes carrying random `U`-over-TGA proteins (both strands, both
  topologies). Per frame, the target and every decoy have equal cell counts, and every decoy
  cell holds `*`. `mirror_cells` is an involution for reverse and the identity for shuffle. The
  reversed-query equality holds after interval mapping.
- **C7** `tests/unit/test_variants_recoded.py`: M7-01 `align_hits` accepts a `*` at a recoded codon in
  `site_keys(annotated_stops(..., recoded=True))` and still raises for a `*` at an internal codon
  outside the set. Example 2's variants row: ALT `U`, `candidate`, flags
  `stop_codon;vcf_unrepresentable`, no VCF record.
- **C8** `tests/unit/test_mapping_outputs_recoded.py`: example 20 (`sec_site`, no `readthrough`,
  `recoded_codon=2578628`); the `pyl_site` row; a row consuming a terminal and a recoded cell
  carries both flags and both GFF3 attributes; example 30 (W over TGA gives `sec_site` and an
  ordinary W ALT).
- **C9** `tests/unit/context/test_placement_kind_rows_recoded.py`: the 14 rule-17 rows, including
  `placement_flag` `sec_site` and `pyl_site` (counts and `n`); the schema `flags` pattern accepts
  `sec_site` and `pyl_site`.
- **C10** `tests/unit/test_te_proxies_recoded.py`: example 31; `te1_result` without `flags`
  equals today's result except for `sec_site` 0; the `per_peptide` fields; `te2_result`
  `alt_stop_codon` and the rule-25 note string.
- **C11** `tests/unit/test_gippasim_selenoproteins.py` (tooling test): example 32's truth, genome
  and annotation checks; the codon is ≥ 10 codons from both ends and overlaps no other feature;
  the expressed protein has `U`; the 0 setting is byte-identical; shortfall and `strain` raise.
- **C12** `tests/integration/test_recoded_sim.py` (integration; `slow` only if over the budget):
  sim100k, `n_selenoproteins = 5`, through `sixframe_index` → `map_tolerant` → `loci` →
  `variants_step`. Every truth site whose spanning true peptide is accepted has a `sec_site` row
  and a `U` ALT site at its codon. No accepted row consumes a recoded cell that is absent from
  the truth file. The sixframe and variants site sets are equal.
- **C13** (`slow`, Nissle recording): example 29 (example 10's 16c column).
- **C14** (`slow`, Slurm): TE1 and TE2 on Nissle with the numbers in "TE2 and TE1 expectations",
  and the 4 peptides' `placement_q` recorded in the result JSON.
- **Not the test-author's:** attack (d), d1–d4 (the statistician; see "Re-runs").

### Re-runs (16c; replaces "fast tier suffices")

Every M3-21 suite runs on gippasim, which has no `U`. With `n_selenoproteins = 0` they are
byte-identical by construction, so they are regression evidence only. Required, with job IDs
in PROGRESS.md, `placement-fdr.md` (16c paragraph) and ADR-0058 Evidence:

1. **Regression, fast tier** (Slurm `test_branch.sbatch`): `test_placement_fdr_sim.py`,
   `test_denovo_precision_sim.py`, attack 5, the stressed arm and the decoy-balance check. The
   FDP tables are identical to main's and the cell counts are unchanged. No full tier.
2. **Attack (d)** in `tests/calibration/attacks/test_attack_readthrough_recoded.py` (a new file), written by the
   statistician once the 16c API exists, with no new tolerance key:
   - **d1:** a recoded-dense gippasim world with ≥ 50 recoded cells per replicate, kind-split
     masses, and criteria (a)–(c) plus the structural check on the recoded part;
   - **d2:** three broken builds that must fail (decoys without recoded cells; recoded cells
     mirrored off by one; target without recoded cells);
   - **d3** (full tier only): Nissle and K-12 structural equality and reverse i.i.d. exactness,
     with the shuffle and tryptic shifts logged and only the ≥ −EQUIV margin asserted;
   - **d4:** `load_stop_cells` sites equal the `variants_step` `site_keys`.

   **(d) passing in both the fast and the full tier is a merge prerequisite.**
3. **TE1 and TE2** (slow integration on Slurm), on Nissle, with the expectations above.
4. **TE8** once on the branch (weekly job, or `te_proxies.py te8` on Slurm). Record the pooled FDR
   (≤ 0.01) and the change in accepted tolerant placements. `k12_exact` does not change.
5. **Example 29** (C13).

No full M3-21 tier, no stressed-arm run beyond item 1, and no 16b observer re-run.

### Existing tests the lead edits (checked on 267d0130)

- `tests/integration/test_te2_proxy.py`: the line-122 comment ("stay unmapped until M7-16"), and
  `COUNT_KEYS` gains `alt_stop_codon`. Optionally assert `alt_stop_codon == 2` on Nissle.
- `tests/unit/test_te2_proxy_result.py`: `COUNT_KEYS` and the exact `counts` dict (line 93) gain
  `alt_stop_codon`. The module docstring's rule 2 ("the Sec sites stay unmapped") changes.
- `tests/unit/test_weekly_te2.py`: checked, it pins no note text or count keys. No edit is
  expected; listed because the brief named it.
- `tests/unit/test_te_proxies.py::test_te1_result_counts`: the exact `counts` dict gains
  `"sec_site": 0`. The other TE1 tests keep `kind = miss` (their `HITS` have no `flags` column).
- `tests/unit/test_mapping_step_readthrough.py` (lines 77 and 96): the exact event dicts gain
  `recoded_skipped`. Their `recoded` value of 0 stands, because the toys have no `U`.
- `tests/integration/test_mapping_readthrough.py`: line 160's event dict gains `recoded_skipped`.
  Example 10 (line 210, `annotated_stops(gm, genome)` with the default) stays 4,689 and needs no
  edit. It changes only if the lead moves it to the step path, and C13 covers the 16c column.
- `tests/unit/context/test_placement_kind_rows.py`: `ORDER` gains `("placement_flag",
  "pyl_site")` and `("placement_flag", "sec_site")`, which makes 14 rows.
- Goldens with the rule-17 rows gain the two `placement_flag` rows (value 0) and are regenerated:
  `tests/data/golden/sim50k/m4/detection_summary.tsv`,
  `tests/data/golden/sim100k_real/m4/detection_summary.tsv`,
  `tests/data/snapshot/sim_plots/OUT/tables/detection_summary.tsv`.
- `tests/unit/test_peptide_hits_schema_flags.py` already allows `sec_site`/`pyl_site`. No edit.

### Size (16c; counted lines: `src`, `workflow`, `scripts`, `validation`; ADR-0007)

| File | Change | Lines |
|---|---|---|
| `src/gippa/mapping/readthrough.py` | recoded scan, merge, `recoded_skips` | 45 |
| `src/gippa/steps/mapping.py` | `recoded=True`, event counts, `recoded_skipped` | 8 |
| `src/gippa/steps/variants.py` | `recoded=True` | 1 |
| `src/gippa/mapping/loci.py` | kind/codon → flag, `recoded_codon` work column | 15 |
| `src/gippa/mapping/outputs.py` | GFF3 `recoded_codon` | 8 |
| `src/gippa/context/summary.py` | two `placement_flag` groups | 4 |
| `src/gippa/schemas/tables/peptide_hits.schema.json` | pattern values | 1 |
| `validation/proxies/te_proxies.py` | TE1 kind, count, `per_peptide`; TE2 note, `alt_stop_codon`; docstring | 30 |
| | **counted total** | **≈ 112** |
| `tests/sim/gippasim/{simulate,layout,proteomics,truth}.py` | `n_selenoproteins` (tooling, not counted) | ≈ 70 |
| tests (C1–C14) | not counted | ≈ 450 |

≈ 112 counted lines fits under the ~400 cap: **one PR, no split.**

### Questions (16c; each has a default, none blocks the test-author)

- **Q21 (lead).** The skip counts go into `mapping.stop_sites` as `recoded_skipped`, which means
  editing two exact-dict tests. **Default:** the event field, as briefed. Alternative: a separate
  WARNING event `mapping.recoded_skipped`, emitted only when a count is > 0, which needs no test
  edit.
- **Q22 (lead).** A merged terminal/recoded codon gets `kind = recoded` and therefore `sec_site`.
  **Default:** recoded, because the protein evidence says the codon is decoded. Alternative:
  `terminal` with `readthrough`.
- **Q23 (Khoa, batch).** No Pyl planting in gippasim. **Default:** Sec only, and `pyl_site` is
  unit-tested only, a stated limitation. Alternative: an `n_pyrrolysine` knob (≈ 15 more
  simulator lines, the same rules with `O`/TAG).
- **Q24 (lead).** `not_a_stop` features (3′-partial) are not scanned for `U`/`O`. **Default:**
  not scanned. Their protein length cannot be checked against a terminal stop. On Nissle the 3
  are plasmid-edge partials with no `U` (statistician count: 3 recoded sites in total).
- **Q25 (lead).** `te1_result` without a `flags` column treats the rows as unflagged.
  **Default:** tolerated, so existing callers and tests keep their meaning. Alternative: require
  the column and edit the TE1 test fixtures.
- **Q26 (lead).** The four Sec peptides' q-values are recorded in the TE1 result JSON through
  `per_peptide`. **Default:** in the JSON, so a test can check it. Alternative: record them by
  hand in the PR body and ADR-0058 Evidence, with no code.
- **Q27 (lead).** The new TE2 count `alt_stop_codon` makes the "two candidate U sites" checkable.
  **Default:** add it (2 lines plus the two `COUNT_KEYS` edits). Alternative: rely on
  `alt_sites` 235.
- **Q28 (lead).** No index format bump. A 16b-built format-2 index has no recoded cells, and
  `-resume` could reuse one across the upgrade. **Default:** no bump, because the index is a
  per-run artifact. Alternative: format 3, or a `recoded` meta key that `load_stop_cells` checks.
- **Q29 (lead).** The #335 decision put 16b-3 (`n_readthrough`) "before 16c". 16c has no code
  dependency on it. **Default:** 16c proceeds in parallel and merges main after 16b-3 if both
  are ready. If not, the lead's ordering stands.

## Public interface and touch points

```python
# src/gippa/mapping/sixframe.py
READTHROUGH: int = 22          # target-only code of a marked (consumable) stop; never in frames
def load_stop_cells(path: Path) -> StopCells: ...            # format 2; format 1 raises
# save_index(path, index, decoys, cells: StopCells | None = None)  -> writes format 2

# src/gippa/mapping/readthrough.py   (new, pure)
@dataclass(frozen=True, slots=True)
class StopSite:
    contig: str; strand: Literal["+", "-"]; codon_start: int
    locus_tags: tuple[str, ...]; codon: Literal["TAA", "TAG", "TGA"]
    kind: Literal["terminal", "recoded"]
@dataclass(frozen=True, slots=True)
class StopCells:
    sites: tuple[StopSite, ...]
    target: tuple[np.ndarray, ...]                       # per frame_id, sorted int64 aa positions
    decoys: Mapping[DecoyType, tuple[np.ndarray, ...]]   # mirrored, same frame order
def annotated_stops(gm: GeneModel, genome: Mapping[str, str]) -> tuple[tuple[StopSite, ...], int]
    # (sites, n_not_a_stop)
def frame_cells(index: SixFrameIndex, sites: Sequence[StopSite]) -> tuple[np.ndarray, ...]
def mirror_cells(kind: DecoyType, frames: tuple[Frame, ...],
                 cells: tuple[np.ndarray, ...]) -> tuple[np.ndarray, ...]
def site_keys(sites: Sequence[StopSite]) -> frozenset[tuple[str, str, int]]
```

| PR | File | Function / change |
|---|---|---|
| 16a | `src/gippa/mapping/readthrough.py` (new) | the above (rules 1–3) |
| 16a | `src/gippa/mapping/sixframe.py` | `READTHROUGH`; `save_index` format 2 + cells; `load_stop_cells` |
| 16a | `src/gippa/mapping/tolerant_dp.py` | `_verify_core` (1, 1) case (rule 4); module docstring |
| 16a | `src/gippa/mapping/tolerant.py` | `Search(frames, marked=None)`; `match_tolerant(..., cells=None)` |
| 16a | `src/gippa/mapping/frameshift.py` | `_FSearch(..., marked)`; `nonfs_overlap` (1, 1) case; `match_frameshift(..., cells=None)` |
| 16a | `src/gippa/steps/mapping.py` | `sixframe_index` input `genes`, event; `map_tolerant` passes cells |
| 16a | `workflow/subworkflows/local/mapping.nf` | `stepIn('sixframe_index', [g, pep, gm[2], gm[0]])` |
| 16a | `docs/methods/placement-fdr.md` | ADR-0058 section; required content in "Methods note (16a)" below |
| 16b | `src/gippa/variants/align.py` | `_residues`, `_dp`, `align_hit`, `align_hits` (`stop_sites`) |
| 16b | `src/gippa/variants/pileup.py` | rule 4c narrowed (U/TGA, O/TAG → ALT) |
| 16b | `src/gippa/variants/vcf_records.py`, `flags.py` | `U → TGA`, `O → TAG` zero-change; no `nt_class`, `vcf_unrepresentable` |
| 16b | `src/gippa/steps/variants.py` | build `site_keys(annotated_stops(...))`, pass to `align_hits` |
| 16b | `src/gippa/novel/evidence.py` | `cterm_extension` rule (rule 12) |
| 16b | `src/gippa/novel/outputs.py` | `patch_features` cterm rows |
| 16b | `src/gippa/makedb/build.py` | patch `X` at the annotated stop (rule 13) |
| 16b | `src/gippa/mapping/readthrough.py` | `consumed_sites` (rule 15) |
| 16b | `src/gippa/mapping/loci.py` | `annotate_peptide_hits(…, consumed=None)`: `readthrough` flag, `readthrough_codon` work column (rule 16) |
| 16b | `src/gippa/steps/mapping.py` | `loci` optional input `sixframe`, `load_stop_cells`, WARNING `mapping.readthrough_unchecked` |
| 16b | `workflow/subworkflows/local/mapping.nf` | `stepIn('loci', [...])` gains `six` |
| 16b | `src/gippa/mapping/outputs.py` | GFF3 attribute `readthrough_codon`; the work column never reaches the TSV |
| 16b | `src/gippa/schemas/tables/peptide_hits.schema.json` | `flags` `pattern` (rule 18) |
| 16b | `src/gippa/context/summary.py`, `src/gippa/steps/context_stats.py` | `placement_kind_rows` (rule 17) as one more summary part |
| 16b | `docs/specs/M3-16-mapping-outputs.md`, `docs/specs/M4-16-detection-summary.md` | amendment notes pointing to rules 15–19 and 17 (written with this amendment) |
| 16b | `docs/methods/placement-fdr.md` (16a section) | a short "16b reporting" paragraph: the flag and the count are counts, not a subset FDR (item 5); the flagged rows are a subset of a list accepted at q and that subset is not FDR-controlled (items 5–6) |
| 16b | `tests/sim/gippasim/{config,features,proteomics,truth}.py` | `n_readthrough`, `truth/readthrough.tsv` (tooling) |
| 16b | `tests/sim/gippasim/{strain,proteomics,truth}.py` | `StrainConfig.n_stop_loss`, `truth/stop_loss.tsv`, C-terminal wrong peptides first (rule 20; tooling) |
| 16b | `validation/proxies/te_proxies.py` | TE2 note text (internal Sec sites are not annotated stops; ADR-0058) |
| 16b | `docs/methods/variants.md`, `docs/methods/novel-orfs.md` | stop_codon reachability, U/O ALT, cterm calls; `variants.md` states that REF `*` and U/O zero-change ALTs pool in `(substitution, mass_equivalent = false)` and that their `class_error_rate` does not measure read-through reads (rule 11) |
| 16c | `src/gippa/mapping/readthrough.py` | `annotated_stops(gm, genome, recoded: bool = False)`; `recoded_skips` (rules 21–22) |
| 16c | `src/gippa/steps/mapping.py`, `steps/variants.py` | both pass `recoded=True`; event `recoded`, cell totals, `recoded_skipped` (rule 22) |
| 16c | `tests/sim/gippasim/{simulate,layout,proteomics,truth}.py` | `n_selenoproteins`, `truth/selenoproteins.tsv` (rule 23; tooling) |
| 16c | `validation/proxies/te_proxies.py` | TE1 `kind = sec_site`, count `sec_site`, `per_peptide` q-values (rule 24); TE2 note, `alt_stop_codon` (rule 25) |
| 16c | `src/gippa/mapping/loci.py`, `outputs.py`, `context/summary.py`, `peptide_hits.schema.json` | `sec_site`/`pyl_site` flags, `recoded_codon`, summary groups, pattern values (rule 19) |

Edits to **existing** tests are the lead's (rules.md). Expected ones: the step-registration test
(`sixframe_index` inputs) and mapping.nf nf-test channel wiring; `test_variants_align` tests that
pin "`*` raises" (it now raises only off the site set); M7-18 tests pinning U-over-TGA `rare`
(examples 7–8 of M7-18 change meaning); M7-09/M7-11b tests pinning "`cterm_extension` never
called"; M7-07 "no cterm patch row"; M3-23 "internal stop raises" for cterm rows; the TE2
note-text test.

Added by the amendment (16b), checked on main c523c1c0:

- **Goldens and snapshots with `detection_summary.tsv`** gain the rule-17 rows and must be
  regenerated: `tests/data/golden/sim50k/m4/detection_summary.tsv`,
  `tests/data/golden/sim100k_real/m4/detection_summary.tsv`,
  `tests/data/snapshot/sim_plots/OUT/tables/detection_summary.tsv`,
  `tests/data/snapshot/sim_empty/OUT/tables/detection_summary.tsv`, and the tests that compare
  against them. `tests/unit/context/test_detection_summary_qc.py::test_qc_rows_shape` pins
  `qc_rows` to group `run`/`all`; it stays valid because the new rows come from a separate
  function. Tests that pin the **whole** assembled summary's row set or count (for example the
  context-stats step tests and `tests/helpers/report_content.py` consumers) may need the new
  rows; the test-author lists the exact failures from one run on the branch.
- **`mappings.gff3` / `peptide_hits.tsv` goldens** (`tests/data/golden/sim50k/mapping/`, the
  snapshot `OUT/mapping/` files) change only if a row there consumes a cell. sim50k is same-strain
  with 905 `same_db` rows and 1 `low_complexity`; expected unchanged, to be confirmed.
- `tests/unit/test_mapping_nf_registry.py` and `tests/unit/test_mapping_steps.py`: they use
  subset checks on `loci` inputs. An optional `sixframe` should not break them; if the registry
  test compares `stepIn` keys with the declared inputs exactly, it needs the new key.
- `tests/unit/test_schemas.py`: no existing assertion on peptide_hits `flags`; the pattern is new
  and gets new tests. Every committed `peptide_hits.tsv` must still validate.

## Methods note (16a): the ADR-0058 section of `docs/methods/placement-fdr.md`

The statistician signs it off before 16a merges (review change 3). It must state all of the following:

1. **The rule:** which stops are consumable (ADR-0058 item 1), at what cost, and that peptides
   shorter than 10 aa over a stop are never placed.
2. **Counts are exactly symmetric.** Each decoy frame has the target's number of marked cells, each
   holding `*`, with the same stop-free run lengths on each side. Origin copies are included.
3. **Context mirroring and the measured ratios.** The reverse decoy is exact under any i.i.d.
   peptide law (D_rt/T_rt = 1 to 1e-16). But the local context is mirrored, not kept: gene | `*` |
   3′ flank becomes reversed flank | `*` | reversed gene. Under a law that is not
   reversal-symmetric (tryptic, C-terminal K/R), the per-genome read-through D_rt/T_rt was
   0.52–1.89 (reverse) and 0.58–1.74 (shuffle) on gippasim 1 Mb, and 0.81–1.15 (reverse) and
   1.14–1.42 (shuffle) on Nissle. The net effect on the total, (D_rt − T_rt)/T_regular, was at
   worst −9.4 × 10⁻⁴ (gippasim) and −1.3 × 10⁻⁴ (Nissle). Cite Slurm jobs 1889319 (Nissle) and
   1889320 (gippasim), and the stressed arm's measured values once it has run.
4. **The added mass is an approximation, not a bound:** ≈ G/(N(1 − s)) at m = 1 and
   ≈ 2G/(N(1 − s)) at m = 2, independent of L, scaling with gene density G/N. The measured realised
   shares were 0.04–0.08 % on Nissle and 0.05–0.20 % on gippasim (i.i.d. and tryptic laws).
5. **The read-through ratio is not a subset FDR.** D_rt/T_rt does not estimate the false discovery
   rate of the read-through placements. The list FDR does not control that subset, and no output,
   table or report may present the ratio as one (R3).
6. **Expected false read-through count per run:** about α · T_acc(m ≥ 1) · ρ, where T_acc(m ≥ 1) is
   the number of accepted target placements in strata 1–2 and ρ ≈ 10⁻³ is the read-through share of
   the chance mass. That count is ≪ 1 for any realistic run, but it is not a per-site guarantee.
7. **Power note.** A true read-through peptide with no other difference lands in stratum 1 with
   π̂ = 1 and the +1. By Appendix B.6, a stratum-1 list with fewer than about 100 targets ahead of
   its first decoy accepts none. So `cterm_extension` has low power on small inputs by design. This
   is a statement, not a change.
8. **What the suites show:** the exercised re-runs and the stressed arm (test plan), with job IDs,
   and the per-replicate count of accepted target rows that consumed a cell. The regression-only
   suites are mentioned only as regression evidence (they cannot see the rule).

## Test plan (test-author; new files unless noted)

16c's cases are C1–C14 in "16c amendment (2026-10-25)" › "Test plan (16c)".

- `tests/unit/test_readthrough_sites.py`: rules 1–3. Cover both strands, linear and circular,
  the origin-crossing codon, the two copies at `[0, 3)`, shared stops, sORF, pseudogene excluded,
  `not_a_stop`, the round-trip assertion, and mirrored cells holding `*` (examples 8, 10, edge
  cases).
- `tests/unit/test_tolerant_dp_readthrough.py`: examples 1, 3, 4, 4b, 4c at `verify` level.
  A marked cell inside every block shape (1↔2, 2↔1, 2↔2) and swap is never consumed. A `NOMATCH`
  query residue (`U`) over a marked cell costs 1. The same inputs with code 19 give the M3-06
  results. `NUMBA_DISABLE_JIT=1` subset agrees.
- **Property: zero false negatives with the new rule** (`tests/property/brute_readthrough.py` +
  `test_readthrough_props.py`). This is an independent brute-force matcher written from this spec.
  It must not import `gippa.mapping.readthrough` or `tolerant_dp`, and the existing oracles are not
  edited. It is the unbanded M3-07 brute force with a set of consumable cells per frame, where a
  cell can be consumed only by a substitution. Random genomes ≤ 50 kb (linear and circular, both
  strands) carry random annotated CDS. Queries come from planting 0–2 differences into windows
  spanning a marked cell, from windows spanning an unmarked stop, and from random queries.
  Criterion (exact, no tolerance): the minimal cost equals the brute force, every brute-force
  cluster is reported, and every reported window costs the minimum. Also: no reported window ever
  consumes an unmarked `*`. Hypothesis `@example`s cover each diff type next to a marked cell, and
  cover the wrap.
- **Property: decoy symmetry.** (a) Per frame, the count of marked cells is equal in the target
  and each decoy, and every decoy cell holds `*`. (b) `mirror_cells` is an involution for
  `reverse` and the identity for `shuffle`. (c) Target search of `reversed(Q)` on the target
  frames, with target cells, equals the decoy search of `Q` on `reverse` decoys with mirrored
  cells, after interval mapping (the analogue of `test_target_search_on_decoy_frames_equals_decoy_search`).
  (d) Example 5's broken build fails (c). (e) **Frameshift route (review change 9):** a unit case
  where the decoy frameshift search builds `_FSearch` with that decoy type's mirrored cells, so
  `nonfs_overlap` suppresses a decoy split exactly as it suppresses the target split of the
  reversed query. Extend (c) to `match_frameshift` output; if the reverse adjacency map makes that
  equality not apply, the test file's docstring says why, and a count-level property (equal
  numbers of suppressed splits per target/decoy pair) is the minimum.
- **Determinism of the decoy cells (review change 7)** (`tests/unit/test_readthrough_sites.py`):
  `mirror_cells` draws no randomness (a spy on `gippa.rng.for_step` and on `numpy.random` sees no
  call, and the shuffle stream of the step yields the same next draw with and without the call).
  For the same seed, the decoy frames saved in index format 2 are byte-identical to those of
  format 1. The existing golden decoy tables do not move.
- **Property: no-op.** With empty cells, `match_tolerant`/`match_frameshift` return frames equal to
  the call without the keyword, over the existing property strategies.
- `tests/unit/test_frameshift_readthrough.py`: rule 7 (a split side never consumes a marked cell;
  precedence suppresses an overlapping split at ≥ cost).
- `tests/integration/test_mapping_readthrough.py`: `sixframe_index` → `map_tolerant` on sim100k with
  `n_readthrough = 5`. Index format 2 round trip, format 1 raises, event counts, cterm rows.
- 16b: `test_variants_readthrough.py` (examples 1, 2 under 16c, 9; rule 11 per codon); M7-01
  equality on every read-through hit; `test_novel_cterm.py` (examples 6, 7, 7b); an e2e local run
  on sim100k with `n_readthrough = 10` (`called` cterm rows, `readthrough_candidate`, patch round
  trip through `makedb --patch`, `stop_codon` ALT sites matching `truth/readthrough.tsv`).
- 16b amendment:
  - `tests/unit/test_readthrough_flag.py`: `consumed_sites` on examples 11–16 and the edge cases
    (both strands, linear and circular, the origin-crossing codon, two cells, frameshift rows, empty
    sites).
  - **Property:** on the hit tables of the existing tolerant property strategies with random
    cells, `consumed_sites(row) != ()` ⇔ `holds_cell(row)` for every target row.
  - `tests/unit/test_mapping_outputs_readthrough.py`: the flag order, the GFF3 attribute
    (1-based, ascending, absent when not flagged, absent on `match_part`), the TSV column tuple
    unchanged, and the work column not written.
  - `tests/unit/context/test_placement_kind_rows.py`: examples 17 and 18, the fixed row set, the
    NA rows without `unmapped`, and the consistency checks of the acceptance criteria.
  - The `loci` step without `sixframe` emits `mapping.readthrough_unchecked`; with a format-1
    file it raises.
  - **e2e, strain stop-loss (D5, rule 20)** (`tests/integration/test_readthrough_strain.py`, not
    `slow` if it fits in the integration budget): sim100k, `related`, `StrainConfig(
    snv_rate_per_kb=5, n_stop_loss=5)`, `wrong_source = "other_strain"`, through `run_chain` (or
    the step calls) → `loci` → `variants_step`. For every `truth/stop_loss.tsv` row whose
    `cterm_peptide` is accepted: its `peptide_hits` row is flagged `readthrough` and ends on the
    site's codon, and `variants.tsv` has a `stop_codon` ALT site at `codon_start` with
    `alt_aa = reference_residue` and the `strain_difference` flag. "Ends on the codon" is in
    the gene direction: `end = codon_start + 3` on `+`, `start = codon_start` on `−` (unwrapped
    on circular replicons). Any accepted row that consumed a cell at a site in neither
    `truth/stop_loss.tsv` (planted or not) nor `truth/readthrough.tsv` fails the test. The seed
    is fixed; a failure of that check is a finding (a false read-through placement) to report,
    never a reason to change the seed. The test asserts that at least 3 of the 5 **planted**
    rows are accepted, so it can never pass vacuously; if the build gives fewer, the
    test-author raises π_wrong or k and reports it, and never drops the assertion.
  - **Calibration observer (D5):** in the `related_other_strain` full-tier observer that counted
    the 0–5 rows (methods note item 8), each counted row is split by the suite's own truth
    scoring into rows placed at their peptide's origin locus and the rest. Every origin row's
    site must have a sample codon that differs from the reference codon at its lifted position
    (a strain stop-loss); under the gippasim strain model this is implied, not sampled (rule 20),
    so it is asserted without a tolerance. Off-origin rows are chance read-through placements:
    they are already false placements in the FDP, so they are counted and reported (expected
    ≈ 0: α · T_acc(m ≥ 1) · ρ, methods note item 6) and not asserted. That confirms the
    strain-model finding on the real suite. It is a diagnostic assertion, not a new gate, and
    no tolerance key is added.
- **Statistical re-run (R3: no tolerance, bound, α or R change).** Run the M3-21 placement suites
  on the 16a branch, fast (`-m calibration`, `scripts/slurm/test_branch.sbatch`) and full
  (`calibration_full.sbatch` placement tasks, `GIPPA_WORKTREE` + `PYTHONPATH` exported), in two
  groups (review change 1):
  - **Regression-only** (byte-identical by construction: they pass `cells=None` or use no
    matcher): `test_tolerant_recall.py` (calls `match_tolerant` without `cells`),
    `test_placement_fdr_genome.py` (`build_sixframe → match_exact`, exact route only),
    `test_placement_fdr_lists.py` (statistics only) and `attacks/test_attack_placement_*` (`pi`,
    `same_db`, `contaminants`, `strain_span`: `placement()` on synthetic cost tables). They must
    pass unchanged. They are **not** evidence for the rule, and ADR-0058 Evidence must not cite
    them as such.
  - **Exercised** (the `run_chain` suites, which reach `sixframe_index → map_tolerant` with cells):
    `test_placement_fdr_sim.py`, `test_denovo_precision_sim.py` and attack 5
    (`attack5_reversed_stretches`). Only this group is evidence. For each, the PR records the job
    IDs, the before → after FDP tables, and the **per-replicate number of accepted target rows
    that consumed a cell** (expected near 0 at 100 kb, where cells are about 5 × 10⁻⁴ of the
    chance mass).
- **Decoy-balance check** (`tests/calibration/test_readthrough_decoy_balance.py`, statistician
  signs off; reuses `replicates_*.placement` and `se_multiplier`, no new tolerance key; review
  change 2). On gippasim 1 Mb null genomes (no read-through planted) with the M3-19/M3-20
  wrong-peptide arms, it counts per replicate and **per decoy type** (`reverse`, `shuffle`) the
  false target placements whose window holds a marked cell (T_rt), the decoy rows holding a
  mirrored cell (D_rt), and the read-through share T_rt / T (all false target placements). It
  reports all three per decoy type. Criterion, one-sided (it fails only in the anti-conservative
  direction): mean(D_rt) ≥ mean(T_rt) − `se_multiplier` · SE(D_rt − T_rt), the ADR-0019
  estimate-check form. It is reported as exercised or vacuous like the M3-21 lists; it is expected
  vacuous at 1 Mb (T_rt ≈ 0). **"Vacuous" is an allowed outcome for this check only.** It is no
  longer the evidence for symmetry; the stressed arm and the statistician's attack are.
- **Stressed arm (required; review change 5; test-author writes it, the statistician signs off)**
  (`tests/calibration/test_readthrough_stressed.py`; tiers `replicates_fast.placement` /
  `replicates_full.placement`; existing keys `levels`, `se_multiplier` only). A sampled arm in
  which T_rt is not vacuous:
  - **Path:** the real chain `sixframe_index → map_tolerant → best_costs → placement` (through
    `run_chain` or the same calls), never a synthetic cost table.
  - **Build:** a synthetic genome with dense short annotated genes (for example 30–40 codons each,
    each ending in a stop, so G/N ≥ 50× the biological density) and a reduced-alphabet coding
    frame (s about 0.2–0.3). Wrong 10–12-mers drawn from the coding-frame law. The build must give
    **≥ 100 expected false target read-through rows per replicate**; the test asserts the realised
    replicate mean of T_rt is ≥ 100 before it evaluates any criterion (otherwise it fails as
    mis-built, never as vacuous).
  - **Arms:** a `reverse` arm (database input) and a `shuffle` arm (de novo input).
  - **Criteria, per arm:** (i) FDP in the m ≥ 1 strata ≤ q + `se_multiplier` · SE at every q in
    `levels`; (ii) mean(D_rt) ≥ mean(T_rt) − `se_multiplier` · SE(D_rt − T_rt).
  - **Negative control:** the same arm with the decoy cells removed (decoy searches without marks)
    must **fail** criterion (ii). The test asserts that failure; a passing control fails the test.
- **Statistician's attack (review change 6)**
  (`tests/calibration/attacks/test_attack_readthrough_decoys.py`). The statistician writes it once
  16a's `frame_cells`/`mirror_cells` exist; it is the exact-mass version of the review table:
  (a) reverse, i.i.d. law: D_rt = T_rt to rtol 1e-9 per replicate; (b) shuffle, i.i.d. law:
  replicate mean D_rt/T_rt ≥ 1 − 0.01 − `se_multiplier` · SE; (c) tryptic law, both decoys: net
  shift (D_rt − T_rt)/T_regular of the replicate mean ≥ −0.01; (d) the recoded cells, d1–d4 as
  in "Re-runs (16c)", written by the statistician once the 16c API exists (the old "(d) … covers
  them unchanged" was vacuous: gippasim has no `U`). The test-author does not write it. **16a's
  merge waited for (a)–(c); 16c's merge waits for (d), fast and full.**

## Acceptance criteria

- [ ] **16a:** examples 1, 3, 4, 4b, 4c, 5, 8, 10 and the edge cases pass as unit tests. The
  zero-false-negative property, the decoy-symmetry property (including (e), the frameshift route)
  and the no-op property pass. `mirror_cells` draws no randomness, and decoy frames are
  byte-identical between index formats 1 and 2 for the same seed. The M3-06 benchmark (10⁶
  `verify` calls ≤ 5 s) still passes.
- [ ] **16a:** the M3-21 placement suites pass fast **and** full on the branch with no edit to
  `tests/calibration/tolerances.yaml`, `ci/gates.yaml` or any criterion. ADR-0058 Evidence and
  PROGRESS.md list the regression-only and exercised groups separately. Only the exercised group
  (`test_placement_fdr_sim.py`, `test_denovo_precision_sim.py`, attack 5) is cited as evidence,
  with Slurm job IDs, before → after FDP tables and the per-replicate count of accepted target
  rows that consumed a cell.
- [ ] **16a:** the decoy-balance check passes or is recorded vacuous, with T_rt, D_rt and the
  read-through share reported per decoy type.
- [ ] **16a:** the stressed arm passes fast **and** full in both the reverse and the shuffle arm:
  realised mean T_rt ≥ 100 per replicate; m ≥ 1 FDP ≤ q + `se_multiplier` · SE at every level;
  mean(D_rt) ≥ mean(T_rt) − `se_multiplier` · SE. Its negative control (decoy cells removed)
  fails the balance criterion.
- [ ] **16a:** the statistician's `tests/calibration/attacks/test_attack_readthrough_decoys.py`
  is on the branch and passes. 16a does not merge before it.
- [ ] **16a:** the `placement-fdr.md` section has every item of "Methods note (16a)", and ADR-0058
  states the added mass as an approximation (scaling formula plus measured shares), not a bound.
- [ ] **16a:** on the Nissle recording `mapping.stop_sites` reports 4,689 terminal sites, 3
  `not_a_stop` and equal target/decoy cell counts. TE1 and TE8 re-run unchanged (TE1 ≥ 0.999,
  the 4 fdnG misses still unmapped; TE8 FDR ≤ 0.01).
- [ ] **16b:** examples 6, 7, 7b, 9 and rule 11 pass. M7-01 recomputes every read-through hit's
  cost. On sim100k with `n_readthrough = 10`, every planted gene with ≥ 2 distinct spanning true
  peptides gives a `called` `cterm_extension` with `readthrough_candidate`. Its patch row
  round-trips through `makedb --patch` with `X` at the old stop. Each planted stop with a placed
  spanning peptide is a `stop_codon` ALT site with the planted residue.
- [ ] **16b:** TE2 re-run passes with 0 violations. Its note says why no `stop_codon` site appears
  on Nissle under items 1–5.
- [ ] **16b:** `variants.md` states that REF `*` sites and U/O zero-change ALTs are in class
  `(substitution, mass_equivalent = false)` and that their `class_error_rate` is that class's
  pooled rate, which does not measure read-through reads, also after a patched `X` search, and
  does not feed `confidence` (rule 11; Q11 answered: pooled).
- [ ] **16b (amendment, rule 15):** `consumed_sites` passes examples 11–16. On the property
  strategies it agrees with `holds_cell` on every target row. Exact, decoy and frameshift rows
  are never flagged.
- [ ] **16b (rule 16):** with `sixframe` wired, every accepted row that consumes a terminal cell
  has `readthrough` in `flags`, and no other row has it. `mappings.gff3` carries
  `readthrough_codon` (1-based codon starts) exactly on those features.
  `PEPTIDE_HITS_COLUMNS` is unchanged. A format-1 index given to `loci` raises. Without
  `sixframe`, `loci` emits `mapping.readthrough_unchecked` and flags nothing. `mapping.nf`'s
  `loci` `stepIn` includes the six-frame index (registry test).
- [ ] **16b (rule 17):** `detection_summary.tsv` has exactly the 12 rule-17 rows (3 `placement_kind` + 4
  `tolerant_diff_type` + 1 `placement_flag` + 4 `unmapped_reason`; 14 after 16c), in that order, with `section = qc` and metric
  `n_peptides`. Examples 17 and 18 hold. Consistency checks:
  - Σ `placement_kind` = distinct `peptide_hits.peptide`;
  - `placement_kind` exact = V05 `peptide/exact/placed`, and tolerant + frameshift = V05
    `peptide/tolerant/placed`;
  - Σ `unmapped_reason` = `unmapped_peptides` + `contaminant_peptides`;
  - `placement_flag` readthrough = distinct peptides with a `readthrough` row.

  `qc_rows`' output is unchanged.
- [ ] **16b (rule 18):** the `peptide_hits` schema `flags` pattern accepts every committed
  `peptide_hits.tsv` and the new flag, and rejects an unknown value. §9.4 is not edited
  (GIPPA_SPEC.md is read-only); the schema description and this spec record the extension.
- [ ] **16b (D5, rule 20):** the strain stop-loss e2e passes. With `n_stop_loss = 0` every
  existing gippasim strain build is byte-identical (a golden or hash test). The related-strain
  observer's assertion holds on the full tier, and the PR records the Slurm job ID.
- [ ] **16c (approved by Khoa, 2026-10-24; amended 2026-10-25):** examples 2 (16c column) and
  21–32 pass as unit/integration tests C1–C13, and every 16c edge case is covered. With
  `recoded=False`, `annotated_stops` equals the 16a output on every input. `mapping.stop_sites`
  on Nissle reports `terminal` 4,689, `recoded` 3, `not_a_stop` 3, `cells_target` 4,692, 4,692
  cells per decoy type and `recoded_skipped` zeros (example 29).
- [ ] **16c (rule 23):** `n_selenoproteins = 0` leaves every gippasim output and committed golden
  byte-identical (no regeneration). With `n > 0` the truth file, the TGA codons, the `U` in the
  annotation and the expressed protein, the ≥ 10-codon margins, the shortfall `ValueError` and the
  `strain` `ValueError` hold (C11). C12 passes.
- [ ] **16c (TE, Slurm, C14):** TE2 on Nissle gives `alt_sites` 235, `high` 80,
  `high_start_codon` 80, `high_stop_codon` 0, `violations` 0, `alt_stop_codon` 2. The two `U`
  sites are JMJJOG_02511 `[2578627, 2578630) −` and JMJJOG_04733 `[4825395, 4825398) +`, each
  `candidate`. fdhF has no row. The note is rule 25's text. TE1 keeps `misses` 4, the share
  0.99998 and `passed`, and its 4 exceptions are `kind = sec_site` with `placement_q`,
  `m_target` and `accepted` recorded. A rejected Sec peptide is reported as a finding, and α is
  not changed (R3).
- [ ] **16c (re-runs):** the five items of "Re-runs (16c)" are done, with job IDs in PROGRESS.md,
  `placement-fdr.md` and ADR-0058 Evidence. The regression FDP tables are identical to main's.
  TE8 pooled FDR ≤ 0.01, with the tolerant delta recorded. **Attack (d), d1–d4, written by the
  statistician after the code exists, passes in the fast and the full tier; 16c does not merge
  before it.**
- [ ] **16c labels (rule 19):** example 20 holds. Rows consuming a recoded TGA/`U` cell carry
  `sec_site` and those consuming a TAG/`O` cell carry `pyl_site`, never `readthrough`, with
  `recoded_codon` in the GFF3. The `placement_flag` rows `sec_site` and `pyl_site` are written.
  On Nissle, the fdnG rows carry `sec_site`; fdhF has 0 rows (no K/R within 30 aa of its
  `U`, so no placeable tryptic peptide spans it). A `W` over a recoded TGA gets `sec_site`
  and a `snv;stop_codon` W ALT (rule 26).
- [ ] Statistician sign-off on this spec and on `placement-fdr.md` **before** 16a code; again on the
  re-run tables before 16a merges.

## Line estimate (counted: `src`, `workflow`, `scripts`, `validation`; ADR-0007)

| PR | Part | Lines |
|---|---|---|
| 16a | `readthrough.py` (sites, cells, mirror, keys) | 95 |
| 16a | `sixframe.py` (code, format 2 save/load) | 40 |
| 16a | `tolerant_dp.py` (one case + docstring) | 12 |
| 16a | `tolerant.py`, `frameshift.py` (keyword, marking, precedence case) | 35 |
| 16a | `steps/mapping.py` + `mapping.nf` | 30 |
| | **16a total** | **≈ 210** (+ ≈ 650 test lines, including the stressed arm and the balance check) |
| 16b | `variants/align.py`, `pileup.py`, `vcf_records.py`, `flags.py`, `steps/variants.py` | 55 |
| 16b | `novel/evidence.py`, `novel/outputs.py` | 30 |
| 16b | `makedb/build.py` | 25 |
| 16b | `te_proxies.py` note | 5 |
| 16b | `readthrough.consumed_sites` (rule 15) | 25 |
| 16b | `loci.py` flag + work column, `steps/mapping.py` optional input + event, `mapping.nf` (rule 16) | 25 |
| 16b | `outputs.py` GFF3 attribute (rule 16) | 10 |
| 16b | `peptide_hits.schema.json` pattern (rule 18) | 2 |
| 16b | `context/summary.py` `placement_kind_rows`, `steps/context_stats.py` part (rule 17) | 45 |
| | **16b total** | **≈ 222** (+ ≈ 115 simulator incl. ≈ 25 for `n_stop_loss`, + ≈ 600 test lines) |
| 16c | `readthrough.py` recoded source, steps (rules 21–22) | 54 |
| 16c | TE1/TE2 (rules 24–25) | 30 |
| 16c | labels: kind → flag, `recoded_codon`, two summary groups, pattern values (rule 19) | 28 |
| | **16c total** | **≈ 112** (+ ≈ 70 simulator, + ≈ 450 test lines; per file in "Size (16c)") |

Total ≈ 544 code lines (16a ≈ 210 as merged, 16b ≈ 222, 16c ≈ 112), against the gated estimate
of 220. The difference is the frameshift precedence (rule 7), the index format and loader, the M7
call paths (16b) and Khoa's reporting decisions (rules 15–19), which the gated estimate did not
count.

**PR split.** 16a (matcher, decoys, index, mapping step, methods note, re-calibration; merged,
#331) → 16b (M7-01 re-alignment, variants U/O rule, cterm call and patch, `makedb` `X`, simulator
planting, the `readthrough` flag and GFF3 attribute, the placement-kind rows, the strain
stop-loss e2e) → 16c (recoded Sec/Pyl codons and their labels; approved, next after 16b). 16b is
≈ 222 counted lines, under the ~400 cap, so it stays one PR. **Fallback split** if the
implementation passes ~400: 16b-1 = variants, novel, `makedb`, the `n_readthrough` simulator and
the TE2 note (≈ 115); 16b-2 = rules 15–18 and 20 (≈ 107 + the `n_stop_loss` simulator), stacked
on 16b-1 because the stop-loss e2e needs the variants path. 16c depends on 16b for the flag
plumbing (rule 19 extends rules 16–18), so it merges main after 16b.

## Questions for the lead (each has a default; none blocks 16a)

1. **(Khoa) Annotated internal Sec/Pyl codons (ADR-0058 item 6).** Option B as worded does not
   reach fdnG, the stated motivation. **Default:** 16a/16b implement items 1–5 only. 16c is fully
   specified but not started until Khoa answers, and PROGRESS.md records "fdnG Sec peptides stay
   unmapped (internal TGA, ADR-0058 item 6 pending)". Batch it with the other Khoa questions.
   **Answered (Khoa, 2026-10-24 evening): approved.** 16c implements item 6 with the labels
   `sec_site`/`pyl_site` (rule 19) and is next after 16b.
2. **Frameshift route.** Should split sides consume a marked stop? **Default:** no; only the
   precedence check sees marked alignments (rule 7).
3. **A window ending on the stop codon** (last residue over `*`, no flank). **Default:** M3-14
   unchanged (`cds_inframe`), so the window gives a `stop_codon` ALT but no `cterm_extension`
   candidate. ADR-0013's "spanning" means extending past the stop.
4. **U over TGA / O over TAG at a marked stop.** **Default:** ALT with `stop_codon` +
   `vcf_unrepresentable`, no `snv`/`multi_nt`, no VCF record (rule 11). This follows §7.5.1 and
   M7-18 decision 1. The alternative (count as REF, M7-15 decision 9) would leave TE2's exception
   vacuous again.
5. **Config switch** (`mapping.tolerant.annotated_stops`)? **Default:** none. The rule is the
   accepted policy; an empty annotation turns it off. A key would need an ADR (defaults.yaml is
   protected).
6. **A hit column marking read-through consumption** (`hits_tolerant`, `peptide_hits.tsv`)?
   **Default:** none (no schema change). Consumers derive it from the interval and the stop sites,
   and `mapping.stop_sites` reports the counts.
   **Answered (Khoa, 2026-10-24 evening): a flag, not a column.** `readthrough` in `flags` plus
   the GFF3 attribute `readthrough_codon` (rules 15–18).
7. **What DIA-NN reports for an `X` in a patched database.** **Default:** no special handling.
   Any letter over a marked cell costs 1, and X/B/Z are skipped in variants (example 9). Record
   it in the methods note, and check it on the first real second-pass search.
8. **Partial and pseudogene features.** **Default:** pseudogenes excluded; a partial CDS is
   included iff its last codon is a stop (rule 1).
9. **Statistical evidence.** **Answered by the statistician (2026-10-24): a stressed arm is
   required.** The added chance-hit mass is an approximation, ≈ G/(N(1 − s)) at m = 1 and
   ≈ 2G/(N(1 − s)) at m = 2 (measured 0.04–0.08 % on Nissle, up to 0.20 % on gippasim under the
   tryptic law), not a bound. The evidence is the exercised M3-21 re-runs, the stressed arm and
   the statistician's attack; the decoy-balance check may be vacuous (test plan).
10. **Index format.** **Default:** bump to 2; format-1 files raise (no silent fallback); Nextflow
    resume re-runs `sixframe_index` anyway because its inputs change.
11. **(Khoa) A separate `stop_codon` error-rate class** (review change 8). REF `*` sites and U/O
    zero-change ALTs pool into `(substitution, mass_equivalent = false)`, whose rate does not measure
    read-through reads. A third class dimension `stop_codon` would change §7.5.4 ("by alt_type and
    by mass_equivalent") and so needs an ADR. **Default (reversible):** pooled, and documented in
    `variants.md` (rule 11). Batch it with the other Khoa questions; it does not block 16a or 16b.
    **Answered (Khoa, 2026-10-24 evening): pooled.** No separate class and no ADR; `variants.md`
    states it (rule 11, D4).
12. **Stressed-arm sizing.** The review sets the target (≥ 100 expected false target read-through
    rows per replicate) and suggests a build. Whether that build runs within the fast-tier time
    budget at `replicates_fast.placement` = 50 is not yet measured. **Default:** the test-author
    sizes the genome for the target and measures one replicate on Slurm. The arm runs in both
    tiers as the acceptance criteria say; if the fast tier is too slow, the test-author reports
    the timing to the lead and the statistician instead of shrinking the arm. No replicate count
    changes (R3).

Added with the amendment (16b, 16c labels). Each has a reversible default, and none blocks the
16b test-author.

13. **GFF3 attribute name and form.** **Default:** `readthrough_codon=<1-based codon start>[,…]`
    (and `recoded_codon` in 16c), following M3-16's `frameshift_reread`. The alternative is a
    bare `readthrough=true`, which adds nothing beyond `flags`.
14. **How `loci` learns the cells.** **Default:** an **optional** `sixframe` input with
    `load_stop_cells`, always wired in `mapping.nf`, plus a WARNING event when it is absent. This
    keeps existing `loci` unit callers valid. Making it **required** is stricter (a missing index
    could not silently drop the flag), but it breaks every test that calls `loci` with a
    hand-built input set, and the lead would have to edit those tests. A third option, computing
    the sites again in `loci` from `annotated_stops(gm, genome)`, is rejected: it could disagree
    with the searched cells (16c's recoded set, and the `features` layer of `_gene_model`).
15. **No ADR for the §9.3/§9.4 additions.** **Default:** none. §9.4 allows extension through the
    schema, and the `detection_summary` rows add no column. ADR-0058's "no new flag value" is
    amended by the lead. If the lead reads §7.2.9's GFF3 attribute set or §9.3's QC list as
    closed, a short ADR (or an ADR-0058 amendment marked accepted-by-khoa, since Khoa made the
    decisions) covers both.
16. **`multi` and `low_complexity` are written to `peptide_hits.flags` but are not in §9.4.**
    **Default:** the new schema pattern lists them as written (rule 18), and a backlog line asks
    whether §9.4 should list them (an ADR, since §9.4 is the controlled vocabulary). The other
    option, dropping the pattern, would leave the vocabulary unchecked, which goes against D1.
17. **Where the counts go.** **Default:** rows in `tables/detection_summary.tsv` (section `qc`,
    rule 17), with no new file. The alternative is a new `mapping/placement_summary.tsv`: a new
    §9.1 file and schema, plus the `lociSpec` stub, written by `loci` next to the other Module 2
    tables. That is cleaner by module, but it is a §9.1/§7.2.9 extension (ADR) and costs about
    40 more lines.
18. **Exclusive `placement_kind` precedence** (exact > tolerant > frameshift) and the
    non-exclusive diff-type groups. **Default:** as in rule 17. An exclusive "dominant diff type"
    would hide `substitution;swap` peptides.
19. **16c: `sec_site`/`pyl_site` replace `readthrough`** on recoded cells. **Default:** they replace
    it, because a recoded codon is decoded, not read through. Both flags appear only when one row
    consumes both kinds. The flags stay placement-scope; variants keep `stop_codon` plus
    `vcf_unrepresentable` (rule 11), and no Variants-scope `sec_site` is added.
20. **D5 e2e floor.** **Default:** at least 3 of 5 planted stop-loss C-terminal peptides
    accepted at sim100k (an exercise floor, not a tolerance; R3 untouched). The test-author
    measures it once and reports the realised count.

## Statistician review (2026-10-24, at e7f7f744)

**APPROVE WITH CHANGES.** 16a code may start; the changes below are required before 16a merges.
Full review: `local/stat/m7-16-review.md`. Exploratory exact-mass runs: Slurm jobs 1889319 (Nissle)
and 1889320 (gippasim 1 Mb × 3), script `local/stat/m7-16/rt_mass.py`.

Findings:
- Decoy symmetry holds exactly in counts. The reverse decoy is exact for i.i.d. laws: the
  read-through D/T was 1 to 1e-16.
- The local context is mirrored, not kept: gene | `*` | flank becomes flank | `*` | gene. Under a
  tryptic law (C-terminal K/R), the read-through D/T per genome is 0.52–1.89 (reverse) and 0.58–1.74
  (shuffle). The net effect on the total chance mass is ≤ 10⁻³, so the list FDR is unaffected.
- D_rt/T_rt is not an FDR for the read-through subset.
- The ADR-0058 added-mass arithmetic is right, but it is an approximation, not a bound. It scales as
  ≈ G/(N(1 − s)) at m = 1 and ≈ 2G/(N(1 − s)) at m = 2. The realised share is 0.04–0.08 % on Nissle
  and 0.05–0.20 % on gippasim.
- The strata, π̂ (m ≥ 1 → 1), E and the §8.1 families are unchanged. The ADR-0019 stratum-0 π̂ is
  untouched, because read-through hits are never at m = 0.
- Item 6 changes none of this: 3 cells, coding context on both sides.

Required changes (numbering as in the review):
1. Split the re-run list into two groups.
   - **Regression-only**, byte-identical by construction: `test_tolerant_recall.py`,
     `test_placement_fdr_genome.py`, `test_placement_fdr_lists.py`, `attacks/test_attack_placement_*`.
   - **Exercised:** the `run_chain` suites `test_placement_fdr_sim.py`, `test_denovo_precision_sim.py`
     and attack 5.
   - Only the exercised group counts as evidence. Record how many accepted target rows consumed a
     cell.
2. Decoy-balance check: use `se_multiplier`, not a literal 2. Report T_rt, D_rt and the read-through
   share per decoy type. Keep it one-sided. "Vacuous" is allowed only for this check.
3. In the `placement-fdr.md` section, state:
   - the context mirroring and the measured ratios;
   - that D_rt/T_rt is not a subset FDR;
   - the expected number of false read-through placements per run, ≈ α · T_acc(m ≥ 1) · ρ with
     ρ ≈ 10⁻³;
   - that power is low on small m ≥ 1 lists (B.6).
4. ADR-0058: replace "≤ 0.11 %" with the scaling and the measured shares.
5. **The stressed arm is required (Q9).** It is a sampled arm through `sixframe_index → map_tolerant
   → best_costs → placement`.
   - Build: dense short annotated genes and a reduced-alphabet coding frame, so there are ≥ 100
     expected false target read-through rows per replicate. Run a reverse arm and a shuffle arm.
   - Criteria: m ≥ 1 FDP ≤ q + `se_multiplier`·SE, and mean(D_rt) ≥ mean(T_rt) −
     `se_multiplier`·SE.
   - A negative control with the decoy cells removed must fail.
6. The statistician adds `tests/calibration/attacks/test_attack_readthrough_decoys.py` (exact mass)
   on top of 16a's `frame_cells`/`mirror_cells`. 16a's merge waits for it.
7. `mirror_cells` draws no randomness. Decoy frames must be byte-identical across index formats 1
   and 2 for the same seed. Add a test.
8. 16b: REF `*` and U/O zero-change ALTs fall into class `(substitution, mass_equivalent)`, and their
   `class_error_rate` is the pooled rate. Document this in `variants.md`. Batch a Khoa question on a
   separate `stop_codon` class (a §7.5.4 change, so it needs an ADR). Until Khoa answers, the default
   is pooled and documented.
9. The frameshift route uses each decoy type's mirrored cells in `_FSearch`/`nonfs_overlap`. Cover
   this with a property or unit test.

**Applied (spec-writer, 2026-10-24):** 1 → test plan "Statistical re-run" and acceptance; 2 →
"Decoy-balance check"; 3 → "Methods note (16a)"; 4 → ADR-0058 Evidence/Consequences and Q9; 5 →
"Stressed arm"; 6 → "Statistician's attack"; 7 → rule 3 and "Determinism of the decoy cells"; 8 →
rule 11, the 16b `variants.md` row and Q11; 9 → rule 7 and symmetry property (e).

### Lead decisions on Q13–Q20 (2026-10-24)
All eight take the stated reversible default: Q13 `readthrough_codon`/`recoded_codon`;
Q14 `sixframe` optional on `loci` with the WARNING event; Q15 no new ADR beyond the ADR-0058
amendment (§9.4 extends through the schema); Q16 `multi` and `low_complexity` listed in the schema
pattern, with a backlog line; Q17 rows in `detection_summary.tsv`, not a new file; Q18 the
precedence exact > tolerant > frameshift, and overlapping diff-type groups; Q19 `sec_site` and
`pyl_site` replace `readthrough` and stay placement flags only; Q20 the e2e floor ≥ 3 of 5. ADR-0058
records item 6 as accepted and Amendment 1.


### Statistician review of the amendment (2026-10-24, at db7d92ab)

**APPROVE WITH CHANGES**, all applied in this commit. Full review: `local/stat/m7-16b-spec-review.md`.

1. Rule 20: `truth/stop_loss.tsv` lists every stop-loss in the build (`planted` column), so a
   natural stop-loss does not fail the e2e's "unexplained row" check; the knob draws from a
   separate generator, picks only genes whose C-terminal peptide is ≥ 10 aa with no other strain
   event, and raises on a shortfall, so the ≥ 3 of 5 floor tests the path, not the digest.
2. The e2e's "ends on the codon" is stated per strand; a failure is a finding, not a reseed.
3. The full-tier observer asserts the stop-loss only for rows at their origin locus (implied by
   the model); off-origin rows are reported, not asserted (they are chance placements already
   in the FDP).
4. Labels: rule 17 and the `placement-fdr.md` row state that the flagged rows are not an
   FDR-controlled subset and that `value / n` is never presented as a rate.
5. Rule 11: the class is written `(substitution, mass_equivalent = false)`; a patched `X`
   database does not lift the degeneracy; the class rate does not feed `confidence` or TE2.

16b re-runs (no tolerance, bound, α or R change): the related-strain placement suite's full tier
with the observer (rule 20 acceptance); fast-tier regression runs of `test_variants_gate.py`,
`test_variant_class_error.py`, `test_novel_orfs_gate.py`, `test_novel_class_fdr_calibration.py`
and `test_detection_summary_calibration.py`, recording "unchanged" plus the count of
`cterm_extension` calls and of REF `*` ALT sites (expected 0 with `n_readthrough = 0`). The
M3-21 placement suites need no re-run: rules 15–19 do not touch placement.

### Lead decision on the rule-20 e2e genome size (2026-10-25, test-author finding)
The sim100k e2e cannot meet its floor: at 100 kb the cost-1 stratum holds T_1 = 27–45 targets, so
the smallest q is (D+1)/T_1 ≥ 0.022 > α 0.01 and no cost-1 row is accepted (measured over seeds 11
and 12, pi_wrong 0.05–0.3 and SNV rates 5–20 /kb). Overriding α would break R3. The e2e therefore
runs on a **1 Mb** gippasim genome (seed 11, `n_stop_loss` = 5, floor ≥ 3 planted rows). There T_1
is about 400, D = 0 and q = 0.0025; the chain takes about 30 s, so the test is not marked slow.
This replaces "sim100k" in rule 20 and its acceptance item. The full-tier observer count (D5) is
reconciled by the statistician (local/stat/m7-16-observer-reconcile.md, Slurm 1891223): the
new observer and the 16a counts agree row for row on the same seeds. Replicates 0-6 of the
`placement.sim_related` stream count 0 rows, so the methods-note numbers stand unchanged.
Simulator (lead, 2026-10-25): `test_planted_cterm_peptide_is_the_reference_c_terminal_product` is
stricter than rule 20 (the whole protein must differ only at the stop). The gippasim knob plants on
such "pure" genes first and falls back to genes that meet rule 20's peptide-level condition, so
both that test and the natural stop-loss test hold at 5–20 SNV/kb. Accepted as tooling.

### Lead decisions in 16b-3 (2026-10-25, implementer findings)
- **The sim100k read-through e2e uses `n_readthrough = 40`, not 10.** At 10 the `cterm_extension`
  class holds T = 12 targets, so its q floor c/T = 0.073 is above α 0.01 and nothing can be called
  (`novel_orfs.class_q_floor_above_alpha`). At 40, T = 96 and the floor is 0.0091. The floors
  (≥ 1 called, ≥ 3 ALT sites) and α are unchanged (R3); only the world grows.
- **`alt_codon_vcf` is the forward-strand allele**, so on `−` the expected codon is the reverse
  complement of the gene-strand codon (a test error, not a pipeline one).
- **A planted sORF is not a `cterm_extension` row.** Rule 12 merges its tolerant C-terminal hits
  with the sORF's own peptides into one ORF, and `CLASS_PRIORITY` ranks `sorf` above
  `cterm_extension`. The called-gene check skips sORF genes; whether a read-through sORF should
  get its own class is a backlog question.

### Lead decisions in the #335 review round (2026-10-25)
- **Rule 14 (`SimConfig.n_readthrough`, `truth/readthrough.tsv`, its sim100k e2e) is deferred to a
  follow-up PR, 16b-3, before 16c.** 16b ships the strain stop-loss knob (rule 20), which already
  drives the e2e; the strain e2e reads `truth/readthrough.tsv` only if it exists. Until 16b-3,
  `cterm_extension` calls are covered by unit tests only (statistician item). Backlog line added.
- **Rule 17 row order: deviation.** `detection_summary.tsv` keeps the M4-16 contract's sort
  (section, metric, level), so the 12 rows appear alphabetically, not in rule 17's listed order.
  The set of rows and their values are as rule 17 says; readers select by key.

### Lead decisions for 16b-3 (rule 14, 2026-10-25, test-author measurement)
- **The rule-14 e2e stays at sim100k, with de novo input** (`DenovoConfig()`, `peptides_denovo.tsv`).
  Database input at sim100k puts T_1 = 7 targets in the cost-1 stratum (q_min = 0.14 > α), the
  16b floor problem. De novo read errors fill stratum 1 (T_1 = 1361, q_min = 7.3e-4), so an
  error-free read over the planted stop is accepted; a read that also carries a de novo error is
  cost 2 (T_2 = 54, q_min = 0.019) and is not. α is unchanged (R3). `gippa run --modules core`
  takes about 35 s, so the test is not marked slow.
- **Readings pinned in `tests/helpers/readthrough_sim_oracle.py`, accepted:**
  - eligible genes are non-pseudo CDS and sORFs (rule 1's feature set);
  - "5–60 codons after" means the new stop is k ∈ [5, 60] codons past the annotated stop;
  - `ext_stop` is the new stop's lowest 0-based coordinate, wrapped;
  - the genome and annotation are not edited (genes are chosen from the existing sequence);
  - `n_true_peptides_spanning` counts distinct detected correct peptides covering the stop.
- **Floors:** at least 1 called `cterm_extension` gene and at least 3 `stop_codon` ALT sites. On a
  shortfall raise `n_readthrough`; never change the seed or drop an assertion.
