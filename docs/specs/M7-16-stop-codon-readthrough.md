# M7-16: Peptides that span a stop codon (Sec, read-through, `cterm_extension`)

Milestone M7 · track R (mapping) · **pass 2 · option B accepted by Khoa (2026-10-24, "B"); ADR-0058**
· depends on: merged M3-06/M3-07/M3-08 tolerant routes, M3-10 placement FDR, M3-14 loci, M3-21
placement suites, M3-23 makedb; M7-01, M7-05, M7-07, M7-09, M7-11b, M7-12, M7-18 · estimate
≈ 370 code lines in three PRs (+ ≈ 130 simulator lines, + ≈ 900 test lines) · **stats**: the
statistician reviews this spec and the `placement-fdr.md` addition **before any code** (§14.1);
reviewed 2026-10-24 at e7f7f744 (APPROVE WITH CHANGES), and the required changes 1–9 are applied
in the rules, test plan and acceptance criteria below

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
codons) is a proposed amendment that needs Khoa (Q1 below), and PR 16c holds it.

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
   "cells_target": n, "cells_decoy": {kind: n}}`. An empty gene model gives zero sites, an
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
      `(substitution, mass_equivalent)`. `mass_equivalent` is false for REF `*` (`_mass("*")` is
      None, and `left_block_read` sees the `*` in its pair). Its `class_error_rate` is therefore
      the pooled substitution rate. That rate is dominated by ordinary SNV sites and does not
      measure read-through reads, and with a reference search on Bakta proteins (which end at the
      stop) every read-through read is in the degenerate case of `variants.md` bias 1: the
      reference cannot contain the peptide. 16b states this in `variants.md` (Q11 asks Khoa
      about a separate class).
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
14. **Simulator** (giapsim, tooling). `SimConfig.n_readthrough: int = 0`. With the default 0, every
    existing scenario and seed is byte-identical. With `n > 0` the simulator picks `n` annotated
    genes (seeded through the generator's stream). Each must have its next in-frame stop 5–60 codons
    after the annotated stop, inside the replicon, wrapping on circular ones. The expressed protein
    reads the stop as `W` (TGA) or `Q` (TAA/TAG) and continues to the next stop. The annotation
    keeps the original stop. Truth goes to `truth/readthrough.tsv`: locus_tag, contig, strand,
    codon_start, codon, residue, ext_stop (0-based start of the new stop codon),
    `n_true_peptides_spanning`. In 16c, `n_selenoproteins: int = 0` plants an internal TGA read as
    `U` with `U` in the annotated protein.

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
| 2 | **Sec peptide over fdnG's UGA**, Nissle, `search_db = same`: `SLGMLAVDNQARVUHGPTVASLAPTFGR` (28 aa); JMJJOG_02511, `−` `[2576167, 2579215)`; frame window `SLGMLAVDNQARV*HGPTVASLAPTFGR` over `[2578585, 2578669) −`, the `*` = internal TGA `[2578627, 2578630)` | **16a/16b (ADR-0058 items 1–5): no target row**: the TGA is not the gene's last codon (`[2576167, 2576170)`), so it is not marked; the peptide stays `no_hit` exactly as on main. **16c (item 6, if Khoa accepts):** one row `[2578585, 2578669) −`, `n_diff` 1, `cds_inframe` JMJJOG_02511; variants: codon 2578627 ALT `U`, `candidate` (one distinct unique peptide), flags `stop_codon;vcf_unrepresentable`, no VCF record (would be POS 2578628) |
| 3 | **Non-annotated in-frame stop:** query `HTLDMRKQAAGLKE` (14 aa, k = 2) over `aa [23, 37)` = `HTLDMRK*AAGLKE` | no target row: p = 30 is not marked, so the `*` is never consumed; no other window within 2; `unmapped` `no_hit` |
| 4 | **Annotated stop + one more mismatch, over budget:** `NEFKWVPAGHA` (11 aa, k = 1) vs `aa [14, 25)` = `NEFK*VPAGHT` | no row (cost 2 > 1) |
| 4b | Same window, `NEFKWVPAGHT` (11 aa, k = 1) | one row `[42, 75) +`, `n_diff` 1 (the stop is the only allowed difference) |
| 4c | `NEFKWVPAGHAL` (12 aa, k = 2) vs `NEFK*VPAGHTL` | one row `[42, 78) +`, `n_diff` 2, `diff_types` `substitution;substitution`, stratum m = 2 |
| 5 | **Decoy symmetry:** `reverse` decoy of the 42-aa frame | mirrored cell 42 − 1 − 18 = 23 holds `*` and is marked; p = 30 → 11 is not. Query `LTHGAPVWKFENLQW` (the reverse of example 1): no target row; one decoy row at decoy `aa [16, 31)`, cost 1, mirrored target interval `[33, 78) +`. With the decoy cells removed (a deliberately broken build) that decoy row disappears: the test that pins symmetry. `shuffle` decoy: cell 18 marked, `*` at 18, 30, 41 in place |
| 6 | **`cterm_extension` past the annotated stop into the 3′ flank:** example 1's peptide plus `KWVPAGHTLDMR` (12 aa, `aa [17, 29)` = `K*VPAGHTLDMR`, cost 1) | both `cterm_extension` on G1. M7-06 ORF: G1's start to the next in-frame stop, `[0, 93) +` (stop TAA at p = 30). M7-09: `spans_annotated_stop` true, label `readthrough_candidate`, class FDR `unestimated`, 2 distinct peptides → `called` true. M7-07 `annotation_patch.gff3`: `c1 GIAP CDS 1 93 . + 0 ID=G1;locus_tag=G1;Note=cterm_extension` |
| 7 | **makedb patch with `X`:** `giap makedb --patch` on example 6's patch | record `>G1_patch …` with sequence `MSDKPLTEAGRWQLNEFKXVPAGHTLDMRK` (p = 18 → `X`, terminal `*` stripped); `db_id_map.tsv` row `G1_patch patch G1 … c1 1 93 +` |
| 7b | Patch row `Note=cterm_extension`, G1, `[0, 126) +` (runs over p = 30) | `ValueError`: "patch CDS G1_patch has an internal stop codon" (only the annotated stop becomes `X`) |
| 8 | **Minus strand:** contig `c2` linear, 300 nt; G2 CDS `−` `[120, 210)`, stop codon `[120, 123)` (gene-strand TAA) | site `(c2, −, 120)`; cell in frame `−0` (`ext_nt` 300, offset 0) at p = (300 − 0 − 3 − 120)/3 = 59; `to_interval(−0, 59, 60)` = `c2:[120, 123) −`. A read-through peptide's flank lies at coordinates < 120 |
| 9 | Second-pass peptide over the patched `X`: the engine reports `…NEFKXVPAG…` | placed as example 1 (any query letter over a marked cell costs 1); variants skip the X op (`ambiguous`, M7-18), so no ALT site; the ORF evidence (example 6) is unchanged |
| 10 | Nissle, ADR-0058 items 1–5, `sixframe_index` | event `mapping.stop_sites`: `terminal` 4,689, `not_a_stop` 3, `recoded` 0, `cells_target` 4,689 (all four replicons are circular in the recording; the 3 are partial genes at plasmid edges, JMJJOG_04967/04968/04975), equal decoy counts |

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

## Public interface and touch points

```python
# src/giap/mapping/sixframe.py
READTHROUGH: int = 22          # target-only code of a marked (consumable) stop; never in frames
def load_stop_cells(path: Path) -> StopCells: ...            # format 2; format 1 raises
# save_index(path, index, decoys, cells: StopCells | None = None)  -> writes format 2

# src/giap/mapping/readthrough.py   (new, pure)
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
| 16a | `src/giap/mapping/readthrough.py` (new) | the above (rules 1–3) |
| 16a | `src/giap/mapping/sixframe.py` | `READTHROUGH`; `save_index` format 2 + cells; `load_stop_cells` |
| 16a | `src/giap/mapping/tolerant_dp.py` | `_verify_core` (1, 1) case (rule 4); module docstring |
| 16a | `src/giap/mapping/tolerant.py` | `Search(frames, marked=None)`; `match_tolerant(..., cells=None)` |
| 16a | `src/giap/mapping/frameshift.py` | `_FSearch(..., marked)`; `nonfs_overlap` (1, 1) case; `match_frameshift(..., cells=None)` |
| 16a | `src/giap/steps/mapping.py` | `sixframe_index` input `genes`, event; `map_tolerant` passes cells |
| 16a | `workflow/subworkflows/local/mapping.nf` | `stepIn('sixframe_index', [g, pep, gm[2], gm[0]])` |
| 16a | `docs/methods/placement-fdr.md` | ADR-0058 section; required content in "Methods note (16a)" below |
| 16b | `src/giap/variants/align.py` | `_residues`, `_dp`, `align_hit`, `align_hits` (`stop_sites`) |
| 16b | `src/giap/variants/pileup.py` | rule 4c narrowed (U/TGA, O/TAG → ALT) |
| 16b | `src/giap/variants/vcf_records.py`, `flags.py` | `U → TGA`, `O → TAG` zero-change; no `nt_class`, `vcf_unrepresentable` |
| 16b | `src/giap/steps/variants.py` | build `site_keys(annotated_stops(...))`, pass to `align_hits` |
| 16b | `src/giap/novel/evidence.py` | `cterm_extension` rule (rule 12) |
| 16b | `src/giap/novel/outputs.py` | `patch_features` cterm rows |
| 16b | `src/giap/makedb/build.py` | patch `X` at the annotated stop (rule 13) |
| 16b | `tests/sim/giapsim/{config,features,proteomics,truth}.py` | `n_readthrough`, `truth/readthrough.tsv` (tooling) |
| 16b | `validation/proxies/te_proxies.py` | TE2 note text (internal Sec sites are not annotated stops; ADR-0058) |
| 16b | `docs/methods/variants.md`, `docs/methods/novel-orfs.md` | stop_codon reachability, U/O ALT, cterm calls; `variants.md` states that REF `*` and U/O zero-change ALTs pool in `(substitution, mass_equivalent)` and that their `class_error_rate` does not measure read-through reads (rule 11) |
| 16c | `src/giap/mapping/readthrough.py` | `annotated_stops(..., recoded=True)`: internal TGA with `U`, TAG with `O` in `Gene.aa` |
| 16c | `src/giap/steps/mapping.py`, `steps/variants.py` | pass `recoded=True`; event `recoded` count |
| 16c | `tests/sim/giapsim/…` | `n_selenoproteins` |
| 16c | `validation/proxies/te_proxies.py` | TE1 exception label `sec_site`; TE2 note |

Edits to **existing** tests are the lead's (rules.md). Expected ones: the step-registration test
(`sixframe_index` inputs) and mapping.nf nf-test channel wiring; `test_variants_align` tests that
pin "`*` raises" (it now raises only off the site set); M7-18 tests pinning U-over-TGA `rare`
(examples 7–8 of M7-18 change meaning); M7-09/M7-11b tests pinning "`cterm_extension` never
called"; M7-07 "no cterm patch row"; M3-23 "internal stop raises" for cterm rows; the TE2
note-text test.

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
   0.52–1.89 (reverse) and 0.58–1.74 (shuffle) on giapsim 1 Mb, and 0.81–1.15 (reverse) and
   1.14–1.42 (shuffle) on Nissle. The net effect on the total, (D_rt − T_rt)/T_regular, was at
   worst −9.4 × 10⁻⁴ (giapsim) and −1.3 × 10⁻⁴ (Nissle). Cite Slurm jobs 1889319 (Nissle) and
   1889320 (giapsim), and the stressed arm's measured values once it has run.
4. **The added mass is an approximation, not a bound:** ≈ G/(N(1 − s)) at m = 1 and
   ≈ 2G/(N(1 − s)) at m = 2, independent of L, scaling with gene density G/N. The measured realised
   shares were 0.04–0.08 % on Nissle and 0.05–0.20 % on giapsim (i.i.d. and tryptic laws).
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
  It must not import `giap.mapping.readthrough` or `tolerant_dp`, and the existing oracles are not
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
  `mirror_cells` draws no randomness (a spy on `giap.rng.for_step` and on `numpy.random` sees no
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
- **Statistical re-run (R3: no tolerance, bound, α or R change).** Run the M3-21 placement suites
  on the 16a branch, fast (`-m calibration`, `scripts/slurm/test_branch.sbatch`) and full
  (`calibration_full.sbatch` placement tasks, `GIAP_WORKTREE` + `PYTHONPATH` exported), in two
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
  change 2). On giapsim 1 Mb null genomes (no read-through planted) with the M3-19/M3-20
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
  shift (D_rt − T_rt)/T_regular of the replicate mean ≥ −0.01; (d) `recoded` cells if ADR-0058
  item 6 is accepted. The test-author does not write it. **16a's merge waits for it.**

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
  `(substitution, mass_equivalent)` and that their `class_error_rate` is the pooled substitution
  rate, which does not measure read-through reads (rule 11; Q11 default).
- [ ] **16c (only after Khoa accepts ADR-0058 item 6):** example 2's 16c column holds. On the
  Nissle TE2 proxy there are two `stop_codon` ALT `U` sites (JMJJOG_02511 `[2578627, 2578630) −`,
  JMJJOG_04733 `[4825395, 4825398) +`), each `candidate`, and `high_stop_codon` is 0 (the 16-mer
  and 27-mer Sec peptides are `shared`). TE1's 4 misses are labelled `sec_site`. The 16a
  statistical re-run is repeated for the 3 extra cells (fast tier suffices: 3 cells against
  4,689).
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
| | **16b total** | **≈ 115** (+ ≈ 90 simulator, + ≈ 350 test lines) |
| 16c | `readthrough.py` recoded source, steps, TE labels | **≈ 45** (+ ≈ 40 simulator, + ≈ 120 test lines) |

Total ≈ 370 code lines, against the gated estimate of 220. The difference is the frameshift
precedence (rule 7), the index format and loader, and the M7 call paths (16b), which the gated
estimate did not count.

**PR split.** 16a (matcher, decoys, index, mapping step, methods note, re-calibration) → 16b
(M7-01 re-alignment, variants U/O rule, cterm call and patch, `makedb` `X`, simulator planting,
e2e) → 16c (recoded Sec/Pyl codons, gated on Khoa). 16b and 16c both depend on 16a, and are
independent of each other except for the shared TE2 note (16c updates it after merging main). Each is under 400 code lines.

## Questions for the lead (each has a default; none blocks 16a)

1. **(Khoa) Annotated internal Sec/Pyl codons (ADR-0058 item 6).** Option B as worded does not
   reach fdnG, the stated motivation. **Default:** 16a/16b implement items 1–5 only. 16c is fully
   specified but not started until Khoa answers, and PROGRESS.md records "fdnG Sec peptides stay
   unmapped (internal TGA, ADR-0058 item 6 pending)". Batch it with the other Khoa questions.
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
7. **What DIA-NN reports for an `X` in a patched database.** **Default:** no special handling.
   Any letter over a marked cell costs 1, and X/B/Z are skipped in variants (example 9). Record
   it in the methods note, and check it on the first real second-pass search.
8. **Partial and pseudogene features.** **Default:** pseudogenes excluded; a partial CDS is
   included iff its last codon is a stop (rule 1).
9. **Statistical evidence.** **Answered by the statistician (2026-10-24): a stressed arm is
   required.** The added chance-hit mass is an approximation, ≈ G/(N(1 − s)) at m = 1 and
   ≈ 2G/(N(1 − s)) at m = 2 (measured 0.04–0.08 % on Nissle, up to 0.20 % on giapsim under the
   tryptic law), not a bound. The evidence is the exercised M3-21 re-runs, the stressed arm and
   the statistician's attack; the decoy-balance check may be vacuous (test plan).
10. **Index format.** **Default:** bump to 2; format-1 files raise (no silent fallback); Nextflow
    resume re-runs `sixframe_index` anyway because its inputs change.
11. **(Khoa) A separate `stop_codon` error-rate class** (review change 8). REF `*` sites and U/O
    zero-change ALTs pool into `(substitution, mass_equivalent)`, whose rate does not measure
    read-through reads. A third class dimension `stop_codon` would change §7.5.4 ("by alt_type and
    by mass_equivalent") and so needs an ADR. **Default (reversible):** pooled, and documented in
    `variants.md` (rule 11). Batch it with the other Khoa questions; it does not block 16a or 16b.
12. **Stressed-arm sizing.** The review sets the target (≥ 100 expected false target read-through
    rows per replicate) and suggests a build. Whether that build runs within the fast-tier time
    budget at `replicates_fast.placement` = 50 is not yet measured. **Default:** the test-author
    sizes the genome for the target and measures one replicate on Slurm. The arm runs in both
    tiers as the acceptance criteria say; if the fast tier is too slow, the test-author reports
    the timing to the lead and the statistician instead of shrinking the arm. No replicate count
    changes (R3).

## Statistician review (2026-10-24, at e7f7f744)

**APPROVE WITH CHANGES.** 16a code may start; the changes below are required before 16a merges.
Full review: `local/stat/m7-16-review.md`. Exploratory exact-mass runs: Slurm jobs 1889319 (Nissle)
and 1889320 (giapsim 1 Mb × 3), script `local/stat/m7-16/rt_mass.py`.

Findings:
- Decoy symmetry holds exactly in counts. The reverse decoy is exact for i.i.d. laws: the
  read-through D/T was 1 to 1e-16.
- The local context is mirrored, not kept: gene | `*` | flank becomes flank | `*` | gene. Under a
  tryptic law (C-terminal K/R), the read-through D/T per genome is 0.52–1.89 (reverse) and 0.58–1.74
  (shuffle). The net effect on the total chance mass is ≤ 10⁻³, so the list FDR is unaffected.
- D_rt/T_rt is not an FDR for the read-through subset.
- The ADR-0058 added-mass arithmetic is right, but it is an approximation, not a bound. It scales as
  ≈ G/(N(1 − s)) at m = 1 and ≈ 2G/(N(1 − s)) at m = 2. The realised share is 0.04–0.08 % on Nissle
  and 0.05–0.20 % on giapsim.
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
