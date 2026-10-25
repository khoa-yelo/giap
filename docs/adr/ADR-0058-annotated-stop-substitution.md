# ADR-0058: The tolerant matcher may consume an annotated stop codon (option B of M7-16)
- Status: accepted-by-khoa (2026-10-24, in session: "B"). Read-through residue in `makedb` patches
  = `X` (the proposed default, taken with B). Item 6 below (annotated internal Sec/Pyl codons) is
  **not** covered by that answer; it is a proposed amendment that needs Khoa.
- Date: 2026-10-24
- Spec sections: §7.2.1 (six-frame index, decoys), §7.2.4 (route 3, the cost table, G4 completeness),
  §7.2.5 (placement FDR, G3), §7.2.6 (locus classes), §7.5.1, §7.5.3 G8 `stop_codon`, §7.6.3,
  §7.9.2 (`makedb --patch`), §11.6 TE1/TE2; ADR-0013 (`cterm_extension`, `readthrough_candidate`)
- Needs Khoa: yes. It changes Module 2's search space and therefore what the placement FDR controls.
  Answered for items 1–5; open for item 6.

## Context
Three accepted parts of the design assume a peptide can be placed across a genome stop codon:
§7.5.3 G8 `stop_codon` ("the REF residue is a stop (e.g., Sec read through UGA)") and §7.5.1 ("Sec
sites … appear as ALT calls"); ADR-0013 (`cterm_extension` evidence is "a peptide spanning the
annotated stop", the `readthrough_candidate` label, and a patch row with a new end); and TE2's
`stop_codon` exception (§11.6).

The merged matcher cannot do any of this. The index writes `*` at every stop (§7.2.1). The exact
automaton never matches across one. The M3-06 DP says "No operation may consume a target `*`"
(`tolerant_dp._verify_core` gates every move on `t1 < 19`), and the frameshift route and the
M7-01 re-alignment (`variants.align._residues` raises on a genome `*`) do the same. So in v1 as
built, `stop_codon` and `cterm_extension` calls are unreachable (M7-05 rule 6, M7-09 rule 4,
M7-11b), and M7-18 lead decision 1 handed the Sec/Pyl ALT rule to this ADR.

M7-16 listed three options (A leave it, B annotated-stop substitution, C any-stop substitution).
Khoa chose B.

**What B reaches, checked on the Nissle recording** (`local/recordings/bakta_full/nissle1917`,
Bakta 1.12.2, full DB). There are 4,689 CDS and sORF features. Each one's last codon is a stop
(TAA 2,943, TGA 1,401, TAG 345). No two features share a stop codon, and none crosses an origin.
The 4 TE1 misses (M3-22) are the selenocysteine peptides of the two fdnG copies (JMJJOG_02511 −,
JMJJOG_04733 +): `VUHGPTVASLAPTFGR`, `VUHGPTVASLAPTFGRGAMTNHWVDIK`,
`SLGMLAVDNQARVUHGPTVASLAPTFGR` and `ALGMLAVDNQARVUHGPTVASLAPTFGR`. Bakta annotates both copies as
full-length selenoproteins (1,015 and 1,016 aa), with `U` at residue 196 over an **internal** TGA:
[2578627, 2578630) − and [4825395, 4825398) +. That TGA is not the last codon of the feature. So
**the option-B rule as worded does not reach the fdnG peptides.** It reaches C-terminal read-through
(`cterm_extension`) and any residue read at an annotated terminal stop. Item 6 states what reaching
fdnG would take.

## Decision
1. **Consumable stops.** The tolerant DP (route 3, M3-06/M3-07) may consume a target `*` **only**
   at the annotated stop codon of a CDS or sORF. That is the last codon of the feature
   (`Gene.interval`, `feature_type` in {`cds`, `sorf`}), on the feature's strand, in the feature's
   frame. It may consume it only by a **substitution** (one query residue against the `*`), at
   **cost 1**, recorded as `diff_types` `substitution` with REF `*`. The rule does not depend on
   the codon (TAA, TAG and TGA alike). A match is never possible: no query code equals a stop.
   Swaps, blocks and near-blocks never consume a stop, consumable or not. Every other `*` and every
   `X` stays unconsumable. Pseudogenes are excluded. So is a feature whose last codon is not a stop
   (for example a partial gene at a contig edge); it is counted in an event. On circular replicons
   the codon is marked at every copy in the origin-extended frames, including a codon that crosses
   the origin.
2. **Decoys.** Each decoy search carries the same rule at the **mirrored** positions. For a
   `reverse` decoy frame of `n` residues, target position `p` becomes `n − 1 − p` (`mirror_aa`).
   For a `shuffle` decoy it stays `p`, because stops keep their positions. Both decoys keep a `*`
   at the mirrored cell, so a target and its decoy expose the same number of consumable cells per
   frame. Decoy hits stay counted only, never placed.
3. **Routes.** The exact route is unchanged. The frameshift route (M3-08) never consumes a stop in
   either side. Its precedence check ("a non-frameshift alignment at cost ≤ the split's overlaps the
   split locus") does see consumable-stop alignments, so a split never coexists with a cheaper or
   equal read-through alignment of the same locus. M7-01 re-alignment uses the same set and rule.
   The recomputed cost must still equal `n_diff`.
4. **Variants.** A substitution op with REF `*` is an ALT site with flag `stop_codon` (G8; M7-05
   rule 6 as written). The `snv`/`multi_nt` class comes from the codon (for example W over TGA →
   `snv`). `U` over TGA and `O` over TAG are the standard decodings (Sec, Pyl). Per §7.5.1 they are
   **ALT** calls with `stop_codon`, zero nucleotide changes, neither `snv` nor `multi_nt`, no VCF
   record (ALT codon = REF codon) and flag `vcf_unrepresentable`. This replaces the M7-18 `rare`
   skip for exactly those two cases. Every other U/O, and every B/J/X/Z, is skipped as M7-18 does
   today.
5. **`makedb --patch` residue: `X`.** A patch CDS with `Note=cterm_extension` translates the
   annotated stop codon of its `locus_tag` gene as `X`. Search engines treat X as any residue, which
   covers Sec (U), Trp (W) and Gln/Tyr read-through. Any other internal stop in any patch CDS still
   raises, as M3-23 does.
6. **Proposed amendment (needs Khoa; not in force).** The consumable set would also include
   **annotated recoded codons**: an internal in-frame codon of a CDS whose Bakta protein has `U`
   over TGA or `O` over TAG at that residue. That is 3 sites in Nissle (fdnG ×2, fdhF). Decoys
   would be mirrored as in item 2. This is what makes the fdnG TE1 misses placeable. Until Khoa
   answers, the code builds only item 1's set (M7-16 PR 16c holds the change).

## Evidence
- **Search-space growth (analytic, Appendix B setting: N = 10⁷ six-frame positions, s = 0.069,
  G = 5,000 consumable stops).** A window can consume a marked stop only if it contains it: L
  windows per stop, with the other L − 1 residues within m − 1 differences. The added chance-hit
  mass is G·L·Σ_{j<m} C(L−1, j)(1−s)^j s^(L−1−j), against E(L, m) of §7.2.5:

  | L, m | E(L, m) | added by B | ratio |
  |---|---|---|---|
  | 10, 1 | 3.3 × 10⁻³ | 1.8 × 10⁻⁶ | 0.05 % |
  | 11, 1 | 2.5 × 10⁻⁴ | 1.4 × 10⁻⁷ | 0.05 % |
  | 12, 2 | 1.4 × 10⁻³ | 1.5 × 10⁻⁶ | 0.11 % |
  | 16, 2 | 5.8 × 10⁻⁸ | 6.2 × 10⁻¹¹ | 0.11 % |

  Option C, with about 5 × 10⁵ stops in 10⁷ positions, would add about 100× this, and that mass
  would sit at every stop of every frame.
- **Symmetry.** By construction (item 2), each decoy frame has exactly the target's number of
  consumable cells. Each such cell holds `*` and is mirrored one to one. The tests in M7-16 make
  this a property test.
- **Calibration.** Before/after numbers come from the re-runs that M7-16 requires (M3-21 placement
  suites fast + full, TE1/TE2/TE8). This section is filled in with the Slurm job IDs when 16a
  merges.

## Consequences
- **Module 2 search space.** Route 3 grows by at most one consumable cell per annotated gene and
  search (≈ 4,700 on Nissle), at ≤ 0.11 % of the chance-hit mass (above). Peptides of 7–9 aa
  never reach route 3 (`min_len` 10), so a stop-spanning peptide shorter than 10 aa is never
  placed. At 10–11 aa (k = 1) the stop is the only allowed difference. With no annotated genes the
  matcher is byte-identical to main.
- **Placement FDR (G3).** The estimator, the strata, π̂ (m ≥ 1 → 1) and the ranking score E are
  unchanged. Read-through placements are tolerant hits at m ≥ 1, so `same_db` exact skipping never
  applies to them. Validity rests on item 2. No tolerance, bound, α or R changes (R3).
- **Re-calibration required.** The M3-21 placement suites must re-run, fast and full, with no
  tolerance change. A failure stops the merge and goes to a new ADR, never to a tolerance change.
  This covers `test_placement_fdr_sim.py`, `test_placement_fdr_genome.py`,
  `test_tolerant_recall.py`, `test_denovo_precision_sim.py`, `attacks/test_attack_placement_*`
  and attack 5. The genome-level M3-21 suites exercise few m ≥ 1 lists at 100 kb to 1 Mb, so
  M7-16 adds a decoy-balance check over the consumable cells for the statistician. TE8 (Nissle
  pooled FDR ≤ 0.01) re-runs weekly as usual.
- **TE1.** Unchanged under items 1–5. The fdnG Sec site is internal, so the 4 misses stay unmapped.
  Under item 6 they would become tolerant placements (cost 1). They would still not be exact, so
  they would stay TE1 exceptions, now labelled as Sec.
- **TE2 and G8.** The `stop_codon` flag becomes reachable. TE2's second exception is exercised
  only when an annotated stop or (item 6) a recoded codon carries a placed peptide. On Nissle
  same-strain peptides, items 1–5 add no site. Item 6 would add two `candidate` U sites at the
  fdnG codons, one unique peptide each (the 16-mer and 27-mer are `shared` between the copies),
  and no `high` site.
- **`cterm_extension`.** It is reachable through route 3: a peptide in frame with a CDS/sORF that
  crosses its annotated stop into the 3′ flank (M3-14 `_cterm`, unchanged). M7-09's class rule
  becomes live (a spanning peptide plus an in-frame continuation; `readthrough_candidate`). Under
  M7-11 decision 2 the class stays `unestimated`, so it is called only with ≥ 2 distinct
  peptides. M7-07 writes its patch row, and `makedb --patch` translates the old stop as `X`.
- **Interfaces.** `sixframe_index` gains the `genes` input and writes index format 2, which
  carries the consumable cells for the target and each decoy type. No config key, no new column,
  no new diff type, no new flag value (`stop_codon`, `vcf_unrepresentable` and
  `readthrough_candidate` already exist).
- **Backlog.** "M7-16: decide `stop_codon` and `mass_equivalent` for blocks and swaps whose REF
  contains `*`" closes: no block or swap ever consumes a stop.

## Alternatives
- **A. Leave it.** No cost. But G8 `stop_codon`, ADR-0013's evidence and patch row, and TE2's
  second exception stay dormant or vacuous, and PROGRESS.md must keep saying so. Rejected by
  Khoa.
- **C. Any-stop substitution.** About 150 lines, but every stop in every frame becomes a cost-1
  chance-hit site: roughly 100× the added mass above, concentrated next to stops. The §7.2.5 /
  Appendix B.2 chance-hit table and the placement-FDR calibration would need re-derivation.
  Every pseudogene and frame-shifted ORF would also become a "read-through". Rejected.
