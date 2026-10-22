# M7-02: Pileup, calls and confidence at every covered codon

Milestone M7 · track V (variants) · **pass 1** · depends on: M7-01 · estimate ≈ 260 code lines
(+ ≈ 450 test lines) · not statistical (counting rules), but the F1 gate (M7-13) measures it

## Spec sections implemented

- §7.5 preamble: "At every genome codon covered by accepted peptides, GIAP reports whether the
  peptides support the genome's residue (REF), a different residue (ALT), or both. It never names
  modifications."
- §7.5.2 steps 2–4: "**Pileup.** For each covered codon, count distinct peptides supporting REF
  (exact residue identity), each ALT residue, and `IL` (the residues differ only by I versus L).
  **Call:** `REF`: only REF support; `ALT`: only ALT support; `mixed`: both; `no_call`: only `IL`
  support … `IL` support is ignored when REF or ALT support exists. **Confidence.** An ALT (or the
  ALT part of a mixed call) is `high` with ≥ 2 distinct peptides supporting the same ALT residue,
  and `candidate` with 1. Every supporting peptide needs placement q ≤ 0.01."
- §7.5.1 rows 1–2 (the exclusions that change which peptides may support an ALT): "Search against
  another strain … with `--genome-relation same`: no ALT calls from these peptides".
  The flags of §7.5.1/§7.5.3 are M7-05.
- §2.1 decision 15 ("from distinct-peptide counts"); config `variants.high_min_peptides` (2),
  `stats.alpha.variants` (0.01); §9.3 `tables/pileup.tsv.gz` and the site/call columns of
  `tables/variants.tsv`; §6.5 sort order.

## Public interface

`src/giap/variants/pileup.py` (pure)
```python
CallKind = Literal["REF", "ALT", "mixed", "no_call"]
AltType = Literal["substitution", "swap", "block", "frameshift"]

@dataclass(frozen=True, slots=True)
class Support:                    # per peptide, joined from C9 peptides.norm + peptide_hits
    peptide: str
    source: Literal["database", "denovo"]
    search_db: Literal["same", "other"]
    placement_q: float
    uniqueness: str

def eligible(aligns: Sequence[HitAlignment], support: Mapping[str, Support],
             cfg: GiapConfig) -> list[HitAlignment]: ...
def pileup(aligns: Sequence[HitAlignment], gm: GeneModel, genome: Mapping[str, str],
           cfg: GiapConfig) -> pd.DataFrame: ...          # one row per covered codon
def sites(pile: pd.DataFrame, cfg: GiapConfig) -> pd.DataFrame: ...   # variants.tsv core rows
```

## Rules

1. **Eligible placements.** A `peptide_hits` row supports anything only if `placement_q ≤
   stats.alpha.variants` and `uniqueness != shared` (decision 3). Under `genome_relation = same`,
   a `database` peptide with `search_db = other` keeps its `match`/`il` support but contributes
   **no ALT** (§7.5.1 row 2). Every other combination contributes fully (flags are M7-05).
2. **Distinct peptides.** Counting is over distinct peptide sequences **after I/L folding**
   (decision 1): two rows that differ only by I/L are one peptide; charges, runs and modified forms
   were already collapsed by preflight (C9 is one row per stripped peptide).
3. **Codon key.** `(contig, strand, codon)` with `codon` the lowest 0-based coordinate, wrapped
   into `[0, length)` on circular replicons (so the same codon reached by an unwrapped hit and a
   plain hit is one key).
4. **What a peptide supports at a codon.** `match` → REF; `il` → IL; `substitution` → ALT
   `(alt, substitution)`; `swap` and `block`/`near_block` → ALT `(ref>alt token, swap|block)` at
   `codons[0]` only, and **nothing** at the block's other codon(s) (decision 4). `near_block` is
   written as `alt_type = block` (the §9.3 enum has no `near_block`; the `mass_equivalent` flag of
   M7-05 tells them apart). A peptide counts once per codon, in one category.
5. **ALT residues and I/L.** A substitution ALT of `I` or `L` is recorded as `J` (I or L, the
   IUPAC–IUBMB ambiguity code), so peptides reading I and L at a V codon support **one** ALT
   (decision 2). A substitution where the genome residue is I or L and the peptide reads L or I is
   `il`, never ALT.
6. **Calls** exactly as §7.5.2 step 3: REF-only → `REF`; ALT-only (one or more ALT keys) → `ALT`;
   REF and ALT → `mixed`; IL-only → `no_call`; IL ignored otherwise.
7. **Confidence** per ALT key: `high` when its distinct peptides ≥ `variants.high_min_peptides`,
   else `candidate`.
8. **Frameshift events.** Keyed `(contig, strand, last_prefix_codon, shift)`. n_alt = distinct
   peptides with that event; n_ref = distinct eligible peptides whose ops cover both the event's
   last prefix codon and the next codon of the same frame without a junction. Call and confidence
   as rules 6–7. They are `variants.tsv` rows with `alt_type = frameshift`, never pileup rows.
9. **Gene and residue position.** `gene` = the C2 CDS, sORF or pseudogene on the same strand whose
   frame contains the codon (M3-14 `in_frame`, origin-aware); with several, the smallest locus tag.
   `aa_pos` = 1-based residue index from the annotated start: (codon − start)/3 + 1 on `+`,
   (end − codon − 3)/3 + 1 on `−`. No gene → both NA. `ref_aa` = the genome residue (table 11 on
   the gene's strand; an annotated non-ATG start translates as V or L here, not M).
10. **`pileup.tsv.gz` columns** (§9.3): `contig`, `codon_start`, `strand`, `gene`, `aa_pos`,
    `ref_aa`, `n_ref`, `n_il`, `alt_counts` (`X:n;Y:m`, keys sorted by count descending, then key),
    `call`. One row per covered codon, §6.5 order (contig, codon_start, strand).
11. **`sites` rows** (the core of `variants.tsv`, columns finished in M7-03/05/10): one row per
    (codon, ALT key) where the call is `ALT` or `mixed`, plus one per frameshift event:
    `variant_id` (`VAR_<n>` over the §6.5-sorted rows, ties by `alt_type`, then `alt_aa`),
    `contig`, `codon_start`, `strand`, `gene`, `aa_pos`, `ref_codon` (gene-strand codon as read),
    `ref_aa`, `alt_aa`, `alt_type`, `call`, `confidence`, `n_ref`, `n_alt`, `source` (sorted
    `;`-list of the ALT peptides' sources).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no accepted peptides | header-only pileup and sites |
| every peptide at a codon is `shared` | the codon is not covered (no pileup row) |
| peptide with placement q above `stats.alpha.variants` | ignored (it cannot be in `peptide_hits.tsv` at the default 0.01, but `stats.alpha.variants` may be set lower) |
| same peptide reaches one codon through two loci (unwrapped and wrapped) | counted once (rule 3) |
| I-peptide and L-peptide at an L codon | `n_il` = 1 (folded, rule 2), call `no_call` |
| two ALT residues, no REF | call `ALT`, two site rows, each with its own confidence |
| codon outside every gene | row with `gene` NA, `aa_pos` NA |
| codon crossing the origin | normal row; VCF handling is M7-03 |
| `search_db = other` peptide, `genome_relation = same`, tolerant hit | REF/IL support kept, ALT dropped; the codon may be `REF` |
| block peptide only | ALT at the first codon; the second codon has no row unless another peptide covers it |
| annotated GTG start covered by `M…` N-terminal peptides | `ref_aa` V, ALT `M` (the `start_codon` flag is M7-05) |

## Examples

Gene G1 `+` [100, 400) on circular `c1` (1,000 nt); codon 112 encodes `N` (aa_pos 5).

| # | Peptides covering codon 112 (all q ≤ 0.01, unique) | Expected pileup / sites |
|---|---|---|
| 1 | P1, P2 read `N` | `n_ref` 2, `alt_counts` empty, call **REF**; no site row |
| 2 | P1 reads `D` | call **ALT**, one site `N→D`, `candidate`, `n_alt` 1 |
| 3 | P1, P2 read `D` (different sequences) | **ALT**, `high`, `n_alt` 2 |
| 4 | P1 reads `D` with charge 2 and 3, two runs | one distinct peptide → `candidate` |
| 5 | P1 reads `N`, P2 and P3 read `D` | **mixed**, site `N→D` `high`, `n_ref` 1, `n_alt` 2 |
| 6 | codon encodes `L`; P1 reads `I`, P2 reads `L` | `n_ref` 1 (P2), `n_il` 1 (P1): call **REF** (IL ignored when REF exists) |
| 7 | codon encodes `L`; only P1 reads `I` | `n_il` 1, call **no_call** |
| 8 | codon encodes `V`; P1 reads `I`, P2 reads `L` | one ALT key `J`, `n_alt` 2, **high** |
| 9 | codon encodes `N`; P1 reads `D`, P2 reads `E` | call **ALT**, two sites, each `candidate`; `alt_counts` `D:1;E:1` |
| 10 | block: P1 reads `GG` at the `N` codon 112 | site at codon 112, `alt_type` block, `alt_aa` `GG`, `ref_aa` `N`; `alt_counts` `N>GG:1`; codon 115 untouched by P1 |
| 11 | P1 `database`, `search_db other`, `genome_relation same`, reads `D` | no ALT; codon 112 not covered by P1 at all (no REF at the differing codon); no site |
| 12 | `−` gene G2 [500, 800); codon 785 | `aa_pos` = (800 − 785 − 3)/3 + 1 = **5** |

## Acceptance criteria

- [ ] `tests/unit/variants/test_pileup.py`: examples 1–12 and every degenerate row; each §7.5.6
  unit item: distinct-peptide counting (charges and runs collapse), IL → `no_call`, `mixed`, the
  confidence thresholds (with `high_min_peptides` 2 and 3).
- [ ] Property: permuting the input rows leaves both frames byte-identical; adding a REF-only
  peptide never turns `REF` into anything else; adding an ALT peptide never lowers a confidence.
- [ ] Property: on a genome with no planted difference and peptides cut from its own genes, every
  covered codon is `REF` (or `no_call` only where the peptide's I/L differs from the genome).
- [ ] Sort order and IDs follow rule 11; `variant_id` is a pure function of the sorted rows.

## Calibration and attacks owed

The F1 of `high` calls on simulated variants is the M7 gate (M7-13). Attack 6 (§8.5, isobaric
swaps everywhere) tests these counting rules together with the M7-05 flags.

## Hardening owed (pass 2)

- Golden pileup and variants tables for the 50-kb fixture (§11.2 "Golden … variant tables"; M7-15).
- Re-examine decision 3 (shared peptides) against the pilot (M10), where IS copies and paralogs are
  real.

## Questions for the lead

1. Are two peptides that differ only by I/L one "distinct peptide"? **Proposed default:** yes.
   Mass spectrometry cannot tell them apart (the reason `IL` exists), so counting them as two would
   let one spectrum's I/L ambiguity make a `high` call.
2. How is an ALT of I or L written? **Proposed default:** `J` (the IUPAC–IUBMB code for "Leu or
   Ile"), so I- and L-reading peptides at a V codon support one ALT; the VCF ALT codon is the
   minimal change to any I or L codon (M7-03). The alternative (two ALTs, `I` and `L`) splits one
   observation into two `candidate` calls.
3. Do `shared` peptides (different genes, §7.2.6) support REF/ALT? **Proposed default:** no — they
   count toward neither detection nor abundance, and a peptide that fits two paralogs equally well
   at its best cost says nothing about which one carries the residue. `multi_identical` peptides
   (identical proteins) count at every copy.
4. A block or swap peptide at its block's second codon: no support (proposed), REF, or ALT?
   **Proposed default:** none. §7.5.2 records blocks "at the first codon of the block"; counting
   the second codon as REF would be false (its residue is not observed) and as ALT would double the
   event.
5. `variant_id` scheme: §6.5 defines none. **Proposed default:** `VAR_<n>` in §6.5 order (rule
   11). `giap compare` (M9) matches sites by coordinates, not by ID.
