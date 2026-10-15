# M3-14: Locus classes, uniqueness, gene assignment and genome relation

Milestone M3 · depends on: M3-01 (`GeneModel`), M3-02 (frame rule), M3-04 (hit schema), M3-10
(accepted placements), M3-12 (protein rows) · creates `src/giap/mapping/loci.py` · pure Python ·
estimate ≈ 380 code lines

## Spec sections implemented

- §7.2.6 "**Locus class** of each placement relative to Bakta features. When several apply, the
  first in this priority order wins: 1. `rna_overlap` (suspicious); 2. `pseudogene`; 3. `sorf`;
  4. `cds_inframe`; 5. `nterm_extension` (same frame, upstream of the annotated start, after the
  previous in-frame stop); 6. `cds_outframe`; 7. `antisense`; 8. `intergenic`. A peptide with loci
  of different classes also gets `multi`."
- "**Uniqueness:** `unique`: one locus; `multi_identical`: several loci in genes with identical
  protein sequences, such as IS copies; counts toward detection of each copy, flagged, never
  toward abundance; `shared`: different genes; counts toward neither detection nor abundance of a
  single gene."
- "**Peptides → genes.** `cds_inframe` and `sorf` placements assign a peptide to that gene.
  `nterm_extension` counts toward the gene's detection (flag `nterm_extension_evidence`) but not
  its abundance, and goes to Module 6."
- "**Proteins → genes** through §7.2.3. A cross-strain protein … carries its abundance to the
  strain gene it maps to. Two input proteins on one gene → flag `two_proteins_one_gene`."
- "**Consistency check.** With same-strain input, a peptide whose engine `protein_ids` disagree
  with the mapped gene is flagged `engine_protein_mismatch`."
- §7.2.7 "**`related`** (default): … Accepted tolerant hits count toward detection and abundance
  and are labelled `strain_difference`. **`same`**: … Tolerant hits are reported but only exact
  hits count." (the counting rule itself is applied by M4, §7.3.1).
- §9.4 placement flags `same_db`, `fs`, `st`, `engine_protein_mismatch`, `origin_spanning`; gene
  flags `two_proteins_one_gene`, `cross_strain_protein`, `strain_difference`,
  `nterm_extension_evidence`, `multi_identical`.

## Definitions (normative)

All tests run on unwrapped intervals; when a hit and a feature are on opposite sides of the origin
of a circular replicon, the hit is shifted by +L before comparing. "In frame" with gene `g`:
plus `(hit.start − g.start) % 3 == 0`; minus `(g.end − hit.end) % 3 == 0`. Frameshift hits use the
prefix block for the frame and the whole span for overlaps.

| class | rule |
|---|---|
| `rna_overlap` | overlaps (≥ 1 bp, either strand) an `rrna`, `trna`, `tmrna` or `ncrna` feature |
| `pseudogene` | overlaps a pseudogene on the same strand |
| `sorf` | inside a sORF, same strand, in frame |
| `cds_inframe` | inside a CDS, same strand, in frame |
| `nterm_extension` | same strand and frame as a CDS, starts upstream of its start (plus: `hit.start < g.start`; minus: `hit.end > g.end`), may reach into the CDS, and no in-frame stop codon lies between the hit and the annotated start |
| `cds_outframe` | overlaps a CDS or sORF on the same strand, not in frame |
| `antisense` | overlaps a CDS or sORF on the opposite strand |
| `intergenic` | none of the above |

## Public interface

```python
LocusClass = Literal["rna_overlap", "pseudogene", "sorf", "cds_inframe", "nterm_extension",
                     "cds_outframe", "antisense", "intergenic"]
Uniqueness = Literal["unique", "multi_identical", "shared"]

def classify_locus(iv: Interval, gm: GeneModel, genome: Mapping[str, str]) -> tuple[LocusClass, str | None]: ...
def annotate_peptide_hits(accepted: pd.DataFrame, gm: GeneModel, genome: Mapping[str, str],
                          peptides: pd.DataFrame, genome_relation: Literal["related", "same"]) -> pd.DataFrame: ...
#   adds: locus_class, gene, uniqueness, flags (+ `multi` per Q1), novel_candidate_id = NA
def annotate_protein_hits(rows: pd.DataFrame, gm: GeneModel) -> pd.DataFrame: ...   # flags
```

- `gene` is set for `cds_inframe`, `sorf` and `nterm_extension` (flag
  `nterm_extension_evidence`), and for `pseudogene` (the pseudogene's locus tag; Q3); NA otherwise.
- **Uniqueness** over the peptide's accepted target loci: one gene-or-locus → `unique` (several
  loci inside one gene count as one; Q2); several loci that all assign to genes with identical
  `aa` → `multi_identical` (flag `multi_identical`); otherwise `shared`.
- **Flags:** `same_db` (from M3-10); `origin_spanning` when the locus crosses the origin;
  `strain_difference` for accepted tolerant hits when `genome_relation = related`;
  `engine_protein_mismatch` when `search_db = same`, the gene is set, and no non-decoy,
  non-contaminant entry of `protein_ids` equals the gene's locus tag.
- **Proteins:** `cross_strain_protein` when `protein_id` is not a locus tag of the model;
  `two_proteins_one_gene` on every row of a gene assigned ≥ 2 distinct input proteins.

## Examples

Model on circular `chr1` (L = 10,000): G1 cds + `[100,400)`; G2 cds − `[1000,1300)`; S1 sorf +
`[2000,2075)`; P1 pseudogene + `[3000,3600)`; R1 rRNA + `[5000,6542)`; G5 cds + `[5100,5400)`
(inside R1, artificial); G0 cds + `[9900,10200)` (crosses the origin).

| # | Hit | Expected (class, gene) |
|---|---|---|
| 1 | `[130,160) +` | `cds_inframe`, G1 |
| 2 | `[131,161) +` | `cds_outframe`, NA |
| 3 | `[130,160) −` | `antisense`, NA |
| 4 | `[70,100) +`, no in-frame stop in `[70,100)` | `nterm_extension`, G1, flag `nterm_extension_evidence` |
| 5 | `[40,70) +` with `TAA` at `[76,79)` | `intergenic`, NA (a stop separates it from G1) |
| 6 | `[1270,1300) −`; `[1300,1330) −` (no in-frame stop) | `cds_inframe` G2; `nterm_extension` G2 |
| 7 | `[5130,5160) +` | `rna_overlap`, NA (priority over `cds_inframe`) |
| 8 | `[2003,2030) +` | `sorf`, S1 |
| 9 | `[3100,3130) +` | `pseudogene`, P1 (Q3) |
| 10 | `[20,50) +` | `cds_inframe`, G0 (shifted to `[10020,10050)`), flag none; `[9990,10020) +` → G0, flag `origin_spanning` |
| 11 | peptide with loci in G3 and G4 whose `aa` are identical | `multi_identical` on both rows |
| 12 | peptide with loci in G1 and intergenic | `shared`; loci classes differ → `multi` (Q1) |
| 13 | accepted tolerant hit, `related` / `same` | flag `strain_difference` / no flag |
| 14 | `search_db = same`, `protein_ids = GSIM_00010`, mapped to G1 (`GSIM_00005`) | `engine_protein_mismatch` |
| 15 | proteins `P0A7V8` and `GSIM_00005` both assigned to G1 | both rows `two_proteins_one_gene`; `P0A7V8` also `cross_strain_protein` |

## Edge cases

- `nterm_extension` scan for the in-frame stop walks codons from the hit to the annotated start on
  the gene's strand, across the origin on circular replicons; on a linear contig, reaching the
  contig end without a stop is fine.
- A hit inside two overlapping same-strand CDSs in different frames: `cds_inframe` for the one in
  frame (priority), gene = that CDS.
- Partial genes are classified like others (M4 excludes them from statistics).

## Acceptance criteria

- [ ] `tests/unit/test_loci.py`: examples 1–15 and the edge cases.
- [ ] Property: on giapsim genomes, every peptide cut in frame from inside a CDS classifies
  `cds_inframe` with that gene; the same window reverse-complemented classifies `antisense`.
- [ ] Property: classification is invariant to the order of features in the model.
- [ ] Benchmark: 500,000 hits on a 5-Mb model ≤ 20 s (interval lookup by sorted arrays).

## Questions for the lead

1. "also gets `multi`": a flag or a class value? `multi` is not in the §9.4 vocabulary.
   **Proposed default:** a placement flag `multi`, added to the flag vocabulary in the schema.
2. Several loci of one peptide inside one gene (internal repeat): **Proposed default:** `unique`
   (gene-level unambiguous).
3. Pseudogene placements: assign `gene` (the pseudogene's locus tag)? §7.2.6 names only
   `cds_inframe` and `sorf` (+ `nterm_extension` for detection), but pseudogenes are rows of
   `genes.tsv` with covariates (§7.3.3), and M4-06 Q3 proposes assigning them. **Proposed
   default:** assign (aligned with M4-06), any frame on the same strand (frameshifted pseudogenes
   have no single frame). If the lead says no, M4 drops pseudogenes from the eligible set.

## Lead decisions (test-author questions)
1. Deterministic output: rows sorted by (peptide, contig, start, strand); flags are `;`-joined in
   alphabetical order; input flags (`same_db`, `low_complexity`) are kept.
2. Timing budget (500k hits ≤ 20 s) lives in benchmarks/; the unit suite keeps at most a loose 10× guard
   (project rule since the flaky-timing fix).
3. Accepted: example 10 uses its own G0-only model (the spec's wrapping G0 would break examples 2 and 5);
   `multi` is a flag on every row of the peptide; one gene hit with two classes = `unique` + `multi`; two
   intergenic loci = `shared`; pseudogenes assigned in any frame on the same strand, opposite strand =
   `intergenic`; decoy/contaminant-prefixed protein_ids never clear `engine_protein_mismatch`; unknown
   genome_relation raises ValueError.
4. Novel-locus clustering is not part of M3-14 (Module 6); `novel_candidate_id` stays NA here.

## Lead decisions after code review (PR #50, changes requested)
1. In-frame same-strand overlap across a CDS/sORF 3′ end → `cterm_extension` (ADR-0013, accepted by Khoa),
   assigned to that gene; `same_strand` (→ `cds_outframe`) only when the hit is OUT of frame. Tests: the
   reviewer's geometries `[370,415)+` on G1, `[2060,2090)+` on the sORF, `[985,1015)−` on G2, and the
   origin-shifted `[180,210)+` over G0's end.
2. Loci are `multi_identical` only when both have a non-empty `aa`; pseudogenes (`aa=""`) on two unrelated
   genes → `shared`.
3. Tests pin every tie-break (nested same-frame CDSs, two overlapping pseudogenes, two genes sharing a start)
   and the reviewer's missing cases: stop immediately before the start (both strands), a hit 3 nt past the
   CDS end, an upstream hit past the end of a short gene, a stop codon split across the origin,
   `origin_spanning` with `end == L`.
4. Remove dead `_overlap_shift`; an unknown contig raises ValueError naming it.
