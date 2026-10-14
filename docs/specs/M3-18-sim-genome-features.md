# M3-18: giapsim, genome features that mapping needs

Milestone M3 · depends on: M0 simulator skeleton · changes `tests/sim/giapsim/` (new
`features.py`, `strain.py`; `simulate.py` extended) · simulator tooling (implementer writes it,
test author tests it; AGENTS.md) · pure Python · estimate ≈ 380 lines (not counted by the
ADR-0007 gate, kept under 400 anyway)

## Spec sections implemented (§11.3)

- Replicons: "a circular chromosome (50 kb – 5 Mb) and optional plasmids".
- Non-coding features: "rRNA operons whose 16S ends in `…GATCACCTCCTTA`; tRNA genes with
  anticodons".
- Planted features (the subset Module 2 classifies): "pseudogenes (internal stop or frameshift);
  sORFs (< 30 aa); unannotated expressed ORFs (intergenic, out-of-frame, antisense); N-terminal
  extensions (annotated start set downstream of the true start)"; plasmids.
- "**Strain differences.** A "sample" genome differs from the "reference" genome by SNVs (ALT
  residues) and indels."
- "**Annotation.** Written as a Bakta-like JSON containing only the fields GIAP parses (validated
  by GIAP's own parser), plus GFF3. All truth tables are written too."
- Speed: "a 100-kb genome in < 1 s, a 5-Mb genome in < 30 s".
- Needed by M3-01 (types), M3-14 (all eight locus classes, origin crossing), M3-02/04 (linear and
  circular contigs), M3-07 (strain differences → tolerant hits). Promoter/RBS/terminator signals,
  prophage/BGC blocks, CRISPR and oriC are **not** in this issue (M2/M4 needs; see Q2).

## New `SimConfig` fields (defaults keep the M0 output byte-identical when all are off)

```python
n_plasmids: int = 0; plasmid_length: tuple[int, int] = (5_000, 60_000)
linear_contig_length: int = 0            # > 0 adds one linear contig (no wrap)
origin_crossing_gene: bool = False       # place one gene across the chromosome origin
n_rrna_operons: int = 0; n_trna: int = 0
n_pseudogenes: int = 0; n_sorfs: int = 0
n_nterm_extensions: int = 0              # annotated start moved 3·k nt downstream (k ∈ [5, 30])
n_novel_orfs: dict[str, int] = {}        # {"intergenic": n, "cds_outframe": n, "antisense": n}
strain: StrainConfig | None = None       # snv_rate_per_kb, indel_rate_per_kb, seed offset
```

## Outputs (added to `SimResult`)

- `bakta_like.json`: feature types `cds` (with `pseudo: true` for pseudogenes), `sorf`, `rRNA`,
  `tRNA`; `sequences` gains plasmids (`type: plasmid`, circular) and the linear contig
  (`topology: linear`); origin-crossing features written Bakta-style (`start > stop`).
- `annotation.gff3` (via M1-05's writer).
- `truth/features.tsv`: `feature_id`, `type`, `contig`, `start`, `end`, `strand`, `true_start`
  (N-terminal extensions), `pseudo_reason` (`internal_stop` | `frameshift`).
- `truth/novel_orfs.tsv`: `orf_id`, `class`, `contig`, `start`, `end`, `strand`, `aa`.
- With `strain`: `reference.fasta` + `reference_bakta_like.json` (the related genome: the
  sample genome with SNVs and indels applied; coordinates lifted), `truth/strain_differences.tsv`
  (`contig`, `pos_sample`, `pos_reference`, `ref`, `alt`, `type` snv|ins|del, `gene`, `aa_change`).

## Examples

| # | Config (seed 1) | Expected |
|---|---|---|
| 1 | M0 defaults | `genome.fasta`, `bakta_like.json`, `truth/genes.tsv` byte-identical to M0's output |
| 2 | `n_sorfs=3` | 3 `sorf` features, each < 30 aa, start codon to stop codon in frame |
| 3 | `n_pseudogenes=2` | 2 `cds` with `pseudo: true`; `pseudo_reason` set; an internal stop or a 1-nt indel present in `nt` |
| 4 | `n_rrna_operons=1` | a 16S feature whose sequence ends with `GATCACCTCCTTA` (on its strand), plus 23S and 5S |
| 5 | `n_nterm_extensions=2` | the annotated `aa` is a suffix of the translation from `true_start`; no in-frame stop between them |
| 6 | `origin_crossing_gene=True` | one gene with `start > stop` in JSON; M3-01 reads it as `end > length` |
| 7 | `linear_contig_length=20_000` | a second contig with `topology: linear`; no feature crosses its end |
| 8 | `n_novel_orfs={"antisense": 2}` | 2 ORFs on the opposite strand of a CDS, absent from the JSON, present in truth |
| 9 | `strain=StrainConfig(snv_rate_per_kb=1)` on 100 kb | ≈ 100 SNVs in `truth/strain_differences.tsv`; `reference.fasta` differs from `genome.fasta` exactly there |
| 10 | 5-Mb genome with every option on | runs in < 30 s; 100 kb in < 1 s |

## Acceptance criteria (test author)

- [ ] `tests/unit/test_sim_features.py`: examples 1–9.
- [ ] Every simulated annotation is read by M3-01's `gene_model_from_bakta_json` without issues
  ("validated by GIAP's own parser").
- [ ] Property: for random seeds, features never overlap in ways the config does not ask for (rRNA
  never overlaps a CDS unless planted), every non-pseudo CDS translates without internal stops.
- [ ] Timing test (`slow`): example 10.

## Questions for the lead

1. Which genome does GIAP receive in a strain-difference run? **Proposed default:** `genome.fasta`
   is the sample's own genome (peptides map exactly); `reference.fasta` is the related genome, so
   `--genome reference.fasta --annotation reference_bakta_like.json` exercises
   `genome_relation = related` and tolerant hits.
2. Should this issue also add the M2-facing signals (promoters, terminators, RBS, oriC, prophage
   and BGC blocks)? **Proposed default:** no; whichever milestone needs them first adds them.

## Lead decisions after code review (PR #25)
1. A lifted CDS is "clean" only with a valid start codon (ATG/GTG/TTG) as well as no frameshift, no stop
   gained and no stop lost.
2. A CDS that is not clean in the reference gets its reference `aa` from its own reference `nt`
   (translated to the first stop, table 11) and `ref_pseudo = true` in truth; no protein is ever listed
   that its `nt` does not encode.
3. Indels are left-normalized and written VCF-style: anchored `pos` (1-based, the base before the event),
   `ref` and `alt` both include the anchor base; `strain_differences.tsv` documents this convention.
4. `aa_change` is computed from the translated proteins (codon 1 reads as M) and accounts for in-frame
   indels in the same gene.
5. `truth/genes.tsv` covers every CDS including origin-crossing, plasmid, linear-replicon and planted ones;
   `features.tsv` `end` is the unwrapped end (> length for origin-crossing features), documented.
6. The origin-crossing gene carries Bakta's `edge: true`; sORFs stay `CDS` in GFF3 (as Bakta writes them)
   with an `sorf=true` attribute.
7. `SimConfig` stays hashable (tuple of pairs instead of a dict).
