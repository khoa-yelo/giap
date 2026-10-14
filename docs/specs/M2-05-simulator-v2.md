# M2-05: Simulator v2 (`giapsim`): 100-kb genome with features, signals and precomputed layers

Milestone M2 · depends on: M2-02 (precomputed formats), M2-03 (Bakta field names), **M3-18** (simulator
plasmids, rRNA/tRNA, pseudogenes, sORFs, origin-crossing genes; whichever merges second rebases) · estimate ≈ 400
code lines under `tests/sim/` (tooling; not gated by ADR-0007, kept small anyway) · capabilities: none

## Spec sections implemented

- §11.3 (subset needed for annotation): "a circular chromosome (50 kb – 5 Mb) and optional plasmids";
  "TU sizes geometric (mean 1.8 genes); intra-operon gaps −20 to +40 bp; inter-TU gaps 50–300 bp";
  "starts ATG 90 %, GTG 8 %, TTG 2 %"; "**Signals with tunable strength:** RBS motifs (AGGAGG at
  spacer 5–10); σ70 −35/−10 motifs upstream of TU starts; GC hairpin + U-tract terminators after TU
  ends"; "**Non-coding features:** rRNA operons whose 16S ends in `…GATCACCTCCTTA`; tRNA genes with
  anticodons; an oriC-like region next to a marker gene"; "**Planted features:** a prophage block with
  shifted GC; a BGC-like block; plasmids; a CRISPR array; pseudogenes (internal stop or frameshift);
  sORFs (< 30 aa)"; "**Annotation.** Written as a Bakta-like JSON containing only the fields GIAP
  parses (validated by GIAP's own parser), plus GFF3. All truth tables are written too."
- §11.3 outputs: "precomputed layers (`operon_pairs.tsv` with UniOP-like probabilities; promoters,
  Tx rates, terminators, RBS and regions, each with true elements plus false positives at known
  rates); `truth/*.tsv`; `config.yaml` with `annotation.precomputed` filled in"; "a 100-kb genome in < 1 s".
- §11.4: "The agent records them by running the tools on the 100-kb synthetic genome on FarmShare".

Out of scope here (M3+): peptides, abundance, strain differences, de novo errors, unannotated ORFs,
N-terminal extensions.

## Public interface

**Overlap with M3-18.** `docs/specs/M3-18-sim-genome-features.md` already adds plasmids
(`n_plasmids`, `plasmid_length`), `n_rrna_operons` (16S ending `…GATCACCTCCTTA`), `n_trna`,
`n_pseudogenes`, `n_sorfs`, `origin_crossing_gene`, the Bakta-style JSON (`start > stop` across the
origin) and `annotation.gff3`, in `tests/sim/giapsim/features.py`. This issue reuses those fields and
that module and adds only what annotation needs: signals, CRISPR, oriC, prophage/BGC blocks,
precomputed layers, annotation truth, `SIM100K`. New modules: `signals.py`, `blocks.py`, `layers.py`.

`tests/sim/giapsim/` (additions)
```python
@dataclass(frozen=True)
class SimConfig:                     # new fields (M0 and M3-18 fields unchanged); defaults = off
    real_rrna_trna: bool = False                   # use bundled E. coli rrnB / tRNA sequences (Q1)
    crispr: bool = False; oric: bool = False; prophage_bp: int = 0; bgc_bp: int = 0
    prophage_gc_shift: float = -0.08
    rbs_strength: float = 0.0; promoter_strength: float = 0.0; terminator_strength: float = 0.0
    layer_false_positive_rate: float = 0.05
    write_precomputed: bool = False

SIM100K = SimConfig(genome_length=100_000, n_plasmids=1, plasmid_length=(8_000, 8_000),
                    n_rrna_operons=1, n_trna=10, real_rrna_trna=True, crispr=True, oric=True,
                    prophage_bp=12_000, bgc_bp=10_000, n_pseudogenes=2, n_sorfs=3,
                    rbs_strength=0.8, promoter_strength=0.8, terminator_strength=0.8,
                    write_precomputed=True)
SIM100K_SEED = 20261014

@dataclass(frozen=True)
class SimResult:  # existing fields plus
    gff3: Path; config_yaml: Path; precomputed: Mapping[str, Path]; truth: Mapping[str, Path]
```
Outputs under `outdir`: `genome.fasta`, `bakta_like.json` (Bakta 1.12.2 field names per M2-03:
`sequences[]` with `id, type, topology, length, nt`; `features[]` with `type, sequence, start, stop,
strand, locus, gene, product, start_type, nt, aa, db_xrefs`, `pseudo`, `edge`; tRNA `anti_codon,
amino_acid`; `version: {bakta: "sim", db: {version: "sim", type: "sim"}}`), `annotation.gff3`
(M1-05 writer), `precomputed/{operon_pairs,promoters,tx_rates,terminators,rbs,regions}.tsv`
(M2-02 formats), `truth/{genes,tus,promoters,terminators,rbs,regions,features}.tsv`, `config.yaml`.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `simulate(SIM100K, SIM100K_SEED, d)` | `genome.fasta` has 2 records: `chr1` (100,000 bp, `[topology=circular]`), `p1` (8,000 bp, circular); runs in < 1 s |
| 2 | same | `parse_bakta_json(d/"bakta_like.json")` succeeds; gene count = rows of `truth/genes.tsv` |
| 3 | same | the 16S feature's `nt` ends with `GATCACCTCCTTA`; M2-04 `anti_sd` → `("ACCTCCTTA", ())` |
| 4 | same | `truth/features.tsv` has ≥ 10 tRNA rows with anticodons; one `crispr` row; one `oric` row adjacent to a gene named `mnmG` |
| 5 | same | every TU in `truth/tus.tsv` with `promoter=true` has `TTGACA` at −35 and `TATAAT` at −10 (17-bp spacer) of its TSS, upstream of the first start codon |
| 6 | `rbs_strength=1.0` | every gene has `AGGAGG` ending 5–10 nt before its start codon (strand-aware) |
| 7 | same | GC content of the prophage block is 0.08 ± 0.02 below the rest of `chr1` |
| 8 | same | each pseudogene has an internal in-frame stop or a frameshift; each sORF is < 30 aa |
| 9 | same | each `precomputed/*.tsv` loads with `load_precomputed` (M2-02); `config.yaml` sets every `annotation.precomputed.*` key to these files |
| 10 | two runs with the same seed | byte-identical directories |

## Edge cases

- Bundled real sequences (`data/ecoli_rrnB.fa`, `data/ecoli_trna.fa`, from NCBI K-12 MG1655,
  < 20 kb total, with a `PROVENANCE.md`) are inserted in an inter-TU gap so that Bakta's own rRNA/tRNA
  search finds them in recordings; no gene overlaps them.
- The prophage and BGC blocks are made of genes from the same codon model (GC shifted for the
  prophage); they carry truth labels only. Real tools may not call them (see Q2).
- Precomputed false positives (`layer_false_positive_rate`) are drawn with the simulator's RNG and
  listed in truth with `true_element=false`.
- `p1` gets 1–3 genes and is marked `type: plasmid` in the JSON.

## Acceptance criteria

- [ ] `tests/unit/test_sim_v2.py` (test-author; AGENTS.md: "the implementer writes it, the test-author tests it"): examples 1–10.
- [ ] Existing M0 and M3-18 simulator tests still pass (all new fields default to off, so the M0
  output stays byte-identical).
- [ ] `SIM100K` + `SIM100K_SEED` generate the genome that M2-10 records; its `genome.fasta.gz` is
  committed by M2-10 with this seed in `PROVENANCE.md`.

## Questions for the lead

1. Bundle real E. coli rRNA/tRNA sequences in the simulator? **Proposed default:** yes (NCBI
   sequences carry no restrictions); without them Bakta finds no rRNA/tRNA and the anti-SD and tAI
   contract paths are untested on real output.
2. Non-empty geNOMAD/antiSMASH fixtures cannot come from synthetic blocks. **Proposed default:** the
   sim100k recordings may be empty for those tools; non-empty contract fixtures come from K-12 and
   Nissle (M2-10).

## Lead decisions after code review (PR #38, changes requested)
1. Truth bug: `signals.get()` on the − strand returns the reverse complement in reading order (complement
   AND reverse is wrong for in-place checks; use the strand's reading order consistently with `put`), so
   the "no second AGGAGG within 5–10 nt" clean-up runs on − strand genes too. Test: over many seeds, no
   − strand (or + strand) gene has a second Shine–Dalgarno core within the spacer window besides the true RBS.
2. `truth/features.tsv` stays byte-identical for M0/M3-18 configs: the `anticodon` column is written only
   when v2 output is on (`write_precomputed` or `real_rrna_trna`). Add `truth/features.tsv` to the M0 and
   M3-18 hash tests.
3. Layer values match what the parse steps will write (the layers are what M2-12…M2-20 tests compare
   against): region ids `<kind>_<n>` (e.g. `prophage_1`, `bgc_1`) as in docs/specs/M2-19-regions.md;
   tool/source names with the official casing of docs/specs/M2-15-terminators.md and M2-18/M2-19
   (`TerminatorNet`, `TransTermHP`, `geNOMAD`, `antiSMASH`, `Bakta`); BGC score NA; false-positive regions
   list the genes they overlap like true ones. Where the M2-12…M2-20 specs define a value, the simulator
   follows them; the implementer reads those specs.
4. Tests kill the surviving mutants: terminator truth start ±1, U-tract change, terminator offset range,
   `tus.tsv` start ±1, region membership boundary, TSS +1 (each checked against the genome sequence).
5. Distances in `layers.py` exponentials are clipped (no overflow warning on gaps > 10 kb).
6. Accepted implementer choices (record in docstrings): plasmid id `plasmid_1` (tests accept `p1`); intra-TU
   gaps of 8–10 bp become 11–13 (no RBS spacer fits); no RBS planted for N-terminal-extension genes;
   signals are planted without truth tables when `write_precomputed=False`.
