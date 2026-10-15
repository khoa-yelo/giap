# M4-07: giapsim layout: sequence-free genomes and contract-table fixtures

Milestone M4 · track F (simulator) · depends on: M4-06 (column names; the test compares against
`TABLES`) · estimate ≈ 380 lines in `tests/sim/giapsim/` (tooling, written by the implementer;
§AGENTS "the implementer writes it, the test-author tests it") · **can start now**

Calibration suites need hundreds of replicate genomes of 1,000–5,000 genes in seconds, and M4 must be
developed before M2/M3 exist. This issue adds a fast, sequence-free layout generator (amino-acid
sequences only, no nucleotide genome) and a writer of the M4-06 contract tables, plus an adapter that
turns the M0 full simulator output into the same tables.

## Spec sections implemented

- §11.3 Genome: "Replicons: a circular chromosome (50 kb – 5 Mb) and optional plasmids." "TU sizes
  geometric (mean 1.8 genes) …; strand runs with set convergent and divergent frequencies; log-normal
  lengths (median about 280 aa)". "Planted features: a prophage block with shifted GC; a BGC-like
  block; plasmids; a CRISPR array; pseudogenes (internal stop or frameshift); sORFs (< 30 aa)".
  "precomputed layers (`operon_pairs.tsv` with UniOP-like probabilities; … regions …) each with true
  elements plus false positives at known rates; `truth/*.tsv`".
- §11.3: "It is fast: a 100-kb genome in < 1 s, a 5-Mb genome in < 30 s."
- §6.5 IDs: regions `<type>_<n>` numbered by contig order then start; TUs `TU_<n>`.
- §7.1.7 prophage function subsets: `integrase`, `repressor`, `structural`, `lysis`, `cargo`.

## Interface (`tests/sim/giapsim/layout.py`, `contract_tables.py`)

```python
@dataclass(frozen=True)
class LayoutConfig:
    n_genes_chromosome: int = 1000
    plasmids: tuple[int, ...] = ()                 # genes per plasmid (circular)
    linear_contigs: tuple[int, ...] = ()           # genes per linear contig (draft assemblies)
    mean_tu_size: float = 1.8
    p_same_strand: float = 0.6
    gene_len_median_aa: float = 280.0; gene_len_sigma: float = 0.45
    frac_sorf: float = 0.02; frac_pseudogene: float = 0.01; frac_hypothetical: float = 0.25
    regions: tuple[RegionSpec, ...] = DEFAULT_REGIONS   # prophage 40 genes, bgc 25, crispr 6, amr 10 scattered
    uniop_beta_in: tuple[float, float] = (8, 2)          # p for true intra-TU pairs
    uniop_beta_out: tuple[float, float] = (2, 8)         # p for adjacent same-strand pairs across TUs
    aa_composition: Mapping[str, float] = ECOLI_AA       # bundled
    short_genes_in: str | None = None                    # attack null 2: region type whose genes are drawn short

@dataclass(frozen=True)
class Layout:
    replicons: pd.DataFrame; genes: pd.DataFrame; tus: pd.DataFrame; operon_pairs: pd.DataFrame
    regions: pd.DataFrame; region_genes: pd.DataFrame; truth: pd.DataFrame   # truth: tu_id, region effects slots

def simulate_layout(cfg: LayoutConfig, rng: np.random.Generator) -> Layout: ...
def layout_from_simulation(sim: SimResult) -> Layout: ...          # adapter for the M0 full simulator
def write_contract_tables(layout: Layout, work_dir: Path) -> None: ...   # C1–C8 (M2-side files) under work_dir
CONTRACT_COLUMNS: Mapping[str, tuple[str, ...]]                    # literal copy; a test compares to giap TABLES
```
- Gene coordinates are synthetic but consistent (3·(len_aa+1) bp, 50–300 bp gaps between TUs,
  −20…+40 within), so `contig/start/end/strand` satisfy the M4-06 invariants; `nt` is a
  back-translation with the bundled codon table (sORFs and pseudogenes included; a pseudogene gets an
  internal stop at a random codon and `aa = NA` with probability 0.5).
- Linear contigs: the first and last gene of each linear contig are `partial = true`.
- `operon_pairs`: every adjacent same-strand pair, intra-TU p ~ Beta(8, 2), inter-TU p ~ Beta(2, 8);
  TUs in the truth are the planted ones (`tus` = planted TUs, not re-derived from p).
- Regions: blocks of consecutive genes at random positions on the chromosome (non-overlapping), AMR
  as a scattered gene set, one `plasmid` region per plasmid listing all its genes; prophage genes get
  subsets with probabilities (0.05, 0.05, 0.4, 0.1, 0.4). `region_genes.core_biosynthetic` true for
  30 % of BGC genes.
- `layers.parquet` (C6): `tu_id`, `tu_position`, `tu_size`, `uniop_p_prev/next` (transcription order,
  as M2-20 defines), `region_ids`, `flags` (`partial`, `crosses_origin`, `pseudogene`); every other
  M2-20 layer column all-NA. `xrefs.parquet` (C3): random `COG####`, `GO:…`, `K#####` for 70 % of genes.
- `layout_from_simulation`: reads `bakta_like.json` and `truth/genes.tsv` of `simulate()`; TUs from
  `tu_id`; no regions (empty tables) until M2 extends the full simulator.

## Examples

| # | Call | Expected |
|---|---|---|
| 1 | `simulate_layout(LayoutConfig(), default_rng(1))` | 1,000 genes on `chr1` (circular); every gene in exactly one TU; mean TU size within 1.8 ± 0.2 |
| 2 | example 1 → `write_contract_tables` → `giap.io.work_tables.load_m4_inputs` | no `error` issues (M4-06 contract) |
| 3 | `LayoutConfig(plasmids=(60,))` | replicon `p1` circular, kind `plasmid`, 60 genes; region `plasmid_1` with 60 members |
| 4 | `LayoutConfig(linear_contigs=(30, 30))` | contigs `ctg1`, `ctg2` linear; 4 genes `partial = true` |
| 5 | default regions | region IDs `prophage_1`, `bgc_1`, `crispr_1`, `amr_1`; `prophage_1` has 40 consecutive gene indices; every prophage member has a non-NA `subset` |
| 6 | same config and seed twice | byte-identical Parquet files (TE7 discipline) |
| 7 | `LayoutConfig(n_genes_chromosome=4500)` | completes in < 2 s (simulator speed budget) |
| 8 | `layout_from_simulation(simulate(SimConfig(), 1, tmp))` | genes and TUs match `truth/genes.tsv`; converted coordinates are 0-based (`start = start1 − 1`) |
| 9 | `LayoutConfig(short_genes_in="prophage")` | prophage genes' median length ≤ 0.5 × the median elsewhere (attack null 2 input) |

## Edge cases

- Regions that cannot fit (more region genes than chromosome genes) → `ValueError`.
- `n_genes_chromosome < 2` → `ValueError`.
- A plasmid of 1 gene: one TU, no operon pairs.
- Locus tags `GSIM_<5·i:05d>` across replicons in contig order (as M0), unique.

## Acceptance criteria

- [ ] `tests/sim/test_layout.py` (test-author): examples 1–9.
- [ ] `CONTRACT_COLUMNS == {k: tuple(c.name for c in v) for k, v in giap.genetable.contract.TABLES.items()}`
  for the tables written (drift check, M4-06).
- [ ] `giapsim` imports nothing from `giap` (independence of truth and code under test); grep test.
- [ ] Speed: 1,000 genes < 0.2 s, 4,500 genes < 2 s (`benchmarks/`).

## Methods note

None (tooling). Documented in `tests/sim/README.md` section "layout" (the implementer may add it).

## Questions for the lead

1. A sequence-free generator separate from the M0 nucleotide simulator? **Proposed default:** yes:
   calibration needs 1,000 replicates per suite; the full simulator stays the source for integration
   and golden fixtures, and M2 extends it (genome signals, regions) without M4 touching it.
2. Gene numbers for calibration genomes: 1,000 (fast) and 4,500 (full, E. coli-like)? **Proposed
   default:** yes; per-PR suites must fit 10 minutes on 4 CPUs.

## Lead decisions after code review (PR #62, changes requested)
1. Circular replicons get the wrap pair (last gene, first gene) in `operon_pairs`, as M2-12 does (layout and
   M0 adapter); first/last genes then have `uniop_p_prev/next` like production.
2. A layout option `origin_crossing: bool` (default False, so existing outputs are unchanged) rotates `chr1`
   so that one gene, one TU and one block region straddle the origin (`crosses_origin` true); M4-11, M4-15
   and M4-17 calibration configs turn it on.
3. Tests kill the surviving mutants: adapter pair p 1 within / 0 across TUs; genes without a truth TU become
   singleton TUs (never merged); region numbering by contig then start (a config with two regions of the same
   type); `_MIN_AA = 20` (no sORF-sized CDS); `region_ids` sorted and complete; a linear contig is never
   `kind = chromosome`.
4. `amr` is written (= membership in an AMR region); `ribosomal_protein` stays unwritten (no M4 consumer).
   The dead `"p"` fallback in `_uniop_lookup` is removed.
5. Fix the M0 simulator crash with `real_rrna_trna=True, n_trna=0` (`blocks.trna_picks` on an empty list) and
   test it.
6. Deferred to M4-17 (noted in its spec): overlapping regions and per-contig region placement (`RegionSpec.contig`,
   e.g. a 60-gene plasmid with its own AMR set).
