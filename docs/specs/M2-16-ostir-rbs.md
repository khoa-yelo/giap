# M2-16: OSTIR input builder, parser and RBS features

Milestone M2 · depends on: M2-04 (genes, anti-SD) · estimate ≈ 180 code lines (tests not counted) ·
capabilities: none (recording needs `needs_conda`)

## Spec sections implemented

- §7.1.1 OSTIR: "OSTIR 1.1.3 (GPL) | own env with ViennaRNA + `pip install ostir==1.1.3` | `ostir -i
  ostir_input.csv -o ostir.csv -j N` | `expression` (TIR), `dG_total`, `dG_rRNA:mRNA`, `dG_mRNA`,
  `dG_spacing`, `dG_standby`, `dG_start_codon`, `RBS_distance_bp`". Gotcha: "`start_position` is
  1-based in the input sequence. Restrict to the annotated start with CSV `start`/`end` columns."
- §7.1.6: "**OSTIR input.** One row per gene: `id` = locus tag; `seq` = −50 to +50 nt around the
  start codon, strand-aware (clipped at contig ends); `start = end = 51`, or the clipped equivalent;
  `anti-Shine-Dalgarno` from §7.1.2." Per-gene features: "`tir_ostir` = `expression`; `sd_dg` =
  `dG_rRNA:mRNA`; `mrna_dg` = `dG_mRNA`; `spacing_dg`, `standby_dg`; `spacer_bp` = `RBS_distance_bp`; `start_codon`."
- §5.2 `ostir: {upstream_nt: 50, downstream_nt: 50, anti_sd: auto}`. §9.3 `genes.tsv` translation
  initiation: `tir_ostir`, `sd_dg`, `mrna_dg`, `spacer_bp`, `start_codon`.

## Public interface

`src/gippa/annotate/ostir.py`, `src/gippa/annotate/ostir_cmd.py`; steps `ostir_input`, `rbs_layer`.
```python
@dataclass(frozen=True, slots=True)
class OstirRow: id: str; seq: str; start: int; end: int; anti_sd: str

def ostir_row(gene: Gene, replicon: Replicon, genome_seq: str, anti_sd: str, up=50, down=50) -> OstirRow
def write_ostir_input(path: Path, rows: Sequence[OstirRow]) -> None   # header: id,seq,start,end,anti-Shine-Dalgarno (Q1)
def ostir_command(inp: Path, out: Path, threads: int) -> list[str]     # ["ostir","-i",inp,"-o",out,"-j",N]

@dataclass(frozen=True, slots=True)
class RbsFeatures:
    tir_ostir: float; sd_dg: float; mrna_dg: float; spacing_dg: float; standby_dg: float
    spacer_bp: int; start_codon: str
def parse_ostir(path: Path, requested: Mapping[str, int]) -> dict[str, RbsFeatures]
    # keeps rows whose start_position equals the requested start; one per id; missing → absent
```
Sequence geometry (S = first start-codon base, E = interval end, 0-based):
- `+`: `seq = genome[S − 50, S + 50)`; `start = end = 51`.
- `−`: `seq = revcomp(genome[E − 50, E + 50))`; `start = end = 51`.
- Circular: wrap. Linear: clip; `start = end = 1 + (bases kept upstream)`.
`start_codon` in `RbsFeatures` is the gene's start codon (from `Gene`), not OSTIR's guess.

## Examples (chr1 circular, 100,000 bp; anti-SD `ACCTCCTTA`)

| # | Input | Expected |
|---|---|---|
| 1 | gene `(1000, 1300, +)` | `OstirRow(id, genome[950:1050], 51, 51, "ACCTCCTTA")`; `seq[50:53]` = start codon |
| 2 | gene `(2000, 2900, −)` | `seq = revcomp(genome[2850:2950])`; `seq[50:53]` = start codon (revcomp of `genome[2897:2900]`) |
| 3 | gene `(20, 300, +)` on a linear contig | `seq = genome[0:70]`, `start = end = 21` |
| 4 | same gene on circular chr1 | `seq = genome[99970:100000] + genome[0:70]`, `start = end = 51` |
| 5 | − gene `(99960, 99990, −)` on linear 100,000-bp contig | `seq = revcomp(genome[99940:100000])` (60 nt: 10 upstream kept), `start = end = 11` |
| 6 | `ostir.csv` row `name=GSIM_00005,start_codon=ATG,start_position=51,expression=1523.4,RBS_distance_bp=7,dG_total=-3.1,dG_rRNA:mRNA=-9.8,dG_mRNA=-4.2,dG_spacing=0.3,dG_standby=0.0,dG_start_codon=-1.19` | `RbsFeatures(1523.4, -9.8, -4.2, 0.3, 0.0, 7, <gene start codon>)` |
| 7 | rows for the same id at positions 51 and 63 | only 51 kept |
| 8 | no row for a gene (e.g. non-canonical start OSTIR skips) | gene absent → `NA` features downstream |

## Edge cases

- Genes shorter than 50 nt downstream (sORFs) still get +50 nt (reads into the next gene); fine.
- Pseudogenes and partial genes get rows too; partial genes on linear contig ends use the clip rule.
- `anti_sd` comes from `annotation_meta.json` (M2-04); with `anti_sd_default` it is `ACCTCCTTA`.
- Rows are written in §6.5 gene order; the file is byte-identical on rerun.

## Acceptance criteria

- [ ] `tests/unit/test_ostir.py`: examples 1–8 (hand-written CSV for 6–8).
- [ ] `tests/property/test_ostir_rows.py`: for random genes on circular replicons `len(seq) == 100`
  and `seq[start-1:start+2] == gene.start_codon`; on linear contigs the same identity holds after clipping.
- [ ] `tests/contract/test_ostir_contract.py` on the recorded sim100k `ostir.csv`; `pending_sample` until recorded.

## Questions for the lead

1. The spec names the id column `id`; OSTIR's documented CSV header may be `name`. **Proposed
   default:** write whatever header OSTIR 1.1.3 accepts (checked once with the env), map it back to
   `id` in GIPPA; the recorded fixture fixes it.

## Lead decisions (after the test-author's pass)
1. The test-author's resolutions are accepted, incl. the M2-index change: `ostir_input` also takes
   `genome` and `replicons` (annotate.nf `stepIn('ostir_input', [... reps, kv('genome', genome)])`);
   registry = the regenerated workflow/tests/data/annotate/step_registry.json; `rbs.parquet` = the
   precomputed `rbs` columns, rows in input order, genes without an OSTIR row left out (yecM in K-12 and
   Nissle is the real example); start codon from the gene's input row.
