# M2-12: UniOP adapter, distance-model fallback, operon pairs, transcription units, TU upstream regions

Milestone M2 · depends on: M2-04 · estimate ≈ 380 code lines (tests not counted) · capabilities: none
(recording the UniOP fixture via the M2-10 harness needs `needs_conda`)

## Spec sections implemented

- §7.1.1 UniOP: "`python <UniOP>/src/UniOP -f cds.gff --gff_source prodigal_gff -t outdir` (confirm
  the entry-point name on the pinned commit)"; keeps "probability for each adjacent same-strand pair (`uniop.pred`)".
- §7.1.3: "**UniOP adapter.** Write the Bakta CDS and sORFs as a Prodigal-style GFF in genome order,
  keep a map back to locus tags, run UniOP, and parse `uniop.pred`." "**Fallback.** When UniOP fails,
  the distance model gives p(operon | d) = 1 / (1 + exp((d − d₀)/s)), where d is the intergenic
  distance in bp (negative for overlaps). Defaults are d₀ = 50, s = 20; … stored in
  `src/giap/data/calibration.yaml`. The manifest records `operon_source = distance_model`." "**TUs**
  are maximal runs of same-strand adjacent genes linked by p ≥ `uniop.pair_threshold`." "**TU upstream
  region.** From −300 to +50 bp around the first start codon, strand-aware, clipped at the end of the
  nearest upstream gene on either strand."
- §5.2 `uniop: {pair_threshold: 0.5, fallback: distance}`, `tu_upstream_bp: [-300, 50]`.
- §9.3 `annotation/operon_pairs.tsv`: "`gene_a`, `gene_b`, `uniop_p`, `source` enum{uniop,
  distance_model}"; `annotation/tus.tsv`: "`tu_id`, `contig`, `start`, `end`, `strand`, `genes`
  (list), `upstream_start`, `upstream_end`". §6.5 "TUs: `TU_<n>`". §9.4 run flag `operon_source_distance_model`.
- Module 1 tests: "an adapter round-trip for UniOP (Bakta → Prodigal-style GFF → back to locus tags)";
  "the distance-model fallback".

## Public interface

`src/giap/annotate/uniop_adapter.py`, `src/giap/annotate/uniop_cmd.py`, `src/giap/annotate/operons.py`,
`src/giap/data/calibration.yaml` (new; `operon_distance: {d0: 50, s: 20, provenance: "GIAP_SPEC §7.1.3 default", calibrated: false}`,
`terminatornet_yes: {value: 0.5, calibrated: false}`), steps `uniop_input`, `operon_pairs`.
```python
def write_prodigal_gff(path: Path, genes: Sequence[Gene], replicons: Sequence[Replicon]) -> list[tuple[str, str]]
    # returns [(prodigal_id, locus_tag)]; IDs "<contig>_<k>" with k = 1-based index per contig in genome order
def parse_uniop_pred(path: Path, id_map: Mapping[str, str]) -> list[PairProb]
def uniop_command(gff: Path, outdir: Path, uniop_home: Path) -> list[str]

@dataclass(frozen=True, slots=True)
class PairProb: gene_a: str; gene_b: str; uniop_p: float; source: Literal["uniop", "distance_model"]

def adjacent_pairs(genes: Sequence[Gene], replicons: Mapping[str, Replicon]) -> list[tuple[Gene, Gene]]
    # consecutive in genome order per contig, same strand; wrap pair (last, first) on circular replicons
def intergenic_distance(a: Gene, b: Gene, replicon: Replicon) -> int   # b.start − a.end (genome order), origin-aware; < 0 = overlap
def distance_model_p(d: int, d0: float, s: float) -> float
def build_pairs(genes, replicons, uniop: Sequence[PairProb] | None, d0: float, s: float) -> tuple[list[PairProb], tuple[str, ...]]

@dataclass(frozen=True, slots=True)
class Tu:
    tu_id: str; contig: str; interval: Interval; strand: str
    genes: tuple[str, ...]          # transcription order
    upstream: Interval

def build_tus(genes, pairs, threshold: float, replicons, upstream_bp=(-300, 50)) -> list[Tu]
def tu_upstream(first: Gene, others: Sequence[Gene], replicon: Replicon, bp=(-300, 50)) -> Interval
def write_operon_pairs_tsv(path, pairs) / write_tus_tsv(path, tus)
def write_tu_regions(outdir: Path, tus: Sequence[Tu], genome: Mapping[str, str]) -> list[Path]  # one FASTA per TU (used by M2-14)
```
Pair orientation: `gene_a` is the earlier gene in genome order (lower start), `gene_b` the next; on
the wrap pair `gene_a` is the last gene of the contig. Genes = CDS, sORF and pseudogene (Q1).
UniOP is used for every pair it reports; pairs it does not report (incl. the wrap pair) get the
distance model with `source = distance_model`. UniOP failure (non-zero exit, missing or empty
`uniop.pred`) → all pairs from the distance model and run flag `operon_source_distance_model`.

## Examples (d₀ = 50, s = 20; chr1 circular, 100,000 bp)

| # | Input | Expected |
|---|---|---|
| 1 | A `(100,400,+)`, B `(420,700,+)` | d = 20; `distance_model_p` = 1/(1+e^−1.5) = **0.817574** |
| 2 | A `(100,400,+)`, B `(396,700,+)` (4-bp overlap) | d = −4; p = 1/(1+e^−2.7) = **0.937027** |
| 3 | d = 50 | p = 0.5; linked at `pair_threshold` 0.5 (≥) |
| 4 | `uniop.pred` row for `chr1_1`/`chr1_2` = 0.91, map `{chr1_1: GSIM_00005, chr1_2: GSIM_00010}` | `PairProb("GSIM_00005","GSIM_00010",0.91,"uniop")` |
| 5 | genes A+ B+ C− D− with p(A,B)=0.9, p(C,D)=0.3 | `TU_1 = (A,B)`, `TU_2 = (C,)`, `TU_3 = (D,)`; C/D TUs list genes in transcription order |
| 6 | TU_1 first gene + start S = 1000 (0-based); nearest upstream gene (either strand) ends at 850 | upstream `(850, 1050, +)`; TSV `upstream_start=851`, `upstream_end=1050` |
| 7 | − TU whose first gene is `(5000,5600,−)`; nearest upstream gene (higher coordinates) starts at 5800 | upstream `(5550, 5800, −)`; TSV `5551`, `5800` |
| 8 | last gene `(99800,99950,+)`, first gene `(30,300,+)` on circular chr1 | wrap pair, d = 80; on a linear contig no pair |
| 9 | UniOP exit 1 | every pair `source=distance_model`; flags `("operon_source_distance_model",)` |
| 10 | round trip: simulator genes → `write_prodigal_gff` → fake `uniop.pred` listing all adjacent same-strand prodigal IDs → `parse_uniop_pred` | the locus-tag pairs equal `adjacent_pairs` minus the wrap pair |

## Edge cases

- A TU spanning the origin on a circular contig has `end > length` (written unwrapped in TSV, §6.5).
- Upstream clip when the upstream gene overlaps the first gene: the region never starts after the
  start codon (minimum region = `[S, S+50)` on +), see Q2.
- Upstream region on a linear contig is clipped at the contig ends; on circular it wraps.
- A contig with one gene → one TU, no pairs. No genes → empty tables with headers.
- TU numbering: `TU_<n>` by contig (FASTA order), then TU start, then strand.
- UniOP IDs must survive UniOP's own parsing: prodigal-style `ID=<contig>_<k>` attributes; contig
  names are already safe IDs.

## Acceptance criteria

- [ ] `tests/unit/test_operons.py`: examples 1–3, 5–9; `tests/property/test_tus.py`: TUs partition
  the genes; every TU is single-strand; consecutive genes in a TU have p ≥ threshold.
- [ ] `tests/unit/test_uniop_adapter.py`: examples 4, 10.
- [ ] `tests/contract/test_uniop_contract.py` on `tests/data/recorded/uniop/<git12>/{sim100k,k12}/`
  (recorded with the M2-10 harness on the Bakta full-DB genes); `pending_sample` until recorded.
- [ ] `calibration.yaml` created; values read through one loader (`giap.data.calibration.load()`).

## Questions for the lead

1. Include pseudogenes in UniOP input and TUs? **Proposed default:** yes (they are rows of
   `genes.tsv` and break or join TUs physically); non-coding RNAs between genes are ignored here
   (M4 flags `rna_between`).
2. "Clipped at the end of the nearest upstream gene" when that gene overlaps the first start codon.
   **Proposed default:** clip to the start codon, so the region is at least `[S, S+50)`.
3. `uniop.pred` column names are not documented in the spec. **Proposed default:** the parser keys
   on the recorded header; until recorded, a doc-derived 3-column file (`gene_a gene_b p`) is used in
   unit tests only.

## Lead decisions (test-author questions)
1. List-valued columns in Parquet work tables are native lists (same rule as M4-06 decision 1): in
   `tus.parquet`, `genes` is a list of locus tags in TU order; `;`-joining (ADR-0006) only in TSV outputs.
2. Accepted conventions: replicons as a mapping in FASTA order; a circular contig with two genes has both
   the forward and the wrap pair; an empty `uniop.pred` is a UniOP failure; nearest-upstream clip =
   `min(gene end, S)` (Q2); unit-test `uniop.pred` uses the 3-column `gene_a gene_b p` format (Q3);
   `uniop_map.tsv` header `prodigal_id locus_tag`; 0-based coordinates in Parquet.
3. With `annotation.precomputed.operon_pairs`, the precomputed layer is read by `operon_pairs` directly
   (its `source` column kept); exact recovery of simulator `truth/tus.tsv` from 0/1 oracle probabilities
   is the accuracy check (the spec sets no calibration threshold).
4. The adapter module is `giap.annotate.uniop_adapter` (not `uniop.py`): the vendoring guard
   (`test_no_tool_source_vendored`) forbids files named after unlicensed tools under src/ (ADR-0008).
5. Accepted implementation notes: a circular TU is capped at the ring length; the distance-model fallback
   is logged as a WARNING (and recorded as `operon_source = distance_model` in the manifest).

## Lead decisions after code review (PR #49, changes requested)
1. Parser follows the REAL `uniop.pred` format (pin 3ed14dfc): a title line `operon prediction scores`, a
   header `Gene A\tGene B\tPrediction`, rows of 1-based gene ROW NUMBERS (not IDs); opposite-strand and
   cross-contig rows have an empty probability and are skipped.
2. Row numbers index UniOP's own gene order: all genes stable-sorted by start across contigs, as UniOP reads
   them from the input GFF/Prodigal file. Row k = sorted gene k; gene_b = the next gene ON THE SAME CONTIG.
   The adapter rebuilds that order from the exact input the builder wrote, requires both genes on the same
   contig and strand, and maps to locus tags. Unit test on a multi-contig genome.
3. Record a real UniOP run on SIM100K as `tests/data/recorded/uniop/3ed14dfc/sim100k/` (built `giap-uniop`
   env; run as a process; PROVENANCE.md with command, env lock and input sha256); the contract test is
   rewritten to the real format and must pass on it (ID-based assumptions removed).
4. `workflow/modules/local/uniop/main.nf` creates UniOP's `-t` output directory (`mkdir -p out`).
5. The distance-model fallback is recorded, not only logged: the operon step declares an output
   `operon_meta.json` (`operon_source`: `uniop` | `distance_model` | `precomputed`, plus the reason), from
   which the manifest takes `operon_source` and the run flag `operon_source_distance_model` (§7.1.3, §9.4).
6. Tests kill the surviving mutants: first gene of a minus-strand multi-gene TU (`run[0]` vs last),
   WARNING level of the fallback, precomputed `source` kept, d0/s read via `calibration.load()`,
   `tu_upstream_bp` honoured, the ring-length cap and the ring self-bound.
Rule from here on (lead): every tool-wrapping issue records a real-tool fixture before its PR is reviewed.

## Test-author pins for the code-review decisions (items 1, 2, 4, 5, 6)
- Items 1-2: `parse_uniop_pred(path: Path, gff: Path, id_map: Mapping[str, str]) -> list[PairProb]`,
  where `gff` is the exact file `write_prodigal_gff` wrote. Rows in file order; empty probability cells
  skipped; the printed `Gene B` number is not used for pairing. `ValueError` on: missing title/header
  (incl. the invented `gene_a gene_b p` file), a row number outside the gene list, a row with a
  probability whose partner is absent or on the other strand, a GFF ID missing from `id_map`. Empty
  file or title+header only → `[]` (a UniOP failure for `build_pairs`). Fixture: real excerpt in
  `tests/data/uniop/` (plasmid gene sorts first, so row 1 is the plasmid pair); ties on start keep
  GFF file order. `operon_pairs` reads `uniop_cds` (required or optional input).
- Item 4: the UNIOP `script:` runs `mkdir -p <the -t dir>` before `src/UniOP`; the `stub:` writes the
  real title and header.
- Item 5: `operon_pairs` outputs `("pairs_raw", "tus", "tu_upstream", "operon_meta")`;
  `operon_meta.json` is a JSON object with `operon_source` ∈ {`uniop`, `distance_model`,
  `precomputed`} and a non-empty `reason` string.
- Item 6: the fallback logs a WARNING whose message contains "distance" (no such record on a UniOP
  run); d₀/s are read by calling `giap.data.calibration.load()` at step run time (monkeypatched in tests).
- Contract: `tests/data/recorded/uniop/3ed14dfc/{sim100k,k12}/` holds one `*.gff`, `uniop.pred`,
  `uniop_map.tsv` and `PROVENANCE.md` containing `3ed14dfc` and the GFF's sha256.
