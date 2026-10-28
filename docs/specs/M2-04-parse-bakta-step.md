# M2-04: `parse_bakta` step: gene and feature tables, 16S and anti-SD, tRNA counts, xrefs, TE6

Milestone M2 · depends on: M2-02, M2-03 · estimate ≈ 300 code lines (tests not counted) · capabilities: none

## Spec sections implemented

- §6.2 `parse_bakta` | py | "Bakta JSON → `genes.parquet`, `features.parquet` (rRNA, tRNA, CRISPR,
  oriC, gaps), `rrna16s.fa`, `trna_counts.tsv`, `xrefs.parquet`".
- §7.1.2: "**Ribosomal-protein genes** (the CAI reference) are those whose product matches
  `^(30S|50S) ribosomal protein` or whose gene name matches `^(rps|rpl|rpm)[A-Z]`." · "**Anti-Shine–
  Dalgarno sequence.** Take the last 9 nt of each 16S rRNA (the most common 9-mer across copies). If
  no 16S exists, or the 9-mer lacks `CCTCC`, use the E. coli default `ACCTCCTTA` and flag
  `anti_sd_default`." · "**tRNA gene copy numbers** are counted per anticodon from Bakta tRNA features, for tAI."
- §5.2 `annotation.ostir.anti_sd: auto  # auto = 3′ end of the strain's 16S rRNA`.
- §9.3 `genes.tsv` `annotation_source` enum{bakta_full, bakta_light, user}; §9.4 run flags
  `anti_sd_default`, `bakta_db_light`; gene flag `cai_reference`.
- §11.6 TE6: "on K-12, Bakta's oriC call sits next to mnmG (gidA): mnmG or mioC is the nearest gene,
  and mnmG is within 2 kb".

## Public interface

`src/gippa/annotate/bakta_layers.py`
```python
RIBOSOMAL_PRODUCT: re.Pattern[str]   # ^(30S|50S) ribosomal protein
RIBOSOMAL_GENE: re.Pattern[str]      # ^(rps|rpl|rpm)[A-Z]
DEFAULT_ANTI_SD = "ACCTCCTTA"

def is_ribosomal_protein(g: Gene) -> bool
def rrna16s(features: Sequence[NcFeature]) -> list[NcFeature]       # rRNA with product containing "16S" or gene "rrs"
def anti_sd(rrna16s_nts: Sequence[str], configured: str = "auto") -> tuple[str, tuple[str, ...]]  # (9-mer, run flags)
def trna_counts(features: Sequence[NcFeature]) -> dict[str, int]   # anticodon → copies (sorted)
def split_xrefs(genes: Sequence[Gene]) -> list[tuple[str, str, str]]  # (locus_tag, db, id)
def annotation_source(db_type: str) -> Literal["bakta_full", "bakta_light", "user"]
def nearest_genes(feature: NcFeature, genes: Sequence[Gene], replicon: Replicon, n: int = 2) -> list[tuple[str, int]]
```
`src/gippa/steps/annotation.py`: `@step("parse_bakta", inputs=["bakta_json", "genome"],
outputs=["genes", "features", "replicons", "rrna16s", "trna_counts", "xrefs", "annotation_meta"],
precomputed=None)`.

Output files (work dir, internal; 0-based half-open in Parquet):
- `genes.parquet`: `locus_tag, gene, product, feature_type, contig, start, end, strand, partial,
  crosses_origin, hypothetical, aa, nt, start_codon, ribosomal_protein, amr`.
- `features.parquet`: `feature_id, type, contig, start, end, strand, gene, product, amino_acid, anticodon`.
- `replicons.parquet`: `contig, length, topology, kind` (FASTA order).
- `rrna16s.fa`: one record per 16S copy, header `<feature_id>`, 5'→3'.
- `trna_counts.tsv`: `anticodon`, `amino_acid`, `n_copies` (sorted by anticodon), via `write_tsv`.
- `xrefs.parquet`: `locus_tag, db, id` (db = text before the first `:`; for `GO` the id keeps the
  `GO:` prefix, e.g. `GO:0005524`; others drop it, e.g. `COG0593`, `K02313`, `3.6.4.12`).
- `annotation_meta.json`: `bakta_version`, `db_version`, `db_type`, `annotation_source`, `anti_sd`,
  `flags` (sorted), `n_genes`, `n_features`.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | seven 16S copies all ending `…GATCACCTCCTTA` | `anti_sd` → `("ACCTCCTTA", ())` |
| 2 | no 16S rRNA | `("ACCTCCTTA", ("anti_sd_default",))` |
| 3 | one 16S ending `…GATCATTTAAG` (9-mer `TCATTTAAG` lacks `CCTCC`) | `("ACCTCCTTA", ("anti_sd_default",))` |
| 4 | copies ending `…ACCTCCTTA` ×2 and `…ACCTCCTTT` ×1 | `("ACCTCCTTA", ())` |
| 5 | `configured="TCCTCCACT"` (user override) | `("TCCTCCACT", ())`, 16S ignored |
| 6 | product `50S ribosomal protein L2`, gene `rplB` / product `DNA-directed RNA polymerase subunit beta`, gene `rpoB` / gene `rpsA`, product `hypothetical protein` | `True` / `False` / `True` |
| 7 | tRNAs with anticodons `CAG, CAG, GCC` | `{"CAG": 2, "GCC": 1}` |
| 8 | db_xrefs `["COG:COG0593","GO:0005524","EC:3.6.4.12","KEGG:K02313","UniRef:UniRef90_P03004"]` on `b3702` | rows `(b3702,COG,COG0593)`, `(b3702,GO,GO:0005524)`, `(b3702,EC,3.6.4.12)`, `(b3702,KEGG,K02313)`, `(b3702,UniRef,UniRef90_P03004)` |
| 9 | `db_type="light"` | `annotation_source="bakta_light"`, flags contain `bakta_db_light` |
| 10 | TE6: recorded K-12 full JSON, oriC feature | `nearest_genes(oriC)[0]` locus has gene `mnmG` or `mioC`; `mnmG` within 2,000 bp |

## Edge cases

- Ties in the most common 9-mer: the 9-mer of the first copy in genome order wins.
- A 16S shorter than 9 nt (truncated at a contig end) is skipped.
- Distances in `nearest_genes` are gap lengths in bp (0 when overlapping), origin-aware on circular
  replicons; linear contigs never wrap.
- A gene can be both `ribosomal_protein` and `pseudogene`; the CAI reference (M2-17) excludes pseudogenes.
- Two writes of the same input give identical bytes (Parquet written without timestamps:
  `pyarrow` with `write_statistics=False` and fixed metadata).

## Acceptance criteria

- [ ] `tests/unit/test_bakta_layers.py`: examples 1–9.
- [ ] `tests/unit/test_parse_bakta_step.py`: run the step on the simulator output; every listed file
  exists; `genes.parquet` row count = simulator gene count; byte-identical on rerun (TE7).
- [ ] `tests/contract/test_te6_oric.py` (example 10) on the recorded K-12 full-DB JSON;
  `pending_sample` until M2-10.

## Questions for the lead

1. Is a user override of the anti-SD allowed? §5.2 has `anti_sd: auto`. **Proposed default:** any
   other value is a 9-nt DNA string used as is (validated `^[ACGT]{9}$`).
2. `annotation_source = user` is undefined. **Proposed default:** `bakta_full`/`bakta_light` from the
   JSON DB type whether Bakta ran or the JSON was supplied; `user` only when the DB type is unknown
   (e.g. the simulator JSON or a non-Bakta source).

## Lead decisions (test-author questions)
1. Q1 accepted as proposed: the anti-SD override must match `^[ACGT]{9}$`, else `ValueError`.
2. Unknown Bakta feature types are counted in `annotation_meta.json` key `skipped_types` (type → count,
   `{}` when none) and each is named in a WARNING. ncRNA `class` is not carried (M2-03 decision 7).

## Lead decisions after code review (PR #36, approved)
1. Tests pin: genes/features ordered by FASTA replicon order, then start; "Putative 30S ribosomal
   protein" is not ribosomal (product match anchored at the start); a tRNA without anticodon writes no
   anticodon row; `nearest_genes` ties keep input order; parquet columns carry no statistics; the 16S
   tail is upper-cased.
2. The deterministic parquet writer moves to `gippa.io` (public, reusable by later steps); steps import
   no private helpers (`_slice` becomes public or moves with it).
