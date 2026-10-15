# M3-23: `giap makedb`, optional strain search database

Milestone M3 · depends on: M3-01, M3-02, M3-14 (`classify_locus`), M1-05 (GFF3 reader), M1-03 ·
creates `src/giap/makedb/{__init__,build}.py`, wires the `makedb` CLI stub · pure Python ·
estimate ≈ 330 code lines

## Spec sections implemented

- §7.9.2 "`giap makedb` (convenience; no FDR is built on it)": "**Writes `DB/strain_db.faa`:**
  Bakta CDS and sORF proteins (`>{locus_tag} {product}`); with `--patch annotation_patch.gff3`, the
  patched proteins; with `--six-frame-orfs`, class-tagged six-frame ORFs ≥ `--min-orf-aa` (20) that
  do not overlap annotated CDS in frame: `>GIAPORF_{class}_{contig}_{start1}_{end1}_{p|m}`, classes
  `intergenic`, `cds_outframe`, `antisense`, `nterm_extension`; contaminants as `>CONT_{id}`
  (unless `--no-contaminants`). **Also writes** `DB/id_map.tsv` and `DB/makedb_manifest.json`.
  **No decoys.** Search engines generate their own; with DIA-NN, use `--report-decoys` so decoy
  rows reach GIAP. If no annotation is given, Bakta runs first."
- §5.1 "optional strain search database for the user's own search (no FDR is built on it)".
  Improvement 3. §15.2 step 2 (pilot A2 database). §6.6 determinism (no timestamps in the FASTA;
  the manifest may carry a date).

## Public interface

```python
@dataclass(frozen=True)
class MakedbOptions:
    annotation: Path; genome: Path | None; outdir: Path
    patch: Path | None = None; six_frame_orfs: bool = False; min_orf_aa: int = 20
    contaminants: Path | None = None        # library FASTA from `giap setup`; None with --no-contaminants
def six_frame_orfs(index: SixFrameIndex, gm: GeneModel, genome: Mapping[str, str],
                   min_aa: int) -> list[OrfRecord]: ...
def build_strain_db(opts: MakedbOptions) -> MakedbResult: ...
```

CLI: `giap makedb --annotation JSON [--genome FASTA] --outdir DB [--patch GFF3] [--six-frame-orfs]
[--min-orf-aa 20] [--no-contaminants]`. The genome comes from the Bakta JSON `sequences` when
`--genome` is absent. Without `--annotation`: exit 2 with "Bakta wrapper required (M2)" until M2's
wrapper exists, then Bakta runs first.

- **ORFs** are maximal stop-to-stop segments of each six-frame translation (M3-02, no I → L in the
  written sequence: translate from the genome directly), ≥ `min_orf_aa`, without `X`. An ORF that
  overlaps an annotated CDS or sORF **in frame** is dropped, except for its upstream part, which
  becomes the `nterm_extension` entry: from the codon after the previous in-frame stop to the
  annotated start, followed by the first 30 aa of the annotated protein (so junction peptides
  exist), kept when the upstream part is ≥ 7 aa (Q1). Other ORFs take the class of
  `classify_locus` (`intergenic`, `cds_outframe`, `antisense`); ORFs whose class would be
  `rna_overlap` or `pseudogene` are omitted.
- **Order:** CDS and sORF proteins in gene-model order; patched proteins; ORFs sorted by contig,
  start, end, strand; contaminants in library order. Wrapped at 60 residues, no `*`.
- **`id_map.tsv`:** `db_id`, `kind` enum{cds, sorf, patch, orf, contaminant}, `locus_tag`,
  `class`, `contig`, `start`, `end`, `strand` (1-based closed; NA for contaminants).
- **`makedb_manifest.json`:** GIAP version, inputs with SHA-256, options, counts per kind, date.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | sim 50-kb JSON, defaults | `strain_db.faa` has one record per CDS/sORF, header `>GSIM_00005 hypothetical protein`, sequence = `aa` |
| 2 | `--no-contaminants` | no `CONT_` record |
| 3 | contaminant library with `P00761` | record `>CONT_P00761` |
| 4 | `--six-frame-orfs` on a genome with an intergenic ORF `[1200, 1290) +` (30 aa) | `>GIAPORF_intergenic_chr1_1201_1290_p` |
| 5 | an antisense ORF on the minus strand `[300, 390)` | `>GIAPORF_antisense_chr1_301_390_m` |
| 6 | an ORF of 19 aa (default `min_orf_aa=20`) | absent |
| 7 | CDS with an in-frame upstream extension of 12 aa | `>GIAPORF_nterm_extension_…` whose sequence is the 12 aa + the first 30 aa of the protein |
| 8 | the same inputs twice | byte-identical `strain_db.faa` and `id_map.tsv` |
| 9 | `--patch` GFF3 with one CDS `ID=NOVEL_chr1_1201_p` | record `>NOVEL_chr1_1201_p` translated from the genome |

## Edge cases

- Origin-crossing ORFs on circular contigs: coordinates `start1 > end1` are not used; the header
  carries `end1 > length` (unwrapped, as in TSV). Linear contigs: ORFs truncated at the contig end
  are kept only if ≥ `min_orf_aa` (they have no stop; flag in `id_map.tsv` `class` unchanged).
- Bakta proteins with a terminal `*` are written without it; pseudogenes are not written.
- Duplicate sequences under different IDs are all written (engines handle them; `id_map.tsv`
  keeps provenance).

## Acceptance criteria

- [ ] `tests/unit/test_makedb.py`: examples 1–9.
- [ ] CLI test: `giap makedb --annotation … --outdir …` exit 0; missing annotation exit 2 with the
  message above.
- [ ] Property: every in-silico tryptic peptide of every `strain_db.faa` record (except
  contaminants and the junction part of `nterm_extension`) maps exactly with M3-04 to the locus
  given in `id_map.tsv`.
- [ ] Benchmark: `--six-frame-orfs` on a 5-Mb genome ≤ 30 s.

## Questions for the lead

1. Shape of the `nterm_extension` entries (not fixed by the spec). **Proposed default:** upstream
   part + first 30 aa of the annotated protein, upstream part ≥ 7 aa.
2. Should A2 (pilot) use `--no-contaminants` because the DIA-NN template already passes
   `--fasta contaminants.fasta`? **Proposed default:** yes, to avoid duplicate entries (M3-24).

## Lead decisions (after the test-author's pass)
1. The test-author's resolutions are accepted (missing product → "hypothetical protein"; patched proteins
   start with M; six-frame ORFs translated codon by codon; db_ids unique, id_map in FASTA order with NA;
   patched records after CDS/sORF and before ORFs; contaminants from `--db-dir`, exit 2 when missing
   unless `--no-contaminants`).
2. nterm_extension entries contain the upstream part plus the first 30 codons of the annotated protein (so
   peptides spanning the junction are searchable); the header and id_map say so.
3. An ORF containing X is split at each X; pieces shorter than `min_orf_aa` are dropped.

## Lead decisions after code review (PR #79, changes requested)
4. nterm_extension is undefined (no entry) when the upstream walk runs off a linear contig or meets an
   X; the oracle and property test skip those tags.
5. The junction part (first 30 codons of the annotated protein) is translated codon by codon from the
   genome (GTG → V, TTG → L inside the entry), not copied from Bakta's `aa` (which starts with M).
6. The borrowed private helpers become public APIs (M3-08 precedent): M3-14 `LocusClassifier(gm,
   genome).classify(...)` (`classify_locus` delegates to it) and `in_frame`; M3-01
   `translate11_text(nt, offset)` and the ambiguity sentinel as a public constant; M3-23 uses only those.
7. Rejected inputs (exit 1, ValueError): a patch CDS whose translation has an internal `*`; an annotation
   JSON with no features (`{}`) — never an empty database with exit 0.
8. The makedb id map is written as `db_id_map.tsv` (avoids the clash with C14 `annotation/id_map.tsv`).
   Resolving `GIAPORF_`/`NOVEL_`/`*_patch` protein IDs through it (instead of `cross_strain_protein`) is
   backlog for M3-17/M4.
9. Tests kill the surviving mutants: nterm interval pinned (upstream + 30 codons only); X in the upstream
   part; `*` stripping for contaminants and the fallback translation; `--genome` length mismatch; duplicate
   ORFs at the origin; in-frame drop across the origin.
