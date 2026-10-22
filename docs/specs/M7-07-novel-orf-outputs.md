# M7-07: `novel_orfs.tsv`, `novel_orfs.gff3`, `novel_orfs.faa` and `annotation_patch.gff3`

Milestone M7 · track N (novel ORFs) · **pass 1** · depends on: M7-06; merged M1-05 GFF3 writer,
M1-03 protein FASTA writer, M3-23 `giap makedb --patch` · estimate ≈ 300 code lines (+ ≈ 400 test
lines) · not statistical

## Spec sections implemented

- §7.6.6 outputs: "`tables/novel_orfs.tsv`; `novel_orfs.gff3`; `novel_orfs.faa` (both start
  variants); `annotation_patch.gff3`: new ORFs as CDS (IDs `NOVEL_*`), and N-terminal extensions as
  the existing locus tag with the new start and `Note=nterm_extension`. `giap makedb --patch` uses
  it to rebuild the search database for a second-pass search." Tests: "GFF3 validity; patch →
  `makedb` round trip".
- ADR-0013 decision 3: "`annotation_patch.gff3` writes it like an N-terminal extension (existing
  locus tag, new end, `Note=cterm_extension`)".
- ADR-0012 decision 3: the four columns after `promoter_score` (`promoter_supported`,
  `terminator_downstream`, `terminator_conf`, `terminator_distance_nt`), and "The novel-ORF GFF3
  carries the same as attributes".
- §9.3 `tables/novel_orfs.tsv` columns; §9.1 (all four files under `tables/`); §6.5 (GFF3 1-based
  closed, origin-crossing `end > length`, `Is_circular=true`); §9.5 (GFF3 round trip, Frictionless).

## Public interface

`src/giap/novel/outputs.py` (pure: rows, features, records)
```python
NOVEL_ORFS_COLUMNS: tuple[str, ...]          # §9.3 order + the ADR-0012 columns after promoter_score
def orf_table(orfs: Sequence[Orf], scored: pd.DataFrame | None) -> pd.DataFrame: ...
def gff_features(table: pd.DataFrame) -> list[GffFeature]: ...
def patch_features(table: pd.DataFrame, gm: GeneModel) -> list[GffFeature]: ...
def protein_records(table: pd.DataFrame, genome: Mapping[str, str],
                    replicons: Mapping[str, Replicon]) -> list[ProteinRecord]: ...
```
`src/giap/io/novel_tables.py` (I/O): `write_novel_orfs`, `read_novel_orfs`.
Schema: `src/giap/schemas/tables/novel_orfs.schema.json`.

## Rules

1. **Columns** (§9.3, then ADR-0012 inserted after `promoter_score`): `orf_id`, `class`, `contig`,
   `strand`, `frame`, `start_longest`, `start_best_ostir`, `end`, `length_aa_longest`,
   `length_aa_best`, `start_codon_best`, `n_peptides`, `peptides`, `best_peptide_q`, `ostir_tir`,
   `promoter_score`, `promoter_supported`, `terminator_downstream`, `terminator_conf`,
   `terminator_distance_nt`, `rbs_supported`, `class_fdr`, `called`, `rank`, `overlapping_gene`,
   `flags`. Columns the pass-1 pipeline cannot fill yet are written as `NA` (second-round scores
   before M7-08; evidence and rank before M7-09) — never omitted.
2. **Coordinates.** `start_*` are the 1-based position of the start codon's **first base in
   transcription direction** (on `−`, the highest coordinate of the codon) and `end` the 1-based
   position of the stop codon's last base in transcription direction (decision 1). `length_aa_*`
   exclude the stop. When `start_best_ostir` is NA (no OSTIR score), `length_aa_best` and
   `start_codon_best` describe `start_longest` and the flag says so (`ostir_missing`, decision 2).
3. **`class_fdr`**: a float (`%.6g`) or the literal `unestimated` (§9.3 "float or `unestimated`");
   the schema declares it as a string with that pattern.
4. **Order.** §6.5: contig (FASTA order), then the ORF's lowest coordinate, end, strand, `orf_id`.
   `rank` (M7-09) is a column, not the row order.
5. **`novel_orfs.gff3`**: every ORF row (called or not), feature type `CDS`, source `GIAP`, phase
   `0` (ADR-0005: CDS carries a phase), interval = best start (or longest) to stop, attributes
   `ID=<orf_id>`, `class`, `called`, `rank`, `n_peptides`, `start_longest`, `start_codon`, the four
   ADR-0012 values, `Note=<class>`. `##sequence-region` for every replicon, `Is_circular=true`
   region lines for circular ones (M1-05).
6. **`novel_orfs.faa`**: two records per ORF with a start (decision 3): `>NOVEL_…|longest` and
   `>NOVEL_…|best_ostir` (the second only when `start_best_ostir` exists and differs; when equal,
   one record `>NOVEL_…|longest,best_ostir`). Translation from the genome with table 11, the
   initiator codon read as `M` (as Bakta and `makedb` do), no `*`. ORFs without a start codon are
   not written.
7. **`annotation_patch.gff3`**: **called** ORFs only. New ORFs (`intergenic`, `cds_outframe`,
   `antisense`) as `CDS` with `ID=<orf_id>`, `product=hypothetical protein (GIAP novel ORF)`, best
   start → stop. `nterm_extension`: `ID=<locus_tag>;locus_tag=<locus_tag>;Note=nterm_extension`
   with the new start and the gene's own end. `cterm_extension`: `Note=cterm_extension`, the gene's
   start and the new end — **written only once M7-16 decides how `makedb` translates the read-through
   stop**, because M3-23 rejects a patch CDS with an internal stop (decision 4). `pseudogene` and
   `sorf` rows are not patched (they are already in Bakta's model).
8. **Round trip.** `giap makedb --patch annotation_patch.gff3` must produce, for every new ORF, a
   protein identical to its `|best_ostir` (or `|longest`) record, and for every N-terminal
   extension a `<locus_tag>_patch` record (M3-23's collision rule).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no ORFs | header-only TSV; GFF3 with pragmas only; empty FASTA; patch GFF3 with pragmas only |
| no called ORF | patch has pragmas only; the other three files list every candidate |
| ORF crossing the origin | GFF3 `end > length` (M1-05 convention); FASTA translated across the origin |
| ORF with `end_3p_open` (linear end) | GFF3 interval to the last whole codon, `partial=true` attribute; FASTA without a stop |
| ORF without a start codon | TSV row, GFF3 feature from the first peptide (`start_codon=NA`), no FASTA record, never patched |
| `start_best_ostir == start_longest` | one FASTA record (rule 6) |
| `cterm_extension` called | TSV/GFF3 rows; no patch row until M7-16 (rule 7) |
| two called ORFs with one ID | impossible (ID from the unique longest interval) — assert |

## Examples

| # | Input (ORF of M7-06 example 1, called, OSTIR best at the GTG 60) | Expected |
|---|---|---|
| 1 | TSV row | `orf_id NOVEL_c1_43_p`, `start_longest 43`, `start_best_ostir 61`, `end 123`, `length_aa_longest 26`, `length_aa_best 20`, `start_codon_best GTG`, `n_peptides 2`, `class_fdr unestimated`, `called true` |
| 2 | FASTA | `>NOVEL_c1_43_p\|longest` 26 aa starting `M`; `>NOVEL_c1_43_p\|best_ostir` 20 aa starting `M` (GTG read as M) |
| 3 | GFF3 | `c1 GIAP CDS 61 123 . + 0 ID=NOVEL_c1_43_p;class=intergenic;called=true;start_longest=43;…` |
| 4 | patch | `c1 GIAP CDS 61 123 . + 0 ID=NOVEL_c1_43_p;product=hypothetical protein (GIAP novel ORF)` |
| 5 | M7-06 example 6 (`nterm_extension` of G at TTG 249, called) | patch `c1 GIAP CDS 250 900 . + 0 ID=G;locus_tag=G;Note=nterm_extension` |
| 6 | M7-06 example 4 (`−` mirror) | `start_longest 958`, `end 878`; GFF3 interval 878–958 on `−` |
| 7 | `giap makedb --patch` on examples 4 and 5 | `strain_db.faa` holds `NOVEL_c1_43_p` (= the `best_ostir` protein) and `G_patch` (the extended G protein, starting `M`); `db_id_map.tsv` rows with kind `patch` |
| 8 | `read_gff3(write(...))` | the same features (round trip) |

## Acceptance criteria

- [ ] `tests/unit/novel/test_outputs.py`: examples 1–6, 8 and every degenerate row.
- [ ] `tests/integration/test_patch_makedb_roundtrip.py`: example 7 on a simulated genome with
  planted ORFs and extensions (§7.6.6 "patch → `makedb` round trip").
- [ ] `frictionless validate` on `novel_orfs.tsv`; both GFF3 files parse with the M1-05 reader and
  pass a GFF3 validity check (§9.5); `schemas/datapackage.json` lists the resource.
- [ ] TE7: the four files are byte-identical across reruns.

## Hardening owed (pass 2)

- Golden `novel_orfs.tsv` for the 50-kb fixture once it plants ORFs with ≥ 2 peptides (M7-12/15).
- The `cterm_extension` patch rows (M7-16).

## Questions for the lead

1. `start_*` and `end` on the minus strand: §9.3 gives names, not a convention. **Proposed
   default:** transcription-direction positions (start = first base of the start codon, end = last
   base of the stop), so `start > end` on `−`; the GFF3 uses the usual lowest/highest pair. The
   alternative (lowest/highest in the TSV) loses which end is the start.
2. When OSTIR has no score (no second-round layer, or OSTIR found no row at any candidate): report
   the longest start as "best" with flag `ostir_missing` (proposed), or leave the `*_best` columns
   NA? **Proposed default:** NA in `start_best_ostir` and `ostir_tir`; `length_aa_best` and
   `start_codon_best` fall back to the longest start and carry `ostir_missing`, so ranking (M7-09)
   still has a start codon.
3. FASTA headers for the two start variants. **Proposed default:** `|longest` / `|best_ostir`
   suffixes, merged when equal; `makedb --patch` never reads this file, so the suffix cannot leak
   into search IDs.
4. `cterm_extension` in the patch: M3-23 raises on an internal stop. **Proposed default:** no patch
   row until M7-16 decides whether `makedb` writes the read-through residue as `X`, `U` or `W`.

## Lead decisions

1. **Questions 1–4:** the proposed defaults are accepted.
   - **Minus-strand positions:** `start_*`/`end` follow the transcription direction (on `−`, start > end). The `novel_orfs` schema field descriptions must say so explicitly, because every other table is lowest/highest. The GFF3 uses lowest/highest.
   - **No OSTIR score:** the `*_best` OSTIR columns are NA; `length_aa_best`/`start_codon_best` fall back to the longest start, with `ostir_missing`.
   - **FASTA header suffixes:** `|longest` and `|best_ostir`, merged when equal.
   - **`cterm_extension`:** no patch row until M7-16.
2. **From the M7-06 review (#161):** `overlapping_gene` is now set for antisense ORFs (the opposite-strand gene). I/L twins share the contaminant ban and the multi-locus count.
3. **The test author's choices** (commit 7c4dc0f) are accepted as the contract.
   - **`orf_table(orfs, scored, *, genome, replicons)`:** `scored` is keyed by `orf_id` and holds the `minimal_call` columns plus optional overrides.
   - **`length_aa_longest` without a start:** counted from the first peptide to the stop. The implementer re-finds the first peptide in the translated frame (I/L folded). `Orf` keeps the M7-06 field set, because M7-06's tests pin its field order.
   - **Cell formats, GFF3 attributes, patch rows and the schema** are as listed in the test docstrings.
   - **`datapackage.json`** is created here with the `novel_orfs` resource only. M7-03 (#160) creates the same file with `pileup` and `variants`, and the lead merges the two resource lists when the second lands.
4. **After the code review (PR #165):**
   - **B1, fixed by the lead.** This revises the length rule of decision 3.
     - Re-finding the first peptide by exact text fails for tolerant-route, frameshift and junction peptides: it gives a short ORF, or raises and stops the step.
     - `Orf` therefore gains a trailing field, `first_peptide`. It holds the lowest 0-based coordinate of the most upstream peptide's first codon, taken from its mapped hit, and `reconstruct` always sets it.
     - `length_aa_longest` without a start is counted from that position.
     - The text search remains only for `Orf`s built by hand with `first_peptide = None`.
     - The M7-06 field-order test now lists the new field.
   - **B2, deferred to pass 2 (M7-15).** A − strand ORF with no start codon, whose stop sits at genome 0–2 on a linear contig, is treated as 3'-open when it is open only at its 5' end. Its GFF3 feature is then one codon short.
     - The table cannot tell the two cases apart without saying which end is open.
     - That needs a vocabulary addition to §9.4 (for example `open_3p`), which is a spec change and so an ADR.
     - The case is rare: a short linear contig spanned by one ORF.
   - The rest of the review findings go to docs/backlog.md.
