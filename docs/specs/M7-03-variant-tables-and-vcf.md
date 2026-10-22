# M7-03: `pileup.tsv.gz`, `variants.tsv` and `variants.vcf` (VCF 4.3, G14)

Milestone M7 · track V (variants) · **pass 1** · depends on: M7-02; merged M1-05 VCF writer
(`giap.io.vcf.write_vcf`), M1-04 TSV writer · estimate ≈ 330 code lines (+ ≈ 400 test lines) ·
not statistical

## Spec sections implemented

- §7.5.5 in full:
  - "`tables/pileup.tsv.gz`: one row per covered codon, with REF, ALT and IL distinct-peptide counts."
  - "`tables/variants.tsv`: sites with an ALT or `mixed` call, plus frameshift events."
  - "`tables/variants.vcf`, VCF 4.3: **Header:** `##fileformat=VCFv4.3`; `##source=giap-<version>`;
    `##reference=` with the genome SHA-256; `##contig` lines with lengths; `##INFO` and `##FORMAT`
    definitions. **Records:** CHROM = contig; POS = the lowest 1-based coordinate of the codon. REF
    = the codon bases on the forward strand. ALT = the minimal-change codon encoding each ALT
    residue, on the forward strand (reverse-complemented for minus-strand genes), multi-allelic
    when needed (G14). **INFO:** `GENE`, `AA_POS`, `AA_REF`, `AA_ALT`, `CALL`, `CONF`,
    `NT_CHANGES`, `ALT_TYPE`, `FLAGS`, `CLASS_ERR`. **FORMAT:** `GT:AD:DP`. GT is haploid: `1`,
    `2`, … for ALT, `.` for `mixed`. AD = REF then each ALT, in distinct peptides. DP = total
    distinct peptides. **Sample column:** `config.sample_name`."
- G14 ("VCF ALT = the minimal-change codon for the ALT residue; the residue itself goes in INFO");
  Appendix B.7 (the record example); §6.5 ("VCF positions are 1-based within `[1, length]`");
  §9.2/§9.3 (`variants.tsv` columns `nt_changes`, `alt_codon_vcf`); §9.5 (Frictionless resources,
  `pysam` check); ADR-0005 item 1 ("a codon crossing the origin raises in the VCF writer (M7 decides
  how to report it)"; strict INFO/FORMAT; equal-length REF/ALT); ADR-0006 (value encoding).

## Public interface

`src/giap/variants/vcf_records.py` (pure)
```python
def minimal_alt_codon(ref_codon: str, alt_aa: str) -> tuple[str, int]: ...   # gene strand
def vcf_records(sites: pd.DataFrame, pile: pd.DataFrame, genome: Mapping[str, str],
                replicons: Mapping[str, Replicon]) -> tuple[list[VcfRecord], pd.DataFrame]: ...
INFO_FIELDS: tuple[VcfHeaderField, ...]
FORMAT_FIELDS: tuple[VcfHeaderField, ...]
```
`src/giap/io/variant_tables.py` (I/O)
```python
def write_pileup(path: Path, pile: pd.DataFrame) -> None: ...        # gzip, deterministic (mtime 0)
def write_variants(path: Path, sites: pd.DataFrame) -> None: ...
def read_variants(path: Path) -> pd.DataFrame: ...
```
Schemas: `src/giap/schemas/tables/pileup.schema.json`, `variants.schema.json`
(`coordinates: 1-based-closed`; enums `call`, `confidence`, `alt_type`; `flags` enum = the §9.4
Variants scope).

## Rules

1. **Minimal-change codon (G14).** Over the table-11 codons encoding `alt_aa` (for `J`: every I or
   L codon), choose the fewest nucleotide differences from `ref_codon`; ties go to the
   lexicographically smallest codon (decision 1). Returns `(codon, n_changes)`; `n_changes` is
   `nt_changes`. For swap and 2↔2 blocks (two codons), apply it per codon to the two ALT residues.
2. **What enters the VCF** (decision 2): `substitution` sites (3-nt REF/ALT) and equal-length
   two-codon events — `swap` and 2↔2 `block`s such as GT↔AS (6-nt REF/ALT at the first codon).
   1↔2 and 2↔1 blocks and `frameshift` events are **not** in the VCF (the M1-05 writer requires
   equal-length alleles, and §7.5.2 already keeps frameshifts out); their `alt_codon_vcf` is NA.
3. **Strand.** REF = the forward-strand bases at `[codon, codon + 3·k)`. On a `−` gene, ALT = the
   reverse complement of the gene-strand minimal codon(s); `alt_codon_vcf` stores the forward-strand
   ALT.
4. **Origin-crossing codons** (circular replicons, a codon at `length − 1` or `length − 2`): not
   representable (ADR-0005); the site stays in `variants.tsv` with `alt_codon_vcf` NA and the flag
   `vcf_unrepresentable` (decision 3; a new Variants flag value).
5. **Records.** One record per codon and VCF length class (3-nt, 6-nt) holding every ALT of that
   class at the codon, alleles ordered by `n_alt` descending, then `alt_aa`. Two records may share
   a POS (a substitution and a swap starting at one codon); that is legal VCF.
6. **INFO.** `GENE` (String; omitted when NA), `AA_POS` (Integer; omitted when NA), `AA_REF`
   (String), `AA_ALT` (String, Number=A), `CALL` (String), `CONF` (String, Number=A), `NT_CHANGES`
   (Integer, Number=A), `ALT_TYPE` (String, Number=A), `FLAGS` (String, Number=., `|`-joined per
   allele as B.7; omitted when no allele has a flag — ADR-0005 forbids empty lists), `CLASS_ERR`
   (String, Number=A: `%.6g` or `unestimated`). Descriptions are fixed strings in `INFO_FIELDS`.
7. **FORMAT.** `AD` = `n_ref`, then each allele's `n_alt`; `DP` = `n_ref + n_il + Σ n_alt` (every
   distinct eligible peptide at the codon, decision 4). `GT`: call `ALT` with one allele → `1`;
   call `ALT` with several alleles → the 1-based index of the allele with the largest `n_alt`, `.`
   on a tie (decision 5); `mixed` → `.`.
8. **FILTER** is `PASS` for every record (decision 6); confidence lives in `CONF`.
9. **Header** as §7.5.5 via `write_vcf` (`##reference=sha256:<hex>`, ADR-0005), sample =
   `config.sample_name` resolved by `giap run` (M7-04).
10. **TSV writers.** 1-based closed `codon_start` (`codon + 1`, wrapped into `[1, length]`), `NA`,
    `%.6g`, `;`-lists, §6.5 order; `pileup.tsv.gz` gzip with `mtime=0` so TE7 (byte-identical
    reruns) holds. `variants.tsv` column order exactly §9.3; columns filled later
    (`class_error_rate`, `class_error_n`: `unestimated`/NA until M7-10; `flags` empty until M7-05)
    are present from the start.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no sites | header-only `variants.tsv`; a VCF with header lines and no records (valid for pysam) |
| empty pileup | header-only gzip TSV |
| ALT residue `*` (cannot come from a peptide) | `minimal_alt_codon` still defined (stop codons), used only by the exhaustive table test |
| `J` ALT | minimal over I and L codons; `AA_ALT=J` in INFO |
| tie between minimal codons | lexicographic (example 3) |
| site on a `−` gene | ALT reverse-complemented (example 4) |
| origin-crossing codon | rule 4 |
| block 1↔2 / frameshift | `variants.tsv` only (rule 2) |
| contig ID not VCF-safe | the M1-05 writer raises; preflight's safe-ID mapping (ADR-0005) makes this unreachable |
| sample name `auto` reaching the writer | `ValueError` (it must have been resolved by `giap run`) |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `+` gene, codon at 0-based 1000 = `GCT` (A), ALT `S`, 3 distinct ALT peptides, no REF; flags `mass_equivalent`, `snv` (M7-05) | the B.7 record: `chr 1001 . GCT TCT . PASS GENE=…;AA_POS=45;AA_REF=A;AA_ALT=S;CALL=ALT;CONF=high;NT_CHANGES=1;ALT_TYPE=substitution;FLAGS=mass_equivalent\|snv;CLASS_ERR=unestimated GT:AD:DP 1:0,3:3` |
| 2 | `minimal_alt_codon("GTT", "J")` (V → I or L) | ATT (I) and CTT (L) are both 1 change → **`("ATT", 1)`** (lexicographic) |
| 3 | `minimal_alt_codon("AAA", "N")` | AAC and AAT both 1 change → **`("AAC", 1)`** |
| 4 | `−` gene, forward bases at the codon `AGC` (gene strand `GCT`, A), ALT `S` | gene-strand ALT `TCT`; VCF REF **`AGC`**, ALT **`AGA`** (revcomp of TCT); POS = codon + 1 |
| 5 | `mixed`: 1 REF, 2 ALT `D` | `GT` `.`, `AD` `1,2`, `DP` 3 |
| 6 | `+` codon `AAT` (N); ALT `D` (2 peptides) and `E` (1), no REF | one record, `ALT=GAT,GAA` (D: GAT 1 change; E: GAA and GAG 2 changes, tie → GAA); `GT` `1`; `AD` `0,2,1` |
| 7 | ALT `D` (1) and `E` (1), no REF | `GT` `.` (tie) |
| 8 | swap `SD → DS` at two `+` codons `TCT GAT` | one 6-nt record: REF `TCTGAT`; codon 1 → D: GAT (2 changes; GAC 3) ; codon 2 → S: AGT and TCT (2 each, tie → AGT); ALT **`GATAGT`**, `NT_CHANGES=4` |
| 9 | block `N → GG` (1↔2) | `variants.tsv` row only, `alt_codon_vcf` NA; no VCF record |
| 10 | codon at 0-based 999 on a circular 1,000-nt contig | `variants.tsv` row with flag `vcf_unrepresentable`; no VCF record |
| 11 | the written VCF | `pysam.VariantFile` reads it without warnings; INFO values round-trip; `bcftools`-free |

Examples 3, 6 and 8 exercise the tie-break; the test-author re-derives them with the brute-force
oracle of the exhaustive table rather than trusting this table.

## Acceptance criteria

- [ ] `tests/unit/variants/test_vcf_records.py`: examples 1–10, every degenerate row.
- [ ] §7.5.6 unit items: "minus-strand codon conversion" and "a VCF round trip with
  `pysam.VariantFile`" (example 11); the exhaustive `minimal_alt_codon` table over 64 codons ×
  21 residues (20 + `*`) plus `J` against a brute-force oracle.
- [ ] Property: for every record, applying ALT to the genome and translating on the gene's strand
  gives `AA_ALT`; REF matches the genome (`write_vcf(..., sequences=genome)` passes).
- [ ] `frictionless validate` passes on both TSVs; `schemas/datapackage.json` lists the two
  resources; `pileup.tsv.gz` is byte-identical across two runs with the same seed (TE7).

## Hardening owed (pass 2)

- Golden `variants.tsv` and `variants.vcf` for the 50-kb fixture (M7-15).
- The backlog item "declare GIAP's VCF FILTER values with descriptions": closed by decision 6
  (PASS only) unless the lead prefers a `candidate` filter; M7-15 removes the generic `##FILTER`
  path from the writer if no FILTER other than PASS is ever written.

## Questions for the lead

1. G14 tie-break among equally minimal ALT codons. **Proposed default:** lexicographically smallest
   codon on the gene's strand — deterministic, codon-usage-free. (The alternative, the strain's most
   used codon, makes the VCF depend on codon statistics.)
2. Which `alt_type`s enter the VCF? **Proposed default:** substitutions and equal-length two-codon
   events (swap, 2↔2 blocks); 1↔2 / 2↔1 blocks stay in `variants.tsv` only, like frameshifts.
   Representing them needs unequal-length alleles, which ADR-0005 forbids and which would read as
   indels to every VCF consumer.
3. Origin-crossing codons (ADR-0005 left this to M7). **Proposed default:** omit from the VCF and
   flag the `variants.tsv` row `vcf_unrepresentable` (new Variants flag value, added through the
   schema). Splitting a codon across the origin into two records would misstate a single residue
   change.
4. Does `DP` count IL-only peptides? **Proposed default:** yes — DP is "total distinct peptides"
   covering the codon, and an I/L-ambiguous peptide covered it.
5. `GT` with two or more ALT alleles and no REF: the spec defines only "`1`, `2`, … for ALT".
   **Proposed default:** the best-supported allele's index, `.` on a tie.
6. VCF FILTER values (backlog): **Proposed default:** every record `PASS`; `CONF` carries
   high/candidate. A `candidate` FILTER would make single-peptide calls disappear from filtered VCF
   views, which hides exactly the calls `class_error_rate` is meant to qualify.
