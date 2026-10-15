# M4-10: Protein covariates: observable peptides, GRAVY, pI, MW, TM segments

Milestone M4 · track G (gene features) · depends on: none (reads C2 columns `aa`, `nt`,
`feature_type`) · estimate ≈ 220 code lines (+ ≈ 250 test lines) · **can start now**

## Spec sections implemented

- §7.3.3 Covariates: "log1p(`n_observable_peptides`): unique tryptic 7–30 aa under the configured
  enzyme; log(`length_aa`); GRAVY; pI (Biopython ProtParam); `tm_segments`: Kyte–Doolittle, window 19,
  mean ≥ 1.6, runs merged, capped at 10 (G6)."
- §7.3.3 Pseudogenes: "Covariates come from Bakta's amino-acid sequence when present, otherwise from
  the in-frame translation up to the first stop."
- §7.3.4: "observable peptides are counted on the strain gene sequence with the configured enzyme, 0
  missed cleavages and length 7–30."
- §5.2 `abundance.enzyme: trypsin_p` "(cleaves before P too; DIA-NN default) | trypsin";
  `abundance.observable_len: [7, 30]`; `annotation.translation_table: 11`.
- §9.3 `genes.tsv` protein columns: `mw_da`, `pi`, `gravy`, `tm_segments`, `n_observable_peptides`;
  identity column `length_aa`.
- §2.4 G6.

## Public interface

`src/giap/genetable/covariates.py` (pure)
```python
Enzyme = Literal["trypsin_p", "trypsin"]

def digest(aa: str, enzyme: Enzyme, missed: int = 0) -> list[str]: ...
def n_observable_peptides(aa: str, enzyme: Enzyme, length: tuple[int, int] = (7, 30)) -> int: ...
def tm_segments(aa: str, window: int = 19, threshold: float = 1.6, cap: int = 10) -> int: ...
def protein_sequence(aa: str | None, nt: str, feature_type: str) -> str: ...   # pseudogene rule
@dataclass(frozen=True, slots=True)
class Covariates:
    length_aa: int; n_observable_peptides: int; gravy: float; pi: float; mw_da: float; tm_segments: int
def covariates(aa: str, enzyme: Enzyme, length: tuple[int, int]) -> Covariates: ...
def covariate_table(genes: pd.DataFrame, cfg: GiapConfig) -> pd.DataFrame: ...   # locus_tag + the six columns
```
- `digest`: cleave after K or R (`trypsin_p`), and not before P for `trypsin`; `missed > 0` adds
  concatenations of up to `missed + 1` consecutive fragments.
- `n_observable_peptides`: **distinct** peptide strings (as written; I and L distinct) of length in the
  range, 0 missed cleavages.
- GRAVY, pI, MW: Biopython `ProteinAnalysis` (`gravy()`, `isoelectric_point()`,
  `molecular_weight()` average masses), on the sequence with a leading/trailing `*` removed and
  letters outside the 20 standard amino acids (X, B, Z, J, U, O) removed; `length_aa` counts the
  sequence before removal (without `*`).
- `tm_segments`: per-position Kyte–Doolittle values (unknown letters 0), sliding mean over windows of
  19; a segment is a maximal run of consecutive windows with mean ≥ 1.6 (overlapping windows merge);
  count capped at 10; sequences shorter than 19 → 0.
- `protein_sequence`: `aa` when present; else (pseudogene) table-11 translation of `nt` from
  the first codon, first codon read as M, stopping before the first stop codon.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `digest("MAGLLKAAGPEPTIDERPAAAGGSTLKAAAAAAAK", "trypsin_p")` | `["MAGLLK", "AAGPEPTIDER", "PAAAGGSTLK", "AAAAAAAK"]` |
| 2 | same, `"trypsin"` | `["MAGLLK", "AAGPEPTIDERPAAAGGSTLK", "AAAAAAAK"]` (no cleavage before P) |
| 3 | `n_observable_peptides` of `"MAGLLKAAGPEPTIDERPAAAGGSTLK" + "W"*31 + "KAAAAAAAK"` | `trypsin_p`: 3 (`AAGPEPTIDER` 11, `PAAAGGSTLK` 10, `AAAAAAAK` 8; `MAGLLK` 6 and the 32-aa W-run excluded); `trypsin`: 2 |
| 4 | `n_observable_peptides("AAAAAAAKAAAAAAAKGGGGGGGR", "trypsin_p")` | 2 (distinct: `AAAAAAAK` twice counts once) |
| 5 | `covariates("MKTAYIAKQR", …)` | `gravy = −0.78`, `pi ≈ 10.289`, `mw_da ≈ 1209.46`, `length_aa = 10` (equal to Biopython 1.88 to 1e-9) |
| 6 | `tm_segments("L"*25)` | 1 |
| 7 | `tm_segments("L"*20 + "K"*20 + "L"*20)` | 2 (no 19-window spans both L runs) |
| 8 | `tm_segments("L"*20 + "K"*3 + "L"*20)` | 1 (windows across the 3-K gap still have mean ≥ 1.6, so the runs merge) |
| 9 | `tm_segments("L"*18)`; `tm_segments(("L"*20 + "K"*20)*12)` | 0; 10 (12 capped) |
| 10 | pseudogene, `aa = NA`, `nt = "GTG" + "GCT"*5 + "TAA" + "GCT"*3 + "TAA"` | `protein_sequence` = `"MAAAAA"` (start read as M; stops at the internal TAA) |

## Edge cases

- Empty protein after cleaning → `length_aa = 0`, other covariates NaN (gene stays in the table;
  M4-13 treats NaN covariates by median imputation, flagged there).
- `aa` with a trailing `*` (some annotators) is accepted and stripped.
- Lower-case letters → `ValueError` (contract says upper case).
- `nt` length not a multiple of 3 (frameshifted pseudogene) → translate whole codons only.
- Enzyme value outside the enum → `ValueError` (config already validates it).

## Acceptance criteria

- [ ] `tests/unit/genetable/test_covariates.py`: examples 1–10.
- [ ] Property: `"".join(digest(s, e)) == s`; `n_observable_peptides` ≤ number of fragments;
  `tm_segments` ∈ [0, 10] and is 0 for any sequence made only of hydrophilic residues (D, E, K, R).
- [ ] `covariate_table` on the M4-07 1,000-gene layout returns one row per gene, sorted like C2,
  in < 2 s.

## Methods note (statistician sign-off)

Part of `docs/methods/detectability.md` (covariate definitions, the pseudogene rule, G6 TM proxy).

## Questions for the lead

1. "unique tryptic peptides": distinct within the gene (as above) or unique in the whole proteome?
   **Proposed default:** distinct within the gene: it is the iBAQ denominator (§7.3.4) and a
   per-protein property; proteome uniqueness is handled by evidence uniqueness (§7.2.6).
2. I/L kept distinct when counting distinct peptides? **Proposed default:** yes (as written, matching
   §5.3 "I/L are kept as written; equivalence is applied only in matching").

## Lead decisions (test-author questions)
1. `length_aa` of an all-`X` protein ("XXX") is 0, as the spec's edge-case example says (explicit example
   wins over the general "count before cleaning" rule).
2. Accepted: empty protein → NaN `gravy`, `pi`, `mw_da` (`n_observable_peptides`/`tm_segments` 0); non-standard
   letters score 0 in TM windows; leading and trailing `*` stripped; `covariate_table` keeps input row order.
