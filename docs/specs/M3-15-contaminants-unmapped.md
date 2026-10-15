# M3-15: Unmapped bucket and contaminants

Milestone M3 · depends on: M3-10 (rejected placements), M3-12/13 (protein leftovers), M1-03,
M1-06 (`is_contaminant`) · creates `src/giap/mapping/contaminants.py` · pure Python for peptides
(`pyahocorasick`) + **external tool BLAST+ `blastp`** for proteins (same env as M3-13) · estimate
≈ 250 code lines

## Spec sections implemented

- §7.2.8 "Sequences that land nowhere are matched against the universal contaminant library
  (Frankenfield 2022, 381 entries), downloaded by `giap setup --contaminants` (vendor it only if
  its licence allows): peptides by exact substring (I = L); proteins by `blastp` at ≥ 90 %
  identity over ≥ 50 % of the query; engine protein IDs carrying a contaminant prefix also count
  as contaminants. Matches are labelled `contaminant`. Everything else is `unmapped`, with a
  reason: `no_hit`, `rejected_placement_fdr`, `below_identity` or `too_short`. Plass is not used."
- §9.3 `mapping/unmapped.tsv`: `sequence_id`, `sequence`, `kind` enum{peptide, protein}, `label`
  enum{contaminant, unmapped}, `reason` enum{no_hit, rejected_placement_fdr, below_identity,
  too_short, contaminant_id}, `contaminant_match`.
- §11.6 TE5 "trypsin and keratin peptides land in the contaminant bucket" (in-silico proxy: digests
  of porcine trypsin and human keratins added to a peptide table).
- Config `input.contaminant_prefixes`.

The library itself is downloaded by `giap setup --contaminants` (M2-09, §7.9.1) to
`DIR/contaminants/contaminants.fasta` with a `db_manifest.json` entry; this issue reads the FASTA path from there (tests use
a small fixture FASTA with a few library-like entries).

## Public interface

```python
Reason = Literal["no_hit", "rejected_placement_fdr", "below_identity", "too_short", "contaminant_id"]

@dataclass(frozen=True, slots=True)
class UnmappedRow:
    sequence_id: str; sequence: str; kind: Literal["peptide", "protein"]
    label: Literal["contaminant", "unmapped"]; reason: Reason; contaminant_match: str | None

def peptide_reason(has_target_hit: bool, accepted: bool, length: int, tolerant_min_len: int) -> Reason | None
def match_contaminant_peptides(peptides: Sequence[str], library: Sequence[ProteinRecord]) -> dict[str, str]: ...
def blastp_cmd(query: Path, library: Path, threads: int) -> list[str]: ...
def match_contaminant_proteins(blastp_tsv: Path, min_identity: float = 0.90,
                               min_coverage: float = 0.50) -> dict[str, str]: ...
def unmapped_table(peptides: pd.DataFrame, placements: pd.DataFrame, protein_state: pd.DataFrame,
                   pep_matches: Mapping[str, str], prot_matches: Mapping[str, str],
                   contaminant_prefixes: Sequence[str]) -> list[UnmappedRow]: ...
```

**Reasons (peptides):** `rejected_placement_fdr` = has a target hit but `placement_q > α`;
`too_short` = no exact target hit and shorter than `mapping.tolerant.min_len` (no tolerant
search was possible); `no_hit` = searched by every applicable route, no target hit within the
budget. **Proteins:** `below_identity` = alignments exist but none passed the thresholds;
`no_hit` = no alignment from either route.
**Labels:** a library match (peptide substring, I = L; protein blastp ≥ 90 % identity over ≥ 50 %
of the query) → `contaminant`, `reason` keeps the mapping reason, `contaminant_match` = library
ID; otherwise, every `protein_ids` entry carrying a contaminant prefix → `contaminant`,
`reason = contaminant_id`, `contaminant_match` = the first such ID; otherwise `unmapped`.
`sequence_id` = the peptide sequence for peptides, the protein ID for proteins. Sorted by
`kind`, then `sequence_id`.

## Examples

| # | Input | Expected row |
|---|---|---|
| 1 | peptide `VATVSLPR` (porcine trypsin), no genome hit, library contains trypsin | `contaminant`, `too_short`, `contaminant_match` = trypsin entry ID |
| 2 | peptide `SLNNQFASFIDK` (human keratin 1), no target hit, library contains KRT1 | `contaminant`, `no_hit`, match = KRT1 ID |
| 3 | same keratin peptide written with `I` vs library `L` | still matched (I = L) |
| 4 | peptide with a target hit at placement q 0.03 | `unmapped`, `rejected_placement_fdr` |
| 5 | peptide with no hit, `protein_ids = CON__P02768` | `contaminant`, `contaminant_id`, match `CON__P02768` |
| 6 | 8-aa peptide with no exact hit, not in the library | `unmapped`, `too_short` |
| 7 | 14-aa peptide with no exact or tolerant hit | `unmapped`, `no_hit` |
| 8 | protein with a 70 %-identity tblastn alignment only, not in the library | `unmapped`, `below_identity` |
| 9 | protein with blastp 92 % identity over 60 % of the query to bovine albumin | `contaminant`, reason from mapping (`no_hit` or `below_identity`), match = albumin ID |
| 10 | protein with blastp 95 % over 40 % | not a contaminant (coverage < 0.5) |

## Edge cases

- Peptides placed and accepted never enter the bucket, even if their `protein_ids` carry a
  contaminant prefix (Q1).
- A peptide matching several library entries: `contaminant_match` = the lowest library ID
  (deterministic).
- Library absent (setup not run): peptides and proteins get `unmapped` with their reasons; the run
  records `layer_missing:contaminants` (§9.4) and preflight already warned.

## Acceptance criteria

- [ ] `tests/unit/test_contaminants.py`: examples 1–10 (fixture library FASTA in
  `tests/data/contaminants_mini.faa` with a trypsin, a keratin and an albumin entry).
- [ ] TE5 proxy (in-silico, `tests/integration/test_te5_proxy.py`): tryptic digests of porcine
  trypsin and human keratins added to a simulated peptide table → every such peptide that does not
  map to the genome is labelled `contaminant`.
- [ ] `tests/contract/test_blastp_contract.py` on a recorded blastp output with `PROVENANCE.md`.

## Questions for the lead

1. Should a peptide that maps to the genome but whose `protein_ids` are all contaminant IDs be
   moved to the bucket? **Proposed default:** no; the spec applies contaminant checks to
   "sequences that land nowhere"; such peptides keep their placement (rare, and visible through
   `is_contaminant` in the work table).
2. `too_short` is not defined by the spec. **Proposed default:** as above (no exact hit and too
   short for tolerant search), which separates "could not be searched tolerantly" from "searched,
   nothing found".

## Lead decisions (test-author questions)
1. `peptides` is the M1-06/M1-12 aggregated frame with `protein_ids` as a native list of strings (as M1-12
   writes it; M4-06 decision 1) — `;` joining only in TSV writers.
2. Several passing contaminant-library hits for one protein: deterministic tie rule — highest bitscore, then
   highest identity, then highest coverage, then lowest subject ID (lexicographic).
3. Accepted: `placements` = M3-10 per-peptide frame (presence = target hit); `protein_state` columns
   `protein_id, sequence, assigned, has_alignment`; `unmapped_table(..., tolerant_min_len=10)` keyword; blastp
   uses tblastn's 15 outfmt-6 columns; thresholds per HSP, inclusive, identity = pident/100, coverage =
   (qend−qstart+1)/qlen; a library match beats `contaminant_id` and keeps the mapping reason; every
   `protein_ids` entry must carry a prefix and the first ID is the match; a protein whose own ID has a
   contaminant prefix → `contaminant_id`; accepted placements never enter the unmapped bucket; decoy-only
   peptides → `no_hit` or `too_short`.
4. The implementer records `tests/data/recorded/blast/2.16.0/contaminants/` with the pinned BLAST image
   (real blastp output; lead rule for tool wrappers).
