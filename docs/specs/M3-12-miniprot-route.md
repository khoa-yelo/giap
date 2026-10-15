# M3-12: Route 2a, proteins > 40 aa with miniprot

Milestone M3 · depends on: M3-01 (`GeneModel`), M1-03 (protein FASTA) · creates
`src/giap/mapping/miniprot.py`, `src/giap/mapping/protein_assign.py`; reuses or creates
`workflow/envs/miniprot.yml`, `workflow/modules/local/miniprot/main.nf` · **external tool:
miniprot 0.18** (conda `miniprot=0.18` or Apptainer image by digest; MIT licence, called as a
separate program, never imported) · estimate ≈ 300 code lines

## Spec sections implemented

- §7.2.3 "**Alignment.** miniprot as in §7.1.1. Parse GFF attributes `Identity`, `Positive`,
  `Rank`, `Target`, `Frameshift` and `StopCodon`. Keep the flags `fs` and `st`."
  "**Assignment to a Bakta CDS** requires all of: identity ≥ 0.80; query coverage ≥ 0.80;
  reciprocal overlap with the CDS ≥ 0.80 (overlap / CDS length and overlap / alignment span)."
  "**Similarity class:** `identical` (identity 1.0 and coverage 1.0), `near` (≥ 0.95) or
  `divergent` (0.80–0.95). Below the thresholds a protein goes to the next route."
  "**Multiple loci.** Secondary alignments (Rank > 1) are kept as additional loci."
- §7.1.1 gotcha: "`-T 11` is an indexing option. Build the index with
  `miniprot -t N -T 11 -d genome.mpi genome.fa`, then align with
  `miniprot -t N -S -N 5 --outs=0.95 --gff genome.mpi proteins.faa`. miniprot indexes ORFs ≥ 30 aa
  only, which is why short and divergent proteins fall to tblastn."
- Decision 5 and 6 ("Proteins are assigned at ≥ 80 % identity over ≥ 80 % of their length").
  Config `mapping.proteins` (`min_identity`, `min_coverage`, `near_identity`, `miniprot_args`).
- §11.2 contract layer: "miniprot … outputs" parsed from a recorded native output.

**Relation to M2.** M2-01 adds `workflow/envs/miniprot.yml` and the container pin; M2 leaves the
wrapper, parser, process module and recording to M3 (M2-index decision 4). If M2-01 has not merged,
this issue adds the env file itself (env build via `sbatch scripts/slurm/build_envs.sbatch
miniprot`, never on the login node, never `pip install`).

## Public interface

```python
def miniprot_index_cmd(genome: Path, mpi: Path, threads: int) -> list[str]: ...
def miniprot_align_cmd(mpi: Path, proteins: Path, threads: int, args: str) -> list[str]: ...

@dataclass(frozen=True, slots=True)
class ProteinAlignment:
    protein_id: str; protein_len_aa: int; route: Literal["miniprot", "tblastn"]
    interval: Interval; rank: int
    identity: float; positive: float | None
    q_start: int; q_end: int                 # 0-based half-open on the protein (from Target)
    frameshift: int; stop_codon: int          # counts; > 0 → flags fs / st
    @property
    def query_coverage(self) -> float: ...   # (q_end − q_start) / protein_len_aa

def parse_miniprot_gff(path: Path, replicons: Mapping[str, Replicon],
                       protein_len: Mapping[str, int]) -> list[ProteinAlignment]: ...
def reciprocal_overlap(aln: Interval, cds: Interval, replicon: Replicon) -> float: ...  # min of both ratios
def similarity_class(identity: float, coverage: float, near: float) -> Literal["identical", "near", "divergent"]
def assign_alignments(alns: Sequence[ProteinAlignment], gm: GeneModel,
                      cfg: ProteinsConfig) -> tuple[pd.DataFrame, set[str]]: ...
#   → (protein_hits rows: §9.3 columns, gene NA where no CDS passes; leftover protein IDs)
```

- Parse `mRNA` lines (one per alignment); `Target=ID qs qe` is 1-based closed on the protein;
  coordinates converted once to `Interval` (minus strand from column 7). `##PAF` lines are ignored.
- A protein is **assigned** when at least one alignment passes identity, coverage and reciprocal
  overlap with a same-strand CDS or sORF (the CDS with the largest overlap wins; ties → lower locus
  tag). All passing alignments (any rank) become `protein_hits` rows. A protein with no passing
  alignment is a **leftover** for tblastn (M3-13).

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `miniprot_index_cmd(g, m, 8)` | `["miniprot", "-t", "8", "-T", "11", "-d", "m", "g"]` |
| 2 | `miniprot_align_cmd(m, p, 8, "-S -N 5 --outs=0.95")` | `["miniprot", "-t", "8", "-S", "-N", "5", "--outs=0.95", "--gff", "m", "p"]` |
| 3 | mRNA line `chr1 miniprot mRNA 101 400 … + . ID=MP000001;Rank=1;Identity=1.0000;Positive=1.0000;Target=GSIM_00005 1 99` with protein length 99 and CDS `[100, 403)` | identity 1.0, coverage 1.0, reciprocal overlap 300/303 = 0.990 → `class=identical`, `gene=GSIM_00005` |
| 4 | identity 0.97, coverage 0.92 on the same CDS | `class=near` |
| 5 | identity 0.85, coverage 0.81 | `class=divergent` |
| 6 | identity 0.79 | not assigned → leftover |
| 7 | alignment covering `[100, 250)` of a CDS `[100, 400)` (overlap/CDS = 0.5) | not assigned (reciprocal overlap < 0.8) → leftover |
| 8 | `Rank=2` alignment on a second copy, passing | an extra row with `rank=2` |
| 9 | `Frameshift=1;StopCodon=1` | flags `fs;st` |
| 10 | contract: recorded output for the 100-kb sim genome (`tests/data/recorded/miniprot/0.18/sim100k/`) | parses; ≥ 99 % of simulated proteins > 40 aa assigned `identical` to their own gene |

## Edge cases

- CDS crossing the origin: overlap computed on unwrapped intervals (shift the alignment by +L when
  needed); miniprot itself reports alignments on the linear sequence, so an origin-spanning protein
  may come back as two partial alignments → coverage < 0.8 → tblastn (whose HSP merge handles the
  wrap, M3-13).
- Proteins ≤ 40 aa never reach this route (preflight routes them to the peptide path).
- Empty GFF (no alignments): every protein is a leftover.
- `Identity` absent on a line → `ValueError` naming the line (contract drift).

## Acceptance criteria

- [ ] `tests/unit/test_miniprot.py`: examples 1–9 on hand-written GFF snippets.
- [ ] `tests/contract/test_miniprot_contract.py`: example 10 on the recorded output, with
  `PROVENANCE.md` (command, env lock or digest, date, host).
- [ ] Integration (`needs_conda` or `needs_apptainer`): live miniprot on the 100-kb sim genome
  gives the same parsed table as the recorded fixture.
- [ ] nf-test stub for the miniprot process writes a minimal valid GFF.

## Questions for the lead

1. An alignment passing identity and coverage on a region without a CDS (unannotated gene): what
   is reported? **Proposed default:** a `protein_hits` row with `gene = NA` (visible to Module 6),
   while the protein still counts as a leftover and goes to tblastn as the spec says.
2. `Identity` in miniprot is computed over aligned columns; coverage comes from `Target`. Confirm
   no extra penalty for frameshifts. **Proposed default:** none; `fs` is a flag only.

## Lead decisions (test-author questions)
1. Protein alignments are assigned to CDS or sORF only (spec). Pseudogenes are not assigned by this route;
   an alignment over a pseudogene stays a leftover (M3-14 decides novel loci). M3-index decision 3
   concerns peptide placements, not protein alignments. Partial CDS can be assigned.
2. Q1 and Q2 take the proposed defaults: a gene-NA row while the protein stays a leftover; frameshifts are
   a flag with no identity penalty.
3. Accepted test-author choices: builders, `ProteinAlignment` and the parser in `mapping/miniprot.py`;
   `reciprocal_overlap`, `similarity_class`, `assign_alignments` in `mapping/protein_assign.py`; 0-based
   half-open start/end; thresholds inclusive (≥), `near` = identity ≥ `near_identity`; largest overlap in bp,
   ties → lower locus tag; flags `"fs;st"`; row order independent of input order; leftovers are aligned
   protein IDs only (the caller adds unaligned proteins); a protein ID missing from `protein_len` raises.

## Lead decisions after code review (PR #47, changes requested)
1. Alternative starts: the parser reads the `cs` tag of the `##PAF` line preceding each alignment; a lone
   `*gtgM` / `*ttgM` (any `*[acgt]tgM`) at query position 0 counts as a match for identity, and
   `ProteinAlignment.alt_start` (bool) records it. Unit test with GTG and TTG starts; example 10 (≥ 99 %
   identical) must hold on the recorded SIM100K run.
2. Record the fixture `tests/data/recorded/miniprot/0.18/sim100k/` with the pinned miniprot image (seconds;
   implementers may add recorded fixtures) incl. PROVENANCE.md (command, image digest, input sha256).
3. A real nf-test stub test for `workflow/modules/local/miniprot/` (workflow/tests/miniprot.nf.test);
   `mapping.proteins.miniprot_args` reaches the process via `ext.args` (no hard-coded default in main.nf).
4. `assign_alignments` validates canonical intervals (`check_interval`: 0 ≤ start < L, unwrapped end ≤ start + L)
   and raises ValueError otherwise; M3-13 normalises chains before calling it.
5. miniprot runs on preflight's `proteins.long.faa`, whose IDs are already made safe (`_safe_protein_ids`);
   the parser raises ValueError if a `Target` ID contains `;` or `=` (defensive).
