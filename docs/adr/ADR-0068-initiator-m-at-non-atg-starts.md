# ADR-0068: Initiator M at annotated non-ATG starts
- Status: accepted-by-khoa (Accepted by Khoa, 2026-10-28, relayed by the lead; Q-K56 / M11-20
  Q-K1: option A, annotation-only, exact route only)
- Date: 2026-10-28
- Spec sections: §2.2 improvement 6 and §15 ("No: … non-ATG starts, Met removal"), §7.2.1
  (six-frame index; storage), §7.2.2 (route 1, exact), §9.4 (Placement flags), §11.6 TE1/TE8
  proxy wording; docs/specs/M11-20-initiator-methionine.md
- Needs Khoa: answered 2026-10-28 (overrides a declined spec item)

## Context
Bakta writes `M` as the first residue of every complete CDS whatever its start codon, while the
six-frame index reads an annotated GTG start as `V` and TTG as `L`. The Met-retaining N-terminal
peptide of a non-ATG gene therefore had no exact hit at its own gene: 7–9 aa peptides were lost,
longer ones became tolerant `strain_difference` rows with an ALT `M` at the start codon, and a
chance exact copy at a real ATG elsewhere became its only locus (an unestimated false novel ORF,
M7-13 decision 14). GIPPA_SPEC §2.2 / §15 decline a start catalogue and the TE1 proxy accepted
the miss.

## Decision
1. A *start site* is the first codon of a `cds` or `sorf` feature whose Bakta `aa` begins with
   `M` and whose genome codon (gene strand, wrapped on a circular replicon) is in
   `annotation.start_codons` and is not `ATG`. Only the run's own annotation is used, as
   ADR-0058 does for annotated stops; no catalogue.
2. Route 1 credits query position 0 (`M`) at a start cell: window `[p, p + k)` in the target and
   `shuffle` decoy, the mirrored `[c − k + 1, c + 1)` with `c = n − 1 − p` in the `reverse`
   decoy, so decoys are searched identically (§7.2.5). Credited rows are ordinary exact rows.
   The tolerant route is unchanged.
3. Credited target rows carry the new §9.4 Placement flag `initiator_m` (`OUTPUT_SCHEMA_VERSION`
   1.2). Variants give the credited codon no category (op `initiator`).
4. The six-frame index is format 4 (start sites and cells); formats 2 and 3 load with no start
   cells and the WARNING event `mapping.index_without_start_cells`.
5. The TE1/TE8 proxies no longer exclude N-terminal products of non-ATG genes from the
   denominator; their count stays as a diagnostic.

## Consequences
- Recall of non-ATG N-terminal peptides at their own gene; they leave the novel candidate set.
- About 0.25 % more chance windows for M-initial peptides, matched by the mirrored decoy cells.
- Follow-ups (backlog): credit in the tolerant route (A′) and novel-ORF non-ATG starts (A″).
