# ADR-0067: Genetic codes 4, 11 and 25
- Status: accepted-by-khoa (Accepted by Khoa, 2026-10-28, relayed by the lead; supports codes 4
  and 25 alongside 11, superseding Q-K30; Bakta gets `--translation-table` only at 4 and 25)
- Date: 2026-10-28
- Spec sections: §4 (config default `translation_table: 11`), §7.0 (preflight), §7.1.1 (Bakta,
  miniprot, BLAST+), §7.2.1 ("translate all six frames of every contig with table 11"), §12.1
  ("table-11 translation"); docs/specs/M11-16-translation-tables-4-25.md
- Needs Khoa: answered 2026-10-28 (an interface and a protected default change:
  `annotation.translation_table` accepts 4 and 25, and D1 changes the
  `mapping.proteins.tblastn_args` default in `defaults.yaml`)

## Context
GIPPA_SPEC §7.2.1 fixes the six-frame translation to table 11, and #413 made preflight refuse any
other `annotation.translation_table` (Q-K30: code 11 only for 1.0). Mollicutes such as
*Mycoplasma* and *Spiroplasma* use code 4 (TGA = Trp), and candidate division SR1 and
Gracilibacteria use code 25 (TGA = Gly). Run at table 11, M. hyorhinis SK76 (PR #411) gave 801 of
1,200 CDSs ending in TGA, a median CDS protein of 80 aa and thousands of unmapped peptides with W.
Bakta 1.12.2 accepts `--translation-table` 4, 11 and 25. In all three tables only TGA differs.

## Decision
1. `annotation.translation_table` may be 4, 11 or 25; 11 stays the default. The pydantic model
   stays `int`; preflight `config.translation_table` refuses any other value with exit
   `INVALID_INPUT`, and the config schema lists the values in its description only (D7).
2. One pure module, `gippa.model.genetic_code`, holds the three tables (built from the NCBI
   TCAG strings, not from Biopython at run time): codon table, stops, minus-strand stops, recoded
   letters (Sec over TGA only under 11, Pyl over TAG in all three) and names. Every step takes
   the table from `annotation.translation_table` and calls it; no module keeps its own stop set.
3. GIPPA sets every tool's code from that one key. Preflight `config.tool_genetic_code` refuses a
   Bakta `--translation-table`, a miniprot `-T`, or a tblastn `-db_gencode` that differs from the
   table or is repeated. Until the D1 default change is made, the shipped `tblastn_args`
   default (which still says `-db_gencode 11`) is exempt; the tblastn command then derives
   `-db_gencode` from the table and the exemption is removed.
4. The `genome.genetic_code` warning (M11-14w) is evaluated at every supported table, at the same
   0.55 cut: ratio < 0.55 warns at table 11 (W11, names codes 4 and 25), ratio ≥ 0.55 warns at
   tables 4 and 25 (W4, the genome looks like code 11). It stays a `warn` (D5). The genome cannot
   tell code 4 from code 25, so the rule never warns about that pair.
5. Bakta gets `--translation-table <t>` only at tables 4 and 25. At table 11 no flag is passed
   (Bakta's default is 11), so code-11 Bakta argv and cached Bakta tasks are unchanged (D2).
6. A six-frame index or a precomputed Bakta JSON made with another table is refused, never
   reused (six-frame index format 3 records its table; format 2 is read as table 11, D3).
7. D1: the default `mapping.proteins.tblastn_args` becomes `"-seg no -evalue 1e-5"` and
   `-db_gencode` is always derived from the table. This changes protected
   `src/gippa/config/defaults.yaml`; the PR that makes it carries
   `Gate-Change-Approved: ADR-0067`.

This ADR overrides the "table 11" sentences of GIPPA_SPEC §7.2.1 and §12.1.

## Evidence (before → after: metrics, tests, simulations)
- `tests/unit/test_m11_16a_genetic_code.py`: the three tables equal Biopython's
  `CodonTable.unambiguous_dna_by_id`; every other table raises.
- `tests/unit/test_m11_16a_preflight.py`: preflight at 4, 11 and 25; `config.tool_genetic_code`
  case table; W11/W4 wording; at 4/25 the decision is the exact complement of 11.
- The M11-14w panel (61 official NCBI assemblies plus draft stress tests) has every code-4/25
  genome at ratio ≤ 0.470 (≤ 0.477 as drafts) and every code-11 genome at ≥ 0.679 (≥ 0.641), so
  the same cut gives no false W4 and no missed code-11 genome on it (D4; statistician review).
- Later PRs: per-table six-frame index and stop sites (16b), simulator codes (16c), tool argv and
  manifest (16d), downstream modules (16e), calibration arms at 4 and 25 with unchanged tolerances
  (16f), and a code-4 re-run of M. hyorhinis (16g).

## Consequences
- Table 11 output stays byte-identical apart from the documented changes (spec AC 1). The Bakta
  argv at table 11 does not change (D2), so cached code-11 Bakta tasks are kept.
- Codes 4 and 25 run end to end only once 16b and 16d are merged; between 16a and 16d, preflight
  accepts them but Bakta and the six-frame index still use table 11.
- One table per run: a plasmid with a different code from its chromosome is not supported.
- Choosing between code 4 and 25 is the user's, from taxonomy; GIPPA cannot check it.
