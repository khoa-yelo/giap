# ADR-0005: Decisions on M1 input handling left open by the spec
- Status: accepted (items 1–3); item 4 rejected by Khoa (2026-10-14)
- Date: 2026-10-14
- Spec sections: §5.3, §5.4, §7.0, §9.2
- Needs Khoa: answered: item 4 rejected

## Context
The M1 issue specs (`docs/specs/M1-*.md`) list "Questions for the lead". The lead accepts the
spec-writer's proposed defaults (R10) and records them here.

## Decision
1. **Formats.** Genome FASTA: `U` and `-` are errors; unsafe ID characters map to `_`, and IDs
   that collide after mapping are duplicate-ID errors. Proteins: group ID = resolved accessions
   joined by `;`; `O`, digits and `-` are counted as "other" with a warning. proBed missing values
   are `.`; a BED score outside 0–1000 raises. VCF header `##reference=sha256:<hex>`; a codon
   crossing the origin raises in the VCF writer (M7 decides how to report it).
   GFF3: every replicon gets a `##sequence-region` pragma; only circular replicons also get a
   `region` feature with `Is_circular=true`; `read_gff3` skips GIPPA's own region lines.
   After the M1-05 review: CDS features must carry a phase (GFF3 3.1.26; overrides the issue's
   "phase None → `.`" edge case); the reader accepts strand `?`, rejects duplicate attribute keys
   and invalid percent-escapes, and accepts only plain digits for coordinates. VCF: INFO/FORMAT
   keys must be declared in the header; NaN → `.`, bool → Flag, empty lists and `|` inside FLAGS
   items are errors; REF/ALT are uppercase `[ACGTN]+`, ALT ≠ REF, same length, GT within the ALT
   count; header text may not contain control characters; contig IDs must be safe IDs. REF is
   checked against the genome only when the caller passes the sequences (optional argument).
   A FILTER value other than `PASS` gets a generic `##FILTER` header line; M7 replaces these with
   explicit descriptions when it defines GIPPA's filters (docs/backlog.md).
2. **Peptide and abundance tables.** Median across runs uses only runs where the peptide was
   quantified. The generic reader warns and drops rows with modification syntax (stripping is
   `gippa convert`'s job) and peptides > 40 aa. Mixed `source` for one peptide: `database` wins,
   with a warning. DIA-NN `PG.MaxLFQ` is collapsed across runs by median; intensity 0 → NA; the
   user passes `--abundance-type maxlfq`.
   After the M1-06 review: rows with any error-level issue (bad type, out of range, `nan`/`inf`
   text in q_value, pep or intensity) are dropped from the evidence frame and `peptide_runs`, and
   the issue is still reported, so a malformed q-value can never pass as evidence. Rows sharing
   (run, spectrum_id, stripped peptide, charge) are one spectrum's evidence: de-duplicated before
   intensities are summed (first row kept), with the warning `duplicate_spectrum_row`; different
   modified forms of one stripped peptide on one spectrum count once. `duplicate_spectrum_id`
   (error; rows dropped) fires only when one (run, spectrum_id) carries different stripped
   peptides. Summing modified forms (§5.3) applies across different spectra. `protein_ids` entries are stripped. Precedence when
   rows disagree: `search_db` same > other > auto; `source` database > denovo > auto. Abundance IDs resolve UniProt accessions as in
   `proteins.faa`; linear zeros are not "missing"; a header row is required.
3. **Preflight.** Every evaluated rule is recorded (`level: ok` when it passes); coordinates in
   `preflight.json` are 1-based closed; exit-code precedence 1 > 2 > 3. The chromosome-size
   threshold reuses `preflight.min_largest_contig_warn`. The abundance join check runs only when
   at least one ID matches a protein; otherwise it is reported as deferred (info). Without
   `--outdir`, preflight prints its summary and writes nothing. A missing `db_manifest.json` is
   info under `gippa preflight` and an error under `gippa run` until M2 lands. Proteins of 7–40 aa
   become peptide rows.
4. **New warning thresholds: rejected by Khoa (2026-10-14).** No new config keys. The decoy share and
   free RAM/disk are reported as Info only; the only resource warning is §7.0's ≥ 24 GB for genome-wide
   Promotech; "no row passes `input.q_value_max`" stays a warning (existing key). Original proposal:
   `preflight.decoy_share_warn: 0.5` (warn when ≥ half the peptide rows are decoys);
   warn when no row passes `input.q_value_max` (no new key);
   `preflight.min_free_ram_gb: 8` (24 GB is already implied for genome-wide Promotech, §7.0);
   `preflight.min_free_disk_gb: 10`.

## Evidence (before → after: metrics, tests, simulations)
Unit tests per issue (T-PRE-* rules); no defaults change for existing keys.

## Consequences
M1-10/M1-11 implement the Info-only behaviour above.
