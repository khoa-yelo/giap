# ADR-0005: Decisions on M1 input handling left open by the spec
- Status: accepted (items 1–3); items 4 proposed, needs Khoa
- Date: 2026-10-14
- Spec sections: §5.3, §5.4, §7.0, §9.2
- Needs Khoa: yes for item 4 only (new warning thresholds, added as new config keys)

## Context
The M1 issue specs (`docs/specs/M1-*.md`) list "Questions for the lead". The lead accepts the
spec-writer's proposed defaults (R10) and records them here.

## Decision
1. **Formats.** Genome FASTA: `U` and `-` are errors; unsafe ID characters map to `_`, and IDs
   that collide after mapping are duplicate-ID errors. Proteins: group ID = resolved accessions
   joined by `;`; `O`, digits and `-` are counted as "other" with a warning. proBed missing values
   are `.`; a BED score outside 0–1000 raises. VCF header `##reference=sha256:<hex>`; a codon
   crossing the origin raises in the VCF writer (M7 decides how to report it).
2. **Peptide and abundance tables.** Median across runs uses only runs where the peptide was
   quantified. The generic reader warns and drops rows with modification syntax (stripping is
   `giap convert`'s job) and peptides > 40 aa. Mixed `source` for one peptide: `database` wins,
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
   info under `giap preflight` and an error under `giap run` until M2 lands. Proteins of 7–40 aa
   become peptide rows.
4. **New warning thresholds (proposed; additive keys in `defaults.yaml`, Khoa to accept):**
   `preflight.decoy_share_warn: 0.5` (warn when ≥ half the peptide rows are decoys);
   warn when no row passes `input.q_value_max` (no new key);
   `preflight.min_free_ram_gb: 8` (24 GB is already implied for genome-wide Promotech, §7.0);
   `preflight.min_free_disk_gb: 10`.

## Evidence (before → after: metrics, tests, simulations)
Unit tests per issue (T-PRE-* rules); no defaults change for existing keys.

## Consequences
Item 4 lands with M1-10/M1-11; those PRs carry the `needs-khoa` label until this ADR's item 4 is
accepted.
