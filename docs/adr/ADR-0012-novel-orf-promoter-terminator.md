# ADR-0012: Promoter and terminator evidence for novel ORFs (Module 6)
- Status: accepted-by-khoa (2026-10-14: "when there is a novel orf, check for promoter and terminator
  nearby"; option "Evidence + ranking" chosen in session)
- Date: 2026-10-14
- Spec sections: §7.6.2, §7.6.3, §7.6.5, §9 `tables/novel_orfs.tsv`, §7.1.5 (terminators), decision 8
  (promoters); config `promotech_yes`, `terminatornet_yes`
- Needs Khoa: yes (adds output columns and a ranking key) — answered

## Context
§7.6.2 already gives every novel-ORF candidate start a targeted Promotech window (−200 to +20) and
`novel_orfs.tsv` reports `promoter_score` (ranking key 5). Nothing checks for a transcription terminator
after a novel ORF. Khoa asked that novel ORFs be checked for a nearby promoter AND terminator.
ORFs inside operons legitimately lack their own promoter or terminator, so neither is a calling requirement.

## Decision
1. **Terminator downstream.** For each novel ORF, look for TerminatorNet calls already made genome-wide
   (§7.1.5; TransTermHP when `--terminators transtermhp`) on the same strand whose start lies between the
   stop codon and +300 nt downstream (circular replicons wrap). Report the best confidence.
2. **Promoter upstream.** Keep §7.6.2's targeted Promotech window at the chosen start (−200 to +20);
   `promoter_supported` = `promoter_score` ≥ `promotech_yes`.
3. **Outputs** (`tables/novel_orfs.tsv`, after `promoter_score`): `promoter_supported` (bool),
   `terminator_downstream` (bool: best confidence ≥ `terminatornet_yes`), `terminator_conf` (float or NA),
   `terminator_distance_nt` (int or NA). Flags: `promoter_supported`, `terminator_supported`. The novel-ORF
   GFF3 carries the same as attributes; locus plots (V3) already draw promoters and terminators.
4. **Ranking (§7.6.5):** insert "terminator confidence, descending" after "promoter score, descending" and
   before "length, descending".
5. **Calling rules (§7.6.3) are unchanged:** promoter/terminator evidence never gates a call.

## Evidence
None needed for validity (reporting and tie-break only). Calibration of the cut-offs follows M2-13/M2-15
(RegulonDB at M6, §15.4).

## Consequences
- Module 6 specs (M6) must implement items 1–4 and test: a terminator 0–300 nt downstream on the same
  strand is found, one on the other strand or beyond 300 nt is not, origin wrap, NA when no call.
- No new config keys (reuses `promotech_yes`, `terminatornet_yes`; 300 nt is a module constant).
