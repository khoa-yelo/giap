# ADR-0013: `cterm_extension` locus class (M3-14, Module 6)
- Status: accepted-by-khoa (2026-10-14, in session: "New class cterm_extension")
- Date: 2026-10-14
- Spec sections: §7.2.6 (locus classes), §7.6.1 (novel candidates), §7.6.3 (classes and evidence), §9
- Needs Khoa: yes (adds a class value) — answered

## Context
A peptide on the same strand and in frame with an annotated CDS (or sORF) can run across its 3′ end, past
the annotated stop codon — stop-codon readthrough or a mis-annotated stop. §7.6.3 has no class for it: it
is not `cds_outframe` (that is "other frame") and not `intergenic` ("no annotated CDS"). The M3-14
implementation called it `cds_outframe` (code review, PR #50).

## Decision
1. New locus class `cterm_extension`: same strand, in frame with an annotated CDS/sORF, overlapping it and
   extending past its annotated stop (3′ end), including across the origin on circular replicons.
2. The placement is assigned to that gene (like `nterm_extension`), and the locus is a Module 6 candidate
   (§7.6.1 list extended by `cterm_extension`). Evidence before calling (§7.6.3 row): the in-frame
   continuation to the next stop and a peptide spanning the annotated stop; labelled
   `readthrough_candidate`. The annotated stop's codon is reported.
3. Downstream ORF reconstruction (§7.6.2) applies unchanged; `annotation_patch.gff3` writes it like an
   N-terminal extension (existing locus tag, new end, `Note=cterm_extension`).

## Consequences
M3-14 classifies it; Module 6 specs implement items 2–3. Output vocabularies (§9 locus classes, novel ORF
classes) gain one value. No thresholds or defaults change.
