# ADR-0070: 1.0.0 ships with TE4's rplL/rplJ pair failing on the pilot
- Status: accepted-by-khoa (2026-10-31, in session: "waive it")
- Date: 2026-10-31
- Spec sections: §11.6 (TE4), §7.3.4, §14.8; M10-06; `scripts/release/readiness.py` row R5
- Needs Khoa: yes (a strict true-example criterion fails)

Waives readiness R5

## Context
TE4 (§11.6) requires log2(RplL/RplJ) in [1, 3] and log2(AtpA/AtpD) in [−1, 1] by iBAQ on pilot arm
A2. It was first evaluated on the report-freeze run at commit `fa49b8f1` (every earlier pass
recorded `pending`): rplL/rplJ is −0.48 (Top3 −0.53), a fail; atpA/atpD is −0.12, a pass.

The range is unchanged. It was fixed in advance by the spec, and it stays the criterion.

The centre comes from ribosome composition (four L7/L12 per L10). The ±1 log2 width has no recorded
derivation, and 4:1 was never measured independently in this sample; the PXD078459 record gives no
protocol or growth condition.

## Decision
Release 1.0.0 with this failure waived and disclosed. It is reported as a fail, not as a pass:
in `REPORT.md`, in the release notes and in `docs/methods/abundance.md`.

**Scope.** This ADR waives only the rplL/rplJ pair of TE4 for version 1.0.0. `apply_waivers` works on
whole rows, so it also turns any other failing or pending part of R5 into `waived`. The lead must
therefore read the R5 value before step 8's PR: every part other than TE4 must say `pass` or
`deferred`. If another part does not, this waiver does not cover it, and it goes back to Khoa.

## Evidence
An independent statistician review (read-only, 2026-10-31):

- The computation is correct: iBAQ and Top3 were recomputed from the inputs and match `genes.tsv`,
  and the observable peptide counts (rplL 7, rplJ 8) reproduce with a separate digest.
- The value is identical at pass 3 (commit `8b6e93d4`), so the code changes since then did not move it.
  A1 reads +0.17, also out of range.
- rplL rests on few peptides: three carry 95 % of its signal in A2. One of them has two states on the
  same raw files, about 6 log2 apart at different retention times, which points to DIA-NN choosing
  different elution peaks. No chromatogram was inspected, so this is not proven.
- Among 48 ribosomal proteins expected at 1:1 in A2, the iBAQ SD is 0.99. If the true ratio were 4:1,
  a value inside [1, 3] would occur only about 60-68 % of the time, and one this low about 1-5 %. So
  noise covers much of the miss, but not comfortably all of it.
- No reanalysis of these peptides reached 4:1. Only a deliberately one-sided one reached +1.

Reviewer's confidence: not a TE4 bug, about 97 %; not a GIPPA quantification bug, about 85-90 %.

Not verified: chromatograms; a re-search allowing acetyl or methyl modifications; published
label-free L7/L12 to L10 values; the sample protocol. A true ratio below 4:1 in this sample is not excluded.

## Consequences
- Do not read iBAQ ratios of a single protein pair as stoichiometry on data like this.
- DIA-NN peak-choice swings reach iBAQ: 53 of 2,766 genes in A2 have a peptide swinging more than
  3 log2 between runs. The limitation is now in `docs/methods/abundance.md`. Handling it is work for
  after 1.0.
- The waiver applies to 1.0.0 only. A later release evaluates TE4 again.
