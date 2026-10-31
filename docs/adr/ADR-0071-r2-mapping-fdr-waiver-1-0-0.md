# ADR-0071: 1.0.0 ships without a measured mapping FDR on pilot arm A2
- Status: accepted-by-khoa (2026-10-31, in session: "Waive and disclose (Recommended)")
- Date: 2026-10-31
- Spec sections: §13.2, §7.2.5, §11.6; M11-03 gap G1; ADR-0010, ADR-0019; `scripts/release/readiness.py` row R2
- Needs Khoa: yes (a release criterion cannot be evaluated)

Waives readiness R2

## Context
R2 reads "the decoy-estimated mapping FDR is ≤ 1 % at defaults", and M11-03 points it at arm A2. On the
report-freeze run at commit `fa49b8f1`, A2 and its replicate have no number to report:

| Group of 28,234 query peptides | Peptides | Handling |
|---|---|---|
| Exact target hit | 28,129 | `same_db`: exempt from both decoy lists (§7.2.5 step 5) |
| Target hit at one difference | 1 | Contaminant, rejected |
| No target hit | 104 | All contaminants |

So 0 accepted peptides are decoy-estimated. A2 searches a database built from its own genome, and a
wrong identification from that database still has an exact target hit, so decoys cannot see it.
This is structural and not a defect in the run. The replicate adds no independent evidence: its
`mapping_fdr.tsv` and `placement_costs.parquet` are byte-identical to A2's.

## Decision
Release 1.0.0 with R2 waived and state plainly that A2's mapping FDR was **not measured**. The ≤ 1 %
threshold is unchanged. The value is not reported as 0 or as a pass, and no `mapping_fdr` key is added.

## Evidence
An independent statistician review (read-only, 2026-10-31) fixed the definition from the spec before
computing it: the accepted-count-weighted FDR over the lists (decoy type × differences), with
`same_db` peptides excluded from both numerator and denominator.

- A2 and replicate: not estimable (0 accepted).
- A1, context only: 7.27e-05 (23,675 accepted). This does not validate A2, and it cannot fail: every
  list is accepted at FDR ≤ α by construction, so a value on any arm shows only that the chain ran.
- Rejected readings of A2: 0.005155 from `te_proxies.pooled_decoy_fdr` and 145 / 28,129 both count
  exempt peptides and describe a chance-hit rate, not an FDR. Counting `same_db` at q = 0 gives 0.

What the data do support is weaker. 144 of 28,129 peptides hit the reversed frames exactly, almost all
7-mers and 8-mers. Identification error is inherited from DIA-NN (input q ≤ 0.01), not estimated by
GIPPA.

Not verified: the argument that decoys cannot see these errors comes from the method and the
database's construction, and was not measured with simulated wrong identifications on A2.

## Consequences
- For `--search-db same`, GIPPA adds no estimate of identification error. Users should treat that
  error rate as coming from their search engine.
- Later releases should either define a measure that can see these errors, such as entrapment on
  a same-database arm, or re-word R2. The waiver applies to 1.0.0 only.
