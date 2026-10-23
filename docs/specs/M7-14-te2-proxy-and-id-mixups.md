# M7-14: The TE2 proxy and the `id_mixups` QC row

Milestone M7 · track V (validation) · **pass 2 — half of the gate's TE2 clause** · depends on:
M7-04, M7-05 (flags `start_codon`/`stop_codon`); merged M3-22 (`validation/proxies/te_proxies.py`:
`run_chain`, `digest`, the Nissle recording), M4-16 (`id_mixups` row, NA + `requires_variants`) ·
estimate ≈ 170 code lines (+ ≈ 250 test lines) · not statistical (a zero-count assertion)

## Spec sections implemented

- §11.6 TE2: "no `high` ALT calls from A2 except at sites flagged `start_codon` or `stop_codon` |
  pilot; in-silico proxy on same-strain peptides"; §13 gate M7 "TE2 proxy passes".
- §7.5.1 row 1: "Search against this strain's Bakta proteins (`search_db = same`) | cannot occur,
  except at Sec sites and non-ATG starts, which appear as ALT calls | consistency check: ALT sites
  without the `start_codon` or `stop_codon` flag are counted as `id_mixups` in
  `detection_summary.tsv`".
- §9.3 `detection_summary.tsv` QC rows (`id_mixups`); backlog "M7: compute the M4-16 QC metric
  id_mixups … M4-16 reports NA + requires_variants until then".

## Public interface

`validation/proxies/te_proxies.py` gains
```python
def te2_result(sites: pd.DataFrame) -> ProxyResult: ...
def te2_proxy(genome: Path, bakta_json: Path) -> ProxyResult: ...   # Nissle, as TE1
```
`src/giap/variants/qc.py`: `id_mixups(sites: pd.DataFrame, support: pd.DataFrame) -> tuple[int, int]`
(count, n same-db ALT sites). The `variants` step gains the output `variants_qc`
(`work/variants_qc.parquet`); `context_stats` gains the optional input `variants_qc`.

## Rules

1. **TE2 proxy.** Peptides = M3-22's TE1 digest of every Nissle Bakta CDS/sORF protein (missed ≤ 1,
   7–30 aa), input as `database`, `search_db = same`; `run_chain` maps them on the Nissle genome;
   then `run_variants` (M7-04). Pass: no `variants.tsv` row with `confidence = high` lacking both
   `start_codon` and `stop_codon`. `ProxyResult` counts: `alt_sites`, `high`, `high_start_codon`,
   `high_stop_codon`, `violations`; up to 50 example violations (locus tag, aa_pos, ref, alt).
2. **Expected exceptions.** The 431 non-ATG N-terminal products of TE1 appear as ALT `M` at GTG/TTG
   starts (flag `start_codon`) when they map through the tolerant route (≥ 10 aa); shorter ones stay
   unmapped. The 4 TE1 misses (selenocysteine, fdnG) stay unmapped until M7-16, so no `stop_codon`
   site is expected yet; the result records 0 and says why.
3. **`id_mixups`.** Number of ALT sites (any confidence) supported by at least one `database`,
   `search_db = same` peptide and flagged neither `start_codon` nor `stop_codon`; `n` = all such
   same-db ALT sites. `context_stats` writes the row with `value` = count and `n`, replacing NA and
   the `requires_variants` flag. Without `variants_qc` (no peptides), the M4-16 behaviour stays.
4. **Wiring.** `variants_qc` is an extra output of the already-wired `variants` step; CONTEXT takes
   it as one more optional channel from VARIANTS (so CONTEXT now waits for VARIANTS). The local
   executor needs no change (optional inputs are picked up when present).
5. **Where it runs.** The proxy maps ≈ 164,000 peptides on a 5.4-Mb genome: a FarmShare job, not
   the session allocation — an entry in `scripts/slurm/weekly.sbatch` beside TE1/TE8, plus
   `tests/integration/test_te2_proxy.py` (marked as TE1's test is).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no `search_db = same` peptides | `id_mixups` 0 with `n` 0 (not NA: the check ran) |
| a same-db ALT at a GTG start with ALT `M` | not a mix-up (flag `start_codon`) |
| a same-db ALT with ALT ≠ M at a start codon | a mix-up |
| `genome_relation = related` | `id_mixups` still counted (it concerns `search_db`, not the relation) |
| TE2 violation | the proxy fails; the gate is red; the example rows name the genes |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | sim `same_strain_insilico` (M7-12) through `run_variants` | zero `high` rows without `start_codon` |
| 2 | Nissle proxy | `passed` true, `violations` 0; `high_start_codon` recorded (expected > 0: non-ATG starts with ≥ 2 distinct ≥ 10-aa N-terminal peptides that pass stratum-1 placement), not asserted |
| 3 | a hand-built same-db ALT `N → D` at an internal codon | `id_mixups` 1 |
| 4 | the same ALT supported only by de novo peptides | `id_mixups` 0 |
| 5 | `detection_summary.tsv` after a sim100k run with peptides | `id_mixups` row numeric, no `requires_variants` |

## Acceptance criteria

- [ ] `tests/unit/test_te_proxies.py` gains `te2_result` cases; `tests/unit/variants/test_qc.py`
  examples 3–4; `tests/unit/steps` example 5.
- [ ] `tests/integration/test_te2_proxy.py` (example 2) runs in the weekly Slurm suite; the job ID
  and the counts in PROGRESS.md; the result JSON under `validation/proxies/` (dates only there).
- [ ] The M4-16 `requires_variants` path is still covered for runs without peptides.
- [ ] The backlog line "M7: compute the M4-16 QC metric id_mixups" is deleted.

## Questions for the lead

1. `id_mixups` counts sites of any confidence (proposed) or only `high`? **Proposed default:** any —
   it is a QC count of contradictions, and a single contradicting same-db peptide is already a
   mix-up signal; TE2 is the `high`-only test.
2. CONTEXT now depends on VARIANTS (one more optional channel). **Proposed default:** accept; the
   alternative (VARIANTS rewriting `detection_summary.tsv`) gives one table two writers.

## Lead decisions (2026-10-23)
1. **Question 1 accepted:** `id_mixups` counts contradicting sites of any confidence. It is a QC
   count, and TE2 is the `high`-only test.
2. The weekly Slurm TE2 proxy follows R14: jobs named `giap-*`, with the CPU total checked before
   submitting. A change to the weekly scrontab, if one is needed, is pending Khoa and not blocking,
   as in MVP-04 decision 3. Until then the job is a one-off sbatch.
3. **Test-author (cbea786c), gaps 1–9: the test-author's readings are accepted.**
   - Gap 1: the lead adds `variants_qc` to the step-output assertion in test_variants_step.py.
   - Gap 2: the tests live in new files.
   - Gap 3: `support` is a long `(variant_id, peptide, source, search_db)` frame.
   - Gap 4: the TE2 JSON mentions M7-16 when there are 0 stop_codon sites.
   - Gap 5: no ALT sites counts as a pass.
   - Gap 6: TE2 runs on weekly.sbatch's Nissle task (not protected). The first run is a one-off
     sbatch; the scrontab is unchanged (decision 2).
   - Gap 7: example 5 runs as the sim100k integration test on Slurm.
   - Gap 8: example 1 activates once #225 merges.
4. **Implementer (c47da49e).** The Nissle TE2 proxy passed (Slurm 1871578): 233 ALT sites, 80
   `high`, 80 `high_start_codon`, 0 violations. Lead edits to existing tests:
   - `variants_qc` (`work/variants_qc.parquet`) was added to OUTPUTS in test_variants_step.py and to
     OUT_SPEC in test_nf_variants_registry.py. It is the only Parquet output. The oracle
     comparisons cover the published tables only.
   - context.nf.test passes `input[11]`.
   - The `sim100k_real` M4 golden was regenerated; only the `id_mixups` row changes, now
     value 0 / n 0.
   - The GTG TE2 test runs at 400 kb: placing N-terminal peptides needs at least 100 in the
     tolerant stratum, and 100 kb gave 44.
   - Accepted: GIAP_STEP publishDir and the local `publish()` skip `work/...` outputs. Work
     intermediates are never published, and no other step had a `work/` output.
   - Accepted: the regenerated step_registry.json mirror.
