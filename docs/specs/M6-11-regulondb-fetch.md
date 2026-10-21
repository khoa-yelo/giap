# M6-11: **BLOCKED** — RegulonDB v12 fetch and the K-12 pair and terminator truth

Milestone M6 · track V (validation) · **blocked, not scheduled** · depends on: network access to
`regulondb.ccg.unam.mx` (see below) · estimate ≈ 180 code lines in `validation/` (+ ≈ 150 test
lines) · not statistical (M6-12 is)

## Status: blocked

| Date | Probe | Result |
|---|---|---|
| 2026-10-14 | PROGRESS.md "Blocked" entry | `regulondb.ccg.unam.mx` unreachable from FarmShare |
| 2026-10-21 (lead) | HTTPS GET | connection refused |
| 2026-10-21 (this issue) | TLS handshake + `GET /` with the official GlobalSign intermediate added to the CA bundle | **HTTP 200**; the server certificate was issued 2026-10-14 by "GlobalSign GCC R46 OV TLS CA 2025" and the deployment **does not send the intermediate**, so a normal client fails with `unable to get local issuer certificate`. The matching root (`GlobalSign_Root_R46.pem`) is already in FarmShare's trust store. |

So the likely cause is a server-side incomplete certificate chain introduced on 2026-10-14, not a
block or an outage. **This needs Khoa, not the lead:** may GIAP fetch the intermediate from
GlobalSign's own CA endpoint (`http://secure.globalsign.com/cacert/gsgccr46ovtlsca2025.crt`) and
pin it for this download, or must we wait for RegulonDB to publish a complete chain? R6 is not at
issue (the data would come from RegulonDB itself; no mirror, no unofficial copy), but trusting a
hand-fetched intermediate is a configuration decision.

Until Khoa answers, this issue and M6-12 stay in PROGRESS.md › Blocked, and M6 proceeds without
them (index, "RegulonDB"). **No substitute source may be used.**

## Spec sections implemented (when unblocked)

- §15.3 table row: "RegulonDB v12 | transcription units and terminators | §15.4. Not committed;
  downloaded to `$GIAP_DATA_DIR` by `validation/fetch_regulondb.py`".
- §15.4 steps 1–2: "Inputs. The K-12 Bakta JSON (recorded), UniOP and TerminatorNet run on the K-12
  genome"; "Pair truth. A same-strand adjacent pair is positive when both genes lie in one
  RegulonDB transcription unit. It is negative when they lie in different TUs and no TU contains
  both. Everything else is ambiguous and excluded."
- §3.9 download discipline (checksums, provenance, no commit of the data), ADR-0014 and ADR-0017
  (unversioned sources).

## Scope

1. `validation/fetch_regulondb.py`: download the v12 TU and terminator datasets to
   `$GIAP_DATA_DIR/regulondb/<version>/`, record URL, date, bytes and SHA-256 in a provenance JSON
   beside them, verify before use, never commit the data. Same shape as
   `validation/fetch_uniprot_k12.py` (merged).
2. `validation/regulondb_truth.py`: parse TUs and terminators; map RegulonDB gene identifiers to
   the recorded K-12 Bakta locus tags through coordinates (± a tolerance) and gene names, reporting
   the mapped share; build
   - `pair_truth.tsv`: `gene_a`, `gene_b`, `label` enum{positive, negative, ambiguous}, `tu_ids`;
   - `terminator_truth.tsv`: `contig`, `start`, `end`, `strand`, `source_id`.
3. The ambiguity rule of §15.4 step 2 is implemented literally and the excluded share is reported.
4. A Slurm job is **not** needed (the download is small); the mapping runs in seconds.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| download unreachable | clear error naming the host and this issue; the PROGRESS.md Blocked entry is the documented state |
| a RegulonDB gene with no Bakta counterpart | counted and reported; the pair is `ambiguous` |
| a TU of one gene | contributes no positive pair |
| two TUs sharing a pair (overlapping TUs) | positive (the spec's "both genes lie in one TU") |
| genes adjacent in Bakta but separated by a RegulonDB-only gene | `ambiguous` |
| version string absent from the download | ADR-0017 applies: record the date and the file hash |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | TU `lacZYA` = (lacZ, lacY, lacA) | pairs (lacZ, lacY) and (lacY, lacA) positive |
| 2 | lacA and the next gene downstream in a different TU, with no TU containing both | negative |
| 3 | a pair whose downstream gene is in no TU at all | ambiguous, excluded |
| 4 | a Bakta locus tag with no RegulonDB gene | every pair touching it ambiguous |
| 5 | the whole K-12 genome | ≥ 2,000 labelled same-strand adjacent pairs, and the mapped share recorded (a share below 0.8 is escalated, not worked around) |
| 6 | fetch run twice | identical hashes; the provenance JSON records both dates |

## Acceptance criteria

- [ ] `validation/tests/test_regulondb_truth.py`: examples 1–4 and the degenerate rows on a small
  hand-written RegulonDB-shaped fixture **committed as a fixture we wrote**, not as RegulonDB data.
- [ ] The fetcher is offline-testable (the parser takes a path; the network call is one function).
- [ ] Provenance JSON validated; no RegulonDB bytes under version control (a CI check).
- [ ] Example 5 recorded in PROGRESS.md with the counts — this is the evidence M6-12 consumes.

## Questions for Khoa (not the lead)

1. The TLS question above. Nothing else unblocks this issue.
2. If the chain is not fixed and the intermediate may not be pinned, M6's gate stays partial
   indefinitely. Khoa should then decide whether M6 is declared done on its simulation half with
   §15.4 deferred to M10 (the pilot already needs K-12 data), or whether the milestone waits.
