# ADR-0014: Checksum sources for the Bakta full DB and the Promotech model (M2-08)
- Status: accepted
- Date: 2026-10-14
- Spec sections: §3.2, §3.9, R6; docs/specs/M2-08-gippa-setup-core.md question 1
- Needs Khoa: no (fills a gap: the spec gives no checksum for these two artefacts; no default, threshold or
  interface changes)

## Context
R6 requires every downloaded artefact to be checksummed. §3.9 gives the Bakta light DB MD5 but not the full
DB's, and Promotech's www.cs.mun.ca directory publishes no checksum.

## Decision
`src/gippa/data/db_sources.yaml` pins:
- Bakta full DB v6.0 `db.tar.xz`: MD5 `4c1115e40abfa2b464ae5dd988bdd88e` from the Zenodo record 14916843
  metadata (the record lists MD5 for every file; `bakta_db` checks it too).
- geNOMAD DB v1.9 `genomad_db_v1.9.tar.gz`: MD5 `67244b528bb8bed464d1ca147136d33e` (Zenodo 14886553).
- Promotech `RF-HOT.zip`: SHA-256 `fcd96ffe5a43b579af965839ee2bb896de94de3601978d88bcabd32c4fb5093f`, recorded on
  the first download (2026-10-14, job 1846845, zip integrity verified). Trust-on-first-use: a later
  download with a different hash is a mismatch, not an update; changing the pin needs a new ADR.

## Evidence
`$GIPPA_DB_DIR/prefetch_manifest.tsv` (prefetch job 1846845): all three verified against these values.

## Consequences
`gippa setup` verifies against these pins and adopts the prefetched copies without re-downloading.
