# ADR-0017: Unversioned sources (KEGG REST) are recorded, not pinned

- Status: accepted (lead)
- Date: 2026-10-15
- Issue: M2-09 (PR #67)

## Context
§3.2 asks every downloaded resource to be checksum-verified. KEGG REST tables (`/link/pathway/ko`,
`/list/pathway`) are unversioned and change without notice, so no checksum can be pinned in advance.
ADR-0014 covers the pins of versioned resources only.

## Decision
Unversioned sources carry `version: ""` and `checksum: ""` in db_sources.yaml. At download, `giap setup`
validates the content (KEGG link/name format checks), records the retrieval date as version and the
sha256 of the bytes in db_manifest.json, and `giap setup --check` verifies against that recorded checksum.
Every run report states the KEGG retrieval date.

## Consequences
KEGG-derived results are reproducible from the recorded manifest, not from the package alone; two installs
on different dates may differ. All other sources stay pinned.
