# M2-09: `giap setup --go --kegg --cog --contaminants`

Milestone M2 · depends on: M2-08 · estimate ≈ 220 code lines (tests not counted) · capabilities: unit
none (fake fetcher); live `needs_net(purl.obolibrary.org)`, `needs_net(rest.kegg.jp)`,
`needs_net(ftp.ncbi.nlm.nih.gov)`, `needs_net(github.com)`

## Spec sections implemented

- §7.9.1 items: "`go-basic.obo` and `goslim_generic.obo`; KEGG KO → pathway links and pathway names;
  COG 2020 tables; the contaminant library (at a pinned commit)." Manifest and exit-2 rules as M2-08.
- §3.2 sources: `purl.obolibrary.org`, `*.geneontology.org` "`go-basic.obo`"; `rest.kegg.jp` "KO →
  pathway links (academic use)"; `ftp.ncbi.nlm.nih.gov` "COG category definitions".
- §7.2.8: "the universal contaminant library (Frankenfield 2022, 381 entries), downloaded by `giap
  setup --contaminants` (vendor it only if its licence allows)".
- R5/NOTICE: "KEGG academic use".

## Public interface

New `DbItem`s in `src/giap/data/db_sources.yaml`, plain HTTPS fetch (`tool: None`):

| Name | URL (pinned) | Files under `db_dir` | Verify |
|---|---|---|---|
| `go_basic` | `https://purl.obolibrary.org/obo/go/releases/<YYYY-MM-DD>/go-basic.obo` | `go/go-basic.obo` | header `data-version: releases/<date>`; sha256 |
| `goslim_generic` | `https://purl.obolibrary.org/obo/go/releases/<YYYY-MM-DD>/subsets/goslim_generic.obo` | `go/goslim_generic.obo` | same |
| `kegg_ko_pathway` | `https://rest.kegg.jp/link/pathway/ko` | `kegg/ko_pathway.tsv` | ≥ 1 line `ko:K\d{5}\tpath:map\d{5}`; sha256 recorded at retrieval |
| `kegg_pathway_names` | `https://rest.kegg.jp/list/pathway` | `kegg/pathway_names.tsv` | lines `map\d{5}\t…` |
| `cog2020_def` | `https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/cog-20.def.tab` | `cog/cog-20.def.tab` | sha256 pin |
| `cog2020_fun` | `https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/fun-20.tab` | `cog/fun-20.tab` | sha256 pin |
| `contaminants` | GitHub raw URL at a pinned commit of the Frankenfield 2022 library | `contaminants/contaminants.fasta` | 381 records; sha256 pin |

```python
def verify_go(path: Path, release: str) -> bool
def verify_kegg_links(path: Path) -> bool
def count_fasta_records(path: Path) -> int
```

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `giap setup --db-dir D --bakta skip --no-genomad --no-promotech-models --go` with a fake fetcher serving an OBO whose header is `data-version: releases/2026-07-01` and pin `2026-07-01` | exit 0; items `go_basic`, `goslim_generic` `ok` |
| 2 | same, header `releases/2026-09-01` | exit 2; `go_basic: version mismatch (expected releases/2026-07-01)` |
| 3 | `--kegg`, fetcher returns `ko:K00001\tpath:map00010\n` | `ok`; checksum = sha256 of the bytes; `version` = retrieval date |
| 4 | `--contaminants`, FASTA with 380 records | exit 2, `mismatch`, `expected 381 records` |
| 5 | `--cog`, `ftp.ncbi.nlm.nih.gov` times out | exit 2; message names the host; nothing left under `D/cog/` |
| 6 | `--check --go --kegg --cog --contaminants` after a full install | exit 0, no network |

## Edge cases

- KEGG REST data are unversioned: the checksum is recorded, not pinned; `--check` verifies the file
  against the recorded checksum, not against the server.
- The contaminant library is downloaded, never committed (licence unknown until checked; R11).
- Two items share `D/go/`; partial-download cleanup removes only the failing item's file.

## Acceptance criteria

- [ ] `tests/unit/test_setup_resources.py`: examples 1–6 with a fake fetcher.
- [ ] `tests/integration/test_setup_resources_live.py` (`needs_net(...)` per host, `slow`): each item
  downloads and verifies.
- [ ] The live run is part of `giap_setup.sbatch` (M2-08); sizes in `PROGRESS.md`.
- [ ] `NOTICE.md`: KEGG academic use; contaminant library licence as found.

## Questions for the lead

1. Which GO release to pin? **Proposed default:** the newest monthly release on the day this PR is
   written, recorded in `db_sources.yaml` (re-pinned only via ADR).
2. Contaminant library source repository and commit are not named. **Proposed default:** the
   repository cited by Frankenfield et al. 2022 (J. Proteome Res.), HEAD on the PR date, pinned by
   commit SHA; if its licence allows, still not vendored in v1.

## Lead decisions (after the test-author's pass)
1. Cross-host redirects: allowed only to hosts on a per-source allowlist of the same organisation's
   official hosts (GO: purl.obolibrary.org → release.geneontology.org, current.geneontology.org); any other
   host is a SourceError (R6). Same-host redirects and the redirect cap stay as in M2-08.
2. A body shorter than Content-Length (truncated transfer) raises SourceError and leaves no file behind.
3. The test-author's choices are accepted: GO version = release date; both `data-version` header forms are
   valid; KEGG rows carry checksum "" and `--check` compares with the checksum recorded in the manifest at
   download (ADR-0014); antiSMASH `url: "image"` is verified by `antismash --check-prereqs` inside the image.

## Lead decisions after code review (PR #67, changes requested)
1. `--check` accepts antiSMASH installed from the image (`url: "image"`): a previous `ok` record with url
   `image` passes after `antismash --check-prereqs` succeeds in the pinned image. Test.
2. Live test tests/integration/test_setup_resources_live.py (`needs_net` per host, `slow`; runs only in the
   weekly Slurm job) downloads each item into a tmp dir and verifies it. Sizes go to PROGRESS.md (lead):
   go-basic.obo 32.2 MB, KEGG links 2.3 MB, other files < 0.4 MB.
3. KEGG checksum rule is ADR-0017 (not ADR-0014); decision 3 above is amended accordingly.
4. A body without Content-Length stops at the item's `size_gib` (or 1 GiB when unset) → SourceError. Test.
5. `_load_sources` runs every `path` through `_safe_name` (no `..`, no absolute paths). Test.
6. db_manifest.json is written with mode 0644 (shared DB directory). Test.
