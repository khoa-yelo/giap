# M2-08: `giap setup` core: database items, `db_manifest.json`, `--check`, Bakta/geNOMAD/antiSMASH/Promotech

Milestone M2 · depends on: M2-01 · estimate ≈ 360 code lines (tests not counted) · capabilities: unit
tests none (fake command runner and fetcher); live runs `needs_slurm`, `needs_net`, `needs_conda`

## Spec sections implemented

- §7.9.1: "Each item is downloaded inside its own tool environment or container where the tool has a
  downloader: `bakta_db download --type light|full`, which also runs `amrfinder_update`; `genomad
  download-database`; antiSMASH `download-databases` (skipped when the container image already
  contains them); the Promotech models zip; …". "`DIR/db_manifest.json` records name, version, URL,
  checksum, size and retrieval date for each. `--check` only verifies. An unreachable source or a
  checksum mismatch gives exit code 2, names the host, and leaves the item absent; it never falls back
  to another source. On FarmShare, `giap setup` runs as a batch job (`--bakta full` downloads about 84 GB)."
- §5.1: `giap setup --db-dir DIR [--bakta light|full|skip] [--genomad/--no-genomad]
  [--antismash/--no-antismash] [--promotech-models/--no-promotech-models] [--go] [--kegg] [--cog]
  [--contaminants] [--check]`.
- §3.2 sources: `zenodo.org` "Bakta DB v6.0 (record 14916843), geNOMAD DB v1.9 (record 14886553)";
  `ftp.ncbi.nlm.nih.gov` "AMRFinderPlus DB (pulled by `bakta_db`)"; `www.cs.mun.ca` "Promotech RF
  models (`/~lourdes/public/PromoTech_models/`)"; `dl.secondarymetabolites.org` antiSMASH DBs.
- §3.9: Bakta light "MD5 `4a6e059ded39e9c5537ef4137d2f5648`"; "delete the archive after its checksum passes".
- R6 (checksums, no mirrors, retry once later), AGENTS.md `sbatch scripts/slurm/giap_setup.sbatch`.

## Public interface

`src/giap/setupdb/`
```python
@dataclass(frozen=True, slots=True)
class DbItem:
    name: str            # "bakta_light", "bakta_full", "genomad", "antismash", "promotech_models"
    version: str; url: str; host: str
    checksum: str        # "md5:<hex>" or "sha256:<hex>" of the archive as published
    path: str            # relative to db_dir: "bakta/db-light", "bakta/db", "genomad_db", "antismash", "promotech"
    tool: str | None     # TOOLS key whose env runs the downloader; None = plain HTTPS fetch

@dataclass(frozen=True, slots=True)
class DbRecord:
    name: str; version: str; url: str; checksum: str; size_bytes: int; retrieved: str  # ISO date
    path: str; status: Literal["ok", "absent", "mismatch"]; message: str = ""

class CommandRunner(Protocol):
    def run(self, tool: str, argv: Sequence[str]) -> int: ...   # micromamba run -n giap-<tool> …
class Fetcher(Protocol):
    def fetch(self, url: str, dest: Path) -> None: ...          # raises SourceError(host, reason)

SOURCES: Mapping[str, DbItem]                                  # from src/giap/data/db_sources.yaml
def plan(opts: SetupOptions) -> list[DbItem]
def install(item: DbItem, db_dir: Path, runner: CommandRunner, fetcher: Fetcher) -> DbRecord
def verify(item: DbItem, db_dir: Path) -> DbRecord             # no network
def write_manifest(db_dir: Path, records: Sequence[DbRecord]) -> None  # merges with an existing manifest
def run_setup(opts: SetupOptions, runner: CommandRunner, fetcher: Fetcher) -> ExitCode
```
`verify` rules: Bakta → `version.json` exists with `major == 6` and `type` matching; geNOMAD →
`genomad_db/version.txt` starts `1.9`; antiSMASH → `antismash --check-prereqs` exit 0 in its env
(or image); Promotech → `RF-HOT.model` present and its SHA-256 equals the pinned value.
`install`: skip when `verify` is ok (adopts the prefetched `local/db/bakta/*`); otherwise run the
downloader into `<path>.partial`, verify, then rename; on failure delete the partial and return
`absent`/`mismatch` with the host. Downloader commands: `bakta_db download --output DIR/bakta --type
light|full`; `genomad download-database DIR`; `download-databases --database-dir DIR/antismash`;
Promotech: HTTPS fetch of the models zip, checksum, unzip, delete archive.

Files: `src/giap/setupdb/{items.py,manifest.py,run.py}`, `src/giap/data/db_sources.yaml`,
`src/giap/cli.py` (`setup` replaces the stub), `scripts/slurm/giap_setup.sbatch` (§3.8 header,
4 CPU, 16 GB, `--time=24:00:00`, runs `giap setup --db-dir $GIAP_DB_DIR --bakta full …` per AGENTS.md).

`db_manifest.json`: `{"schema_version": "1.0", "items": [DbRecord …]}` sorted by `name`.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `giap setup --db-dir D --bakta light --no-genomad --no-promotech-models` with a fake runner that creates `D/bakta/db-light/version.json` `{"major":6,"minor":0,"type":"light"}` | exit 0; manifest item `bakta_light`, `status: ok`, `version: "6.0"` |
| 2 | same, fake runner exits 1 with `zenodo.org` unreachable | exit 2; stderr `bakta_light: source zenodo.org unreachable`; no `D/bakta/db-light`, no `.partial` left |
| 3 | `--check` with nothing installed | exit 2; every planned item `absent`; no network call (fetcher/runner never invoked) |
| 4 | `--check` after example 1 | exit 0 |
| 5 | Promotech zip whose SHA-256 differs from the pin | exit 2; `status: mismatch`; archive deleted |
| 6 | prefetched `D/bakta/db` (full) present, `--bakta full` | no download; manifest `bakta_full` `ok` |
| 7 | `--bakta skip --no-genomad --no-promotech-models` | exit 0; manifest unchanged |
| 8 | two setups adding different items | manifest contains both, sorted |

## Edge cases

- Free disk below the item's expected size (`--bakta full` ≈ 84 GB unpacked) → exit 2 before
  downloading, with the measured free space.
- Running `giap setup` outside a Slurm job for `--bakta full` prints a warning (R14) but proceeds.
- antiSMASH: when `--antismash` and the pinned image already contains the DBs, record `ok` with
  `url: "image"`.
- Concurrent setups on the same `db_dir`: a lock file `db_dir/.giap_setup.lock`; the second exits 2.

## Acceptance criteria

- [ ] `tests/unit/test_setupdb.py`: examples 1–8 with fake runner/fetcher; manifest validated
  against `src/giap/schemas/db_manifest.schema.json` (new).
- [ ] `tests/integration/test_setup_live.py` (`needs_net`, `needs_conda`, `slow`): `--promotech-models` live.
- [ ] `sbatch scripts/slurm/giap_setup.sbatch` completed on FarmShare: Bakta light + full adopted,
  geNOMAD v1.9 and Promotech models downloaded; job ID, wall time, sizes in `PROGRESS.md` (§3.9).
- [ ] `scripts/probe_env.py` reads `db_manifest.json` (no new probe logic beyond reading it).

## Questions for the lead

1. The full-DB MD5 and the Promotech zip checksum are not in the spec. **Proposed default:** take the
   MD5 from the Zenodo record metadata (bakta_db verifies it itself) and the Promotech SHA-256 from
   the first download, both recorded in `db_sources.yaml` with an ADR.
2. Does `amrfinder_update` (pulled by `bakta_db`) need its own manifest row? **Proposed default:**
   yes, `amrfinderplus` with the version from `bakta/db/amrfinderplus-db/latest/version.txt`.

## Lead decisions (questions and test-author notes)
1. Q1: proposed default accepted — ADR-0014 pins the Bakta full DB MD5 (Zenodo metadata), the geNOMAD MD5
   and the Promotech RF-HOT.zip first-download SHA-256 in `db_sources.yaml`.
2. Q2: proposed default accepted — an `amrfinderplus` manifest row with the version from
   `bakta/db/amrfinderplus-db/latest/version.txt` (currently `2026-08-07.1`), status from its presence.
3. geNOMAD lives at `DB/genomad_db` (spec §7.1.1); the prefetch data was moved there and the prefetch
   script fixed. `SOURCES["genomad"].path == "genomad_db"`; no special case for `genomad/genomad_db`.
4. Accepted test-author choices: `giap.setupdb` module; CLI-named `SetupOptions`; checksum mismatch recorded
   as `status: "mismatch"` (example 5), never `ok`; adopted prefetched items keep the prefetch checksum and
   date; preflight maps bakta_* → bakta_db, genomad → genomad_db, antismash → antismash_db, counting `ok`
   rows only; free disk via `shutil.disk_usage`. The kept `promotech/RF-HOT.zip` is what the SHA-256 pin
   verifies. antiSMASH `verify` and `--go/--kegg/--cog/--contaminants` are covered by M2-09.
5. Tool commands inside `giap setup` (`bakta_db`, `amrfinder_update`, `download-antismash-databases`) run
   through the image pinned in `workflow/conf/containers.config` (Apptainer) for tools without a conda env
   (M2-01 decision 4: container-run tools have none) — the same images the pipeline uses. A missing image
   is an environment error (exit 2) naming `scripts/slurm/build_envs.sbatch`.
6. Tests never reach the network: tests that only exercise option parsing skip every resource (lead fixed
   two M1 CLI tests that downloaded the 1 GB Promotech zip once `giap setup` became real).

## Lead decisions after code review (PR #57, changes requested)
1. R6: geNOMAD is fetched by GIAP itself from the pinned Zenodo URL (record 14886553, HTTPS), MD5-checked
   against ADR-0014, and unpacked with a safe tar extractor (no absolute paths, no `..`, no links outside);
   `tool: null` — `genomad download-database` (portal.nersc.gov, unchecked) is never used. Any tool-driven
   download must point at a §3.2 host and be checksummed, or it is not used.
2. R6: prefetch provenance is adopted only when its checksum equals the pinned checksum AND its host equals
   the pinned host; otherwise the pin is used for verification and a mismatch is recorded as `mismatch`.
3. Tests can never reach the network or start tools: an autouse fixture (tests/unit setup tests and CLI
   tests) makes `ToolRunner.run`, `subprocess.run`/`Popen` and `HttpsFetcher.fetch` raise unless a test
   injects a fake explicitly.
4. TLS: redirects are followed only to https URLs on the item's official host; anything else is an error.
5. Tests kill the surviving mutants: zip-slip check, `ToolRunner` binds / env-vs-image choice / missing image
   (`ToolRunner.command` tests), https-only check, post-download verification, remove-before-rename order,
   reuse of previous `ok` provenance; the prefetch-provenance fixture uses an old date (not today).
6. `amrfinderplus` row also for `bakta_light`. The antiSMASH `url: "image"` case is M2-09's (backlog).
