# M1-11: Environment rules, `preflight.json` schema and `preflight.md`

Milestone M1 · depends on: M1-09, M1-10 · estimate ≈ 360 changed lines (schema JSON ≈ 90)

## Spec sections implemented

- §7.0 Environment rows: "tools and databases — binaries, environments and DB versions for enabled
  modules (Bakta DB v6 for Bakta 1.12) — Error (exit 2)"; "resources — free disk and RAM vs. module
  needs (genome-wide Promotech needs ≥ 24 GB) — Warn".
- §7.0 Outputs: "`preflight.json`: checks with severity, value and rule; derived routing; the evidence
  tier; replicon table." · "`preflight.md`: a human-readable summary." · Test: "The JSON report
  validates against its schema."
- §9.5: JSON Schemas live in `src/giap/schemas/`. §6.6: no timestamps or absolute paths in tables
  (applied here too, for byte-identical reports). M1 gate: "`preflight.json` validated against its schema".

## Public interface

`src/giap/preflight/environment.py`
```python
@dataclass(frozen=True, slots=True)
class Requirement:
    name: str; module: str; kind: Literal["binary", "env", "db", "model"]; min_version: str | None

@dataclass(frozen=True, slots=True)
class EnvProbe:                          # injected; real probe reads .giap/capabilities.json + db_manifest.json
    found: Mapping[str, str | None]      # requirement name → version found, None if missing
    free_disk_bytes: int; total_ram_bytes: int

def requirements(modules: Sequence[str], promoters_mode: str) -> list[Requirement]: ...
def environment_checks(reqs: Sequence[Requirement], probe: EnvProbe, need_disk_bytes: int,
                       need_ram_bytes: int) -> list[Check]: ...
```
`src/giap/preflight/report.py`
```python
@dataclass(frozen=True, slots=True)
class PreflightReport:
    giap_version: str; inputs: list[InputRef]          # role, file name (no directory), sha256, size
    checks: list[Check]; routing: Mapping[str, JsonValue]; evidence_tier: str
    replicons: list[Replicon]; contig_map: list[tuple[str, str]]
    @property
    def exit_code(self) -> int: ...     # see M1-12 for --strict

def to_json_dict(r: PreflightReport) -> dict[str, JsonValue]: ...
def write_preflight_json(path: Path, r: PreflightReport) -> None: ...    # sorted keys, indent 2, "\n"
def render_preflight_md(r: PreflightReport) -> str: ...
```
`src/giap/schemas/preflight.schema.json` (Draft 2020-12), required top-level keys:
`schema_version`, `giap_version`, `inputs`, `checks`, `summary` (`n_error`, `n_warn`, `n_info`,
`n_ok`), `evidence_tier` (enum of the four tiers), `routing`, `replicons` (`id`, `length`,
`topology`, `kind`), `contig_map`. Each check: `rule_id` (pattern `^[a-z]+\.[a-z_]+$`), `input`
(enum), `level` (enum `ok|info|warn|error`), `rule`, `value`, `message`, `exit_code` (0 or 2).
`routing`: `n_proteins_long`, `n_proteins_peptide_path`, `n_proteins_dropped`, `n_peptides`,
`abundance_used` (bool), `abundance_scale`. `additionalProperties: false` at the top level.

## Examples

| # | Input | Expected output |
|---|---|---|
| 1 | T-PRE-23 pass: modules `core`, probe finds `bakta=1.12.2`, `bakta_db=6.0`, all core tools | every `environment.tools` check `ok` |
| 2 | T-PRE-23 fail: probe has `bakta_db=None` | `environment.tools` `error`, `exit_code=2`, value `{"missing":["bakta_db"]}` |
| 3 | T-PRE-23 fail: `bakta_db=5.1` with Bakta 1.12 | `error`, value `{"wrong_version":{"bakta_db":"5.1 (need 6.x)"}}` |
| 4 | T-PRE-24 pass: `promoters: genome`, RAM 64 GB | `environment.resources` `ok` |
| 5 | T-PRE-24 fail: `promoters: genome`, RAM 16 GB | `warn`, value `{"ram_gb":16,"need_ram_gb":24}` |
| 6 | T-PRE-25: report from a passing tiny genome | `jsonschema.validate(to_json_dict(r), schema)` passes; `summary.n_error == 0` |
| 7 | T-PRE-25 negative: same dict with `level: "fatal"` | schema validation fails |
| 8 | same report written twice | byte-identical files; no absolute path or timestamp in the text |
| 9 | `render_preflight_md` with one warn and one error | Markdown with a summary line `1 error, 1 warning`, then a table `Level \| Rule \| Value \| Message` ordered error → warn → info → ok |

## Edge cases

- Optional modules off (`modules: [core]`): geNOMAD/antiSMASH requirements are not checked.
- `function` module on: GO/KEGG/COG resources are requirements (`kind="db"`).
- Probe raises or `capabilities.json` missing: one `environment.tools` error naming the probe failure.
- `value` objects must be JSON-serialisable; numpy scalars are converted before writing.

## Acceptance criteria

- [ ] `tests/unit/preflight/test_environment_rules.py`: `test_T_PRE_23_pass/_fail`,
  `test_T_PRE_24_pass/_fail` with fake `EnvProbe`s (no real tools).
- [ ] `tests/unit/test_schemas.py` gains `preflight.schema.json` in `SCHEMAS` (metaschema check).
- [ ] `tests/unit/preflight/test_report.py`: examples 6–9.

## Questions for the lead

1. `giap preflight` has no `--db-dir` (§5.1). Where does it find DB versions? **Proposed default:**
   no new option or config key; the probe reads `$GIAP_DB_DIR/db_manifest.json` (env var from
   `farmshare_env.sh`) and `.giap/capabilities.json`. Until M2 writes the DB manifest, a missing
   manifest is an environment **error only under `giap run`** and **info** under `giap preflight`.
2. RAM/disk needs per module are only given for genome-wide Promotech (24 GB). **Decision (Khoa
   rejected new thresholds, ADR-0005 item 4):** preflight warns only when `promoters: genome` and RAM
   < 24 GB; free RAM and free disk are otherwise recorded as Info values and never warn.
3. Should `inputs` store the path given on the command line? **Proposed default:** file name only
   (the full path goes to `run_manifest.json`), keeping `preflight.json` reproducible.
4. `jsonschema` is only a transitive dependency (via `frictionless`, dev extra). **Proposed default:**
   validation happens in tests only (dev extra); the runtime writer does not import `jsonschema`.
