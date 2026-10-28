# M1-12: `run_preflight`, normalised inputs and `gippa preflight` CLI with exit codes

Milestone M1 · depends on: M1-02 … M1-11 · estimate ≈ 360 changed lines (≈ 160 code, ≈ 200 tests)

## Spec sections implemented

- §5.1: `gippa preflight --genome G.fasta [--proteins P.faa] [--peptides T.tsv] [--abundance A.tsv]
  [--annotation X] [--outdir DIR] [--strict]`; exit codes "0 success · 1 invalid input (preflight
  errors) · 2 environment error · 3 warnings were raised and `--strict` was given".
- §7.0: "Errors stop the run (exit 1)." · "The exit codes are 0, 1 and 3 (`--strict`)." ·
  "Normalised inputs in the work directory (§6.2)."
- §6.2 `preflight` step: "raw inputs → `preflight.json`, `genome.norm.fa`, `contig_map.tsv`,
  `replicons.tsv`, `peptides.norm.parquet`, `peptide_runs.parquet`, `proteins.long.faa`,
  `abundance.norm.parquet`". §9.1: `OUT/preflight/` holds `preflight.json`, `preflight.md`.
- §5.3: "Unsafe ID characters are mapped …; the mapping is saved." §7.1.2 `replicons.tsv` for Bakta.
- §6.7 logging: one `events.jsonl` line per check is out of scope here (runner, M2+).

## Public interface

`src/gippa/preflight/run.py`
```python
@dataclass(frozen=True, slots=True)
class PreflightInputs:
    genome: Path; proteins: Path | None = None; peptides: Path | None = None
    abundance: Path | None = None; annotation: Path | None = None

def run_preflight(inp: PreflightInputs, cfg: GippaConfig, probe: EnvProbe | None,
                  mode: Literal["preflight", "run"]) -> tuple[PreflightReport, NormalisedInputs]: ...
def exit_code(report: PreflightReport, strict: bool) -> ExitCode: ...
def write_outputs(outdir: Path, report: PreflightReport, norm: NormalisedInputs) -> None: ...
```
`NormalisedInputs` holds the in-memory genome records (safe IDs), aggregated peptides and
`peptide_runs` frames, long proteins, and the abundance frame. `write_outputs` writes
`OUT/preflight/preflight.json`, `OUT/preflight/preflight.md` always, and the normalised files
under `OUT/work/preflight/` only when the report has no error. File formats:
- `contig_map.tsv`: columns `original_id`, `contig` (FASTA order), via `write_tsv`.
- `replicons.tsv` (Bakta `--replicons`, no header): `contig`, `contig`, type
  `chromosome|plasmid|contig`, topology `c|l`, name `-`.
- `genome.norm.fa`: safe IDs, upper case (M1-02 `write_fasta`).
- `proteins.long.faa`: records > `peptide_len[1]` aa, header = group ID.
- Proteins of 7–40 aa join `peptides.norm.parquet` as rows with `protein_ids=[id]`.

`exit_code` precedence: any `error` with `exit_code=1` (input) → 1; else any environment error → 2;
else `strict` and any `warn` → 3; else 0. CLI: `preflight` in `src/gippa/cli.py` calls these, prints
`preflight.md` to stdout, and exits with the code. Without `--outdir`, nothing is written.

## Examples (tiny files written in `tmp_path`; thresholds lowered with `--set`)

| # | Command | Expected |
|---|---|---|
| 1 | `gippa preflight --genome g.fa` with `>c circular=true` + 120 kb random ACGT, `--set preflight.min_n50_warn=1 --set preflight.min_largest_contig_warn=1` | exit 0; stdout contains `0 errors` |
| 2 | same genome without the threshold overrides (largest < 1 Mb → contiguity warn) | exit 0 |
| 3 | #2 with `--strict` | exit 3 |
| 4 | genome file containing `>c\nACGU\n` | exit 1 |
| 5 | #4 with `--strict` | exit 1 (errors win over strict) |
| 6 | #1 plus `--peptides p.tsv` lacking a `peptide` column, `--outdir out` | exit 1; `out/preflight/preflight.json` exists, validates against the schema, has `peptides.schema` `error`; `out/work/preflight/` absent |
| 7 | #1 with `--outdir out` and a valid peptides file | exit 0; `out/work/preflight/` has `genome.norm.fa`, `contig_map.tsv`, `replicons.tsv`, `peptides.norm.parquet`, `peptide_runs.parquet` |
| 8 | header `>chr:1 circular=true` | `contig_map.tsv` = `original_id\tcontig\nchr:1\tchr_1\n`; `replicons.tsv` line `chr_1\tchr_1\tchromosome\tc\t-` (with `min_largest` lowered) |
| 9 | #7 run twice into two directories | `preflight.json`, `contig_map.tsv`, `replicons.tsv` byte-identical |

## Edge cases

- `--annotation` given: M1 records its role, size and SHA-256 and checks only that it is a file of a
  known kind by extension (`.json`, `.gff3`/`.gff`, `.gbff`/`.gbk`, optionally `.gz`); parsing is M2.
- Exception inside a rule → exit 5 (internal error) with the "please report" message, never 1.
- `--strict` with only `info`/`ok` checks → exit 0.
- Environment probe under `mode="preflight"` follows M1-11 question 1.

## Acceptance criteria

- [ ] `tests/unit/test_cli_preflight.py` (CliRunner): examples 1–8 assert exit codes and files.
- [ ] Example 9 is a determinism test (R7).
- [ ] `tests/unit/test_cli_exit_codes.py` no longer expects `preflight` to be a stub.
- [ ] M1 gate evidence: every `T-PRE-01 … T-PRE-25` test exists and passes; every
  `preflight.json` written in this file's tests validates against `preflight.schema.json`.

## Questions for the lead

1. Without `--outdir`, should `gippa preflight` write `./preflight/`? **Proposed default:** no; print
   the summary only (the command is "input validation only (seconds)").
2. Bakta's replicon table: original and new IDs both the safe ID, name `-`, no header. **Proposed
   default:** as written; the M2 Bakta contract test confirms the format.
3. Should proteins of 7–40 aa really become peptide rows (with `source` from config)? **Proposed
   default:** yes, per §5.3 routing; their `protein_ids` keeps the protein ID so gene pairing works.

## Lead decisions after code review
- Outputs are atomic: every run first removes a stale `work/preflight/`; files are written into a
  temporary directory under `--outdir` and moved into place with `os.replace`; `preflight/preflight.json`
  is written last as the completion marker. A run with errors leaves no normalised inputs behind.
- Resources are probed live (`shutil.disk_usage(outdir)`; Slurm/cgroup memory limit, else total RAM),
  not from `.gippa/capabilities.json`. `preflight.json` keeps only the outcome (`level`, and for
  `environment.resources` the needed vs. available verdict without raw numbers), so it is identical
  across nodes (§6.6). The raw numbers go to `logs/events.jsonl` (§6.7).
- Protein FASTA headers are written as safe IDs (`write_fasta` default); `work/preflight/protein_map.tsv`
  maps safe IDs to the original IDs and group members (§5.3 "mapping saved"). `safe_ids=False` is removed.
- `--outdir` pointing at an existing file exits 1 (invalid input), not 5.
- `db_manifest.json` (provisional): `{"<name>": {"version": "<v>", "sha256": "<hex>", ...}}`; M2's
  `gippa setup` pins it with `src/gippa/schemas/db_manifest.schema.json`.
- The probe code lives in `gippa.preflight.environment`; `route_by_length` runs once.
