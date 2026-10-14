# M2-02: Step registry, `giap _step`, local executor, precomputed layer formats

Milestone M2 · depends on: M0, M1-04 · estimate ≈ 330 code lines (tests not counted) · capabilities: none

## Spec sections implemented

- §6.2: "Every step is a Python function or an external process with declared file inputs and
  outputs. Python steps live in `src/giap/steps/` and are registered with a decorator that declares
  names, inputs, outputs and the module flag. Intermediate files are Parquet under the work directory."
- §5.1: "Hidden commands: `giap _step NAME --config resolved.yaml --in key=path ... --out key=path
  ...` (used by Nextflow, §6.3)"; exit codes "4 a pipeline step failed · 5 internal error".
- §6.4: "`giap run --executor local` runs the Python steps in topological order in one process.
  Every external layer must come from `annotation.precomputed` (Bakta-like JSON, operon pairs,
  promoters, Tx rates, terminators, RBS, regions), and a missing layer is an error naming the key."
- §6.6: "Each step draws from `numpy.random.default_rng(SeedSequence([config.seed,
  crc32(step_name)]))`" (already `giap.rng.for_step`).
- §5.2 `annotation.precomputed: {bakta_json, operon_pairs, promoters, tx_rates, terminators, rbs, regions}`.

## Public interface

`src/giap/steps/registry.py`
```python
@dataclass(frozen=True, slots=True)
class StepContext:
    cfg: GiapConfig; inputs: Mapping[str, Path]; outputs: Mapping[str, Path]
    rng: np.random.Generator; log: logging.Logger

StepFn = Callable[[StepContext], None]

@dataclass(frozen=True, slots=True)
class StepSpec:
    name: str; kind: Literal["py", "ext"]; inputs: tuple[str, ...]; outputs: tuple[str, ...]
    optional_inputs: tuple[str, ...]; module: Literal["core", "mobile", "bgc", "function", "opt_in"]
    precomputed: str | None          # annotation.precomputed key that replaces this step
    func: StepFn | None              # None for ext

def step(name: str, *, inputs: Sequence[str], outputs: Sequence[str], optional_inputs=(),
         module="core", precomputed: str | None = None) -> Callable[[StepFn], StepFn]: ...
def external(name: str, *, inputs, outputs, module="core", precomputed=None) -> None: ...
def get(name: str) -> StepSpec: ...
def all_steps() -> list[StepSpec]: ...                  # sorted by name
def topo_order(names: Iterable[str]) -> list[str]: ...  # output key → input key edges; ties by name
def run_step(name: str, cfg: GiapConfig, inputs: Mapping[str, Path], outputs: Mapping[str, Path]) -> None
class StepError(Exception): ...          # raised step → exit 4
class WiringError(Exception): ...        # unknown step / key / missing output → exit 5
```
`run_step`: unknown name or a key not declared → `WiringError`; a declared required input missing
or not existing → `WiringError`; creates output parent dirs; wraps any exception from `func` in
`StepError(name, cause)`; after success every declared output must exist, else `WiringError`.

`src/giap/cli.py` `_step` replaces the M0 stub: parses repeatable `--in k=path` / `--out k=path`,
loads `resolved.yaml`, calls `run_step`, maps `StepError` → exit 4 (stderr names the step and
`logs/<step>.log`), `WiringError` → exit 5.

`src/giap/runner/local.py`
```python
class LayerMissing(Exception): key: str     # "annotation.precomputed.terminators"
def plan_local(cfg: GiapConfig, targets: Sequence[str]) -> list[str]: ...  # py steps only; ext steps resolved via precomputed
def run_local(cfg: GiapConfig, workdir: Path, targets: Sequence[str]) -> None: ...
```
`src/giap/annotate/precomputed.py`
```python
PRECOMPUTED_COLUMNS: Mapping[str, tuple[str, ...]]
def load_precomputed(key: str, path: Path, replicons: Mapping[str, Replicon]) -> pd.DataFrame  # TSV → validated frame (0-based internally)
```

| Key | File | Columns (1-based closed coordinates, §9.2) | Replaces steps |
|---|---|---|---|
| `bakta_json` | Bakta JSON | Bakta 1.12.2 fields (M2-03) | `bakta` |
| `operon_pairs` | TSV | `gene_a`, `gene_b`, `uniop_p`, `source` (= §9.3 `operon_pairs.tsv`) | `uniop_input`, `uniop` |
| `promoters` | TSV | `contig`, `strand`, `end3p`, `score` (window scores; absent window = not scored) | `promoter_windows`, `promotech`, `promoter_layer` |
| `tx_rates` | TSV | `contig`, `strand`, `tss`, `tx_rate` | `promoter_calculator`, `tx_layer` |
| `terminators` | TSV | `contig`, `start`, `end`, `strand`, `confidence`, `source` | `terminatornet`, `transterm_input`, `transtermhp`, `bactermfinder`, `terminator_layer` |
| `rbs` | TSV | `locus_tag`, `tir_ostir`, `sd_dg`, `mrna_dg`, `spacing_dg`, `standby_dg`, `spacer_bp`, `start_codon` | `ostir_input`, `ostir`, `rbs_layer` |
| `regions` | TSV | `region_id`, `region_type`, `contig`, `start`, `end`, `source_tool`, `score`, `genes`, `notes` (= §9.3 `regions.tsv`) | `genomad`, `antismash`, `region_layer` |

## Examples

| # | Input | Expected output |
|---|---|---|
| 1 | `@step("toy", inputs=["a"], outputs=["b"])` on a function; `get("toy")` | `StepSpec(name="toy", kind="py", inputs=("a",), outputs=("b",), …)` |
| 2 | registering `"toy"` twice | `ValueError("step 'toy' already registered")` |
| 3 | steps x: a→b, y: b→c, z: a→d; `topo_order(["y","z","x"])` | `["x", "y", "z"]` (x before y; z tie-broken by name) |
| 4 | cycle x: a→b, y: b→a | `ValueError` naming both steps |
| 5 | `giap _step toy --config r.yaml --in a=in.txt --out b=out.txt` where func copies | exit 0, `out.txt` exists |
| 6 | same, func raises `RuntimeError("boom")` | exit 4; stderr contains `step 'toy' failed` and `boom` |
| 7 | `giap _step nope …` / `--in zz=…` / func forgets to write `b` | exit 5 each, message names the step and key |
| 8 | `run_local` with `annotation.precomputed.terminators: null` | `LayerMissing`; CLI exit 1, message `annotation.precomputed.terminators is required with --executor local` |
| 9 | `load_precomputed("terminators", f)` with a row `chr1 100 160 - 0.82 TerminatorNet` | frame row `contig=chr1, start=99, end=160, strand="-", confidence=0.82` |
| 10 | `load_precomputed("rbs", f)` missing column `sd_dg` | `ValueError("precomputed rbs: missing column sd_dg")` |

## Edge cases

- `optional_inputs` may be absent on the command line; the function sees them missing from `inputs`.
- A step whose module is off (e.g. `bgc`) is not planned; downstream optional inputs stay absent.
- `topo_order` is deterministic regardless of registration order (T: shuffle registrations).
- Precomputed coordinates are validated with `check_interval` (M1-01); origin-crossing rows keep
  `end > length` on circular replicons only.
- `--in` paths that are directories (e.g. `genomad/`) are allowed.

## Acceptance criteria

- [ ] `tests/unit/test_step_registry.py` (examples 1–4, 7), `tests/unit/test_cli_step.py` (5–7),
  `tests/unit/test_local_executor.py` (8), `tests/unit/test_precomputed.py` (9–10, one test per key).
- [ ] Every annotation step in the index contract table can be registered without conflicts
  (a test registers placeholders under a temporary registry).
- [ ] RNG in `StepContext` equals `giap.rng.for_step(name, cfg.seed)` (same first draw).

## Questions for the lead

1. Precomputed formats are not given in the spec. **Proposed default:** the table above (each equals
   the output of the layer's parse step); the simulator (M2-05) writes exactly these.
2. Exit code for a missing precomputed layer. **Proposed default:** 1 (invalid input/config), since
   no tool is missing; `giap run` validation catches it before any step runs.
3. New step names beyond §6.2 (`uniop_input`, `promoter_layer`, …). **Proposed default:** allowed
   without ADR (internal plumbing; no output changes); listed in `docs/methods/annotation_steps.md`.

## Lead decisions after code review
1. `giap run --executor local` imports the step modules through a loader (every module under
   `giap.steps`); with no steps registered it exits 5 ("not implemented"), never 0 with empty output.
2. `regions.tsv` rows with `region_type = amr` may have NA `start`/`end` (M2 index decision 6).
3. Declared outputs are written to temporary paths and renamed into place on success (atomic); outputs
   are removed before a step runs, so a stale file never passes the output check.
4. Opt-in steps follow configuration (e.g. `--add-bactermfinder`, `annotation.terminators.bactermfinder`),
   not the target list.
5. A configured precomputed layer always wins over a stale `workdir/key` on rerun.
6. Registration rejects two steps declaring the same output key.
7. A failing step writes `logs/<step>.log` with the traceback; stderr points to it.
8. A repeated `--in key=path` collects a list of paths (needed for Promotech chunks).
9. Parse errors in precomputed layers name the key, column and 1-based row.
