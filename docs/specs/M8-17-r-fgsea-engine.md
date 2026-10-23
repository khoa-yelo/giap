# M8-17: the R fgsea engine for `enrichment_gsea` (ADR-0038)

Milestone M8 · track T (tool integration) · **pass 2** · depends on: M8-05 (`gsea_table`), M8-06a
(#181, the `function` step, merged), M8-06b (#187, `function.nf`), M8-12 (#192, the golden inputs
and outputs, the R pins) · estimate ≈ 470 code lines (+ ≈ 450 test lines), so it is split into
**M8-17a** (the R side and its wrapper, ≈ 290) and **M8-17b** (the engine switch in `gsea_table`,
the step, Nextflow and preflight wiring, ≈ 180) · statistical (label `stats`) for 14b only

> **Numbering.** `docs/specs/M8-17-sim-function-determinants.md` already holds the number M8-17
> (giapsim scenarios). The lead named this file; see lead decision 1.

## Spec sections implemented

- **ADR-0038** (accepted-by-khoa, 2026-10-22), which overrides §2.4 G15 and the §7.7.4 "Method"
  bullet. `enrichment_gsea.tsv` comes from `fgsea::fgseaMultilevel` 1.36.2 in R 4.5, run as a
  separate process in its own environment.
- §7.7.4: TU ranking unit, sizes in TUs within `function.gsea_size`, ES, NES, p, `log2err`, BH q
  within each collection, the leading edge in genes, and the GO `redundant` re-test. All of these
  stay in Python and are unchanged. Only the engine call changes.
- §7.7.5 and §11.5: the fgsea golden data on #192 become a **contract test of the wrapper**
  (ADR-0038 item 7).
- R5: R (GPL-2) and fgsea (MIT) run only as a separate process. Nothing under `src/` mentions
  `Rscript`, so `scripts/ci/check_no_rscript.py` still holds. No `rpy2`.
- R6 and §3.2: packages come only from `conda.anaconda.org` (conda-forge and bioconda). The image
  comes only from `quay.io/biocontainers`, pulled by digest. **Nothing is downloaded from
  bioconductor.org**, because the environment has no Bioconductor data package.
- R7 and §6.6: one seed integer from `giap.rng.for_step("function.gsea", seed)`, and `nproc = 1`.
  The same inputs and seed give a byte-identical `tables/enrichment_gsea.tsv`.
- §6.4: both executors call the same step function. The R call happens inside the `function` step,
  so local and Nextflow runs stay byte-identical (lead decision 2).
- §11.4: recorded R outputs under `tests/data/recorded/fgsea/1.36.2/`. CI never runs R.
- ADR-0038 item 5 (accepted): **the detection ORA stays on `scipy.stats.nchypergeom_wallenius`.
  The pipeline never calls goseq**, and goseq is not in the engine environment. goseq runs only
  in `giap-golden-r` to produce the §11.5 golden data.

Not in scope: removing `giap.stats.fgsea` and its kernels. That belongs to the hardening pass
(ADR-0038 item 6). PR #193 stays a draft. M8-13's fgsea half is superseded by the contract test
here; its Wallenius half is unchanged.

## Design in one paragraph

`gsea_table` stays pure. It builds the TU ranking, the in-range sets per collection and, after
BH, the GO re-test sets, exactly as it does today. It then hands them to an **injected engine**,
a callable that takes one `EngineCall` (a ranking plus a list of jobs) and returns rows per job.
The step injects `RFgseaEngine` from `giap.io.fgsea_r`. That engine writes three TSVs to a
temporary directory, runs `workflow/bin/giap_fgsea.R` as a subprocess in the `giap-fgsea` env or
image, and parses one result TSV. Unit tests inject fakes or a replay of recorded R outputs. A
step makes **at most two R calls**: one for every collection's main run, and one for every GO
re-test, batched across namespaces, made only when at least one re-test exists.

---

## M8-17a: the environment, `giap_fgsea.R` and the Python wrapper

### 1. Environment `workflow/envs/fgsea.yml` and its lock

```yaml
# R fgsea engine for the function step's GSEA (ADR-0038; GIAP_SPEC.md §7.7.4, R5).
# fgsea is MIT, R is GPL-2: run only as a separate process. Exact solve: locks/fgsea.explicit.txt.
name: giap-fgsea
channels: [conda-forge, bioconda]
dependencies:
  - r-base=4.5.3
  - bioconductor-fgsea=1.36.2
  - r-data.table=1.18.6.1
```

- These pins match #192's `golden_r.explicit.txt`: `r-base-4.5.3-h502d0c9_4`,
  `bioconductor-fgsea-1.36.2-r45ha27e39d_0` and `r-data.table-1.18.6.1-r45hbd9402d_0`. The same
  build of fgsea is what makes the contract test (T4) meaningful.
- The file passes `giap.annotate.tools.env_problems`: the channels are exactly conda-forge and
  bioconda, every dependency is pinned exactly, and the name is `giap-<stem>`.
- `workflow/envs/locks/fgsea.explicit.txt` is written by `sbatch scripts/slurm/build_envs.sbatch
  fgsea` (`micromamba env export --explicit`), the same way the ostir and uniop locks are. There
  is no pip lock.
- **No data packages.** The lock must contain none of `bioconductor-goseq`,
  `bioconductor-genelendatabase`, `bioconductor-go.db`, `bioconductor-data-packages`, any `org.*`,
  `TxDb.*` or `BSgenome.*`, and none of `r-biasedurn`. Those are the packages whose post-link
  scripts download from bioconductor.org (ADR-0038 item 5), and goseq is never part of this env.
  A unit test greps the committed lock for these names (A3).
- The build log of the env job is kept beside the lock's PROVENANCE line (job ID, host, date).

### 2. BioContainers image (apptainer profile)

- The tag is `quay.io/biocontainers/bioconductor-fgsea:1.36.2--r45ha27e39d_0`, the bioconda build
  in the lock. The implementer resolves it to a digest with `python scripts/resolve_digest.py`,
  which reads quay.io, an allowed §3.2 source.
- The pin lives in one place: `giap.annotate.tools.FGSEA_IMAGE =
  "docker://quay.io/biocontainers/bioconductor-fgsea@sha256:<digest>"`, plus a dated
  "resolved from" comment in `workflow/conf/containers.config`, in the same style as the
  existing tags.
- `scripts/slurm/build_envs.sbatch fgsea` also pulls the image to
  `$GIAP_ROOT/containers/<mangled ref>.img`. That is the path `missing_tools` already derives
  with `_image_file_name`.
- If the tag does not exist on quay.io, the implementer records it under PROGRESS.md › Blocked
  (R6) and fgsea stays conda-only under apptainer, as uniop and ostir are today. The image is
  never built or mirrored by us.

### 3. `workflow/bin/giap_fgsea.R`

The script lives outside `src/` (R5). Its shebang is `#!/usr/bin/env -S Rscript --vanilla`. Nextflow
puts `workflow/bin` on the task `PATH`; the wrapper calls the script by absolute path.

**CLI** (every argument is required; there is no default inside R):

```
giap_fgsea.R --stats stats.tsv --sets sets.tsv --jobs jobs.tsv \
             --sample-size 101 --eps 0 --n-perm-simple 1000 --seed <int> \
             --out results.tsv --versions versions.txt
```

**Inputs.** All inputs are UTF-8, tab-separated and uncompressed, with a header and `\n` line ends.

| File | Columns | Rules |
|---|---|---|
| `stats.tsv` | `tu_id`, `score` | one ranking per call. Rows are in the **seeded tie order** (§5). `score` is written with `repr(float)`, which round-trips in R's `read.delim` as 17 significant digits. Every score is finite and every `tu_id` is unique |
| `jobs.tsv` | `job_id`, `min_size`, `max_size` | jobs are run in file order. `job_id` is `main:<collection>` or `retest:<collection>:<term>` |
| `sets.tsv` | `job_id`, `set_id`, `member` | long form. Sets keep file order within a job, and every member is a `tu_id` in `stats.tsv` |

**Per job:**

```r
set.seed(seed)
fg <- fgseaMultilevel(pathways, stats, sampleSize = sample_size, minSize = min_size,
                      maxSize = max_size, eps = eps, nPermSimple = n_perm_simple,
                      scoreType = "std", nproc = 1)
```

- `stats` is `setNames(score, tu_id)` in file order. `pathways` is
  `split(member, factor(set_id, levels = unique(set_id)))`. This is exactly #192's construction,
  so a one-job call on the golden inputs is `make_golden_r.R`'s call.
- `set.seed(seed)` runs **before every job**. A job's numbers therefore do not depend on which
  other jobs share the call. With one job this is the golden script's single `set.seed`.
- Rows are re-ordered to the input set order:
  `fg[match(intersect(names(pathways), fg$pathway), fg$pathway), ]`.

**Output `results.tsv`.** It is written with `writeLines` and `useBytes = TRUE`, the same way as
`make_golden_r.R`:

| Column | Format |
|---|---|
| `job_id` | as given |
| `pathway` | the set ID |
| `size` | integer |
| `ES`, `NES`, `pval`, `log2err` | `format(x, digits = 17)`, or `NA` |
| `leadingEdge` | TU IDs joined by `;` in fgsea's order. The cell is empty when the edge is empty |

- A job whose sets all fall outside `[min_size, max_size]` adds no rows. It is not an error.
- With no job at all, the script writes the header only. The wrapper never makes such a call
  (§4), but the script stays total.

**Provenance `versions.txt`.** One `key<TAB>value` per line:

- `R` = `R.version.string`;
- `fgsea` = `packageVersion("fgsea")`;
- `data.table`, `BiocParallel` = their `packageVersion`;
- `seed`;
- `jobs` = the number of jobs.

Then one `warning<TAB><job_id><TAB><message>` line per distinct warning. Warnings are caught with
`withCallingHandlers` and `muffleWarning`, as in #192. This catches the ties warning, the
"unbalanced … not calculated properly" warning and the "likely overestimated" warning.

**Validation.** The script writes a one-line message to stderr and exits with a code from this
table. It never writes a partial `results.tsv`: it writes to `results.tsv.tmp` and renames.

| Exit code | Meaning |
|---|---|
| 0 | ok |
| 2 | usage error (missing or unknown argument, non-integer seed) |
| 3 | invalid input: non-finite or NA score, duplicate `tu_id`, a member not in `stats`, a job in `sets.tsv` missing from `jobs.tsv`, `min_size < 1` or `max_size < min_size` |
| 1 | any R error (R's default), including errors from `library(fgsea)` |

### 4. Python wrapper `src/giap/io/fgsea_r.py`

Process I/O is allowed only in `io/` and `steps/` (AGENTS.md). This module is that I/O. It never
spells `Rscript`.

```python
FGSEA_VERSION = "1.36.2"
SCRIPT = "giap_fgsea.R"               # resolved as workflow_main().parent / "bin" / SCRIPT
TIMEOUT_S = 3600                      # per call (lead decision 6)

class FgseaEngineError(StepError): ...      # exit code 4, like ToolUnavailable / bad resources
class FgseaEngineUnavailable(FgseaEngineError): ...

@dataclass(frozen=True, slots=True)
class Launcher:                       # how to start the script
    kind: Literal["prefix", "image"]
    path: Path                        # env prefix, or the .img file

def resolve_launcher(environ: Mapping[str, str]) -> Launcher: ...
def r_fgsea_engine(launcher: Launcher, *, timeout_s: float = TIMEOUT_S,
                   tmp_root: Path | None = None) -> GseaEngine: ...
def write_call(call: EngineCall, d: Path) -> None: ...            # stats/jobs/sets TSVs
def read_results(path: Path, call: EngineCall) -> EngineResult: ...
def read_versions(path: Path) -> EngineProvenance: ...
```

**Launcher resolution.** The first match wins:

1. `GIAP_FGSEA_PREFIX`: a micromamba env prefix. The command is `[script, *args]`, with `PATH`
   set to `<prefix>/bin` followed by the inherited `PATH`.
2. `GIAP_FGSEA_IMAGE`: an `.img` file. The command is `["apptainer", "exec", "--cleanenv",
   "--bind", tmp, "--bind", script_dir, img, script, *args]`.
3. `$GIAP_ROOT/micromamba/envs/giap-fgsea`, if its `.giap-built` matches the sha256 of
   `workflow/envs/fgsea.yml`. This is the same test as `preflight.environment._env_built`, and
   the function is reused, not copied.
4. Otherwise `FgseaEngineUnavailable("the R fgsea environment giap-fgsea is not built: sbatch
   scripts/slurm/build_envs.sbatch fgsea")`.

The child environment always sets `R_LIBS_USER`, `R_LIBS` and `R_PROFILE_USER` to empty or
`/dev/null` and `LC_ALL=C`, so a user's R library can never shadow the pinned fgsea.

**One call:**

1. Make a fresh `tempfile.mkdtemp(dir=tmp_root)`. In Nextflow `tmp_root` is the task directory;
   locally it is the system temp directory.
2. `write_call`.
3. `subprocess.run(cmd, cwd=tmp, env=..., capture_output=True, timeout=timeout_s, check=False)`.
4. Check the outcome:
   - A non-zero exit raises `FgseaEngineError`. The message carries the exit code, the job count
     and the last 20 lines of stderr.
   - `subprocess.TimeoutExpired` kills the process group and raises `FgseaEngineError("R fgsea
     timed out after <s> s")`.
   - `FileNotFoundError` or `PermissionError` on start raises `FgseaEngineUnavailable`.
5. `read_versions`. Raise `FgseaEngineError` unless `fgsea == FGSEA_VERSION`.
6. `read_results`, which validates the output:
   - the header is exact;
   - every `job_id` is one that was sent, and a `main:` job's rows are a subset of its in-range
     sets in the sent order;
   - `size` equals the size Python sent. A mismatch means the names did not match, so the call
     fails;
   - every leading-edge TU is a member of its set.

   Any violation raises `FgseaEngineError`.
7. Remove the temporary directory, also when the call fails.

The engine never retries. A failed R call fails the step (lead decision 5).

**Provenance.** The step logs one §6.7 event per call, `function.gsea_engine`, with `{"call":
"main" | "retest", "jobs": n, "r": ..., "fgsea": "1.36.2", "launcher": "prefix" | "image",
"warnings": [...]}`. The events log is not under `tables/`, so R7 is unaffected. No timing goes
into any table.

### 5. The seed handoff (spec §6.6)

`gsea_table` keeps its single generator `rng = for_step("function.gsea", cfg.seed)`, which is
already a sub-stream that the step passes in. The order of draws is:

1. `perm = rng.permutation(n)`. This is the seeded tie order and is unchanged (M8-02 lead
   decision 7).
2. `r_seed = int(rng.integers(1, 2**31 - 1))`. This is one integer in `[1, 2147483646]`, inside
   R's integer range and never `NA_integer_`.

Both R calls of the step use `r_seed` (ADR-0038 item 3: one integer). `stats.tsv` lists the TUs in
`perm` order. fgsea sorts the stats decreasing with R's `order(method = "radix")`, which is stable,
so tied scores keep the seeded order and never fall back to TU-ID order. T5 checks this on the real
engine. The fakes do not consume `rng`, so the draws of anything after the engine call do not
depend on the engine.

### 6. Tests (14a)

All of the following run in CI without R unless marked otherwise.

- **T1 env and lock (unit).**
  - `fgsea.yml` passes `env_problems`.
  - The pins equal `FGSEA_VERSION`, r-base 4.5.3 and data.table 1.18.6.1.
  - The lock exists and contains the same three builds as #192's lock.
  - The lock contains none of the forbidden names in §1.
- **T2 wrapper I/O (unit).**
  - `write_call` produces byte-exact TSVs for a hand-built `EngineCall`, including `repr`
    round-trip of the scores.
  - `read_results` parses each recorded `results.tsv` and rejects each mutation: a wrong header,
    an unknown job, a size mismatch, a leading-edge TU outside its set, or a non-numeric ES.
  - `read_versions` rejects fgsea 1.38.0.
- **T3 subprocess behaviour (unit, no R).** A fake script stands in for `giap_fgsea.R`: a small
  `sh` script in `tmp_path` reached through `GIAP_FGSEA_PREFIX`. It covers:
  - exit 1 with stderr, which becomes `FgseaEngineError` with the stderr tail;
  - `sleep` past a 1 s timeout, which becomes a timeout error and leaves no process behind;
  - a missing prefix or a missing `.giap-built`, which becomes `FgseaEngineUnavailable`;
  - a success, after which the temporary directory is removed;
  - launcher resolution order 1 to 4, and the `apptainer exec` argv;
  - the child environment, with `R_LIBS_USER` empty and `LC_ALL=C`.
- **T4 golden contract (slow, Slurm only; `@pytest.mark.slow`, `@pytest.mark.requires_r`).**
  - Turn #192's `inputs/gsea_stats.tsv.gz` and `gsea_sets.tsv.gz` into one job: `job_id =
    main:golden`, min 10, max 500, `--seed` GOLDEN_SEED from `tolerances.yaml`, and the
    `function.gsea` defaults.
  - Run `giap_fgsea.R` under the prefix launcher. Columns `pathway`, `pval`, `log2err`, `ES`,
    `NES`, `size` and `leadingEdge` must equal `outputs/fgsea.tsv.gz` **as strings, row for
    row** (exact; ADR-0038 item 7).
  - Repeat under the image launcher when the image is present. Without it, that test is
    deselected by marker, not skipped.
  - Submitted by `scripts/slurm/test_branch.sbatch -m requires_r`, at most 2 CPUs.
- **T5 tie order (slow, `requires_r`).**
  - On a ranking where 40 % of the scores tie, two calls that differ only in the order of the
    tied TUs give leading edges that follow that order.
  - Two identical calls give byte-identical `results.tsv`.
- **Recorded outputs (§11.4).** `scripts/slurm/record_fgsea.sbatch` (giap-record-fgsea, 2 CPUs,
  30 min) records `tests/data/recorded/fgsea/1.36.2/<input>/`, each with `stats.tsv`, `jobs.tsv`,
  `sets.tsv`, `command.txt`, `results.tsv`, `versions.txt` and `PROVENANCE.md` (the command, the
  lock sha256, the date, the host and the job ID). The inputs are:
  - `sim_main`: the main call of the M8-06 step fixture (`tests/helpers/function_step_data.py`),
    with collections cog, kegg, go_bp and go_slim;
  - `sim_retest`: that run's re-test call;
  - `ties`: the T5 ranking;
  - `unbalanced`: a ranking with 1 negative score and 499 positive ones, so that one set returns
    `pval` NA;
  - `tiny_p`: a set packed at the top, so that p is below 1e-20 and `log2err` is finite. Also, if
    the implementer can produce one, a case with `log2err` NA from `isCpGeHalf`. Otherwise the
    implementer states that it could not be produced, and lead decision 7 applies.

  Every input is synthetic and every file is under 1 MB.

---

## M8-17b: the engine switch, the step, Nextflow and preflight

### 7. The engine interface (pure, in `src/giap/function/gsea.py`)

```python
@dataclass(frozen=True, slots=True)
class GseaJob:
    job_id: str                                   # "main:<coll>" | "retest:<coll>:<term>"
    sets: Mapping[str, tuple[str, ...]]           # term -> member TU IDs, natural order
    min_size: int
    max_size: int

@dataclass(frozen=True, slots=True)
class EngineCall:
    tu_ids: tuple[str, ...]                       # perm order (§5)
    scores: F64                                   # same order, all finite
    jobs: tuple[GseaJob, ...]
    params: GseaConfig
    seed: int

@dataclass(frozen=True, slots=True)
class EngineRow:
    term: str; size: int; es: float; nes: float; p: float; log2err: float
    leading_edge: tuple[str, ...]                 # TU IDs, fgsea's order

EngineResult: TypeAlias = Mapping[str, tuple[EngineRow, ...]]   # job_id -> rows

class GseaEngine(Protocol):
    method: str                                   # written to the method column
    def __call__(self, call: EngineCall) -> EngineResult: ...

def gsea_table(genes, colls, size, cfg, alpha, rng, replicons=None, *,
               engine: GseaEngine = SIMPLE_ENGINE) -> pd.DataFrame: ...
```

- **`SIMPLE_ENGINE`** wraps the merged `fgsea_simple` port, with `method = "fgsea_simple"`, so the
  existing M8-05 tests keep running without R until the hardening pass (ADR-0038 item 6). It
  draws its permutations from `np.random.default_rng(call.seed)`. Exact numbers pinned by the old
  tests therefore move. The test-author re-pins them (lead decision 4).
- **`RFgseaEngine`** comes from §4, with `method = "fgsea_multilevel"`. It is the only engine the
  step uses.
- The module still imports nothing from `io/`. Injection keeps it pure.

### 8. What `gsea_table` does with the engine

**Main call.** Build one `GseaJob` per collection, `main:<collection>`. Each job holds the
collection's terms with `min_size ≤ size_tus ≤ max_size` (`function.gsea_size`) and
`size_tus < n`. A set covering every ranked TU has no complement, so it cannot be tested (lead
decision 3). A collection with no such term gets no job. If there is no job at all, **the engine
is not called**: the table is header-only and gets no new run flag. This keeps today's behaviour
for an empty collection.

**BH.** BH runs within each collection over the returned rows, exactly as today, through
`giap.stats.multitest.bh`. A NaN p keeps a NaN q and is excluded from the BH count, as today.

**Re-tests (one batched call).**
- For every GO collection and every significant term with significant descendants, compute
  `tus = members(set − ∪ desc sets)` as today.
- `k = 0`: no job, `redundant = True`. This follows the M8-05 lead decision 7.
- `k ≥ 1`: a job `retest:<coll>:<term>` with that single set, `min_size = 1` and `max_size = n`.
- All re-test jobs of all GO collections go into **one** second `EngineCall` with the same
  ranking and seed. With no re-test job there is no second call.
- `size_filtered` (`k < min_size`), `redundant = bh(p with p_re substituted)[i] > alpha` and
  "NaN p_re gives NA" keep their current rules.

**Row mapping.**

| Column | Source |
|---|---|
| `es`, `stat` | `EngineRow.es` |
| `nes`, `effect` | `EngineRow.nes` |
| `p_value` | `EngineRow.p` |
| `log2err` | `EngineRow.log2err` |
| `size_tus` | Python's count. The engine's `size` must equal it, otherwise the call fails |
| `method` | `engine.method` |
| `null_model` | `tu_set_permutation` |
| `n_perm` | `cfg.n_perm_simple` |
| `ci_low`, `ci_high` | NA, with flag `ci_unestimated` |

**Flags.**
- `unestimated` when `p_value` is NA. This is fgsea's unbalanced-signs case, where NES and q are
  NA too.
- `log2err_unestimated` when p is finite and `log2err` is NA (the M8-index row "GSEA `log2err`").
  This is a **new flag value**, so the `enrichment_gsea` schema's flags pattern gains it (§10).
- `size_filtered`, `genomic_cluster` and `detection_confounded` are unchanged.

**Leading edge back to genes.** For each TU in `EngineRow.leading_edge`, in fgsea's order, take
that TU's scoring genes in `tu_position` order (`rank.tu_genes[tu]`) and keep those in the term's
gene set. The cell is empty (written `NA`) when `es == 0` or the edge is empty. This is today's
rule; only the source of the TU list changes.

**Edge cases.**
- **NA scores.** These cannot reach the engine. `tu_ranking` keeps only genes with a non-NA
  `abundance_log2`, and `EngineCall` asserts that every score is finite. If one is not, Python
  raises `ValueError` before any R start. The R exit code 3 is a second guard.
- **Ties.** The seeded order goes into `stats.tsv` (§5). fgsea's ties warning goes into the event
  and is never an error. Imputed values tie by construction (M8-05 decision 9).

### 9. The step and the executors

`steps/function.py` builds `engine = r_fgsea_engine(resolve_launcher(os.environ),
tmp_root=Path.cwd())` **only when** `gsea_table` would call it. The engine is built lazily, so a
run with no abundance or with no testable set never needs R. It passes the engine as
`engine=engine`. `run_function` gains the same keyword and a default of `SIMPLE_ENGINE`, so its
existing unit tests are unchanged. Exceptions from the engine propagate and the step exits with
code 4. No table is written: the step writes after computing, as it does today.

- **Local executor.** It calls the same function. The launcher comes from `$GIAP_ROOT` (rule 3) or
  `GIAP_FGSEA_PREFIX`.
- **Nextflow.**
  - `workflow/conf/prebuilt_envs.config` adds `env { GIAP_FGSEA_PREFIX =
    "${GIAP_ROOT ?: projectDir/../local}/micromamba/envs/giap-fgsea" }`.
  - The apptainer profile adds `env { GIAP_FGSEA_IMAGE = "<GIAP_ROOT>/containers/<mangled>.img" }`.
    This runs after `prebuilt_envs` and takes precedence, because rule 2 is only consulted when
    rule 1 is unset. The apptainer profile therefore also unsets `GIAP_FGSEA_PREFIX`.
  - `GIAP_STEP` itself stays uncontainerised, so `apptainer exec` runs on the host, as Nextflow's
    own tasks do.
- **Resources.**
  - `function.nf` keeps its `FUNCTION_STEP` alias.
  - `base.config` gives `withName: 'FUNCTION_STEP'` 1 CPU, 6 GB and 4 h. That is about 4 GB for
    R on roughly 4,000 TUs × 3,000 GO sets plus the Python step, and the 4 h is above the
    per-call 1 h timeout times 2 calls.
  - On FarmShare that comes to 2 CPUs (4 GB per CPU, R14). There is no `maxForks`.
- **nf-test.**
  - `workflow/tests/function.nf.test` keeps its stub test. The stub does not touch R.
  - A new real-run test, `function_rfgsea.nf.test` with tag `requires_r`, runs the M8-06 fixture
    with `-profile farmshare` on Slurm only. It asserts a non-empty `enrichment_gsea.tsv` whose
    `method` is `fgsea_multilevel` on every row, and a byte-identical file on a second run.

### 10. Preflight and `giap check`

- `giap.annotate.tools.TOOLS` gains `ToolSpec("fgsea", "1.36.2", "MIT", env_file
  "workflow/envs/fgsea.yml", env_name "giap-fgsea", version_cmd ("R", "--vanilla", "--slave",
  "-e", "cat(as.character(packageVersion('fgsea')))"), module "function")`.
  - `ToolSpec` gains `image: str | None = None`, set to `FGSEA_IMAGE`.
  - `version_cmd` must not mention `Rscript` (check_no_rscript).
- `required_tools`: `"fgsea": "function" in mods`. See lead decision 8 for whether abundance
  should also be required.
- `missing_tools`: use `ref = pinned.get(tool.name) or tool.image`. The rest is unchanged. Under
  apptainer the tool is missing when the image file is absent; otherwise it is missing when
  `.giap-built` does not match the yml's sha256.
- A missing fgsea is reported exactly like a missing ostir: the same check ID, message shape and
  severity, and `giap run` refuses to start. One unit test per profile (A9).
- Schema `enrichment_gsea.schema.json`:
  - the flags pattern gains `log2err_unestimated`;
  - the `method` description drops "reserved for M8-11";
  - the enum is unchanged.
- `scripts/ci/check_no_rscript.py`: the docstring's "GIAP never shells out to R" becomes "src/
  starts R only through workflow/bin/giap_fgsea.R (ADR-0038)". The check itself is unchanged.

### 11. Tests (14b)

- **T6 engine use (unit, fakes).** A `RecordingEngine` stores every `EngineCall` and returns
  scripted rows.
  - The number of calls is 0, 1 or 2 as in examples 2, 6, 7 and 8.
  - The job IDs, sets, min and max sizes and the seed are as in §5 and §8.
  - The `perm` order is the one in `stats`.
  - The seed equals the second draw of `for_step("function.gsea", seed)`.
- **T7 replay (unit).** `ReplayEngine` reads `tests/data/recorded/fgsea/1.36.2/<input>/`. It
  checks that `write_call(call)` equals the recorded input TSVs byte for byte, then returns
  `read_results` of the recorded output. Fed with the M8-06 step fixture, `gsea_table` gives a
  table equal to a committed expected TSV. This is how CI exercises real R numbers without R.
- **T8 mapping and flags (unit).** Examples 3, 4 and 5; the size-mismatch error; and the
  whole-ranking set (example 10).
- **T9 step (unit).**
  - `function_step` with `GIAP_FGSEA_PREFIX` pointing to a missing directory exits with code 4
    when there is abundance, and runs fine without abundance.
  - The event `function.gsea_engine` is logged once per call.
- **T10 preflight (unit).** fgsea is required when `function` is in the modules and missing
  under each profile kind. It is not required under `test` or `stub`.
- **T11 calibration (Slurm only).** M8-15's `gsea` suite and attack 3 gain an `engine` parameter:
  - CI keeps running them on `SIMPLE_ENGINE` (fast replicates);
  - `-m "calibration and requires_r"` runs them on `RFgseaEngine` in
    `scripts/slurm/calibration_full.sbatch gsea_r`: 50 / 300 replicates, at most 2 R calls of a
    few seconds each, xdist with 8 workers, 8 CPUs, under 2 h for the full suite;
  - the criterion is unchanged: level, pooled, from `tolerances.yaml`, which is protected and not
    edited.

  This replaces M8-11's multilevel calibration (ADR-0038 Consequences).

---

## Examples (input → expected output)

1. **Golden contract.** #192's frozen `gsea_stats`/`gsea_sets`, one job `main:golden` with
   `[10, 500]`, seed 50000, 101 / 0 / 1000, through the prefix launcher → `results.tsv` columns
   `pathway…leadingEdge` equal `outputs/fgsea.tsv.gz` string for string. `versions.txt` has
   `fgsea\t1.36.2`.
2. **Two collections, one call.** Collections cog (14 in-range terms) and kegg (9), n = 300 ranked
   TUs, seed 7:
   - exactly one engine call;
   - `jobs.tsv` is `main:cog\t10\t500` then `main:kegg\t10\t500`;
   - `stats.tsv` has 300 rows in `for_step("function.gsea", 7).permutation(300)` order;
   - `--seed` is the next `integers(1, 2**31 - 1)` draw of that generator;
   - 23 rows, all with `method = fgsea_multilevel` and `n_perm = 1000`.
3. **Unbalanced set.** The engine returns `pval = NA, NES = NA, log2err = NA, ES = 0.81` for
   `K00001` → row `p_value`, `q_value`, `nes`, `effect` and `log2err` NA, `es = 0.81`, flags
   `ci_unestimated;unestimated`. The other kegg q values equal `bh` over the non-NA p values.
4. **log2err not estimated.** `pval = 3.2e-15, log2err = NA` → `p_value = 3.2e-15`, `log2err` NA,
   flags `ci_unestimated;log2err_unestimated`. q is computed normally.
5. **Leading edge to genes.**
   - The engine returns `leadingEdge = ("TU_7", "TU_2")`. TU_7's scoring genes are
     `[b0007a, b0007b]` and TU_2's are `[b0002]`. The term's genes are `{b0007b, b0002, b0099}`.
   - Result: `leading_edge = "b0007b;b0002"`. b0099 is not in a leading-edge TU, and b0007a is
     not in the set.
6. **Batched redundancy.**
   - go_bp parent `GO:P` has q = 0.001; its significant child `GO:C` has q = 0.002;
     `set(GO:P) − set(GO:C)` covers 12 TUs. go_mf parent `GO:Q` has `set − desc` covering 0 TUs.
   - Result: exactly two calls. The second has the single job `retest:go_bp:GO:P` with
     `[1, n]` and the same seed. `GO:Q` gets `redundant = True` without a job.
   - The re-test returns p = 0.2. The p vector with 0.2 substituted gives q > 0.05, so `GO:P` gets
     `redundant = True`.
7. **No abundance.** `abundance_scale` and `abundance_log2` are all NA → zero engine calls, the
   engine is never constructed, the table is header-only, and the run flag is
   `layer_missing:abundance`. This holds even when `GIAP_FGSEA_PREFIX` points nowhere.
8. **Every set out of range.** Every term has 3 to 8 TUs against `[10, 500]`, or the collection is
   empty → zero calls, the table is header-only, and there is no run flag.
9. **R failure.**
   - The script exits 1 with stderr `Error in fgseaMultilevel(...) : ...` →
     `FgseaEngineError("R fgsea failed (exit 1, 2 jobs): …Error in fgseaMultilevel…")`. The step
     exits with code 4 and no `tables/enrichment_*.tsv` is written.
   - A 1 s timeout → `FgseaEngineError("R fgsea timed out after 1 s")` and no child process is
     left.
10. **Set covering the whole ranking.** n = 40 ranked TUs, and term `C:J` has all 40 → no job for
    it. Its row has `size_tus = 40`, NA es, nes, p and q, and flags `ci_unestimated;unestimated`
    (lead decision 3).

## Acceptance criteria (each checked by a test)

- **A1** (T4) On #192's inputs with GOLDEN_SEED, `giap_fgsea.R` reproduces `fgsea.tsv.gz`
  exactly. This runs as a slow test on Slurm, and its job ID goes in the PR.
- **A2** (T5, nf-test real) Two runs with the same seed give byte-identical `results.tsv` and
  `tables/enrichment_gsea.tsv`.
- **A3** (T1) `fgsea.yml` is clean under `env_problems`, has the pins in §1 and a committed
  explicit lock with no forbidden data package. `check_no_rscript.py src` and `check_imports.py`
  pass.
- **A4** (T2, T3) The wrapper turns every failure (non-zero exit, timeout, missing env, wrong
  fgsea version, malformed output, size mismatch) into `FgseaEngineError` or
  `FgseaEngineUnavailable`, with the messages above. The temporary directory is always removed.
- **A5** (T6) The number of engine calls is 0, 1 or 2 as in examples 2 and 6–8. The seed is
  `for_step`'s second draw. `stats.tsv` is in `perm` order.
- **A6** (T7) The replay of `sim_main` and `sim_retest` through `gsea_table` equals the committed
  expected table.
- **A7** (T8) Examples 3, 4, 5 and 10 hold, and `log2err_unestimated` validates against the
  schema.
- **A8** (T9) The step exits with code 4 on an engine error, and runs without R when there is no
  abundance.
- **A9** (T10) preflight and `giap check` list `fgsea` as missing under farmshare and under
  apptainer when it is absent, and never under `test` or `stub`.
- **A10** (T11) The `gsea` level suite and attack 3 pass on the R engine on Slurm. The job ID goes
  in PROGRESS.md.
- **A11** The ORA path is untouched: `enrichment_detection.tsv` is byte-identical before and after
  this change on the M8-06 fixture, and no goseq appears anywhere in the pipeline.

## Line estimate

| Part | Files | Code lines |
|---|---|---|
| 14a | `fgsea.yml` (12), `giap_fgsea.R` (≈ 120), `io/fgsea_r.py` (≈ 150), `record_fgsea.sbatch` (≈ 25), `build_envs` entry (≈ 5) | **≈ 290** (+ ≈ 230 test lines; the lock and recorded files are data) |
| 14b | `gsea.py` engine types and switch (≈ 100 changed), `steps/function.py` (≈ 15), `tools.py` + `environment.py` (≈ 15), configs (≈ 12), schema (2), `function_rfgsea.nf.test` (≈ 30), calibration parameter (≈ 10) | **≈ 180** (+ ≈ 220 test lines) |

14a lands first, and its T4 job is the evidence that the wrapper is right. 14b depends on 14a.

## Lead decisions needed

1. **Issue number.** M8-17 is taken by `M8-17-sim-function-determinants.md`. Should this be
   M8-17 (17a and 17b), with the M8 index row added, or does the lead want to renumber the
   simulator issue? The filename stays as instructed until then.
2. **R inside the step (recommended), or a separate Nextflow process.**
   - The recommendation is to start R from the Python `function` step, through the wrapper. That
     keeps §6.4 (both executors, the same step function, byte-identical) and keeps the
     interactive re-test round as a plain function call.
   - The alternative is the OSTIR pattern: a `gsea_jobs` Python step, an `FGSEA` process, a
     `gsea_retest_jobs` step, `FGSEA_RETEST` and then `function`. That is five tasks. Nextflow
     would own the container, but the local executor could not run GSEA at all (§6.4 forbids
     external steps there).
   - Note also that the brief said the engine "is injectable". On main it is not: tests
     monkeypatch `fgsea_simple`. §7 adds the injection.
3. **A set covering every ranked TU.** This spec writes an `unestimated` row with NA stats,
   mirroring the ORA's whole-universe rule. Should such a set get no row instead? Under default
   sizes it can only happen when n ≤ 500.
4. **Re-pinning old tests.** `SIMPLE_ENGINE` seeds itself from `call.seed`, so exact numbers that
   current M8-05 and M8-06 tests pin for `fgsea_simple` move. Should the test-author re-pin them in
   14b? The other option is to keep the old rng-threading inside `SIMPLE_ENGINE` and avoid the
   churn, at the cost of a second, divergent seed path.
5. **R failure policy.** Any R error fails the step with exit code 4, and ORA and the category
   summary are lost with it. Alternatively, GSEA could fall back to a header-only table with a run
   flag such as `gsea_engine_failed`. That would need a new run-flag value, and it hides a real
   failure.
6. **Timeout.** A fixed 3,600 s per call in code, because `defaults.yaml` is protected and adding
   a key needs an ADR. Is a constant acceptable, or should `function.gsea.timeout_s` go into a
   small ADR?
7. **`isCpGeHalf`.** If the `log2err`-NA case cannot be produced on synthetic data for the
   recording, should example 4 stay a fake-only test?
8. **When fgsea is required.** Should it be required whenever `function` is in the modules
   (simple, as specified), or only when abundance is also given (`required_tools` would then need
   inputs, not just cfg)? As specified, a no-abundance `function` run still demands the env at
   preflight.
9. **Apptainer under a non-container GIAP_STEP.** Running `apptainer exec` from inside a Python
   task works on FarmShare, but it is not Nextflow-managed. Is that acceptable, or should fgsea
   stay conda-only under every profile (like uniop and ostir) and drop the image pin from
   ADR-0038 item 1?
10. **M8-13.** Should the fgsea half of M8-13 (port vs golden) be closed as superseded by A1 here,
    with its Wallenius half kept as is?

## Lead decisions (2026-10-22)

1. **Number.** The issue is M8-17 (M8-14 is the simulator issue), split into M8-17a and M8-17b.
2. **Where R runs.** R is started from inside the Python `function` step through `giap.io.fgsea_r`. Both executors share one code path.
3. **A set covering every ranked TU.** It gets an `unestimated` row with NA statistics, matching the ORA's whole-universe rule.
4. **Old test numbers.** The M8-05/M8-06 tests that pin port numbers are re-pinned by the lead after M8-17b's implementation. The test author writes new tests only (role guard).
5. **R failure.** Any R error fails the step with exit 4, the same as every other tool failure. No partial GSEA table.
6. **Timeout.** The 3,600 s code constant is accepted; it moves to config in the hardening pass if it is ever needed.
7. **The `log2err` NA case.** Tested with fakes and recorded outputs only.
8. **When fgsea is required.** Preflight requires the fgsea env only when `function` is enabled **and** abundance is given, which is when GSEA runs.
9. **Apptainer.** fgsea is conda-only under every profile for now. The BioContainers pin in ADR-0038 item 1 is deferred to the hardening pass (backlog). The step never calls `apptainer exec` itself.
10. **M8-13.** Its fgsea half is superseded by the M8-17a contract test; its Wallenius half stays.

## Lead decisions on the M8-17a tests (2026-10-22, after 9d1b517)

1. **Golden contract.** Python restricts each set to the ranked TUs before handing it to R, in file order, as `gsea.py` already does. The R script's exit 3 for unranked members stays as an input check. So the contract compares 119 rows: SET_027 drops below minSize.
2. The §7 types `GseaJob`, `EngineCall`, `EngineRow` and `EngineProvenance(r, fgsea, warnings)` are added to `giap.function.gsea` in M8-17a, without the engine switch, which stays in 17b.
3. **Results header.** `job_id pathway size ES NES pval log2err leadingEdge`, as pinned.
4. **Warnings.** The R script replaces embedded newlines in warnings with a space.
5. **Provenance.** The engine object keeps the `EngineProvenance` of its last call as `engine.provenance`, which is `None` before the first call. 17b reads it for the `function.gsea_engine` event.
6. **Missing R.** The wrapper checks for `<prefix>/bin/R` and raises `FgseaEngineUnavailable` when it is absent.
7. **The `requires_r` marker.** It stays registered in `tests/contract/conftest.py` for 17a and moves to `pyproject.toml` in 17b with the calibration tests.
8. **EXPECTED_ENVS.** The lead added `fgsea` to `EXPECTED_ENVS` and `CONDA_ONLY` in `tests/unit/annotate/test_tool_env_files.py`.
9. **check_no_rscript.** Its test passes once #192 is merged into this branch.
10. **No PROVENANCE line in the lock.** Provenance lives in the recordings' README and the versions file.
11. **Fixtures.** `by_hand/` stays as a permanent format fixture. The five real recordings are made on Slurm by `scripts/slurm/record_fgsea.sbatch` once the env is built.

## Lead notes for M8-17b (after the 17a implementation)
- Move fgsea into `TOOLS`. Then remove the lead's `{"fgsea"}` exemption in `test_every_env_file_has_a_registry_entry`, and add fgsea to `ALL_TOOLS` and `ALL_BUILDS`.
- The step fixture yields six collections (cog, kegg, go_bp, go_mf, go_cc, go_slim), while `sim_main` was recorded for the four §6 collections. T7 must either replay only those four or re-record with six.
- `load_call` rebuilds `params` from `GseaConfig()` defaults, but sim_main and sim_retest used `n_perm_simple` 200. T7 must pass the recorded params.
