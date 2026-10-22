# M8-12: The R golden environment, the frozen inputs, `make_golden_r.R`, and the committed outputs

Milestone M8 · track G (golden data) · **pass 2** (the env build is submitted in wave A) · depends
on: nothing in `src/` (inputs come from giapsim and plain numpy); M8-02's `prepare` order is
mirrored, not imported · estimate ≈ 220 code lines (R ≈ 110, Python ≈ 70, shell/YAML ≈ 40)
(+ ≈ 120 test lines) · not statistical · **Slurm**

## Spec sections implemented

- §11.5 "`scripts/golden/make_golden_r.R` runs in a conda environment (`r-base`,
  `bioconductor-fgsea`, `bioconductor-goseq`; conda-forge and bioconda). It writes: fgsea
  multilevel results (ES, NES, p, log2err) for frozen ranked TU scores and gene sets; goseq
  Wallenius p-values for frozen sets with a supplied PWF equal to p_detect. Outputs go to
  `tests/data/golden/r/` with the R and package versions. `golden-r.yml` (manual) can regenerate
  them and open a pull request."
- §3.7 "R golden references (fgsea, goseq) | yes, conda R env | or a short batch job | on demand";
  §12 `golden-r.yml` "manual | regenerate the R golden data and open a PR (may also be done on
  FarmShare with the conda R env)"; §3.2 sources (conda-forge and bioconda only, never `defaults`);
  §3.8 batch scripts; R5, R6, R11 (≤ 5 MB per file, `tests/data/` ≤ 50 MB, gzip text fixtures).

## Licences and the boundary (index §Licences)

fgsea (MIT), goseq (LGPL ≥ 2), BiasedUrn (GPL-3) and R (GPL-2+) are installed **only** in
`giap-golden-r`. No `src/` file imports, vendors or shells out to R; `scripts/ci/check_imports.py`
adds `rpy2` to FORBIDDEN, and a new CI check fails if any file under `src/` contains `Rscript`.
The committed outputs are numbers (not covered by those licences). The golden **inputs are
synthetic** — random gene IDs and random set memberships — so they contain no KEGG, GO or COG
content (KEGG's academic-use terms, index §KEGG).

## What is added

| File | Content |
|---|---|
| `scripts/golden/golden_r.yml` | env `giap-golden-r`, channels `[conda-forge, bioconda]`: `r-base=4.5.3`, `bioconductor-fgsea=1.36.2`, `bioconductor-goseq=1.62.0`, `r-biasedurn=2.0.12`, `r-data.table` (fgsea dep, solver-pinned) — exact pins confirmed by the first solve and frozen in the lock |
| `scripts/golden/golden_r.explicit.txt` | `micromamba env export --explicit` lock (the reproducibility record) |
| `scripts/golden/make_golden_inputs.py` | writes the frozen inputs from `GOLDEN_SEED` (50000, `tolerances.yaml`) with numpy only |
| `scripts/golden/make_golden_r.R` | reads the inputs, runs fgsea and goseq, writes the outputs and `VERSIONS.txt` |
| `scripts/slurm/golden_r.sbatch` | `giap-golden-r` job: builds the env if its marker is stale (the `build_envs.sh` marker pattern: sha256 of the yml), then runs both scripts |
| `.github/workflows/golden-r.yml` | manual (`workflow_dispatch`): setup-micromamba from the lock, run both scripts, open a PR with the diff (never pushes to main) |
| `tests/data/golden/r/` | `inputs/*.tsv.gz`, `outputs/*.tsv.gz`, `VERSIONS.txt`, `PROVENANCE.md` |

1. **Inputs** (`make_golden_inputs.py`, deterministic from `GOLDEN_SEED`):
   - **GSEA**: 2,000 "TU" scores, named `TU_0001…`, a mixture designed to exercise the port —
     N(0, 1) background, 5 % exact ties (two values repeated), a few zeros; already centred
     (M8-05 decision 8), in TU-ID order (M8-05 step 2). 120 sets of sizes 10–500 (log-uniform),
     of which 20 are planted toward the top or bottom with shifts chosen so that true p spans
     10⁻² … 10⁻³⁰, 5 are tiny-but-legal (size 10), and 2 contain names absent from the scores.
   - **ORA**: a universe of 3,000 genes with p_detect from Beta(2, 2), a few at exactly 0 and 1;
     detected ~ Bernoulli(p_detect); 150 categories of sizes 5–500, including: ω ≫ 1 (members
     drawn from high p_detect), ω ≪ 1, k = n, k = 0, a category equal to the universe minus one
     gene, strong planted enrichments giving p down to ~10⁻⁴⁰, and 30 % of genes in no category.
2. **R** (`make_golden_r.R`, `set.seed(GOLDEN_SEED)` before each call):
   - `fgseaMultilevel(pathways, stats, sampleSize = 101, minSize = 10, maxSize = 500, eps = 0,
     nPermSimple = 1000, scoreType = "std", nproc = 1)` → `pathway, pval, log2err, ES, NES, size,
     leadingEdge` (`;`-joined); plus `calcGseaStat` ES for every set (deterministic reference).
   - `goseq(pwf, gene2cat = <long table>, method = "Wallenius", use_genes_without_cat = TRUE)`
     with `pwf = data.frame(DEgenes, bias.data = 1, pwf = p_detect, row.names = genes)` →
     `category, over_represented_pvalue, under_represented_pvalue, numDEInCat, numInCat`.
   - Numbers written with `format(x, digits = 17)` (round-trip exact), gzip, ≤ 1 MB total.
   - `VERSIONS.txt` = `sessionInfo()` + `packageVersion()` of fgsea, goseq, BiasedUrn, + the lock's
     sha256.
3. **`PROVENANCE.md`** (the §11.4 recorded-output template): command, Slurm job ID, host, date,
   env lock sha256, git commit of the scripts, `GOLDEN_SEED`.
4. **Regeneration** is the only way the outputs change: a PR that edits `tests/data/golden/r/`
   without a matching `PROVENANCE.md` update fails a CI check. Re-generation needs a reason
   (package upgrade → ADR, as for tool upgrades in §11.4).

## Degenerate cases (what the scripts must refuse)

| Situation | Behaviour |
|---|---|
| env solve picks a different fgsea/goseq version than pinned | build fails (the `build_envs.sh` version-check pattern: `packageVersion` must equal the pin) |
| `conda.anaconda.org` unreachable | job fails after one retry; PROGRESS.md › Blocked (R6); no other source |
| outputs > 5 MB | script fails (R11) |
| an input set with < 10 members after name matching | kept in the inputs; R drops it (`minSize`), and the Python side must drop the same set (tested in M8-13) |
| goseq warns about genes without categories | expected (`use_genes_without_cat = TRUE` keeps them); the warning text is captured in `VERSIONS.txt` |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `sbatch scripts/slurm/golden_r.sbatch` on FarmShare | env built (first time, ≤ 1 h, ≤ 4 CPUs, 16 GB), outputs written, `VERSIONS.txt` lists fgsea 1.36.2, goseq 1.62.0, R 4.5.x |
| 2 | the job run twice | byte-identical outputs (fixed seeds, `nproc = 1`) |
| 3 | `make_golden_inputs.py` run twice | byte-identical inputs |
| 4 | `golden-r.yml` dispatched | a PR whose diff in `outputs/` is empty (same env lock) |
| 5 | `python scripts/ci/check_imports.py src` on a file with `import rpy2` | exit 1 |

## Acceptance criteria

- [ ] The env builds on FarmShare from the yml in a `giap-golden-r` job (job ID in PROGRESS.md);
  the lock is committed.
- [ ] `tests/golden/test_golden_r_files.py`: inputs and outputs exist, parse, have the expected
  set/category counts, total < 5 MB, and `VERSIONS.txt` names the pinned versions.
- [ ] Re-running inputs and the R script gives identical bytes (example 2, 3), checked once in the
  job and recorded.
- [ ] `check_imports.py` and the new `Rscript` guard extended and tested.
- [ ] `golden-r.yml` exists, is `workflow_dispatch` only, and opens a PR (never pushes to main).

## Hardening owed

- When bioconda publishes fgsea 1.38 / goseq 1.64 for R 4.6, an upgrade ADR re-generates the
  golden data and re-runs M8-13 (the port must then be re-checked against the new version).

## Questions for the lead

1. Env location: `scripts/golden/golden_r.yml` built by its own sbatch (proposed), or
   `workflow/envs/golden_r.yml` with a `TOOLS` registry entry so `build_envs.sbatch golden_r` works?
   The registry is for pipeline tools; `build_envs.sh all` would start building R for every user.
2. Create `golden-r.yml` now (proposed; §12 lists it) or keep FarmShare as the only route?
3. Pins: bioconda's fgsea 1.36.2 / goseq 1.62.0 on R 4.5 (proposed; the only versions bioconda has
   on 2026-10-21), not Bioconductor's current 1.38.0 / 1.64.0.
