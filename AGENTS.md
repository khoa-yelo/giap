# AGENTS.md: how agents build, test and review GIAP

## Source of truth
- GIAP_SPEC.md (read-only). Accepted ADRs in docs/adr/ override it.
- Rules R1–R14 (GIAP_SPEC.md §0.2) bind every agent and subagent.

## Setup (FarmShare; inside a Slurm allocation, never on a login node)
source scripts/farmshare_env.sh               # every shell and every sbatch script
bash scripts/bootstrap_farmshare.sh           # once: uv, Nextflow, nf-test, gh (user space)
uv sync --frozen --all-extras
uv run python scripts/probe_env.py            # writes .giap/capabilities.json
sbatch scripts/slurm/build_envs.sbatch <tool> [...]   # micromamba envs (workflow/envs/*.yml) and Apptainer pulls
sbatch scripts/slurm/giap_setup.sbatch        # giap setup --bakta light,full --genomad --promotech-models --go --kegg --cog --contaminants

## Test commands
- every commit:   uv run pytest -m "not slow and not calibration_full and not e2e" -n auto
- calibration:    uv run pytest tests/calibration -m calibration -n auto      (full: -m calibration_full)
- snapshots:      uv run pytest tests/snapshot --mpl
- integration:    uv run pytest tests/integration ; nf-test test workflow/tests
- end-to-end:     sbatch scripts/slurm/weekly.sbatch    (never in the session allocation)
- pipeline run:   sbatch scripts/slurm/giap_run.sbatch <args>   (-profile farmshare,apptainer)
- benchmarks:     uv run pytest benchmarks --benchmark-only
- lint and types: uv run ruff check . && uv run ruff format --check . && uv run mypy src/giap

## Conventions
- src layout, type hints everywhere, pure functions for statistics; file I/O only in io/ and steps/.
- 0-based half-open coordinates internally; convert only in writers (spec §6.5).
- Randomness only through giap.rng.for_step(step_name, seed) (spec §6.6).
- Never import GPL/AGPL/unlicensed tools (scripts/ci/check_imports.py).
- Commit trailer `Agent-Role: <role>`; conventional commits; PRs ≤ ~400 lines.
- Shared cluster (R14): jobs named `giap-*`, ≤ 128 CPUs / 4 GPUs / 20 jobs, large files only under $GIAP_ROOT.

## Definition of done (module): spec §14.9
## Where things live: specs docs/specs/, methods docs/methods/, ADRs docs/adr/, backlog docs/backlog.md, state PROGRESS.md
