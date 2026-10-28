# Fix: isolate tool runs from user-site Python packages (PYTHONNOUSERSITE)

## Problem
The Python 3.10 tool environments (gippa-ostir, gippa-promoter_calculator, gippa-terminatornet,
gippa-cascadia) import packages from `~/.local/lib/python3.10/site-packages` (ViennaRNA 2.7.1,
biopython 1.83, ...) because nothing sets `PYTHONNOUSERSITE`. OSTIR results differed on 6 of
~9,000 genes depending on the user's `~/.local`. Apptainer runs have no `--cleanenv` and see
`$HOME` as well, so containers are affected too.

## Required behaviour
1. `workflow/nextflow.config` has a top-level `env { }` scope (outside `profiles`, so it applies
   under every profile) that sets `PYTHONNOUSERSITE = '1'`.
2. Under the `docker` profile the variable reaches the container: `docker.envWhitelist` lists
   `PYTHONNOUSERSITE` (comma-separated, next to `TMPDIR`).
3. `gippa.setupdb.run.ToolRunner.run` starts the child with an explicit `env=` mapping that is a
   copy of `os.environ` (every other variable passes through unchanged) with
   `PYTHONNOUSERSITE=1`; a pre-existing value (e.g. `0`) is overridden. `os.environ` itself is
   not modified. The return value is still the child's return code.
4. `scripts/farmshare_env.sh` exports `PYTHONNOUSERSITE=1`, so every session shell and every
   sbatch job (which sources it) runs with the user site disabled.
5. The OSTIR module (`workflow/modules/local/ostir/main.nf`) records the ViennaRNA Python binding
   version in `versions.yml` as a `viennarna:` line read from `RNA.__version__` in the task's
   Python; the stub writes a `viennarna:` line as well.

## Out of scope
Changing which ViennaRNA version is pinned; `--cleanenv` for Apptainer.

## Tests
`tests/unit/test_fix_pythonnousersite.py` (static config/module text checks, a monkeypatched
`subprocess.run` for ToolRunner, and `farmshare_env.sh` sourced in a bash subprocess with a fake
`HOME`).
