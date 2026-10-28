# Fix: BacTermFinder genome_scan arguments and working directory

## Problem
`workflow/modules/local/bactermfinder/main.nf` runs `genome_scan.py genome.fa bactermfinder`. At
the pinned commit `024870447f6288ebc9d29b314e957a6dfdfeb5f6` (ADR-0008) the script takes four
positional arguments, `<fasta> <step> <prefix> <batch>` (`int(sys.argv[2])`, `int(sys.argv[4])`),
so the process fails. It also reads its four `deep-bactermfinder_*_saved_model.h5` models and runs
`python iLearnPlus/util/FileProcessing.py` relative to the working directory, and it writes
`Sample.csv`, `df_sample.fasta`, `output_sample/` and its CSVs there. The source lives in
`$CONDA_PREFIX/share/bactermfinder` (envs/bactermfinder.yml) or `/opt/conda/share/bactermfinder`
(containers/bactermfinder/Dockerfile). That directory is shared by every task and read-only in
an Apptainer image. Smoke job 1895829 ran the 4-argument form in a source dir (106 windows).

## Required behaviour
1. Arguments: `genome.fa 3 out 10000`, which is the README example at the pinned commit. GIPPA_SPEC.md
   §7.1.1 says "per its README" and no GIPPA spec or methods document sets another step.
   The prefix `out` has no `.`, because the script uses `prefix.split(".")[0]`.
2. Working directory: the script makes a fresh `btf_work/` in the task dir and `cd`s into it.
   It symlinks `genome_scan.py`, `iLearnPlus` and the four `*.h5` models from
   `$BACTERMFINDER_HOME`, and `genome.fa` from the task dir. Then it runs
   `python genome_scan.py genome.fa 3 out 10000 > log.out`. It never `cd`s into
   `$BACTERMFINDER_HOME` and never writes there. The genome argument is the bare name
   `genome.fa`, because the script embeds it in its output names.
3. Outputs: the result table is `out` + `genome.fa` + `_mean.csv` = `outgenome.fa_mean.csv`.
   Its columns are `SampleName,probability_binary,probability_ENAC,probability_PS2,probability_NCP,probability_mean`,
   as in the smoke log. The table is copied to `bactermfinder/outgenome.fa_mean.csv`, and
   `bactermfinder/` holds no other `.csv`. The step reader (`gippa.steps.terminators._bactermfinder`,
   M2-15 review item 5) rejects a directory with several CSVs. The intermediate CSVs stay in
   `btf_work/`, outside the `bactermfinder` output glob. `output:` keeps `path 'bactermfinder', emit: dir`
   and `versions.yml`. The `BACTERMFINDER_HOME` default and `GIPPA_COMMIT` version line stay the same.
4. The stub keeps writing an empty `bactermfinder/` (nf-test annotate example 5).
## Out of scope (follow-up)
`parse_bactermfinder` (provisional, M2-15) expects contig/start/end/strand columns. The real table
encodes these in `SampleName` (`<n>_<contig>_<start>_<end>_<strand>`), so the parser needs its own
recorded-fixture PR (docs/backlog.md M2-15).

## Tests
`tests/unit/test_fix_bactermfinder_argv.py`: static checks on the `script:` block. It also runs
the rendered block in bash against a fake, read-only `BACTERMFINDER_HOME` whose `genome_scan.py`
mimics the pinned script's file I/O.
