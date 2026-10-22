# OSTIR 1.1.3 on three sim100k novel candidates (M7-08 rule 9)

- Tool: OSTIR 1.1.3 with ViennaRNA 2.6.4 (`run.log`), from the project's official micromamba env
  `giap-ostir` built from `workflow/envs/ostir.yml` and its lock (`versions.txt`: sha256 of
  `ostir.explicit.txt` and `ostir.pip.txt`, the same env as the `../sim100k` recording).
- Run: Slurm job 1865755 (`giap-record-ostir-novel`, 4 CPUs, node wheat-04, 2 s wall),
  `sbatch scripts/slurm/record_ostir_novel.sbatch <worktree>` at giap commit 336cc1b
  (branch m7/novel-steps). Command as the OSTIR module runs it: `command.txt`.
- Input (`ostir_input.csv`): `scripts/dev/record_ostir_novel_input.py` simulates giapsim
  `SIM100K` (seed 20261014) with the M7-08 test configuration (M3-19 proteomics, planted novel
  ORFs 3 intergenic / 2 cds_outframe / 2 antisense, 2 N-terminal extensions;
  `tests/helpers/novel_steps_fixture.py`), runs `giap run --executor local` and keeps three of
  the 9 rows `novel_candidates` wrote to `ostir_input_novel.csv`: the first multi-start row
  (`NOVEL_chr1_3872_p`, start 51, end 150) and the first two others
  (`NOVEL_chr1_2635_p`, single start; `NOVEL_chr1_11220_m`, start 51, end 153), in file order.
- Output (`ostir.csv`): as OSTIR wrote it, unedited. OSTIR also reports start codons that are
  not in frame with `start` (e.g. 66, 106, 118); `parse_ostir_starts` drops them.
- Raw copy with `SHA256SUMS` and `runtime.txt`: `$GIAP_ROOT/recordings/ostir/sim100k_novel/`.
