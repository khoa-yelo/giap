# Recorded UniOP run on giapsim SIM100K (M2-12, lead decision 3 after the PR #49 review)

- Tool: UniOP, git commit `3ed14dfc2a5d69b993d99a26e50c1f8e9723c087` (pin 3ed14dfc; ADR-0008), from
  the built `giap-uniop` micromamba env (`$CONDA_PREFIX/share/uniop/GIAP_COMMIT`). Run as a separate
  process; UniOP is never imported or vendored (only its output is recorded here).
- Env lock: `workflow/envs/locks/uniop.explicit.txt` (sha256
  `a16e1731537f901edfff5c1b7604b40e05de40ba8a062c7d4647ef54b7ff463b`) and
  `workflow/envs/locks/uniop.pip.txt` (sha256
  `3a2d491cf84fd08caf79b73d9ca48ea06442975cb41f3f5cdc17a3f8b1b5eaff`); spec `workflow/envs/uniop.yml`.
- Input: `uniop_cds.gff` and `uniop_map.tsv` written by the `uniop_input` step
  (`giap.annotate.uniop_adapter.write_prodigal_gff`) after `parse_bakta` on
  `giapsim.simulate(giapsim.SIM100K, giapsim.SIM100K_SEED)`; 87 genes on `chr1` and `plasmid_1`.
  - `uniop_cds.gff` sha256 `1d9c91e955d816df7b3e23890661af3e430f09317465dc6f479d233b7049053e`
  - `uniop_map.tsv` sha256 `1b7c6475770e931f4054fa8b7cf061c14ea48759c4e52f88233704e507f3d3fe`
- Command (`giap.annotate.uniop_cmd.uniop_command`, as in `workflow/modules/local/uniop/main.nf`), run
  in the env with 2 CPUs (`OMP_NUM_THREADS=2`) after `mkdir -p out`:

      python $CONDA_PREFIX/share/uniop/src/UniOP -f uniop_cds.gff --gff_source prodigal_gff -t out

  Exit 0; stderr only a SciPy/NumPy version UserWarning. UniOP wrote `out/uniop.pred` and
  `out/uniop.operon`; `uniop.pred` is copied byte-for-byte (sha256 `ee4fcc9d09db8a34471c4ad03a2a632d4b9eb29048fbd3ddf7663e1b51b82853`).
- Recorded 2026-10-14 on Stanford FarmShare by the M2-12 implementer.
