Tool runs no longer import Python packages from the user site (`~/.local`): every workflow task,
including Docker and Apptainer tasks, and every `gippa setup` tool call run with
`PYTHONNOUSERSITE=1`, so a user-installed ViennaRNA or biopython can no longer change OSTIR and
other tool results. OSTIR's `versions.yml` now records the ViennaRNA version it actually used.
