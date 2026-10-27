`gippa setup` no longer needs `--db-dir`: it uses `$GIPPA_DB_DIR`, else
`$XDG_DATA_HOME/gippa/db` (`~/.local/share/gippa/db`), and records the directory in
`~/.config/gippa/db_dir`, which `gippa run`, `gippa preflight --profile` and `gippa makedb` read
when neither `--db-dir` nor `$GIPPA_DB_DIR` is given. Before downloading, setup prints the plan
(each item and image with its size, the apptainer build space, the total, the free space and
user quota of each filesystem), refuses with exit 2 when less than 1.1 x the need is free, warns
about a default DB directory in a shared home, and asks `Proceed? [Y/n]` on a terminal (`--yes`
skips it). A full disk or quota during a download or pull marks the item `failed` and leaves no
partial file; `apptainer pull` runs with `--disable-cache` unless `APPTAINER_CACHEDIR` is set.
