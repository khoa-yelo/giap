`gippa setup --uniop` builds UniOP for every runtime and asks first:
- the terms notice is always printed; without `--accept-uniop-terms` setup asks
  `Build UniOP for your own use? [y/N]` on a terminal and builds nothing otherwise (no is the
  default; `--uniop` then exits 2, `--images` warns and exits 0);
- apptainer builds `uniop.def` into a `0600` `.sif` (`--fakeroot`, then
  `--fakeroot --ignore-fakeroot-command`); plain `conda` creates `<cache-dir>/envs/gippa-uniop`
  from the lock; `farmshare*` profiles check and record the prebuilt env;
- an explicit `setup --uniop` installs UniOP only, and a rerun with a current build does nothing;
- `db_manifest.json` records the build as `local_build`; preflight accepts it only when its
  runtime, key and image ID, `.sif` sha256 or env commit still match, and otherwise names the
  setup command;
- the run manifest's `uniop` entry carries `local_build`, `runtime`, `source_sha256` and
  `container_digest`.

Stub runs (`-stub-run`) write `operon_meta.json` with `operon_source: uniop` when UNIOP ran, so
the stub run's manifest records the operon source.
