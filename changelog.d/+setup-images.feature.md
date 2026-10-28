`gippa setup --images --profile apptainer|docker` pulls the pinned container images a config
needs (`--images all`: every pinned image), checks each pulled digest against its pin and records
the images in `db_manifest.json`; `--cache-dir` chooses the Apptainer image cache. `gippa setup`
itself now runs the Bakta and antiSMASH downloaders through Docker when no environment or cached
image is available, so it needs neither `$GIPPA_ROOT` nor a source checkout. `gippa setup --check`
reports a database installed from another release as `installed <old>, this release pins <new>`.
