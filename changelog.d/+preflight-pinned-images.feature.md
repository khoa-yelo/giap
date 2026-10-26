Preflight now treats a pinned container image as the install of its tool under the `docker` and
`apptainer` profiles: under docker no local environment is needed, and under apptainer the image
must be in the image cache (`NXF_APPTAINER_CACHEDIR`, else `$GIAP_ROOT/containers`, else
`<db_dir>/containers`), so a pip install without `$GIAP_ROOT` works. A tool that has only a conda
environment is reported missing under `docker` or plain `apptainer`, with the reason, instead of
running on the host without its environment; a missing image names `giap setup --images`.
