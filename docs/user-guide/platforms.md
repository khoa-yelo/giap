# Platforms

Where GIAP is tested and which execution profiles are supported.

!!! note "Stub"
    Filled by M11-09 (platforms).

## Image references and registries

`workflow/conf/containers.config` pins each image as a bare `<image>@sha256:<digest>` reference,
with no `docker://` prefix: Docker takes it as is, and Nextflow adds `docker://` for Apptainer.
The Bakta and antiSMASH images use Docker Hub short names (`oschwengers/bakta`,
`antismash/standalone`). If your own configuration sets `docker.registry` or `apptainer.registry`
(nf-core configs often set `quay.io`), Nextflow prefixes those two names with that registry, where
they do not exist. Leave both settings unset, or override the two processes' `container` with
`docker.io/…` references.
