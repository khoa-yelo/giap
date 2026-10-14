# giap-ostir image (GIAP_SPEC.md §4, §6.3). Built by containers.yml from the repository root:
#   docker build -f workflow/containers/ostir/Dockerfile -t ghcr.io/<owner>/giap-ostir .
# and referenced by digest in workflow/conf/containers.config. Pins: ADR-0008.
# Upstream OSTIR 1.1.3 (GPL-3.0) from PyPI via the env file.
FROM mambaorg/micromamba@sha256:f90ebbd10a6dabe1ab87497babd33085d6cf743da6176cf6bf119b9431cdf251
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/ostir.yml /tmp/env.yml
RUN micromamba install -y -n base -c conda-forge -c bioconda --override-channels -f /tmp/env.yml \
    && micromamba clean -a -y
LABEL org.opencontainers.image.title="giap-ostir"
