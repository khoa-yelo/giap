# giap-terminatornet image (GIAP_SPEC.md §4, §6.3). Built by containers.yml from the repository root:
#   docker build -f workflow/containers/terminatornet/Dockerfile -t ghcr.io/<owner>/giap-terminatornet .
# and referenced by digest in workflow/conf/containers.config. Pins: ADR-0008.
# Upstream https://github.com/btjaden/TerminatorNet at commit 72180309545988dba63b8bd5db5d0b09bf4dba70 (MIT); fetched at build time, never vendored (R5).
FROM mambaorg/micromamba@sha256:f90ebbd10a6dabe1ab87497babd33085d6cf743da6176cf6bf119b9431cdf251
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/terminatornet.yml /tmp/env.yml
RUN micromamba install -y -n base -c conda-forge -c bioconda --override-channels -f /tmp/env.yml \
    && micromamba clean -a -y
ARG GIT_SHA=72180309545988dba63b8bd5db5d0b09bf4dba70
ADD --chown=$MAMBA_USER:$MAMBA_USER https://github.com/btjaden/TerminatorNet/archive/${GIT_SHA}.tar.gz /tmp/src.tar.gz
RUN mkdir -p /opt/conda/share/terminatornet \
    && tar -xzf /tmp/src.tar.gz --strip-components=1 -C /opt/conda/share/terminatornet \
    && echo "$GIT_SHA" > /opt/conda/share/terminatornet/GIAP_COMMIT && rm /tmp/src.tar.gz
WORKDIR /opt/conda/share/terminatornet
LABEL org.opencontainers.image.title="giap-terminatornet"
