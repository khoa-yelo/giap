# giap-promotech image (GIAP_SPEC.md §4, §6.3). Built by containers.yml from the repository root:
#   docker build -f workflow/containers/promotech/Dockerfile -t ghcr.io/<owner>/giap-promotech .
# and referenced by digest in workflow/conf/containers.config. Pins: ADR-0008.
# Upstream https://github.com/BioinformaticsLabAtMUN/PromoTech at commit 56251ad9b883ef831b4753fc623d5ec970fe65e0 (GPL-3.0); fetched at build time, never vendored (R5).
FROM mambaorg/micromamba@sha256:f90ebbd10a6dabe1ab87497babd33085d6cf743da6176cf6bf119b9431cdf251
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/promotech.yml /tmp/env.yml
RUN micromamba install -y -n base -c conda-forge -c bioconda --override-channels -f /tmp/env.yml \
    && micromamba clean -a -y
ARG GIT_SHA=56251ad9b883ef831b4753fc623d5ec970fe65e0
ADD --chown=$MAMBA_USER:$MAMBA_USER https://github.com/BioinformaticsLabAtMUN/PromoTech/archive/${GIT_SHA}.tar.gz /tmp/src.tar.gz
RUN mkdir -p /opt/conda/share/promotech \
    && tar -xzf /tmp/src.tar.gz --strip-components=1 -C /opt/conda/share/promotech \
    && echo "$GIT_SHA" > /opt/conda/share/promotech/GIAP_COMMIT && rm /tmp/src.tar.gz
WORKDIR /opt/conda/share/promotech
LABEL org.opencontainers.image.title="giap-promotech"
