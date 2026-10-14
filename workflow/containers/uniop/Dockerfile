# giap-uniop image (GIAP_SPEC.md §4, §6.3). Built by containers.yml from the repository root:
#   docker build -f workflow/containers/uniop/Dockerfile -t ghcr.io/<owner>/giap-uniop .
# and referenced by digest in workflow/conf/containers.config. Pins: ADR-0008.
# Upstream https://github.com/hongsua/UniOP at commit 3ed14dfc2a5d69b993d99a26e50c1f8e9723c087 (unlicensed: built for local use only, never pushed publicly without the authors' permission); fetched at build time, never vendored (R5).
FROM mambaorg/micromamba@sha256:f90ebbd10a6dabe1ab87497babd33085d6cf743da6176cf6bf119b9431cdf251
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/uniop.yml /tmp/env.yml
RUN micromamba install -y -n base -c conda-forge -c bioconda --override-channels -f /tmp/env.yml \
    && micromamba clean -a -y
ARG GIT_SHA=3ed14dfc2a5d69b993d99a26e50c1f8e9723c087
ADD --chown=$MAMBA_USER:$MAMBA_USER https://github.com/hongsua/UniOP/archive/${GIT_SHA}.tar.gz /tmp/src.tar.gz
RUN mkdir -p /opt/conda/share/uniop \
    && tar -xzf /tmp/src.tar.gz --strip-components=1 -C /opt/conda/share/uniop \
    && echo "$GIT_SHA" > /opt/conda/share/uniop/GIAP_COMMIT && rm /tmp/src.tar.gz
WORKDIR /opt/conda/share/uniop
LABEL org.opencontainers.image.title="giap-uniop"
