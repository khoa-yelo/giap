# giap-uniop image (GIAP_SPEC.md §4, §6.3; M11-14 item 1b). Pins: ADR-0008.
# This image is built only locally by `giap setup --uniop` from the authors' archive; never pushed
# (unlicensed).
# Upstream https://github.com/hongsua/UniOP at commit 3ed14dfc2a5d69b993d99a26e50c1f8e9723c087 has
# no licence file, so it is never vendored (R5) and containers.yml never builds this image. Setup
# fetches the archive, checks its pinned sha256 and stages it as uniop-src.tar.gz in the build
# context (with this file, uniop.def and workflow/envs/uniop.yml), then runs:
#   docker build -f <context>/Dockerfile -t localhost/giap-uniop:<version>-<key> <context>
FROM mambaorg/micromamba@sha256:f90ebbd10a6dabe1ab87497babd33085d6cf743da6176cf6bf119b9431cdf251
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/uniop.yml /tmp/env.yml
RUN micromamba install -y -n base -c conda-forge -c bioconda --override-channels -f /tmp/env.yml \
    && micromamba clean -a -y
ARG GIT_SHA=3ed14dfc2a5d69b993d99a26e50c1f8e9723c087
COPY --chown=$MAMBA_USER:$MAMBA_USER uniop-src.tar.gz /tmp/src.tar.gz
RUN mkdir -p /opt/conda/share/uniop \
    && tar -xzf /tmp/src.tar.gz --strip-components=1 -C /opt/conda/share/uniop \
    && echo "$GIT_SHA" > /opt/conda/share/uniop/GIAP_COMMIT && rm /tmp/src.tar.gz
WORKDIR /opt/conda/share/uniop
LABEL org.opencontainers.image.title="giap-uniop"
