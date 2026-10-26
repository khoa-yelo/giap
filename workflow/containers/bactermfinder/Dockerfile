# giap-bactermfinder image (GIAP_SPEC.md §4, §6.3). Built by containers.yml from the repository root:
#   docker build -f workflow/containers/bactermfinder/Dockerfile -t ghcr.io/<owner>/giap-bactermfinder .
# and referenced by digest in workflow/conf/containers.config. Pins: ADR-0008.
# Upstream https://github.com/BioinformaticsLabAtMUN/BacTermFinder at commit 024870447f6288ebc9d29b314e957a6dfdfeb5f6 (GPL-3.0); fetched at build time, never vendored (R5).
FROM mambaorg/micromamba@sha256:f90ebbd10a6dabe1ab87497babd33085d6cf743da6176cf6bf119b9431cdf251
# Installed from the committed locks (M11-14 item 2): the explicit conda lock, then the pip lock
# with --no-deps, so nothing is solved at build time. pip lines `name @ file://…` are conda's own
# packages (already installed by the explicit lock) and are filtered out.
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/locks/bactermfinder.explicit.txt /tmp/explicit.txt
RUN micromamba install -y -n base --override-channels -f /tmp/explicit.txt \
    && micromamba clean -a -y && rm /tmp/explicit.txt
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/locks/bactermfinder.pip.txt /tmp/pip.txt
RUN { grep -v ' @ file://' /tmp/pip.txt || true; } > /tmp/pip-pypi.txt \
    && micromamba run -n base python -m pip install --no-cache-dir --no-deps -r /tmp/pip-pypi.txt \
    && rm /tmp/pip.txt /tmp/pip-pypi.txt
# sha256 of the GitHub archive of that commit (computed on FarmShare from source_repo; R6)
ARG ARCHIVE_SHA256=a6951130bad0fbb04e4d6d912f177c111f95fadda65cdb8ba218af813936937a
ARG GIT_SHA=024870447f6288ebc9d29b314e957a6dfdfeb5f6
ADD --checksum=sha256:${ARCHIVE_SHA256} --chown=$MAMBA_USER:$MAMBA_USER https://github.com/BioinformaticsLabAtMUN/BacTermFinder/archive/${GIT_SHA}.tar.gz /tmp/src.tar.gz
RUN mkdir -p /opt/conda/share/bactermfinder \
    && tar -xzf /tmp/src.tar.gz --strip-components=1 -C /opt/conda/share/bactermfinder \
    && echo "$GIT_SHA" > /opt/conda/share/bactermfinder/GIAP_COMMIT && rm /tmp/src.tar.gz
WORKDIR /opt/conda/share/bactermfinder
LABEL org.opencontainers.image.title="giap-bactermfinder"
