# gippa-terminatornet image (GIPPA_SPEC.md §4, §6.3). Built by containers.yml from the repository root:
#   docker build -f workflow/containers/terminatornet/Dockerfile -t ghcr.io/<owner>/gippa-terminatornet .
# and referenced by digest in workflow/conf/containers.config. Pins: ADR-0008.
# Upstream https://github.com/btjaden/TerminatorNet at commit 72180309545988dba63b8bd5db5d0b09bf4dba70 (MIT); fetched at build time, never vendored (R5).
FROM mambaorg/micromamba@sha256:f90ebbd10a6dabe1ab87497babd33085d6cf743da6176cf6bf119b9431cdf251
# Installed from the committed locks (M11-14 item 2): the explicit conda lock, then the pip lock
# with --no-deps, so nothing is solved at build time. pip lines `name @ file://…` are conda's own
# packages (already installed by the explicit lock) and are filtered out.
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/locks/terminatornet.explicit.txt /tmp/explicit.txt
RUN micromamba install -y -n base --override-channels -f /tmp/explicit.txt \
    && micromamba clean -a -y && rm /tmp/explicit.txt
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/locks/terminatornet.pip.txt /tmp/pip.txt
RUN { grep -v ' @ file://' /tmp/pip.txt || true; } > /tmp/pip-pypi.txt \
    && micromamba run -n base python -m pip install --no-cache-dir --no-deps -r /tmp/pip-pypi.txt \
    && rm /tmp/pip.txt /tmp/pip-pypi.txt
# sha256 of the GitHub archive of that commit (computed on FarmShare from source_repo; R6)
ARG ARCHIVE_SHA256=b3716eb26ca2a5b1835383ac85a30448a058f1249ac5095c69e85d1ea9eeb33c
ARG GIT_SHA=72180309545988dba63b8bd5db5d0b09bf4dba70
ADD --checksum=sha256:${ARCHIVE_SHA256} --chown=$MAMBA_USER:$MAMBA_USER https://github.com/btjaden/TerminatorNet/archive/${GIT_SHA}.tar.gz /tmp/src.tar.gz
RUN mkdir -p /opt/conda/share/terminatornet \
    && tar -xzf /tmp/src.tar.gz --strip-components=1 -C /opt/conda/share/terminatornet \
    && echo "$GIT_SHA" > /opt/conda/share/terminatornet/GIPPA_COMMIT && rm /tmp/src.tar.gz
WORKDIR /opt/conda/share/terminatornet
# apptainer exec does not run the entrypoint that activates the base env.
ENV PATH="/opt/conda/bin:$PATH" CONDA_PREFIX=/opt/conda
LABEL org.opencontainers.image.title="gippa-terminatornet"
