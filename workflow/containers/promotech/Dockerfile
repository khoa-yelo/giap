# gippa-promotech image (GIPPA_SPEC.md §4, §6.3). Built by containers.yml from the repository root:
#   docker build -f workflow/containers/promotech/Dockerfile -t ghcr.io/<owner>/gippa-promotech .
# and referenced by digest in workflow/conf/containers.config. Pins: ADR-0008.
# Upstream https://github.com/BioinformaticsLabAtMUN/PromoTech at commit 56251ad9b883ef831b4753fc623d5ec970fe65e0 (GPL-3.0); fetched at build time, never vendored (R5).
FROM mambaorg/micromamba@sha256:f90ebbd10a6dabe1ab87497babd33085d6cf743da6176cf6bf119b9431cdf251
# Installed from the committed locks (M11-14 item 2): the explicit conda lock, then the pip lock
# with --no-deps, so nothing is solved at build time. pip lines `name @ file://…` are conda's own
# packages (already installed by the explicit lock) and are filtered out.
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/locks/promotech.explicit.txt /tmp/explicit.txt
RUN micromamba install -y -n base --override-channels -f /tmp/explicit.txt \
    && micromamba clean -a -y && rm /tmp/explicit.txt
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/locks/promotech.pip.txt /tmp/pip.txt
RUN { grep -v ' @ file://' /tmp/pip.txt || true; } > /tmp/pip-pypi.txt \
    && micromamba run -n base python -m pip install --no-cache-dir --no-deps -r /tmp/pip-pypi.txt \
    && rm /tmp/pip.txt /tmp/pip-pypi.txt
# sha256 of the GitHub archive of that commit (computed on FarmShare from source_repo; R6)
ARG ARCHIVE_SHA256=186d12df99ea83d02355ea29cc6edb5e6d4a0c0d87c5241df0455692910e5c9f
ARG GIT_SHA=56251ad9b883ef831b4753fc623d5ec970fe65e0
ADD --checksum=sha256:${ARCHIVE_SHA256} --chown=$MAMBA_USER:$MAMBA_USER https://github.com/BioinformaticsLabAtMUN/PromoTech/archive/${GIT_SHA}.tar.gz /tmp/src.tar.gz
RUN mkdir -p /opt/conda/share/promotech \
    && tar -xzf /tmp/src.tar.gz --strip-components=1 -C /opt/conda/share/promotech \
    && echo "$GIT_SHA" > /opt/conda/share/promotech/GIPPA_COMMIT && rm /tmp/src.tar.gz
WORKDIR /opt/conda/share/promotech
# promotech.py imports PyQt5 at the top, and Qt links libGL (absent from the conda lock).
USER root
RUN apt-get update && apt-get install -y --no-install-recommends libgl1 \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER
# apptainer exec does not run the entrypoint that activates the base env.
ENV PATH="/opt/conda/bin:$PATH" CONDA_PREFIX=/opt/conda
LABEL org.opencontainers.image.title="gippa-promotech"
