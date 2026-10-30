# gippa-promoter_calculator image (GIPPA_SPEC.md §4, §6.3). Built by containers.yml from the repository root:
#   docker build -f workflow/containers/promoter_calculator/Dockerfile -t ghcr.io/<owner>/gippa-promoter_calculator .
# and referenced by digest in workflow/conf/containers.config. Pins: ADR-0008.
# Upstream promotercalculator 1.2.4 (GPL-3.0-or-later) from PyPI via the pip lock.
FROM mambaorg/micromamba@sha256:f90ebbd10a6dabe1ab87497babd33085d6cf743da6176cf6bf119b9431cdf251
# Installed from the committed locks (M11-14 item 2): the explicit conda lock, then the pip lock
# with --no-deps, so nothing is solved at build time. pip lines `name @ file://…` are conda's own
# packages (already installed by the explicit lock) and are filtered out.
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/locks/promoter_calculator.explicit.txt /tmp/explicit.txt
RUN micromamba install -y -n base --override-channels -f /tmp/explicit.txt \
    && micromamba clean -a -y && rm /tmp/explicit.txt
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/locks/promoter_calculator.pip.txt /tmp/pip.txt
RUN { grep -v ' @ file://' /tmp/pip.txt || true; } > /tmp/pip-pypi.txt \
    && micromamba run -n base python -m pip install --no-cache-dir --no-deps -r /tmp/pip-pypi.txt \
    && rm /tmp/pip.txt /tmp/pip-pypi.txt
# apptainer exec does not run the entrypoint that activates the base env.
ENV PATH="/opt/conda/bin:$PATH" CONDA_PREFIX=/opt/conda
LABEL org.opencontainers.image.title="gippa-promoter_calculator"
