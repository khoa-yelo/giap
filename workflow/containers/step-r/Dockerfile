# gippa-step-r image (M11-15 item 8a, slice 15g; the Q-K39 ADR, ADR-0065): the gippa-step contents
# plus the fgsea image's locked R env, for FUNCTION_STEP and SENSITIVITY (they start R fgsea, and a
# container cannot start a sibling container). R (GPL-2) and fgsea (MIT) are separate programs
# started by subprocess, never imported (R5): an aggregate, as in the gippa-fgsea image (Q-K37).
# Built from the repository root after `uv build --wheel`, with the build args of step/Dockerfile:
#   docker build -f workflow/containers/step-r/Dockerfile --build-arg WHEEL=<file>.whl ... .
# Stage 1: the fgsea explicit lock (M11-14 item 2) installed at /opt/conda; nothing is solved.
FROM mambaorg/micromamba@sha256:f90ebbd10a6dabe1ab87497babd33085d6cf743da6176cf6bf119b9431cdf251 AS fgsea
RUN --mount=type=bind,source=workflow/envs/locks/fgsea.explicit.txt,target=/tmp/explicit.txt \
    micromamba install -y -n base --override-channels -f /tmp/explicit.txt \
    && micromamba clean -a -y
# Stage 2: the same python base and steps as step/Dockerfile, then the R env copied in.
FROM python:3.12-slim-bookworm@sha256:6dd8b6dab6b3adc7812bab69eca5dbdc1a5cc690ad28306d234ffb4ef8ac9753
ARG WHEEL
ARG WHEEL_SHA256=unknown
ARG VERSION=unknown
ARG REVISION=unknown
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps ca-certificates \
    && rm -rf /var/lib/apt/lists/*
# Bind mounts (BuildKit), not COPY, so neither file stays in a layer.
RUN --mount=type=bind,source=workflow/containers/step/requirements.lock,target=/tmp/requirements.lock \
    pip install --no-cache-dir --require-hashes -r /tmp/requirements.lock
RUN --mount=type=bind,source=dist/${WHEEL},target=/tmp/${WHEEL} \
    pip install --no-cache-dir --no-deps /tmp/${WHEEL}
COPY --from=fgsea /opt/conda /opt/conda
# The marker read by gippa.io.fgsea_r: this image's fgsea prefix (it overrides GIPPA_FGSEA_PREFIX).
RUN mkdir -p /opt/gippa && printf '%s\n' /opt/conda > /opt/gippa/STEP_IMAGE
# R is found by its prefix, not PATH: python stays the image's own (/usr/local/bin first).
ENV HOME=/tmp PYTHONDONTWRITEBYTECODE=1
LABEL org.opencontainers.image.title="gippa-step-r" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.source="https://github.com/khoa-yelo/gippa" \
      org.opencontainers.image.revision="${REVISION}" \
      org.opencontainers.image.licenses="MIT AND GPL-2.0-or-later" \
      org.gippa.wheel.sha256="${WHEEL_SHA256}" \
      org.opencontainers.image.description="gippa with R fgsea 1.36.2 (MIT) in R 4.5.3 (GPL-2)"
