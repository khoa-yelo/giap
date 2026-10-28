# gippa-step image (M11-15 item 8a, slice 15g; the Q-K39 ADR, ADR-0065): the package and its
# permissive Python dependencies only, so GIPPA_STEP tasks need no Python env on compute nodes.
# No bioinformatics tool, no conda, no R. Built from the repository root after `uv build --wheel`:
#   docker build -f workflow/containers/step/Dockerfile --build-arg WHEEL=<file>.whl \
#     --build-arg WHEEL_SHA256=<sha256> --build-arg VERSION=<version> \
#     --build-arg REVISION=<git sha> -t ghcr.io/<owner>/gippa-step:<version> .
# scripts/dev/build_step_image.sh does this for dev builds; scripts/release/check_image_licences.py
# runs inside the built image (licences.allow).
FROM python:3.12-slim-bookworm@sha256:6dd8b6dab6b3adc7812bab69eca5dbdc1a5cc690ad28306d234ffb4ef8ac9753
ARG WHEEL
ARG WHEEL_SHA256=unknown
ARG VERSION=unknown
ARG REVISION=unknown
# Nextflow's task metrics need ps (procps); the slim base has none. System software, aggregated.
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps ca-certificates \
    && rm -rf /var/lib/apt/lists/*
# The hashed setuptools, then the hashed runtime lock (uv export --frozen --no-dev
# --no-emit-project) without build isolation: pandarallel and upsetplot are sdists, and an
# isolated build would fetch an unhashed setuptools. Then the wheel alone. Bind mounts
# (BuildKit), not COPY, so no file stays in a layer.
RUN --mount=type=bind,source=workflow/containers/step/build-requirements.lock,target=/tmp/build-requirements.lock \
    --mount=type=bind,source=workflow/containers/step/requirements.lock,target=/tmp/requirements.lock \
    pip install --no-cache-dir --require-hashes -r /tmp/build-requirements.lock \
    && pip install --no-cache-dir --require-hashes --no-build-isolation -r /tmp/requirements.lock
# An empty WHEEL would mount dist/ itself; the wheel must match WHEEL_SHA256 (the label).
RUN test -n "${WHEEL}" || { echo "build arg WHEEL (the dist/ wheel file name) is required" >&2; exit 1; } \
    && test "${WHEEL_SHA256}" != unknown || { echo "build arg WHEEL_SHA256 is required" >&2; exit 1; }
RUN --mount=type=bind,source=dist/${WHEEL},target=/tmp/${WHEEL} \
    echo "${WHEEL_SHA256}  /tmp/${WHEEL}" | sha256sum -c - \
    && pip install --no-cache-dir --no-deps /tmp/${WHEEL}
# The marker read by gippa.io.fgsea_r: no fgsea prefix in this image (empty content).
RUN mkdir -p /opt/gippa && touch /opt/gippa/STEP_IMAGE
# Runs as any UID (HOME is writable by everyone).
ENV HOME=/tmp PYTHONDONTWRITEBYTECODE=1
LABEL org.opencontainers.image.title="gippa-step" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.source="https://github.com/khoa-yelo/gippa" \
      org.opencontainers.image.revision="${REVISION}" \
      org.opencontainers.image.licenses="MIT" \
      org.gippa.wheel.sha256="${WHEEL_SHA256}" \
      org.opencontainers.image.description="gippa and its permissive Python dependencies"
