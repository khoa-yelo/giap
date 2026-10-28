# gippa-fgsea image (M11-14 item 1, 14c; ADR-0038, ADR-0060 §1). Built by containers.yml from the
# repository root:
#   docker build -f workflow/containers/fgsea/Dockerfile -t ghcr.io/<owner>/gippa-fgsea .
# Its ref lives only in workflow/conf/images.json ("processes": []): no Nextflow process runs it.
# gippa.io.fgsea_r starts it (docker run / apptainer exec) when no GIPPA_FGSEA_PREFIX and no prebuilt
# gippa-fgsea env is found. R (GPL-2) and fgsea (MIT) run only as a separate process (R5).
FROM mambaorg/micromamba@sha256:f90ebbd10a6dabe1ab87497babd33085d6cf743da6176cf6bf119b9431cdf251
# Installed from the committed explicit lock (M11-14 item 2): nothing is solved at build time.
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/locks/fgsea.explicit.txt /tmp/explicit.txt
RUN micromamba install -y -n base --override-channels -f /tmp/explicit.txt \
    && micromamba clean -a -y && rm /tmp/explicit.txt
# The script at the path the image argv names (Rscript /opt/gippa/gippa_fgsea.R <options>).
COPY workflow/bin/gippa_fgsea.R /opt/gippa/gippa_fgsea.R
# apptainer exec does not run the entrypoint that activates the base env: Rscript by name.
ENV PATH="/opt/conda/bin:$PATH" LC_ALL=C
LABEL org.opencontainers.image.title="gippa-fgsea" \
      org.opencontainers.image.licenses="MIT AND GPL-2.0-or-later" \
      org.opencontainers.image.description="R fgsea 1.36.2 (MIT) in R 4.5.3 (GPL-2) with gippa_fgsea.R"
