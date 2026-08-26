FROM mambaorg/micromamba:2.3.2

COPY --chown=$MAMBA_USER:$MAMBA_USER docker/environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml \
    && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
WORKDIR /opt/dockmate-vs
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/dockmate-vs
RUN python -m pip install --no-cache-dir openmm==8.6.0 pdbfixer==1.12.0 \
    && python -m pip install --no-deps . \
    && test -x /opt/conda/bin/vina \
    && test -x /opt/conda/bin/smina \
    && test -x /opt/conda/bin/rbdock \
    && test -x /opt/conda/bin/obabel

# DockMate-VS treats the environment prefix as the rDock installation root.
# Bioconda may place rDock data below share/, so expose it at the conventional
# <root>/data location expected by rDock and the existing workflow.
USER root
RUN if [ ! -e /opt/conda/data/scripts/dock.prm ]; then \
      dock_prm="$(find /opt/conda/share -path '*/data/scripts/dock.prm' -print -quit)"; \
      test -n "$dock_prm"; \
      data_dir="${dock_prm%/scripts/dock.prm}"; \
      ln -s "$data_dir" /opt/conda/data; \
    fi
USER $MAMBA_USER

ENV PATH=/opt/conda/bin:$PATH \
    RBT_ROOT=/opt/conda \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg

WORKDIR /workspace
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "dockmate-vs"]
CMD ["doctor", "--strict"]
