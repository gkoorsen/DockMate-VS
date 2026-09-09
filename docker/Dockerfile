FROM mambaorg/micromamba:2.3.2

COPY --chown=$MAMBA_USER:$MAMBA_USER docker/environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml \
    && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# WORKDIR creates the directory owned by root even though the image runs as
# $MAMBA_USER, and COPY --chown only sets ownership on the entries it copies,
# not on the pre-existing parent. pip builds in-tree, so setuptools needs to
# write dockmate_vs.egg-info into this directory: create it and hand it over
# before switching back to the unprivileged user.
USER root
RUN mkdir -p /opt/dockmate-vs && chown $MAMBA_USER:$MAMBA_USER /opt/dockmate-vs
USER $MAMBA_USER

WORKDIR /opt/dockmate-vs
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/dockmate-vs
RUN python -m pip install --no-cache-dir openmm==8.6.0 pdbfixer==1.12.0 \
    && python -m pip install --no-deps . \
    && rm -rf /opt/dockmate-vs/dockmate_vs.egg-info \
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
