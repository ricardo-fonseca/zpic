# Docker image to use ZPIC inside jupyter (Python) notebooks
# This version downloads and compiles the latest ZPIC version
#
# Create image with:
#
# > docker build --force-rm --tag zamb/zpic .
#
# --force-rm        : Always remove intermediate containers
# --tag name        : Name for image (in this case zamb/zpic)
#
# Run image with:
#
# > docker run -p 8888:8888 -t --rm zamb/zpic
#
# -p host:container : Publish a container's port(s) to the host. In this case
#                     port 8888 is published as 8888
# -t                : Allocate a pseudo-TTY
# --rm              : Automatically remove the container when it exits
#
# Then just open a web browser on https://locahost:8888
#
# Upload image to Docker hub
# > docker push zamb/zpic

# our base image
FROM jupyter/base-notebook

USER root

# Install additional required system packages
RUN apt-get update && apt-get install -yq --no-install-recommends \
    ffmpeg \
    && rm -rf /var/lib/apt/lists/*


# Install Python 3 packages
# Remove pyqt and qt pulled in for matplotlib since we're only ever going to
# use notebook-friendly backends in these images
RUN conda install --quiet --yes \
    'ipywidgets=7.6*' \
    'ipympl' \
    'matplotlib=3.5*' \
    'cython=0.29*' \
    'scipy' \
    && \
    conda remove --quiet --yes --force qt pyqt && \
    conda clean -tipsy && \
    # Enable passwordless login to Jupyter
    echo "c.NotebookApp.token = ''" >> /home/$NB_USER/.jupyter/jupyter_notebook_config.py && \
    # Activate ipywidgets extension in the environment that runs the notebook server
    jupyter nbextension enable --py widgetsnbextension --sys-prefix && \
    jupyter labextension install @jupyter-widgets/jupyterlab-manager jupyter-matplotlib && \
    jupyter labextension enable @jupyter-widgets/jupyterlab-manager jupyter-matplotlib && \
    npm cache clean --force && \
    rm -rf $CONDA_DIR/share/jupyter/lab/staging && \
    rm -rf /home/$NB_USER/.cache/yarn && \
    rm -rf /home/$NB_USER/.node-gyp && \
    fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER

# Import matplotlib the first time to build the font cache.
ENV XDG_CACHE_HOME /home/$NB_USER/.cache/
RUN MPLBACKEND=Agg python -c "import matplotlib.pyplot" && \
    fix-permissions /home/$NB_USER

# Get ZPIC python modules and install them
# The required build utilities are installed and uninstalled in this step
RUN apt-get update && apt-get install -yq --no-install-recommends \
    gcc \
    libc6-dev \
    make \
    git && \
    rm -rf /var/lib/apt/lists/* && \
    cd /tmp && \
    git clone https://github.com/zambzamb/zpic.git && \
    cd zpic/python && \
    make CC='gcc' && \
    mkdir -p /opt/zpic/lib && \
    chmod a+x /opt/zpic /opt/zpic/lib && \
    cp -v lib/*.so lib/zdf.py /opt/zpic/lib && \
    cp -rv notebooks/* /home/$NB_USER/ && \
    chown $NB_UID /home/$NB_USER/* && \
    rm -rf /tmp/zpic && \
    apt-get --purge remove -yq \
    gcc \
    libc6-dev \
    make \
    git && \
    rm -rf /var/lib/apt/lists/*

# Add additional iPython configuration
COPY ipython_config.py /home/$NB_USER/.ipython/profile_default/
RUN fix-permissions /home/$NB_USER/.ipython

COPY plugin.jupyterlab-settings /home/$NB_USER/.jupyter/lab/user-settings/@jupyterlab/docmanager-extension/
RUN fix-permissions /home/$NB_USER/.jupyter

# Switch back to jovyan to avoid accidental container runs as root
USER $NB_UID

# Add ZPIC modules to Python searchpath
ENV PYTHONPATH /opt/zpic/lib


