# Running ZPIC on mybinder.org

To run ZPIC on [mybinder.org](https://mybinder.org/) just click the badge below:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ricardo-fonseca/zpic/HEAD?urlpath=/lab/tree/python/notebooks/README.ipynb)

And choose any of the available notebooks.

If you prefer you can copy the URL below to any browser you like:

`https://mybinder.org/v2/gh/ricardo-fonseca/zpic/HEAD?urlpath=/lab/tree/python/notebooks/README.ipynb`

## More information

The files on this directory set the required information for [mybinder.org](https://mybinder.org/):

* The `requirements.txt` file specifies the required python packages (on top of a standard conda distribution).
* The `apt.txt` file specifies packages to be installed using the `apt` package manager.
* The `postBuild` file instructs [mybinder.org](https://mybinder.org/) to compile the the python modules after creating the container, and to add the required extensions to the Jupyter Lab environment.
