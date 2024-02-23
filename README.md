# statistics-resampling-online

Try the statistics-resampling (and statistics) package online in a Notebook with an Octave kernel within Jupyter Lab using Binder by clicking the launch-binder button:

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/acpennlab/statistics-resampling-online/master?labpath=statistics-resampling.ipynb)

Note that the first time you access statistics-resampling online (since the last commit at the GitHub repository) it may take a while to build the docker image, but frequent access to statistics-resampling-online thereon should give load times less than a minute.

The documentation for the statistics-resampling package can be found [here](https://gnu-octave.github.io/statistics-resampling/index.html). If you wish to download and use the package offline, you can find the package source code on [GitHub](https://github.com/gnu-octave/statistics-resampling/).

Data files (.tsv and .csv) can be conveniently modified using the [jupyterlab-spreadsheet-editor](https://jupyterlab-contrib.github.io/jupyterlab-spreadsheet-editor.html) included in the environment. The environment also includes kernels to run workbooks in R and Python (in which additional packages and modules can be installed from CRAN and with pip respectively, or [conda-forge](https://conda-forge.org/packages/)). You could consider loading Jupyter notebooks and data files in your own GitHub repositories using this Binder environment by using the [nbgitpuller link generator](https://nbgitpuller.readthedocs.io/en/latest/link.html?tab=binder). This Binder repository is `nbgitpuller` enabled, so you can set it use it's environment (and image if it already) and use content (e.g. notebooks and data), and install further dependencies (in an environment.yaml or apt.txt file), in your own GitHub repository, potentially speeding up load times.
