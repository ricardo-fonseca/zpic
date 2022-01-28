---
title: Getting started

#usemathjax: true
permalink: /start/

layout: single
toc: true
toc_label: Getting started

sidebar:
  nav: "docs"
---

## Quickstart

The quickest way to start using ZPIC and learn about PIC simulations is to use a [Docker](https://www.docker.com/) Image. Check the [Docker image](/start/docker) page for instruction on how to get the image running, and start with the tutorial notebook available in `notebooks/tutorial/ZPIC`.

## Using ZPIC

ZPIC includes a suite of 1D/2D fully relativistic electromagnetic PIC codes, as well as 1D electrostatic. These codes can all be run from a Python environment, and in particular inside Jupyter notebooks, which is the recommended way of using ZPIC. To use ZPIC with Python two options are available:

1. [Download](/start/download) ZPIC and [compile](/start/compile) it on your computer. This option is recommended for users wanting to create their own simulations and notebooks, and requires a working C99 compiler and python environment.
2. Use a [Docker image](/start/docker). This is usually the best option for use in a classroom type environment, as it requires the least amount of configuration (i.e., have Docker installed).

For both options, you should check that your ZPIC installation by running the tutorial `ZPIC.ipynb` notebook, available in the `python/notebooks/tutorial/` folder.

Once you have ZPIC up and running the next steps are to read the documentation, available [here](/documentation/python), and to to check the examples library, in the `python/notebooks` folder. The [examples](/examples) page has a detailed description of the notebooks available.

## ZPIC as a C code

The core components of ZPIC are written in C (C99). This is not only for performance reasons but also to more easily allow the use of ZPIC codes in other projects and libraries. As a result, ZPIC codes can also be run as pure C codes, whith no dependencies, making them extremely portable.

For most users this will not be the best choice, as the codes cannot be run interactively this way, and it will be necessary to save simulation output to disk for visualization and post-processing. However, if users want to avoid Python and/or want to run very long simulations, this may be a good choice.

To run the ZPIC codes as pure C codes begin by [downloading](/start/download) the code to your computer, and then follow the instructions in the [Using ZPIC with C](/documentation/c) page.
