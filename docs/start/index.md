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

The quickest way to start using ZPIC and learn about PIC simulations is to use `mybinder.org`, just click the button below:

[<i class="fas fa-cloud" aria-hidden="true"></i>  Launch ZPIC binder](https://mybinder.org/v2/gh/ricardo-fonseca/zpic/HEAD?urlpath=/lab/tree/python/notebooks/README.ipynb){: .btn .btn--info}

(It takes a few minutes to load)

## Using ZPIC

ZPIC includes a suite of 1D/2D fully relativistic electromagnetic PIC codes, as well as 1D electrostatic. These codes can all be run from a Python environment, and in particular inside Jupyter notebooks, which is the recommended way of using ZPIC. To use ZPIC with Python three options are available:

1. Launch ZPIC in [mybinder.org](https://mybinder.org/v2/gh/ricardo-fonseca/zpic/HEAD?urlpath=/lab/tree/python/notebooks/README.ipynb). This is the easiest way to use ZPIC but you will be limited to short sessions. (By the way, mybinder.org rules!)

2. Use a [Docker image](docker). This is usually the best option for use in a classroom type environment, as it requires the least amount of configuration (i.e., have Docker installed).

3. [Download](download) ZPIC and [compile](compile) it on your computer. This option is recommended for users wanting to create their own simulations and notebooks, and requires a working C99 compiler and python environment.

For all options, you should check that your ZPIC installation is ok by running the tutorial `ZPIC.ipynb` notebook, available in the `python/notebooks/tutorial/` folder.

Once you have ZPIC up and running the next steps are to read the documentation, available [here](../documentation/python), and to to check the examples library, in the `python/notebooks` folder. The [examples](../examples) page has a detailed description of the notebooks available.

## ZPIC as a C code

The core components of ZPIC are written in C (C99). This is not only for performance reasons but also to more easily allow the use of ZPIC codes in other projects and libraries. As a result, ZPIC codes can also be run as pure C codes, whith no dependencies, making them extremely portable.

For most users this will not be the best choice, as the codes cannot be run interactively this way, and it will be necessary to save simulation output to disk for visualization and post-processing. However, if users want to avoid Python and/or want to run very long simulations, this may be a good choice.

To run the ZPIC codes as pure C codes begin by [downloading](download) the code to your computer, and then follow the instructions in the [Using ZPIC with C](../documentation/c) page.
