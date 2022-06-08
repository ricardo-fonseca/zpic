---
title: Compiling ZPIC
permalink: /start/compile/

layout: single
toc: true
toc_label: Compiling ZPIC

sidebar:
  nav: "docs"
---

## Quickstart

__ZPIC__ requires only a C99 compiler (gcc works great) and (GNU) make. If you have these tools available on a Linux/Mac OS X system you can compile the code simply by navigating into one of the source folders and running `make`, e.g.:

```bash
$ cd em2d
$ make
```

Building the Python modules requires [Cython](https://cython.org/) and [Numpy](https://numpy.org/). To build the Python modules just navigate into the `python` source folder and run `make`:

```bash
$ cd python
$ make
```

## Compiling ZPIC

### MacOS X

For Mac OS X systems all you will need to install are the Xcode Command Line Tools and the additional Python packages.

You can find detailed instructions [here](macosx).

### Linux

For Linux (and other *NIX systems) you will need a set of build tools installed and a working Python3 installation. It should be straightforward to reproduce the procedure on any Linux distro.

You can find detailed instructions for Ubuntu [here](linux).

### Windows (Windows Subsystem for Linux)

The Microsoft Visual C++ compiler (MSVC) does not support the required C99 features (in particular the Complex type support which is crucial for the spectral versions of ZPIC) and so it cannot be used. If you plan to use only the C versions of ZPIC then installing `gcc` from [MinGW-w64](https://www.mingw-w64.org) or some other source is sufficient. Running the Python versions, however, is a little more complicated because of incompatibilities between the required compilers and standard Windows python binaries.

The recommended way of compiling ZPIC on Windows systems is using the Windows Subsystem for Linux (WSL). This will install a Linux OS (usually Ubuntu) that is run alongside your normal Windows applications, and that you will use for compiling/running ZPIC.

You can find detailed instructions [here](win64).

### Windows (MSYS2)

If you are unable (or unwilling) to use WSL, then you can alternatively use MSYS2. MSYS2 is a collection of tools and libraries providing you with an easy-to-use environment for building, installing and running native Windows software. It is more lightweight than WSL, requiring fewer resources and working on older systems, but it is also less powerful and not as versatile as WSL.

You can find detailed instructions [here](msys2).
