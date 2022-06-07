---
title: Compiling ZPIC for Mac OS X
permalink: /start/compile/macosx/

layout: single
toc: true
toc_label: Compiling for Mac OS X

sidebar:
  nav: "docs"
---

## Pre-requisites - All

To compile the ZPIC C codes on Mac OS X all you need are the Xcode Command Line Tools. You don't need the full [Xcode](https://developer.apple.com/xcode/) package, but if you have it installed just make sure you have it up to date.

To install the Xcode Command Line tools just open the Terminal and type the following command:

```bash
$ xcode-select --install
```

If you have the tools already installed you should see something like this:

```bash
$ xcode-select --install
xcode-select: error: command line tools are already
installed, use "Software Update" to install updates
```

Otherwise a dialog box should pop up. Just click "Install" to download and install the Xcode Command Line Tools.

## Pre-requisites - Python

Compiling and using the Python modules also requires the [Cython](https://cython.org/) and [Numpy](https://numpy.org/) Python packages, just install them using whatever package manager you prefer (e.g. `pip3 install Cython`). The modules can be used from any Python script, but using them in [Jupyter](https://jupyter.org) notebooks is recommended.

Additional recommended python packages are:

* ipython
* scipy
* matplotlib
* jupyterlab
* Cython
* ipympl

## Compiling the code

### Individual C code versions

Each of the `C` code versions, available in different directories, is a complete standalone version, and can be compiled separately. Each of these directories contains a local `Makefile` that you can edit if you need to change any of the compiler parameters (or want to specify a different compiler).

To compile any of the code versions, just open a terminal window and navigate to the directory with the code that you want to use. You can then compile the code using the `make` command. You can clean up all the generated object files and executable by issuing the `make clean` command.

In this example we will compile the `em2d` code. We assume you have the code on inside a directory named `source` in your home folder. Just go into that directory and run `make`:

```bash
$ cd source/zpic/em2d
$ make
gcc -c -Ofast -std=c99 -pedantic current.c -o current.o
gcc -c -Ofast -std=c99 -pedantic emf.c -o emf.o
(...)
gcc -Ofast -std=c99 -pedantic  current.o emf.o particles.o random.o timer.o main.o simulation.o zdf.o -o zpic
$ 
```

Once compilation has finished there will be a `zpic` binary file in the same directory that you can run.

### Python modules

To compile the python modules navigate to the `python` directory on the ZPIC sourcetree and run `make`:

```shell
$ cd source/zpic/python
$ make
cd source && /Library/Developer/CommandLineTools/usr/bin/make
python3 setup.py build_ext --build-lib=../lib
Compiling em1d.pyx because it changed.
Compiling em2d.pyx because it changed.
(...)
clang -bundle -undefined dynamic_lookup -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX12.sdk build/temp.macosx-12-x86_64-3.9/Users/zamb/Source/zpic/em2ds/charge.o build/temp.macosx-12-x86_64-3.9/Users/zamb/Source/zpic/em2ds/current.o build/temp.macosx-12-x86_64-3.9/Users/zamb/Source/zpic/em2ds/emf.o build/temp.macosx-12-x86_64-3.9/Users/zamb/Source/zpic/em2ds/fft.o build/temp.macosx-12-x86_64-3.9/Users/zamb/Source/zpic/em2ds/filter.o build/temp.macosx-12-x86_64-3.9/Users/zamb/Source/zpic/em2ds/grid2d.o build/temp.macosx-12-x86_64-3.9/Users/zamb/Source/zpic/em2ds/particles.o build/temp.macosx-12-x86_64-3.9/Users/zamb/Source/zpic/em2ds/random.o build/temp.macosx-12-x86_64-3.9/Users/zamb/Source/zpic/em2ds/simulation.o build/temp.macosx-12-x86_64-3.9/Users/zamb/Source/zpic/em2ds/timer.o build/temp.macosx-12-x86_64-3.9/Users/zamb/Source/zpic/em2ds/zdf.o build/temp.macosx-12-x86_64-3.9/em2ds.o -o ../lib/em2ds.cpython-39-darwin.so
$ 
```

Libraries will be created in the `python/lib` directory, you need to add that library to your Python environment so you can use the modules.

## Launching ZPIC notebooks

Once you have compiled the Python modules, you shoud be able to launch the ZPIC notebooks by navigating into the `notebooks` directory and launching Jupyter:

```bash
$ cd python/notebooks
$ jupyter lab
[I 2022-06-07 14:05:36.220 ServerApp] jupyterlab | extension was successfully linked.
(...)
[I 2022-06-07 14:05:36.673 ServerApp] Jupyter Server 1.13.3 is running at:
[I 2022-06-07 14:05:36.673 ServerApp] http://localhost:8888/lab?token=a13a9886d3948335ac4a920cdbd24d9c498d9b9f493c700a
(...)
```

Depending on your configuration, a browser should have opened with your Jupyter session. If not, just copy the URL shown in your terminal (including the token) into the URL bar of any browser.

Be sure to open the `README.ipynb` file first.