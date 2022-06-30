---
title: Compiling ZPIC for Linux
permalink: /start/compile/linux/

layout: single
toc: true
toc_label: Compiling for Linux

sidebar:
  nav: "docs"
---

## Pre-requisites - All

The examples given here are for the Ubuntu distribution and the `apt` package manager. Other Linux distros (e.g. CentOS, Alpine) are known to work with similar instructions.

To compile the ZPIC C codes on Linux all you need are the basic developer packages. On Ubuntu you can install them using:

```bash
$ sudo apt install build-essential
```

## Pre-requisites - Python

Compiling and using the Python module also requires the [Cython](https://cython.org/) and [Numpy](https://numpy.org/) Python packages. Ubuntu comes with a recent enough `python3` installation so 
just install these packages using whatever package manager you prefer. One good choice is using `pip`:

```bash
$ sudo apt install python3-pip
```

Numpy and Cython can now be installed using:

```bash
$ sudo pip3 install Cython Numpy
```

The modules can be used from any Python script, but using them in [Jupyter](https://jupyter.org) notebooks is recommended. Additional recommended python packages are:

* ipython
* scipy
* matplotlib
* jupyterlab
* Cython
* ipympl

You can install all of the recommended packages doing:

```bash
$ sudo pip3 install Jinja2 ipython scipy matplotlib jupyterlab Cython ipympl
```

Installing the `ffmpeg` library is also useful for generating animations:

```bash
$ sudo apt install ffmpeg
```

### Optional Jupyter extensions

If you want to use (`ipywidgets`) HTML widgets inside your notebooks you will also need to install a recent version of Node.js (≥ 12.0). With the official Ubuntu 20.04 package repositories you'll get Node.js 10.x, so we need to do a little more work:

```bash
$ curl -sL https://deb.nodesource.com/setup_12.x | sudo -E bash -
$ sudo apt-get install -y nodejs
```

This will update the package repositories to look for Node.js 12.x and install it. You should now be able to install the required additional Jupyter extensions:

```bash
$ sudo jupyter labextension install @jupyter-widgets/jupyterlab-manager jupyter-matplotlib
```

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

```bash
$ cd source/zpic/python
$ make
cd source && make
make[1]: Entering directory '/home/zamb/zpic/python/source'
python3 setup.py build_ext --build-lib=../lib
Compiling em1d.pyx because it changed.
Compiling em2d.pyx because it changed.
(...)
x86_64-linux-gnu-gcc -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro -g -fwrapv -O2 -Wl,-Bsymbolic-functions -Wl,-z,relro -g -fwrapv -O2 -g -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 build/temp.linux-x86_64-3.8/em2ds.o build/temp.linux-x86_64-3.8/home/zamb/zpic/em2ds/charge.o build/temp.linux-x86_64-3.8/home/zamb/zpic/em2ds/current.o build/temp.linux-x86_64-3.8/home/zamb/zpic/em2ds/emf.o build/temp.linux-x86_64-3.8/home/zamb/zpic/em2ds/fft.o build/temp.linux-x86_64-3.8/home/zamb/zpic/em2ds/filter.o build/temp.linux-x86_64-3.8/home/zamb/zpic/em2ds/grid2d.o build/temp.linux-x86_64-3.8/home/zamb/zpic/em2ds/particles.o build/temp.linux-x86_64-3.8/home/zamb/zpic/em2ds/random.o build/temp.linux-x86_64-3.8/home/zamb/zpic/em2ds/simulation.o build/temp.linux-x86_64-3.8/home/zamb/zpic/em2ds/timer.o build/temp.linux-x86_64-3.8/home/zamb/zpic/em2ds/zdf.o -o ../lib/em2ds.cpython-38-x86_64-linux-gnu.so
make[1]: Leaving directory '/home/zamb/zpic/python/source'
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

## CentOS 9

To use the code on CentOS 9 you follow the same instructions but use the `dnf` installer instead, and install a couple of extra packages (e.g. `python3-devel`)

```bash
$ sudo dnf group install "Development Tools"
$ sudo dnf install python3-pip python3-devel nodejs
```

Installing `ffmpeg` is a little trickier because it is not included in the standard repos:

```bash
$ sudo dnf install epel-release
$ sudo dnf install https://mirrors.rpmfusion.org/free/el/rpmfusion-free-release-9.noarch.rpm
$ sudo dnf install ffmpeg ffmpeg-devel
```

The python modules should now be installed with `pip3`:

```bash
$ sudo pip3 install Jinja2 ipython scipy matplotlib jupyterlab Cython ipympl
$ sudo /usr/local/bin/jupyter labextension install @jupyter-widgets/jupyterlab-manager jupyter-matplotlib
```

And you should now be able to compile the code and launch the notebooks as described above.
