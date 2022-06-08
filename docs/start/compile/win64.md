---
title: Compiling ZPIC for Windows (MSL)
permalink: /start/compile/win64/

layout: single
toc: true
toc_label: Compiling for Windows (MSL)

sidebar:
  nav: "docs"
---

## ZPIC on Windows using WSL

The [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about) (WSL) lets developers run a GNU/Linux environment directly on Windows, unmodified, without the overhead of a traditional virtual machine or dual-boot setup. In this guide we will be using the Ubuntu Linux distribution that is installed with WSL. Other distributions are also known to work.

The instructions provided here were tested with Microsoft __Windows 10 Pro__, Version __20H2__ (OS Build 19042.1706) on June 2022.

## Installing the Windows Subsystem for Linux + Ubuntu

To install WSL + Ubuntu you should just follow the instructions provided by Microsoft here:

[https://docs.microsoft.com/en-us/windows/wsl/install](https://docs.microsoft.com/en-us/windows/wsl/install)

We also recommend that you install the [Windows Terminal](https://apps.microsoft.com/store/detail/windows-terminal/9N0DX20HK701) application from the Microsoft Store and use it as your terminal application.

### Initial Configuration

Once you have installed WSL + Ubuntu and the Windows Terminal you can simply open the latter and choose "Ubuntu". If you have not done so already begin by updating all your software packages:

```bash
$ sudo apt-get update
$ sudo apt-get upgrade
```

After updating all your software packages python should be (as of June 2022) up to version 3.8.10, which is fine. We just need some additional dependencies:

```bash
$ sudo apt install build-essential ffmpeg python3-pip
```

Now we need to add some additional python packages. We will be doing it using `pip3` (that we just installed):

```bash
$ sudo pip3 install Jinja2 ipython scipy matplotlib jupyterlab Cython ipympl
```

_Note_: Although it is possible to run the above `pip3` command without `sudo` this would install the packages in a local user directory not in the path.

At this point you should be ready to launch Jupyter:

```bash
$ jupyter lab
[I 2022-06-07 11:38:12.245 ServerApp] jupyterlab | extension was successfully linked.
[I 2022-06-07 11:38:12.255 ServerApp] nbclassic | extension was successfully linked.
(...)
```

If you get an error like this:

```text
(...)
  File "/usr/lib/python3/dist-packages/jinja2/utils.py", line 656, in <module>
    from markupsafe import Markup, escape, soft_unicode
ImportError: cannot import name 'soft_unicode' from 'markupsafe' (/usr/local/lib/python3.8/dist-packages/markupsafe/__init__.py)
```

You need to update your Jinja2 library like this:

```bash
$ sudo pip3 install -U Jinja2
```

And things should be fine.

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

### Cloning (Downloading) ZPIC

The Ubuntu distribution comes with `git` pre-installed so if you have not done so already, you can clone ZPIC from GitHub by doing:

```bash
$ mkdir zpic
$ cd zpic
$ git clone https://github.com/ricardo-fonseca/zpic.git .
```

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
$ cd source/zpic/python/notebooks
$ jupyter lab
[I 2022-06-07 14:05:36.220 ServerApp] jupyterlab | extension was successfully linked.
(...)
[I 2022-06-07 14:05:36.673 ServerApp] Jupyter Server 1.13.3 is running at:
[I 2022-06-07 14:05:36.673 ServerApp] http://localhost:8888/lab?token=a13a9886d3948335ac4a920cdbd24d9c498d9b9f493c700a
(...)
```

You will probably need to copy the URL shown in your terminal (including the token) into the URL bar of any browser.

Be sure to open the `README.ipynb` file first.
