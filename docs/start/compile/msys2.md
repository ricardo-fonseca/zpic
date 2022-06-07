---
title: Compiling ZPIC for Windows (MSYS2)
permalink: /start/compile/msys2

layout: single
toc: true
toc_label: Compiling for Windows (MSYS2)

sidebar:
  nav: "docs"
---

## ZPIC on Windows using MSYS2

For using ZPIC on Windows we recommend using the Windows Subsystem for Linux (WSL) instead. You can find instructions for running ZPIC on Windows using WSL [here](../win64).

If you are unable (or unwilling) to use WSL you can try MSYS2 instead. MSYS2 is a collection of tools and libraries providing you with an easy-to-use environment for building, installing and running native Windows software. For more information on MSYS2 please check their website:

[https://www.msys2.org/](https://www.msys2.org/)

The instructions provided here allow you to set up a working MSYS2 environment on windows for using the ZPIC code notebooks. This was tested with:

* Microsoft Windows 10 Pro, Version 10.0.19042 Build 19042
* MSYS2 version msys2-x86_64-20220319

On April 2022.

## Installing MSYS2

1. Download installer from [https://www.msys2.org](https://www.msys2.org/)
2. Double-click and Install

### Initial Configuration

(Follows the instructions from the MSYS2 website)

Open MSYS2 and run the following command to update the package database and base packages:

```bash
$ pacman -Syu
```

This will close the MSYS2 window once it is complete. Now open the `MSYS2 MinGW x64` application and run the following command to update the rest of the packages:

```bash
$ pacman -Su
```

### Additional modules

To install the additional packages required by ZPIC notebooks run the following commands:

```bash
$ pacman -S mingw-w64-x86_64-python-ipython
$ pacman -S mingw-w64-x86_64-python-scipy
$ pacman -S mingw-w64-x86_64-python-matplotlib
$ pacman -S mingw-w64-x86_64-python-jupyter_console
$ pacman -S mingw-w64-x86_64-python-argon2_cffi
$ pacman -S mingw-w64-x86_64-python-jupyterlab-pygments
$ pacman -S git
$ pacman -S mingw-w64-x86_64-ffmpeg
```

### Compiler toolchain

You now need to add a compiler toolchaing for building the ZPIC C/Python modules:

```bash
$ pacman -S --needed base-devel mingw-w64-x86_64-toolchain
```

### Additional Python modules

Finally you will need a few additional modules that must be installed from within python:

```bash
$ pacman -S mingw-w64-x86_64-python-pip
$ python -m pip install jupyterlab
$ python -m pip install ipympl
$ python -m pip install cython 
```

## Build ZPIC Python modules

From now on your MSYS2 environment will work like a Linux/Mac OS X environment. Just navigate to your `zpic/python` directory and type `mingw32-make`:

```bash
$ cd zpic/python
$ mingw32-make
$ cd source && C:/msys64/mingw64/bin/mingw32-make.exe
mingw32-make[1]: Entering directory 'C:/msys64/home/zpic_user/zpic/python/source'
python3 setup.py build_ext --build-lib=../lib
Compiling em1d.pyx because it changed.
Compiling em2d.pyx because it changed.
Compiling es1d.pyx because it changed.
Compiling em1ds.pyx because it changed.
Compiling em2ds.pyx because it changed.
[1/5] Cythonizing em1d.pyx
[2/5] Cythonizing em1ds.pyx

(...)

writing build/temp.mingw_x86_64-3.9/msys64/home/zpic_user/zpic/em2ds/em2ds.cp39-mingw_x86_64.def
gcc -shared -Wl,--enable-auto-image-base -pipe -pipe -s build/temp.mingw_x86_64-3.9/msys64/home/zpic_user/zpic/em2ds/charge.o build/temp.mingw_x86_64-3.9/msys64/home/zpic_user/zpic/em2ds/current.o build/temp.mingw_x86_64-3.9/msys64/home/zpic_user/zpic/em2ds/emf.o build/temp.mingw_x86_64-3.9/msys64/home/zpic_user/zpic/em2ds/fft.o build/temp.mingw_x86_64-3.9/msys64/home/zpic_user/zpic/em2ds/filter.o build/temp.mingw_x86_64-3.9/msys64/home/zpic_user/zpic/em2ds/grid2d.o build/temp.mingw_x86_64-3.9/msys64/home/zpic_user/zpic/em2ds/particles.o build/temp.mingw_x86_64-3.9/msys64/home/zpic_user/zpic/em2ds/random.o build/temp.mingw_x86_64-3.9/msys64/home/zpic_user/zpic/em2ds/simulation.o build/temp.mingw_x86_64-3.9/msys64/home/zpic_user/zpic/em2ds/timer.o build/temp.mingw_x86_64-3.9/msys64/home/zpic_user/zpic/em2ds/zdf.o build/temp.mingw_x86_64-3.9/em2ds.o build/temp.mingw_x86_64-3.9/msys64/home/zpic_user/zpic/em2ds/em2ds.cp39-mingw_x86_64.def -LC:/msys64/mingw64/lib/python3.9/config-3.9 -LC:/msys64/mingw64/lib -lpython3.9 -lm -lversion -lshlwapi -o ../lib/em2ds.cp39-mingw_x86_64.pyd
mingw32-make[1]: Leaving directory 'C:/msys64/home/zpic_user/zpic/python/source'

```

## Launching notebooks

Before you can launch your notebooks there is an additional step that must be performed: you must specify the path to the Jupyter python kernel. To this end edit the file:

 `/mingw64/share/jupyter/kernels/python3/kernel.json`

 Change line 3 from:

```text
   "D:/a/_temp/msys/msys64/mingw64/bin/python.exe",
```

To the path to your python binary, usually:

```text
   "C:/msys64/mingw64/bin/python.exe",
```

At this point you should be able to run the ZPIC notebooks by doing:

```
$ jupyter lab
[I 2022-04-26 17:07:06.013 ServerApp] jupyterlab | extension was successfully linked.
[I 2022-04-26 17:07:06.023 ServerApp] nbclassic | extension was successfully linked.
[I 2022-04-26 17:07:06.277 ServerApp] notebook_shim | extension was successfully linked.
[I 2022-04-26 17:07:06.307 ServerApp] notebook_shim | extension was successfully loaded.
[I 2022-04-26 17:07:06.307 LabApp] JupyterLab extension loaded from C:/msys64/mingw64/lib/python3.9/site-packages/jupyterlab
[I 2022-04-26 17:07:06.307 LabApp] JupyterLab application directory is C:/msys64/mingw64/share/jupyter/lab
[I 2022-04-26 17:07:06.307 ServerApp] jupyterlab | extension was successfully loaded.
[I 2022-04-26 17:07:06.317 ServerApp] nbclassic | extension was successfully loaded.
[I 2022-04-26 17:07:06.317 ServerApp] Serving notebooks from local directory: C:/msys64
[I 2022-04-26 17:07:06.317 ServerApp] Jupyter Server 1.16.0 is running at:
[I 2022-04-26 17:07:06.317 ServerApp] http://localhost:8888/lab?token=ffad91c015129aadca9715ddb4d1943007a2063397e88ab5
[I 2022-04-26 17:07:06.317 ServerApp]  or gcc
[I 2022-04-26 17:07:06.317 ServerApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 2022-04-26 17:07:06.357 ServerApp]

    To access the server, open this file in a browser:
        file:///C://Users/zpic_user/AppData/Roaming/jupyter/runtime/jpserver-4712-open.html
    Or copy and paste one of these URLs:
        http://localhost:8888/lab?token=ffad91c015129aadca9715ddb4d1943007a2063397e88ab5
     or http://127.0.0.1:8888/lab?token=ffad91c015129aadca9715ddb4d1943007a2063397e88ab5

```

Just access the Jupyter server by one of the methods specified above, and open the `README.ipynb` notebook to verify that everything is running.
