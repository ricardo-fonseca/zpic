---
layout: page
title: Compiling
permalink: /compile/
---

## Pre-requisites

The only pre-requisites for ZPIC are a C99 compliant compiler and (GNU) make (check the [Getting a C99 compiler](#Getting-a-C99-compiler) below). Our recommendation is to use the GCC compiler, but other compilers such as Clang, Intel icc, and PGI pgcc are known to work. No additional / external libraries are required.

### Python module / notebooks

The Python modules are designed for Python 3. Python 2 is not supported. Compiling and using the Python module also requires the `Cython` and `numpy` Python packages, just install them using whatever package manager you prefer (e.g. `pip3 install Cython`). The modules can be used from any Python script, but using them in `Jupyter` notebooks is recommended.

## Compiling the code

### Individual code versions

Each of the code versions, available in different directories, is a complete standalone version, and can be compiled separately. Each of these directories contains a local `Makefile` that you can edit if you need to change any of the compiler parameters (or want to specify a different compiler).

To compile any of the code versions, just open a terminal window and navigate to the directory with the code that you want to use. You can then compile the code using the `make` command (`mingw32-make` on Windows). You can clean up all the generated object files and executable by issuing the `make clean` command.

#### Linux / Mac OS X example - em2d code

This example assumes you have the code on your home folder, inside a directory named "source". Just go into that directory and run type `make`:

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
cd source && /Library/Developer/CommandLineTools/usr/bin/make
python3 setup.py build_ext --build-lib=../lib
Compiling em1d.pyx because it changed.
Compiling em2d.pyx because it changed.
(...)
/temp.macosx-12-x86_64-3.9/Users/zamb/Source/zpic/em2ds/zdf.o build/temp.macosx-12-x86_64-3.9/em2ds.o -L/usr/local/lib -L/usr/local/opt/openssl@1.1/lib -L/usr/local/opt/sqlite/lib -o ../lib/em2ds.cpython-39-darwin.so
$ 
```

Libraries will be created in the `python/lib` directory, you need to add that library to your Python environment so you can use the modules.

***

## Getting a C99 compiler

Compiling ZPIC requires a standards compliant C99 compiler supporting C `float complex` types (GCC works great!) and GNU make.

### Linux

Most Linux distributions will let you install GCC using a package manager. For example, if your distribution uses the APT package manager, you can simply use `apt-get` to install gcc and make. Just open a terminal and do:

```bash
$ sudo apt-get install gcc make
```

This will be similar for most Linux distributions, just google for the instructions for your specific distribution.

### Mac OS X

For Mac OS X all you need are the Xcode Command Line Tools. You don't need the full [Xcode](https://developer.apple.com/xcode/) package, but if you have it installed just make sure you have it up to date.

To install the Xcode Command Line tools just open the Terminal and type the following command:

```text
$ xcode-select --install
```

If you have the tools already installed you should see something like this:

```text
$ xcode-select --install
xcode-select: error: command line tools are already
installed, use "Software Update" to install updates
```

Otherwise a dialog box should pop up. Just click "Install" to download and install the Xcode Command Line Tools.

### Windows

There are many ways to get GCC working under Windows. The recommended way is to use the mingw-w64 project: [http://mingw-w64.org/](http://mingw-w64.org/). Just download and run the installer from the Downloads page. We recommend the Mingw-builds package: [http://mingw-w64.org/doku.php/download/mingw-builds](http://mingw-w64.org/doku.php/download/mingw-builds).

The default settings work well. To launch a terminal with the appropriate path just click the Start button, and choose "Run terminal" under the "MinGW-W64 project" folder.

Please note that the `make` command will be installed as a rather criptic `mingw32-make`, but this is actually GNU make:

```
C:\>mingw32-make --version
GNU Make 4.2.1
Built for i686-w64-mingw32
Copyright (C) 1988-2016 Free Software Foundation, Inc.
Licence GPLv3+: GNU GPL version 3 or later <http://gnu.org/licences/gpl.html>
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.
```

#### Windows example - 1D electrostatic (es1d)

Launch a terminal by clicking the Start button, and choosing “Run terminal” under the “MinGW-W64 project” folder. This example assumes you have the source code on your root folder.

```
C:\>cd zpic\es1d
C:\zpic\es1d>mingw32-make
gcc -Ofast -g -std=c99 -pedantic -c charge.c -o charge.o
gcc -Ofast -g -std=c99 -pedantic -c field.c -o field.o
(...)
gcc -Ofast -g -std=c99 -pedantic  charge.o field.o particles.o grid.o fft.o random.o timer.o main.o simulation.o zdf.o -o zpic.exe
C:\zpic\es1d>
```

#### Note regarding MSVC

MSVC does not support the required C99 features (in particular the Complex type support which is crucial for the spectral versions of ZPIC) and so it cannot be used.
