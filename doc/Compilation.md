# Compiling and Running ZPIC

## Pre-requisites

The only pre-requisites for ZPIC are a C99 compliant compiler and (GNU) make. Our recommendation is to use the GCC compiler, but other compilers such as Clang, Intel icc, and PGI pgcc are known to work. No additional / external libraries are required.

## Getting a C99 compiler

If you already have a C99 compiler on your system you can skip to the [Compiling the Code](#Compiling-the-Code) section below.

### Linux

Most Linux distributions will let you install GCC using a package manager. For example, if your distribution uses the APT package manager, you can simply use apt-get to install gcc and make. Just open a terminal and do:

```
$ sudo apt-get install gcc make
```

This will be similar for most Linux distributions, just google for the instructions for your specific distribution.

### Mac OS X

For Mac OS X all you need are the Xcode Command Line Tools. You don't need the full [Xcode](https://developer.apple.com/xcode/) package, but if you have it installed just make sure you have it up to date.

To install the Xcode Command Line tools just open the Terminal and type the following command:

```
$ xcode-select --install
```

If you have the tools already installed you should see something like this:

```
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

#### Note

MSVC does not support the required C99 features (in particular the Complex type support which is crucial for the spectral versions of ZPIC) and so it cannot be used.

## Compiling the code

Each of the code versions, available in different directories, is a complete standalone version, and can be compiled separately. Each of these directories contains a local `Makefile` that you can edit if you need to change any of the compiler parameters (or want to specify a different compiler).

To compile any of the code versions, just open a terminal window and navigate to the directory with the code that you want to use. You can then compile the code using the `make` command (`mingw32-make` on Windows). You can clean up all the generated object files and executable by issuing the `make clean` command.

#### Linux / Mac OS X example - 2D electromagnetic (em2d)

This example assumes you have the code on your home folder, inside a directory named "source":

```
$ cd source/zpic/em2d
$ make
gcc -c -Ofast -std=c99 -pedantic current.c -o current.o
gcc -c -Ofast -std=c99 -pedantic emf.c -o emf.o
(...)
gcc -Ofast -std=c99 -pedantic  current.o emf.o particles.o random.o timer.o main.o simulation.o zdf.o -o zpic
$
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

## Running the code

After compiling the code, you can run it simply by calling the `zpic` program that you have just created. Output files will be written in the same directory as the code is run.

We recommend that you run the code in a different directory from the source directory. For example on Linux / Mac OS X:

```
$ mkdir test
$ cd test
$ ../zpic
$ ../zpic
Starting simulation ...

n = 0, t = 0.000000
n = 1, t = 0.250000
(...)
```

Or on Windows:

```
C:\zpic\es1d>mkdir test
C:\zpic\es1d>cd test
C:\zpic\es1d\test>..\zpic
Starting simulation ...

n = 0, t = 0.000000
n = 1, t = 0.250000
(...)
```
