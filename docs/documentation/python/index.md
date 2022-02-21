---
title: Using ZPIC with python
permalink: /documentation/python/

layout: single
toc: true
toc_label: Using ZPIC with python

sidebar:
  nav: "docs"
---

ZPIC simulations can be run directly from a python environment, allowing detailed control over the simulation parameter, and easy access to simulation results. To get started with running ZPIC simulations be sure to check the [Getting started](start) page.

## Documentation

* [Getting started](start)
* [Initializing particle species](particles)
* [Accessing simulation data](data)
* [Saving simulation output](save)
* [Using laser pulses](laser)
* [Setting a moving simulation window](movwindow)
* [Using electric current smoothing](smooth)
* [Setting initial EM fields](eminit)
* [Setting external EM fields](emext)
* [Setting the EM field solver type in spectral codes](spectralsolver)

## Example notebooks

ZPIC also includes several Jupyter notebooks exemplifying code usage, be sure to check the [examples](../../examples) section of the website.

## Python API Reference

The ZPIC python modules are extensively documented using docstrings, you can acess the documentation using the python `help()` command:

```text
>>> help( em2d.Density )
Help on class Density in module em2d:

class Density(builtins.object)
 |  Density(type='uniform', start=0.0, end=0.0, n=1.0, custom_x=None, custom_y=None)
 |  
 |  Class representing charge density profiles for particle species
 |  initialization
(...)
```

The complete API can also be consulted here:

* [em1d](api/em1d.html) - 1D EM code (finite-difference)
* [em2d](api/em2d.html) - 2D EM code (finite-difference)
* [em1ds](api/em1ds.html)- 1D EM code (spectral)
* [em2ds](api/em2ds.html) - 2D EM code (spectral)
* [es1d](api/es1d.html) - 1D electrostatic code (spectral)
