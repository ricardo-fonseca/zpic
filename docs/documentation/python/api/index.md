---
title: Python API reference
layout : single
#usemathjax: true
permalink: /documentation/python/api/

sidebar:
  nav: "docs"
---

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

The complete API can also be viewed here:

* [em1d](em1d.html) - 1D EM code (finite-difference)
* [em2d](em2d.html) - 2D EM code (finite-difference)
* [em1ds](em1ds.html)- 1D EM code (spectral)
* [em2ds](em2ds.html) - 2D EM code (spectral)
* [es1d](es1d.html) - 1D electrostatic code (spectral)