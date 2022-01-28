---
title: Particle species
description: Initializing particle data in ZPIC simulations
permalink: /documentation/python/particles
usemathjax: true

layout: single
toc: true
toc_label: Particle species

sidebar:
  nav: "docs"
---

## Introduction

Adding particle species to the simulation is done using the `species` keyword of the `Simulation` class.  Users may use an an arbitrary number of particle species, defining the required plasma properties, such as density, temperature, fluid velocity, etc. Particle species must be created before creating the simulation object, e.g.:

```python
right = zpic.Species( "right", -1.0, ppc, ufl = [ 0.2,0,0], uth = [0.001,0.001,0.001] )
left  = zpic.Species( "left",  -1.0, ppc, ufl = [-0.2,0,0], uth = [0.001,0.001,0.001] )

sim = zpic.Simulation( nx, box, dt, species = [right,left] )
```

## The Species objects

When creating a particle species object, the `Species` class constructor takes the following parameters:

* `name` - Name of the particle species
* `m_q` - Mass over charge ratio for the particle species, in normalized units. For electrons this would be -1.
* `ppc` - Number of particles per cell. In 2D this should be a 2 element integer list specifying the number of particles per cell in each direction.

Additionally it may also take the following optional keyworkds:

* `ufl` - Initial fluid (generalized) velocity for the particles, defaults to [0,0,0] (no fluid velocity).
* `uth` - Initial thermal velocity for the particles, defaults to [0,0,0] (no thermal velocity).
* `density` - Density profile for the particle species, defaults to using a uniform density $n = 1$, see below for details.
* `n_sort` - Particle sorting frequency, defaults to 16.

The particles per cell parameter (`ppc`) of the species object specifies the number of particles per cell that correspond to the reference density (1.0). If you require a large range of density values, this parameter should be set to a reasonably high value, otherwise you may have cells with zero particles in the low density regions.

The particle sorting frequency (`n_sort`) refers to a performance optimization that improves memory cache use. To this end, particles are reordered in the particle buffer every `n_sort` iterations so that particles in the same simulatio cell are together in the buffer. If you don't want this behavior you may safely disable it by doing `n_sort = 0`.

## 1D Density profiles

ZPIC defaults to using a uniform density profile for particle species. To use a different density profile the user must create a `Density` object and set it as the `density` keyword when creating the particle species. Different density values are obtained by placing more or less particles per cell. As mentioned above, the number of particles per cell corresponding to a density value of 1.0 is set by the parameter `ppc`.

The following example creates an electron species using 128 particles per cell starting at position 17.5.

```python
density = zpic.Density( type = "step", start = 17.5 )
electrons = zpic.Species( "electrons", -1.0, 128, density = density )
```

The specific type of density profile is set through the `type` keyword when creating the `Density` object. Available options are:

* __"uniform"__ - (default) Uniform density
* __"step"__ - Step like density profile (0 before the step, 1 afterwards), the keyword `start` specifies the position of the step
* __"slab"__ - Slab like density profile (1 inside the slab, 0 outside), the keywords `start`/`end` specify the beginning/end of the slab
* __"ramp"__ - Linear density ramp, with 0 density outside the ramp region. The keywords `start`/`end` specify the beginning/end of the ramp, and the keyword `ramp` specifies the required density at these two points.
* __"custom"__ - Custom density profile. Density is specified through a custom function specified through the `custom` keyword. This function takes a single parameter specifying the position and returns the required density at that point.

## 2D Density profiles

Setting the density profile in 2D works the same way as in 1D. Available `Density` options are:

* __"uniform"__ - (default) Uniform density
* __"step"__ - Step like density profile along $x$ (0 before the step, 1 afterwards), the keyword `start` specifies the $x$ position of the step
* __"slab"__ - Slab like density profile along $x$ (1 inside the slab, 0 outside), the keywords `start`/`end` specify the $x$ poitison of the beginning/end of the slab
* __"custom"__ - Custom density profile. Density is specified using 2 custom functions specified through the `custom_x` and `custom_y` keywords. The density at each point will be `custom_x(x) * custom_y(y)`.

## Diagnostics

The `Species` object includes special methods / properties to simplify getting diagnostic information:

* The `Species.energy` property has the kinetic energy of the species. Bare in mind that this energy calculated during the particle advance, so it will be 0 before the simulation starts.
* The `Species.charge()` method returns a NumPy array with the same shape as the simulation grid with the charge density of the species.
* The `Species.phasespace(quants, pha_nx, pha_range)` method returns a NumPy array with the phasespace density of the species. Method parameters are:
  * `quants` - list of 2 quantities that make the phasespace axis, possible values are `'x1'`, `'x2'` (2D codes only), `'u1'`, `'u2'`, and `'u3'` specifying positions and generalized velocities, respectively.
  * `pha_nx` - Dimensions of phasespace grid
  * `pha_range`- Physical limits of the phasespace grid.
