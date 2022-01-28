---
title: Moving simulation window
description: Simulating phenomena travelling close to the speed of light
usemathjax: true
permalink: /documentation/python/movwindow

layout: single
toc: true
toc_label: Moving simulation window

sidebar:
  nav: "docs"
---

## Introduction

The finite difference models can be run using a moving simulation window that moves at the speed of light along the _x_ direction. Please note that simulation is done in the lab reference frame; it is the simulation box that moves and follows relevant phenomena moving at, or close to, this speed, such as laser pulses or relativistic particle beams.

To perform this motion the simulation window will be shifted left by 1 cell, (optionally) injecting new simulation particles in the rightmost cell of the simulation box, and setting the EM fields to 0 in that cell. These shifts occur whenever the simulation time $$t$$ divided by the cell size $$\Delta x$$ is larger than the number of simulation window shifts $$N$$ plus 1, i.e., whenever $$t / \Delta x \ge N + 1$$.

## Using a moving simulation window

To use a moving window, you must call the `Simulation.set_moving_window()` method after creating the simulation object, e.g.:

```python
import em1d

sim = em1d.Simulation( 1000, 100.0, 0.09, species = ... )

sim.set_moving_window()
```

## Checking window motion

The `n_move` property of the simulation object keeps track of the number of grid cell shifts $$N$$ that have occured so far in the simulation. For example, to get the current position of the lower _x_ boundary of the simulation we could do:

```python
import em1d

sim = em1d.Simulation( 1000, 100.0, 0.09, species = ... )
sim.set_moving_window()
sim.run( 10.0 )

print( "lower boundary = {:g}".format( sim.n_move * sim.dx ) )

```
