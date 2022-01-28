---
title: Electric Current Smoothing
description: Units used for ZPIC simulations
permalink: /documentation/python/smooth
usemathjax: true

layout: single
toc: true
toc_label: Electric Current Smoothing

sidebar:
  nav: "docs"
---

## Introduction

In the finite difference ZPIC codes (`em1d`, `em2d`) you may specify a digital filter that will be applied to the electric current density after this has been deposited by particle motion and before using it for advancing the EM fields.

The ZPIC codes implement 2 type of low pass spatial filters:

* __Binomial__ - Multi-pass binomial (gaussian) filtering, i.e., do a convolution with a [1,2,1] stencil multiple times.
* __Compensated__ - Same as binomial filtering followed by an appropriate compensator that eliminates $k^2$ dependency of the transfer function near $k=0$.

In 2D you may define different smoothing types along $x$ and $y$.

## Setting Current Smoothing Parameters

The electric current density smoothing (filtering) parameters are defined in a `Smooth` object, that is then used with the `Simulation.set_smooth()` method. Valid parameters for the `Smooth` class constructor are:

* `xtype` - Type of digital filtering to apply along the x direction, must be one of "none" (default), "binomial"  or "compensated"
* `xlevel` - Number of filtering passes to apply
* `ytype` / `ylevel` (2D only) - Same as `xtype` / `xlevel` for the y direction

The following example sets a 4 level compensated smoother for a `em1d` simulation:

```python
sim = em1d.Simulation(...)
sim.set_smooth( em1d.Smooth(xtype = "compensated", xlevel = 4) )
```
