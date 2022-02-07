---
title: Laser Pulses
usemathjax: true
permalink: /documentation/python/laser

layout: single
toc: true
toc_label: Laser pulses

sidebar:
  nav: "docs"
---

## Introduction

To simplify the use of laser pulses in ZPIC simulations, the code includes the ability to launch laser pulses. The code allows the user to choose the frequency, amplitude, temporal envelope and polarization of the laser pulse, and will calculate the self-consistent electric and magnetic fields:

* The laser pulses are created at once: the full laser fields are super-imposed (added) on the simulation fields.
* The temporal envelope of the laser pulses is defined as a $\sin^2$ function; different durations for the rising and falling edges may be specified, and the user may also specify an interval where the laser pulse holds its intensity constant.
* All laser pulses propagate along the positive direction of the $x$ axis. Other directions are not currently supported.

## 1D Laser Pulses

To add a laser pulse to a simulation object `sim` we use the `add_laser()` method providing the laser pulse parameters as a `Laser` object, e.g.:

```python
sim = zpic.Simulation( nx, box, dt, species = electrons )
sim.add_laser( zpic.Laser( start = 17.0, fwhm = 2.0, a0 = 1.0, omega0 = 10.0 ))
```

In 1D the `Laser` class constructor accepts the following parameters:

* `fwhm` - Full width at half-max of the laser pulse. If set it overrides the `rise`, `flat`, and `fall` parameters;
* `rise`, `flat`, `fall` - Rise time (`rise`), flat time (`flat`) and fall time (`fall`) of the temporal envelope;
* `start` - Position of the starting point (right) of the laser pulse;
* `a0` - Normalized vector potential value at peak intensity of the laser pulse;
* `omega0` - Laser frequency in simulation units;
* `polarization` - Laser polarization in radians measured in reference to the $ y $ direction.

## 2D Laser Pulses

Using laser pulses in 2D simulations works the same way as in 1D, but the `Laser` class also allows for Gaussian laser pulses. Besides the parameters available in 1D, the `Laser` class constructor in 2D also accepts the following parameters:

* `type` - Specifies the laser pulse type, plane wave ("plane", the default) or Gaussian pulse ("gaussian")
* `W0` - Gaussian pulse waist size in simulation units
* `focus` - Position of focal plane for Gaussian pulse in simulation units, may be set outside of the simulation box
* `axis` - $y$ position of the propagation axis for Gaussian pulses in simulation units

