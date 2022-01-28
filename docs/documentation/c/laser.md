---
title: Using Laser pulses
permalink: /documentation/c/laser

usemathjax: true
layout: single

toc: true
toc_label: Laser pulses
sidebar:
  nav: "docs"
---

## Introduction

Users may optionally define an arbitrary number of laser pulses to be added to the simulation. Each laser pulse is added through a call to the `sim_add_laser()` routine:

```c
void sim_add_laser( t_simulation* sim, t_emf_laser* laser)
```

This routine should be called inside `sim_init()`, somewhere after the call to `sim_new()`, e.g.:

```c
void sim_init( t_simulation* sim ){

  //...
  // Initialize Simulation data
  sim_new( sim, nx, box, dt, tmax, ndump, species, n_species );

  // Add laser pulse (this must come after sim_new)
  t_emf_laser laser = {
    .start = 17.0,
    .fwhm  = 2.0,
    .a0 = 2.0,
    .omega0 = 10.0,
    .polarization = M_PI_2
    };
  sim_add_laser( sim, &laser );

  // ...
}
```

## Laser pulse parameters - 1D

Laser pulse parameters are specified by a `t_emf_laser` variable supplied to the `sim_add_laser()` function. In 1D the `t_emf_laser` structure has the following fields:

| Laser parameters| Description|
|---|---|
| type | PLANE for a plane wave or GAUSSIAN for a gaussian beam |
| start | Front edge of the laser pulse |
| fwhm  | FWHM of the laser pulse duration |
| rise, flat, fall  | Rise / flat / fall time of the laser pulse |
| a0  | Normalized peak vector potential of the pulse |
| omega0 | Laser frequency |
| polarization | Laser polarization angle, 0 aligns $E$ field along $x_2$ |
| W0 | Gaussian beam waist |
| focus | Focal plane position |
| axis | Position of the optical axis |

All parameters are in normalized simulation units except for the normalized peak vector potential, $a_0$, which is adimensional.

Using the `fwhm` parameter will override the `rise`, `flat` and `fall` parameters. Specifically, it sets `rise = fwhm/2`, `flat = 0`, and `fall = fwhm/2`.

The following example launches a laser starting at position 17.0, with a (temporal) full width at half max of 2.0. The peak normalized vector potential is 2.0, and the laser frequency is 10.0. The polarization degree is $\pi/2$, which aligns the $E$ field along the $x_3$ direction. All values are in normalized simulation units.

```c
t_emf_laser laser = {
  .start = 17.0,
  .fwhm  = 2.0,
  .a0 = 2.0,
  .omega0 = 10.0,
  .polarization = M_PI_2
};
sim_add_laser( sim, &laser );
```

## Laser pulse parameters - 2D

In 2D ZPIC also supports Gaussian beam laser pulses (besides plane wave pulses). The `t_emf_laser` structure has the additional fields:

| Laser parameters| Description|
|---|---|
| type | `PLANE` (default) for a plane wave or `GAUSSIAN` for a gaussian beam |
| W0 | Gaussian beam waist |
| focus | Focal plane position |
| axis | Position of the optical axis |

All parameters are in normalized simulation units. The following is a 2D example of a gaussian laser pulse. It uses the same parameters as the previous example, set the beam focus waist to $W_0=4.0$, the focal plane position to $x=20.0$, and the propagation axis to $y=12.8$.

```c
t_emf_laser laser = {
  .type = GAUSSIAN,
  .start = 17.0,
  .fwhm  = 2.0,
  .a0 = 2.0,
  .omega0 = 10.0,
  .W0 = 4.0,
  .focus = 20.0,
  .axis = 12.8,
  .polarization = M_PI_2
};
sim_add_laser( sim, &laser );
```
