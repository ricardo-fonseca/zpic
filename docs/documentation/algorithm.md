---
title: Algorithm details
description: Implementation details of the ZPIC codes
layout: single
usemathjax: true
permalink: /documentation/algorithm/

sidebar:
  nav: "docs"

layout: single
toc: true
toc_label: Algorithm details
---

## Units

The current implementation uses the same normalized units as the OSIRIS code. This requires choosing a reference normalization frequency, $\omega_n$. Time is normalized to $1/\omega_n$. (Proper) velocities are normalized to the speed of light, $c$. Space is normalized to $c/\omega_n$. The fields are then normalized appropriately.

The density is normalized to $\omega_n^2$ (the normalization frequency squared). So if the density is 1 at a given location then the normalization frequency is the plasma frequency at that location.

|       Normalized units                   |
| ---------------------------------------- |
| $x' = \frac{\omega_n}{c} x$|
| $v' = \frac{v}{c}$ |
| $u' = \frac{u}{c} = \frac{\gamma v}{c}$ |
| $E' = e \frac{c / \omega_n }{m_e c^2} E$ |
| $B' = e \frac{c / \omega_n }{m_e c^2} B$ |

### Laser pulses

When using laser pulses in ZPIC simulations the laser frequency is set using the normalized units. This means that if the laser frequency is chosen to be 1, then the normalization frequency ($\omega_n$) is the laser frequency and the density is normalized to the critical densify (for that laser frequency)

## Time centering of quantities

Time integration of quantities in these codes is done using a leap frog method. For this purpose particle positions are defined at integral time steps, $t_i$, $t_{i+1}$, etc., whereas velocities (or generalized velocities in the relativistic versions) are defined at half time steps $t_{i-1/2}$, $t_{i+1/2}$, etc.

As a consequence, the electromagnetic fields, $E$ and $B$, as well as charge density, $\rho$, are also defined at integral timesteps $t_i$. Energy will also be defined at integral time steps. For this purpose particle energy is calculated during the particle advance, using time centered values of the velocity.

Current density, $j$, is defined at half time steps $t_{i+1/2}$.
