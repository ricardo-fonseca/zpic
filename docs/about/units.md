---
title: ZPIC units
layout: single
usemathjax: true

permalink: /about/units/

sidebar:
  nav: "about"
---

The current implementation uses the same normalized units as the OSIRIS code. This requires choosing a reference normalization frequency, $\omega_n$. Time is normalized to $1/\omega_n$. (Proper) velocities are normalized to the speed of light, $c$. Space is normalized to $c/\omega_n$. The fields are then normalized appropriately.

The density is normalized to $\omega_n^2$ (the normalization frequency squared). So if the density is 1 at a given location then the normalization frequency is the plasma frequency at that location.

## Normalized units

The normalized ($'$) quantities refer to the original quantities using the following expressions:

$ x' = \frac{\omega_n}{c} x $

$ v' = \frac{v}{c} $

$ u' = \frac{u}{c} = \frac{\gamma v}{c} $

$ E' = e \frac{c / \omega_n }{m_e c^2} E $

$ B' = e \frac{c / \omega_n }{m_e c^2} B $

## A note on laser pulses

When using laser pulses in ZPIC simulations the laser frequency is set using the normalized units. This means that if the laser frequency is chosen to be 1, then the normalization frequency ($\omega_n$) is the laser frequency and the density is normalized to the critical density (for that laser frequency).
