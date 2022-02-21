---
title: Time centering of quantities
layout: single
usemathjax: true

permalink: /about/timecenter/

sidebar:
  nav: "about"
---

Time integration of particle quantities in these codes is done using a leap frog method. For this purpose particle positions are defined at integral time steps, $t_i$, $t_{i+1}$, etc., whereas velocities (or generalized velocities in the relativistic versions) are defined at half time steps $t_{i-1/2}$, $t_{i+1/2}$, etc.

As a consequence, the electromagnetic fields, $E$ and $B$, as well as charge density, $\rho$, are also defined at integral timesteps $t_i$. Energy will also be defined at integral time steps. To this end, particle energy is calculated during the particle advance, using time centered values of the velocity.

Current density, $j$, is defined at half time steps $t_{i+1/2}$.
