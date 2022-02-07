---
layout: single
title: Simulation particles
permalink: /about/particles

usemathjax: true

sidebar:
  nav: "about"
---

## Particle advance

In a PIC code, after interpolating the fields at particle positions, particle velocity and positions must then be advanced in time. In the ZPIC codes this is acomplished through a leap-frog scheme, that first advances velocity and then positions, which provides 2nd order accuracy in time. This requires that particles positions and velocities are not time-centered, with particle positions being known at time $t$, and particle velocity being known at $t−\Delta t/2$.  As a result, field quantities will also have to be known at time $t$, and the values interpolated at particle positions will be used to calculate the Lorentz force and advance particle momenta from $t− \Delta t/2$ to $t+ \Delta t/2$.

In the ZPIC EM codes the velocity advance is done using a (relativistic) Boris pusher [Boris 1970, Birdsall and Langdon 1991]: This method considers the effects of the electric and magnetic forces separately, splitting the particle advance into electric field impulses and magnetic rotations. This allows for easy calculation of the time centered relativistic gamma factor that is not affected by the magnetic rotation but is required to calculate the magnetic force, and results in a second order accurate momentum integration. After completing this integration, the momentum is known at time $t + \Delta t/2$, which is then used to advance particle positions in time from $t$ to $t + \Delta t$.

In the electrostatic codes the magnetic field is absent, so the velocity is straightforwardly advanced from $t− \Delta t/2$ to $t+ \Delta t/2$.

## References

* Boris JP 1970, _Proceedings of the 4th Conference on the Numerical Simulation of Plasmas_, Naval Research Laboratory, Washington, D. C., pp. 3-67, 1970.
* Birdsall C K and Langdon  A B 1991, _Plasma Physics via Computer Simulations_, IoP Publishing, Bristol, UK
