---
layout: single
title: Field interpolation
permalink: /about/interpolation

sidebar:
  nav: "about"
---

In a PIC algorithm fields and charge/current density are discretized on a grid. Simulation particles, however, are free to move into any position inside the simulation domain. To advance particles we need to calculate the forces acting on them, so we interpolate the fields at the particle positions.

In ZPIC this is achieved using linear interpolation (bi-linear in 2D codes). Using the particle cell index the code determines which grid values need to be used, namely the neares grid points to the left/right of particle position, and calculates the value using the particle distance to these points. Note that this requires that the fields define an extra (guard) cell at the upper boundary of the simulation box for a particle sitting close to this boundary.

In the finite difference codes this operation must account for the staggered grid used for the Yee scheme (see the [fields](fields) section). For the spectral codes (including the electrostatic code) all quantities are defined in the lower left corner of the cell.
