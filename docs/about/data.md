---
layout: single
title: Data structures in ZPIC
permalink: /about/data

usemathjax: true

sidebar:
  nav: "about"

toc: true
toc_label: Data structures
---

## Particle data

### Positions

Given that all grid related operations for particles, namely field interpolation and charge/current deposition, require the particle position inside the simulation grid cell, particle positions in ZPIC are stored using two values: an integer value (`ix`) specifying the cell where the particle and a floating-point value (`x`) specifying the position inside the cell. This choice, as oposed to say storing only the particle position in reference to the simulation box boundary, increases significantly the accuracy of the interpolation calculation and ensures uniform numerical properties over the whole simulation grid. It also allows us to use single precision for particle positions, which would not be feasible otherwise. Position values are normalized to cell size, and are always in the $[ 0, 1 )$ range (meaning that a value $x = 1$ is oustside of the current cell). In 2D we use 2 integers (`ix`, `iy`) and 2 floating point values (`x`, `y`).

### Velocities

All EM codes are relativistic, meaning that they use relativistic equations of motion for advancing the particles. These codes are expected to deal with multi-GeV particle energies, and as such, in ZPIC velocities are store as proper velocities, $\mathbf{u} = \gamma \mathbf{v}$, where $\gamma$ is the relativistic Lorentz factor, normalized to the speed of light, $c$. This avoids issues with loss of precision as the particle velocity approaches $c$. Also note that all EM codes use 3 velocity components, `ux`, `uy` and `uz`. Components not in the simulation domain (`y`, `z` in 1D, and `z` in 2D simulations) will not result in particle motion, but will contribute to electric current deposition.

Electro-static codes are not relativistic, using classical equations of motion for advancing the particles, so only normal velocities are stored. For consistency in the choice of units, we also normalize these to the speed of light. Only velocity components in the simulation domain (`vx` in 1D) are stored.

## Electro-Magnetic fields

Field data in ZPIC is stored on a simple floating point grid. This grid includes additional cells before and after the simulation domain, that we refer to as "guard cells". For example, in a 1D simulation grid with `nx` simulation cells, cell indices `0` to `nx-1` correspond to actual simulation cells, whereas cell indices `-1` or `nx` would correspond to guard cells.

The use of guard cells greatly simplifies the implementation because the avoid the need to write specific code for fringe case, such as advancing the fields in the 1st cell, or interpolating the fields in the last one. After each field advance a separate routine updates the values in the guard cells according to the boundary conditions used.

## Electric current/charge density

The electric current/charge density is stored in a manner similar to the EM-field data. These grids also use guard cells, but the number of guard cells may be different from other grids.
