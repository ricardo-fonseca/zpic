---
layout: single
title: Charge and current deposit
permalink: /about/deposit

usemathjax: true

sidebar:
  nav: "about"
---

After advancing particle quantities, the resulting current density is calculated on the grid at time $t + \Delta t/2$ to self-consistently evolve the electromagnetic field. This requires the current to be deposited on the grid using the same interpolation order (linear) as the interpolation scheme used to calculate the fields for the forces, otherwise non-physical self-forces can arise. This step will connect particle and grid quantities, closing the loop on particle/grid interaction. The charge / current density deposition scheme must also be chosen to match the electromagnetic field solver.

## Finite-difference EM codes

In the finite difference EM codes, the field solver only requires the electric current density for advancing the fields. Our electric current deposition scheme must account for the staggered Yee grid that we are using in our field solver; each electric current component must be defined at the same position as the corresponing electric current component. Furthermore, we must also insure that the charge continuity equation, as calculated using the same differential operators used in the field solver, is obeyed. This means that simply depositing $q v$ for each particle in not sufficient, and a more advanced method must be used. To this end, in ZPIC we use the method which is analytically equivalent to that proposed by [Villasenor 1992]. This method splits the particle trajectories into smaller segments falling inside a single cell; for each of these segments the algorithm essentially solves the continuity equation for the charge density evolution of the single particle segment, ensuring exact charge conservation.

## Spectral EM codes

In the spectral EM codes both the electric current and charge densities are required. The electric current needs to be deposited at $t + \Delta t/2$, like in the previous algorithm, but simply depositing $q v$ is enough, given that the details of the spectral solvers (in particular the handling of the longitudinal field) will ensure correct charge conservation. The electric charge density needs to be deposited at $t + \Delta t$, so the charge deposit is done after the particle push is complete.

Also, since there is no need for a staggered grid, all quantities are defined in the lower left corner of the cell.

## References

* Villasenor J and Buneman O 1992 _Comput. Phys. Commun._ __69__ 306–16
* Dawson J M 1983, _Rev. Mod. Physs_, __55__(2), 403-445
