---
layout: single
title: Electromagnetic Fields
permalink: /about/fields

usemathjax: true

toc: true
toc_label: EM Fields

sidebar:
  nav: "about"
---

## Field Advance

After advancing the particles and depositing the charge/current density, the field quantities are then advanced from time $t$ to time $t + \Delta t$. In ZPIC codes this can either be achieved using a finite difference method (`em1d`, `em2d`) or a spectral method (`em1ds`, `em2ds`).

### Finite-difference field solver

The finite-difference field solver in ZPIC is a local Maxwell solver that is equivalent to the Finite Difference Time Domain (FDTD) method using a Yee mesh [Yee 1966], adding the contribution of the plasma currents to the electrical field advance equation. This method approximates the spatial derivatives using finite differences. Electric and magnetic quantities are defined in different in a staggered grid allowing for second order accuracy in space for this operator. In ZPIC we chose to place the charge point at the lower left corner of the simulation cell

The standard FDTD method requires that the electric and magnetic field to be also non-time centered. However, as we saw above, this is a requirement for the particle push. To overcome this issue, we follow the method proposed by Boris that splits the magnetic field integration into 2 steps. The algorithm first advances the magnetic field from time $t$ to time $t + \Delta t/2$ using the electric field values at time $t$, then performs a standard electric field advance from $t$ to $t + \Delta t$ using the magnetic field and current density defined at time $t + \Delta t/2$, and finally advances the magnetic field from $t + \Delta t/2$ to $t + \Delta t$. This method has been demonstrated to have the same accuracy as the FDTD method, while allowing the electric and magnetic fields to be time centered, and not requiring any additional storage for the magnetic field.

### Pseudo-spectral field solvers

Alternatively to (spatial) finite difference methods, we can use spectral methods to represent the spatial derivates [Dawson 1983 and references therein]. In these methods the field equations are advanced in Fourier space, making use of fast fourier transforms (FFT). To facilitate this we split $\mathbf{E}$ and $\mathbf{j}$ into transverse and longitudinal components with respect to the wave vector $\mathbf{k}$; $\mathbf{B}$ has only transverse components. The transverse parts of $\mathbf{E}$ and $\mathbf{B}$ are then advanced in time using the transverse part of $\mathbf{j}$. The longitudinal part of $\mathbf{E}$ is then calculated from the charge density by solving Poisson's equations. The pseudo-spectral methods therefore require that the charge density for the particles, defined at $t + \Delta t$ as also been deposited.

The simplest method for advancing the (transverse) field equations in time is similar to the strategy described for the finite difference solver described above: transverse $\mathbf{E}$ and $\mathbf{B}$ fields are advanced using a leap-frog scheme; to maintain time-centering of field quantities, we first advance $\mathbf{B}$ half a time step, then advance $\mathbf{E}$ a full time step, and finally advance $\mathbf{B}$ another half time step. This is known as the "Pseudo-spectral time domain" (__PSTD__) method [Liu 1997].

Alternatively, instead of using a finite difference in time approach, we may advance the field equations in time using an exact analytical solution in Fourier space. The method proposed by [Haber 1973] is exact for plasma currents that are constant during the time step, is free of electromagnetic wave numerical dispersion for wave numbers that are properly resolved by the simulation grid. This is known as the "Pseudo-spectral analytical time domain" (__PSATD__) method. Both transverse $\mathbf{E}$ and $\mathbf{B}$ are advanced simultaneously so the time-centering is preserved.

### Which field solver should I use?

All the methods available have different properties, making them more suitable in some scenarios. The main advantages of the __FDTD__ solver are the posibility of using a moving simulation window and/or absorbing boundary conditions. Pseudo-spectral methods have better dispersion relation properties and allow more freedom in the choice of simulation time steps, but are limited to periodic boundaries.

## Dispersion relation

All 3 solvers yield different dispersion relations for EM waves:

* __FDTD__ - The dispersion relation for EM waves is always sub-luminal; best results are obtained by choosing the largest possible time step (just below the CFL condition).
* __PSTD__ - The dispersion relation for EM waves is always super-luminal; best results are obtained by choosing smaller time steps.
* __PSATD__ - The dispersion relation for EM waves is always correct.

Please see the example [Dispersion relation]({{page.nbroot}}/classroom/Field%20solver%20dispersion.ipynb) notebook for a comparison of the 3 types of field solver.

## Neutralizing background

When using the finite difference models, the code will always behave as if the total charge density was initially zero, even when this is not the case. This effect is generally described as having the code add a neutralizing background; however, for finite difference models, this is simply a consequence of the way the field solver is implemented, and no actual addition takes place. So to model a plasma with a fixed ion background, it is enough to add the electron species, and no additional particle species are required.

For the spectral models this is no longer the case, and if required a neutralizing background must be explicitly added.

## References

* Dawson J M 1983, _Rev. Mod. Physs_, __55__(2), 403-445
* Haber I et al 1973 _Proceedings of the 6th conference on the Numerical Simulation of Plasmas_, Berkeley, CA
* Liu QH 1997, _Micro. Opt. Tech. Lett._, __15__(3), 158-165
* Yee K 1966 _IEEE Trans. Ant. Prop._  __AP14__ 302-307
