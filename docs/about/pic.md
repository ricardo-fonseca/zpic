---
layout: single
title: Particle-in-Cell algorithm
usemathjax: true
permalink: /about/pic

sidebar:
  nav: "about"
---

Based on the highly nonlinear and kinetic processes that occur in many plasma physics scenarios, we must use a kinetic description of the plasma in our numerical model. Solving the plasma kinetic equations (e.g. Vlasov or Focker-Planck equations) numerically is not achievable in realistic geometries with nowadays computers, and presently the best choice is that of “particle simulations”. This follows the pioneering work of Buneman and Dawson in the early 1960’s [Buneman 1959, Dawson 1962] that showed we can accurately model the collective behavior of plasmas by following the evolution of a relatively small number (when compared to the number of particles in a real plasma) of simulation particles. This lead to the concept of a macro particle, that represents many plasma particles, and allows for the simulation of large spatial regions on the order of several tens of collisionless plasma skin depths.

As a first approach we might consider calculating the interaction between these particles directly, implementing a so called particle-particle algorithm. However, even using macro particles, realistic simulations would still lead to enormous requirements in computing power, requiring months of computing time in present day state of the art machines, since the number of operations required for such algorithm will scale roughly as the square of the number of particles used.

## The particle-in-cell algorithm (PIC)

The solution to this problem comes from the fact that we are generally interested in situations where the relevant physics is dominated by the collective effects of the plasma, meaning that short-range collisions can be neglected. In this sense, rather than interacting particles directly we will do so through the electromagnetic fields, which we will sample on a grid. The field values are interpolated at particle positions to calculate forces acting on particles, and particle motion will be used to determine the electrical current density used to solve the field equations. This is a so called particle-mesh algorithm, which greatly reduces the computing power requirements. With this algorithm, the number of operations for a simulation with $N_P$ particles and $N_{cells}$ grid cells will scale as $\sim \alpha \times N_P + \beta(N_{cells})$, where $\alpha$ is an integer constant (on the order of $\sim 10^3$) and $\beta$ is a function of the number of cells, usually scaling as $\sim \mathcal{O}(N_{cells})$ or $\sim \mathcal{O}( N_{cells} \log N_{cells})$ depending on the field solver type. The first term generally dominates, and the number of operations now scales linearly with the number of particles, and an extraordinary performance boost. However, this is gained at the expense of spatial resolution, given that we are now discretizing our fields on a grid thus limiting the shortest wavelengths that we can study.

{% include figure image_path="/assets/images/pic_loop.svg" alt="PIC loop" caption="__PIC loop__ The _i_ subscripts indicate particle quantities while the _j_ subscripts indicate grid quantities." %}

There are a variety of PIC codes in common use, differentiated by the kinds of forces retained in the model. The simplest is the electrostatic force, described by a Poisson equation. More complex are the Darwin (non-radiative electromagnetics) and the fully electromagnetic models. PIC codes generally have four important procedures in the main iteration loop. The first is the deposit, where some particle quantity, such as a charge, is accumulated on a grid via interpolation to produce a source density. Various other quantities can also be deposited, such as current densities, depending on the model. The second important procedure is the field solver, which solves Maxwells equations or a subset to obtain the electric and/or magnetic fields from the source densities. Once the fields are obtained, the particle forces are found by interpolation from the grid, and finally the particle coordinates are updated, using Newtons second law and the Lorentz force.

The image above illustrates a typical iteration loop in a PIC code [Hockney 1981, Dawson 1983, Birdsall 1991]. After setting the initial conditions for the problem being simulated, the main simulation loop will advance particles advanced using the Lorentz force resulting from the self-consistent electromagnetic field in the simulation. These fields are defined on a grid, and need to be interpolated at the particle positions using some form of weighing function, effectively connecting grid and particles quantities. Traditionally PIC codes use linear interpolation, that combines the values of the two nearest grid points in each dimension, which is acceptable for many simulation problems, with iteration counts up to $\sim 10^4$ [Villasenor 1992].

### References

* Buneman O 1959 _Phys. Rev._ __115__ 503-517.
* Dawson J M 1962, _Phys. Fluids_ __5__, 445-459.
* Dawson J M 1983, _Rev. Mod. Phys._, __55__(2), 403-445
* Hockney R W and Eastwood J W 1981, _Computer Simulation Using Particles_, McGraw-Hill, New York.
* Birdsall C K and Langdon  A B 1991, _Plasma Physics via Computer Simulations_, IoP Publishing, Bristol, UK.
