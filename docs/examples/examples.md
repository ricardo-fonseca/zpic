---
title: ZPIC Examples
layout: single
#usemathjax: true
permalink: /examples/

toc: true
toc_label: "ZPIC Examples"

sidebar:
  nav: "docs"

nbroot : https://github.com/ricardo-fonseca/zpic/blob/master/python/notebooks
---

ZPIC makes available a set of well documented Jupyter notebooks, with example problems of textbook and advanced plasma mechanisms, ranging from Debye shielding, to EM waves in density ramps and kinetic instabilities, as well as notebooks that reproduce and extend the work done in seminal plasma physics papers.

If you haven't, make sure you start with the [ZPIC]({{page.nbroot}}/tutorial/ZPIC.ipynb) notebook that covers the basic usage of ZPIC in these environments.

## Contributing

We invite contributions to the __ZPIC__ codes and repository of test problems, to be made freely available to the community. Come find us on [GitHub](https://github.com/ricardo-fonseca/zpic).

## Tutorials

Notebooks introducing code usage and functionalities:

* [ZPIC]({{page.nbroot}}/tutorial/ZPIC.ipynb) - Main tutorial file, **be sure to start here**;
* [Saving Results]({{page.nbroot}}/tutorial/Saving%20results.ipynb) - Saving simulation results to disk;
* [Animation]({{page.nbroot}}/tutorial/Animation.ipynb) - Ceating animations from your simulations;
* [External Fields]({{page.nbroot}}/tutorial/External%20Fields.ipynb) - Using external (constant) EM fields in your simulation;
* [Initial Fields]({{page.nbroot}}/tutorial/Initial%20Fields.ipynb) - Setting the initial EM fields for your simulation;
* [Density]({{page.nbroot}}/tutorial/Density.ipynb) - Setting the density profile of particle species;
* [Thermal distribution]({{page.nbroot}}/tutorial/Thermal%20Distribution.ipynb) - Setting the thermal distribution of particle species;
* [Cathode]({{page.nbroot}}/tutorial/Cathode.ipynb) - Implementing a cathode (particle injection from simulation wall);
* [Custom velocity distribution]({{page.nbroot}}/tutorial/Custom%20velocity%20distribution.ipynb) - Implementing custom (arbitrary) intial velocity distributions;
* [Laser Pulses]({{page.nbroot}}/tutorial/Laser%20Pulses.ipynb) - Launching laser pulses in ZPIC simulations.

## Classroom examples

Examples of ZPIC notebooks that can be used in a classroom to showcase some of the most fundamental plasma physics phenomena:

* [Field solver dispersion]({{page.nbroot}}/classroom/Field%20solver%20dispersion.ipynb) - Analyses the different dispersion properties of the field solvers used in finite difference and spectral codes;
* [Electron Plasma Waves]({{page.nbroot}}/classroom/Electron%20Plasma%20Waves.ipynb) - Focuses on electrostatic and electromagnetic plasma waves in unmagnetized plasmas;
* [O-X Waves]({{page.nbroot}}/classroom/O-X%20Waves.ipynb) - Studies electromagnetic waves in a magnetized plasma, in particular polarized either along the applied magnetic fields (O-waves) or perpendicular to it (X-waves);
* [R-L Waves]({{page.nbroot}}/classroom/R-L%20Waves.ipynb) - Studies electromagnetic waves in a magnetized plasma, in particular waves propagating along the applied magnetic field;
* [Debye Shielding]({{page.nbroot}}/classroom/Debye%20Shielding.ipynb) - Focuses on Debye shielding;
* [Faraday Rotation]({{page.nbroot}}/classroom/Faraday%20Rotation.ipynb) - Observation of the Faraday rotation effect;
* [Two-Stream]({{page.nbroot}}/classroom/Two-Stream.ipynb) - Numerical simulations of the Two-Stream instability using an electromagnetic code;
* [Electrostatic Two-Stream]({{page.nbroot}}/classroom/Electrostatic%20Two-Stream.ipynb) - Numerical simulations of the Two-Stream instability using an electrostatic code;
* [LWFA 1D]({{page.nbroot}}/classroom/LWFA%201D.ipynb) - Numerical simulations of a laser wakefield accelerator in 1D;
* [LWFA 2D]({{page.nbroot}}/classroom/LWFA%202D.ipynb) - Numerical simulations of a laser wakefield accelerator in 2D;
* [PWFA 1D]({{page.nbroot}}/classroom/PWFA%201D.ipynb) - Numerical simulations of a plasma wakefield accelerator in 1D;
* [PBWA 1D]({{page.nbroot}}/classroom/PBWA%201D.ipynb) - Numerical simulations of a plasma beat wave accelerator in 1D;
* [Weibel]({{page.nbroot}}/classroom/Weibel.ipynb) - Demonstration of the Weibel (electromagnetic filamentation) instability in the collision of neutral electron/positron plasma clouds;
* [Coulomb Collisions]({{page.nbroot}}/classroom/Coulomb%20Collisions.ipynb) - Study of Coulomb Collisions / Rutherford scattering;
* [Diamagnetic Drift]({{page.nbroot}}/classroom/Diamagnetic%20Drift.ipynb) - Focuses on the diamagnetic drift (or diamagnetic current), a fluid drift in which a current arrises in a magnetized plasma with a given density gradient;
* [E×B Drift]({{page.nbroot}}/classroom/ExB%20Drift.ipynb) - Studies of single particle **E**×**B** drifts;
* [Magnetic Bottle]({{page.nbroot}}/classroom/Magnetic%20Bottle.ipynb) - Numerical studies of the magnetic bottle confinement scheme;

## Landmark papers

Notebooks reproducing and extending the work done in seminal plasma physics papers:

* [Tajima and Dawson (1979)]({{page.nbroot}}/papers/Tajima%20and%20Dawson%201979.ipynb) - Laser Electron Accelerator, Physical Review Letters, Volume 43, Number 4, July 1979, [DOI: 10.1103/PhysRevLett.43.267](https://doi.org/10.1103/PhysRevLett.43.267);
* [Morse and Nielsen (1971)]({{page.nbroot}}/papers/Morse%20and%20Nielsen%201971.ipynb) - Numerical Simulation of the Weibel Instability in One and Two Dimensions, The Physics of Fluids, Volume 14, Number 4, April 1971, [DOI: 10.1063/1.1693518](https://doi.org/10.1063/1.1693518).
