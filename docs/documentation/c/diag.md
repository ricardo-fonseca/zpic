---
title: Simulation diagnostics
permalink: /documentation/c/diag

usemathjax: true

layout: single
toc: true
toc_label: Simulation diagnostics

sidebar:
  nav: "docs"
---

## Introduction

Simulation diagnostics are defined in the `sim_report()` function of the input file. This function will be called from the simulation main loop at a frequency defined by the `ndump` parameter of the `sim_init()` function. All input files must define this function.

Here's an example from a 2D Weibel instability simulation:

```c
void sim_init( t_simulation* sim ){
  //...
  
  // Diagnostic frequency - call sim_report() every 50 iterations
  int ndump = 50;

  // Initialize Simulation data
  sim_new( sim, nx, box, dt, tmax, ndump, species, n_species );

  //...
}

void sim_report( t_simulation* sim ){
  // Bx, By
  emf_report( &sim -> emf, BFLD, 0 );
  emf_report( &sim -> emf, BFLD, 1 );

  // Jz
  current_report( &sim -> current, 2 );

  // electron and positron density
  spec_report( &sim -> species[0], CHARGE, NULL, NULL );
  spec_report( &sim -> species[1], CHARGE, NULL, NULL );
}
```

## Electromagnetic Field Diagnostics

Electromagnetic field diagnostics are done through a call to the `emf_report()` function:
```c
void emf_report( t_emf *emf, char field, char fc )
```

This function may be called multiple times inside `sim_report()` with the following parameters:

| EMF report parameters| Description|
|---|---|
| field | Field to save. Must be one of `BFLD` or `EFLD` |
| fc | Field component to save. Must be one of 0(x), 1(y), or 2(z).  |

## Electric Current Diagnostics

Electric Current diagnostics are done through a call to the `current_report()` function:

```c
void current_report( const t_current *emf, const char field, const char fc )
```

This function may be called multiple times inside `sim_report()` with the following parameters:

| Current report parameters| Description|
|---|---|
| fc | Electric current component to save. Must be one of 0(x), 1(y), or 2(z). |

## Electric Charge Diagnostics

Global electric Charge diagnostics are done through a call to the `charge_report()` function:

```c
void charge_report( t_current *charge )
```

This function is only available in spectral code (es1d, em1ds, em2ds) as this quantity is not required for the finit difference models.

This function should only be called once inside `sim_report()`.

## Particle diagnostics

Particle species diagnostics are done through a call to the `spec_report()` function:

```c
void spec_report( t_species *spec, int rep_type, int pha_nx[], float pha_range[][2] )
```

This function may be called multiple times inside `sim_report()` with the following parameters:

| EMF report parameters| Description|
|---|---|
| rep_type | Type of diagnostic. Can be set to `CHARGE`, `PARTICLES`, or the `PHASESPACE()` macro |
| pha_nx | Dimensions (number of grid points) of the grid for phasespace density diagnostics  |
| pha_range | Physical grid size for phasespace density diagnostics |

This function is meant to be called per simulation species (set by the `spec` parameter), and acting only on the specified species.

The following types of reports are available:

* `CHARGE` - Saves the charge density of the species. This density is calculated in the simulation grid. The parameters `pha_nx` and `pha_range` should be set to `NULL`.
* `PARTICLES` -  Saves the complete particle dataset: position, (generalized) velocity, charge for all particles in the species. The parameters `pha_nx` and `pha_range` should be set to `NULL`.
* `PHASESPACE(ax1,ax2)` - Saves the specified phasespace density for this plot. The phasespace parameters are controlled through the following parameters:
  * `ax1`, `ax2` - Define the phasespace axis. Can be set to `X1`, `X2` (2D only), `U1`, `U2`, `U3` (EM codes) or `V1` (ES codes).
  * `pha_nx` - Defines the number of grid points in each direction for the phasespace density grid.
  * `pha_range` - Defines the physical limits (in simulation units) of the phasespace density grid in each direction.

In the following example, a simulation with two species saves diagnostic information for the charge density and $(x_1,u_1)$ phasespace density of species 0, and the complete particle dataset for species 1:

```c
void sim_report( t_simulation* sim ){

  spec_report( &sim -> species[0], CHARGE, NULL, NULL );
  spec_report( &sim -> species[0], PARTICLES, NULL, NULL );

  const int pha_nx[] = {1024,512};
  const float pha_range[][2] = { {0.0,20.0}, {-2.0,+2.0} };
  spec_report( &sim -> species[1], PHASESPACE(X1,U1), pha_nx, pha_range);
}
```
