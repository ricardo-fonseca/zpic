---
title: Particle species initialization
permalink: /documentation/c/particles

usemathjax: true
layout: single
toc: true
sidebar:
  nav: "docs"
---

## Introduction

Adding particle data to a ZPIC simulation is done using the `species` and `n_species` parameters of the `sim_new()` function. The `species` parameter is an array of `t_species` structure. Before calling `sim_new()` each item in this array must be initialized using the `spec_new()` function:

```c
void spec_new( t_species* spec, char name[], float m_q, int ppc[],
    float ufl[], float uth[], int nx[], float box[], float dt,
    t_density* density );
```

The parameters for this function are as follows:

* `spec` - Pointer to `t_species` variable being initialized.
* `name` - Name for the species. This is used for diagnostic information only.
* `m_q` - Mass over charge ratio for particles in the species in simulation units (e.g. -1 for electrons)
* `ppc` - Reference number of particles per cell in each direction (scalar `int` in 1D codes)
* `ufl` - Initial fluid velocity for particle species (`float[3]`), set to NULL to use a 0 fluid velocity
* `uth` - Initial thermal velocity for particle species (`float[3]`), set to NULL to use a 0 thermal velocity
* `nx` - Simulation grid size
* `box` - Simulation physical box size in simulation units
* `dt` - Simulation time step in simulation units
* `density` - Pointer to `t_density` variable describing initial charge density, set to `NULL` to use a uniform density.

The following example adds a particle species with "step like" density profile:

```c
void sim_init( t_simulation* sim ){

    // Simulation box
    int   nx[2]  = { 1500, 128 };
    float box[2] = { 30.0, 25.6 };

    // Initialize particles
    const int n_species = 1;

    // Use 4x2 particles per cell
    int ppc[] = {4,2};

    // Density profile
    t_density density = { .type = STEP, .start = 30.0 };

    t_species* species = (t_species *) malloc( n_species * sizeof( t_species ));
    spec_new( &species[0], "electrons", -1.0, ppc, NULL, NULL, nx, box, dt, &density );

    //...
    // Initialize Simulation data
    sim_new( sim, nx, box, dt, tmax, ndump, species, n_species );
}
```

All particles in the same species share these parameters; if you require any of these parameters to change (e.g. initial fluid velocity) you will need to use additional species.

## Particle injection - 1D

In ZPIC all particles in a species have the same charge. Initializing a given density profile is done by varying the number of particles placed inside each simulation cell. The `ppc` parameter of the `spec_new()` function discussed above sets the reference number of particles for a density of 1 (in simulation units). The individual particle charge $q_p$ will respect:

$$
n = \frac{ ppc \, q_p }{\Delta x}
$$

with $\Delta x$ being the cell size, meaning that a cell with $ppc$ particles will have a charge density of $n$.

To choose a different initial density profile we need to reate a `t_density` variable specifying the density profile parameters and use it in the `density` parameter of our `spec_new()` call. The `t_density` structure has the following fields:


| Density parameters ||
|---|---|
| n | Reference density (default 1.0) |
| type | Density profile type: UNIFORM (default), STEP, SLAB, RAMP, or CUSTOM |
| start | Start of the particle injection region (STEP, SLAB, RAMP) |
| end  | End of the particle injection region (SLAB,RAMP) |
| ramp[2]  | Density at the start (ramp[0]) and end (ramp[1]) of the particle injection region for type RAMP, normalized to $n$ |
| custom | Pointer to a function defining an arbitrary density profile, normalized to $n$ |

After defining a `t_density` variable, this must be supplied (by reference) as the last parameter of the `spec_new()` function:

```c
// Uniform density example
t_density density = { .type = UNIFORM };

// (...)
spec_new( &species[0], "electrons", -1.0, ppc, NULL, NULL, nx, box, dt, &density );
```

### Density profile types

#### Uniform

Injects a uniform density `n`.

```c
// Uniform density example
t_density density = { .n = 1, .type = UNIFORM };
```


#### Step

Injects a step like density profile starting at position `start`, such that the density before this position is set to 0, and after this position it is set to `n`.

```c
// Step density example
t_density density = { .type = STEP, .start = 17.5 };
```

#### Slab

Injects a density profile that has a value of `n` between `start` and `end` and 0 otherwise.

```c
// Slab density example
t_density density = { .type = SLAB, .start = 17.5, .end = 22.5 };
```

#### Ramp

Injects a density profile that varies linearly from `ramp[0]` at position `start` to `ramp[1]` at position `end`. The values of the `ramp` parameter are normalized to the reference density `n`, meaning that `ramp = 1.0`, corresponds to a density `n`.

```c
// Ramp density example
t_density density = {
   .type = RAMP,
   .start = 17.5,
   .end = 22.5,
   .ramp = {1.0, 2.0}
};
```

#### Custom

The custom type injects a density profile set by a user defined function. The function must accept two parameters: one parameter of type `float` (the position at which the density is to be evaluated in simulation units) and a `void *` pointer (which can be used to send additional data to the function). It must return the density value as a value of type `float`. Note that density values are normalized to `n`.

```c
// Custom density example
// The density will oscillate between 0 and 0.1
// note the n parameter below

float custom_n0( float x, void *data ) {
    return sin(x/M_PI)*sin(x/M_PI);
}

//(...)

t_density density = {
  .type = CUSTOM,
  .n = 0.1,
  .custom = &custom_n0,
  .custom_data = NULL,
};
```

## Particle Injection - 2D

Particle injection in 2D is similar to 1D, with a couple of differences:

* The `ppc` parameter is a two element integer array specifying the number of particles per cell in both directions. The total number of particles per cell will be the product of these 2 numbers.
* The `RAMP` density profile is not available (use the `CUSTOM` profile instead)
* The `CUSTOM` density profile requires that the density can be described by a separable function in the $x$ and $y$ coordinates, i.e., $n(x,y) = n_x(x) \times n_y(y)$, and these two functions are defined separably.

The `t_density` structure in 2D has the following fields:

| Density parameters ||
|---|---|
| n | Reference density (default 1.0) |
| type | Density profile type: `UNIFORM` (default), `STEP`, `SLAB`, or `CUSTOM` |
| start | Start of the particle injection region (`STEP`, `SLAB`) |
| end  | End of the particle injection region (`SLAB`) |
| custom_x | Pointer to a function defining the $n_x(x)$ function of a custom density profile $n(x,y) = n_x(x) \times n_y(y)$, normalized to `n` |
| custom_y | Pointer to a function defining the $n_y(y)$ function of a custom density profile $n(x,y) = n_x(x) \times n_y(y)$, normalized to `n`  |

### Density profile types

The types `UNIFORM`, `STEP`, and `SLAB` work exactly the same way as in the 1D version.

#### Custom

The custom type injects a density profile set by two user defined functions. The functions must accept two parameters: one parameter of type `float` (the position $x$ or $y$ at which the density is to be evaluated in simulation units) and a `void *` pointer (which can be used to send additional data to the function). They must return the density value as a value of type `float`. Note that density values are normalized to `n`.

```c
// 2D Custom density example
// The 2D density will be the product of nx * ny
// It will oscillate between 0 and 0.1 - note the n parameter below

float nx( float x, void *data ) {
  return sin(x/M_PI)*sin(x/M_PI);
}

float ny( float y, void *data ) {
  return cos(y/M_PI)*cos(y/M_PI);
}

(...)

t_density density = {
  .type = CUSTOM,
  .n = 0.1,
  .custom_x = &nx,
  .custom_y = &ny,
};
```
