---
title: Moving simulation window
permalink: /documentation/c/movwindow

layout: single
sidebar:
  nav: "docs"
---

The finite difference models can be run using a moving simulation window that moves at the speed of light. Please note that simulation is done in the lab reference frame; it is the simulation box that moves and follows relevant phenomena moving at, or close to, this speed, such as laser pulses or relativistic particle beams.

Using a moving window requires calling the `sim_set_moving_window()` routine:

```c
void sim_set_moving_window( t_simulation* sim )
```

This routine should be called inside `sim_init()`, somewhere after the call to `sim_new()`, e.g.:

```c
void sim_init( t_simulation* sim ){
    // (...)
    // Initialize Simulation data
    sim_new( sim, nx, box, dt, tmax, ndump, species, n_species );
    // (...)
    // Set moving window (this must come after sim_new)
    sim_set_moving_window( sim ); 
    // (...)
}
```

Spectral models currently do not support this feature.
