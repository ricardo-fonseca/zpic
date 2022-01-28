---
title: Electric current smoothing
permalink: /documentation/c/smooth

usemathjax: true

layout: single
sidebar:
  nav: "docs"
---

For the finite difference models (em1d, em2d), the input parameters may also optionally define some smoothing (filtering) parameters to be applied to the electric current following the deposition of this quantity from particle motion. This is achieved through a call to the `sim_set_smooth()` routine:

```c
void sim_set_smooth( t_simulation* sim, t_smooth* smooth)
```

This routine should be called inside `sim_init()`, somewhere after the call to `sim_new()`. Smoothing parameters are defined in the using a `t_smooth` structure:

| Smoothing parameters| Description|
|---|---|
| xtype | Smoothing type along x. Must be one of `NONE`, `BINOMIAL` or `COMPENSATED` |
| xlevel | Number of passes of the smoothing kernel along x |
| ytype  | (2D only) Same as xtype for the y direction |
| ylevel  | (2D only) Same as xlevel for the y direction |

For `BINOMIAL` smoothing, filtering is performed through a convolution of the grid data with a \[1,2,1\] kernel. This operation is repeated `n` times, where `n` is defined by the `xlevel | ylevel` parameters.

For `COMPENSATED` smoothing, the code will first apply `BINOMIAL` smoothing as described above, and will perform an additional pass with a specially calculated kernel that ensures that the filter transfer function near $k=0$ will be on the order of $1 + O(k^3)$, zeroing out $k^2$ dependencies. In pratice this leads to a filter that is flatter near $k=0$ and has a sharper drop near the cutoff frequency.

The following example sets up a 4th order `COMPENSATED` smoothing along the _x_ direction:

```c
void sim_init( t_simulation* sim ){

  (...)

  // Initialize Simulation data
  sim_new( sim, nx, box, dt, tmax, ndump, species, n_species );

  (...)

  // Set current smoothing (this must come after sim_new)
  t_smooth smooth = {
    .xtype = COMPENSATED,
    .xlevel = 4
  };

  sim_set_smooth( sim, &smooth );

  (...)
}
```
