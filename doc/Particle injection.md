# Particle injection - 1D

Particle injection is controlled through a _t\_density_ structure that is supplied to the _spec\_new()_ routine when initializing the particle species. It accepts the following parameters:


| Density parameters ||
|---|---|
| n | Reference density (default 1.0) |
| type | Density profile type: UNIFORM (default), STEP, SLAB, RAMP, or CUSTOM |
| start | Start of the particle injection region (STEP, SLAB, RAMP) |
| end  | End of the particle injection region (SLAB,RAMP) |
| ramp[2]  | Density at the start (ramp[0]) and end (ramp[1]) of the particle injection region for type RAMP, normalized to _n_ |
| custom | Pointer to a function defining an arbitrary density profile, normalized to _n_ |

Particle injection also relies on the _ppc_ (particles per cell) parameter supplied to the _spec\_new()_ routine. The individual particle charge $q_p$is calculated so that:

$n = \frac{ ppc \, q_p }{\Delta x}$

with $\Delta x$ being the cell size, meaning that a cell with _ppc_ particles will have a charge density of _n_.

After defining a _t\_density_ variable, this must be supplied (by reference) as the last parameter of the _spec\_new()_ function:

```C
// Uniform density example
t_density density = { .type = UNIFORM };

// (...)
spec_new( &species[0], "electrons", -1.0, ppc, NULL, NULL, nx, box, dt, &density );
```



## Density profile types

### Uniform

Injects a uniform density _n_.

```C
// Uniform density example
t_density density = { .n = 1, .type = UNIFORM };
```


### Step

Injects a step like density profile starting at position _start_, such that the density before this position is set to 0, and after this position it is set to _n_.

```C
// Step density example
t_density density = { .type = STEP, .start = 17.5 };
```

### Slab

Injects a density profile that has a value of _n_ between _start_ and _end_ and 0 otherwise.

```C
// Slab density example
t_density density = { .type = SLAB, .start = 17.5, .end = 22.5 };
```

### Ramp

Injects a density profile that varies linearly from _ramp[0]_ at position _start_ to _ramp[1]_ at position _end_. The values of the _ramp_ parameter are normalized to the reference density _n_, meaning that _ramp_ = 1.0, corresponds to a density _n_.

```C
// Ramp density example
t_density density = {
   .type = RAMP,
   .start = 17.5,
   .end = 22.5,
   .ramp = {1.0, 2.0}
};
```


### Custom

The custom type injects a density profile set by a user defined function. The function must accept two parameters: one parameter of type _float_ (the position at which the density is to be evaluated in simulation units) and a void pointer (which can be used to send additional data to the function). It must return the density value as a value of type _float_. Note that density values are normalized to _n_.

```C
// Custom density example
// The density will oscillate between 0 and 0.1
// note the n parameter below

float custom_n0( float x, void *data ) {
	return sin(x/M_PI)*sin(x/M_PI);
}

(...)

t_density density = {
  .type = CUSTOM,
  .n = 0.1,
  .custom = &custom_n0,
  .custom_data = NULL,
};
```

# Particle Injection - 2D

Particle injection in 2D is similar to 1D, with two differences:

* The `ramp` density profile is not available, and
* The `custom` density profile requires that the density can be described by a separable function in the $x$ and $y$ coordinates, i.e., $n(x,y) = n_x(x) \times n_y(y)$, and these two functions are defined separably.

The `t_density` structure in 2D has the following options:

| Density parameters ||
|---|---|
| n | Reference density (default 1.0) |
| type | Density profile type: UNIFORM (default), STEP, SLAB, or CUSTOM |
| start | Start of the particle injection region (STEP, SLAB) |
| end  | End of the particle injection region (SLAB) |
| custom_x | Pointer to a function defining the $n_x(x)$ function of a custom density profile $n(x,y) = n_x(x) \times n_y(y)$, normalized to _n_ |
| custom_y | Pointer to a function defining the $n_y(y)$ function of a custom density profile $n(x,y) = n_x(x) \times n_y(y)$, normalized to _n_  |

## Density profile types

The types `UNIFORM`, `STEP`, and `SLAB` work exactly the same way as in the 1D version. 

### Custom

The custom type injects a density profile set by two user defined functions. The functions must accept two parameters: one parameter of type _float_ (the position x or y at which the density is to be evaluated in simulation units) and a void pointer (which can be used to send additional data to the function). They must return the density value as a value of type _float_. Note that density values are normalized to _n_.

```C
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
