---
title: Saving simulation data
description: Saving simulation results to disk for post-processing
permalink: /documentation/python/save
usemathjax: true

layout: single
toc: true
toc_label: Saving simulation data

sidebar:
  nav: "docs"
---

## Introduction

You can run ZPIC from Python and access all the simulation data directly in memory which, for most situations, is the recommended way of using the code. However, if you your simulation takes a long time to compute, you may want to write diagnostic information to disk for post-processing later.

To this effect, the `EMF`, `Current` and `Species` classes (as well as the `Charge` class on spectral codes) include a `report()` method that can be used to save specific dianostic information to disk. Data will be saved in the `ZDF` format, as in the rest of the ZPIC framework. The includes `zdf` python module allows reading ZDF file data and metadata, check the documentation for details.

## Saving EM Fields

To save the current value of any of the EM fields we would use the `EMF.report( type, fc )` method where:

* `type` - selects which field is to be saved. Valid values are `'E'` for the electric field and `'B'` for the magnetic field. If the simulation uses external fields then the `'Ep'` and `'Bp'` options are also available, specifying the total (self consistent + external) electric/magnetic fields.
* `fc` - selects which field component is saved, `0` - _x_, `1` - _y_ or `2` - _z_.

Files will be saved in the `EMF` directory. The following example runs a simulation up to $$t = 10.0$$ and then saves $$E_z$$.

```python
sim = em1ds.Simulation( ... )

sim.run(10)

# Save Ez at t = 10.0
sim.emf.report("E",2 )
```

## Saving Current/Charge density

To save the electric current density we use the `Current.report(fc)` method, where `fc` specifies which component to save, e.g.:

```python
# Save Jz
sim.current.report(1)
```

Files are saved to the `CURRENT` directory.

For the charge density (spectral codes only) the `Charge.report()` method is used, taking no parameters:

```python
# Save rho
sim.charge.report()
```

Files are saved to the `CHARGE` directory.

## Particle data

To save diagnostic particle data we use the `Species.report( type, quants = [], pha_nx = [], pha_range = [] )` method:

* `type` specifies the actual type of diagnostic data to be saved. Must be one of:
  * `'particles'` - Raw particle data, including position and velocity of all particles
  * `'charge'`- Species electric charge density
  * `'pha'` - Species phasespace density, see parameters `quants`, `pha_nx` and `pha_range` below.

When choosing `type = 'pha'` (phasespace density) the following parameters must also be specified:

* `quants` - list of 2 quantities that make the phasespace axis, possible values are `'x1'`, `'x2'` (2D codes only), `'u1'`, `'u2'`, and `'u3'` specifying positions and generalized velocities, respectively.
* `pha_nx` - Dimensions of phasespace grid
* `pha_range`- Physical limits of the phasespace grid.

The next example saves a $$ x_1 - u_1 $$ phasespace:

```python
right.report("pha",                       # Phasespace density
    quants = ["x1","u1"],                 # Use 'x1' in the x axis and 'u1' in the y axis
    pha_nx = [120,256],                   # Use a 120 x 256 grid
    pha_range = [[0,4*np.pi],             # Use [0, 4*pi] range for the x1 axis
                 [-1.5,1.5]]              # and [-1.5, 1.5] for the u1 axis
)
```

## Using the report function

To simplify the process of writing diagnostic data regularly over the course of a simulation, ZPIC codes allow defining a `report` function that will be called automatically at the end of every iteration when running the simulation using the `Simulation.run()` method. This function must accept as a single argument a simulation object, allowing access to all simulation data and the use of the `save` methods discussed above.

The following example saves the $E_x$ field at every 10 timesteps until the simulation reaches `tmax`:

```python
ndump = 10

def rep( sim ):
  if (sim.n % ndump == 0 ):
    sim.emf.report("E", 0 )

# ...
sim = zpic.Simulation( ..., report = rep )
sim. run( tmax )
```

Alternatively, we could create a special simulation loop in Python, and advance the simulation using the `Simulation.iter()` method:

```python
sim = zpic.Simulation(...)

while sim.t <= tmax :
  print( "Now at t = {:g}".format(sim.t))
  if ( sim.n % ndump == 0 ):
    sim.emf.report("E", 0 )
  sim.iter()
```
