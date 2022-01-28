---
title: Getting started
description: Getting started with ZPIC simulations in Python
permalink: /documentation/python/start
usemathjax: true

layout: single
toc: true
toc_label: Getting started

sidebar:
  nav: "docs"
---

## Choosing the code version

To run a ZPIC simulation in Python (including Jupyter notebooks) users must first choose which code version is to be used and import the appropriate module. Available versions are:

* `em1d`, `em1ds` - 1D electromagnetic code using a finite difference field solver (`em1d`) or a spectral field solver (`em1ds`)
* `em2d`, `em2ds` - 2D electromagnetic code using a finite difference field solver (`em2d`) or a spectral field solver (`em2ds`)
* `es1d` - 1D electrostatic code using a spectral field solver

To run the simulation users must start by creating a `Simulation` object from the appropriate module, e.g.:

```python
import em1d
sim = em1d.Simulation( ... )
```

## Setting the simulation parameters

Simulation parameters are set at the creation of the `Simulation` object. The `Simulation` class constructor has the following syntax:

```python
Simulation( nx, box, dt, species = [list], report = function )
```

The parameters are the following:

* `nx` - Number of grid cells to use
* `box`- Simulation box (phyiscal) size, in simulation units
* `dt` - Simulation time step, in simulation units
* `species` (optional) - List of particle species to use in the simulation
* `report` (optional) - Python function to call at every time-step when using the `run()` method.

The input parameters must define the number of grid points in each direction of the simulation grid, and the physical dimensions (in simulation units) of this box. These dimensions are measured from the lower boundary of the first cell to the upper boundary of the last cell. The input parameters must also define the time step (in simulation units)

### Example

The following example sets the minimal parameters for a 2D ZPIC simulation with no particles:

```python
import em2d

nx = [100,100]          # 100 x 100 grid
box = [10.0, 10.0]      # 10.0 x 10.0 simulation box
dt = 0.07               # 0.07 time step

sim = em2d.Simulation( nx, box, dt )
```

## Adding particles

Adding particle species to the simulation is done using the `species` keyword of the `Simulation` class. Users may define an an arbitrary number of particle species, defining the required plasma properties, such as density, temperature, fluid velocity, etc. See the [particles](particles) section for details.

The following example initializes a 1D simulation with a warm plasma:

```python
import em1ds

plasma = em1ds.Species("plasma", -1.0, 500, uth = [0.001, 0.001, 0.001] )
sim = em1ds.Simulation( 120, 4 * np.pi, 0.1, species = plasma )
```

## Running the simulation

The `Simulation` class provides 2 methods for advancing the simulation in time:

* `run(tmax)` - Advances ths simulation up to time `tmax`. If the simulation is already at `tmax` or later a warning message is issued. If the user specified a `report` function when creating the simulation object, this function will be called once for every time step, mimicking the behaviour of the C codes.
* `iter()` - Advances the simulation 1 iteration. This gives more control over the simulation, giving access to simulation data at every time-step, and allowing for additional actions to be performed.

The 2 following examples run the simulation up to `t = 10.0`

```python
import em1ds

sim = em1ds.Simulation( 120, 4 * np.pi, 0.1 )

# Option 1 - use the run() method
# sim.run(10.0)

# Option 2 - use the iter() method
while sim.t < 10.0:
    print('n = {:d}, t = {:g}'.format(sim.n,sim.t) )
    sim.iter()
```
