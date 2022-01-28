---
title: Simulation data
description: Accessing simulation data from Python
usemathjax: true
permalink: /documentation/python/data

layout: single
toc: true
toc_label: Accessing simulation data

sidebar:
  nav: "docs"
---

The Python interface allows users to access simulation data directly, without requiring file output. All simulation data is exposed as data members of the `Simulation` object that was used.

## Simulation members

* `n` (read only) - Current iteration number
* `t` (read only)- Current simulation time
* `emf` - EM fields (EMF object)
* `current` - Electric current density (Current object)
* `charge` (spectral codes only) - Charge density (Charge object)
* `species` - Particle species data (list of Species objects)

The iteration and simulation time properties are read-only, and are updated when calling the `iter()` or `run()` methods of the `Simulation` object. 

## Electric and magnetic fields

The raw electric field and magnetic fields of a `Simulation` object `sim` are available through the `sim.emf.E[x|y|z]`, `sim.emf.B[x|y|z]` properties, respectively (see below for `es1d`). Each of these properties will be a [nx] ([nx,ny] in 2D) NumPy float32 array that can be used as usual:

```python
sim = em1d.Simulation( 100, 10.0, 0.09, species = ...)

# Run up to t = 20
sim.run(20.0)

# Get Ex field at the center of the box
print( 'Ex[5.0] = {:g}'.format( sim.emf.Ex[50] ) )
```

### Electrostatic code

The electrostatic code (`es1d`) only uses electric field along the $x$ direction. Accessing this field is done through the `sim.field.E` property:

```python
sim = es1d.Simulation( ... )
# (...)
print( 'Ex[5.0] = {:g}'.format(sim.field.E[50]))
```

## Electric current density

Similarly the electric current density can be accessed through the `sim.current.J[x|y|z]` property, e.g.:

```python
sim = em2d.Simulation( [100,100], ... )
# (...)

# Get Jz current density at the center of the box
print( 'Jz[5.0, 5.0] = {:g})'.format( sim.current.Jz[50,50] ) )
```

## Charge density (Spectral codes only)

The spectral codes (`em1ds`, `em2ds` and `es1d`) also define the total electric charge density. This can be accessed through the `sim.charge.rho` property:

```python
sim = em2ds.Simulation( [100,100], ... )
# (...)

# Get charge density at the center of the box
print( 'rho[5.0, 5.0] = {:g})'.format( sim.charge.rho[50,50] ) )
```

## Particle data

Accessing raw particle data using the `particles` property of each Species object. This property is a NumPy array of `t_part` structures containing:

* `ix` - the particle cell
* `x` - the particle position inside the cell normalized to the cell size ( 0 <= x < 1 )
* `ux`, `uy`, `uz` - the particle generalized velocity in each direction

Note that you can use this properties to change the simulation data, e.g., allowing you to implement arbitrary velocity distributions.

```python
#(...)
right = em1ds.Species( "right", m_q, ppc, ufl = ufl, uth = uth )

# Particle data is initialized when initializing the simulation object
sim = em1ds.Simulation( ..., species = right )

#(...)

# Generalized x velocity of the first particle
print( 'ux[0] = {:g}'.format( right.particles['ux'][0]))
```
