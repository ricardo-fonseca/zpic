---
title: Initial EM fields
description: Setting the initial electro-magnetic fields for your simulation
usemathjax: true
permalink: /documentation/python/eminit

layout: single
toc: true
toc_label: Initial EM fields

sidebar:
  nav: "docs"
---

## Introduction

The EM codes in ZPIC work by advancing the EM fields in time using the electric current density from particle motion as a source term. The algorithm implcitly assumes that the initial fields in the simulation are self-consistent, and they will be evolved accordingly.

## Setting the initial EM fields

By default, the initial EM fields in ZPIC codes are set to 0. Other initial values can be set using `InitialField` objects and the `EMF.init_fld()`. The following types of initial fields are supported:

* "**none**" - no initial field
* "**uniform**" - spatially uniform field
* "**custom**" - field defined by a custom function

The initial electric and magnetic fields can be set separately. Please note that the method can only be called before the simulation starts. Currently only EM1D and EM2D support initial fields, adding this feature to the remaining codes is in the pipeline.

The `InitialField` class constructor uses the following parameters:

* `E_type` / `B_type` - Type of initial electric / magnetic field, must be one of:
  * "none" (default) - No initial field
  * "uniform" - Uniform initial field
  * "custom" - Custrom initial field
* `E_0` / `B_0` - Value of the initial electric / magnetic field when using the "uniform" type
* `E_custom` / `B_custom` - Function defining the initial electric / magnetic field when using the "custom" type

The following example sets an initial uniform magnetic field $B_z = 1$ (in simulation units):

```python
sim = em1d.Simulation( ... )
init = em1d.InitialField(B_type = 'uniform', B_0 = [0.0,0.0,1.0])
sim.emf.init_fld( init )
```

## Setting "custom" initial fields

The "custom" functions defining initial fields must take exactly 2 arguments in 1D. These arguments are:

* `ix` - the cell index for which the field is being calculated
* `dx` - the cell size in simulation units

This function must return a 3 value list of all the field components. Also note that the position of all field components is not the same inside the cell, and that should be taken into account when calculating the field values.

The following example creates an initial sinusoidal field $B_z = \sin(x)^2$:

```python
def sin2( ix, dx ):
    # Bz is located at the center of the cell
    x = (ix+0.5)*dx
    return [0,0,np.sin(x)**2]

init = em1d.InitialField(B_type = 'custom', B_custom = sin2 )
```

In 2D the function takes exactly 4 arguments:

* `ix` - the x cell index for which the field is being calculated
* `dx` - the x cell size in simulation units
* `iy` - the y cell index for which the field is being calculated
* `dy` - the y cell size in simulation units

Just like in 1D, this function must return a 3 value list of all the field components, and the position of field components inside the cell should be taken into account when calculating the field values.

The following example creates a circular magnetic field around the center of the box, falling off with $r^2$:

```python
def initB(ix,dx,iy,dy):
    x0 = 6.4
    y0 = 6.4
    
    x = ix*dx       - x0;
    y = (iy+0.5)*dy - y0;

    r2 = x*x+y*y;
    bx = -y/r2;

    x = (ix+0.5)*dx - x0;
    y = iy*dy       - y0;

    r2 = x*x+y*y;
    by = x/r2;

    return [bx,by,0] 

init = em2d.InitialField(B_type = 'custom', B_custom = initB )
```
