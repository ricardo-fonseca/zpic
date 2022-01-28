---
title: External EM fields
description: Setting external electro-magnetic fields for you simulation
usemathjax: true
permalink: /documentation/python/emext

layout: single
toc: true
toc_label: External EM fields

sidebar:
  nav: "docs"
---

## Introduction

The ZPIC `em1d` and `em2d` codes allow the user to specify external EM fields to be used in the simulations. These fields, that are constant, will be superimposed on the self-consistent simulation EM fields and used for particle advance only; they do not contribute directly to the evolution of the simulation fields.

Please note that the usual field diagnostics will only report on the self consistent fields. If you want to the look at the full (self consistent + external) fields use the `sim.emf.*_part` properties, e.g. `sim.emf.Bz_part`. See below for more details.

## Setting the external EM fields

Adding an external field is done using `ExternalField` objects and the `emf.set_ext_fld()` method as exemplified below. The following types of external fields are supported:

* "**none**" - no external field
* "**uniform**" - spatially uniform field
* "**custom**" - field defined by a custom function

The external electric and magnetic fields can be set separately.

The `ExternalField` class constructor uses the following parameters:

* `E_type` / `B_type` - Type of initial electric / magnetic field, must be one of:
  * "none" (default) - No initial field
  * "uniform" - Uniform initial field
  * "custom" - Custrom initial field
* `E_0` / `B_0` - Value of the initial electric / magnetic field when using the "uniform" type
* `E_custom` / `B_custom` - Function defining the initial electric / magnetic field when using the "custom" type

The following example sets an external uniform magnetic field $B_z = 1$ (in simulation units):

```python
sim = em1d.Simulation( ... )
ext = em1d.ExternalField(B_type = 'uniform', B_0 = [0.0,0.0,1.0])
sim.emf.set_ext_fld( init )
```

## Setting "custom" external EM fields

The "custom" functions defining initial fields must take exactly 2 arguments in 1D. These arguments are:

* `ix` - the cell index for which the field is being calculated
* `dx` - the cell size in simulation units

This function must return a 3 value list of all the field components. Also note that the position of all field components is not the same inside the cell, and that should be taken into account when calculating the field values.

The following example creates an external sinusoidal field $B_z = \sin(x)^2$:

```python
def sin2( ix, dx ):
    # Bz is located at the center of the cell
    x = (ix+0.5)*dx
    return [0,0,np.sin(x)**2]

ext = em1d.ExternalField(B_type = 'custom', B_custom = sin2 )
```

In 2D the function takes exactly 4 arguments:

* `ix` - the x cell index for which the field is being calculated
* `dx` - the x cell size in simulation units
* `iy` - the y cell index for which the field is being calculated
* `dy` - the y cell size in simulation units

Just like in 1D, this function must return a 3 value list of all the field components, and the position of field components inside the cell should be taken into account when calculating the field values.

The following example creates a circular magnetic field around the center of the box, falling off with $r^2$:

```python
def extB(ix,dx,iy,dy):
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

init = em2d.ExternalField(B_type = 'custom', B_custom = initB )
```

## Self-consistent EM fields and total EM fields

As mentioned above, the external fields are super-imposed (added) to the simulation self-consistent fields at particle push time, directly affecting particle motion. The field advance does not use the external fields, and the self-consistent fields are kept separately.

To simplify diagnostics, ZPIC allows access to both the self-consistent fields and the total fields as seen by the particles. These can be accessed through the following properties:

* `EMF.Ex`, `EMF.Ey`, `EMF.Ez` - Self-consistent simulation electric fields
* `EMF.Bx`, `EMF.By`, `EMF.Bz` - Self-consistent simulation magnetic fields
* `EMF.Ex_part`, `EMF.Ey_part`, `EMF.Ez_part` - Total electric field (self-consistent plus external)
* `EMF.Bx_part`, `EMF.By_part`, `EMF.Bz_part` - Total magnetic field (self-consistent plus external)
