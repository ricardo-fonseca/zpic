---
layout: single
title: ZPIC codes
permalink: /about/implementation

sidebar:
  nav: "about"
---

The ZPIC project implements a suite of 1D/2D fully relativistic electromagnetic PIC codes, as well as 1D electrostatic, geared towards plasma physics education and research. For a brief discussion about the particle-in-cell algorithm, check the [PIC algorithm section](pic).

## The ZPIC codes

The ZPIC project consists of 5 different PIC codes, implementing different variants of the PIC algorithm:

* __em1d__ and __em2d__ - Implement fully-relativistic electro-magnetic particle-in-cell codes in 1D and 2D, using a finite difference field solver.
* __em1d__ and __em2ds__ - Implement fully-relativistic electro-magnetic particle-in-cell codes in 1D and 2D, using a spectral field solver.
* __es1d__ - Implements an electro-static particle-in-cell code in 1D.

All codes work in cartesian geometry, meaning that particles will behave as infinite planes (1D codes) or infinite wires (2D codes). No other geometry is currently supported.

## The ZPIC implementation

As in any PIC code, the main algorithm of any ZPIC code can essencially be divided into 4 parts. See below for the ZPIC implementation details on each of these parts:

1. [Particle advance](particles)
2. [Current/charge deposition](deposit)
3. [Field advance](fields)
4. [Field interpolation](interpolation)

Also check the following sections for details on the units used in ZPIC codes, as well as aditional implementation details.

* [ZPIC units](units)
* [Time-centering of quantities](timecenter)
* [Data structures](data)
