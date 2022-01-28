---
title: Spectral field solver
description: Choosing the spectral field solver type
permalink: /documentation/python/spectralsolver
usemathjax: true

layout: single
toc: true
toc_label: Spectral field solver

sidebar:
  nav: "docs"
---

## Introduction

The ZPIC `em1ds` and `em2ds` use spectral based methods to avance EM fields in time. In these codes, users can choose between the so called Pseudo Spectral Time Domain (PSTD) and the Pseudo Spectral Analytical Time Domain (PSATD).

## Choosing the spectral field solver

Users can choose the EM field solver used in spectral codes by setting the `EMF.solver_type` property. Valid values are:

* `"PSTD"` (default) - Use the Pseudo Spectral Time Domain method
* `"PSATD"` - Use the Pseudo Spectral Analytic Time Domain method

The solver type should be set after the simulation object has been created e.g.:

```python
sim = em1ds.Simulation( ... )

# Use PSATD field solver
sim.emf.solver_type = "PSATD"

# Run the simulation
sim.run(tmax)
```
