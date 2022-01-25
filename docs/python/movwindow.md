---
layout: page
title: Moving simulation window
permalink: /python/mov_window
---

The finite difference models can be run using a moving simulation window that moves at the speed of light. Please note that simulation is done in the lab reference frame; it is the simulation box that moves and follows relevant phenomena moving at, or close to, this speed, such as laser pulses or relativistic particle beams.

To use a moving window, you must call the `Simulation.set_moving_window()` method after creating the simulation object, e.g.:

```python
import em1d

sim = em1d.Simulation( 1000, 100.0, 0.09, species = ... )

sim.set_moving_window()
```
