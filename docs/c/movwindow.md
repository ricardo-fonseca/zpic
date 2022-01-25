## Moving simulation window

The finite difference models can be run using a moving simulation window that moves at the speed of light. Please note that simulation is done in the lab reference frame; it is the simulation box that moves and follows relevant phenomena moving at, or close to, this speed, such as laser pulses or relativistic particle beams.

Using a moving window requires calling the `sim_set_moving_window()` routine:

```C
void sim_set_moving_window( t_simulation* sim )
```

This routine should be called inside `sim_init()`, somewhere after the call to `sim_new()`. Spectral models currently do not support this feature.

See for example the Laser Wakefield input files (e.g. [lwfa.c](../em1d/input/lwfa.c)).
