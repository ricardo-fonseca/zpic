## Adding laser pulses

The input parameters may optionally define an arbitrary number of laser pulses to be added to the simulation. Each laser pulse is added through a call to the `sim_add_laser()` routine:

```C
void sim_add_laser( t_simulation* sim, t_emf_laser* laser)
```

This routine should be called inside `sim_init()`, somewhere after the call to `sim_new()`. Laser parameters are defined in the supplied _t\_emf\_laser_ structure.

| Laser parameters| Description|
|---|---|
| type | PLANE for a plane wave or GAUSSIAN for a gaussian beam |
| start | Front edge of the laser pulse |
| fwhm  | FWHM of the laser pulse duration |
| rise, flat, fall  | Rise / flat / fall time of the laser pulse |
| a0  | Normalized peak vector potential of the pulse |
| omega0 | Laser frequency |
| polarization | Laser polarization angle, 0 aligns $E$ field along $x_2$ |
| W0 | Gaussian beam waist |
| focus | Focal plane position |
| axis | Position of the optical axis |

All parameters are in normalized simulation units except for the normalized peak vector potential, $a_0$, which is adimensional. In 1D only plane waves exist, so the _type_, _W0_, _focus_ and _axis_ parameters are not supported.

Using the _fwhm_ parameter will override the _rise_, _flat_ and _fall_ parameters. Specifically, it sets _rise = fwhm/2_, _flat = 0_, and _fall = fwhm/2_.

The following example launches a laser starting at position 17.0, with a (temporal) full width at half max of 2.0. The peak normalized vector potential is 2.0, and the laser frequency is 10.0. The polarization degree is $\pi/2$, which aligns the $E$ field along the $x_3$ direction. All values are in normalized simulation units.

```C
t_emf_laser laser = {
	.start = 17.0,
	.fwhm  = 2.0,
	.a0 = 2.0,
	.omega0 = 10.0,
	.polarization = M_PI_2
};
sim_add_laser( sim, &laser );
```

The following is a 2D example of a gaussian laser pulse. It uses the same parameters as the previous example, set the beam focus waist to 4.0, the focal plane position to x=20.0, and the propagation axis to y=12.8.

```C
t_emf_laser laser = {
	.type = GAUSSIAN,
	.start = 17.0,
	.fwhm  = 2.0,
	.a0 = 2.0,
	.omega0 = 10.0,
	.W0 = 4.0,
	.focus = 20.0,
	.axis = 12.8,
	.polarization = M_PI_2
};
sim_add_laser( sim, &laser );
```

See for example the Laser Wakefield input files (e.g. [lwfa.c](../em1d/input/lwfa.c)).
