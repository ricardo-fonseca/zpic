# The ZPIC educational code suite


## Units

The current implementation uses the same normalized units as the OSIRIS code. This requires choosing a reference normalization frequency, $\omega_n$. Time is normalized to $1/\omega_n$. (Proper) velocities are normalized to the speed of light, $c$. Space is normalized to $c/\omega_n$. The fields are then normalized appropriately.

The density is normalized to $\omega_n^2$ (the normalization frequency squared). So if the density is 1 at a given location then the normalization frequency is the plasma frequency at that location. 

If the laser frequency is 1, then the normalization frequency is the laser frequency and the density is normalized to the critical densify (for that laser frequency)

|       Normalized units                   |
| ---------------------------------------- |
| $x' = \frac{\omega_n}{c} x$|
| $v' = \frac{v}{c}$ |
| $u' = \frac{u}{c} = \frac{\gamma v}{c}$ |
| $E' = e \frac{c / \omega_n }{m_e c^2} E$ |
| $B' = e \frac{c / \omega_n }{m_e c^2} B$ |

## Initializing simulation particles

...

## Moving simulation window

The finite difference models can be run using a moving simulation window that moves at the speed of light. Please note that simulation is done in the lab reference frame; it is the simulation box that moves and follows relevant phenomena moving at, or close to, this speed, such as laser pulses or relativistic particle beams.

Using a moving window requires calling the sim\_set\_moving\_window() routine:

```C
void sim_set_moving_window( t_simulation* sim )
```

This routine should be called inside sim\_init(), somewhere after the call to sim\_new(). Spectral models currently do not support this feature.

## Adding laser pulses

Adding laser pulses to the simulation is done through a call to the sim\_add\_laser() routine:

```C
void sim_add_laser( t_simulation* sim, t_emf_laser* laser)
```

This routine should be called inside sim\_init(), somewhere after the call to sim\_new(). Laser parameters are defined in the supplied _t\_emf\_laser_ structure.

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
 
All parameters are in normalized simulation units except for the normalized peak vector potential, $a_0$. In 1D only plane waves exist, so the _type_, _W0_, _focus_ and _axis_ parameters are not supported.

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

## Neutralizing backgound

### Finite difference models

When using the finite difference models, the code will always behave as if the total charge density was initially zero, even when this is not the case. This effect is generally described as having the code add a neutralizing background; however, for finite difference models, this is simply a consequence of the way the field solver is implemented, and no actual addition takes place. So to model a plasma with a fixed ion background, it would suffice to add the electron species, and no additional changes are required.

### Spectral models

For the spectral models this is no longer the case, and if required a neutralizing background must be explicitly added. This is achieved through a call to the sim\_add\_neutral\_bkg() routine:

```C
void sim_add_neutral_bkg( t_simulation* sim )
```

This routine should be called inside sim\_init(), somewhere after the call to sim\_new().

The code will calculate the charge density resulting from the initial charge distribution of the multiple species and store that value on an auxiliary buffer. This buffer will then be subtracted from the simulation charge density calculated at each time step, before calling the field advance.

# Algorithm details

## Time centering of quantities

Time integration of quantities in these codes is done using a leap frog method. For this purpose particle positions are defined at integral time steps, $t_i$, $t_{i+1}$, etc., whereas velocities (or generalized velocities in the relativistic versions) are defined at half time steps $t_{i-1/2}$, $t_{i+1/2}$, etc.

As a consequence, the electromagnetic fields, $E$ and $B$, as well as charge density, $\rho$, are also defined at integral timesteps $t_i$. Energy will also be defined at integral time steps. For this purpose particle energy is calculated during the particle advance, using time centered values of the velocity.

Current density, $j$, is defined at half time steps $t_{i-1/2}$, $t_{i+1/2}$.
