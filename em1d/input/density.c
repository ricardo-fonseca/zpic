/**
 * ZPIC - em1d
 *
 * Initializing different density profiles
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "../simulation.h"

float custom_n0( float x, void *data ) {

	return 1.0 + 0.5*sin(x/M_PI)*sin(x/M_PI);

}

void sim_init( t_simulation* sim ){

	// Time step
	float dt = 0.019;
	float tmax = 10.00;

	// Simulation box
	int   nx  = 64;
	float box = 20.0;

	// Diagnostic frequency
	int ndump = 100;

    // Initialize particles
	const int n_species = 1;

	// Use 128 particles per cell
	int ppc = 128;

	// Density profile
//	t_density density = { .type = UNIFORM };
//	t_density density = { .type = STEP, .start = 17.5 };
//	t_density density = { .type = SLAB, .start = 17.5, .end = 22.5 };
//	t_density density = { .type = RAMP, .start = 17.5, .end = 22.5, .ramp = { 1.0, 2.0 } };
	t_density density = { .type = CUSTOM, .custom = &custom_n0 };

	t_species* species = (t_species *) malloc( n_species * sizeof( t_species ));
	spec_new( &species[0], "electrons", -1.0, ppc, NULL, NULL, nx, box, dt, &density );

	// Initialize Simulation data
	sim_new( sim, nx, box, dt, tmax, ndump, species, n_species );

	// Set moving window (this must come after sim_new)
	sim_set_moving_window( sim );

}


void sim_report( t_simulation* sim ){

	// Charge density
	spec_report( &sim->species[0], CHARGE, NULL, NULL );

	// RAW dump
	spec_report( &sim->species[0], PARTICLES, NULL, NULL );

}
