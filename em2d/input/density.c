/**
 * ZPIC - em2d
 *
 * Particle density initialization options
 */

#include <stdlib.h>
#include "../simulation.h"

float custom_x( float x, void *data ) {
	return exp( -(x-6.4)*(x-6.4)/3.0);
}

float custom_y( float y, void *data ) {
	return exp( -(y-6.4)*(y-6.4)/1.0);
}

void sim_init( t_simulation* sim ){

	// Time step
	float dt = 0.07;
	float tmax = 35.0;

	// Simulation box
	int   nx[2]  = { 128, 128 };
	float box[2] = { 12.8, 12.8 };

	// Diagnostic frequency
	int ndump = 100;

    // Initialize particles
	const int n_species = 1;

	// Use 2x2 particles per cell
	int ppc[] = {8,8};

    // Density profile
//	t_density density = { .type = UNIFORM };
//	t_density density = { .type = STEP, .start = 10.0 };
//	t_density density = { .type = SLAB, .start = 10.0, .end = 15.6 };
    t_density density = { 
        .type = CUSTOM, 
        .custom_x = &custom_x,
        .custom_y = &custom_y };

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
