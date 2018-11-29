/**
 * ZPIC - em2ds
 *
 * Weibel instability
 */

#include <stdlib.h>
#include <math.h>

#include "../simulation.h"

void sim_init( t_simulation* sim ){

	// Time step
	float dt = 0.07;
	float tmax = 0.7;

	// Simulation box
	int   nx[2]  = { 128, 128 };
	float box[2] = { 12.8, 12.8 };

	// Diagnostic frequency
	int ndump = 1;

    // Initialize particles
	const int n_species = 1;
	t_species* species = (t_species *) malloc( n_species * sizeof( t_species ));

	// Use 2x2 particles per cell
	int ppc[] = {2,2};

	// Initial fluid and thermal velocities
	float ufl[] = { 10000.0, 10000.0, 10000.0 };


	spec_new( &species[0], "positrons", +1.0, ppc, ufl, NULL, nx, box, dt, NULL );

	// Initialize Simulation data
	sim_new( sim, nx, box, dt, tmax, ndump, species, n_species );

}


void sim_report( t_simulation* sim ){

	current_report( &sim ->current, 0 );
	current_report( &sim ->current, 1 );
	current_report( &sim ->current, 2 );

	charge_report( &sim ->charge );

}
