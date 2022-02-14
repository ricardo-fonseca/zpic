/**
 * ZPIC - em1ds
 *
 * Laser Wakefield Acceleration
 */

#include <stdlib.h>
#include <math.h>

#include "../simulation.h"

void sim_init( t_simulation* sim ){

	// Time step
	float dt = 0.1;
	float tmax = 50.0;

	// Simulation box
	int   nx  = 120;
	float box = 4*M_PI;

	// Diagnostic frequency
	int ndump = 10;

    // Initialize particles
	const int n_species = 2;

	// Use 1000 particles per cell
	int ppc = 500;

	t_species* species = (t_species *) malloc( n_species * sizeof( t_species ));


	// Initial fluid and thermal velocities
	float ufl[] = { 0.2 , 0.0 , 0.0 };
	float uth[] = { 0.001 , 0.001 , 0.001 };


	spec_new( &species[0], "right", -1.0, ppc,
		      ufl, uth, nx, box, dt, NULL );

	ufl[0] = -ufl[0];
	spec_new( &species[1], "left", -1.0, ppc,
		      ufl, uth, nx, box, dt, NULL );

	// Initialize Simulation data
	sim_new( sim, nx, box, dt, tmax, ndump, species, n_species );

}


void sim_report( t_simulation* sim ){

	// Charge report - global
	charge_report( &sim->charge );

	spec_report( &sim->species[0], PARTICLES, NULL, NULL );
	spec_report( &sim->species[1], PARTICLES, NULL, NULL );

	spec_report( &sim->species[0], CHARGE, NULL, NULL );

    // x1 electric field component
	emf_report( &sim->emf, EFLD, 0 );

}
